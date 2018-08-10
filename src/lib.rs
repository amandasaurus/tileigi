#![allow(dead_code,unused_imports,unused_variables,unused_mut,unused_assignments)]
#![feature(extern_prelude)]

#[macro_use(to_sql_checked)]
extern crate postgres;
extern crate slippy_map_tiles;
extern crate yaml_rust;
extern crate mapbox_vector_tile;
extern crate wkb;
extern crate geo;
#[macro_use]
extern crate serde_json;

#[macro_use]
extern crate log;

extern crate users;
extern crate num_traits;
extern crate rusqlite;
extern crate md5;
extern crate byteorder;
extern crate separator;
extern crate procinfo;

use std::fs::File;
use std::fs;
use std::io::prelude::*;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::PathBuf;
use std::collections::{HashSet, HashMap};
use std::time::Instant;
use std::borrow::{Cow, Borrow};
use std::rc::Rc;
use std::fmt::Write;

use std::thread;
use std::sync::mpsc::{channel, sync_channel, Sender, SyncSender};
use std::sync::{Arc, Mutex};

use yaml_rust::{YamlLoader, Yaml};

use postgres::{Connection, TlsMode};
use postgres::params::ConnectParams;
use postgres::types::{Type, ToSql, IsNull};

use slippy_map_tiles::{BBox, Metatile, MetatilesIterator};

use geo::*;
use geo::algorithm::map_coords::MapCoords;
use geo::algorithm::map_coords::MapCoordsInplace;
use geo::algorithm::boundingbox::BoundingBox;
use geo::algorithm::contains::Contains;
use separator::Separatable;

mod clip;
use clip::{clip_to_bbox,clip_geometry_to_tiles};

mod validity;
use validity::{is_valid, is_valid_skip_expensive};

macro_rules! memory {
    () => (
        use separator::Separatable;
        let total_memory = procinfo::pid::status_self().ok().map(|s| s.vm_rss.separated_string()).unwrap_or("N/A".to_string());
        info!("total memory {} KiB", total_memory);
    );
    ($msg: expr) => (
        use separator::Separatable;
        let total_memory = procinfo::pid::status_self().ok().map(|s| s.vm_rss.separated_string()).unwrap_or("N/A".to_string());
        info!("total memory {} KiB: {}", total_memory, $msg);
    );
    ($msg: expr, $($arg:tt)*) => (
        use separator::Separatable;
        let total_memory = procinfo::pid::status_self().ok().map(|s| s.vm_rss.separated_string()).unwrap_or("N/A".to_string());
        info!("total memory {} KiB: {}", total_memory, format!($msg, $($arg)*));
    );
}

mod printer;
mod fileio;
mod simplify;

use fileio::{FileIOMessage,TileDestination};

mod stringstore;
use stringstore::StringStore;

#[cfg(test)]
mod test;

#[derive(Clone)]
pub enum TileDestinationType {
    TileStashDirectory(PathBuf),
    MBTiles(PathBuf),
    ModTileDirectory(PathBuf),
}

pub struct ConnectionPool {
    connections: HashMap<ConnectParams, Connection>,
    layer_to_param: HashMap<String, ConnectParams>,

}

impl ConnectionPool {
    pub fn new(params_to_layers: HashMap<ConnectParams, Vec<String>>) -> Self {
        let mut layer_to_param = HashMap::new();
        for (cp, ls) in params_to_layers.iter() {
            for l in ls.iter() {
                layer_to_param.insert(l.clone(), cp.clone());
            }
        }

        let mut connections = HashMap::with_capacity(params_to_layers.len());
        for (cp, layers) in params_to_layers.into_iter() {
            let connection = Connection::connect(cp.clone(), TlsMode::None).unwrap();
            connections.insert(cp, connection);
        }

        ConnectionPool{ connections: connections, layer_to_param: layer_to_param }
    }

    fn connection_for_layer<'a>(&'a self, layer_id: &str) -> &'a Connection {
        let cp = &self.layer_to_param[layer_id];
        
        &self.connections[cp]
    }

}

#[derive(Clone,Debug)]
struct TableSQL {
    query: String,
    has_pixel_width: bool,
    has_pixel_height: bool,
    has_scale_denominator: bool,
        
}

impl TableSQL {
    fn new(query: String) -> Self {
        let has_pixel_width = query.contains("!pixel_width!");
        let has_pixel_height = query.contains("!pixel_height!");
        let has_scale_denominator = query.contains("!scale_denominator!");

        let mut query = query;

        query = query.replace("!bbox!", "$1");

        let mut param_num = 2;
        if has_pixel_width {
            query = query.replace("!pixel_width!", &format!("${}", param_num));
            param_num += 1;
        }
        if has_pixel_height {
            query = query.replace("!pixel_height!", &format!("${}", param_num));
            param_num += 1;
        }
        if has_scale_denominator {
            query = query.replace("!scale_denominator!", &format!("${}", param_num));
            param_num += 1;
        }
                
        let query = format!("SELECT ST_AsBinary(way), * from {} where way && $1", query);
        TableSQL{
            query, has_pixel_width, has_pixel_height, has_scale_denominator,
        }
    }

    fn params<'a, T: num_traits::Float+Into<f64>+'a+std::fmt::Debug>(&self, bbox: &'a LocalBBox<T>, pixel_width: &'a f32, pixel_height: &'a f32, scale_denominator: &'a f32) -> Vec<&'a postgres::types::ToSql> {
        // we always have bbox
        let mut results: Vec<&postgres::types::ToSql> = Vec::with_capacity(4);
        results.push(bbox);    // bbox
        if self.has_pixel_width {
            results.push(pixel_width);
        }
        if self.has_pixel_height {
            results.push(pixel_height);
        }
        if self.has_scale_denominator {
            results.push(scale_denominator);
        }

        results
    }
}

#[derive(Clone,Debug)]
pub struct Layers {
    layers: Vec<Layer>,
    global_maxzoom: u8,
    global_minzoom: u8,
}

#[derive(Clone,Debug)]
struct Layer {
    minzoom: u8,
    maxzoom: u8,
    buffer: u16,
    id: String,
    table: TableSQL,
    dbname: Option<String>,
}

impl Layers {
    fn from_file(filename: &str) -> Self {
        let mut file = File::open(filename).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();

        let mut data_yml = YamlLoader::load_from_str(&contents).unwrap();
        let data_yml = data_yml.remove(0);

        let global_maxzoom = data_yml["maxzoom"].as_i64().unwrap() as u8;
        let global_minzoom = data_yml["minzoom"].as_i64().unwrap() as u8;

        // rust-yaml really needs an into_hash (etc)
        // clone all the things
        let layers = data_yml.as_hash().unwrap().clone().remove(&Yaml::String("Layer".to_string())).unwrap();
        let layers = layers.as_vec().unwrap();

        let layers: Vec<Layer> = layers.into_iter().filter(
            |l| l["Datasource"]
                .as_hash()
                .and_then(
                    |d| d.get(&Yaml::String("type".to_string()))
                        .and_then(|t| Some(t.as_str() == Some("postgis")))
                )
                .unwrap_or(false)
            )
            .map(|layer| {
                let table = layer["Datasource"]["table"].as_str().unwrap();
                let table = TableSQL::new(table.to_owned());
                
                Layer {
                    id: layer["id"].as_str().unwrap().to_owned(),
                    dbname: layer["Datasource"]["dbname"].as_str().map(|x| x.to_owned()),
                    minzoom: layer["properties"]["minzoom"].as_i64().map(|x| x as u8).unwrap_or(global_minzoom) as u8,
                    maxzoom: layer["properties"]["maxzoom"].as_i64().map(|x| x as u8).unwrap_or(global_maxzoom) as u8,
                    buffer: layer["properties"]["buffer-size"].as_i64().map(|x| x as u16).unwrap_or(0) as u16,
                    table: table,
                }
            })
            .map(|layer| {
                layer
            })
            .collect();

        Layers{ layers: layers, global_minzoom: global_minzoom, global_maxzoom: global_maxzoom }

    }

    fn get_all_connections(&self) -> HashMap<ConnectParams, Vec<String>> {
        let mut conns = HashMap::new();
        for layer in self.layers.iter() {

            let mut conn_params = postgres::params::Builder::new();
            // TODO do others
            if let Some(ref dbname) = layer.dbname {
                conn_params.database(&dbname);
            }

            if let Some(username) = users::get_current_username() {
                conn_params.user(&username, None);
            }
            // TODO read hosts
            let conn_params = conn_params.build(postgres::params::Host::Tcp("localhost".to_string()));

            if ! conns.contains_key(&conn_params) {
                conns.insert(conn_params.clone(), Vec::new());
            }

            let layer_id: String = layer.id.clone();
            conns.get_mut(&conn_params).unwrap().push(layer_id);
            
        }

        conns
    }
}

#[inline]
fn fmt_duration(dur: &std::time::Duration) -> String {
    format!("{:.2}s", duration_to_float_secs(dur))
}

#[inline]
fn duration_to_float_secs(dur: &std::time::Duration) -> f64 {
    (dur.as_secs() as f64) + (dur.subsec_nanos() as f64 / 1e9)
}

fn scale_denominator_for_zoom(zoom: u8) -> f32 {
    match zoom {
        0 => 250000000000.,
        1 => 500000000.,
        2 => 200000000.,
        3 => 100000000.,
        4 => 50000000.,
        5 => 25000000.,
        6 => 12500000.,
        7 => 6500000.,
        8 => 3000000.,
        9 => 1500000.,
        10 => 750000.,
        11 => 400000.,
        12 => 200000.,
        13 => 100000.,
        14 => 50000.,
        15 => 25000.,
        16 => 12500.,
        17 => 5000.,
        18 => 2500.,
        _ => {
            eprintln!("Unsupported zoom ({})", zoom);
            unimplemented!();
        }
    }
}

pub fn generate_all(filename: &str, min_zoom: u8, max_zoom: u8, bbox: &Option<BBox>, dest: &TileDestinationType, if_not_exists: bool, compress: bool, metatile_scale: u8, num_threads: usize, tile_list: Option<String>, file_writer_buffer: usize, quiet: bool) {
    let layers = Layers::from_file(filename);

    let connection_pool = ConnectionPool::new(layers.get_all_connections());

    let (metatile_iterator, total_num_of_metatiles) = match tile_list {
        None => {
            let total_num_of_metatiles: Option<usize> = (min_zoom..max_zoom+1).map(|z| {
                match *bbox {
                    None => {
                        let scale = (metatile_scale.trailing_zeros()) as u32;
                        let z = z as u32;
                        if scale >= z {
                            Some(1)
                        } else {
                            Some(2_u64.pow(2u32 * (z-scale)) as usize)
                        }
                    },
                    Some(ref bbox) => {
                        let this = slippy_map_tiles::size_bbox_zoom_metatiles(&bbox, z, metatile_scale);
                        this
                    },
                }})
                .fold(Some(0_usize), |acc, on_this_zoom| {
                    match (acc, on_this_zoom) {
                        (Some(x), Some(y)) => x.checked_add(y),
                        _ => None,
                    }
                });

            let metatile_iterator = MetatilesIterator::new_for_bbox_zoom(metatile_scale, &bbox, min_zoom, max_zoom);

            (metatile_iterator, total_num_of_metatiles)
        },
        Some(tile_list) => {
            let mt_file_list = MetatilesIterator::new_from_filelist(tile_list);
            let total_num_of_metatiles = mt_file_list.total();
            (mt_file_list, total_num_of_metatiles)
        }
    };

    let metatile_iterator = Arc::new(Mutex::new(metatile_iterator));


    let (printer_tx, printer_rx) = channel();
    let new_bbox: Option<BBox> = bbox.clone();
    let mut printer_thread = if quiet {
        thread::spawn(move || { printer::quiet_printer(printer_rx, total_num_of_metatiles) })
    } else {
        thread::spawn(move || { printer::printer(printer_rx, total_num_of_metatiles) })
    };

    let (fileio_tx, fileio_rx) = sync_channel(file_writer_buffer);

    let mut fileio_thread = match dest {
        &TileDestinationType::TileStashDirectory(ref path) => {
            let tile_dest = fileio::TileStashDirectory::new(&path);
            write_tilejson(&layers, &connection_pool, &path);
            thread::spawn(move || { fileio::fileio_thread(fileio_rx, Box::new(tile_dest)) })
        },
        &TileDestinationType::MBTiles(ref path) => {
            let mut tile_dest = fileio::MBTiles::new(&path);
            tile_dest.set_tilejson_vector_layers(tilejson_vector_layers(&layers, &connection_pool));
            thread::spawn(move || { fileio::fileio_thread(fileio_rx, Box::new(tile_dest)) })
        },
        &TileDestinationType::ModTileDirectory(ref path) => {
            write_tilejson(&layers, &connection_pool, &path);
            let tile_dest = fileio::ModTileMetatileDirectory::new(&path);
            thread::spawn(move || { fileio::fileio_thread(fileio_rx, Box::new(tile_dest)) })
        },
    };



    let mut workers = Vec::with_capacity(num_threads);
    for _ in 0..num_threads {
        // TODO do I need all these clones?
        let my_connection_pool = ConnectionPool::new(layers.get_all_connections());
        let my_printer_tx = printer_tx.clone();
        let my_fileio_tx = fileio_tx.clone();
        let my_metatile_iterator = Arc::clone(&metatile_iterator);
        let my_layers = layers.clone();
        let my_dest = dest.clone();

        let should_do_metatile = move |mt: &slippy_map_tiles::Metatile| {
            if if_not_exists {
                match my_dest {
                    TileDestinationType::TileStashDirectory(ref path) => {
                        !fileio::TileStashDirectory::does_metatile_exist(&path, &mt)
                    },
                    TileDestinationType::ModTileDirectory(ref path) => {
                        !fileio::ModTileMetatileDirectory::does_metatile_exist(&path, &mt)
                    }
                    TileDestinationType::MBTiles(ref path) => {
                        unimplemented!();
                    },
                }
            } else {
                true
            }
        };

        let handle = thread::spawn(move || {
            worker_all_layers(my_printer_tx, my_fileio_tx, my_metatile_iterator, &my_connection_pool, &my_layers, should_do_metatile);
        });
        workers.push(handle);
    }

    for worker in workers {
        // If one of our worker threads has panic'ed, then this main programme should fail too
        worker.join().ok();
    }

    printer_tx.send(printer::PrinterMessage::Quit).unwrap();
    printer_thread.join().unwrap();
    fileio_tx.send(FileIOMessage::Quit).unwrap();

    if ! quiet {
        println!("All tiles generated. Waiting for all to be written to disk...");
    }

    fileio_thread.join().unwrap();

    memory!("Finished");
    if ! quiet {
        println!("Finished.");
    }

}

fn worker_all_layers<F>(printer_tx: Sender<printer::PrinterMessage>, fileio_tx: SyncSender<FileIOMessage>, mut metatile_iterator: Arc<Mutex<Iterator<Item=Metatile>>>, connection_pool: &ConnectionPool, layers: &Layers, should_do_metatile: F)
    where F: Fn(&slippy_map_tiles::Metatile) -> bool,
{
    loop {
        let metatile = metatile_iterator.lock().unwrap().next();
        if let None = metatile {
            // The iterator is finished.
            break;
        }
        let metatile = metatile.unwrap();

        if ! should_do_metatile(&metatile) {
            continue;
        }

        let tiles = single_metatile(&layers, &metatile, &connection_pool);
        let num_tiles = tiles.len();

        let tiles: Vec<_> = tiles.into_iter().map(|(tile, mvt)| (tile, mvt.to_compressed_bytes())).collect();

        printer_tx.send(printer::PrinterMessage::DoneTiles(metatile.zoom(), 1, num_tiles)).unwrap();

        fileio_tx.send(FileIOMessage::SaveMetaTile(metatile, tiles)).unwrap();

    }

}

pub fn generate_by_layer(filename: &str, min_zoom: u8, max_zoom: u8, bbox: &Option<BBox>, dest: &TileDestinationType, if_not_exists: bool, compress: bool, metatile_scale: u8, num_threads: usize, tile_list: Option<String>, file_writer_buffer: usize) {
    if tile_list.is_some() {
        unimplemented!();
    }
    let layers = Layers::from_file(filename);

    let connection_pool = ConnectionPool::new(layers.get_all_connections());

    let (fileio_tx, fileio_rx) = sync_channel(file_writer_buffer);


    let mut fileio_thread = match dest {
        &TileDestinationType::TileStashDirectory(ref path) => {
            let tile_dest = fileio::TileStashDirectory::new(&path);
            write_tilejson(&layers, &connection_pool, &path);
            thread::spawn(move || { fileio::fileio_thread(fileio_rx, Box::new(tile_dest)) })
        },
        &TileDestinationType::MBTiles(ref path) => {
            let mut tile_dest = fileio::MBTiles::new(&path);
            tile_dest.set_tilejson_vector_layers(tilejson_vector_layers(&layers, &connection_pool));
            thread::spawn(move || { fileio::fileio_thread(fileio_rx, Box::new(tile_dest)) })
        },
        &TileDestinationType::ModTileDirectory(ref path) => {
            write_tilejson(&layers, &connection_pool, &path);
            let tile_dest = fileio::ModTileMetatileDirectory::new(&path);
            thread::spawn(move || { fileio::fileio_thread(fileio_rx, Box::new(tile_dest)) })
        },
    };
    let global_maxzoom = layers.global_maxzoom;

    for layer in layers.layers.iter() {

        if min_zoom > layer.maxzoom || max_zoom < layer.minzoom {
            // This layer doesn't apply.
            continue;
        }
        println!("\nLayer {}", layer.id);


        let total_num_of_metatiles: Option<usize> = (min_zoom..max_zoom+1).map(|z| {
            match *bbox {
                None => {
                    let scale = (metatile_scale.trailing_zeros()) as u32;
                    let z = z as u32;
                    if scale >= z {
                        Some(1)
                    } else {
                        Some(2_u64.pow(2u32 * (z-scale)) as usize)
                    }
                },
                Some(ref bbox) => {
                    let this = slippy_map_tiles::size_bbox_zoom_metatiles(&bbox, z, metatile_scale);
                    this
                },
            }})
            .fold(Some(0_usize), |acc, on_this_zoom| {
                match (acc, on_this_zoom) {
                    (Some(x), Some(y)) => x.checked_add(y),
                    _ => None,
                }
            });

        let (printer_tx, printer_rx) = channel();
        let new_bbox: Option<BBox> = bbox.clone();
        let mut printer_thread = thread::spawn(move || { printer::printer(printer_rx, total_num_of_metatiles) });

        let mut metatile_iterator = MetatilesIterator::new_for_bbox_zoom(metatile_scale, &bbox, min_zoom, max_zoom);
        let mut metatile_iterator = Arc::new(Mutex::new(metatile_iterator));

        let mut workers = Vec::with_capacity(num_threads);
        for _ in 0..num_threads {
            // TODO do I need all these clones?
            let my_connection_pool = ConnectionPool::new(layers.get_all_connections());
            let my_printer_tx = printer_tx.clone();
            let my_fileio_tx = fileio_tx.clone();
            let my_metatile_iterator = Arc::clone(&metatile_iterator);
            let my_layer = layer.clone();
            let my_dest = dest.clone();

            // TODO 'should do metatile'

            let handle = thread::spawn(move || {
                worker_one_layer(my_printer_tx, my_fileio_tx, my_metatile_iterator, &my_connection_pool, &my_layer, global_maxzoom);
            });
            workers.push(handle);
        }

        for worker in workers {
            // If one of our worker threads has panic'ed, then this main programme should fail too
            worker.join().ok();
        }

        printer_tx.send(printer::PrinterMessage::Quit).unwrap();
        printer_thread.join().unwrap();
    }


    fileio_tx.send(FileIOMessage::Quit).unwrap();

    println!("All tiles generated. Waiting for all to be written to disk...");

    fileio_thread.join().unwrap();

    println!("Finished.");

}

fn worker_one_layer(printer_tx: Sender<printer::PrinterMessage>, fileio_tx: SyncSender<FileIOMessage>, mut metatile_iterator: Arc<Mutex<MetatilesIterator>>, connection_pool: &ConnectionPool, layer: &Layer, global_maxzoom: u8)
{
    loop {
        let metatile = metatile_iterator.lock().unwrap().next();
        if let None = metatile {
            // The iterator is finished.
            break;
        }
        let metatile = metatile.unwrap();
        let scale = metatile.size() as u32;

        let mut string_store = StringStore::new();
        let tiles = single_layer(&layer, global_maxzoom, &metatile, &connection_pool, &mut string_store);

        let num_tiles = tiles.len();

        for (i, mvt_layer) in tiles.into_iter().enumerate() {
            let i = i as u32;
            let x = i / scale + metatile.x();
            let y = i % scale + metatile.y();
            let tile = slippy_map_tiles::Tile::new(metatile.zoom(), x, y).unwrap();
            
            fileio_tx.send(FileIOMessage::AppendToTile(tile, mvt_layer.to_bytes())).unwrap();
        }

        printer_tx.send(printer::PrinterMessage::DoneTiles(metatile.zoom(), 1, num_tiles)).unwrap();

    }

}

fn tilejson_vector_layers(layers: &Layers, connection_pool: &ConnectionPool) -> serde_json::Value {
    json!({
        "vector_layers": layers.layers.iter().map(|layer| {
            let layer_name = &layer.id;
            let columns = columns_for_layer(layer, connection_pool);
            let minzoom = layer.minzoom;
            let maxzoom = layer.maxzoom;
            let maxzoom = if maxzoom > layers.global_maxzoom { layers.global_maxzoom } else { maxzoom };
            json!({
                "id": layer_name,
                "description": "",
                "minzoom": minzoom,
                "maxzoom": maxzoom,
                "fields": columns.into_iter().collect::<HashMap<_, _>>(),
            })
        }).collect::<Vec<_>>(),
    })
}


fn write_tilejson(layers: &Layers, connection_pool: &ConnectionPool, dest: &PathBuf) {
    let tilejson = json!({
        "tilejson": "2.2.0",
        "tiles": [
            "http://www.example.com/{z}/{x}/{y}.pbf"
        ],
        "minzoom": 0,
        "maxzoom": 14,
        "format": "pbf",
        "vector_layers": tilejson_vector_layers(layers, connection_pool),
    });

    fs::create_dir_all(&dest).unwrap();
    let mut tilejson_file = File::create(dest.join("index.json")).unwrap();
    serde_json::to_writer_pretty(tilejson_file, &tilejson).unwrap();
    
}

fn columns_for_layer(layer: &Layer, connection_pool: &ConnectionPool) -> Vec<(String, String)> {
    let layer_name = &layer.id;

    let conn = connection_pool.connection_for_layer(&layer_name);
    
    let bbox = LocalBBox(0., 0., 0., 0.);
    let res = conn.query(&layer.table.query, &layer.table.params(&bbox, &0., &0., &0.)).unwrap();

    res.columns().iter()
        .filter_map(|column| {
            let name = column.name();
            if name == "way" {
                return None;
            }

            // Sometimes a NULL value can be returned, hence the dance with Option<Value>
            let column_type: Option<String> = match column.type_().name() {
                "float4" => Some("Number".to_string()),
                "float8" => Some("Number".to_string()),
                "text" => Some("String".to_string()),
                "int4" => Some("Number".to_string()),
                "int8" => Some("Number".to_string()),
                "numeric" => Some("Number".to_string()),
                "varchar" => Some("String".to_string()),
                
                // why is there unknown?
                "unknown" => None,
                // Should this be Vec<u8>??
                "bytea" => None,
                x => {
                    eprintln!("Postgres type {:?} not known for layer {:?}", x, layer);
                    unimplemented!()
                },
            };

            if let Some(column_type) = column_type {
                Some((name.to_owned(), column_type))
            } else {
                None
            }
        })
        .collect()
}

pub fn single_metatile(layers: &Layers, metatile: &slippy_map_tiles::Metatile, connection_pool: &ConnectionPool) -> Vec<(slippy_map_tiles::Tile, mapbox_vector_tile::Tile)> {
    let empty_tile = mapbox_vector_tile::Tile::new();
    let scale = metatile.size() as u32;

    let mut results: Vec<mapbox_vector_tile::Tile> = vec![empty_tile; (scale*scale) as usize];

    let mut string_store = stringstore::StringStore::new();

    for layer in layers.layers.iter() {
        let minzoom = layer.minzoom;
        let maxzoom = layer.maxzoom;
        let maxzoom = if maxzoom > layers.global_maxzoom { layers.global_maxzoom } else { maxzoom };
        // Skip layers which are not on this zoom
        if metatile.zoom() < minzoom || metatile.zoom() > maxzoom {
            continue;
        }

        let mvt_layers = single_layer(layer, layers.global_maxzoom, metatile, connection_pool, &mut string_store);
        for (mvt_tile, mvt_layer) in results.iter_mut().zip(mvt_layers.into_iter()) {
            mvt_tile.add_layer(mvt_layer);
        }
        //memory!("Done layer {}", layer.id);

    }


    memory!("Metatile {:?} finished", metatile);
    results.into_iter().enumerate().map(|(i, mvt_tile)| {
        let i = i as u32;
        let x = i / scale + metatile.x();
        let y = i % scale + metatile.y();
        (slippy_map_tiles::Tile::new(metatile.zoom(), x, y).unwrap(), mvt_tile)
    }).collect()
}

fn single_layer(layer: &Layer, global_maxzoom: u8, metatile: &slippy_map_tiles::Metatile, connection_pool: &ConnectionPool, mut string_store: &mut StringStore) -> Vec<mapbox_vector_tile::Layer> {
    let scale = metatile.size() as u32;
    let layer_name = &layer.id;

    let new_layer = mapbox_vector_tile::Layer::new(layer_name.to_string());
    let mut results: Vec<mapbox_vector_tile::Layer> = vec![mapbox_vector_tile::Layer::new(layer_name.to_string()); (scale*scale) as usize];

    // One 'pixel' of buffer space is actually 16 pixels of space now
    let buffer = (layer.buffer as i32) * 16;

    let layer_name = &layer.id;

    let conn = connection_pool.connection_for_layer(&layer_name);
    
    let table = &layer.table;
    // TODO should this be 4096??
    // TODO not confident about this calculation.
    let canvas_size = 256.*(metatile.size() as f64);
    let ll = metatile.sw_corner().to_3857();
    let ur = metatile.ne_corner().to_3857();

    // TODO vtiles have y positive going down, is this correct??
    let tile_width = (ur.0 - ll.0) as f64;
    let tile_height = (ur.1 - ll.1) as f64;

    // calculate how much to expand the bbox to get the buffer.
    // Similar calculation for the pixel size
    let buffer_width = (tile_width / canvas_size)*(buffer as f64);
    let buffer_height = (tile_height / canvas_size)*(buffer as f64);

    let minx = ll.0 as f64;
    let miny = ll.1 as f64;
    let maxx = ur.0 as f64;
    let maxy = ur.1 as f64;

    let bbox = LocalBBox(minx-buffer_width, miny-buffer_height, maxx+buffer_width, maxy+buffer_height);
    assert!(tile_height > 0.);
    assert!(tile_width > 0.);

    let pixel_width = (tile_width / canvas_size) as f32;
    let pixel_height = (tile_height / canvas_size) as f32;

    let scale_denom = scale_denominator_for_zoom(metatile.zoom());
    let stmt = conn.prepare_cached(&layer.table.query).unwrap();
    let res = stmt.query(&table.params(&bbox, &pixel_width, &pixel_height, &scale_denom)).unwrap();

    if res.is_empty() {
        return results;
    }


    let extent = (new_layer.extent as f64)*(metatile.size() as f64);

    let columns: Vec<_> = res.columns().iter().skip(1).filter(|c| c.name() != "way").collect();

    let mut res = res.iter().enumerate();

    let mut num_objects = 0;

    for (i, row) in res {
        num_objects += 1;
        let bad_obj = false && metatile.zoom() == 3 && i == 4_579;
        if i % 5_000 == 0 {
            //memory!("layer {} have done {} objects", layer_name, i.separated_string());
        }

        // First object is the ST_AsBinary
        // TODO Does this do any copies that we don't want?
        let wkb_bytes: Vec<u8> = row.get(0);

        //println!("\nL {} bytes {:?}", line!(), wkb_bytes);

        let geom: geo::Geometry<f64> = match wkb::wkb_to_geom(&mut wkb_bytes.as_slice()) {
            Err(e) => {
                // TODO investigate this more
                //eprintln!("{}:{} Metatile: {:?} WKB reading error {:?}, layer {} row {:?}", file!(), line!(), metatile, e, layer_name, row);
                continue;
            },
            Ok(g) => g,
        };
        drop(wkb_bytes);

        // TODO not sure about this
        let pixel_size: f64 = tile_width/extent;

        // TODO there are a lot of calls to `is_valid`, which is computationally expensive, but
        // removing them makes it slower, probably because of the lots of invalid geoms

        //if bad_obj {
        //    println!("\nL {} geom {:100}", line!(), format!("{:?}", geom));
        //    println!("\nL {} minx {} maxx {} miny {} maxy {} extent {}", line!(), minx, maxx, miny, maxy, extent);
        //}
        let mut geom = match remap_geometry(geom, minx, maxx, miny, maxy, extent) {
            None => { continue; }
            Some(g) => g,
        };

        let geom = match simplify::remove_unneeded_points(geom) {
            None => { continue; },
            Some(g) => g,
        };
        //if bad_obj { println!("{}:{} geom {:100}", file!(), line!(), format!("{:?}", geom)); }

        //let geom = validity::make_valid(geom);

        //debug_assert!(is_valid(&geom), "L {} Geometry is invalid after remap: {:100}", line!(), format!("{:?}", geom));
        //validity::ensure_polygon_orientation(&mut geom);
        //if ! is_valid_skip_expensive(&geom) {
        //    continue;
        //}

        // Only do the simplification if we're not at maxzoom. We've already removed extra
        //
        // points in remove_unneeded_points above
        //println!("{} L {}", file!(), line!());
        let geom = if metatile.zoom() < global_maxzoom {
                match simplify::simplify(geom, 8) {
                    None => {
                        continue;
                    },
                    Some(g) => g,
                }
        } else { geom };
        //println!("{} L {}", file!(), line!());
        //debug_assert!(is_valid(&geom), "L {} Geometry is invalid after remap: {:100}", line!(), format!("{:?}", geom));
        //if bad_obj { println!("{}:{} geom {:102}", file!(), line!(), format!("{:?}", geom)); }
        //if bad_obj {
        //    let mut g2 = geom.clone();
        //    validity::ensure_polygon_orientation(&mut g2);
        //    print_geom_as_geojson(&g2, 4096.*8.);
        //}
        
        // After simplifying a geometry, it's possible it becomes invalid. So we just skip the
        // geometries in that case.
        //if ! is_valid(&geom) {
        //    continue;
        //}

        // clip geometry, so no part of it goes outside the bbox. PostgreSQL will return
        // anything that overlaps.
        let geom = match clip_to_bbox(Cow::Owned(geom), &geo::Bbox{ xmin: -(buffer as i32), xmax: extent as i32 + buffer as i32, ymin: -(buffer as i32), ymax: extent as i32 + buffer as i32 }) {
            None => {
                // geometry is outside the bbox, so skip
                continue;
            },
            Some(g) => g,
        };

        //let geom = validity::make_valid(geom);
        //debug_assert!(is_valid(&geom), "L {} Geometry is invalid after clip_to_bbox: {:100}", line!(), format!("{:?}", geom));
                
        let mut properties = mapbox_vector_tile::Properties::new();

        for column in columns.iter() {
            let name = column.name();

            // Sometimes a NULL value can be returned, hence the dance with Option<Value>
            let value: Option<mapbox_vector_tile::Value> = match column.type_().name() {
                "float4" => row.get_opt(name).map(|x| x.ok().map(mapbox_vector_tile::Value::Float)).unwrap_or(None),
                "float8" => row.get_opt(name).map(|x| x.ok().map(mapbox_vector_tile::Value::Double)).unwrap_or(None),

                "text" =>  row.get_opt(name).map(|x| x.ok().map(|s: String| mapbox_vector_tile::Value::String(string_store.get_string(s)))).unwrap_or(None),
                "varchar" =>  row.get_opt(name).map(|x| x.ok().map(|s: String| mapbox_vector_tile::Value::String(string_store.get_string(s)))).unwrap_or(None),

                "int4" => row.get_opt(name).map(|x| x.ok().map(|y| { let val: i32 = y; mapbox_vector_tile::Value::Int(val as i64) })).unwrap_or(None),
                "int8" => row.get_opt(name).map(|x| x.ok().map(|y| { let val: i64 = y; mapbox_vector_tile::Value::Int(val as i64) })).unwrap_or(None),
                
                // TODO not 100% sure numeric is correct here
                "numeric" => row.get_opt(name).map(|x| x.ok().map(mapbox_vector_tile::Value::Double)).unwrap_or(None),
                
                // why is there unknown?
                "unknown" => None,
                x => {
                    eprintln!("Postgres type {:?} not known", x);
                    unimplemented!()
                },
            };


            if let Some(v) = value {
                properties.insert(string_store.get_string(name.to_string()), v);
            }
        }
        properties.0.shrink_to_fit();

        let mut geoms: Vec<_> = clip_geometry_to_tiles(&metatile, geom, buffer).into_iter().filter_map(
            |(t, g)| match g {
                None => None,
                // TODO probably could use .map/.and_then here
                Some(mut g) => {
                    //debug_assert!(is_valid(&g), "L {} Geometry is invalid after clip_geometry_to_tiles: {:?}", line!(), g);

                    match validity::make_valid(g) {
                        None => None,
                        Some(mut g) => {
                            if is_valid(&g) {
                                validity::ensure_polygon_orientation(&mut g);
                                Some((t, g))
                            } else {
                                warn!("make_valid returned an invalid geometry: {:?}", g);
                                None
                            }
                        },
                    }
                },
            }).collect();

        // If there are >1 tiles, then we don't want to clone the properties everytime. So share
        // the data between all mapbox_vector_tile::Features using a Rc.
        // This is only a small speed up.
        let properties = Rc::new(properties);

        for (tile, mut geom) in geoms.into_iter() {

            let i = (tile.x() - metatile.x()) as i32;
            let j = (tile.y() - metatile.y()) as i32;

            let xoff = i*4096;
            let yoff = j*4096;

            geom.map_coords_inplace(&|&(x, y)| ( (x - xoff), (y - yoff)));

            let feature = mapbox_vector_tile::Feature::new(geom, properties.clone());
            let n = (i*(scale as i32) + j) as usize;
            results.get_mut(n).unwrap().add_feature(feature);

        };

    }
    memory!("Finished layer {}, there were {} object", layer_name, num_objects.separated_string());

    results

}


fn remap_linestring(ls: LineString<f64>, minx: f64, maxx: f64, miny: f64, maxy: f64, size: f64, should_be_ring: bool) -> Option<LineString<i32>> {
    
    let remap_xy = |x: f64, y: f64| -> (i32, i32) {
        let x: f64 = ((x - minx) / (maxx - minx))*size;
        let x = x.round();
        debug_assert!(x <= i32::max_value() as f64);
        debug_assert!(x >= i32::min_value() as f64);
        let x: i32 = x as i32;

        // y axies goes down, hence different ordering for y
        let y: f64 = ((maxy - y) / (maxy - miny))*size;
        let y = y.round();
        debug_assert!(y <= i32::max_value() as f64);
        debug_assert!(y >= i32::min_value() as f64);
        let y: i32 = y as i32;


        (x, y)
    };

    let mut new_points: Vec<Point<_>> = Vec::with_capacity(ls.0.len());
    let last_xy = remap_xy(ls.0[0].x(), ls.0[0].y());
    let mut last_x = last_xy.0;
    let mut last_y = last_xy.1;
    new_points.push(Point::new(last_x, last_y));

    // Remap all points, but don't add a point if it's the same location as the last point.
    for p in ls.0.into_iter().skip(1) {
        let new_xy = remap_xy(p.x(), p.y());
        if ! ( new_xy.0 == last_x && new_xy.1 == last_y ) {
            last_x = new_xy.0;
            last_y = new_xy.1;
            new_points.push(Point::new(last_x, last_y));
        }
    }

    if should_be_ring {
        if new_points.len() >= 4 && new_points[0] == new_points[new_points.len()-1] { 
            Some(LineString(new_points))
        } else {
            None
        }
    } else {
        if new_points.len() >= 2 {
            Some(LineString(new_points))
        } else {
            None
        }
    }
}

fn remap_geometry(geom: Geometry<f64>, minx: f64, maxx: f64, miny: f64, maxy: f64, size: f64) -> Option<Geometry<i32>> {

    fn conv(x: f64) -> i32 {
        let x = x.round();
        debug_assert!(x <= i32::max_value() as f64);
        debug_assert!(x >= i32::min_value() as f64);

        x as i32
    }
    
    let remap_xy = |x: f64, y: f64| -> (i32, i32) {
        (
            conv(((x - minx) / (maxx - minx))*size),
            conv(((maxy - y) / (maxy - miny))*size)
        )
    };

    match geom {
        Geometry::Point(p) => {
            let xy = remap_xy(p.x(), p.y());
            Some(Geometry::Point(Point::new(xy.0, xy.1)))
        },
        Geometry::MultiPoint(mp) => {
            if mp.0.is_empty() {
                None
            } else {
                Some(Geometry::MultiPoint(MultiPoint(mp.0.into_iter().map(|p| {
                        let xy = remap_xy(p.x(), p.y());
                        Point::new(xy.0, xy.1)
                    }
                    ).collect::<Vec<Point<_>>>())))
            }
        },
        Geometry::LineString(ls) => {
            remap_linestring(ls, minx, maxx, miny, maxy, size, false).and_then(|ls| Some(Geometry::LineString(ls)))
        },
        Geometry::MultiLineString(mls) => {
            let mut res: Vec<LineString<_>> = mls.0.into_iter().filter_map(|ls| remap_linestring(ls, minx, maxx, miny, maxy, size, false)).collect();
            match res.len() {
                0 => None,
                1 => Some(Geometry::LineString(res.remove(0))),
                _ => Some(Geometry::MultiLineString(MultiLineString(res))),
            }
        },
        Geometry::Polygon(p) => {
            let Polygon{ exterior, interiors } = p;
            match remap_linestring(exterior, minx, maxx, miny, maxy, size, true) {
                // Exterior gets simplified away
                None => None,
                Some(exterior) => {
                    let interiors: Vec<LineString<_>> = interiors.into_iter().filter_map(|int| remap_linestring(int, minx, maxx, miny, maxy, size, true)).collect();
                    Some(Geometry::Polygon(Polygon::new(exterior, interiors)))
                }
            }
        }
        Geometry::MultiPolygon(mp) => {
            let mut res: Vec<Polygon<_>> = mp.0.into_iter().filter_map(|p| {
                let Polygon{ exterior, interiors } = p;
                match remap_linestring(exterior, minx, maxx, miny, maxy, size, true) {
                        // Exterior gets simplified away
                    None => None,
                    Some(exterior) => {
                        let interiors: Vec<LineString<_>> = interiors.into_iter().filter_map(|int| remap_linestring(int, minx, maxx, miny, maxy, size, true)).collect();
                        Some(Polygon::new(exterior, interiors))
                    }
                }
            }).collect();

            match res.len() {
                0 => None,
                1 => Some(Geometry::Polygon(res.remove(0))),
                _ => Some(Geometry::MultiPolygon(MultiPolygon(res))),
            }

        }

        _ => unimplemented!()
    }
}

fn x_to_lon<T: CoordinateType+Into<f64>>(x: T, extent: f64) -> f64 {
    let earth_radius = 6378137.;
    let x: f64 = x.into();
    let x = (x/extent) * (2.*20037508.34) - 20037508.34;
    //let x = self.lon() * 20037508.34 / 180.;

    (x/earth_radius).to_degrees()
}

fn y_to_lat<T: CoordinateType+Into<f64>>(y: T, extent: f64) -> f64 {
    let old_y = y;
    let y: f64 = y.into();
    let y = y/extent;
    
    let pi = std::f64::consts::PI;

    ((1. - 2.*y) * pi).sinh().atan().to_degrees()
}

fn print_geom_as_geojson<T: CoordinateType+Into<f64>>(geom: &Geometry<T>, extent: f64) {
    println!("{}", geom_as_geojson(geom, extent));
}

fn geom_as_geojson<T: CoordinateType+Into<f64>>(geom: &Geometry<T>, extent: f64) -> String {
    let mut output = String::new();

    let geojson = |ls: &LineString<T>| -> String {
        format!("[{}]", ls.0.iter().map(|p| format!("[{}, {}]", x_to_lon(p.x(), extent), y_to_lat(p.y(), extent))).collect::<Vec<_>>().join(", "))
    };

    let points_numbered = |p: &Vec<Point<T>>| -> String {
        p.iter().enumerate().filter_map(|(i, p)| {
            let col = match i / 100 {
                0 => "#f00",
                1 => "#0f0",
                2 => "#00f",
                3 => "#ff0",
                4 => "#f0f",
                5 => "#fff",
                6 => "#000",
                7 => "#800",
                8 => "#080",
                9 => "#008",
                10 => "#880",
                11 => "#808",
                12 => "#400",
                13 => "#040",
                14 => "#004",
                15 => "#440",
                _ => { return None; }
            };
            let num = i % 100;
            Some(format!("{{\"type\":\"Feature\", \"properties\": {{\"marker-symbol\": \"{num}\", \"marker-color\": \"{col}\", \"title\":\"{i}\"}}, \"geometry\": {{ \"type\":\"Point\", \"coordinates\":[{x:?}, {y:?}]}}}}", x=x_to_lon(p.x(), extent), y=y_to_lat(p.y(), extent), num=num, col=col, i=i))
        }).collect::<Vec<_>>().join(", ")
    };

    write!(output, "\n{{\"type\":\"FeatureCollection\", \"features\":[").unwrap();
    match *geom {
        Geometry::Polygon(ref poly) => {
            write!(output, "{{\"type\": \"Feature\", \"properties\":{{}}, \"geometry\":{{ \"type\":\"Polygon\", \"coordinates\": [").unwrap();
            write!(output, "{}", geojson(&poly.exterior)).unwrap();
            if !poly.interiors.is_empty() {
                write!(output, ", ").unwrap();
                write!(output, "{}", poly.interiors.iter().map(|int| geojson(&int)).collect::<Vec<_>>().join(", ")).unwrap();
            }
            write!(output, "]}}}}").unwrap();
            write!(output, ", {}", points_numbered(&poly.exterior.0)).unwrap();
        },
        Geometry::MultiPolygon(ref mp) => {
            write!(output, "{{\"type\": \"Feature\", \"properties\":{{}}, \"geometry\":{{ \"type\":\"MultiPolygon\", \"coordinates\": [").unwrap();
            let mut first = true;
            for poly in mp.0.iter() {
                if !first { write!(output, ", ").unwrap(); }
                write!(output, "[").unwrap();
                write!(output, "{}", geojson(&poly.exterior)).unwrap();
                if !poly.interiors.is_empty() {
                    write!(output, ", ").unwrap();
                    write!(output, "{}", poly.interiors.iter().map(|int| geojson(&int)).collect::<Vec<_>>().join(", ")).unwrap();
                }
                write!(output, "]").unwrap();
                first = false;
            }
            write!(output, "]}}}}").unwrap();
        },
        Geometry::LineString(ref ls) => {
            write!(output, "{{\"type\": \"Feature\", \"properties\":{{}}, \"geometry\":{{ \"type\":\"LineString\", \"coordinates\": ").unwrap();
            write!(output, "{}", geojson(&ls)).unwrap();
            write!(output, "}}}}").unwrap();
            write!(output, ", {}", points_numbered(&ls.0)).unwrap();
        },
        Geometry::MultiLineString(ref mls) => {
            write!(output, "{{\"type\": \"Feature\", \"properties\":{{}}, \"geometry\":{{ \"type\":\"MultiLineString\", \"coordinates\": [").unwrap();
            write!(output, "{}", mls.0.iter().map(geojson).collect::<Vec<_>>().join(", ")).unwrap();
            write!(output, "]}}}}").unwrap();
        },

        _ => unimplemented!(),
    }
    write!(output, "]}}\n").unwrap();

    output
}

#[derive(Debug)]
struct LocalBBox<T: num_traits::Float+Into<f64>>(T, T, T, T);
//let bbox = format!("ST_SetSRID(ST_MakeBox2D(ST_Point({llx}, {lly}), ST_Point({urx}, {ury})), 3857)", llx=(minx-buffer_width), lly=(miny-buffer_height), urx=(maxx+buffer_width), ury=(maxy+buffer_height));

impl<T: num_traits::Float+Into<f64>+std::fmt::Debug> ToSql for LocalBBox<T> {
    fn accepts(ty: &Type) -> bool {
        ty.name() == "geometry"
    }

    fn to_sql(&self, ty: &Type, mut out: &mut Vec<u8>) -> Result<IsNull, Box<::std::error::Error+Sync+Send>> {
        let minx = self.0.into();
        let miny = self.1.into();
        let maxx = self.2.into();
        let maxy = self.3.into();
        // a--b
        // |  |
        // d--c
        let a = Point::new(minx, miny);
        let b = Point::new(maxx, miny);
        let c = Point::new(maxx, maxy);
        let d = Point::new(minx, maxy);
        let polygon: Geometry<f64> = Polygon::new(vec![a, b, c, d, a].into(), vec![]).into();
        wkb::write_geom_to_wkb(&polygon, &mut out);
        Ok(IsNull::No)
    }

    to_sql_checked!();
}
