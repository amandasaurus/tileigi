#![allow(dead_code,unused_imports,unused_variables,unused_mut,unused_assignments)]

extern crate postgres;
extern crate slippy_map_tiles;
extern crate yaml_rust;
extern crate mapbox_vector_tile;
extern crate wkb;
extern crate geo;
#[macro_use]
extern crate serde_json;

extern crate users;
extern crate num_traits;
extern crate rusqlite;
extern crate md5;

use std::fs::File;
use std::fs;
use std::io::prelude::*;
use std::path::PathBuf;
use std::collections::{HashSet, HashMap};
use std::time::Instant;
use std::borrow::{Cow, Borrow};

use std::thread;
use std::sync::mpsc::{channel, Sender};
use std::sync::{Arc, Mutex};

use yaml_rust::{YamlLoader, Yaml};

use postgres::{Connection, TlsMode};
use postgres::params::ConnectParams;

use slippy_map_tiles::{BBox, Metatile, MetatilesIterator};

use geo::*;
use geo::algorithm::map_coords::MapCoords;
use geo::algorithm::map_coords::MapCoordsInplace;
use geo::algorithm::boundingbox::BoundingBox;
use geo::algorithm::contains::Contains;

mod clip;
use clip::{clip_to_bbox,clip_geometry_to_tiles};

mod validity;
use validity::{is_valid, is_valid_skip_expensive};

mod printer;
mod fileio;
mod simplify;
mod fraction;

use fileio::FileIOMessage;

#[cfg(test)]
mod test;

pub enum TileDestinationType {
    TileStashDirectory(String),
    MBTiles(String),
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
    table: String,
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
            .map(|layer| Layer{
                id: layer["id"].as_str().unwrap().to_owned(),
                dbname: layer["Datasource"]["dbname"].as_str().map(|x| x.to_owned()),
                minzoom: layer["properties"]["minzoom"].as_i64().map(|x| x as u8).unwrap_or(global_minzoom) as u8,
                maxzoom: layer["properties"]["maxzoom"].as_i64().map(|x| x as u8).unwrap_or(global_maxzoom) as u8,
                buffer: layer["properties"]["buffer-size"].as_i64().map(|x| x as u16).unwrap_or(0) as u16,
                table: layer["Datasource"]["table"].as_str().unwrap().to_owned(),
            })
            // TODO finish this thing
            //.map(|layer| {
            //    layer["Datasource"]["table"].as_str().unwrap()
            //        .replace("!bbox!", "$1")
            //        .replace("!pixel_width!", "$2").replace("!pixel_height!", "$3")
            //        .replace("!scale_denominator!", "$4");
            //    layer
            //})
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

fn scale_denominator_for_zoom(zoom: u8) -> &'static str {
    match zoom {
        0 => "250000000000",
        1 => "500000000",
        2 => "200000000",
        3 => "100000000",
        4 => "50000000",
        5 => "25000000",
        6 => "12500000",
        7 => "6500000",
        8 => "3000000",
        9 => "1500000",
        10 => "750000",
        11 => "400000",
        12 => "200000",
        13 => "100000",
        14 => "50000",
        15 => "25000",
        16 => "12500",
        17 => "5000",
        18 => "2500",
        _ => {
            eprintln!("Unsupported zoom ({})", zoom);
            unimplemented!();
        }
    }
}

pub fn generate_all(filename: &str, min_zoom: u8, max_zoom: u8, bbox: &Option<BBox>, dest: &TileDestinationType, if_not_exists: bool, compress: bool, metatile_scale: u8, num_threads: usize) {
    let layers = Layers::from_file(filename);

    let connection_pool = ConnectionPool::new(layers.get_all_connections());


    let (printer_tx, printer_rx) = channel();
    let mut printer_thread = thread::spawn(move || { printer::printer(printer_rx) });

    let (fileio_tx, fileio_rx) = channel();

    let mut fileio_thread = match dest {
        &TileDestinationType::TileStashDirectory(ref path) => {
            let path = PathBuf::from(path);
            fs::create_dir_all(&path).unwrap();
            write_tilejson(&layers, &connection_pool, &path);
            let tile_dest = fileio::TileStashDirectory::new(&path);
            thread::spawn(move || { fileio::fileio_thread(fileio_rx, Box::new(tile_dest)) })
        },
        &TileDestinationType::MBTiles(ref path) => {
            let path = PathBuf::from(path);
            //fs::create_dir_all(&path).unwrap();
            //write_tilejson(&layers, &connection_pool, &path);
            let tile_dest = fileio::MBTiles::new(&path);
            thread::spawn(move || { fileio::fileio_thread(fileio_rx, Box::new(tile_dest)) })
        },
    };

    let mut metatile_iterator = MetatilesIterator::new_for_bbox_zoom(metatile_scale, &bbox, min_zoom, max_zoom);
    let mut metatile_iterator = Arc::new(Mutex::new(metatile_iterator));

    let mut workers = Vec::with_capacity(num_threads);
    for _ in 0..num_threads {
        // TODO do I need all these clones?
        let my_connection_pool = ConnectionPool::new(layers.get_all_connections());
        let my_printer_tx = printer_tx.clone();
        let my_fileio_tx = fileio_tx.clone();
        let my_metatile_iterator = Arc::clone(&metatile_iterator);
        let my_layers = layers.clone();
        let handle = thread::spawn(move || {
            worker(my_printer_tx, my_fileio_tx, my_metatile_iterator, &my_connection_pool, &my_layers);
        });
        workers.push(handle);
    }

    for worker in workers {
        // If one of our worker threads has panic'ed, then this main programme should fail too
        worker.join().ok();
    }

    fileio_tx.send(FileIOMessage::Quit).unwrap();
    fileio_thread.join().unwrap();

    printer_tx.send(printer::PrinterMessage::Quit).unwrap();
    printer_thread.join().unwrap();
}

fn worker(printer_tx: Sender<printer::PrinterMessage>, fileio_tx: Sender<FileIOMessage>, mut metatile_iterator: Arc<Mutex<MetatilesIterator>>, connection_pool: &ConnectionPool, layers: &Layers) {
    loop {
        let metatile = metatile_iterator.lock().unwrap().next();
        if let None = metatile {
            // The iterator is finished.
            break;
        }
        let metatile = metatile.unwrap();

        let tiles = single_metatile(&layers, &metatile, &connection_pool);
        let num_tiles = tiles.len();

        // TODO the tile exists check needs to work with metatiles

        for (tile, pbf) in tiles.into_iter() {
            let bytes = pbf.to_compressed_bytes();
            fileio_tx.send(FileIOMessage::SaveTile(tile, bytes)).unwrap();
        }

        printer_tx.send(printer::PrinterMessage::DoneTiles(metatile.zoom(), 1, metatile.scale() as usize)).unwrap();

    }

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
    });

    let mut tilejson_file = File::create(dest.join("index.json")).unwrap();
    serde_json::to_writer_pretty(tilejson_file, &tilejson).unwrap();
    
}

fn columns_for_layer(layer: &Layer, connection_pool: &ConnectionPool) -> Vec<(String, String)> {
    let layer_name = &layer.id;

    let conn = connection_pool.connection_for_layer(&layer_name);
    
    let table = &layer.table;
    let table = table
        .replace("!pixel_width!", "0").replace("!pixel_height!","0")
        .replace("!bbox!", "ST_Point(0, 0)")
        .replace("!scale_denominator!", "0");

    let query = format!("SELECT * from {table} LIMIT 0", table=table);
    let res = conn.query(&query, &[]).unwrap();

    res.columns().iter().filter(|c| c.name() != "way")
        .filter_map(|column| {
            let name = column.name().to_owned();

            // Sometimes a NULL value can be returned, hence the dance with Option<Value>
            let column_type: Option<String> = match column.type_().name() {
                "float4" => Some("Number".to_string()),
                "float8" => Some("Number".to_string()),
                "text" => Some("String".to_string()),
                "int4" => Some("Number".to_string()),
                "int8" => Some("Number".to_string()),
                "numeric" => Some("Number".to_string()),
                
                // why is there unknown?
                "unknown" => None,
                x => {
                    eprintln!("Postgres type {:?} not known", x);
                    unimplemented!()
                },
            };

            if let Some(column_type) = column_type {
                Some((name, column_type))
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

    for layer in layers.layers.iter() {
        let minzoom = layer.minzoom;
        let maxzoom = layer.maxzoom;
        let maxzoom = if maxzoom > layers.global_maxzoom { layers.global_maxzoom } else { maxzoom };

        // Skip layers which are not on this zoom
        if metatile.zoom() < minzoom || metatile.zoom() > maxzoom {
            continue;
        }

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

        let bbox = format!("ST_SetSRID(ST_MakeBox2D(ST_Point({llx}, {lly}), ST_Point({urx}, {ury})), 3857)", llx=(minx-buffer_width), lly=(miny-buffer_height), urx=(maxx+buffer_width), ury=(maxy+buffer_height));
        assert!(tile_height > 0.);
        assert!(tile_width > 0.);

        let pixel_width = format!("{}", tile_width / canvas_size);
        let pixel_height = format!("{}", tile_height / canvas_size);

        // Would it be faster to have a prepared statement which we then execute many times,
        // with these !params! being $1 etc?
        // TODO use prepared statements
        let table = table.replace("!pixel_width!", &pixel_width).replace("!pixel_height!", &pixel_height).replace("!bbox!", &bbox).replace("!scale_denominator!", scale_denominator_for_zoom(metatile.zoom()));

        let query = format!("SELECT ST_AsBinary(way), * from {table} where way && {bbox}", table=table, bbox=bbox);
        let res = conn.query(&query, &[]).unwrap();

        if res.is_empty() {
            continue;
        }

        // Ensure all tiles in this metatile have this layer
        let new_layer = mapbox_vector_tile::Layer::new(layer_name.to_string());
        for mvt in results.iter_mut() {
            mvt.add_layer(new_layer.clone());
        }

        let extent = (new_layer.extent as f64)*(metatile.size() as f64);

        let columns: Vec<_> = res.columns().iter().skip(1).filter(|c| c.name() != "way").collect();

        let mut res = res.iter().enumerate();

        let mut num_objects = 0;

        for (i, row) in res {
            num_objects += 1;


            // First object is the ST_AsBinary
            // TODO Does this do any copies that we don't want?
            let wkb_bytes: Vec<u8> = row.get(0);

            //println!("\nL {} bytes {:?}", line!(), wkb_bytes);

            let geom: geo::Geometry<f64> = match wkb::wkb_to_geom(&mut wkb_bytes.as_slice()) {
                Err(e) => {
                    // TODO investigate this more
                    eprintln!("Metatile: {:?} WKB reading error {:?}, first few bytes geom: {:?}", metatile, e, wkb_bytes.into_iter().take(20).collect::<Vec<u8>>());
                    continue;
                },
                Ok(g) => g,
            };

            let mut bad_obj = false;

            if bad_obj {
                println!("\nL {} starting tile {:?}", line!(), metatile);
            }

            // TODO not sure about this
            let pixel_size: f64 = tile_width/extent;

            // TODO there are a lot of calls to `is_valid`, which is computationally expensive, but
            // removing them makes it slower, probably because of the lots of invalid geoms

            if bad_obj {
                println!("\nL {} geom {:?}", line!(), geom);
                println!("\nL {} minx {} maxx {} miny {} maxy {} extent {}", line!(), minx, maxx, miny, maxy, extent);
            }
            let mut geom = remap_geometry(geom, minx, maxx, miny, maxy, extent);
            if geom.is_none() {
                if bad_obj {
                    println!("\nL {} none after remap", line!());
                }
                continue;
            }
            if bad_obj {
                println!("\nL {} geom {:?}", line!(), geom);
            }
            let mut geom = geom.unwrap();

            simplify::remove_unneeded_points(&mut geom);

            //let geom = validity::make_valid(geom);

            //debug_assert!(is_valid(&geom), "L {} Geometry is invalid after remap: {:?}", line!(), geom);
            if bad_obj {
                println!("\nL {} geom {:?}", line!(), geom);
            }
            //validity::ensure_polygon_orientation(&mut geom);
            //if ! is_valid_skip_expensive(&geom) {
            //    continue;
            //}

            // Only do the simplification if we're not at maxzoom. We've already removed extra
            //
            // points in remove_unneeded_points above
            let geom = if metatile.zoom() < layers.global_maxzoom {
                    match simplify::simplify(geom, 8) {
                        None => { continue; },
                        Some(g) => g,
                    }
            } else { geom };
            //debug_assert!(is_valid(&geom), "L {} Geometry is invalid after remap: {:?}", line!(), geom);
            
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

            let geom = validity::make_valid(geom);
            debug_assert!(is_valid(&geom), "L {} Geometry is invalid after clip_to_bbox: {:?}", line!(), geom);
            if bad_obj {
                println!("\nL {} geom {:?}", line!(), geom);
            }
                    
            let mut properties = mapbox_vector_tile::Properties::new();

            for column in columns.iter() {
                let name = column.name();

                // Sometimes a NULL value can be returned, hence the dance with Option<Value>
                let value: Option<mapbox_vector_tile::Value> = match column.type_().name() {
                    "float4" => row.get_opt(name).map(|x| x.ok().map(mapbox_vector_tile::Value::Float)).unwrap_or(None),
                    "float8" => row.get_opt(name).map(|x| x.ok().map(mapbox_vector_tile::Value::Double)).unwrap_or(None),
                    "text" => row.get_opt(name).map(|x| x.ok().map(mapbox_vector_tile::Value::String)).unwrap_or(None),
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
                    properties.insert(name, v);
                }
            }

            if bad_obj {
                println!("\nL {} {} properties {:?}", line!(), layer_name, properties);
            }


            if bad_obj {
                println!("\nL {} geom {:?}", line!(), geom);
            }
            let mut geoms: Vec<_> = clip_geometry_to_tiles(&metatile, geom, buffer).into_iter().filter_map(
                |(t, g)| match g {
                    Some(mut g) => {
                        //debug_assert!(is_valid(&g), "L {} Geometry is invalid after clip_geometry_to_tiles: {:?}", line!(), g);
                        if is_valid(&g) {
                            validity::ensure_polygon_orientation(&mut g);
                            Some((t, g))
                        } else {
                            None
                        }
                    },
                    None => None,
                }).collect();
            geoms.reverse();


            let mut save_single_tile = |tile: slippy_map_tiles::Tile, mut geom: Geometry<i32>| {

                let i = (tile.x() - metatile.x()) as i32;
                let j = (tile.y() - metatile.y()) as i32;

                geom.map_coords_inplace(&|&(x, y)| ( (x - (4096*i)), (y - (4096*j))));

                debug_assert!(is_valid(&geom));

                if bad_obj {
                    println!("\nL {} geom {:?}", line!(), geom);
                }

                let feature = mapbox_vector_tile::Feature::new(geom, properties.clone());
                let i = ((tile.x() - metatile.x())*scale + (tile.y() - metatile.y())) as usize;
                let mvt_tile = results.get_mut(i).unwrap();
                mvt_tile.add_feature(&layer_name, feature);

            };

            // In cases where there is only one geometry here, we don't want to clone the
            // `properties`, and instead move it. If there are N geometries, we want to do N-1
            // clones, and 1 move. Hence the duplication with the loop.
            //
            // We want to do this in order, so we reverse the vec, and the pop from the end
            // (which is the original front).
            // TODO rather than filtering out invalid geoms here, prevent the clipping code from
            // generating invalid geoms in the first place
            // One error was creating a linestring with 2 points, both the same
            loop {
                if geoms.len() <= 1 { break; }
                if let Some((tile, geom)) = geoms.pop() {
                    save_single_tile(tile, geom);
                }
            }

            if geoms.is_empty() { continue; }
            if let Some((tile, geom)) = geoms.pop() {
                save_single_tile(tile, geom);
            }

        }

        
    }

    results.into_iter().enumerate().map(|(i, mvt_tile)| {
        let i = i as u32;
        let x = i / scale + metatile.x();
        let y = i % scale + metatile.y();
        (slippy_map_tiles::Tile::new(metatile.zoom(), x, y).unwrap(), mvt_tile)
    }).collect()

}


fn remap_linestring(ls: LineString<f64>, minx: f64, maxx: f64, miny: f64, maxy: f64, size: f64, should_be_ring: bool) -> Option<LineString<i32>> {
    fn conv(x: f64) -> i32 {
        let x = x.round();
        debug_assert!(x <= i32::max_value() as f64);
        debug_assert!(x >= i32::min_value() as f64);

        x as i32
    }
    
    let remap_xy = |x: f64, y: f64| -> (i32, i32) {
        (
            conv(((x - minx) / (maxx - minx))*size),

            // y axies goes down, hence different ordering for y
            conv(((maxy - y) / (maxy - miny))*size)
        )
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

fn x_to_lon(x: i32, extent: f64) -> f64 {
    let earth_radius = 6378137.;
    let x = x as f64;
    let x = (x/extent) * (2.*20037508.34) - 20037508.34;
    //let x = self.lon() * 20037508.34 / 180.;

    (x/earth_radius).to_degrees()
}

fn y_to_lat(y: i32, extent: f64) -> f64 {
    let old_y = y;
    let y = y as f64;
    let y = y/extent;
    
    let pi = std::f64::consts::PI;

    ((1. - 2.*y) * pi).sinh().atan().to_degrees()
}

fn print_geom_as_geojson(geom: &Geometry<i32>, extent: f64) {
    let geojson = |ls: &LineString<i32>| -> String {
        format!("[{}]", ls.0.iter().map(|p| format!("[{}, {}]", x_to_lon(p.x(), extent), y_to_lat(p.y(), extent))).collect::<Vec<_>>().join(", "))
    };

    println!("\n");
    match *geom {
        Geometry::Polygon(ref poly) => {
            print!("{{\"type\": \"Polygon\", \"coordinates\": [");
            print!("{}", geojson(&poly.exterior));
            if !poly.interiors.is_empty() {
                print!(", ");
                print!("{}", poly.interiors.iter().map(|int| geojson(&int)).collect::<Vec<_>>().join(", "));
            }
            print!("]}}");
        },
        Geometry::MultiPolygon(ref mp) => {
            print!("{{\"type\": \"MultiPolygon\", \"coordinates\": [");
            let mut first = true;
            for poly in mp.0.iter() {
                if !first { print!(", "); }
                print!("[");
                print!("{}", geojson(&poly.exterior));
                if !poly.interiors.is_empty() {
                    print!(", ");
                    print!("{}", poly.interiors.iter().map(|int| geojson(&int)).collect::<Vec<_>>().join(", "));
                }
                print!("]");
                first = false;
            }


            print!("]}}");
        },
        _ => unimplemented!(),
    }
    println!("\n");
}
