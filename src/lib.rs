#![allow(dead_code,unused_imports,unused_variables,unused_mut)]

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

use std::fs::File;
use std::fs;
use std::io::prelude::*;
use std::path::Path;
use std::collections::{HashSet, HashMap};
use std::time::Instant;
use std::borrow::{Cow, Borrow};

use yaml_rust::{YamlLoader, Yaml};

use postgres::{Connection, TlsMode};
use postgres::params::ConnectParams;

use slippy_map_tiles::{BBox, Metatile, MetatilesIterator};

use geo::*;
use geo::algorithm::simplify::Simplify;
use geo::algorithm::map_coords::MapCoords;
use geo::algorithm::map_coords::MapCoordsInplace;
use geo::algorithm::boundingbox::BoundingBox;
use geo::algorithm::contains::Contains;

mod clip;
use clip::{clip_to_bbox,clip_geometry_to_tiles};

mod validity;
use validity::{is_valid, is_valid_skip_expensive};

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

pub struct Layers {
    layers: Vec<Layer>,
    global_maxzoom: u8,
    global_minzoom: u8,
}

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
                // FIXME buffer is making different layers offset at different amounts
                buffer: 0, //layer["properties"]["buffer-size"].as_i64().map(|x| x as u16).unwrap_or(0) as u16,
                table: layer["Datasource"]["table"].as_str().unwrap().to_owned(),
            })
            // FIXME finish this thing
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

pub fn generate_all(filename: &str, min_zoom: u8, max_zoom: u8, bbox: &BBox, dest_dir: &str, if_not_exists: bool, compress: bool, metatile_scale: u8) {
    let layers = Layers::from_file(filename);
    let dest_dir = Path::new(dest_dir);
    fs::create_dir_all(dest_dir).unwrap();

    let connection_pool = ConnectionPool::new(layers.get_all_connections());

    let mut started_current_zoom: Option<Instant> = None;

    write_tilejson(&layers, &connection_pool, dest_dir);

    let mut last_zoom = 255;
    let mut num_tiles_done: u64 = 0;
    let logging_freq: u32 = 8;
    let log_every: u64 = 2u64.pow(logging_freq);
    let mask: u64 = (1u64<<logging_freq) - 1;
    for metatile in MetatilesIterator::new_for_bbox_zoom(metatile_scale, bbox, min_zoom, max_zoom) {
        if metatile.zoom() < min_zoom { continue; }

        if num_tiles_done & mask == 0 && num_tiles_done > 0 {
            if let Some(t) = started_current_zoom {
                let duration = duration_to_float_secs(&t.elapsed());
                println!("    Zoom {:2}, done {:4} ({:3}*{}) metatiles in {:>8} ({:9.4} metatiles/sec, {:9.4} tiles/sec)", last_zoom, num_tiles_done, num_tiles_done>>logging_freq, log_every, fmt_duration(&t.elapsed()), (num_tiles_done as f64)/duration, (num_tiles_done*(metatile_scale as u64)) as f64/duration );
            }
        }
        if metatile.zoom() != last_zoom {
            if let Some(t) = started_current_zoom {
                let duration = duration_to_float_secs(&t.elapsed());
                println!("Zoom {:2}, {:4} metatile(s), done in {:>8} ({:9.4} metatiles/sec, {:9.4} tiles/sec)", last_zoom, num_tiles_done, fmt_duration(&t.elapsed()), (num_tiles_done as f64)/duration, (num_tiles_done*(metatile_scale as u64)) as f64/duration );
            }
            started_current_zoom = Some(Instant::now());
            last_zoom = metatile.zoom();
            num_tiles_done = 0;
        }
        if metatile.zoom() > max_zoom {
            break;
        }


        let tiles = single_metatile(&layers, &metatile, &connection_pool);

        assert!(!if_not_exists);
        // FIXME the tile exists check needs to work with metatiles
        //if if_not_exists && filename.exists() {
        //    continue;
        //}

        for (tile, pbf) in tiles.into_iter() {
            let filename = dest_dir.join(tile.ts_path("pbf"));
            fs::create_dir_all(filename.parent().unwrap()).unwrap();

            pbf.write_to_file(filename.to_str().unwrap());
        }
        num_tiles_done += 1;


    }
}

fn write_tilejson(layers: &Layers, connection_pool: &ConnectionPool, dest: &Path) {
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

fn simplify_geom(geom: Geometry<i64>, tolerance: i64) -> Geometry<i64> {
    match geom {
        // Can't simplify a Point. let's hope this doesn't do a copy or memory move or something
        Geometry::Point(p) => Geometry::Point(p),
        Geometry::MultiPoint(p) => Geometry::MultiPoint(p),

        Geometry::LineString(p) => Geometry::LineString(p.simplify(&tolerance)),
        Geometry::MultiLineString(p) => Geometry::MultiLineString(p.simplify(&tolerance)),
        Geometry::Polygon(p) => Geometry::Polygon(p.simplify(&tolerance)),
        Geometry::MultiPolygon(p) => Geometry::MultiPolygon(p.simplify(&tolerance)),

        Geometry::GeometryCollection(_) => unimplemented!(),
    }
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
        let buffer = (layer.buffer as i64) * 16;

        let layer_name = &layer.id;

        let conn = connection_pool.connection_for_layer(&layer_name);
        
        let table = &layer.table;
        // FIXME should this be 4096??
        // FIXME not confident about this calculation.
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

        let minx = ll.0 as f64 - buffer_width;
        let miny = ll.1 as f64 - buffer_height;
        let maxx = ur.0 as f64 + buffer_width;
        let maxy = ur.1 as f64 + buffer_height;

        let bbox = format!("ST_SetSRID(ST_MakeBox2D(ST_Point({llx}, {lly}), ST_Point({urx}, {ury})), 3857)", llx=minx, lly=miny, urx=maxx, ury=maxy);
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

            let geom: geo::Geometry<f64> = match wkb::wkb_to_geom(&mut wkb_bytes.as_slice()) {
                Err(e) => {
                    // TODO investigate this more
                    eprintln!("Metatile: {:?} WKB reading error {:?}, first few bytes geom: {:?}", metatile, e, wkb_bytes.into_iter().take(20).collect::<Vec<u8>>());
                    continue;
                },
                Ok(g) => g,
            };

            let mut bad_obj = false;
            //if let Geometry::Polygon(ref x) = geom {
            //    // at z10 southside dublin
            //    // 39_469 good
            //    // 39_470 bad
            //    // the bad obj is 39_469
            //    if metatile.zoom() == 10 && num_objects >= 39_470 {
            //        continue;
            //    }
            //    bad_obj = metatile.zoom() == 10 && num_objects == 39_469;
            //}
            //if let Geometry::MultiPolygon(_) = geom {
            ////    //println!("Is mulitpolygon");
            ////    //println!("Skipping multipolygon");
            //    continue;
            //}

            // TODO not sure about this
            let pixel_size: f64 = tile_width/extent;

            // TODO there are a lot of calls to `is_valid`, which is computationally expensive, but
            // removing them makes it slower, probably because of the lots of invalid geoms

            // Convert to integer geom. MVT tiles have integers in i32, but I am getting i32
            // overflows when simplifying, so initally convert it to i64, and then convert it back.
            // it's a little poor since there are duplicate data.
            // FIXME if this is a point, maybe don't do the double change.
            let mut geom = remap_geometry(geom, minx, maxx, miny, maxy, extent);
            if geom.is_none() { continue; }
            let mut geom = geom.unwrap();
            validity::ensure_polygon_orientation(&mut geom);
            if ! is_valid_skip_expensive(&geom) {
                continue;
            }

            // TODO after converting to integer, maybe run a simple algorithm that removes points
            // which are on the line, i.e. A-B-C is straight line, so remove B. This keeps the
            // shape the same, but could reduce the number of points, which means other algorithms
            // will have less work to do.

            let simplification: i64 = if metatile.zoom() == layers.global_maxzoom { 1 } else { 8 };
            let geom = simplify_geom(geom, simplification);


            if ! is_valid(&geom) {
                continue;
            }

            // clip geometry, so no part of it goes outside the bbox. PostgreSQL will return
            // anything that overlaps.
            let geom = match clip_to_bbox(Cow::Owned(geom), &geo::Bbox{ xmin: -(buffer as i64), xmax: extent as i64 + buffer as i64, ymin: -(buffer as i64), ymax: extent as i64 + buffer as i64 }) {
                None => {
                    // geometry is outside the bbox, so skip
                    continue;
                },
                Some(g) => g,
            };
                    
            let mut properties = mapbox_vector_tile::Properties::new();

            for column in columns.iter() {
                let name = column.name();

                // Sometimes a NULL value can be returned, hence the dance with Option<Value>
                let value: Option<mapbox_vector_tile::Value> = match column.type_().name() {
                    "float4" => row.get_opt(name).map(|x| x.ok().map(mapbox_vector_tile::Value::Float)).unwrap_or(None),
                    "float8" => row.get_opt(name).map(|x| x.ok().map(mapbox_vector_tile::Value::Double)).unwrap_or(None),
                    "text" => row.get_opt(name).map(|x| x.ok().map(mapbox_vector_tile::Value::String)).unwrap_or(None),
                    "int4" => row.get_opt(name).map(|x| x.ok().map(|y| { let val: i32 = y; mapbox_vector_tile::Value::Int(val as i64) })).unwrap_or(None),
                    
                    // FIXME not 100% sure numeric is correct here
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

            // In cases where there is only one geometry here, we don't want to clone the
            // `properties`, and instead move it. If there are N geometries, we want to do N-1
            // clones, and 1 move. Hence the duplication with the loop.
            //
            // We want to do this in order, so we reverse the vec, and the pop from the end
            // (which is the original front).
            // FIXME rather than filtering out invalid geoms here, prevent the clipping code from
            // generating invalid geoms in the first place
            // One error was creating a linestring with 2 points, both the same
            let mut geoms: Vec<_> = clip_geometry_to_tiles(&metatile, geom, buffer).into_iter().filter_map(
                |(t, g)| match g {
                    Some(mut g) => {
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

            loop {
                if geoms.len() <= 1 { break; }
                if let Some((tile, geom)) = geoms.pop() {

                    let i = (tile.x() - metatile.x()) as i64;
                    let j = (tile.y() - metatile.y()) as i64;

                    // Here we convert it back to i32
                    let geom: Geometry<i32> = geom.map_coords(&|&(x, y)| ( (x - (4096*i)) as i32, (y - (4096*j)) as i32 ));

                    if is_valid(&geom) {

                        let feature = mapbox_vector_tile::Feature::new(geom, properties.clone());
                        let i = ((tile.x() - metatile.x())*scale + (tile.y() - metatile.y())) as usize;
                        let mvt_tile = results.get_mut(i).unwrap();
                        mvt_tile.add_feature(&layer_name, feature);

                    }
                }
            }

            if geoms.is_empty() { continue; }
            if let Some((tile, geom)) = geoms.pop() {
                let i = (tile.x() - metatile.x()) as i64;
                let j = (tile.y() - metatile.y()) as i64;

                // Here we convert it back to i32
                let mut geom: Geometry<i32> = geom.map_coords(&|&(x, y)| ( (x - (4096*i)) as i32, (y - (4096*j)) as i32 ));
                //validity::ensure_polygon_orientation(&mut geom);

                if is_valid(&geom) {
                    if bad_obj {
                        println!("geom {:?}", geom);
                    }
                    let feature = mapbox_vector_tile::Feature::new(geom, properties);
                    let i = ((tile.x() - metatile.x())*scale + (tile.y() - metatile.y())) as usize;
                    let mvt_tile = results.get_mut(i).unwrap();
                    mvt_tile.add_feature(&layer_name, feature);
                }
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

fn remap_linestring(ls: LineString<f64>, minx: f64, maxx: f64, miny: f64, maxy: f64, size: f64, should_be_ring: bool) -> Option<LineString<i64>> {
    let remap_xy = |x: f64, y: f64| -> (i64, i64) {
        (
            (((x - minx) / (maxx - minx))*size).round() as i64,
            (((maxy - y) / (maxy - miny))*size).round() as i64,
        )
    };

    let mut new_points: Vec<Point<i64>> = Vec::with_capacity(ls.0.len());
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

fn remap_geometry(geom: Geometry<f64>, minx: f64, maxx: f64, miny: f64, maxy: f64, size: f64) -> Option<Geometry<i64>> {
    let remap_xy = |x: f64, y: f64| -> (i64, i64) {
        (
            (((x - minx) / (maxx - minx))*size).round() as i64,
            (((maxy - y) / (maxy - miny))*size).round() as i64,
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
                    ).collect::<Vec<Point<i64>>>())))
            }
        },
        Geometry::LineString(ls) => {
            remap_linestring(ls, minx, maxx, miny, maxy, size, false).and_then(|ls| Some(Geometry::LineString(ls)))
        },
        Geometry::MultiLineString(mls) => {
            let mut res: Vec<LineString<i64>> = mls.0.into_iter().filter_map(|ls| remap_linestring(ls, minx, maxx, miny, maxy, size, false)).collect();
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
                    let interiors: Vec<LineString<i64>> = interiors.into_iter().filter_map(|int| remap_linestring(int, minx, maxx, miny, maxy, size, true)).collect();
                    Some(Geometry::Polygon(Polygon::new(exterior, interiors)))
                }
            }
        }
        Geometry::MultiPolygon(mp) => {
            let mut res: Vec<Polygon<i64>> = mp.0.into_iter().filter_map(|p| {
                let Polygon{ exterior, interiors } = p;
                match remap_linestring(exterior, minx, maxx, miny, maxy, size, true) {
                        // Exterior gets simplified away
                    None => None,
                    Some(exterior) => {
                        let interiors: Vec<LineString<i64>> = interiors.into_iter().filter_map(|int| remap_linestring(int, minx, maxx, miny, maxy, size, true)).collect();
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

