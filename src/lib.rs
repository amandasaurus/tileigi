#![allow(dead_code,unused_imports,unused_variables,unused_mut)]

extern crate postgres;
extern crate slippy_map_tiles;
extern crate yaml_rust;
extern crate mapbox_vector_tile;
extern crate wkb;
extern crate geo;

use std::fs::File;
use std::fs;
use std::io::prelude::*;
use std::path::Path;
use std::collections::{HashSet, HashMap};
use std::time::Instant;

use yaml_rust::{YamlLoader, Yaml};

use postgres::{Connection, TlsMode};
use postgres::params::ConnectParams;

use slippy_map_tiles::BBox;

use geo::*;
use geo::algorithm::simplify::Simplify;
use geo::algorithm::map_coords::MapCoords;
use geo::algorithm::boundingbox::BoundingBox;
use geo::algorithm::contains::Contains;

// local stuff
mod cohen_sutherland;
mod sutherland_hodgeman;
//mod weiler_atherton;

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
    layers: Vec<Yaml>,
    global_maxzoom: u8,
    global_minzoom: u8,
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

        let layers: Vec<_> = layers.into_iter().filter(
            |l| l["Datasource"]
                .as_hash()
                .and_then(
                    |d| d.get(&Yaml::String("type".to_string()))
                        .and_then(|t| Some(t.as_str() == Some("postgis")))
                )
                .unwrap_or(false)
            ).cloned().collect();

        Layers{ layers: layers, global_minzoom: global_minzoom, global_maxzoom: global_maxzoom }

    }

    fn get_all_connections(&self) -> HashMap<ConnectParams, Vec<String>> {
        let mut conns = HashMap::new();
        for layer in self.layers.iter() {

            let mut conn_params = postgres::params::Builder::new();
            // TODO do others
            if let Some(dbname) = layer["Datasource"]["dbname"].as_str() {
                conn_params.database(dbname);
            }

            // TODO fix user
            conn_params.user("rory", None);
            // TODO read hosts
            let conn_params = conn_params.build(postgres::params::Host::Tcp("localhost".to_string()));

            if ! conns.contains_key(&conn_params) {
                conns.insert(conn_params.clone(), Vec::new());
            }

            let layer_id: String = layer["id"].as_str().unwrap().to_string();
            conns.get_mut(&conn_params).unwrap().push(layer_id);
            
        }

        conns
    }
}

fn fmt_duration(dur: &std::time::Duration) -> String {
    format!("{:.2}s", (dur.as_secs() as f64) + (dur.subsec_nanos() as f64 / 1e9))
}

pub fn generate_all(filename: &str, max_zoom: u8, bbox: &BBox, dest_dir: &str) {
    let layers = Layers::from_file(filename);
    let dest_dir = Path::new(dest_dir);

    let connection_pool = ConnectionPool::new(layers.get_all_connections());

    let mut started_current_zoom: Option<Instant> = None;

    let mut last_zoom = 255;
    let mut num_tiles_done = 0;
    for tile in bbox.tiles() {
        if tile.zoom() != last_zoom {
            if let Some(t) = started_current_zoom {
                println!("Zoom {}, {} tile(s), done in {}", last_zoom, num_tiles_done, fmt_duration(&t.elapsed()));
            }
            started_current_zoom = Some(Instant::now());
            last_zoom = tile.zoom();
            num_tiles_done = 0;
        }
        if tile.zoom() > max_zoom {
            break;
        }

        let filename = dest_dir.join(tile.ts_path("pbf"));

        fs::create_dir_all(filename.parent().unwrap()).unwrap();

        let pbf = single_tile(&layers, &tile, &connection_pool);
        pbf.write_to_file(filename.to_str().unwrap());
        num_tiles_done += 1;
    }

}

fn simplify_geom(geom: Geometry<f64>, tolerance: f64) -> Geometry<f64> {
    match geom {
        // Can't simplify a Point. let's hope this doesn't do a copy or memory move or something
        Geometry::Point(p) => Geometry::Point(p),
        Geometry::MultiPoint(p) => Geometry::MultiPoint(p),

        Geometry::LineString(ref p) => Geometry::LineString(p.simplify(&tolerance)),
        Geometry::MultiLineString(ref p) => Geometry::MultiLineString(p.simplify(&tolerance)),
        Geometry::Polygon(ref p) => Geometry::Polygon(p.simplify(&tolerance)),
        Geometry::MultiPolygon(ref p) => Geometry::MultiPolygon(p.simplify(&tolerance)),
        Geometry::GeometryCollection(_) => unimplemented!(),
    }
}

fn clip_line<T: CoordinateType>(line: LineString<T>, bbox: &Bbox<T>) -> Option<Geometry<T>> {
    // This is a very simple approach, just split the line into many little segments and clip them
    // all.
    let mut segments: Vec<LineString<T>> = line.0.windows(2).filter_map(|points|
        match cohen_sutherland::clip(&points[0], &points[1], bbox) {
            None => None,
            Some((new_p1, new_p2)) =>
                Some(LineString(vec![Point::new(new_p1.0, new_p1.1), Point::new(new_p2.0, new_p2.1)]))
    }).collect();

    if segments.len() == 0 {
        None
    } else if segments.len() == 1 {
        Some(segments.remove(0).into())
    } else {
        Some(MultiLineString(segments).into())
    }
}


fn clip_to_bbox<T: CoordinateType>(geom: Geometry<T>, bbox: &Bbox<T>) -> Option<Geometry<T>> {
    match geom {
        Geometry::Point(p) => {
            if bbox.contains(&p) {
                Some(Geometry::Point(p))
            } else {
                None
            }
        },
        Geometry::LineString(l) => clip_line(l, bbox),
        Geometry::MultiLineString(ml) => {
            let mut lines = Vec::with_capacity(ml.0.len());
            for clipped_line in ml.0.into_iter().filter_map(|l| clip_line(l, bbox)) {
                match clipped_line {
                    Geometry::LineString(l) => lines.push(l),
                    Geometry::MultiLineString(mut mls) => lines.append(&mut mls.0),
                    _ => unreachable!(),
                }
            }

            if lines.len() == 0 {
                None
            } else {
                Some(MultiLineString(lines).into())
            }
        }
        _ => unimplemented!(),
    }
}



pub fn single_tile(layers: &Layers, tile: &slippy_map_tiles::Tile, connection_pool: &ConnectionPool) -> mapbox_vector_tile::Tile {
    //println!("Tile {:?}", tile);

    let mut result_tile = mapbox_vector_tile::Tile::new();
    for layer in layers.layers.iter() {
        let minzoom = layer["properties"]["minzoom"].as_i64().map(|x| x as u8).unwrap_or(layers.global_minzoom) as u8;
        let maxzoom = layer["properties"]["maxzoom"].as_i64().map(|x| x as u8).unwrap_or(layers.global_maxzoom) as u8;
        let maxzoom = if maxzoom > layers.global_maxzoom { layers.global_maxzoom } else { maxzoom };

        // Skip layers which are not on this zoom
        if tile.zoom() < minzoom || tile.zoom() > maxzoom {
            continue;
        }

        let layer_name = layer["id"].as_str().unwrap().clone();

        let conn = connection_pool.connection_for_layer(layer_name);
        
        let table = layer["Datasource"]["table"].as_str().unwrap();
        // FIXME should this be 4096??
        // FIXME not confident about this calculation.
        let canvas_size = 256.;
        let ll = tile.sw_corner().to_3857();
        let ur = tile.ne_corner().to_3857();

        // TODO vtiles have y positive going down, is this correct??
        let minx = ll.0 as f64;
        let miny = ll.1 as f64;
        let maxx = ur.0 as f64;
        let maxy = ur.1 as f64;

        // FIXME include buffer in bbox
        let bbox = format!("ST_SetSRID(ST_MakeBox2D(ST_Point({llx}, {lly}), ST_Point({urx}, {ury})), 3857)", llx=ll.0, lly=ll.1, urx=ur.0, ury=ur.1);
        let tile_width = (ur.0 - ll.0) as f64;
        let tile_height = (ur.1 - ll.1) as f64;
        assert!(tile_height > 0.);
        assert!(tile_width > 0.);

        let pixel_width = format!("{}", tile_width / canvas_size);
        let pixel_height = format!("{}", tile_height / canvas_size);
        let table = table.replace("!pixel_width!", &pixel_width).replace("!pixel_height!", &pixel_height).replace("!bbox!", &bbox);
        if table.contains("!scale_denominator!") {
            // FIXME implement this
            //println!("This query has scale_denominator, skipping.");
            continue;
        }
        let mut new_layer = mapbox_vector_tile::Layer::new(layer_name.to_string());
        let query = format!("SELECT ST_AsBinary(way), * from {table} where way && {bbox}", table=table, bbox=bbox);
        let res = conn.query(&query, &[]).unwrap();

        if res.is_empty() {
            continue;
        }

        let extent = new_layer.extent as f64;

        let columns: Vec<_> = res.columns().iter().skip(1).filter(|c| c.name() != "way").collect();

        for row in res.iter() {
            // First object is the ST_AsBinary
            // TODO Does this do any copies that we don't want?
            let wkb_bytes: Vec<u8> = row.get(0);

            let geom: geo::Geometry<f64> = wkb::wkb_to_geom(wkb_bytes.as_slice());


            // TODO not sure abou this
            let pixel_size: f64 = tile_width/extent;
            let simplification: f64 = if tile.zoom() == layers.global_maxzoom { pixel_size } else { 8. * pixel_size };
            let geom = simplify_geom(geom, simplification);

            // TODO there might be a problem with y direction here??
            let geom: geo::Geometry<i32> = geom.map_coords(&|&(x, y)|
                (
                    (((x - minx) / (maxx - minx))*extent).round() as i32,
                    (((y - miny) / (maxy - miny))*extent).round() as i32,
                ));

            // clip geometry
            let geom = match clip_to_bbox(geom, &geo::Bbox{ xmin: 0, xmax: new_layer.extent as i32, ymin: 0, ymax: new_layer.extent as i32 }) {
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

            let feature = mapbox_vector_tile::Feature::new(geom, properties);
            new_layer.add_feature(feature);

        }
        
        result_tile.layers.push(new_layer);

    }

    result_tile

}

mod test {
    use super::*;

    #[test]
    fn clip_point() {
        let bbox = Bbox{ xmin: 0., ymin: 0., xmax: 4096., ymax: 4096.};

        assert_eq!(clip_to_bbox(Point::new(0., 0.).into(), &bbox), Some(Point::new(0., 0.).into()));
        assert_eq!(clip_to_bbox(Point::new(-1., -1.).into(), &bbox), None);
        assert_eq!(clip_to_bbox(Point::new(4000., 0.).into(), &bbox), Some(Point::new(4000., 0.).into()));
        assert_eq!(clip_to_bbox(Point::new(5000., 0.).into(), &bbox), None);
    }

    #[test]
    fn clip_linestring() {
        let bbox = Bbox{ xmin: 0., ymin: 0., xmax: 4096., ymax: 4096.};

        let p1 = Point::new(10., 10.);
        let p2 = Point::new(20., 20.);
        let p3 = Point::new(10., 5000.);
        let p4 = Point::new(5000., 5000.);
        let p5 = Point::new(10., 1000.);

        assert_eq!(clip_to_bbox(LineString(vec![p1, p2]).into(), &bbox), Some(LineString(vec![p1, p2]).into()));
        assert_eq!(clip_to_bbox(LineString(vec![p4, p4]).into(), &bbox), None);
        assert_eq!(clip_to_bbox(LineString(vec![p1, p5]).into(), &bbox), Some(LineString(vec![p1, p5]).into()));
        assert_eq!(clip_to_bbox(LineString(vec![p1, p3]).into(), &bbox), Some(LineString(vec![p1, Point::new(10., 4096.)]).into()));

        assert_eq!(clip_to_bbox(LineString(vec![p1, p2, p5]).into(), &bbox),
            Some(MultiLineString(vec![
                LineString(vec![p1, p2]),
                LineString(vec![p2, p5]),
                                ]).into()));

        assert_eq!(clip_to_bbox(LineString(vec![p1, p2, p5, p4]).into(), &bbox),
            Some(MultiLineString(vec![
                LineString(vec![p1, p2]),
                LineString(vec![p2, p5]),
                LineString(vec![p5, Point::new(3872.26, 4096.)]),
                                ]).into()));
    }
    
    #[test]
    fn clip_multilinestring() {
        let bbox = Bbox{ xmin: 0., ymin: 0., xmax: 4096., ymax: 4096.};

        let p1 = Point::new(10., 10.);
        let p2 = Point::new(20., 20.);
        let p3 = Point::new(10., 5000.);
        let p4 = Point::new(5000., 5000.);
        let p5 = Point::new(10., 1000.);

        assert_eq!(clip_to_bbox(
                MultiLineString(vec![
                    LineString(vec![p1, p2, p5]),
                    LineString(vec![p1, p2, p5, p4])
                ]).into(), &bbox),
            Some(MultiLineString(vec![
                LineString(vec![p1, p2]),
                LineString(vec![p2, p5]),
                LineString(vec![p1, p2]),
                LineString(vec![p2, p5]),
                LineString(vec![p5, Point::new(3872.26, 4096.)]),
                                ]).into()));
    }

}
