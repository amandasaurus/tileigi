use postgres::params::ConnectParams;
use postgres::types::{Type, ToSql, IsNull};
use std::collections::HashMap;
use yaml_rust::{YamlLoader, Yaml};
use std::fs::File;
use std::io::prelude::*;
use std::fs;

use LocalBBox;

#[derive(Clone,Debug)]
pub struct Layers {
    pub layers: Vec<Layer>,
    pub global_maxzoom: u8,
    pub global_minzoom: u8,
    pub bounds: [f64; 4],
    pub center: [f64; 3],
    pub name: String,
    pub description: String,
}

#[derive(Clone,Debug)]
pub struct Layer {
    pub minzoom: u8,
    pub maxzoom: u8,
    pub buffer: u16,
    pub id: String,
    pub table: TableSQL,
    pub dbname: Option<String>,
}

impl Layers {
    pub fn from_file(filename: &str) -> Self {
        let mut file = File::open(filename).unwrap();
        let mut contents = String::new();
        file.read_to_string(&mut contents).unwrap();

        let mut data_yml = YamlLoader::load_from_str(&contents).unwrap();
        let data_yml = data_yml.remove(0);

        let global_maxzoom = data_yml["maxzoom"].as_i64().unwrap() as u8;
        let global_minzoom = data_yml["minzoom"].as_i64().unwrap() as u8;
        let bounds = data_yml["bounds"].as_vec().unwrap();
        // '90' would be parsed as an int in yaml
        let bounds: [f64; 4] = [
            bounds[0].as_f64().or_else(|| bounds[0].as_i64().map(|i| i as f64)).unwrap(),
            bounds[1].as_f64().or_else(|| bounds[1].as_i64().map(|i| i as f64)).unwrap(),
            bounds[2].as_f64().or_else(|| bounds[2].as_i64().map(|i| i as f64)).unwrap(),
            bounds[3].as_f64().or_else(|| bounds[3].as_i64().map(|i| i as f64)).unwrap(),
        ];

        let center = data_yml["center"].as_vec().unwrap();
        let center: [f64; 3] = [
            center[0].as_f64().or_else(|| center[0].as_i64().map(|i| i as f64)).unwrap(),
            center[1].as_f64().or_else(|| center[1].as_i64().map(|i| i as f64)).unwrap(),
            center[2].as_f64().or_else(|| center[2].as_i64().map(|i| i as f64)).unwrap(),
        ];

        let name = data_yml["name"].as_str().unwrap().to_string();
        let description = data_yml["description"].as_str().unwrap().to_string();

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

        Layers{ layers: layers, global_minzoom: global_minzoom, global_maxzoom: global_maxzoom, bounds: bounds, center: center, name: name, description: description }

    }

    pub fn get_all_connections(&self) -> HashMap<ConnectParams, Vec<String>> {
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

#[derive(Clone,Debug)]
pub struct TableSQL {
    pub query: String,
    pub has_pixel_width: bool,
    pub has_pixel_height: bool,
    pub has_scale_denominator: bool,
        
}

impl TableSQL {
    pub fn new(query: String) -> Self {
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

    pub fn params<'a, T: num_traits::Float+Into<f64>+'a+std::fmt::Debug>(&self, bbox: &'a LocalBBox<T>, pixel_width: &'a f32, pixel_height: &'a f32, scale_denominator: &'a f32) -> Vec<&'a postgres::types::ToSql> {
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

