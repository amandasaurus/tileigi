use postgres::params::ConnectParams;
use postgres::types::{Type, ToSql, IsNull};
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::fs;

use LocalBBox;

type Result<T> = std::result::Result<T, failure::Error>;

mod tmsource;

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
    pub fn from_file(filename: &str) -> Result<Self> {
        tmsource::layers_from_file(filename)
    }

    pub fn from_tmsource_file(filename: &str) -> Result<Self> {
        tmsource::layers_from_file(filename)
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
    pub has_zoom: bool,
        
}

impl TableSQL {
    pub fn new(query: String) -> Self {
        let has_pixel_width = query.contains("!pixel_width!");
        let has_pixel_height = query.contains("!pixel_height!");
        let has_scale_denominator = query.contains("!scale_denominator!");
        let has_zoom = query.contains("!zoom!") || query.contains("!ZOOM!");

        let mut query = query;

        query = query.replace("!bbox!", "$1");
        query = query.replace("!BBOX!", "$1");

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
        if has_zoom {
            query = query.replace("!zoom!", &format!("${}", param_num));
            query = query.replace("!ZOOM!", &format!("${}", param_num));
            param_num += 1;
        }
                
        let query = format!("SELECT ST_AsBinary(way), * from {} where way && $1", query);
        TableSQL{
            query, has_pixel_width, has_pixel_height, has_scale_denominator, has_zoom,
        }
    }

    pub fn params<'a, T: num_traits::Float+Into<f64>+'a+std::fmt::Debug>(&self, bbox: &'a LocalBBox<T>, pixel_width: &'a f32, pixel_height: &'a f32, zoom: &'a i32, scale_denominator: &'a f32) -> Vec<&'a postgres::types::ToSql> {
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
        if self.has_zoom {
            results.push(zoom);
        }

        results
    }
}
