//! Tegola input format

use std::fs::File;
use std::io::prelude::*;
use std::fs;

use super::{Layers, Layer, TableSQL};

type Result<T> = std::result::Result<T, failure::Error>;

#[derive(Deserialize, Debug)]
struct TegolaConfig {
    cache: Cache,
    websever: Option<WebserverConfig>,
    providers: Vec<Provider>,
    maps: Vec<Map>,
}

#[derive(Deserialize, Debug)]
struct WebserverConfig {
    port: String,
    hostname: String,
    cors_allowed_origin: String,
}

#[derive(Deserialize, Debug)]
struct Provider {
    name: String,

    #[serde(rename="type")]
    type_: String,

    host: String,
    port: u16,
    database: String,
    user: String,
    password: String,

    layers: Vec<TegolaProviderLayer>,
}

#[derive(Deserialize, Debug)]
struct TegolaProviderLayer {
    name: String,
    id_fieldname: String,
    geometry_fieldname: String,
    sql: String,
}

#[derive(Deserialize, Debug)]
struct Map {
    name: String,
    attribution: Option<String>,
    bounds: Option<[f64; 4]>,
    center: Option<[f64; 3]>,
    layers: Vec<TegolaMapLayer>,
}

#[derive(Deserialize, Debug)]
struct TegolaMapLayer {
    provider_layer: String,
    name: Option<String>,
    min_zoom: Option<u8>,
    max_zoom: Option<u8>,
    default_tags: Option<toml::Value>,
    dont_simplify: Option<bool>,
    geometry_type: Option<String>,
}

#[derive(Deserialize, Debug)]
struct Cache {
    #[serde(rename="type")]
    type_: String,

    max_zoom: Option<u8>,
    basepath: String,
}

pub fn layers_from_file(filename: &str) -> Result<Layers> {
    let mut file = File::open(filename)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    let tegola_config: TegolaConfig =  toml::from_str(&contents)?;

    if tegola_config.maps.len() != 1 {
        return Err(format_err!("Invalid number of maps. Must be exactly one, not {}", tegola_config.maps.len()));
    }

    let map = &tegola_config.maps[0];

    Ok(Layers{ 
        global_maxzoom: 14,
        global_minzoom: 0,
        bounds: None,
        name: None,
        description: None,
        center: map.center,
        layers: map.layers.iter().map(|l| {
            let (provider_name, provider_layer) = {
                let x = l.provider_layer.split(".").take(2).collect::<Vec<_>>();
                (x[0], x[1])
            };
            let provider = tegola_config.providers.iter().filter(|p| p.name == provider_name).nth(0).ok_or(format_err!("Missing provider: {}", l.provider_layer))?;
            let sql = &provider.layers.iter().filter(|l| l.name == provider_layer).nth(0).ok_or(format_err!("missing layer {}", l.provider_layer))?.sql;
            Ok(Layer {
                id: (l.name.to_owned()).ok_or(format_err!("Missing name"))?,
                dbname: Some(provider.database.to_owned()),
                minzoom: l.min_zoom.unwrap_or(0),
                maxzoom: l.max_zoom.unwrap_or(22),
                buffer: 0,
                table: TableSQL::new(format!("({}) as t", sql.to_owned())),
            })
        }).collect::<Result<Vec<Layer>>>()?,
    })

}
