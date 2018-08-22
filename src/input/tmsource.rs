//! Read tmsource (aka tm2source) data.yml files

use yaml_rust::{YamlLoader, Yaml};
use std::fs::File;
use std::io::prelude::*;
use std::fs;

use super::{Layers, Layer, TableSQL};

type Result<T> = std::result::Result<T, failure::Error>;

pub fn layers_from_file(filename: &str) -> Result<Layers> {
    let mut file = File::open(filename)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;

    let mut data_yml = YamlLoader::load_from_str(&contents)?;
    let data_yml = data_yml.remove(0);

    let global_maxzoom = data_yml["maxzoom"].as_i64().ok_or(format_err!("maxzoom is not a i64"))? as u8;
    let global_minzoom = data_yml["minzoom"].as_i64().ok_or(format_err!("minzoom is not a i64"))? as u8;
    let bounds = data_yml["bounds"].as_vec().ok_or(format_err!("bounds is not an array"))?;
    // '90' would be parsed as an int in yaml
    let bounds: [f64; 4] = [
        bounds[0].as_f64().or_else(|| bounds[0].as_i64().map(|i| i as f64)).ok_or(format_err!("bounds invalid"))?,
        bounds[1].as_f64().or_else(|| bounds[1].as_i64().map(|i| i as f64)).ok_or(format_err!("bounds       invalid"))?,
        bounds[2].as_f64().or_else(|| bounds[2].as_i64().map(|i| i as f64)).ok_or(format_err!("bounds       invalid"))?,
        bounds[3].as_f64().or_else(|| bounds[3].as_i64().map(|i| i as f64)).ok_or(format_err!("bounds       invalid"))?,
    ];

    let center = data_yml["center"].as_vec().ok_or(format_err!("center is not an array"))?;
    let center: [f64; 3] = [
        center[0].as_f64().or_else(|| center[0].as_i64().map(|i| i as f64)).ok_or(format_err!("center invalid"))?,
        center[1].as_f64().or_else(|| center[1].as_i64().map(|i| i as f64)).ok_or(format_err!("center       invalid"))?,
        center[2].as_f64().or_else(|| center[2].as_i64().map(|i| i as f64)).ok_or(format_err!("center       invalid"))?,
    ];

    let name = data_yml["name"].as_str().ok_or(format_err!("name is not str"))?.to_string();
    let description = data_yml["description"].as_str().ok_or(format_err!("description is not str"))?.to_string();

    // rust-yaml really needs an into_hash (etc)
    // clone all the things
    let layers = data_yml.as_hash().ok_or(format_err!("yml file is not a hash"))?.clone().remove(&Yaml::String("Layer".to_string())).ok_or(format_err!("No Layer key"))?;
    let layers = layers.as_vec().ok_or(format_err!("layers is not an array"))?;

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
            let table = layer["Datasource"]["table"].as_str().ok_or(format_err!("table is not a str"))?;
            let table = TableSQL::new(table.to_owned());
            
            Ok(Layer {
                id: layer["id"].as_str().ok_or(format_err!("id for layer is not a str"))?.to_owned(),
                dbname: layer["Datasource"]["dbname"].as_str().map(|x| x.to_owned()),
                minzoom: layer["properties"]["minzoom"].as_i64().map(|x| x as u8).unwrap_or(global_minzoom) as u8,
                maxzoom: layer["properties"]["maxzoom"].as_i64().map(|x| x as u8).unwrap_or(global_maxzoom) as u8,
                buffer: layer["properties"]["buffer-size"].as_i64().map(|x| x as u16).unwrap_or(0) as u16,
                table: table,
            })
        })
        .collect::<Result<Vec<Layer>>>()?;

    Ok(Layers{ layers: layers, global_minzoom: global_minzoom, global_maxzoom: global_maxzoom, bounds: Some(bounds), center: Some(center), name: Some(name), description: Some(description) })

}
