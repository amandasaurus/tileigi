extern crate slippy_map_tiles;
extern crate clap;

extern crate tilegen;

use std::path::PathBuf;
use clap::{Arg, App, AppSettings};
use slippy_map_tiles::BBox;

use tilegen::*;

fn main() {

    let matches = App::new("test")
        .setting(AppSettings::AllowLeadingHyphen)
        .arg(Arg::with_name("data_yml").long("data-yml").takes_value(true).required(true).help("Filename of the .yml file"))

        .arg(Arg::with_name("dest_dir").long("dest-dir").takes_value(true).conflicts_with("dest-mbtiles").help("Save tiles to this mbtiles file"))
        .arg(Arg::with_name("dest_mbtiles").long("dest-mbtiles").takes_value(true).conflicts_with("dest_dir").help("Save tiles to this TileStash directory path"))
        .arg(Arg::with_name("dest_modtile").long("dest-modtile").takes_value(true).conflicts_with("dest-mbtiles").help("Save tiles to this mod_tile directory path"))

        .arg(Arg::with_name("minzoom").long("minzoom").default_value("0").help("Minimum zoom to generate"))
        .arg(Arg::with_name("maxzoom").long("maxzoom").default_value("14").help("Maximum zoom to generate"))

        .arg(Arg::with_name("bbox").long("bbox").takes_value(true).conflicts_with_all(&["bbox-bottom", "bbox-top", "bbox-left", "bbox-right"]).help("Only generate tiles inside this bbox. 'planet' for planet, or minlon,minlat,maxlon,maxlat"))

        .arg(Arg::with_name("bbox-bottom").long("bbox-bottom").takes_value(true).conflicts_with("bbox").requires_all(&["bbox-bottom", "bbox-top", "bbox-left", "bbox-right"]).help("BBox, bottom"))
        .arg(Arg::with_name("bbox-top").long("bbox-top").takes_value(true).conflicts_with("bbox").requires_all(&["bbox-bottom", "bbox-top", "bbox-left", "bbox-right"]).help("BBox, top"))
        .arg(Arg::with_name("bbox-left").long("bbox-left").takes_value(true).conflicts_with("bbox").requires_all(&["bbox-bottom", "bbox-top", "bbox-left", "bbox-right"]).help("BBox, left"))
        .arg(Arg::with_name("bbox-right").long("bbox-right").takes_value(true).conflicts_with("bbox").requires_all(&["bbox-bottom", "bbox-top", "bbox-left", "bbox-right"]).help("BBox, right"))

        .arg(Arg::with_name("if_not_exists").long("if-not-exists").help("Do not generate a tile if the file already exists. Doesn't work with mbtiles (yet)"))
        .arg(Arg::with_name("no_compress").long("no-compress").help("Do not compress the pbf files"))
        .arg(Arg::with_name("metatile-scale").long("metatile-scale").default_value("8").help("Size of metatile to use (8x8 default)"))
        .arg(Arg::with_name("threads").long("threads").default_value("1").help("Number of concurrent generation threads to run"))
        .get_matches();

    let data_yml = matches.value_of("data_yml").unwrap();

    let dest = match (matches.value_of("dest_dir"), matches.value_of("dest_mbtiles"), matches.value_of("dest_modtile")) {
        (Some(dest_dir), None, None) => TileDestinationType::TileStashDirectory(PathBuf::from(dest_dir)),
        (None, Some(mbtiles_filename), None) => TileDestinationType::MBTiles(PathBuf::from(mbtiles_filename)),
        (None, None, Some(modtile_dir)) => TileDestinationType::ModTileDirectory(PathBuf::from(modtile_dir)),
        (None, None, None) => panic!("Must provide a destination"),
        _ => panic!("Can't provide >1 dest"),
    };

    let minzoom: u8 = matches.value_of("minzoom").unwrap().parse().unwrap();
    let maxzoom: u8 = matches.value_of("maxzoom").unwrap().parse().unwrap();
    let if_not_exists = matches.is_present("if_not_exists");
    let compress = ! matches.is_present("no_compress");
    let metatile_scale: u8 = matches.value_of("metatile-scale").unwrap().parse().unwrap();
    let num_threads: usize = matches.value_of("threads").unwrap().parse().unwrap();

    let bbox: Option<BBox> = match matches.value_of("bbox") {
        Some("planet") => None,
        Some(bbox_string) => Some(BBox::new_from_string(bbox_string).expect("Invalid bbox")),
        None => None,
    };

    generate_all(&data_yml, minzoom, maxzoom, &bbox, &dest, if_not_exists, compress, metatile_scale, num_threads);
}
