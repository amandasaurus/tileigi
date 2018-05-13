extern crate slippy_map_tiles;
#[macro_use]
extern crate clap;

extern crate log;
extern crate env_logger;

extern crate tileigi;

use std::path::{PathBuf, Path};
use std::io::Write;

use clap::{Arg, App, AppSettings, ArgGroup};
use slippy_map_tiles::BBox;
use env_logger::{Color};
use log::Level;

use tileigi::*;

fn main() {
    env_logger::Builder::from_default_env()
        .format(|buf, record| {
            let level = record.level();
            let mut level_style = buf.style();

            match level {
                Level::Trace => level_style.set_color(Color::White),
                Level::Debug => level_style.set_color(Color::Blue),
                Level::Info => level_style.set_color(Color::Green),
                Level::Warn => level_style.set_color(Color::Yellow),
                Level::Error => level_style.set_color(Color::Red).set_bold(true),
            };

            write!(buf, "{:>5} ", level_style.value(level)).unwrap();

            let ts = buf.timestamp();
            write!(buf, "{}: ", ts).unwrap();
            write!(buf, "{}: ", record.module_path().unwrap_or("UNKNOWN_MOD")).unwrap();
            write!(buf, "{}:{} ", record.file().unwrap_or("UNKNOWN_FILE"), record.line().map(|l| format!("{}", l)).unwrap_or("UNKNOWN_LINE".to_string())).unwrap();

            writeln!(buf, "{}", record.args()).unwrap();

            Ok(())
        })
        .init();

    let matches = App::new("tileigi")
        .version(crate_version!())
        .about("Generate vector tiles from a yml file")
        .setting(AppSettings::AllowLeadingHyphen)
        .arg(Arg::with_name("data_yml").long("data-yml").takes_value(true).value_name("FILENAME").required(true).help("Filename of the .yml file"))

        .arg(Arg::with_name("dest_dir").long("dest-dir").takes_value(true).value_name("DIR").help("Save tiles to this mbtiles file"))
        .arg(Arg::with_name("dest_mbtiles").long("dest-mbtiles").takes_value(true).value_name("FILENAME").help("Save tiles to this TileStash directory path"))
        .arg(Arg::with_name("dest_modtile").long("dest-modtile").takes_value(true).value_name("DIR").help("Save tiles to this mod_tile directory path"))
        .group(ArgGroup::with_name("dest").args(&["dest_dir", "dest_mbtiles", "dest_modtile"]).required(true))

        .arg(Arg::with_name("minzoom").long("minzoom").value_name("ZOOM").default_value("0").help("Minimum zoom to generate"))
        .arg(Arg::with_name("maxzoom").long("maxzoom").value_name("ZOOM").default_value("14").help("Maximum zoom to generate"))

        .arg(Arg::with_name("zoom").long("zoom").value_name("ZOOM").conflicts_with_all(&["minzoom", "maxzoom"]).help("Only generate for this zoom"))

        .arg(Arg::with_name("bbox").long("bbox").takes_value(true).value_name("MINLON,MINLAT,MAXLON,MAXLAT").help("Only generate tiles inside this bbox. 'planet' for planet, or minlon,minlat,maxlon,maxlat"))

        .arg(Arg::with_name("bbox-bottom").long("bbox-bottom").takes_value(true).value_name("DEGREES").help("BBox, bottom"))
        .arg(Arg::with_name("bbox-top").long("bbox-top").takes_value(true).value_name("DEGREES").help("BBox, top"))
        .arg(Arg::with_name("bbox-left").long("bbox-left").takes_value(true).value_name("DEGREES").help("BBox, left"))
        .arg(Arg::with_name("bbox-right").long("bbox-right").takes_value(true).value_name("DEGREES").help("BBox, right"))
        .group(ArgGroup::with_name("bbox_individual").args(&["bbox-bottom", "bbox-top", "bbox-left", "bbox-right"]).conflicts_with("bbox").multiple(true))

        .arg(Arg::with_name("metatile-scale").long("metatile-scale").default_value("8").value_name("NUMBER").help("Size of metatile to use (8x8 default)"))
        .arg(Arg::with_name("threads").long("threads").default_value("1").value_name("NUBMER").help("Number of concurrent generation threads to run"))

        .arg(Arg::with_name("iter_mode").long("mode").default_value("tile-then-layer").possible_values(&["tile-then-layer", "layer-then-tile"]))

        .arg(Arg::with_name("if_not_exists").long("if-not-exists").help("Do not generate a tile if the file already exists. Doesn't work with mbtiles (yet)"))
        .arg(Arg::with_name("no_compress").long("no-compress").help("Do not compress the pbf files"))

        .arg(Arg::with_name("file-writer-buffer").long("file-writer-buffer").help("Size of buffer for the file writer thread").takes_value(true))

        .arg(Arg::with_name("tile_list")
             .long("tile-list").alias("list")
             .takes_value(true).required(false).value_name("FILENAME")
             .validator(|s| { if Path::new(&s).exists() { Ok(()) } else { Err(format!("File {} not found", s)) }})
             .help("Generate tiles from a list of tiles, one metatile per line 'SCALE Z/X/Y'"))
        .get_matches();

    let data_yml = matches.value_of("data_yml").unwrap();

    let dest = match (matches.value_of("dest_dir"), matches.value_of("dest_mbtiles"), matches.value_of("dest_modtile")) {
        (Some(dest_dir), None, None) => TileDestinationType::TileStashDirectory(PathBuf::from(dest_dir)),
        (None, Some(mbtiles_filename), None) => TileDestinationType::MBTiles(PathBuf::from(mbtiles_filename)),
        (None, None, Some(modtile_dir)) => TileDestinationType::ModTileDirectory(PathBuf::from(modtile_dir)),
        (None, None, None) => panic!("Must provide a destination"),
        _ => panic!("Can't provide >1 dest"),
    };

    let (minzoom, maxzoom): (u8, u8) = if let Some(z) = matches.value_of("zoom") {
        let z: u8 = z.parse().unwrap();
        (z, z)
    } else {
        ( matches.value_of("minzoom").unwrap().parse().unwrap(),
          matches.value_of("maxzoom").unwrap().parse().unwrap() )
    };

    let if_not_exists = matches.is_present("if_not_exists");
    let compress = ! matches.is_present("no_compress");
    let metatile_scale: u8 = matches.value_of("metatile-scale").unwrap().parse().unwrap();
    let num_threads: usize = matches.value_of("threads").unwrap().parse().unwrap();

    let bbox: Option<BBox> = match matches.value_of("bbox") {
        Some("planet") => None,
        Some(bbox_string) => Some(BBox::new_from_string(bbox_string).expect("Invalid bbox")),
        None => {
            if matches.is_present("bbox-top") {
                // if we have one, we presume we have them all
                Some(BBox::new(
                        matches.value_of("bbox-top").unwrap_or("90.0").parse().unwrap(),
                        matches.value_of("bbox-left").unwrap_or("-180.0").parse().unwrap(),
                        matches.value_of("bbox-bottom").unwrap_or("-90.0").parse().unwrap(),
                        matches.value_of("bbox-right").unwrap_or("180.0").parse().unwrap(),
                    ).unwrap())
            } else {
                None
            }
        },
    };

    let tile_list: Option<String> = matches.value_of("tile_list").map(|s| s.to_string());

    let file_writer_buffer: usize = matches.value_of("file-writer-buffer").map(|s| s.parse().unwrap()).unwrap_or(1_000_000);

    match matches.value_of("iter_mode") {
        Some("tile-then-layer") => {
            generate_all(&data_yml, minzoom, maxzoom, &bbox, &dest, if_not_exists, compress, metatile_scale, num_threads, tile_list, file_writer_buffer);
        },
        Some("layer-then-tile") => {
            generate_by_layer(&data_yml, minzoom, maxzoom, &bbox, &dest, if_not_exists, compress, metatile_scale, num_threads, tile_list, file_writer_buffer);
        },
        _ => panic!(),
    }
}
