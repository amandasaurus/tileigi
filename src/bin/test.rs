extern crate slippy_map_tiles;
extern crate clap;

extern crate tilegen;

use clap::{Arg, App, AppSettings};
use slippy_map_tiles::BBox;

use tilegen::*;

fn main() {

    let matches = App::new("test")
        .setting(AppSettings::AllowLeadingHyphen)
        .arg(Arg::with_name("data_yml").long("data-yml").takes_value(true).required(true))
        .arg(Arg::with_name("dest_dir").long("dest").takes_value(true).required(true))
        .arg(Arg::with_name("minzoom").long("minzoom").default_value("0"))
        .arg(Arg::with_name("maxzoom").long("maxzoom").default_value("14"))
        .arg(Arg::with_name("bbox").long("bbox").default_value("90,-180,-90,180"))
        .arg(Arg::with_name("if_not_exists").long("if-not-exists"))
        .arg(Arg::with_name("no_compress").long("no-compress"))
        .arg(Arg::with_name("metatile-scale").long("metatile-scale").default_value("8"))
        .arg(Arg::with_name("threads").long("threads").default_value("1"))
        .get_matches();

    let data_yml = matches.value_of("data_yml").unwrap();
    let dest_dir = matches.value_of("dest_dir").unwrap();
    let minzoom: u8 = matches.value_of("minzoom").unwrap().parse().unwrap();
    let maxzoom: u8 = matches.value_of("maxzoom").unwrap().parse().unwrap();
    let if_not_exists = matches.is_present("if_not_exists");
    let compress = ! matches.is_present("no_compress");
    let metatile_scale: u8 = matches.value_of("metatile-scale").unwrap().parse().unwrap();
    let num_threads: usize = matches.value_of("threads").unwrap().parse().unwrap();

    let bbox = BBox::new_from_string(matches.value_of("bbox").expect("bbox not provided")).expect("Invalid bbox");

    generate_all(&data_yml, minzoom, maxzoom, &bbox, &dest_dir, if_not_exists, compress, metatile_scale, num_threads);
}
