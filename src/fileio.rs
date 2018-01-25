use std::{thread, time};
use std::sync::mpsc::Receiver;
use std::io::{Write, BufWriter};
use std::path::PathBuf;
use slippy_map_tiles;
use std::fs;
use std::fs::File;
use std::io::prelude::*;
use rusqlite;
use md5;

#[derive(Debug,Eq,PartialEq)]
pub enum FileIOMessage {
    // Time to quit
    Quit,

    SaveTile(slippy_map_tiles::Tile, Vec<u8>),

    // For when we do per layer tiles
    // AppendToTile(slippy_map_tiles::Tile, Vec<u8>),
    // EnsureAllCompressed,
}

pub fn fileio_thread<D: TileDestination+Sized>(rx: Receiver<FileIOMessage>, mut dest: Box<D>)
{
    for msg in rx.iter() {
        //println!("{}:{} msg {:?}", file!(), line!(), msg);
        match msg {
            FileIOMessage::Quit => { break; },
            FileIOMessage::SaveTile(tile, bytes) => {
                dest.save_tile(tile, bytes);
            },
        }
    }

    dest.finish();
}

pub trait TileDestination {
    fn save_tile(&mut self, tile: slippy_map_tiles::Tile, bytes: Vec<u8>);
    fn finish(&mut self) {}
}

pub struct TileStashDirectory {
    dest_dir: PathBuf,
}

impl TileStashDirectory {
    pub fn new(dest_dir: &PathBuf) -> Self {
        fs::create_dir_all(&dest_dir).unwrap();
        TileStashDirectory{ dest_dir: dest_dir.clone() }
    }
}

impl TileDestination for TileStashDirectory {
    fn save_tile(&mut self, tile: slippy_map_tiles::Tile, bytes: Vec<u8>) {
        let filename = self.dest_dir.join(tile.ts_path("pbf"));
        fs::create_dir_all(filename.parent().unwrap()).unwrap();

        let mut file = BufWriter::new(File::create(filename).unwrap());
        file.write_all(&bytes).unwrap();
    }
}

pub struct MBTiles {
    conn: rusqlite::Connection,
}

impl MBTiles {
    pub fn new(filename: &PathBuf) -> Self {
        let path = filename.clone();
        fs::create_dir_all(&path.parent().unwrap()).unwrap();
        if path.is_file() {
            panic!("mbtiles filename already exists");
        }
        let conn = rusqlite::Connection::open(path).unwrap();
        let creation_schema = include_str!("mbtiles-schema.sql");
        conn.execute_batch(creation_schema).unwrap();
        MBTiles{ conn: conn }
    }
}

impl TileDestination for MBTiles {
    fn save_tile(&mut self, tile: slippy_map_tiles::Tile, bytes: Vec<u8>) {
        //let digest = format!("{:x}", md5::compute(&bytes));
        let digest = format!("{}/{}/{}", tile.zoom(), tile.x(), tile.y());
        self.conn.execute(
            "INSERT INTO map (zoom_level, tile_column, tile_row, tile_id) VALUES (?1, ?2, ?3, ?4);",
            &[&tile.zoom(), &tile.x(), &tile.y(), &digest]
            ).unwrap();

        self.conn.execute(
            "INSERT OR IGNORE INTO images (tile_id, tile_data) VALUES (?1, ?2);",
            &[&digest, &bytes]
            ).unwrap();

        // TODO COMMIT?

    }
}
