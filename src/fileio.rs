use std::{thread, time};
use std::sync::mpsc::Receiver;
use std::io::{Write, BufWriter};
use std::path::PathBuf;
use slippy_map_tiles;
use std::fs;
use std::fs::File;
use std::io::prelude::*;

#[derive(Debug,Eq,PartialEq)]
pub enum FileIOMessage {
    // Time to quit
    Quit,

    SaveTile(slippy_map_tiles::Tile, Vec<u8>),

    // For when we do per layer tiles
    // AppendToTile(slippy_map_tiles::Tile, Vec<u8>),
    // EnsureAllCompressed,
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

pub fn fileio_thread<D: TileDestination>(rx: Receiver<FileIOMessage>, mut dest: D)
{
    let mut should_quit = false;

    for msg in rx.iter() {
        match msg {
            FileIOMessage::Quit => { break; },
            FileIOMessage::SaveTile(tile, bytes) => {
                dest.save_tile(tile, bytes);
            },
        }
    }

    dest.finish();
}
