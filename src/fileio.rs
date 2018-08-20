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
use byteorder::{LittleEndian, WriteBytesExt};
use serde_json;

#[derive(Debug,Eq,PartialEq)]
pub enum FileIOMessage {
    // Time to quit
    Quit,

    SaveTile(slippy_map_tiles::Tile, Vec<u8>),

    SaveMetaTile(slippy_map_tiles::Metatile, Vec<(slippy_map_tiles::Tile, Vec<u8>)>),

    AppendToTile(slippy_map_tiles::Tile, Vec<u8>),
    // EnsureAllCompressed,
}

pub fn fileio_thread<D: TileDestination+Sized>(rx: Receiver<FileIOMessage>, mut dest: Box<D>)
{
    for msg in rx.iter() {
        match msg {
            FileIOMessage::Quit => { break; },
            FileIOMessage::SaveTile(tile, bytes) => {
                //println!("{}:{} got savetile {:?} bytes.len {}", file!(), line!(), tile, bytes.len());
                dest.save_tile(tile, bytes);
            },
            FileIOMessage::SaveMetaTile(metatile, tiles) => {
                dest.save_metatile(metatile, tiles);
                memory!("Wrote a metatile");
            },
            FileIOMessage::AppendToTile(tile, bytes) => {
                dest.append_bytes_to_tile(tile, bytes);
            },
        }
    }

    dest.finish();
}

pub trait TileDestination {
    fn new(dest_dir: &PathBuf) -> Self;

    fn save_tile(&mut self, tile: slippy_map_tiles::Tile, bytes: Vec<u8>);

    fn save_metatile(&mut self, metatile: slippy_map_tiles::Metatile, tiles: Vec<(slippy_map_tiles::Tile, Vec<u8>)>) {
        for (tile, bytes) in tiles.into_iter() {
            self.save_tile(tile, bytes);
        }
    }

    fn finish(&mut self) {}

    fn does_metatile_exist(dest: &PathBuf, metatile: &slippy_map_tiles::Metatile) -> bool {
        metatile.tiles().iter().all(|t| Self::does_tile_exist(dest, t))
    }

    fn does_tile_exist(dest: &PathBuf, tile: &slippy_map_tiles::Tile) -> bool;

    fn append_bytes_to_tile(&mut self, tile: slippy_map_tiles::Tile, bytes: Vec<u8>) {
        unimplemented!();
    }
}

pub struct TileStashDirectory {
    dest_dir: PathBuf,
}

impl TileDestination for TileStashDirectory {
    fn new(dest_dir: &PathBuf) -> Self {
        fs::create_dir_all(&dest_dir).unwrap();
        TileStashDirectory{ dest_dir: dest_dir.clone() }
    }

    fn save_tile(&mut self, tile: slippy_map_tiles::Tile, bytes: Vec<u8>) {
        let filename = self.dest_dir.join(tile.ts_path("pbf"));
        fs::create_dir_all(filename.parent().unwrap()).unwrap();

        let mut file = BufWriter::new(File::create(filename).unwrap());
        file.write_all(&bytes).unwrap();
    }

    fn does_tile_exist(dest: &PathBuf, tile: &slippy_map_tiles::Tile) -> bool {
        dest.join(tile.ts_path("pbf")).exists()
    }
}

pub struct MBTiles {
    conn: rusqlite::Connection,
}

impl TileDestination for MBTiles {
    fn new(filename: &PathBuf) -> Self {
        let path = filename.clone();
        fs::create_dir_all(&path.parent().unwrap()).unwrap();
        let file_exists = path.is_file();

        let conn = rusqlite::Connection::open(path).unwrap();

        if ! file_exists {
            let creation_schema = include_str!("mbtiles-schema.sql");
            conn.execute_batch(creation_schema).unwrap();
        }

        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('name', 'Vector Tiles');", &[]).unwrap();
        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('description', 'Vector Tiles');", &[]).unwrap();
        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('version', '1');", &[]).unwrap();
        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('format', 'pbf');", &[]).unwrap();

        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('minzoom', '0');", &[]).unwrap();
        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('maxzoom', '14');", &[]).unwrap();
        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('center', '0.0,0.0,0');", &[]).unwrap();
        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('bounds', '-180.0,-85,180,85');", &[]).unwrap();
        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('type', 'basemap');", &[]).unwrap();
        conn.execute( "INSERT OR REPLACE INTO metadata (name, value) VALUES ('scheme', 'tms');", &[]).unwrap();

        // sqlite is *much* faster if everything happens in one transaction, rather than each
        // statement each in it's own transaction. But we can't easily carry a
        // rusqlite::Transaction around due to life times. This is a hack. It will open a
        // transaction (in this connection), at the start
        conn.execute("BEGIN TRANSACTION;", &[]).unwrap();

        MBTiles{ conn: conn }
    }


    fn save_tile(&mut self, tile: slippy_map_tiles::Tile, bytes: Vec<u8>) {
        let digest = format!("{:x}", md5::compute(&bytes));

        let row: u32 = 2u32.pow(tile.zoom() as u32) - tile.y() - 1;

        self.conn.execute(
            "INSERT OR REPLACE INTO map (zoom_level, tile_column, tile_row, tile_id) VALUES (?1, ?2, ?3, ?4);",
            &[&tile.zoom(), &tile.x(), &row, &digest]
            ).unwrap();

        self.conn.execute(
            "INSERT OR IGNORE INTO images (tile_id, tile_data) VALUES (?1, ?2);",
            &[&digest, &bytes]
            ).unwrap();

        // TODO COMMIT?

    }

    fn finish(&mut self) {
        //self.txn().commit();
        self.conn.execute("COMMIT;", &[]).unwrap();
    }

    fn does_tile_exist(dest: &PathBuf, tile: &slippy_map_tiles::Tile) -> bool {
        // Not implemented yet
        false
    }

    fn append_bytes_to_tile(&mut self, tile: slippy_map_tiles::Tile, bytes: Vec<u8>) {
        let digest = format!("{}/{}/{}", tile.zoom(), tile.x(), tile.y());

        let row: u32 = 2u32.pow(tile.zoom() as u32) - tile.y() - 1;

        let num_changed = self.conn.execute(
            "INSERT OR IGNORE INTO map (zoom_level, tile_column, tile_row, tile_id) VALUES (?1, ?2, ?3, ?4);",
            &[&tile.zoom(), &tile.x(), &row, &digest]
            ).unwrap();

        if num_changed == 1 {
            self.conn.execute(
                "INSERT OR IGNORE INTO images (tile_id, tile_data) VALUES (?1, ?2);",
                    &[&digest, &bytes]
                ).unwrap();
        } else {
            self.conn.execute(
                "UPDATE images SET tile_data = tile_data||?1 WHERE tile_id = ?2;",
                    &[&bytes, &digest]
                ).unwrap();
        }
    }

}

impl MBTiles {
    pub fn set_tilejson_vector_layers(&mut self, vector_layers: serde_json::Value) {
        let vector_layers_string = vector_layers.to_string();
        self.conn.execute(
            "INSERT OR REPLACE INTO metadata (name, value) VALUES ('json', ?1);",
            &[&vector_layers_string]
            ).unwrap();
    }
}

pub struct ModTileMetatileDirectory {
    dest_dir: PathBuf,
}

impl TileDestination for ModTileMetatileDirectory {
    fn new(dest_dir: &PathBuf) -> Self {
        fs::create_dir_all(&dest_dir).unwrap();
        ModTileMetatileDirectory{ dest_dir: dest_dir.clone() }
    }

    fn does_tile_exist(dest: &PathBuf, tile: &slippy_map_tiles::Tile) -> bool {
        // Not implemented yet
        false
    }

    fn save_metatile(&mut self, metatile: slippy_map_tiles::Metatile, tiles: Vec<(slippy_map_tiles::Tile, Vec<u8>)>) {
        let size = metatile.size() as usize;
        let x = metatile.x();
        let y = metatile.y();
        let zoom = metatile.zoom();

        if size > 8 {
            // Split this metatile into smaller 8x8 metatiles.
            // We still require a whole number multiple of 8
            assert_eq!(size % 8, 0);
            let multiple = (size / 8) as usize;
            let mut new_metatiles = vec![Vec::new(); multiple*multiple];
            for (tile, bytes) in tiles.into_iter() {
                let new_mt_x = ((tile.x() - x) / 8) as usize;
                let new_mt_y = ((tile.y() - y) / 8) as usize;
                new_metatiles[new_mt_x*multiple+new_mt_y].push((tile, bytes));
            }

            for (offset, tiles) in new_metatiles.into_iter().enumerate() {
                let offset = offset as u32;
                self.save_metatile(slippy_map_tiles::Metatile::new(8, zoom, x+(offset/8), y+(offset%8)).unwrap(), tiles);
            }

            return;
        }

        // TODO suspect I can optimize this...
        // Are there unnecessary memory copies?
        assert!(size <= 8);
        let mut tiles_array: Vec<Vec<u8>> = vec![Vec::new(); size*size];
        for (tile, bytes) in tiles.into_iter() {
            let i = ((tile.x()-x)*size as u32 + (tile.y()-y)) as usize;
            tiles_array[i] = bytes;
        }
        let tiles = tiles_array;

        let mut offsets = vec![0; size*size];
        let mut sizes = vec![0; size*size];

        // offset is from the start of the file, so it starts with
        // 'META' = 4 bytes. count = 4 bytes. 12 bytes for x, y, z
        // 8 * count for offsets
        let mut curr_offset = 4 + 4 + 12 + 8*(size*size);
        for i in 0..(size*size) {
            offsets[i] = curr_offset;
            let this_size = tiles[i].len();
            sizes[i] = this_size;
            curr_offset += this_size;
        }


        let filename = self.dest_dir.join(xyz_to_mt(metatile.zoom(), metatile.x(), metatile.y(), "meta"));
        fs::create_dir_all(filename.parent().unwrap()).unwrap();

        let mut file = BufWriter::new(File::create(filename).unwrap());
        file.write_all(&[0x4d, 0x45, 0x54, 0x41]).unwrap(); // 'META' magic string
        file.write_u32::<LittleEndian>((size*size) as u32).unwrap();
        file.write_u32::<LittleEndian>(x).unwrap();
        file.write_u32::<LittleEndian>(y).unwrap();
        file.write_u32::<LittleEndian>(metatile.zoom() as u32).unwrap();

        for i in 0..(size*size) {
            file.write_u32::<LittleEndian>(offsets[i] as u32).unwrap();
            file.write_u32::<LittleEndian>(sizes[i] as u32).unwrap();
        }

        for i in 0..(size*size) {
            file.write_all(&tiles[i]).unwrap();
        }

    }

    fn save_tile(&mut self, tile: slippy_map_tiles::Tile, bytes: Vec<u8>) {
        panic!("Use save_metatile instead");
    }

    fn does_metatile_exist(dest: &PathBuf, metatile: &slippy_map_tiles::Metatile) -> bool {
        dest.join(xyz_to_mt(metatile.zoom(), metatile.x(), metatile.y(), "meta")).exists()
        
    }
}

/// Convert x & y to a ModTile metatile directory parts
fn xyz_to_mt(z: u8, x: u32, y: u32, ext: &str) -> String {
    // /[Z]/[xxxxyyyy]/[xxxxyyyy]/[xxxxyyyy]/[xxxxyyyy]/[xxxxyyyy].png
    // i.e. /[Z]/a/b/c/d/e.png

    let mut x = x;
    let mut y = y;

    let e = (((x & 0x0f) << 4) | (y & 0x0f)) as u8;
    x >>= 4;
    y >>= 4;

    let d = (((x & 0x0f) << 4) | (y & 0x0f)) as u8;
    x >>= 4;
    y >>= 4;

    let c = (((x & 0x0f) << 4) | (y & 0x0f)) as u8;
    x >>= 4;
    y >>= 4;

    let b = (((x & 0x0f) << 4) | (y & 0x0f)) as u8;
    x >>= 4;
    y >>= 4;

    let a = (((x & 0x0f) << 4) | (y & 0x0f)) as u8;
    //x >>= 4;
    //y >>= 4;

    format!("{}/{}/{}/{}/{}/{}.{}", z, a, b, c, d, e, ext)
}

