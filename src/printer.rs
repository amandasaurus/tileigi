use std::{thread, time};
use std::sync::mpsc::Receiver;
use std::io::Write;
use std::cmp::max;

use separator::Separatable;

#[derive(Debug,Eq,PartialEq)]
pub enum PrinterMessage {
    // Time to quit
    Quit,

    // has done this many tiles, (zoom, num_metatiles, num_tiles)
    DoneTiles(u8, usize, usize),
}

fn fmt_duration(dur: &time::Duration) -> String {
    format!("{:.2}s", duration_to_float_secs(dur))
}

fn duration_to_float_secs(dur: &time::Duration) -> f64 {
    (dur.as_secs() as f64) + (dur.subsec_nanos() as f64 / 1e9)
}

fn round(x: f64, places: i32) -> f64 {
    let mult: f64 = 10_f64.powi(places);
    (x * mult).round()/mult
}

pub fn printer(rx: Receiver<PrinterMessage>) {
    let mut num_metatiles_done = 0;
    let mut num_tiles_done = 0;
    let mut current_zoom = 0;
    let mut should_quit = false;

    let start = time::Instant::now();

    let mut stdout = ::std::io::stdout();

    loop {
        let mut num_this_sec = 0;

        // collect all the pending messages
        for msg in rx.try_iter() {
            match msg {
                PrinterMessage::Quit => {
                    should_quit = true;
                },
                PrinterMessage::DoneTiles(zoom, num_meta, num_tiles) => {
                    // If we finish some z8s and then a z7 is finished, we don't want the
                    // current_zoom to go bck to 7, so keep it at the max we've seen.
                    current_zoom = max(current_zoom, zoom);

                    num_metatiles_done += num_meta;
                    num_tiles_done += num_tiles;
                    num_this_sec += num_tiles;
                }
            }
        }

        let duration = duration_to_float_secs(&start.elapsed());
        write!(stdout, "\r[{:>6}s] z{:>2}, {:>4} tiles ({:>9} tiles/s, {:>9} metatiles/s, last sec: {:>5} tiles)   ", start.elapsed().as_secs(), current_zoom, num_tiles_done.separated_string(), round((num_tiles_done as f64)/duration, 4).separated_string(), round((num_metatiles_done as f64)/duration, 4).separated_string(), num_this_sec.separated_string()).ok();
        stdout.flush().ok();

        if should_quit {
            break;
        } else {
            thread::sleep(time::Duration::from_secs(1));
        }
    }

    let duration = duration_to_float_secs(&start.elapsed());
    println!("\nFinished. {} tiles ({} metatiles), done in {} ({:>9} metatiles/sec, {:>9} tiles/sec)", num_tiles_done.separated_string(), num_metatiles_done.separated_string(), fmt_duration(&start.elapsed()), round((num_metatiles_done as f64)/duration, 4).separated_string(), round((num_tiles_done as f64)/duration, 4).separated_string() );
}
