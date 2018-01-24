use std::{thread, time};
use std::sync::mpsc::Receiver;
use std::io::Write;
use std::cmp::max;

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
        //write!(stdout, "\r[{:>6}s] z{:>2}, {:4} metatile(s) ({:9.4} metatiles/s, {:9.4} tiles/s, {:5} tiles in last s)", start.elapsed().as_secs(), current_zoom, num_metatiles_done, (num_metatiles_done as f64)/duration, (num_tiles_done as f64)/duration, num_this_sec).ok();
        //stdout.flush().ok();

        if should_quit {
            break;
        } else {
            thread::sleep(time::Duration::from_secs(1));
        }
    }

    let duration = duration_to_float_secs(&start.elapsed());
    println!("\nFinished. {:4} metatile(s), done in {:>8} ({:9.4} metatiles/sec, {:9.4} tiles/sec)", num_metatiles_done, fmt_duration(&start.elapsed()), (num_metatiles_done as f64)/duration, (num_tiles_done as f64)/duration );
}
