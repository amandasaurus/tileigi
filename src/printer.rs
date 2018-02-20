use std::{thread, time};
use std::sync::mpsc::Receiver;
use std::io::Write;
use std::cmp::max;

use slippy_map_tiles;
use separator::Separatable;

#[derive(Debug,Eq,PartialEq)]
pub enum PrinterMessage {
    // Time to quit
    Quit,

    // has done this many tiles, (zoom, num_metatiles, num_tiles)
    DoneTiles(u8, usize, usize),
}

fn fmt_rate(num: f64) -> String {
    if num > 1000. {
        (num.trunc() as i64).separated_string()
    } else {
        format!("{:.4}", num)
    }
}

fn fmt_duration(sec: u64) -> String {
    let (min, sec) = (sec / 60, sec % 60);
    let (hr, min) = (min / 60, min % 60);
    let (day, hr) = (hr / 24, hr % 24);
    if day > 0 {
        format!("{}d{:>02}h{:>02}m{:>02}s", day, hr, min, sec)
    } else if hr > 0 {
        format!("{}h{:>02}m{:>02}s", hr, min, sec)
    } else if min > 0 {
        format!("{}m{:>02}s", min, sec)
    } else {
        format!("{}s", sec)
    }
}

fn duration_to_float_secs(dur: &time::Duration) -> f64 {
    (dur.as_secs() as f64) + (dur.subsec_nanos() as f64 / 1e9)
}

fn round(x: f64, places: i32) -> f64 {
    let mult: f64 = 10_f64.powi(places);
    (x * mult).round()/mult
}

pub fn printer(rx: Receiver<PrinterMessage>, bbox: Option<slippy_map_tiles::BBox>, metatile_scale: u8, minzoom: u8, maxzoom: u8) {
    let mut num_metatiles_done = 0;
    let mut num_tiles_done = 0;
    let mut current_zoom = 0;
    let mut should_quit = false;

    let start = time::Instant::now();

    let mut stdout = ::std::io::stdout();
    
    let total_num_of_metatiles: Option<usize> = (minzoom..maxzoom+1).map(|z| {
        match bbox {
            None => {
                let scale = (metatile_scale.trailing_zeros()) as u32;
                Some(2_u64.pow(2 * (z as u32-scale) as u32) as usize)
            },
            Some(ref bbox) => {
                let this = slippy_map_tiles::size_bbox_zoom_metatiles(&bbox, z, metatile_scale);
                //println!("z {} {:?}", z, this);
                this
            },
        }})
        .fold(Some(0_usize), |acc, on_this_zoom| {
            match (acc, on_this_zoom) {
                (Some(x), Some(y)) => x.checked_add(y),
                _ => None,
            }
        });

        //println!("total_num_of_metatiles {:?}", total_num_of_metatiles);
                   


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
        let percent_done = match total_num_of_metatiles {
            None => "% N/A".to_string(),
            Some(total_num_of_metatiles) => format!("{:>.3}%", (num_metatiles_done as f32*100.)/(total_num_of_metatiles as f32))
        };

        let eta = match total_num_of_metatiles {
            None => "N/A".to_string(),
            Some(total_num_of_metatiles) => {
                let fraction = (num_metatiles_done as f32)/(total_num_of_metatiles as f32);
                let est_total = (start.elapsed().as_secs() as f32/fraction).round() as u64;
                if start.elapsed().as_secs() < 3 || start.elapsed().as_secs() > est_total {
                    // Can happen at the start with rounding and few samples
                    format!("N/A")
                } else {
                    let eta_secs = est_total - start.elapsed().as_secs();
                    fmt_duration(eta_secs)
                }
            }
        };

        let duration = duration_to_float_secs(&start.elapsed());
        write!(stdout, "\r{percent} [{duration:>6}] ETA: {eta} z{zoom:>2}, {num_tiles:>4} tiles ({tiles_p_s:>9} tiles/s, last sec: {last_num:>5} tiles)   ",
            duration=fmt_duration(start.elapsed().as_secs()),
            eta=eta,
            zoom=current_zoom,
            num_tiles=num_tiles_done.separated_string(),
            percent=percent_done,
            //num_metatiles=num_metatiles_done.separated_string(),
            tiles_p_s=fmt_rate((num_tiles_done as f64)/duration),
            //mt_p_s=round((num_metatiles_done as f64)/duration, 4).separated_string(),
            last_num=num_this_sec.separated_string()
        ).ok();

        stdout.flush().ok();

        if should_quit {
            break;
        } else {
            thread::sleep(time::Duration::from_secs(1));
        }
    }

    let duration = duration_to_float_secs(&start.elapsed());
    println!("\nFinished. {} tiles ({} metatiles), done in {} ({:>9} metatiles/sec, {:>9} tiles/sec)", num_tiles_done.separated_string(), num_metatiles_done.separated_string(), fmt_duration(start.elapsed().as_secs()), fmt_rate((num_metatiles_done as f64)/duration), fmt_rate((num_tiles_done as f64)/duration) );
}
