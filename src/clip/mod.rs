use std::collections::{HashSet, HashMap};
use std::fmt::Debug;
use std::borrow::{Cow, Borrow};

use yaml_rust::{YamlLoader, Yaml};

use slippy_map_tiles;
use slippy_map_tiles::{BBox, Metatile, MetatilesIterator};

use geo::*;
use geo::algorithm::simplify::Simplify;
use geo::algorithm::map_coords::MapCoords;
use geo::algorithm::map_coords::MapCoordsInplace;
use geo::algorithm::boundingbox::BoundingBox;
use geo::algorithm::contains::Contains;

use ::simplify::remove_duplicate_points;

// local stuff
mod cohen_sutherland;
mod sutherland_hodgeman;
//mod weiler_atherton;

#[cfg(test)]
mod test;

// A border we want
#[derive(Debug,Clone,Copy)]
pub enum Border<T: CoordinateType> {
    XMin(T),
    XMax(T),
    YMin(T),
    YMax(T),
}

fn is_inside<T: CoordinateType>(p: &Point<T>, border: &Border<T>) -> bool {
    match *border {
        Border::XMin(xmin) => p.x() >= xmin,
        Border::XMax(xmax) => p.x() <= xmax,
        Border::YMin(ymin) => p.y() >= ymin,
        Border::YMax(ymax) => p.y() <= ymax,
    }
}

fn is_on_border<T: CoordinateType>(p: &Point<T>, border: &Border<T>) -> bool {
    match *border {
        Border::XMin(xmin) => p.x() == xmin,
        Border::XMax(xmax) => p.x() == xmax,
        Border::YMin(ymin) => p.y() == ymin,
        Border::YMax(ymax) => p.y() == ymax,
    }
}


#[derive(Clone, Debug, PartialEq, PartialOrd)]
enum LineBorderIntersection<T: CoordinateType> {

    /// All points are inside the border
    AllInside,

    /// All points are outside the border
    AllOutside,

    /// Some are inside, some outside
    Intersections(Vec<IntersectionOption<T>>),
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
enum IntersectionOption<T: CoordinateType> {
    /// This point is inside the border
    Inside,

    /// This point is outside the border
    Outside,

    /// This point is outside the border and the next point is inside, and this (x, y) is where it
    /// crosses the border.
    Entry((T, T)),

    /// This point is inside the border and the next point is outside, and this (x, y) is where it
    /// crosses the border.
    Exit((T, T)),
}


/// given 2 points, and a border, return the (x, y) where the line between the two points crosses
/// the border. This assumes that both points are on different sizes of the border
fn intersection<T: CoordinateType>(p1: &Point<T>, p2: &Point<T>, border: &Border<T>) -> (T, T) {
    let x1 = p1.x();
    let y1 = p1.y();
    let x2 = p2.x();
    let y2 = p2.y();

    // The algorithms are rearranged to reduce the risk of an i32 integer overflow
    match *border {
        Border::XMin(xmin) => {
            let x = xmin;
            assert!(x1 != x2);

            let a = y2 - y1;
            let b = x - x1;
            let c = x2 - x1;
            let mut delta_y = T::zero();
            delta_y = (a/c)*b;
            if delta_y == T::zero() {
                delta_y = (b/c)*a;
                if delta_y == T::zero() {
                    // might be overflow, but we presume not.
                    delta_y = (a*b)/c;
                }
            }
            let y = y1 + delta_y;

            (x, y)
        },
        Border::XMax(xmax) => {
            let x = xmax;
            assert!(x1 != x2);

            let a = y2 - y1;
            let b = x - x1;
            let c = x2 - x1;
            let mut delta_y = T::zero();
            delta_y = (a/c)*b;
            if delta_y == T::zero() {
                delta_y = (b/c)*a;
                if delta_y == T::zero() {
                    // might be overflow, but we presume not.
                    delta_y = (a*b)/c;
                }
            }
            let y = y1 + delta_y;

            (x, y)
        },
        Border::YMin(ymin) => {
            let y = ymin;
            assert!(y1 != y2);

            let a = x2 - x1;
            let b = y - y1;
            let c = y2 - y1;
            let mut delta_x = T::zero();
            delta_x = (a/c)*b;
            if delta_x == T::zero() {
                delta_x = (b/c)*a;
                if delta_x == T::zero() {
                    // might be overflow, but we presume not.
                    delta_x = (a*b)/c;
                }
            }
            let x = x1 + delta_x;

            (x, y)
        },
        Border::YMax(ymax) => {
            let y = ymax;
            assert!(y1 != y2);

            let a = x2 - x1;
            let b = y - y1;
            let c = y2 - y1;
            let mut delta_x = T::zero();
            delta_x = (a/c)*b;
            if delta_x == T::zero() {
                delta_x = (b/c)*a;
                if delta_x == T::zero() {
                    // might be overflow, but we presume not.
                    delta_x = (a*b)/c;
                }
            }
            let x = x1 + delta_x;

            (x, y)
        },
    }
}


fn calculate_intersections<T: CoordinateType>(linestring: &LineString<T>, border: &Border<T>) -> LineBorderIntersection<T> {
    // First we look if everything is all inside or all outside, and early return then, Then we
    // don't have to allocate a vec for the intersections
    let mut all_inside = true;
    let mut not_all_outside = false;

    let point_inside = is_inside(&linestring.0[0], border);
    all_inside &= point_inside;
    not_all_outside |= point_inside;

    for (idx, point) in linestring.0.iter().skip(1).enumerate() {
        let point_inside = is_inside(point, border);

        all_inside &= point_inside;
        not_all_outside |= point_inside;
    }

    let all_outside = ! not_all_outside;
    if all_inside {
        return LineBorderIntersection::AllInside;
    } else if all_outside {
        return LineBorderIntersection::AllOutside;
    }

    assert!(!all_inside && !all_outside);
    drop(all_inside);
    drop(all_outside);

    // otherwise, we need to look at the intersections

    // I'm not happy with this, since you need to iterate over the linestring, twice and I wish
    // there was a way to do it only once.

    #[allow(unused_assignments)]
    let mut last_inside: bool = false;

    let mut intersections: Vec<IntersectionOption<T>> = Vec::with_capacity(linestring.0.len() - 1);

    // first point
    let point_inside = is_inside(&linestring.0[0], border);
    last_inside = point_inside;
    let mut last_point = &linestring.0[0];

    for (idx, point) in linestring.0.iter().skip(1).enumerate() {
        let point_inside = is_inside(point, border);

        let intersection_option = if point_inside {
            if last_inside {
                IntersectionOption::Inside
            } else {
                IntersectionOption::Entry(intersection(last_point, point, border))
            }
        } else {
            if last_inside {
                if is_on_border(last_point, border) {
                    IntersectionOption::Inside
                } else {
                    IntersectionOption::Exit(intersection(last_point, point, border))
                }
            } else {
                IntersectionOption::Outside
            }
        };

        // This is a little confusing. we want intersections[n] to refer to the segment from point
        // n to point n+1. But in the loop, the point is actually the point n+1, but since we
        // didn't push anything when doing the first point, it shouldn't be a problem.
        intersections.push(intersection_option);

        last_inside = point_inside;
        last_point = point;
    }

    // just for the last point
    intersections.push(if last_inside { IntersectionOption::Inside } else { IntersectionOption::Outside });

    LineBorderIntersection::Intersections(intersections)
}


fn clip_point_to_border<T: CoordinateType>(p: Cow<Point<T>>, border: &Border<T>) -> Option<Geometry<T>> {
    if is_inside(&p, border) {
        Some(Geometry::Point(p.into_owned()))
    } else {
        None
    }
}

fn clip_linestring_to_border<T: CoordinateType>(ls: Cow<LineString<T>>, border: &Border<T>) -> Option<Geometry<T>> {
    match calculate_intersections(&ls, border) {
        LineBorderIntersection::AllInside => Some(Geometry::LineString(ls.into_owned())),
        LineBorderIntersection::AllOutside => None,
        LineBorderIntersection::Intersections(intersections) => {
            let mut lines: Vec<LineString<T>> = Vec::new();
            let mut curr_line: Vec<Point<T>> = Vec::new();

            for (point, intersection) in ls.0.iter().zip(intersections.into_iter()) {
                match intersection {
                    IntersectionOption::Inside => {
                        curr_line.push(point.clone());
                    },
                    IntersectionOption::Outside => {
                        if ! curr_line.is_empty() {
                            if curr_line.len() > 1 {
                                lines.push(LineString(curr_line));
                            }
                            curr_line = Vec::new();
                        }
                    },
                    IntersectionOption::Entry((x, y)) => {
                        curr_line.push(Point::new(x, y));
                    },
                    IntersectionOption::Exit((x, y)) => {
                        curr_line.push(point.clone());
                        curr_line.push(Point::new(x, y));
                        if curr_line.len() > 1 {
                            lines.push(LineString(curr_line));
                        }
                        curr_line = Vec::new();
                    }
                }
            }

            if curr_line.len() > 1 {
                lines.push(LineString(curr_line));
            }

            match lines.len() {
                // When there would only be one point returned.
                0 => None,

                1 => Some(Geometry::LineString(lines.remove(0))),
                _ => Some(Geometry::MultiLineString(MultiLineString(lines))),
            }
        }
    }
}


fn clip_multipoint_to_border<T: CoordinateType>(mp: Cow<MultiPoint<T>>, border: &Border<T>) -> Option<Geometry<T>> {
    let points: Vec<Point<T>> = match mp {
        Cow::Owned(mp) => mp.0.into_iter().filter_map(|p| clip_point_to_border(Cow::Owned(p), border).and_then(|g| g.as_point())).collect(),
        Cow::Borrowed(mp) => mp.0.iter().filter_map(|p| clip_point_to_border(Cow::Borrowed(p), border).and_then(|g| g.as_point())).collect(),
    };
    if points.is_empty() {
        None
    } else {
        Some(Geometry::MultiPoint(MultiPoint(points)))
    }
}

fn clip_multilinestring_to_border<T: CoordinateType>(mls: Cow<MultiLineString<T>>, border: &Border<T>) -> Option<Geometry<T>> {
    let mut lines: Vec<LineString<T>> = Vec::with_capacity(mls.0.len());
    match mls {
        Cow::Borrowed(mls) => {
            for ls in mls.0.iter() {
                match clip_linestring_to_border(Cow::Borrowed(ls), border) {
                    None => {},
                    Some(Geometry::LineString(new_ls)) => { lines.push(new_ls); },
                    Some(Geometry::MultiLineString(new_mls)) => { lines.extend(new_mls); }
                    _ => { unreachable!(); },
                }
            }
        },
        Cow::Owned(mls) => {
            for ls in mls.0.into_iter() {
                match clip_linestring_to_border(Cow::Owned(ls), border) {
                    None => {},
                    Some(Geometry::LineString(new_ls)) => { lines.push(new_ls); },
                    Some(Geometry::MultiLineString(new_mls)) => { lines.extend(new_mls); }
                    _ => { unreachable!(); },
                }
            }
        }
    }

    match lines.len() {
        0 => None,
        1 => Some(Geometry::LineString(lines.remove(0))),
        _ => Some(Geometry::MultiLineString(MultiLineString(lines))),
    }
}

/// Clip a geometry to a border. None iff the geometry is entirely outside it, otherwise Some(g)
/// with the new geometry clipped to that border.
fn clip_to_border<T: CoordinateType>(geom: Cow<Geometry<T>>, border: &Border<T>) -> Option<Geometry<T>> {
    match geom {
        Cow::Owned(geom) => {
            match geom {
                Geometry::Point(p) => clip_point_to_border(Cow::Owned(p), border),
                Geometry::LineString(ls) => clip_linestring_to_border(Cow::Owned(ls), border),
                Geometry::Polygon(p) => sutherland_hodgeman::clip_polygon_to_border(Cow::Owned(p), border).map(Geometry::Polygon),
                Geometry::MultiPoint(mp) => clip_multipoint_to_border(Cow::Owned(mp), border),
                Geometry::MultiLineString(mls) => clip_multilinestring_to_border(Cow::Owned(mls), border),
                Geometry::MultiPolygon(p) => sutherland_hodgeman::clip_multipolygon_to_border(Cow::Owned(p), border).map(Geometry::MultiPolygon),
                Geometry::GeometryCollection(_) => unimplemented!(),
            }
        },
        Cow::Borrowed(geom) => {
            match *geom {
                Geometry::Point(ref p) => clip_point_to_border(Cow::Borrowed(p), border),
                Geometry::LineString(ref ls) => clip_linestring_to_border(Cow::Borrowed(ls), border),
                Geometry::Polygon(ref p) => sutherland_hodgeman::clip_polygon_to_border(Cow::Borrowed(p), border).map(Geometry::Polygon),

                Geometry::MultiPoint(ref mp) => clip_multipoint_to_border(Cow::Borrowed(mp), border),
                Geometry::MultiLineString(ref mls) => clip_multilinestring_to_border(Cow::Borrowed(mls), border),
                Geometry::MultiPolygon(ref p) => sutherland_hodgeman::clip_multipolygon_to_border(Cow::Borrowed(p), border).map(Geometry::MultiPolygon),
                Geometry::GeometryCollection(_) => unimplemented!(),
            }
        }
    }
}


pub fn clip_to_bbox<T: CoordinateType+::std::fmt::Debug>(geom: Cow<Geometry<T>>, bbox: &Bbox<T>) -> Option<Geometry<T>> {
    clip_to_border(geom, &Border::XMin(bbox.xmin))
       .and_then(|geom| clip_to_border(Cow::Owned(geom), &Border::XMax(bbox.xmax)))
       .and_then(|geom| clip_to_border(Cow::Owned(geom), &Border::YMin(bbox.ymin)))
       .and_then(|geom| clip_to_border(Cow::Owned(geom), &Border::YMax(bbox.ymax)))
}

/// geom - The geometry
/// metatile_scale - e.g. 8 for an 8x8 tile. This could be a 2x2 metatile, so that'd be 2
/// zoom - of the tile
/// tile_x0/tile_y0 - the x/y value of the top left tile of the metatile
/// x0/y0 - The x/y value of the top left corner of the topleft tile. At the start this will be (0, 0)
/// size - the width (& height) of the metatile. A regular tile is 4096. So a 2x2 is 8192 etc
fn slice_box(geom: Cow<Geometry<i32>>, metatile_scale: u8, zoom: u8, tile_x0: u32, tile_y0: u32, x0: i32, y0: i32, size: i32, buffer: i32) -> Vec<(slippy_map_tiles::Tile, Option<Geometry<i32>>)> {
    if metatile_scale == 1 {
        return vec![(slippy_map_tiles::Tile::new(zoom, tile_x0, tile_y0).unwrap(), Some(geom.into_owned()))];
    }

    let mut results = Vec::new();

    let half = size / 2;
    let tile_half = (metatile_scale/2) as u32;

    if let Some(left) = clip_to_border(Cow::Borrowed(geom.borrow()), &Border::XMax(x0+half+buffer as i32)) {
        if let Some(topleft) = clip_to_border(Cow::Borrowed(&left), &Border::YMax(y0+half+buffer as i32)) {
            let mut tiles = slice_box(Cow::Owned(topleft), metatile_scale/2, zoom, tile_x0, tile_y0, x0, y0, size/2, buffer);
            if metatile_scale == 2 {
                // this is only one tile
                match tiles.len() {
                    0 => {},
                    1 => { results.push(tiles.remove(0)); },
                    _ => { unreachable!(); },
                }
            } else {
                results.append(&mut tiles);
            }
        }

        if let Some(bottomleft) = clip_to_border(Cow::Owned(left), &Border::YMin(y0+half-buffer as i32)) {
            let mut tiles = slice_box(Cow::Owned(bottomleft), metatile_scale/2, zoom, tile_x0, tile_y0+tile_half, x0, y0+half, size/2, buffer);
            if metatile_scale == 2 {
                // this is only one tile
                match tiles.len() {
                    0 => {},
                    1 => { results.push(tiles.remove(0)); },
                    _ => { unreachable!(); },
                }
            } else {
                results.append(&mut tiles);
            }
        }
    }

    if let Some(right) = clip_to_border(geom, &Border::XMin(x0+half-buffer as i32)) {
        if let Some(topright) = clip_to_border(Cow::Borrowed(&right), &Border::YMax(y0+half+buffer as i32)) {
            let mut tiles = slice_box(Cow::Owned(topright), metatile_scale/2, zoom, tile_x0+tile_half, tile_y0, x0+half, y0, size/2, buffer);
            if metatile_scale == 2 {
                // this is only one tile
                match tiles.len() {
                    0 => {},
                    1 => { results.push(tiles.remove(0)); },
                    _ => { unreachable!(); },
                }
            } else {
                results.append(&mut tiles);
            }
        }

        if let Some(bottomright) = clip_to_border(Cow::Owned(right), &Border::YMin(y0+half-buffer as i32)) {
            let mut tiles = slice_box(Cow::Owned(bottomright), metatile_scale/2, zoom, tile_x0+tile_half, tile_y0+tile_half, x0+half, y0+half, size/2, buffer);
            if metatile_scale == 2 {
                // this is only one tile
                match tiles.len() {
                    0 => {},
                    1 => { results.push(tiles.remove(0)); },
                    _ => { unreachable!(); },
                }
            } else {
                results.append(&mut tiles);
            }
        }
    }

    results
}

fn all_points_in_one_tile(metatile: &Metatile, ls: &LineString<i32>, buffer: i32) -> Option<(u32, u32)> {
    assert_eq!(buffer, 0);
    let buffer = buffer as i32;

    let initial_tile_x = ls.0[0].x() / 4096;
    let initial_tile_y = ls.0[0].y() / 4096;

    if initial_tile_x < buffer && initial_tile_x > 0 {
        // point is in buffer to the left, so is in this tile and one to left
        None
    } else if initial_tile_x > (4096 - buffer) && initial_tile_x < 4096 {
        // point is in buffer to the right, so is in this tile and one to right
        None
    } else if initial_tile_y < buffer && initial_tile_y > 0 {
        // point is in buffer to the top, so is in this tile and one aboce
        None
    } else if initial_tile_y > (4096 - buffer) && initial_tile_y < 4096 {
        // point is in buffer to the bottom, so is in this tile and one below
        None
    } else if ls.0.iter().skip(1).all(|&p| (p.x() / 4096 == initial_tile_x && p.y() / 4096 == initial_tile_y)) {
        // FIXME finish this
        Some((initial_tile_x as u32, initial_tile_y as u32))
    } else {
        None
    }
}

pub fn clip_linestring_to_tiles(metatile: &Metatile, mut ls: LineString<i32>, buffer: i32) -> Vec<(slippy_map_tiles::Tile, Option<Geometry<i32>>)> {
    assert_eq!(buffer, 0);
    // Are all the points in the linestring in the same tile? If so, we can skip a lot of steps
    // TODO this is not buffer aware.
    let size = metatile.size() as u32;

    let single_tile = all_points_in_one_tile(metatile, &ls, buffer);

    match single_tile {
        None => {
            // Not all on one tile, so usual approach
            slice_box(Cow::Owned(Geometry::LineString(ls)), metatile.size(), metatile.zoom(), metatile.x(), metatile.y(), 0, 0, metatile.size() as i32*4096, buffer)
        },
        Some((tile_x, tile_y)) => {
            if tile_x < size && tile_y < size {
                vec![(slippy_map_tiles::Tile::new(metatile.zoom(), (tile_x as u32)+metatile.x(), (tile_y as u32)+metatile.y()).unwrap(), Some(Geometry::LineString(ls)))]
            } else {
                // Dunno why
                // TODO why is this needed
                slice_box(Cow::Owned(Geometry::LineString(ls)), metatile.size(), metatile.zoom(), metatile.x(), metatile.y(), 0, 0, metatile.size() as i32*4096, buffer)
            }
        },
    }
    
}

pub fn clip_point_to_tiles(metatile: &Metatile, point: Point<i32>, buffer: i32) -> Vec<(slippy_map_tiles::Tile, Option<Geometry<i32>>)> {
    // TODO support buffer
    assert_eq!(buffer, 0);
    let metatile_scale = metatile.size() as u32;

    let x = point.x();
    let y = point.x();

    let tile_x = (x / 4096) as u32 + metatile.x();
    let tile_y = (y / 4096) as u32 + metatile.y();

    let tile_x_offset = x % 4096;
    let tile_y_offset = y % 4096;


    // TODO fill in
    if tile_x > 0 && tile_y > 0 && tile_x < metatile_scale && tile_y < metatile_scale {
        vec![(slippy_map_tiles::Tile::new(metatile.zoom(), tile_x, tile_y).unwrap(), Some(Geometry::Point(Point::new(tile_x_offset, tile_y_offset))))]
    } else {
        vec![]
    }
}

pub fn clip_geometry_to_tiles(metatile: &Metatile, geom: Geometry<i32>, buffer: i32) -> Vec<(slippy_map_tiles::Tile, Option<Geometry<i32>>)> {
    // TODO somehow in this method, it's making invalid polygons where the last point != first
    // point
    // Simple approach for now
    let mut res = slice_box(Cow::Owned(geom), metatile.size(), metatile.zoom(), metatile.x(), metatile.y(), 0, 0, metatile.size() as i32*4096, buffer);

    // TODO make the slice_box etc not produce geoms with this result
    for &mut (tile, ref mut geom_opt) in res.iter_mut() {
        if let &mut Some(ref mut geom) = geom_opt {
            // TODO replace with remove_unneeded_points ?
            remove_duplicate_points(geom);
        }
    }

    res


    // re-enable this later, when all functions support buffer
    //match geom {
    //    Geometry::Point(p) => clip_point_to_tiles(metatile, p, buffer),
    //    Geometry::LineString(ls) => clip_linestring_to_tiles(metatile, ls, buffer),
    //    _ => slice_box(Cow::Owned(geom), metatile.size(), metatile.zoom(), metatile.x(), metatile.y(), 0, 0, metatile.size() as i32*4096, buffer),
    //}
}
