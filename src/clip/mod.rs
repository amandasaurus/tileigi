use std::collections::{HashSet, HashMap};
use std::fmt::Debug;

use yaml_rust::{YamlLoader, Yaml};

use slippy_map_tiles;
use slippy_map_tiles::{BBox, Metatile, MetatilesIterator};

use geo::*;
use geo::algorithm::simplify::Simplify;
use geo::algorithm::map_coords::MapCoords;
use geo::algorithm::boundingbox::BoundingBox;
use geo::algorithm::contains::Contains;

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
            let y = ((y2-y1)/(x2 - x1))*(x - x1) + y1;
            (x, y)
        },
        Border::XMax(xmax) => {
            let x = xmax;
            assert!(x1 != x2);
            let y = ((y2-y1)/(x2-x1))*(x - x1) + y1;
            (x, y)
        },
        Border::YMin(ymin) => {
            let y = ymin;
            assert!(y1 != y2);
            let x = ((x2-x1)/(y2-y1))*(y - y1) + x1;
            (x, y)
        },
        Border::YMax(ymax) => {
            let y = ymax;
            assert!(y1 != y2);
            let x = ((x2-x1)/(y2-y1))*(y - y1) + x1;
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

    // otherwise, we need to look at the intersections
    // I'm not happy with this, since you need to iterate over the linestring, twice and I wish
    // there was a way to do it only once.

    // FIXME little bit of code duplicating calculaing the AllInside and AllOutside

    #[allow(unused_assignments)]
    let mut last_inside: bool = false;

    let mut intersections: Vec<IntersectionOption<T>> = Vec::with_capacity(linestring.0.len() - 1);

    // first point
    let point_inside = is_inside(&linestring.0[0], border);
    all_inside &= point_inside;
    not_all_outside |= point_inside;
    last_inside = point_inside;
    let mut last_point = &linestring.0[0];

    for (idx, point) in linestring.0.iter().skip(1).enumerate() {
        let point_inside = is_inside(point, border);

        all_inside &= point_inside;
        not_all_outside |= point_inside;

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

    let all_outside = ! not_all_outside;
    if all_inside {
        LineBorderIntersection::AllInside
    } else if all_outside {
        LineBorderIntersection::AllOutside
    } else {
        LineBorderIntersection::Intersections(intersections)
    }
}


fn clip_point_to_border<T: CoordinateType>(p: &Point<T>, border: &Border<T>) -> Option<Geometry<T>> {
    if is_inside(p, border) {
        Some(Geometry::Point(p.clone()))
    } else {
        None
    }
}

fn clip_linestring_to_border<T: CoordinateType>(ls: &LineString<T>, border: &Border<T>) -> Option<Geometry<T>> {
    match calculate_intersections(ls, border) {
        LineBorderIntersection::AllInside => Some(Geometry::LineString(ls.clone())),
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


fn clip_multipoint_to_border<T: CoordinateType>(mp: &MultiPoint<T>, border: &Border<T>) -> Option<Geometry<T>> {
    let points: Vec<Point<T>> = mp.0.iter().filter_map(|p| clip_point_to_border(p, border).and_then(|g| g.as_point())).collect();
    if points.is_empty() {
        None
    } else {
        Some(Geometry::MultiPoint(MultiPoint(points)))
    }
}

fn clip_multilinestring_to_border<T: CoordinateType>(mls: &MultiLineString<T>, border: &Border<T>) -> Option<Geometry<T>> {
    let mut lines: Vec<LineString<T>> = Vec::with_capacity(mls.0.len());
    for ls in mls.0.iter() {
        match clip_linestring_to_border(ls, border) {
            None => {},
            Some(Geometry::LineString(new_ls)) => { lines.push(new_ls); },
            Some(Geometry::MultiLineString(new_mls)) => { lines.extend(new_mls); }
            _ => { unreachable!(); },
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
fn clip_to_border<T: CoordinateType>(geom: &Geometry<T>, border: &Border<T>) -> Option<Geometry<T>> {
    match *geom {
        Geometry::Point(ref p) => clip_point_to_border(p, border),
        Geometry::LineString(ref ls) => clip_linestring_to_border(ls, border),
        Geometry::Polygon(ref p) => sutherland_hodgeman::clip_polygon_to_border(p, border).map(Geometry::Polygon),
        Geometry::MultiPoint(ref mp) => clip_multipoint_to_border(mp, border),
        Geometry::MultiLineString(ref mls) => clip_multilinestring_to_border(mls, border),
        Geometry::MultiPolygon(ref p) => sutherland_hodgeman::clip_multipolygon_to_border(p, border).map(Geometry::MultiPolygon),
        Geometry::GeometryCollection(_) => unimplemented!(),
    }
}


pub fn clip_to_bbox<T: CoordinateType+::std::fmt::Debug>(geom: &Geometry<T>, bbox: &Bbox<T>) -> Option<Geometry<T>> {
    clip_to_border(geom, &Border::XMin(bbox.xmin))
       .and_then(|geom| clip_to_border(&geom, &Border::XMax(bbox.xmax)))
       .and_then(|geom| clip_to_border(&geom, &Border::YMin(bbox.ymin)))
       .and_then(|geom| clip_to_border(&geom, &Border::YMax(bbox.ymax)))
}

fn slice_box(geom: &Geometry<i32>, metatile_scale: u8, zoom: u8, tile_x0: u32, tile_y0: u32, x0: i32, y0: i32, size: i32) -> Vec<(slippy_map_tiles::Tile, Option<Geometry<i32>>)> {
    if metatile_scale == 1 {
        // FIXME do this at the caller and save a .clone()
        // Should do this in the caller, where metatile_scale == 2, where we own the geometry
        // I have a feeling there are too many clones here
        return vec![(slippy_map_tiles::Tile::new(zoom, tile_x0, tile_y0).unwrap(), Some(geom.clone()))];
    }

    let mut results = Vec::with_capacity((metatile_scale*metatile_scale) as usize);

    let half = size / 2;
    let tile_half = (metatile_scale/2) as u32;

    if let Some(left) = clip_to_border(geom, &Border::XMax(x0+half as i32)) {
        if let Some(topleft) = clip_to_border(geom, &Border::YMax(y0+half as i32)) {
            let mut tiles = slice_box(&topleft, metatile_scale/2, zoom, tile_x0, tile_y0, x0, y0, size/2);
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

        if let Some(bottomleft) = clip_to_border(geom, &Border::YMin(y0+half as i32)) {
            let mut tiles = slice_box(&bottomleft, metatile_scale/2, zoom, tile_x0, tile_y0+tile_half, x0, y0+half, size/2);
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

    if let Some(right) = clip_to_border(geom, &Border::XMin(x0+half as i32)) {
        if let Some(topright) = clip_to_border(geom, &Border::YMax(y0+half as i32)) {
            let mut tiles = slice_box(&topright, metatile_scale/2, zoom, tile_x0+tile_half, tile_y0, x0+half, y0, size/2);
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

        if let Some(bottomright) = clip_to_border(geom, &Border::YMin(y0+half as i32)) {
            let mut tiles = slice_box(&bottomright, metatile_scale/2, zoom, tile_x0+tile_half, tile_y0+tile_half, x0+half, y0+half, size/2);
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

pub fn clip_geometry_to_tiles(metatile: &Metatile, geom: &Geometry<i32>) -> Vec<(slippy_map_tiles::Tile, Option<Geometry<i32>>)> {
    slice_box(geom, metatile.scale(), metatile.zoom(), metatile.x(), metatile.y(), 0, 0, metatile.size() as i32*4096)
}
