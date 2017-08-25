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
    AllInside,
    AllOutside,
    Intersections(Vec<IntersectionOption<T>>),
}

#[derive(Clone, Debug, PartialEq, PartialOrd)]
enum IntersectionOption<T: CoordinateType> {
    Inside,
    Outside,
    Entry((T, T)),
    Exit((T, T)),
}


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
    let mut all_inside = true;
    let mut not_all_outside = false;

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

pub fn clip_geometry_to_tiles(metatile: &Metatile, geom: &Geometry<i32>) -> Vec<(slippy_map_tiles::Tile, Option<Geometry<i32>>)> {
    // this is a very simple solution, there are much better ways to do it:
    // Faster would be to do horizontal and vertical slices one at a time. So slice all the geoms
    // on the left and that's ½ the tiles done. there will be much less geometry stuff then.
    metatile.tiles().into_iter().map(|t| {
        let i = t.x() - metatile.x();
        let j = t.y() - metatile.y();
        // FIXME is the y the right way around?
        let bbox = Bbox{ xmin: (i*4096) as i32, xmax: ((i+1)*4096) as i32, ymin: (j*4096) as i32, ymax: ((j+1)*4096) as i32 };
        //println!("clip: tile {:?} bbox {:?} geom {:?}", t, bbox, geom);
        (t, clip_to_bbox(geom, &bbox))
    }).collect()
}
