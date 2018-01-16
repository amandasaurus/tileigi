use geo::*;
use geo::map_coords::MapCoords;
use geo::intersects::Intersects;
use std::cmp::{min, max, Ord, Ordering};
use std::ops::{Add, Sub};
use std::collections::HashMap;
use num_traits::Signed;
use std::fmt::Debug;
use std::hash::Hash;

use fraction::Fraction;

pub fn is_valid<T: CoordinateType+Signed+Debug+Ord>(geom: &Geometry<T>) -> bool {
    match *geom {
        Geometry::LineString(ref ls) => is_linestring_valid(ls),
        Geometry::Polygon(ref p) => is_polygon_valid(p),
        Geometry::MultiPolygon(ref mp) => mp.0.iter().all(|p| is_polygon_valid(p)),
        Geometry::MultiLineString(ref mls) => mls.0.iter().all(|ls| is_linestring_valid(ls)),
        _ => true,
    }
}

pub fn is_valid_skip_expensive<T: CoordinateType+Signed+Debug+Ord>(geom: &Geometry<T>) -> bool {
    match *geom {
        Geometry::LineString(ref ls) => is_linestring_valid(ls),
        Geometry::Polygon(ref p) => is_polygon_valid_skip_expensive(p),
        Geometry::MultiPolygon(ref mp) => mp.0.iter().all(|p| is_polygon_valid_skip_expensive(p)),
        Geometry::MultiLineString(ref mls) => mls.0.iter().all(|ls| is_linestring_valid(ls)),
        _ => true,
    }
}

fn is_linestring_valid<T: CoordinateType>(ls: &LineString<T>) -> bool {
    if ls.0.len() < 2 {
        return false;
    }

    if ls.0.len() == 2 && ls.0[0] == ls.0[1] {
        return false;
    }

    true
}

pub fn is_polygon_valid<T: CoordinateType+Signed+Debug+Ord>(p: &Polygon<T>) -> bool {
    is_polygon_valid_skip_expensive(p) && is_polygon_valid_do_expensive(p)
}

fn is_polygon_valid_skip_expensive<T: CoordinateType+Signed+Debug+Ord>(p: &Polygon<T>) -> bool {
    if p.exterior.0.len() < 4 {
        return false;
    }

    if p.exterior.0[0] != p.exterior.0[p.exterior.0.len()-1] {
        // first != last point
        return false;
    }

    // Sometimes there are duplicate points, e.g. A-A-B-A. If we remove all dupes, we can see if
    // there are <4 points
    if num_points_excl_duplicates(&p.exterior) < 4 {
        return false;
    }
    // TODO fix clipping code etc to not make linestrings with duplicated points

    if p.exterior.0.iter().skip(1).all(|&pt| pt == p.exterior.0[0]) {
        // All points the same
        // Shouldn't this be caught by the num_points_excl_duplicates ?
        return false;
    }


    for i in p.interiors.iter() {
        if num_points_excl_duplicates(i) < 4 {
            return false;
        }

        if i.0[0] != i.0[i.0.len()-1] {
            // first != last point
            return false;
        }

        if i.0.iter().skip(1).all(|&pt| pt == i.0[0]) {
            // All points the same
            return false;
        }

    }


    true
}

fn is_polygon_valid_do_expensive<T: CoordinateType+Signed+Debug+Ord>(p: &Polygon<T>) -> bool {
    if has_self_intersections(&p.exterior) {
        return false;
    }

    if p.interiors.iter().any(|i| has_self_intersections(i)) {
        return false;
    }

    // In theory this is backwards. Ext rings should be CCW, and int rings CW. But in vector tiles
    // the y goes down, so it's flipped.
    if p.exterior.is_ccw() || p.interiors.iter().any(|i| i.is_cw()) {
        return false;
    }

    true
}

fn linestring_has_duplicate_points<T: CoordinateType>(ls: &LineString<T>) -> bool {
    (0..ls.0.len()-1).into_iter().any(|i| ls.0[i] == ls.0[i+1])
}

fn has_duplicate_points<T: CoordinateType>(geom: &Geometry<T>) -> bool {
    match *geom {
        Geometry::Point(_) => false,
        Geometry::MultiPoint(_) => false,
        Geometry::LineString(ref ls) => linestring_has_duplicate_points(ls),
        Geometry::MultiLineString(ref mls) => mls.0.iter().any(|l| linestring_has_duplicate_points(l)),
        Geometry::Polygon(ref p) => {
            linestring_has_duplicate_points(&p.exterior) || p.interiors.iter().any(|l| linestring_has_duplicate_points(l))
        },
        Geometry::MultiPolygon(ref mp) => {
            // This is silly why can't we call this function 
            mp.0.iter().any(|p| linestring_has_duplicate_points(&p.exterior) || p.interiors.iter().any(|l| linestring_has_duplicate_points(l)) )
        },
        Geometry::GeometryCollection(ref gc) => gc.0.iter().any(|g| has_duplicate_points(g)),
    }
}

/// Returns the number of points in this line if you were to remove all consequetive duplicate
/// points. If this is <4 then it's not valid for a ring.
fn num_points_excl_duplicates<T: CoordinateType>(ls: &LineString<T>) -> usize {
    if ls.0.len() <= 1 { return ls.0.len(); }

    let mut curr_point_idx = 0;
    let mut num = 1;
    for i in 1..ls.0.len() {
        if ls.0[i] != ls.0[curr_point_idx] {
            curr_point_idx = i;
            num += 1;
        }
    }

    num

}




pub fn ensure_polygon_orientation<T: CoordinateType>(geom: &mut Geometry<T>) {
    match *geom {
        Geometry::Polygon(ref mut p) => {
            // This is stupid, this is supposed to be the other way around!!
            // FIXME check the geo code for winding, make sure it's not the wrong way around
            // Is this because in MVT the Y goes down? Hence the winding is the other way?
            // Investigate
            p.exterior.make_clockwise_winding();
            for i in p.interiors.iter_mut() {
                i.make_counterclockwise_winding();
            }
        },
        Geometry::MultiPolygon(ref mut mp) => {
            for p in mp.0.iter_mut() {
                p.exterior.make_clockwise_winding();
                for i in p.interiors.iter_mut() {
                    i.make_counterclockwise_winding();
                }
            }
        },
        _ => {},
    }
}

fn has_self_intersections<T: CoordinateType+Signed+Debug+Ord>(ls: &LineString<T>) -> bool {
    if ls.0.len() <= 4 {
        // cannot have a self intersection with this few members. (There shouldn't be <4 anyway)
        // With 4 points, it's a orientation, not self-intersection thing really
        return false;
    }
    //println!("\n\nXXX\nls {:?}\n", ls);

    for (i, points12) in ls.0.windows(2).enumerate() {
        let (p1, p2) = (points12[0], points12[1]);
        
        for points34 in ls.0[i+1..].windows(2).take(ls.0.len()-i-1) {
            let (p3, p4) = (points34[0], points34[1]);
            //println!("looking at i {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", i, p1, p2, p3, p4);

            match intersection(p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y(), p4.x(), p4.y()) {
                Intersection::Crossing(_) | Intersection::Overlapping(_, _)  => { return true; },
                Intersection::Touching(_) => { return true; },
                Intersection::None | Intersection::EndToEnd => {},
            }
            //println!("no intersection");
        }
    }


    false
}

#[derive(PartialEq,Eq,Clone,Copy,Debug)]
enum Intersection<T> {
    // They don't intersect/touch at all
    None,

    // One is wholly, or partially, on top of another, ie infinite number of intersecting points,
    // the intersection is a line
    // These points are the points where the overlap starts and ends
    Overlapping((T, T), (T, T)),

    // The end point of one is the same as the end point of another,
    EndToEnd,

    // The end of one touches the other, but not at it's end, at point (T, T)
    Touching((T, T)),

    // real crossing, at point (T, T)
    Crossing((T, T))
}

fn intersect_incl_end<T: CoordinateType+Signed+Debug+Ord>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T, x4: T, y4: T) -> bool {
    intersection(x1, y1, x2, y2, x3, y3, x4, y4) == Intersection::None
}


/// True iff the segments |p1p2| and |p3p4| intersect at any point, and the intersection point is
/// not on both end points. i.e. 2 lines can join end-to-end in this, but not touch anywhere else.
fn intersection<T: CoordinateType+Signed+Debug+Ord>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T, x4: T, y4: T) -> Intersection<T> {
    if max(x1, x2) < min(x3, x4) || min(x1, x2) > max(x3, x4)
        || max(y1, y2) < min(y3, y4) || min(y1, y2) > max(y3, y4)
    {
        return Intersection::None;
    }
    
    //println!("\nline12 ({:?}, {:?}) - ({:?}, {:?})", x1, y1, x2, y2);
    //println!("line34 ({:?}, {:?}) - ({:?}, {:?})", x3, y3, x4, y4);

    debug_assert!((x1, y1) != (x2, y2));
    debug_assert!((x3, y3) != (x4, y4));

    let a = x2 - x1;
    let b = x3 - x4;
    let c = y2 - y1;
    let d = y3 - y4;

    let determinate = a*d - b*c;
    if determinate == T::zero() {
        // TODO should probably profile & optimize this bit
        // Slope of line12 is a/c, slope of line34 is b/d. Lines are parallel/colinear if a/c =
        // b/d, i.e.  a*d - b*c == 0
        // This branch is when the slopes are the same

        // The lines are the same (if we ignore direction). One lies totally on top of the
        // other
        if ((x1, y1) == (x3, y3) && (x2, y2) == (x4, y4)) || ((x1, y1) == (x4, y4) && (x2, y2) == (x3, y3)) {
            return Intersection::Overlapping((x1, y1), (x2, y2));
        }

        fn delta<T: Ord+Sub<Output=T>>(a: T, b: T) -> T {
            if a > b { a - b } else { b - a }
        }

        let delta_x = delta(x1, x2)+delta(x3, x4);
        let delta_y = delta(y1, y2)+delta(y3, y4);
        if    (delta_x == delta(x1, x4) && delta_y == delta(y1, y4)) // 1-23-4
           || (delta_x == delta(x2, x4) && delta_y == delta(y2, y4)) // 2-13-4
           || (delta_x == delta(x1, x3) && delta_y == delta(y1, y3)) // 1-24-3
           || (delta_x == delta(x2, x3) && delta_y == delta(y2, y3)) // 2-13-4
            {
            // One after the other. We know they have the same slope, so this shortcut calculation
            // works.
            return Intersection::EndToEnd;
        }

        let slope_12 = Fraction::new(a, c);
        let slope_13 = Fraction::new(x3-x1, y3-y1);

        if !slope_12.same_slope(&slope_13) {
            return Intersection::None;
        }

        let slope_34 = Fraction::new(b, d);

        let slope_23 = Fraction::new(x3-x1, y3-y1);
        let slope_14 = Fraction::new(x4-x1, y4-y1);
        let slope_24 = Fraction::new(x4-x2, y4-y2);

        // TODO Maybe remove Copy & use reference &U? Would that mean less memory copies?
        fn in_bounds<U: Ord+Copy>(z: U, a: U, b: U) -> bool {
            z >= min(a,b) && z <= max(a,b)
        }

        let p3_on_end = (x1, y1) == (x3, y3) || (x2, y2) == (x3, y3);
        let p4_on_end = (x1, y1) == (x4, y4) || (x2, y2) == (x4, y4);
        let p1_on_34 = in_bounds(x1, x3, x4) && in_bounds(y1, y3, y4) && slope_13.same_slope(&slope_34);
        let p2_on_34 = in_bounds(x2, x3, x4) && in_bounds(y2, y3, y4) && slope_23.same_slope(&slope_34);
        let p3_on_12 = in_bounds(x3, x1, x2) && in_bounds(y3, y1, y2) && slope_13.same_slope(&slope_12);
        let p4_on_12 = in_bounds(x4, x1, x2) && in_bounds(y4, y1, y2) && slope_14.same_slope(&slope_12);

        //println!("line12 ({:?}, {:?}) - ({:?}, {:?})", x1, y1, x2, y2);
        //println!("line34 ({:?}, {:?}) - ({:?}, {:?})", x3, y3, x4, y4);
        //println!("p3_on_end {} p3_on_12 {}\np4_on_end {} p4_on_12 {}", p3_on_end, p3_on_12, p4_on_end, p4_on_12);
        //println!("p1_on_34 {} p2_on_34 {}", p1_on_34, p2_on_34);
        //println!("p3_on_12 {} p4_on_12 {}", p3_on_12, p4_on_12);

        if ((x1, y1) == (x3, y3) && (x2, y2) == (x4, y4)) || ((x1, y1) == (x4, y4) && (x2, y2) == (x3, y3)) {
            return Intersection::Overlapping((x1, y1), (x2, y2));
        } else {
            let first_point = if (x1,y1) == (x3, y3) || p1_on_34 { (x1, y1) } else if p3_on_12 { (x3,y3) } else { return Intersection::None; };
            let last_point = if (x2,y2) == (x4, y4) || p2_on_34 { (x2, y2) } else if p4_on_12 { (x4,y4) } else { return Intersection::None; };
            debug_assert!(first_point != last_point, "{:?}", first_point);
            return Intersection::Overlapping(first_point, last_point);
        }
    }

    let e = x3 - x1;
    let f = y3 - y1;

    // we know it's not zero
    let signum = determinate.signum();
    let determinate = determinate.abs();

    let sd = signum * (a*f - c*e);
    if sd > determinate || sd < T::zero() {
        return Intersection::None;
    }

    let td = signum*(d*e - b*f);
    if td > determinate || td < T::zero() {
        return Intersection::None;
    }

    if (td == determinate || td == T::zero()) && (sd == T::zero() || sd == determinate) {
        // endpoints overlap
        return Intersection::EndToEnd;
    } else if (td == determinate || td == T::zero()) && (sd > T::zero() || sd < determinate) {
        if td == T::zero() {
            return Intersection::Touching((x1, y1));
        } else if td == determinate {
            return Intersection::Touching((x2, y2));
        } else {
            unreachable!();
        }
    } else if (td < determinate || td > T::zero()) && (sd == T::zero() || sd == determinate) {
        if sd == T::zero() {
            return Intersection::Touching((x3, y3));
        } else if sd == determinate {
            return Intersection::Touching((x4, y4));
        } else {
            unreachable!();
        }
    } else if td > T::zero() && td < determinate && sd > T::zero() && sd < determinate {
        // This will do some roundingin on integers
        let xd = determinate*x1 + td*(x2 - x1);
        let mut x = xd/determinate;
        let yd = determinate*y1 + td*(y2 - y1);
        let mut y = yd/determinate;

        //println!("td {:?} sd {:?} determinate {:?}", td, sd, determinate);
        //println!("xd {:?} yd {:?}", xd, yd);

        // Do regular rounding on the integers (i.e. [0,0.5) is rounded down, [0.5, 1) is rounded
        // up.
        // Look at the remained from *d/determinate, and if it's more than half the value of
        // determinate (or twice it is more than determinate), then the first decimal place would
        // be above 5, ergo we should round up. i.e. we add one to the current numbers
        let two = T::one() + T::one();
        let twice_x_remainder = two*(xd % determinate);
        if twice_x_remainder >= determinate {
            x = x + T::one();
        }

        let twice_y_remainder = two*(yd % determinate);
        if twice_y_remainder >= determinate {
            y = y + T::one();
        }
        //println!("twice_x_remainder {:?} twice_y_remainder {:?}", twice_x_remainder, twice_y_remainder);

        return Intersection::Crossing((x, y));
    }

    // Should have been caught above.
    eprintln!("points {:?} {:?} {:?} {:?}", (x1, y1), (x2, y2), (x3, y3), (x4, y4));
    eprintln!("det {:?} sd {:?} td {:?}", determinate, sd, td);
    unreachable!();
}

pub fn make_valid<T: CoordinateType+Debug+Ord+Signed+Hash>(geom: Geometry<T>) -> Geometry<T> {
    match geom {
        Geometry::Polygon(p) => Geometry::MultiPolygon(make_polygon_valid(p)),
        Geometry::MultiPolygon(mp) => Geometry::MultiPolygon(make_multipolygon_valid(mp)),
        x => x,
    }
}

fn make_multipolygon_valid<T: CoordinateType+Debug+Ord+Signed+Hash>(mut mp: MultiPolygon<T>) -> MultiPolygon<T> {
    let MultiPolygon( polygons ) = mp;

    let rings: Vec<LineString<T>> = polygons.into_iter().flat_map(|p| {
        let Polygon{ exterior, interiors } = p;
        let mut these_rings = interiors;
        these_rings.insert(0, exterior);
        these_rings.into_iter()
    }).collect();


    make_rings_valid(rings)
}

fn make_polygon_valid<T: CoordinateType+Debug+Ord+Signed+Hash>(mut p: Polygon<T>) -> MultiPolygon<T> {
    let Polygon{ exterior, interiors } = p;
    let mut rings = interiors;
    rings.insert(0, exterior);

    make_rings_valid(rings)

}

fn make_rings_valid<T: CoordinateType+Debug+Ord+Signed+Hash>(mut rings: Vec<LineString<T>>) -> MultiPolygon<T> {

    let rings: Vec<LineString<T>> = rings.into_iter().flat_map(|mut ring| {
        add_points_for_all_crossings(&mut ring);
        let these_rings = dissolve_into_rings(ring);
        these_rings.into_iter()
    }).collect();

    
    let result = convert_rings_to_polygons(rings);

    // This takes a geom, so we do a dance
    let mut result = Geometry::MultiPolygon(result);
    ensure_polygon_orientation(&mut result);

    if let Geometry::MultiPolygon(mp) = result {
        return mp;
    } else {
        unreachable!()
    }
}


/// Modify the LineString, so that at all self-intersection places there is a node. i.e. if 2
/// segments cross, add a node in the middle of each segment where they cross. After this all
/// self-intersections will be of the EndToEnd type
fn add_points_for_all_crossings<T: CoordinateType+Debug+Signed+Ord>(ls: &mut LineString<T>) {
    if ls.0.len() <= 3 {
        return;
    }
    //println!("\n\nXXX\nls {:?}\n", ls);

    loop {
        //println!("\nStart of loop\n{:?}", ls.0);
        let mut coords_to_insert = HashMap::new();
        // Keys are the point indexes.
        // Values are a Vec of new points to add after the point with that index.
        // So vec![(0, 0), (1, 0)] for key #3, means to insert those 2 points after ls.0[3]
        // They are initially stored in the order they appear in, but they need to be sorted
        // afterwards

        for (i, points12) in ls.0.windows(2).enumerate() {
            let (p1, p2) = (points12[0], points12[1]);
            
            for (j, points34) in ls.0[i+1..].windows(2).enumerate().take(ls.0.len()-i-1) {
                let j = j + i + 1;
                let (p3, p4) = (points34[0], points34[1]);
                let x1 = p1.x(); let y1 = p1.y();
                let x2 = p2.x(); let y2 = p2.y();
                let x3 = p3.x(); let y3 = p3.y();
                let x4 = p4.x(); let y4 = p4.y();

                match intersection(x1, y1, x2, y2, x3, y3, x4, y4) {
                    Intersection::None | Intersection::EndToEnd => {},

                    Intersection::Crossing((x0, y0)) => {
                        //println!("looking at i {} j {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", i, j, p1, p2, p3, p4);
                        //println!("i {} j {} crossing {:?},{:?}", i, j, x0, y0);
                        coords_to_insert.entry(i).or_insert(vec![]).push((x0, y0));
                        coords_to_insert.entry(j).or_insert(vec![]).push((x0, y0));
                    },

                    Intersection::Overlapping(overlap1, overlap2)  => {
                        //println!("looking at i {} j {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", i, j, p1, p2, p3, p4);
                        //println!("i {} j {} overlapping {:?},{:?}", i, j, overlap1, overlap2);
                        debug_assert!(overlap1 != overlap2);

                        if (x1, y1) != overlap1 && (x2, y2) != overlap1 {
                            coords_to_insert.entry(i).or_insert(vec![]).push(overlap1);
                        }
                        if (x1, y1) != overlap2 && (x2, y2) != overlap2 {
                            coords_to_insert.entry(i).or_insert(vec![]).push(overlap2);
                        }

                        if (x3, y3) != overlap1 && (x4, y4) != overlap1 {
                            coords_to_insert.entry(j).or_insert(vec![]).push(overlap1);
                        }
                        if (x3, y3) != overlap2 && (x4, y4) != overlap2 {
                            coords_to_insert.entry(j).or_insert(vec![]).push(overlap2);
                        }
                    },

                    Intersection::Touching((x0, y0)) => {
                        //println!("looking at i {} j {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", i, j, p1, p2, p3, p4);
                        //println!("i {} j {} touching {:?},{:?}", i, j, x0, y0);
                        // (x0, y0) is the point where they touch
                        if (x1,y1) == (x0,y0) || (x2,y2) == (x0,y0) {
                            // touching point is at end of line12, ergo it's in the middle of line34
                            coords_to_insert.entry(j).or_insert(vec![]).push((x0, y0));
                        } else if (x3,y3) == (x0,y0) || (x4,y4) == (x0,y0) {
                            coords_to_insert.entry(i).or_insert(vec![]).push((x0, y0));
                        } else {
                            unreachable!();
                        }
                    }
                }
            }
        }


        if coords_to_insert.is_empty() {
            break;
        } else {
            // When we insert a point into the vec, it'll push all after that along. Keep track of
            // how many we've inserted.
            let mut offset = 0;

            // Turn hashmap into a sorted vec, sorted by index to add
            let coords_to_insert = ls.0.windows(2).enumerate().filter_map(|(idx, points)| {
                let (point1, point2) = (points[0], points[1]);
                if let Some(mut new_points) = coords_to_insert.remove(&idx) {
                    new_points.sort_by(|&new_coord1, &new_coord2| order_points(((point1.x(), point1.y()), (point2.x(), point2.y())), new_coord1, new_coord2));
                    new_points.dedup();
                    Some((idx, new_points))
                } else {
                    None
                }
            }).collect::<Vec<_>>();

            //println!("line {:?}", ls);
            //println!("coords_to_insert {:?}", coords_to_insert);

            for (point_idx, new_points) in coords_to_insert.into_iter() {
                for new_point in new_points.into_iter() {
                    //println!("Adding {:?} after index {}", new_point, point_idx+offset);
                    // +1 because we want the new point to be *after* the current point we're
                    // looking at
                    ls.0.insert(point_idx+offset+1, Point::new(new_point.0, new_point.1));
                    offset += 1;
                }
            }

            //println!("We added {} new points to the line", offset);
        }
    }

    //println!("finished");

}

fn dissolve_into_rings<T: CoordinateType+Debug+Hash+Eq>(ls: LineString<T>) -> Vec<LineString<T>> {
    let LineString( points ) = ls;
    if points.len() <= 3 {
        // Not enough points for a proper ring
        return vec![];
    }

    let mut outgoing_segments = HashMap::with_capacity(points.len());

    for (i, p) in points.iter().enumerate() {
        // TODO here we could assert that the existing vec is <=2, and generate the loops vec here,
        // rather than do a loop later
        outgoing_segments.entry((p.x(), p.y())).or_insert(vec![]).push(i);
    }

    //println!("outgoing_segments {:?}", outgoing_segments);
    let mut loops = outgoing_segments.iter().filter(|&(_, v)| v.len() > 1).map(|(_, v)| v).collect::<Vec<_>>();

    if loops.len() == 1 {
        if loops[0][0] == 0 && loops[0][1] == points.len()-1 {
            // start & end
            return vec![LineString(points)];
        } else {
            eprintln!("points {:?}", points);
            eprintln!("loops {:?}", loops);
            unreachable!();
        }
    }

    let mut point_unassigned = vec![true; points.len()];
    let mut results: Vec<LineString<T>> = vec![];

    debug_assert!(loops.iter().all(|idxes| idxes.len() == 2));

    //println!("points {:?}", points);
    //println!("outgoing_segments {:?}", outgoing_segments);
    // sort loops where the smaller length (in terms of number of points) are to the front.
    // Ideal: Sort them so that if a loop is a subset of a larger loop, then the smaller is ahead,
    // so the smaller, "inner" loop will be removed first. Unless something really strange is going
    // on, this sort-by-length should do it (since an outer loop will be longer than the inner one
    // it contains)
    loops.sort_by_key(|i| (i[1]-i[0], i[0]));
    //println!("loops {:?}", loops);

    // FIXME fix this
    //assert!(loops.len() == 2);

    for loop_indexes in loops {
        assert!(loop_indexes.len() == 2);
        let (start, end) = (loop_indexes[0], loop_indexes[1]);
        if !point_unassigned[start] {
            // this has already been removed earlier in another loop
            continue;
        }
        //println!("loop from {:?} to {:?}", start, end);

        if start + 2 == end {
            // This is only 3 points, so it's a little spike
            // Don't include it, and ensure the points are skipped
            point_unassigned[start] = false;
            point_unassigned[start+1] = false;
            continue;
        }
        
        let mut new_ls = vec![];
        point_unassigned[start] = false;
        //points[start..end].iter_mut().map( set to true here? )
        new_ls.push(points[start].clone());
        for i in start+1..end {
            //println!("looking at point i {} {:?}, point_unassigned[i] {:?}", i, points[i], point_unassigned[i]);
            if point_unassigned[i] {
                //println!("\tadding point");
                new_ls.push(points[i].clone());
                point_unassigned[i] = false;
            }
        }
        if new_ls.len() > 2 {
            // Any outer loops need at least one point at this, so don't save it
            //point_unassigned[end] = false;
            new_ls.push(points[end].clone());
            //println!("adding ls {:?}", new_ls);
            results.push(LineString(new_ls));
        } else {
            //println!("too short");
        }
    }

    //println!("point_unassigned {:?}", point_unassigned);
    // There will always be the last/first point unassigned since we keep the end around, which
    // means the endpoint of the outer ring is kept. So they should all be false, except the last
    // which is true
    debug_assert!(point_unassigned.iter().take(point_unassigned.len()-1).all(|x| !x));
    debug_assert!(point_unassigned[point_unassigned.len()-1]);

    results
}

/// Possible return values from does_ray_cross
#[derive(PartialEq,Eq,Debug)]
enum Crossing {
    /// Definitly no overlap
    No,

    /// There is a specific overlap
    Yes,

    /// The point is on the line segment, i.e it's a horizontal line and the ray passes
    /// along/through it
    Touches,

    /// The ray goes through the start, or end, point of the line
    StartPoint,
    EndPoint,
}

/// An infinite line from point to the left (ie negative infitity in the x direction), does that
/// line intersect with the line segment from p1-p2?
fn does_ray_cross<T: CoordinateType+Debug+Ord>(point: &Point<T>, p1: &Point<T>, p2: &Point<T>) -> Crossing {
    let (x, y) = (point.x(), point.y());
    assert!(p1 != p2);
    let (x1, y1) = (p1.x(), p1.y());
    let (x2, y2) = (p2.x(), p2.y());

    //println!("x {:?} y {:?} x1 {:?} y1 {:?} x2 {:?} y2 {:?}", x, y, x1, y1, x2, y2);
    if (x == x1 && y == y1) || (x == x2 && y == y2) {
        return Crossing::Touches;
    }

    if y1 > y && y2 > y {
        // Line segment is above the point
        return Crossing::No;
    } else if y1 < y && y2 < y {
        // Line segment is below the point
        return Crossing::No;
    } else if x1 > x && x2 > x {
        // Line segment is to the right of the point
        return Crossing::No;
    } else if (y1 > y && y2 < y) || (y1 < y && y2 > y) {
        // proper crossing
        return Crossing::Yes;
    } else if x1 < x && x2 < x {
        // This linesegment is to the left of the poing

        if y1 != y2 && y1 == y {
            return Crossing::StartPoint;
        } else if y1 != y2 && y2 == y {
            return Crossing::EndPoint;
        } else if (y1 == y && y2 != y) || (y1 != y && y2 == y) {
            return Crossing::Yes;
        } else if y1 == y && y2 == y {
            // This line lies on the ray
            return Crossing::Touches;
        }
    }

    unreachable!();
}


#[derive(PartialEq,Eq,Debug)]
enum RingType { Exterior, Interior }

/// ring is at index `ring_type` in `all_rings`
fn is_ring_ext_int<T: CoordinateType+Debug+Ord>(ring: &LineString<T>, ring_index: usize, all_rings: &Vec<LineString<T>>) -> RingType {
    // Do an even/odd check on a point in `ring` on all rings in all_rings. except this one (that's
    // why we need ring_index. If the point is inside, then this is an interior ring, else
    // exterior.
    // We pick the first point in ring, but if we get a "touch" relation against any other ring, we
    // just move on to another point.
    // We assume that a ring is either entirely inside, or entirely outside another ring. There are
    // no "partially overlapping" rings.
    let point = ring.0[0];
    let mut num_crossings = 0;
    //println!("{} is_ring_ext_int\nring {:?}\nring_index {}\nall_rings{:?}", line!(), ring, ring_index, all_rings);

    'start_point: for point in ring.0.iter() {
        num_crossings = 0;

        
        // loop over all the rings
        for (i, ring) in all_rings.iter().enumerate() {
            if i == ring_index { continue; }
            //println!("i {} point {:?}", i, point);

            // then all the segments in this ring
            for other_points in ring.0.windows(2) {
                debug_assert!(other_points.len() == 2);

                //println!("other_points {:?}, does_ray_cross {:?}", other_points, does_ray_cross(&point, &other_points[0], &other_points[1]));
                match does_ray_cross(&point, &other_points[0], &other_points[1]) {
                    // Choose that when it goes through the start, it's a cross. Otherwise we could
                    // double count the crossings when the point is at the same y value as a point
                    // on this ring.
                    Crossing::StartPoint => { num_crossings += 1 },
                    Crossing::EndPoint => {},

                    Crossing::Yes => { num_crossings += 1 },
                    Crossing::No => {},
                    Crossing::Touches => {
                        //println!("Touches, so try again");
                        // Go back and choose a new start point
                        continue 'start_point;
                    }
                }
            }
        }

        // If we've gotten to here, this start point is good.
        break 'start_point;
    }

    if num_crossings % 2 == 0 {
        RingType::Exterior
    } else {
        RingType::Interior
    }

}

fn calc_rings_ext_int<T: CoordinateType+Debug+Ord>(rings: Vec<LineString<T>>) -> Vec<(LineString<T>, RingType)> {
    let ring_types: Vec<RingType> = rings.iter().enumerate().map(|(i, r)| is_ring_ext_int(&r, i, &rings) ).collect();

    rings.into_iter().zip(ring_types.into_iter()).collect()


}

/// This will look at what rings are inside other rings.
fn convert_rings_to_polygons<T: CoordinateType+Debug+Ord>(mut rings: Vec<LineString<T>>) -> MultiPolygon<T> {
    if rings.is_empty() {
        return MultiPolygon(vec![]);
    }
    if rings.len() == 1 {
        return MultiPolygon(vec![Polygon::new(rings.remove(0), vec![])]);
    }


    let rings_with_type = calc_rings_ext_int(rings);

    //println!("{} rings_with_type {:?}", line!(), rings_with_type);

    // Do a simple case when there are only 2 rings?
    let mut exteriors = Vec::new();
    let mut interiors = Vec::new();

    for (ring, ring_type) in rings_with_type.into_iter() {
        match ring_type {
            RingType::Exterior => {
                exteriors.push(ring);
            },
            RingType::Interior => {
                interiors.push(ring);
            },
        }
    }
    assert!(!(exteriors.is_empty() && interiors.is_empty()));

    // All interiours?!
    assert!(!exteriors.is_empty());

    let mut polygons: Vec<_> = exteriors.into_iter().map(|p| Polygon::new(p, vec![])).collect();

    // we need to calculate the what exterior that each interior is in
    
    if polygons.len() == 1 {
        // There is only one exterior ring, so take a simple approach of assuming all the
        // interiors are part of that
        ::std::mem::swap(&mut polygons[0].interiors, &mut interiors);
        
    } else {
        if interiors.is_empty() {
            // nothing to do
        } else {
            // we need to figure out which exterior each interior is in.
            eprintln!("polygons {:?}", polygons);
            unimplemented!()
        }
    }


    MultiPolygon(polygons)
}

/// Given a line defined by 2 points, and 2 other points (p1 & p2) which were assume are on the
/// line, return where those 2 points are in order when going along the line, or not.
/// Returns Ordering::Equal when the points are the same.
/// Returns Ordering::Less when p1 comes before p2 when moving from the line's start to it's end,
/// i.e. they are sorta in order already.
/// Returns Ordering::Greater when p1 comes after p2 when moving from the line's start to it's end
fn order_points<T: CoordinateType+Debug+Sub<Output=T>+Ord>(line: ((T, T), (T, T)), p1: (T, T), p2: (T, T)) -> Ordering {
    if p1 == p2 {
        return Ordering::Equal;
    }
    assert!(line.0 != line.1);

    fn sub<T: CoordinateType+Ord+Sub<Output=T>>(a: (T, T), b: (T, T)) -> (T, T) {
        (
            match a.0.cmp(&b.0) {
                Ordering::Equal => T::zero(),
                Ordering::Greater => (a.0 - b.0),
                Ordering::Less => (b.0 - a.0),
            },
            match a.1.cmp(&b.1) {
                Ordering::Equal => T::zero(),
                Ordering::Greater => (a.1 - b.1),
                Ordering::Less => (b.1 - a.1),
            },

        )
    }

    fn add<T: Add<Output=T>>(a: (T, T), b: (T, T), c: (T, T)) -> (T, T) {
        (a.0+b.0+c.0, a.1+b.1+c.1)
    }

    // we 'abs' all the slopes so that the line is entirely in the first quarter
    // (delta x, delta y)
    let slope_line = sub(line.1, line.0);

    // slope from the start point to p1
    let slope_start_1 = sub(p1, line.0);//(p1.0 - (line.0).0, p1.1 - (line.0).1);

    // slope from the start point to p2
    let slope_start_2 = sub(p2, line.0);//(p2.0 - (line.0).0, p2.1 - (line.0).1);
    
    // slope from p1 to p2
    let slope_1_2 = sub(p2, p1);//(p2.0 - p1.0, p2.1 - p1.1);

    // slope from p2 to p1
    let slope_2_1 = sub(p1, p2);

    // slope from p2 to the end
    let slope_2_end = sub(line.1, p2);//((line.1).0 - p2.0, (line.1).1 - p2.1);
    let slope_1_end = sub(line.1, p1);

    if add(slope_start_1, slope_1_2, slope_2_end) == slope_line {
        Ordering::Less
        // (p1, p2)
    } else if add(slope_start_2, slope_2_1, slope_1_end) == slope_line {
        Ordering::Greater
        // (p2, p1)
    } else {
        // this shouldn't happen
        // Probably happens when p1 and/or p2 aren't on the line
        eprintln!("line {:?} p1 {:?} p2 {:?}", line, p1, p2);
        eprintln!("slone_line {:?}", slope_line);
        eprintln!("slope_start_1 {:?} slope_start_2 {:?}", slope_start_1, slope_start_2);
        eprintln!("slope_1_2 {:?} slope_2_1 {:?}", slope_1_2, slope_2_1);
        eprintln!("slope_2_end {:?}", slope_2_end);
        unreachable!();
    }

}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn intersect1() {

        assert_eq!(intersection(0, 0,  0, 10,  5, 1,  5, 2), Intersection::None);
        assert_eq!(intersection(0, 0,  0, 10,  0, 5,  5, 5), Intersection::Touching((0, 5)));

        assert_eq!(intersection(0, 0,  0, 10,  0, 0,  0, 10), Intersection::Overlapping((0, 0), (0, 10)));
        assert_eq!(intersection(0, 0,  0, 10,  0, 5,  0, 10), Intersection::Overlapping((0, 5), (0, 10)));
        assert_eq!(intersection(0, 0,  0, 10,  0, 5,  0, 15), Intersection::Overlapping((0, 5), (0, 10)));
        assert_eq!(intersection(0, 0,  0, 10,  0, 0,  0, 5), Intersection::Overlapping((0, 0), (0, 5)));
        assert_eq!(intersection(0, 0,  0, 10,  0, 2,  0, 8), Intersection::Overlapping((0, 2), (0, 8)));
        assert_eq!(intersection(0,2, 0,8,  0,0, 0,10), Intersection::Overlapping((0, 2), (0, 8)));

        assert_eq!(intersection(0, 0, 0, 10,  0, 10,  1, 20), Intersection::EndToEnd);
        assert_eq!(intersection(0, 0, 0, 10,  0, 10,  0, 20), Intersection::EndToEnd);
        assert_eq!(intersection(0, 0, 0, 10,  1, 10,  0, 10), Intersection::EndToEnd);


        //   C
        // A-B
        // test all combinations
        // ABBC
        assert_eq!(intersection(0,0, 0,1,  0,1, 1,1), Intersection::EndToEnd);
        // ABCB
        assert_eq!(intersection(0,0, 0,1,  1,1, 0,1), Intersection::EndToEnd);
        // BABC
        assert_eq!(intersection(0,1, 0,0,  0,1, 1,1), Intersection::EndToEnd);
        // BACB
        assert_eq!(intersection(0,1, 0,0,  1,1, 0,1), Intersection::EndToEnd);
        

        assert_eq!(intersection(0, 0, 0, 10,  1, 10,  1, 20), Intersection::None);
        assert_eq!(intersection(0, 0, 0, 10,  1, 20,  1, 40), Intersection::None);

        assert_eq!(intersection(0, 0,  0, 10,  -5, 5,  5, 5), Intersection::Crossing((0, 5)));
        assert_eq!(intersection(0, 0,  0, 10,  -5, 1,  5, 1), Intersection::Crossing((0, 1)));

        assert_eq!(intersection(0, 0,  10, 0,  10, 0,  10, 10), Intersection::EndToEnd);
        assert_eq!(intersection(-5, 5,  5, 5,  0, 0,  0, 10), Intersection::Crossing((0, 5)));
        assert_eq!(intersection(0, 0,  10, 0,  5, 10,  5, -10), Intersection::Crossing((5, 0)));

        // Rounded down to the nearest whole number in the coordinate system
        // They meet at (0.5, 0.5), so that's rounded to (1, 1)
        assert_eq!(intersection(0,0, 1,1,  1,0, 0,1), Intersection::Crossing((1, 1)));
    }

    #[test]
    fn intersect2() {
        assert!(!has_self_intersections(&vec![(0, 0), (1, 0)].into()));
        assert!(!has_self_intersections(&vec![(0, 0), (1, 0), (2, 0)].into()));
        assert!(!has_self_intersections(&vec![(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0)].into()));

        // should be a problem
        assert!(has_self_intersections(&vec![(0, 0), (10, 0), (10, 10), (5, 10), (5, -10)].into()));

        // closed ring, should be OK
        assert!(!has_self_intersections(&vec![(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)].into()));
    }

    #[test]
    fn intersect3() { assert_eq!(intersection(4,0, 2,-1,  2,1, 0,0), Intersection::None); }

    #[test]
    fn intersect4() {
        assert_eq!(intersection(0,0, 4,0,  2,-1, 2,0), Intersection::Touching((2, 0)));
        assert_eq!(intersection(0,0, 4,0,  2,0, 2,1), Intersection::Touching((2, 0)));

        assert_eq!(intersection(2,-1, 2,0,  0,0, 4,0), Intersection::Touching((2, 0)));
        assert_eq!(intersection(2,0, 2,1,  0,0, 4,0), Intersection::Touching((2, 0)));
    }

    #[test]
    fn intersect5() {
        assert_eq!(intersection(0,0, 4,0,  1,-1, 1,1), Intersection::Crossing((1, 0)));
        assert_eq!(intersection(0,0, 4,0,  2,-1, 2,1), Intersection::Crossing((2, 0)));
    }

    #[test]
    fn intersect6() {
        // bbox overlaps, and they have the same slope, but they don't touch
        assert_eq!(intersection(0,0, 10,10,  1,2, 6,7), Intersection::None);
        assert_eq!(intersection(1,2, 6,7,  0,0, 10,10), Intersection::None);
        assert_eq!(intersection(10,10, 0,0,  6,7, 1,2), Intersection::None);
        assert_eq!(intersection(6,7, 1,2,  10,10, 0,0), Intersection::None);
    }

    #[test]
    fn intersect7() {
        // bbox overlaps, but they don't have the same slope
        assert_eq!(intersection(0,0, 10,10,  1,2, 1,5), Intersection::None);
    }

    #[test]
    fn intersect8() {
        assert_eq!(intersection(1,2, 1,1,  1,3, 1,2), Intersection::EndToEnd);
        assert_eq!(intersection(1,1, 1,2,  1,3, 1,2), Intersection::EndToEnd);

        assert_eq!(intersection(1,2, 1,1,  1,2, 1,3), Intersection::EndToEnd);
        assert_eq!(intersection(1,1, 1,2,  1,2, 1,3), Intersection::EndToEnd);
    }

    #[test]
    fn validity_checks() {
        let geom: LineString<i32> = LineString(vec![]);
        assert!(!is_linestring_valid(&geom));

        assert!(!is_linestring_valid(&LineString(vec![(0i32, 0i32).into()])));

        // Linestrings can self-intersect
        assert!(is_linestring_valid(&vec![(0, 0), (4, 0), (2, -1), (2, 1)].into()));

        assert!(has_self_intersections(&vec![(0, 0), (4, 0), (2, -1), (2, 1), (0,0)].into()));
        assert!(has_self_intersections(&vec![(0, 0), (4, 0), (2, -1), (2, 0), (2, 1), (0,0)].into()));

        // Simple square - valid
        assert!(is_polygon_valid(&Polygon::new(vec![(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)].into(), vec![])));
        // Unclosed - invalid
        assert!(!is_polygon_valid(&Polygon::new(vec![(0, 0), (0, 1), (1, 1), (1, 0)].into(), vec![])));

        // Has a touching inner
        let geom: Polygon<i32> = Polygon::new(vec![(0, 0), (0, 2), (1, 2), (1, 1), (2, 1), (2, 3), (1, 3), (1, 2), (0, 2), (0, 4), (3, 4), (3, 0), (0, 0)].into(), vec![]);
        assert!(!is_polygon_valid(&geom));
    }
    
    #[test]
    fn make_valid1() {
        let unit_square = vec![(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)];
        let geom: Polygon<i32> = Polygon::new(unit_square.clone().into(), vec![]);
        
        let mut new_geom = make_polygon_valid(geom);
        assert_eq!(new_geom.0.len(), 1);
        let new_geom: Polygon<_> = new_geom.0.remove(0);
        assert!(is_polygon_valid(&new_geom));
        assert_eq!(new_geom.exterior, unit_square.into());
    }

    #[test]
    fn make_valid2() {
        // a-----b
        // | g-h |
        // e-f | |
        // | j-i |
        // d-----c
        let a = Point::new(0, 0); let b = Point::new(6, 0);
        let c = Point::new(6, 4); let d = Point::new(0, 4);
        let e = Point::new(0, 2); let f = Point::new(2, 2);
        let g = Point::new(2, 1); let h = Point::new(4, 1);
        let i = Point::new(4, 3); let j = Point::new(2, 3);

        let geom = Polygon::new(vec![a, b, c, d, e, f, j, i, h, g, f, e, a].into(), vec![]);
        assert!(!is_polygon_valid(&geom));
        
        let mut new_geom = make_polygon_valid(geom);
        assert_eq!(new_geom.0.len(), 1);
        let new_geom: Polygon<_> = new_geom.0.remove(0);
        assert!(is_polygon_valid(&new_geom));
        assert_eq!(new_geom, Polygon::new(vec![a, e, d, c, b, a].into(), vec![vec![f, g, h, i, j, f].into()]));
    }

    #[test]
    fn make_valid3() {
        // a-----b
        // | g-h |
        // | | | |
        // | j-i |
        // d-----c

        let a = Point::new(0, 0); let b = Point::new(6, 0);
        let c = Point::new(6, 4); let d = Point::new(0, 4);
        let g = Point::new(2, 1); let h = Point::new(4, 1);
        let i = Point::new(4, 3); let j = Point::new(2, 3);
        
        let p = Polygon::new(vec![a, d, c, b, a].into(), vec![vec![g, h, i, j, g].into()]);
        assert!(is_polygon_valid(&p));
        let original = p.clone();

        let mut p: MultiPolygon<_> = make_polygon_valid(p);
        assert_eq!(p.0.len(), 1);
        let p: Polygon<_> = p.0.remove(0);
        assert!(is_polygon_valid(&p));
        assert_eq!(p, original);
    }

    #[test]
    fn make_valid4() {
        // a-----b
        // | g-h |
        // | | | |
        // | j-i |
        // d-----c

        let a = Point::new(0, 0); let b = Point::new(6, 0);
        let c = Point::new(6, 4); let d = Point::new(0, 4);
        let g = Point::new(2, 1); let h = Point::new(4, 1);
        let i = Point::new(4, 3); let j = Point::new(2, 3);
        

        // Same but the inner 
        let p_outer = Polygon::new(vec![a, d, c, b, a].into(), vec![]);
        assert!(is_polygon_valid(&p_outer));
        let p_inner = Polygon::new(vec![g, j, i, h, g].into(), vec![]);
        assert!(is_polygon_valid(&p_inner));
        let mp = MultiPolygon(vec![p_outer.clone(), p_inner.clone()]);

        let mut new_mp = match make_valid(mp.into()) {
            Geometry::MultiPolygon(x) => x,
            _ => unreachable!(),
        };

        //println!("{:?}", new_mp);
        assert_eq!(new_mp.0.len(), 1);
        let poly = new_mp.0.remove(0);
        assert_eq!(poly.exterior, vec![a, d, c, b, a].into());
        assert_eq!(poly.interiors.len(), 1);
        assert_eq!(poly.interiors[0], vec![g, h, i, j, g].into());

    }

    #[test]
    fn make_valid5() {
        // This polygon touches at a point (d). it should be 2 polygons
        //   a-b
        //   | |
        // g-d-c
        // | |
        // f-e
        let a = Point::new(2, 0); let b = Point::new(4, 0); let c = Point::new(4, 6);
        let d = Point::new(2, 4);
        let e = Point::new(2, 6); let f = Point::new(0, 6); let g = Point::new(0, 4);
        // sanity check
        assert!(is_polygon_valid(&Polygon::new(vec![a, d, c, b, a].into(), vec![])));
        assert!(is_polygon_valid(&Polygon::new(vec![d, g, f, e, d].into(), vec![])));

        let poly = Polygon::new(vec![a, d, g, f, e, d, c, b, a].into(), vec![]);
        //assert!(!is_polygon_valid(&poly));

        let new_mp: MultiPolygon<_> = make_polygon_valid(poly);

        assert_eq!(new_mp.0.len(), 2);
        assert_eq!(new_mp.0[0], Polygon::new(vec![d, g, f, e, d].into(), vec![]));
        assert!(is_polygon_valid(&new_mp.0[0]));
        assert_eq!(new_mp.0[1], Polygon::new(vec![a, d, c, b, a].into(), vec![]));
        assert!(is_polygon_valid(&new_mp.0[1]));

    }

    // Helper function that tests that applying func to in_obj doesn't result in in_obj changing
    fn test_no_change<T, F>(func: F, mut in_obj: T)
        where F: Fn(&mut T), T: Clone+Debug+PartialEq
    {
        let out_obj = in_obj.clone();
        expected_results(func, in_obj, out_obj);
    }

    fn expected_results<T, F>(func: F, mut in_obj: T, out_obj: T)
        where F: Fn(&mut T), T: Debug+PartialEq
    {
        func(&mut in_obj);
        assert_eq!(in_obj, out_obj);
    }

    fn test_no_change_own_vec<T, F>(func: F, mut in_obj: T)
        where F: Fn(T)->Vec<T>, T: Clone+Debug+PartialEq
    {
        let expected_out = vec![in_obj.clone()];
        let out = func(in_obj);
        assert_eq!(out, expected_out);
    }


    #[test]
    fn test_add_points_for_all_crossings1() {
        test_no_change(add_points_for_all_crossings, LineString(vec![(0i32, 0i32).into()]));
        test_no_change(add_points_for_all_crossings, vec![(0, 0), (4, 0), (2, -1)].into());
        test_no_change(add_points_for_all_crossings, vec![(0, 0), (2, 0), (4, 0), (2, -1), (2, 0), (2, 1), (0,0)].into());

        expected_results(add_points_for_all_crossings, vec![(0, 0), (4, 0), (2, -1), (2, 0), (2, 1), (0,0)].into(), vec![(0, 0), (2, 0), (4, 0), (2, -1), (2, 0), (2, 1), (0,0)].into());
        expected_results(add_points_for_all_crossings, vec![(0, 0), (4, 0), (2, -1), (2, 1)].into(), vec![(0, 0), (2, 0), (4, 0), (2, -1), (2, 0), (2, 1)].into());
    }

    #[test]
    fn test_add_points_for_all_crossings2() {
        expected_results(add_points_for_all_crossings, vec![(0, 0), (10, 0), (5, 0), (5, 10), (0, 0)].into(), vec![(0, 0), (5, 0), (10, 0), (5, 0), (5, 10), (0, 0)].into());
    }
    #[test]
    fn test_add_points_for_all_crossings3() {
        expected_results(add_points_for_all_crossings, vec![(0, 0), (10, 0), (-2, 0), (-2, 10), (0, 0)].into(), vec![(0, 0), (10, 0), (0, 0), (-2, 0), (-2, 10), (0, 0)].into());
    }
    #[test]
    fn test_add_points_for_all_crossings4() {
        expected_results(add_points_for_all_crossings,
                         vec![(0, 0), (100, 0), (100, 100), (70, 0), (50, 0), (0, 100), (0, 0)].into(),
                         vec![(0, 0), (50, 0), (70, 0), (100, 0), (100, 100), (70, 0), (50, 0), (0, 100), (0, 0)].into() );
    }
    #[test]
    fn test_add_points_for_all_crossings5() {
        expected_results(add_points_for_all_crossings,
                         vec![(0, 0), (100, 0), (110, 100), (110, 0), (50, 0), (0, 100), (0, 0)].into(),
                         vec![(0, 0), (50, 0), (100, 0), (110, 100), (110, 0), (100, 0), (50, 0), (0, 100), (0, 0)].into() );
    }

    #[test]
    fn test_dissolve_into_rings1() {
        test_no_change_own_vec(dissolve_into_rings, vec![(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)].into());

        // This polygon touches at a point (d). it should be 2 polygons
        //   a-b
        //   | |
        // g-d-c
        // | |
        // f-e
        let a = Point::new(2, 0); let b = Point::new(4, 0); let c = Point::new(4, 6);
        let d = Point::new(2, 4);
        let e = Point::new(2, 6); let f = Point::new(0, 6); let g = Point::new(0, 4);
        // sanity check
        assert!(is_polygon_valid(&Polygon::new(vec![a, d, e, b, a].into(), vec![])));
        assert!(is_polygon_valid(&Polygon::new(vec![d, g, f, e, d].into(), vec![])));

        let ls = vec![a, d, g, f, e, d, c, b, a].into();

        let result = dissolve_into_rings(ls);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], vec![d, g, f, e, d].into());
        assert_eq!(result[1], vec![a, d, c, b, a].into());
        
        assert_eq!(result, vec![
                   vec![d, g, f, e, d].into(),
                   vec![a, d, c, b, a].into(),
                   ] );

    }

    #[test]
    fn test_dissolve_into_rings2() {
        let a = Point::new(0, 0); let b = Point::new(2, 0); let c = Point::new(3, 0);
        let d = Point::new(1, 1);

        // a---b--c
        //  \ /
        //   d
        
        // a-b-a Just the spike
        assert_eq!(dissolve_into_rings(LineString(vec![a, b, a])), vec![]);

        // Triangle (a-b-d-a) is kept, the little spike (b-c-b) is removed.
        assert_eq!(dissolve_into_rings(LineString(vec![a, b, c, b, d, a])), vec![ vec![a, b, d, a].into(), ]);

    }

    #[test]
    fn test_dissolve_into_rings3() {
        let a = Point::new(0, 0); let c = Point::new(2, 0);
        let b = Point::new(1, 1); let d = Point::new(2, 1);
        let e = Point::new(1, 2); let f = Point::new(2, 2);

        // a----c
        //  \ / |
        //   b--d
        //   |  |
        //   e--f
        // Triangle abc is filled in, bcd isn't. cdef is a square that's filled in.
        // Ideally we want 2 rings abdefca and bcdb

        let result = dissolve_into_rings(LineString(vec![a, b, c, d, b, e, f, d, c, a]));
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], vec![b, c, d, b].into());
        assert_eq!(result[1], vec![a, b, e, f, d, c, a].into());
    }

    #[test]
    fn test_dissolve_into_rings4() {
        // a-----b
        // | g-h |
        // e-f | |
        // | j-i |
        // d-----c

        let a = Point::new(0, 0); let b = Point::new(6, 0);
        let c = Point::new(6, 4); let d = Point::new(0, 4);
        let e = Point::new(0, 2); let f = Point::new(2, 2);
        let g = Point::new(2, 1); let h = Point::new(4, 1);
        let i = Point::new(4, 3); let j = Point::new(2, 3);


        let result = dissolve_into_rings(LineString(vec![a, b, c, d, e, f, g, h, i, j, f, e, a]));
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], vec![f, g, h, i, j, f].into());
        assert_eq!(result[1], vec![a, b, c, d, e, a].into());
    }

    #[test]
    fn convert_rings_to_polygons1() {
        assert_eq!(convert_rings_to_polygons(Vec::<LineString<i32>>::new()), MultiPolygon(vec![]));

        let unit_square = vec![(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)];
        assert_eq!(convert_rings_to_polygons(vec![unit_square.clone().into()]), MultiPolygon(vec![Polygon::new(unit_square.into(), vec![])]));
    }

    #[test]
    fn convert_rings_to_polygons2() {
        // a-----b
        // | g-h |
        // e f | |
        // | j-i |
        // d-----c
        let a = Point::new(0, 0); let b = Point::new(6, 0);
        let c = Point::new(6, 4); let d = Point::new(0, 4);
        let e = Point::new(0, 2); let f = Point::new(2, 2);
        let g = Point::new(2, 1); let h = Point::new(4, 1);
        let i = Point::new(4, 3); let j = Point::new(2, 3);


        let outer: LineString<_> = vec![a, b, c, d, e, a].into();
        let inner: LineString<_> = vec![g, h, i, j, f, g].into();
        let rings = vec![ outer.clone(), inner.clone() ];

        assert_eq!(convert_rings_to_polygons(rings), MultiPolygon(vec![Polygon::new(outer, vec![inner])]));
    }

    #[test]
    fn does_ray_cross1() {
        fn know_answer((x1, y1): (i32, i32), (x2, y2): (i32, i32), res: Crossing) {
            assert_eq!(does_ray_cross(&(0,0).into(), &(x1, y1).into(), &(x2, y2).into()), res, "({:?}, {:?}), ({:?}, {:?}) {:?}", x1, y1, x2, y2, res);
        }

        know_answer((1, 1), (10, 10), Crossing::No);
        know_answer((1, 0), (2, 0), Crossing::No);
        know_answer((-10, 10), (-10, 20), Crossing::No);
        know_answer((-10, -10), (-10, -20), Crossing::No);

        know_answer((0, 0), (10, 10), Crossing::Touches);
        know_answer((10, 1), (0, 0), Crossing::Touches);
        know_answer((-10, 0), (-5, 0), Crossing::Touches);

        know_answer((-10, 10), (-10, -10), Crossing::Yes);
        know_answer((-10, 10), (-10, -10), Crossing::Yes);

    }

    #[test]
    fn does_ray_cross2() {
        assert_eq!(does_ray_cross(&(1,2).into(), &(0, 0).into(), &(0, 2).into()), Crossing::EndPoint);
        assert_eq!(does_ray_cross(&(1,2).into(), &(0, 2).into(), &(0, 0).into()), Crossing::StartPoint);
    }

    #[test]
    fn does_ray_cross3() {
        assert_eq!(does_ray_cross(&(1,2).into(), &(0, 0).into(), &(0, 2).into()), Crossing::EndPoint);
        assert_eq!(does_ray_cross(&(1,2).into(), &(0, 2).into(), &(0, 4).into()), Crossing::StartPoint);
    }

    #[test]
    fn calc_rings_ext_int1() {
        // a-----b
        // | g-h |
        // e-f | |
        // | j-i |
        // d-----c

        let a = Point::new(0, 0); let b = Point::new(6, 0);
        let c = Point::new(6, 4); let d = Point::new(0, 4);
        let e = Point::new(0, 2); let f = Point::new(2, 2);
        let g = Point::new(2, 1); let h = Point::new(4, 1);
        let i = Point::new(4, 3); let j = Point::new(2, 3);
        
        let unit_square: LineString<_> = vec![a, b, c, d, a].into();
        let result = calc_rings_ext_int(vec![unit_square.clone()]);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, unit_square);
        assert_eq!(result[0].1, RingType::Exterior);

        let inner_square: LineString<_> = vec![g, h, i, j, g].into();
        let result = calc_rings_ext_int(vec![unit_square.clone(), inner_square.clone()]);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].0, unit_square);
        assert_eq!(result[0].1, RingType::Exterior);
        assert_eq!(result[1].0, inner_square);
        assert_eq!(result[1].1, RingType::Interior);

        // same but with other order
        let result = calc_rings_ext_int(vec![inner_square.clone(), unit_square.clone()]);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].0, inner_square);
        assert_eq!(result[0].1, RingType::Interior);
        assert_eq!(result[1].0, unit_square);
        assert_eq!(result[1].1, RingType::Exterior);

    }
    #[test]
    fn calc_rings_ext_int2() {
        let ring1: LineString<_> = vec![(1, 2), (1, 1), (2, 1), (2, 3), (1, 3), (1, 2)].into();
        let ring2: LineString<_> = vec![(0, 0), (0, 2), (0, 4), (3, 4), (3, 0), (0, 0)].into();

        let result = calc_rings_ext_int(vec![ring1.clone(), ring2.clone()]);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].0, ring1);
        assert_eq!(result[0].1, RingType::Interior);
        assert_eq!(result[1].0, ring2);
        assert_eq!(result[1].1, RingType::Exterior);
    }

    #[test]
    fn order_points1() {
        assert_eq!(order_points( ((0,0), (10, 0)), (5, 0), (1, 0) ),  Ordering::Greater );
        assert_eq!(order_points( ((0,0), (10, 0)), (1, 0), (5, 0) ),  Ordering::Less );
        assert_eq!(order_points( ((10,0), (0, 0)), (1, 0), (5, 0) ),  Ordering::Greater );
        assert_eq!(order_points( ((10,0), (0, 0)), (5, 0), (1, 0) ),  Ordering::Less );

        assert_eq!(order_points( ((0,0), (10, 0)), (0, 0), (10, 0) ), Ordering::Less );
        assert_eq!(order_points( ((0,0), (10, 0)), (10, 0), (0, 0) ), Ordering::Greater );

        assert_eq!(order_points( ((0,0), (10, 0)), (0, 0), (5, 0) ), Ordering::Less  );
        assert_eq!(order_points( ((0,0), (10, 0)), (5, 0), (0, 0) ), Ordering::Greater );

        assert_eq!(order_points( ((0,0), (10, 0)), (5, 0), (10, 0) ), Ordering::Less );
        assert_eq!(order_points( ((0,0), (10, 0)), (10, 0), (5, 0) ), Ordering::Greater );
    }

}

