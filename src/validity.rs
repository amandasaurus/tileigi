use geo::*;
use geo::map_coords::MapCoords;
use geo::intersects::Intersects;
use std::cmp::{min, max, Ord, Ordering};
use std::ops::{Add, Sub, DivAssign,Rem,Mul,AddAssign};
use std::collections::HashMap;
use num_traits::Signed;
use std::fmt::Debug;
use std::hash::Hash;

use ::simplify;

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
            // Y goes positive down, ergo the winding order is 'wrong way around' since the winding
            // order code works with y up
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

    for (i, points12) in ls.0.windows(2).enumerate() {
        let (p1, p2) = (points12[0], points12[1]);
        
        for points34 in ls.0[i+1..].windows(2).take(ls.0.len()-i-1) {

            if max(p1.x(), p2.x()) < min(points34[0].x(), points34[1].x()) || min(p1.x(), p2.x()) > max(points34[0].x(), points34[1].x())
                || max(p1.y(), p2.y()) < min(points34[0].y(), points34[1].y()) || min(p1.y(), p2.y()) > max(points34[0].y(), points34[1].y())
            {
                continue;
            }
            // For some reason it's a little faster to do this here after the check
            let (p3, p4) = (points34[0], points34[1]);

            match intersection(p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y(), p4.x(), p4.y()) {
                Intersection::Crossing(_) | Intersection::Overlapping(_, _)  => { return true; },
                Intersection::Touching(_) => { return true; },
                Intersection::None | Intersection::EndToEnd => {},
            }
        }
    }


    false
}

fn in_bounds<U: Ord+Copy>(z: U, a: U, b: U) -> bool {
    z >= min(a,b) && z <= max(a,b)
}

/// True iff point p is collinear with the line ab (NB: true if it goes beyond the ends)
fn collinear<T: CoordinateType>(a: (T, T), b: (T, T), p: (T, T)) -> bool {
    // (x2 - x1)(y - y1) == (y2 - y1)(x - x1)
    (b.0 - a.0)*(p.1 - a.1) == (b.1 - a.1)*(p.0 - a.0)
}

/// True iff p lies on the line segment ab, i.e. between the two, incl a and not b
/// Assumes that p is already collinear with ab
fn point_on_line_incl_end<T: CoordinateType+Ord+Copy>(a: (T, T), b: (T, T), p: (T, T)) -> bool {
    debug_assert!(collinear(a, b, p));
    in_bounds(p.0, a.0, b.0) && in_bounds(p.1, a.1, b.1)
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

    assert!((x1, y1) != (x2, y2), "(x1, y2) == (x2, y2) == {:?}", (x1, y1));
    assert!((x3, y3) != (x4, y4), "(x3, y3) == (x4, y4) == {:?}", (x3, y3));

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
        

        let p1_collinear_34 = collinear((x3, y3), (x4, y4), (x1, y1));
        let p2_collinear_34 = collinear((x3, y3), (x4, y4), (x2, y2));

        /// True iff p lies on the line segment ab, i.e. between the two, and is not a and not b
        /// (i.e. is on the line, but is not at the end points)
        /// Assumes that p is already collinear with ab
        fn point_on_line<T: CoordinateType+Ord+Copy>(a: (T, T), b: (T, T), p: (T, T)) -> bool {
            debug_assert!(collinear(a, b, p));
            (p != a) && (p != b) && in_bounds(p.0, a.0, b.0) && in_bounds(p.1, a.1, b.1)
        }


        fn delta<T: Ord+Sub<Output=T>>(a: T, b: T) -> T {
            if a > b { a - b } else { b - a }
        }

        match (p1_collinear_34, p2_collinear_34) {
            (false, false) => {
                // neither points are collinear, so no match
                return Intersection::None;
            },

            (true, false) | (false, true) => {
                // One point is on the line and the other is? But the slopes are the same, so this
                // should be impossible
                unreachable!();
            },

            (true, true) => {
                // These lines totally overlap
                let delta_x = delta(x1, x2)+delta(x3, x4);
                let delta_y = delta(y1, y2)+delta(y3, y4);
                if    (delta_x == delta(x1, x4) && delta_y == delta(y1, y4)) // 1-23-4
                   || (delta_x == delta(x2, x4) && delta_y == delta(y2, y4)) // 2-13-4
                   || (delta_x == delta(x1, x3) && delta_y == delta(y1, y3)) // 1-24-3
                   || (delta_x == delta(x2, x3) && delta_y == delta(y2, y3)) // 2-13-4
                    {
                    // One after the other. We know they have the same slope, so this shortcut
                    // calculation works.
                    return Intersection::EndToEnd;
                }
                let p3_on_12 = point_on_line((x1, y1), (x2, y2), (x3, y3));
                let p4_on_12 = point_on_line((x1, y1), (x2, y2), (x4, y4));
                match (p3_on_12, p4_on_12) {
                    (true, true) => {
                        // both on the line
                        return Intersection::Overlapping((x3, y3), (x4, y4));
                    },
                    (true, false) => {
                        // p3 is on the line 12, but which of p1 & p2 is the other point
                        // either p1 or p2 is on the line 34
                        debug_assert!(point_on_line_incl_end((x3, y3), (x4, y4), (x1, y1)) || point_on_line_incl_end((x3, y3), (x4, y4), (x2, y2)));
                        let other_point = if point_on_line_incl_end((x3, y3), (x4, y4), (x1, y1)) {
                            (x1, y1)
                        } else {
                            debug_assert!(point_on_line_incl_end((x3, y3), (x4, y4), (x2, y2)));
                            (x2, y2)
                        };
                        return Intersection::Overlapping((x3, y3), other_point);
                    },
                    (false, true) => {
                        // p4 is on the line 12, but which of p1 & p2 is the other point
                        // either p1 or p2 is on the line 34
                        debug_assert!(point_on_line_incl_end((x3, y3), (x4, y4), (x1, y1)) || point_on_line_incl_end((x3, y3), (x4, y4), (x2, y2)));
                        let other_point = if point_on_line_incl_end((x3, y3), (x4, y4), (x1, y1)) {
                            (x1, y1)
                        } else {
                            debug_assert!(point_on_line_incl_end((x3, y3), (x4, y4), (x2, y2)));
                            (x2, y2)
                        };
                        return Intersection::Overlapping((x4, y4), other_point);
                    },
                    (false, false) => {
                        // This can happen when 12 is a subset of 34
                        debug_assert!(point_on_line_incl_end((x3, y3), (x4, y4), (x1, y1)) && point_on_line_incl_end((x3, y3), (x4, y4), (x2, y2)));
                        return Intersection::Overlapping((x1, y1), (x2, y2));
                    }
                }
            },

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

pub fn make_valid<T: CoordinateType+Debug+Ord+Signed+Hash+DivAssign+Rem+AddAssign+Into<f64>>(mut geom: Geometry<T>) -> Geometry<T> {
    //println!("{} L {}", file!(), line!());
    //simplify::remove_unneeded_points(&mut geom);
    if is_valid(&geom) {
        return geom;
    }

    //println!("{} L {}", file!(), line!());
    match geom {
        Geometry::Polygon(p) => Geometry::MultiPolygon(make_polygon_valid(p)),
        Geometry::MultiPolygon(mp) => Geometry::MultiPolygon(make_multipolygon_valid(mp)),
        x => x,
    }
}

fn make_multipolygon_valid<T: CoordinateType+Debug+Ord+Signed+Hash+Into<f64>>(mut mp: MultiPolygon<T>) -> MultiPolygon<T> {
    //println!("{} L {}", file!(), line!());
    let MultiPolygon( polygons ) = mp;

    let rings: Vec<LineString<T>> = polygons.into_iter().flat_map(|p| {
        let Polygon{ exterior, interiors } = p;
        let mut these_rings = interiors;
        these_rings.insert(0, exterior);
        these_rings.into_iter()
    }).collect();


    make_rings_valid(rings)
}

fn make_polygon_valid<T: CoordinateType+Debug+Ord+Signed+Hash+Into<f64>>(mut p: Polygon<T>) -> MultiPolygon<T> {
    //println!("{} L {}", file!(), line!());
    let Polygon{ exterior, interiors } = p;
    let mut rings = interiors;
    rings.insert(0, exterior);

    make_rings_valid(rings)

}

fn make_rings_valid<T: CoordinateType+Debug+Ord+Signed+Hash+Into<f64>>(mut rings: Vec<LineString<T>>) -> MultiPolygon<T> {
    //println!("{} L {}", file!(), line!());

    let mut new_rings: Vec<LineString<T>> = Vec::with_capacity(rings.len());
    for mut ring in rings.into_iter() {
        //println!("{} L {}", file!(), line!());
        add_points_for_all_crossings(&mut ring);
        //println!("{} L {}", file!(), line!());
        let mut these_rings = dissolve_into_rings(ring);
        new_rings.append(&mut these_rings);
    }
    //println!("{} L {}", file!(), line!());

    let rings = new_rings;
    //println!("{} L {}", file!(), line!());
    
    let result = convert_rings_to_polygons(rings);

    // This takes a geom, so we do a dance
    let mut result = Geometry::MultiPolygon(result);
    ensure_polygon_orientation(&mut result);

    //println!("{} L {}", file!(), line!());
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
    //println!("{}:{}", file!(), line!());

    loop {
        //println!("{}:{} ls.0.len() {}", file!(), line!(), ls.0.len());
        //println!("\nStart of loop\n{:?}", ls.0);
        let mut coords_to_insert = HashMap::new();
        // Keys are the point indexes.
        // Values are a Vec of new points to add after the point with that index.
        // So vec![(0, 0), (1, 0)] for key #3, means to insert those 2 points after ls.0[3]
        // They are initially stored in the order they appear in, but they need to be sorted
        // afterwards

        for (i, points12) in ls.0.windows(2).enumerate() {
            //println!("{} L {}", file!(), line!());
            let (p1, p2) = (points12[0], points12[1]);
            //if bad {
            //    println!("{}:{} p1 {:?} p2 {:?}", file!(), line!(), p1, p2);
            //}
            
            for (j, points34) in ls.0[i+1..].windows(2).enumerate().take(ls.0.len()-i-1) {
                let j = j + i + 1;
                let (p3, p4) = (points34[0], points34[1]);
                let x1 = p1.x(); let y1 = p1.y();
                let x2 = p2.x(); let y2 = p2.y();
                let x3 = p3.x(); let y3 = p3.y();
                let x4 = p4.x(); let y4 = p4.y();
                //println!("looking at i {} j {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", i, j, p1, p2, p3, p4);

                if max(x1, x2) < min(x3, x4) || min(x1, x2) > max(x3, x4)
                    || max(y1, y2) < min(y3, y4) || min(y1, y2) > max(y3, y4)
                {
                    continue;
                }

                //if bad {
                //    println!("{}:{} i {} j {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", file!(), line!(), i, j, p1, p2, p3, p4);
                //}

                //println!("{} L {}", file!(), line!());
                match intersection(x1, y1, x2, y2, x3, y3, x4, y4) {
                    Intersection::None | Intersection::EndToEnd => {},

                    Intersection::Crossing(crosspoint) => {
                        //println!("looking at i {} j {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", i, j, p1, p2, p3, p4);
                        //println!("i {} j {} crossing {:?}", i, j, crosspoint);
                        // A "unit square" can cause a crossing. ie. (0,0)-(1,1) and (1,0)-(0,1)
                        // (diagonal). That's returned as Crossing((1, 1)).
                        // So don't add a point if it would cause a duplicate
                        // We basically never want 2 identical points, one after the other

                        // In cases of a diagonol crossing, the 3 points won't be collinear.
                        //debug_assert!(collinear((x1, y1), (x2, y2), crosspoint), "L {} !collinear {:?} {:?} - {:?} {:?} point {:?}", line!(), (x1,y1), (x2, y2), (x3, y3), (x4, y4), crosspoint);
                        //debug_assert!(point_on_line_incl_end((x1, y1), (x2, y2), crosspoint));

                        if (x1, y1) != crosspoint && (x2, y2) != crosspoint {
                            coords_to_insert.entry(i).or_insert(vec![]).push(crosspoint);
                        }
                        if (x3, y3) != crosspoint && (x4, y4) != crosspoint {
                            coords_to_insert.entry(j).or_insert(vec![]).push(crosspoint);
                        }
                    },

                    Intersection::Overlapping(overlap1, overlap2)  => {
                        //println!("looking at i {} j {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", i, j, p1, p2, p3, p4);
                        //println!("i {} j {} overlapping {:?},{:?}", i, j, overlap1, overlap2);
                        debug_assert!(overlap1 != overlap2);
                        //debug_assert!(collinear((x1, y1), (x2, y2), overlap1));
                        //debug_assert!(point_on_line_incl_end((x1, y1), (x2, y2), overlap1));
                        //debug_assert!(collinear((x1, y1), (x2, y2), overlap2));
                        //debug_assert!(point_on_line_incl_end((x1, y1), (x2, y2), overlap2));

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
                        debug_assert!(collinear((x1, y1), (x2, y2), (x0, y0)));
                        debug_assert!(point_on_line_incl_end((x1, y1), (x2, y2), (x0, y0)));
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


        //println!("{}:{}", file!(), line!());
        if coords_to_insert.is_empty() {
            break;
        } else {
            //println!("{}:{}", file!(), line!());
            // When we insert a point into the vec, it'll push all after that along. Keep track of
            // how many we've inserted, so we know the correct place to push the later ones
            let mut offset = 0;

            // Turn hashmap into a sorted vec, sorted by index to add
            let coords_to_insert = ls.0.windows(2).enumerate().filter_map(|(idx, points)| {
                let (point1, point2) = (points[0], points[1]);
                if let Some(mut new_points) = coords_to_insert.remove(&idx) {
                    //println!("index {:?} point1 {:?} point2 {:?} new_points {:?}", idx, point1, point2, new_points);
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
            //println!("{} L {}", file!(), line!());

            //println!("We added {} new points to the line", offset);
        }
    }
    //println!("{}:{} finished", file!(), line!());

}

fn dissolve_into_rings<T: CoordinateType+Debug+Hash+Eq+Into<f64>>(ls: LineString<T>) -> Vec<LineString<T>> {
    let LineString( points ) = ls;
    if points.len() <= 3 {
        // Not enough points for a proper ring
        return vec![];
    }

    // Key: a point (x,y) values are a vec of usizes representing the index (in the linestring)
    // where this point is the first point of that segment. e.g. (0, 0): [1, 5], means that the
    // line segment [1, 2] starts at point (0,0), i.e. points[1] == (0,0), likewise for segment [5,
    // 6]
    let mut outgoing_segments = HashMap::with_capacity(points.len());

    for (i, p) in points.iter().enumerate() {
        // TODO here we could assert that the existing vec is <=2, and generate the loops vec here,
        // rather than do a loop later
        outgoing_segments.entry((p.x(), p.y())).or_insert(vec![]).push(i);
    }
    //println!("{}:{} outgoing_segments {:?}", file!(), line!(), outgoing_segments);
    //println!("{}:{} outgoing_segments w/ > 1 segment: {:?}", file!(), line!(), outgoing_segments.iter().filter_map(|(k, v)| if v.len() > 1 { Some((k, v)) } else { None }).collect::<Vec<_>>());

    // loops: a Vec of Vec's. Each inner vec is 2 point indexes, and means 'there is a loop from
    // this point to that point'
    let mut loops: Vec<Vec<usize>> = outgoing_segments.into_iter().filter_map(|(_, v)| if v.len() > 1 { Some(v) } else { None }).collect();
    //println!("{}:{} ls {:?}", file!(), line!(), points);
    //println!("{}:{} loops {:?}", file!(), line!(), loops);

    // If the linestring touches the first/last point, then this algoritm will result in a loop
    // with 3 numbers, [0, touching point idx, last]. Make this into 2 loops
    //
    // This is a list of indices in loops where these problems occur
    let mut loop_with_3_points = loops.iter().enumerate()
            .filter_map(|(i, l)| if l.len() == 3 && l[0] == 0 && l[2] == points.len()-1 { Some(i) } else { None })
            .collect::<Vec<usize>>();
    // 0 = nothing is strange, 1 = simple case of first/last point being 'touched'. >1 who knows?
    assert!(loop_with_3_points.len() == 0 || loop_with_3_points.len() == 1);
    if loop_with_3_points.len() == 1 {
        let single_loop = loops.swap_remove(loop_with_3_points.remove(0));
        loops.push(vec![single_loop[0], single_loop[1]]);
        loops.push(vec![single_loop[1], single_loop[2]]);
    }

    if loops.len() == 1 {
        if loops[0].len() == 2 && loops[0][0] == 0 && loops[0][1] == points.len()-1 {
            // start & end the same, ergo only one loop here
            return vec![LineString(points)];
        } else {
            //return Vec::new();
            // FIXME do something here
            // There is only one loop, and it is not a simple outer loop.
            //eprintln!("outgoing_segments {:?}", outgoing_segments);
            eprintln!("points {:?}", points);
            eprintln!("loops {:?}", loops);
            for (i, p) in points.iter().enumerate() {
                eprintln!("{:03} {:?},{:?}", i, p.x(), p.y());
            }
            unreachable!();
        }
    }

    let mut point_unassigned = vec![true; points.len()];
    let mut results: Vec<LineString<T>> = vec![];


    debug_assert!(loops.iter().all(|idxes| idxes.len() == 2), "There is a point with != 2 segments: {:?}\npoints: {:?}", loops, points);
    if loops.iter().any(|idxes| idxes.len() != 2) {
        //eprintln!("{} {} loops {:?}", file!(), line!(), loops);
        // FIXME do something here
        return vec![];
    }
    // Remove weird loops.
    // This is not ideal, but it seems to make valid geometries for MbSC
    //let mut loops = loops.into_iter().filter(|idxes| idxes.len() == 2).collect::<Vec<_>>();

    //println!("points {:?}", points);
    //println!("outgoing_segments {:?}", outgoing_segments);
    // sort loops where the smaller length (in terms of number of points) are to the front.
    // Ideal: Sort them so that if a loop is a subset of a larger loop, then the smaller is ahead,
    // so the smaller, "inner" loop will be removed first. Unless something really strange is going
    // on, this sort-by-length should do it (since an outer loop will be longer than the inner one
    // it contains)
    loops.sort_by_key(|i| (i[1]-i[0], i[0]));
    //println!("{}:{} loops {:?}", file!(), line!(), loops);

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

    //println!("{}:{} point_unassigned {:?}", file!(), line!(), point_unassigned);
    // There will always be the last/first point unassigned since we keep the end around, which
    // means the endpoint of the outer ring is kept. So they should all be false, except the last
    // which is true
    debug_assert!(point_unassigned.iter().take(point_unassigned.len()-1).all(|x| !x), "{:?}", point_unassigned);
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

    /// The ray passes though the point entirely
    /// (i) The start or end point is the point
    /// (ii) it's a horizontal line and the ray passes along/through it
    /// (iii) The point is part of the line
    Touches,

    /// The ray goes through the start, or end, point of the line, but the point is not the
    /// start/end
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

    if ( y1 > y && y2 > y ) || ( y1 < y && y2 < y ) || (x1 > x && x2 > x ){
        // segment is entirely above, below, or to the right of, the point.
        return Crossing::No;
    } else if (x == x1 && y == y1)  // point is start point
       || (x == x2 && y == y2) // point is end point
       || ( (x2-x1)*(y - y1) == (x-x1)*(y2-y1) )  // point is on the line of a-b
       || ( y1 == y2 && y1 == y && ( x1 <= x || x2 <= x2 )  )  // the ray goes through all, or part of, the line segment
    {
        return Crossing::Touches;
    } else if x1 < x && y1 == y {
        return Crossing::StartPoint;
    } else if x2 < x && y2 == y {
        return Crossing::EndPoint;
    } else if (x1 < x || x2 < x) && ( (y1>y && y2<y) || (y1<y && y2>y) ) {
        return Crossing::Yes;
    } else {
        // I don't like this and would like to have all "No" cases explicity covered
        return Crossing::No;
    }

    //eprintln!("x {:?} y {:?} x1 {:?} y1 {:?} x2 {:?} y2 {:?}", x, y, x1, y1, x2, y2);
    //unreachable!();
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
        let point_x = point.x();
        let point_y = point.y();

        
        // loop over all the rings
        for (i, ring) in all_rings.iter().enumerate() {
            if i == ring_index { continue; }
            //println!("i {} point {:?}", i, point);

            // then all the segments in this ring
            for other_points in ring.0.windows(2) {
                debug_assert!(other_points.len() == 2);

                if ( other_points[0].y() > point_y && other_points[1].y() > point_y ) || ( other_points[0].y() < point_y && other_points[1].y() < point_y ) || (other_points[0].x() > point_x && other_points[1].x() > point_x ) {
                    // line is entirely above, below or to right of point.
                    continue;
                }

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
            RingType::Exterior => { exteriors.push(ring); },
            RingType::Interior => { interiors.push(ring); },
        }
    }
    assert!(!(exteriors.is_empty() && interiors.is_empty()));

    if exteriors.is_empty() {
        // FIXME implement this properly, esp if there are interiors
        return MultiPolygon(vec![]);
    }

    // this is a simple hack, take largest poly
    exteriors.sort_unstable_by_key(|p| -1 * p.0.len() as i64);
    let ring = exteriors.remove(0);
    //println!("ring {}", ring.0.len());
    return MultiPolygon(vec![Polygon::new(ring, vec![])]);

    
    // FIXME implement this
    // All interiors?!
    //assert!(!(exteriors.is_empty() && !interiors.is_empty()), "convert_rings_to_polygons {}\ninteriors {:?}\nexteriors {:?}\nno exteriors, but interiors", line!(), interiors, exteriors);

    //let mut polygons: Vec<_> = exteriors.into_iter().map(|p| Polygon::new(p, vec![])).collect();

    //// we need to calculate the what exterior that each interior is in
    //
    //if polygons.len() == 1 {
    //    // There is only one exterior ring, so take a simple approach of assuming all the
    //    // interiors are part of that
    //    ::std::mem::swap(&mut polygons[0].interiors, &mut interiors);
    //    
    //} else {
    //    if interiors.is_empty() {
    //        // nothing to do
    //    } else {
    //        // we need to figure out which exterior each interior is in.
    //        //eprintln!("polygons {:?}", polygons);
    //        //eprintln!("interiors {:?}", interiors);
    //        //unimplemented!()

    //        // Just skip it for now
    //        // FIXME implement this
    //    }
    //}


    //MultiPolygon(polygons)
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

    //debug_assert!(collinear(line.0, line.1, p1), "line {:?} p1 {:?}", line, p1);
    //debug_assert!(collinear(line.0, line.1, p2), "line {:?} p2 {:?}", line, p2);
    //debug_assert!(point_on_line_incl_end(line.0, line.1, p1));
    //debug_assert!(point_on_line_incl_end(line.0, line.1, p2));

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
    let slope_start_1 = sub(p1, line.0);

    // slope from the start point to p2
    let slope_start_2 = sub(p2, line.0);
    
    // slope from p1 to p2
    let slope_1_2 = sub(p2, p1);

    // slope from p2 to p1
    let slope_2_1 = sub(p1, p2);

    // slope from p2 to the end
    let slope_2_end = sub(line.1, p2);
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

        // Gonna presume they are equal, and if we do a stable sort then the order won't change
        // TODO Should this be a PartialOrd instead?
        Ordering::Equal

        //eprintln!("line {:?} p1 {:?} p2 {:?}", line, p1, p2);
        //eprintln!("slone_line {:?}", slope_line);
        //eprintln!("slope_start_1 {:?} slope_start_2 {:?}", slope_start_1, slope_start_2);
        //eprintln!("slope_1_2 {:?} slope_2_1 {:?}", slope_1_2, slope_2_1);
        //eprintln!("slope_2_end {:?}", slope_2_end);
        //unreachable!();
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
        assert_eq!(intersection(0, 0,  0, 10,  0, 0,  0, 5), Intersection::Overlapping((0, 5), (0, 0)));
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

    fn test_overlapping(p1: (i32, i32), p2: (i32, i32), p3: (i32, i32), p4: (i32, i32), res1: (i32, i32), res2: (i32, i32)) {
        let res = intersection(p1.0, p1.1, p2.0, p2.1, p3.0, p3.1, p4.0, p4.1);
        assert!(res == Intersection::Overlapping(res1, res2) || res == Intersection::Overlapping(res2, res1));

        let res = intersection(p2.0, p2.1, p1.0, p1.1, p3.0, p3.1, p4.0, p4.1);
        assert!(res == Intersection::Overlapping(res1, res2) || res == Intersection::Overlapping(res2, res1), "res {:?}", res);

        let res = intersection(p1.0, p1.1, p2.0, p2.1, p4.0, p4.1, p3.0, p3.1);
        assert!(res == Intersection::Overlapping(res1, res2) || res == Intersection::Overlapping(res2, res1));

        let res = intersection(p2.0, p2.1, p1.0, p1.1, p4.0, p4.1, p3.0, p3.1);
        assert!(res == Intersection::Overlapping(res1, res2) || res == Intersection::Overlapping(res2, res1));
    }

    #[test]
    fn intersect10() {
        test_overlapping((0,2), (0,0), (0,0), (0,1),   (0, 0), (0, 1));

        test_overlapping((2,0), (0,0), (0,0), (1,0),   (0, 0), (1, 0));

        test_overlapping((0,0), (5,0), (-5,0), (1,0),   (0, 0), (1, 0));

        test_overlapping((0,0), (0,5), (0,-5), (0,1),   (0, 0), (0, 1));

        test_overlapping((-10,-10), (10,10), (0,0), (5,5),   (0, 0), (5, 5));

        test_overlapping((0,0), (10,10), (0,0), (5,5),   (0, 0), (5, 5));
    }

    #[test]
    fn intersect11() {
        test_overlapping((0,0), (10,0), (10,0), (-2,0),   (0, 0), (10, 0));
    }

    #[test]
    fn intersect12() {
        assert_eq!(intersection(0,0, 1,1,  1,0, 0,1), Intersection::Crossing((1, 1)));
        assert_eq!(intersection(1,1, 0,0,  1,0, 0,1), Intersection::Crossing((1, 1)));
        assert_eq!(intersection(0,0, 1,1,  0,1, 1,0), Intersection::Crossing((1, 1)));
        assert_eq!(intersection(1,1, 0,0,  0,1, 1,0), Intersection::Crossing((1, 1)));

        assert_eq!(intersection(3,1, 4,0,  3,0, 4,1), Intersection::Crossing((4, 1)));
        assert_eq!(intersection(75,43, 76,42,  75,42, 76,43), Intersection::Crossing((76, 43)));
        assert_eq!(intersection(1975,1243, 1976,1242,  1975,1242, 1976,1243), Intersection::Crossing((1976, 1243)));
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
    fn add_points_for_all_crossings1() {
        test_no_change(add_points_for_all_crossings, LineString(vec![(0i32, 0i32).into()]));
        test_no_change(add_points_for_all_crossings, vec![(0, 0), (4, 0), (2, -1)].into());
        test_no_change(add_points_for_all_crossings, vec![(0, 0), (2, 0), (4, 0), (2, -1), (2, 0), (2, 1), (0,0)].into());

        expected_results(add_points_for_all_crossings, vec![(0, 0), (4, 0), (2, -1), (2, 0), (2, 1), (0,0)].into(), vec![(0, 0), (2, 0), (4, 0), (2, -1), (2, 0), (2, 1), (0,0)].into());
        expected_results(add_points_for_all_crossings, vec![(0, 0), (4, 0), (2, -1), (2, 1)].into(), vec![(0, 0), (2, 0), (4, 0), (2, -1), (2, 0), (2, 1)].into());
    }

    #[test]
    fn add_points_for_all_crossings2() {
        expected_results(add_points_for_all_crossings, vec![(0, 0), (10, 0), (5, 0), (5, 10), (0, 0)].into(), vec![(0, 0), (5, 0), (10, 0), (5, 0), (5, 10), (0, 0)].into());
    }
    #[test]
    fn add_points_for_all_crossings3() {
        expected_results(add_points_for_all_crossings, vec![(0, 0), (10, 0), (-2, 0), (-2, 10), (0, 0)].into(), vec![(0, 0), (10, 0), (0, 0), (-2, 0), (-2, 10), (0, 0)].into());
    }
    #[test]
    fn add_points_for_all_crossings4() {
        expected_results(add_points_for_all_crossings,
                         vec![(0, 0), (100, 0), (100, 100), (70, 0), (50, 0), (0, 100), (0, 0)].into(),
                         vec![(0, 0), (50, 0), (70, 0), (100, 0), (100, 100), (70, 0), (50, 0), (0, 100), (0, 0)].into() );
    }
    #[test]
    fn add_points_for_all_crossings5() {
        expected_results(add_points_for_all_crossings,
                         vec![(0, 0), (100, 0), (110, 100), (110, 0), (50, 0), (0, 100), (0, 0)].into(),
                         vec![(0, 0), (50, 0), (100, 0), (110, 100), (110, 0), (100, 0), (50, 0), (0, 100), (0, 0)].into() );
    }

    #[test]
    fn dissolve_into_rings1() {
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
    fn dissolve_into_rings2() {
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
    fn dissolve_into_rings3() {
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
    fn dissolve_into_rings4() {
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
    fn dissolve_into_rings5() {
        // This is a complicated real world example.
        let ls = LineString(vec![
            Point::new(31071, 21260),
            Point::new(31071, 21259),
            Point::new(31071, 21258),
            Point::new(31072, 21258),
            Point::new(31072, 21259),
            Point::new(31071, 21259),
            Point::new(31071, 21260),
            Point::new(31072, 21260),
            Point::new(31072, 21262),
            Point::new(31073, 21262),
            Point::new(31073, 21264),
            Point::new(31074, 21264),
            Point::new(31074, 21265),
            Point::new(31073, 21265),
            Point::new(31073, 21264),
            Point::new(31072, 21264),
            Point::new(31072, 21262),
            Point::new(31071, 21262),
            Point::new(31071, 21260)]);

        let result = dissolve_into_rings(ls);
        assert_eq!(result.len(), 4);
        assert_eq!(result[0], LineString(vec![Point::new(31071, 21259), Point::new(31071, 21258), Point::new(31072, 21258), Point::new(31072, 21259), Point::new(31071, 21259)]));
        assert_eq!(result[1], LineString(vec![Point::new(31073, 21264), Point::new(31074, 21264), Point::new(31074, 21265), Point::new(31073, 21265), Point::new(31073, 21264)]));
        assert_eq!(result[2], LineString(vec![Point::new(31072, 21262), Point::new(31073, 21262), Point::new(31073, 21264), Point::new(31072, 21264), Point::new(31072, 21262)]));
        assert_eq!(result[3], LineString(vec![Point::new(31071, 21260), Point::new(31072, 21260), Point::new(31072, 21262), Point::new(31071, 21262), Point::new(31071, 21260)]));
    }

    #[test]
    fn dissolve_into_rings6() {
        // b--c
        // | /
        // a
        // | \
        // e--d
        let b = Point::new(0, 0); let c = Point::new(5, 0);
        let a = Point::new(0, 5);
        let e = Point::new(0, 10); let d = Point::new(5, 10);
        let ls: LineString<_> = vec![a, b, c, a, d, e, a].into();
        let result = dissolve_into_rings(ls);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], vec![a, b, c, a].into());
        assert_eq!(result[1], vec![a, d, e, a].into());
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
    fn does_ray_cross4() {
        assert_eq!(does_ray_cross(&(50, 3).into(), &(50, 2).into(), &(49, 3).into()), Crossing::EndPoint);
        assert_eq!(does_ray_cross(&(50, 3).into(), &(49, 3).into(), &(50, 2).into()), Crossing::StartPoint);
    }

    #[test]
    fn does_ray_cross5() {
        assert_eq!(does_ray_cross(&(0, 0).into(), &(1, 0).into(), &(0, 1).into()), Crossing::No);
        assert_eq!(does_ray_cross(&(0, 0).into(), &(0, 1).into(), &(1, 0).into()), Crossing::No);
        assert_eq!(does_ray_cross(&(0, 0).into(), &(-1, 0).into(), &(0, -1).into()), Crossing::StartPoint);
        assert_eq!(does_ray_cross(&(0, 0).into(), &(0, -1).into(), &(-1, 0).into()), Crossing::EndPoint);

        assert_eq!(does_ray_cross(&(0, 0).into(), &(0, -1).into(), &(1, 0).into()), Crossing::No);
        assert_eq!(does_ray_cross(&(0, 0).into(), &(1, 0).into(), &(0, -1).into()), Crossing::No);
    }

    #[test]
    fn does_ray_cross6() {
        assert_eq!(does_ray_cross(&(0, 0).into(), &(-5, 5).into(), &(0, 5).into()), Crossing::No);
        assert_eq!(does_ray_cross(&(0, 0).into(), &(-5, 5).into(), &(3, 1).into()), Crossing::No);
    }

    #[test]
    fn does_ray_cross7() {
        assert_eq!(does_ray_cross(&(0, 0).into(), &(0, 0).into(), &(0, 5).into()), Crossing::Touches);
        assert_eq!(does_ray_cross(&(0, 0).into(), &(0, 5).into(), &(0, 0).into()), Crossing::Touches);

        assert_eq!(does_ray_cross(&(0, 0).into(), &(0, 5).into(), &(0, -5).into()), Crossing::Touches);
        assert_eq!(does_ray_cross(&(0, 0).into(), &(-1, 1).into(), &(1, -1).into()), Crossing::Touches);
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

    #[test]
    fn order_points2() {
        assert_eq!(order_points( ((29147, 10518), (17365, 10520)), (-16552, 10518), (-4238, 10518) ), Ordering::Greater );
    }

}

