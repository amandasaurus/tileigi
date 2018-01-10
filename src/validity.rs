use geo::*;
use geo::map_coords::MapCoords;
use geo::intersects::Intersects;
use std::cmp::{min, max, Ord};
use num_traits::Signed;
use std::fmt::Debug;

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

fn is_polygon_valid<T: CoordinateType+Signed+Debug+Ord>(p: &Polygon<T>) -> bool {
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
            //println!("looking at i {} j {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", i, j, p1, p2, p3, p4);

            match intersection(p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y(), p4.x(), p4.y()) {
                Intersection::Crossing | Intersection::Overlapping  => { return true; },
                Intersection::Touching => { return true; },
                Intersection::None | Intersection::EndToEnd => {},
            }
            //println!("no intersection");
        }
    }


    false
}

#[derive(PartialEq,Eq,Clone,Copy,Debug)]
enum Intersection {
    // They don't intersect/touch at all
    None,

    // One is wholly, or partially, on top of another, ie infinite number of intersecting points,
    // the intersection is a line
    Overlapping,

    // The end point of one is the same as the end point of another,
    EndToEnd,

    // The end of one touches the other, but not at it's end
    Touching,

    // real crossing
    Crossing,
}

fn intersect_incl_end<T: CoordinateType+Signed+Debug+Ord>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T, x4: T, y4: T) -> bool {
    intersection(x1, y1, x2, y2, x3, y3, x4, y4) == Intersection::None
}


/// True iff the segments |p1p2| and |p3p4| intersect at any point, and the intersection point is
/// not on both end points. i.e. 2 lines can join end-to-end in this, but not touch anywhere else.
fn intersection<T: CoordinateType+Signed+Debug+Ord>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T, x4: T, y4: T) -> Intersection {
    if max(x1, x2) < min(x3, x4) || min(x1, x2) > max(x3, x4)
        || max(y1, y2) < min(y3, y4) || min(y1, y2) > max(y3, y4)
    {
        return Intersection::None;
    }
    
    debug_assert!((x1, y1) != (x2, y2));
    debug_assert!((x3, y3) != (x4, y4));

    let a = x2 - x1;
    let b = x3 - x4;
    let c = y2 - y1;
    let d = y3 - y4;

    let determinate = a*d - b*c;
    if determinate == T::zero() {
        let slope_12 = Fraction::new(a, c);
        let slope_34 = Fraction::new(b, d);
        let neg_slope_34 = Fraction::new(-b, d);
        assert!((slope_12 == slope_34) || (slope_12 == neg_slope_34)); // if this is false, then I don't know what's going on.

        let delta_12 = (a, c);
        let delta_13 = (x3 - x1, y3 - y1);
        let delta_14 = (x4 - x1, y4 - y1);
        let delta_23 = (x3 - x2, y3 - y2);
        let delta_24 = (x4 - x2, y4 - y2);

        let zero = (T::zero(), T::zero());
        let p3_on_end = (x1, y1) == (x3, y3) || (x2, y2) == (x3, y3);
        let p4_on_end = (x1, y1) == (x4, y4) || (x2, y2) == (x4, y4);
        let p3_on_12_excl_end = delta_13 > zero && delta_13 < delta_12;
        let p4_on_12_excl_end = delta_14 > zero && delta_14 < delta_12;

        if ((x1, y1) == (x3, y3) && (x2, y2) == (x4, y4)) || ((x1, y1) == (x4, y4) && (x2, y2) == (x3, y3)) {
            return Intersection::Overlapping;
        } else if (p3_on_end && !p4_on_12_excl_end) || (p4_on_end && !p3_on_12_excl_end) {
            return Intersection::EndToEnd;
        } else if (p3_on_end && p4_on_12_excl_end) || (p4_on_end && p3_on_12_excl_end) {
            return Intersection::Overlapping;
        } else if (!p3_on_end && p4_on_12_excl_end) || (!p4_on_end && p3_on_12_excl_end) {
            return Intersection::Overlapping;
        }

        println!("({:?}, {:?}) - ({:?}, {:?})", x1, y1, x2, y2);
        println!("({:?}, {:?}) - ({:?}, {:?})", x3, y3, x4, y4);
        println!("delta_12 {:?}\ndelta_13 {:?} delta_14 {:?}\ndelta_23 {:?} delta_24 {:?}", delta_12, delta_13, delta_14, delta_23, delta_24);
        println!("p3_on_end {} p3_on_12_excl_end {}\np4_on_end {} p4_on_12_excl_end {}", p3_on_end, p3_on_12_excl_end, p4_on_end, p4_on_12_excl_end);
        unreachable!();
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
        return Intersection::Touching
    } else if (td < determinate || td > T::zero()) && (sd == T::zero() || sd == determinate) {
        return Intersection::Touching
    } else if td > T::zero() && td < determinate && sd > T::zero() && sd < determinate {
        return Intersection::Crossing;
    }

    // Should have been caught above.
    println!("points {:?} {:?} {:?} {:?}", (x1, y1), (x2, y2), (x3, y3), (x4, y4));
    println!("det {:?} sd {:?} td {:?}", determinate, sd, td);
    unreachable!();
}

pub fn make_valid<T: CoordinateType>(geom: Geometry<T>) -> Geometry<T> {
    match geom {
        Geometry::Polygon(p) => Geometry::MultiPolygon(make_polygon_valid(p)),
        Geometry::MultiPolygon(mp) => {
            unimplemented!();
        },
        x => x,
    }
}

fn make_polygon_valid<T: CoordinateType>(p: Polygon<T>) -> MultiPolygon<T> {
    MultiPolygon(vec![])

}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn intersect1() {

        assert_eq!(intersection(0, 0,  0, 10,  5, 1,  5, 2), Intersection::None);
        assert_eq!(intersection(0, 0,  0, 10,  0, 5,  5, 5), Intersection::Touching);

        assert_eq!(intersection(0, 0,  0, 10,  0, 0,  0, 10), Intersection::Overlapping);
        assert_eq!(intersection(0, 0,  0, 10,  0, 5,  0, 10), Intersection::Overlapping);
        assert_eq!(intersection(0, 0,  0, 10,  0, 5,  0, 15), Intersection::Overlapping);
        assert_eq!(intersection(0, 0,  0, 10,  0, 0,  0, 5), Intersection::Overlapping);
        assert_eq!(intersection(0, 0,  0, 10,  0, 2,  0, 8), Intersection::Overlapping);

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

        assert_eq!(intersection(0, 0,  0, 10,  -5, 5,  5, 5), Intersection::Crossing);
        assert_eq!(intersection(0, 0,  10, 0,  10, 0,  10, 10), Intersection::EndToEnd);
        assert_eq!(intersection(-5, 5,  5, 5,  0, 0,  0, 10), Intersection::Crossing);
        assert_eq!(intersection(0, 0,  10, 0,  5, 10,  5, -10), Intersection::Crossing);

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
    fn validity_checks() {
        let geom: LineString<i32> = LineString(vec![]);
        assert!(!is_linestring_valid(&geom));

        assert!(!is_linestring_valid(&LineString(vec![(0i32, 0i32).into()])));

        assert!(!is_linestring_valid(&vec![(0, 0), (4, 0), (2, -1), (2, 1)].into()));
        // TODO fix
        //assert!(!is_linestring_valid(&vec![(0, 0), (4, 0), (2, -1), (2, 0), (2, 1)].into()));

        // Simple square - valid
        assert!(is_polygon_valid(&Polygon::new(vec![(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)].into(), vec![])));
        // Unclosed - invalid
        assert!(!is_polygon_valid(&Polygon::new(vec![(0, 0), (0, 1), (1, 1), (1, 0)].into(), vec![])));

        // Has a touching inner
        let geom: Polygon<i32> = Polygon::new(vec![(0, 0), (0, 2), (1, 2), (1, 1), (2, 1), (2, 3), (1, 3), (1, 2), (0, 2), (0, 4), (3, 4), (3, 0), (0, 0)].into(), vec![]);
        assert!(!is_polygon_valid(&geom));
    }
    
    #[test]
    fn test_make_valid() {
        let geom: Polygon<i32> = Polygon::new(vec![(0, 0), (0, 2), (1, 2), (1, 1), (2, 1), (2, 3), (1, 3), (1, 2), (0, 2), (0, 4), (3, 4), (3, 0), (0, 0)].into(), vec![]);
        assert!(!is_polygon_valid(&geom));
        
        let mut new_geom = make_polygon_valid(geom);
        assert_eq!(new_geom.0.len(), 1);
        let new_geom: Polygon<_> = new_geom.0.remove(0);
        assert!(is_polygon_valid(&new_geom));
        assert_eq!(new_geom, Polygon::new(vec![(0, 0), (0, 4), (3, 4), (3, 0), (0, 0)].into(), vec![vec![(1, 1), (2, 1), (2, 3), (1, 3), (1, 1)].into()]));
    }
}

