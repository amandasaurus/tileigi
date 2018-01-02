use geo::*;
use geo::map_coords::MapCoords;
use geo::intersects::Intersects;
use std::cmp::{min, max, Ord};
use num_traits::Signed;
use std::fmt::Debug;

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

        if has_self_intersections(&i) {
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

    for (i, points12) in ls.0.windows(2).enumerate() {
        let (p1, p2) = (points12[0], points12[1]);
        for points34 in ls.0.windows(2).skip(i+1) {
            let (p3, p4) = (points34[0], points34[1]);

            if intersect_excl_end(p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y(), p4.x(), p4.y()) {
                //println!("\n\nintersection i {} p1 {:?} p2 {:?} p3 {:?} p4 {:?}", i, p1, p2, p3, p4);
                return true;
            }
        }
    }


    false
}

/// True iff the segments |p1p2| and |p3p4| intersect at any point except their endpoints
fn intersect_excl_end<T: CoordinateType+Signed+Debug+Ord>(x1: T, y1: T, x2: T, y2: T, x3: T, y3: T, x4: T, y4: T) -> bool {
    // FIXME add initiall bbox check which should speed it up
    
    if max(x1, x2) < min(x3, x4) || min(x1, x2) > max(x3, x4)
        || max(y1, y2) < min(y3, y4) || min(y1, y2) > max(y3, y4)
    {
        return false;
    }
    
    debug_assert!((x1, y1) != (x2, y2));
    debug_assert!((x3, y3) != (x4, y4));

    let a = x2 - x1;
    let b = x3 - x4;
    let c = y2 - y1;
    let d = y3 - y4;

    let determinate = a*d - b*c;
    if determinate == T::zero() {
        return false;
    }

    let e = x3 - x1;
    let f = y3 - y1;

    // we know it's not zero
    let signum = determinate.signum();
    let determinate = determinate.abs();

    let sd = signum * (a*f - c*e);
    if sd > determinate || sd < T::zero() {
        return false
    }

    let td = signum*(d*e - b*f);
    if td > determinate || td < T::zero() {
        return false
    }

    if (td == determinate || td == T::zero()) && (sd == T::zero() || sd == determinate) {
        // endpoints overlap
        return false;
    }

    if td >= T::zero() && td <= determinate && sd >= T::zero() && sd <= determinate {
        return true;
    }

    // Should have been caught above.
    println!("points {:?} {:?} {:?} {:?}", (x1, y1), (x2, y2), (x3, y3), (x4, y4));
    println!("det {:?} ds {:?} dt {:?}", determinate, sd, td);
    unreachable!();
}

#[cfg(test)]
mod test {
    use super::*;
    use super::super::validity::is_valid;

    //#[test]
    //fn intersect1() {
    //    assert!(!intersect_excl_end(0, 0,  0, 10,  5, 1,  5, 2));
    //    assert!(intersect_excl_end((0, 0), (0, 10), (0, 5), (5, 5)));

    //    assert!(!intersect_excl_end((0, 0), (0, 10), (0, 10), (0, 20)));
    //    assert!(intersect_excl_end((0, 0), (0, 10), (-5, 5), (5, 5)));
    //    assert!(!intersect_excl_end((0, 0), (10, 0), (10, 0), (10, 10)));
    //    assert!(intersect_excl_end((-5, 5), (5, 5), (0, 0), (0, 10)));

    //}
    
}

