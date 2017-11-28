use geo::*;
use geo::map_coords::MapCoords;
use geo::intersects::Intersects;
use std::cmp::{min, max};

pub fn is_valid<T: CoordinateType>(geom: &Geometry<T>) -> bool {
    match *geom {
        Geometry::LineString(ref ls) => is_linestring_valid(ls),
        Geometry::Polygon(ref p) => is_polygon_valid(p),
        Geometry::MultiPolygon(ref mp) => mp.0.iter().all(|p| is_polygon_valid(p)),
        _ => true,
    }
}

pub fn is_linestring_valid<T: CoordinateType>(ls: &LineString<T>) -> bool {
    if ls.0.len() < 2 {
        return false;
    }

    if ls.0.len() == 2 && ls.0[0] == ls.0[1] {
        return false;
    }

    true
}

pub fn is_polygon_valid<T: CoordinateType>(p: &Polygon<T>) -> bool {
    // Sometimes there are duplicate points, e.g. A-A-B-A. If we remove all dupes, we can see if
    // there are <4 points
    // Gah clones!
    let mut ext = p.exterior.clone();
    remove_duplicate_points_linestring(&mut ext);
    // TODO fix clipping code etc to not make linestrings with duplicated points

    if ext.0.len() < 4 {
        return false;
    }

    if ext.0.iter().skip(1).all(|&pt| pt == ext.0[0]) {
        // All points the same
        return false;
    }

    for i in p.interiors.iter() {
        let mut i = i.clone();
        remove_duplicate_points_linestring(&mut i);

        if i.0.len() < 4 {
            return false;
        }

        if i.0.iter().skip(1).all(|&pt| pt == i.0[0]) {
            // All points the same
            return false;
        }
    }

    // In theory this is backwards. Ext rings should be CCW, and int rings CW. But in vector tiles
    // the y goes down, so it's flipped.
    if ext.is_ccw() || p.interiors.iter().any(|i| i.is_cw()) {
        return false;
    }

    true
}


fn remove_duplicate_points_linestring<T: CoordinateType>(ls: &mut LineString<T>) {
    let mut i = 0;

    // This could be more effecient for cases of many duplicate points in a row.
    loop {
        if i >= ls.0.len()-1 {
            break;
        }

        if ls.0[i] == ls.0[i+1] {
            ls.0.remove(i+1);
        } else {
            i += 1;
        }

    }
}


pub fn remove_duplicate_points<T: CoordinateType>(geom: &mut Geometry<T>) {
    match *geom {
        Geometry::LineString(ref mut ls) => remove_duplicate_points_linestring(ls),
        _ => {},
    }
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
        _ => {},
    }
}

fn has_self_intersections<T: CoordinateType>(ls: &LineString<T>) -> bool {
    false
}

/// True iff the segments ab intersect at any point except their endpoints
fn intersect<P: Into<Point<T>>, T: CoordinateType>(a: P, b: P, c: P, d: P) -> bool {
    // This really need improving
    // FIXME add initiall bbox check which should speed it up
    let a: Point<T> = a.into();
    let b: Point<T> = b.into();
    let c: Point<T> = c.into();
    let d: Point<T> = d.into();

    assert!(a != b);
    assert!(c != d);

    let (x1, y1) = (a.x(), a.y());
    let (x2, y2) = (b.x(), b.y());
    let (x3, y3) = (c.x(), c.y());
    let (x4, y4) = (d.x(), d.y());

    let a = x2 - x1;
    let b = x3 - x4;
    let c = y2 - y1;
    let d = y3 - y4;

    let e = x3 - x1;
    let f = y3 - y1;

    let determinate = a*d - b*c;
    if determinate == T::zero() {
        return false;
    }
    // FIXME what if determinate < 0 ?
    assert!(determinate >= T::zero());

    let sd = d*e - b*e;
    if sd > determinate || sd < T::zero() {
        return false
    }

    let td = e*f - a*f;
    if td > determinate || td < T::zero() {
        return false
    }

    if (td == determinate || td == T::zero()) && (sd == T::zero() || sd == determinate) {
        // endpoints overlap
        return false;
    }

    true
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn intersect1() {
        assert!(!intersect((0, 0), (0, 10), (5, 1), (5, 2)));
        assert!(intersect((0, 0), (0, 10), (0, 5), (5, 5)));

        assert!(!intersect((0, 0), (0, 10), (0, 10), (0, 20)));
        assert!(intersect((0, 0), (0, 10), (-5, 5), (5, 5)));
    }


}

