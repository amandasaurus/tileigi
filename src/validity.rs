use geo::*;

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
    if p.exterior.0.len() < 4 {
        return false;
    }

    if p.exterior.0.iter().skip(1).all(|&pt| pt == p.exterior.0[0]) {
        // All points the same
        return false;
    }

    for i in p.interiors.iter() {
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
    if p.exterior.is_ccw() || p.interiors.iter().any(|i| i.is_cw()) {
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


