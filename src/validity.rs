use geo::*;

pub fn is_valid<T: CoordinateType>(geom: &Geometry<T>) -> bool {
    match *geom {
        Geometry::LineString(ref ls) => is_linestring_valid(ls),
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

