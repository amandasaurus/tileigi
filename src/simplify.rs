use std::cmp::Ordering;
use std::fmt::Debug;

use fraction::Fraction;

use geo::*;

fn distance_sqr(a: &Point<i32>, b: &Point<i32>) -> i64 {
    let delta_x = (a.x() - b.x()).abs();
    let delta_y = (a.y() - b.y()).abs();
    let delta_x = delta_x as i64;
    let delta_y = delta_y as i64;
    assert!(delta_x >= 0);
    assert!(delta_y >= 0);

    (delta_x*delta_x + delta_y*delta_y)
}

// perpendicular distance from a point to a line
fn point_line_distance_sqr(point: &Point<i32>, start: &Point<i32>, end: &Point<i32>) -> Fraction<i64>
{
    if start == end {
        Fraction::new(distance_sqr(point, start), 1)
    } else {
        let start_x = start.y() as i64;
        let start_y = start.y() as i64;
        let end_x = end.y() as i64;
        let end_y = end.y() as i64;
        let point_x = point.x() as i64;
        let point_y = point.y() as i64;
        // point = (x0, y0)
        // |(y2 - y1)x0 - (x2 - x1)y0 + x2y1 - y2x1|
        // |a - b + c - d|
        // |a+c - (b+d)|
        // which we turn into a + c compared with b+d to get the diff
        let ac = (end_y - start_y)*point_x + end_x*start_y;
        let bd = (end_x - start_x)*point_y + end_y*start_x;
        let numerator = if ac > bd { ac - bd } else { bd - ac };
        let numerator = numerator*numerator;

        let denominator = distance_sqr(start, end);
        assert!(denominator != 0);
        if denominator < 0 {
            eprintln!("point {:?} start {:?} end {:?}", point, start, end);
            eprintln!("denominator {:?}", denominator);
            panic!();
        }

        assert!(denominator > 0);

        Fraction::new(numerator, denominator)
    }
}

// Ramer-Douglas-Peucker line simplification algorithm
fn rdp(mut points: Vec<Point<i32>>, epsilon: i32) -> Vec<Point<i32>>
{
    let initial_num_points = points.len();
    if points.is_empty() || initial_num_points <= 2 {
        return points;
    }

    let mut points_to_keep: Vec<bool> = vec![true; points.len()];
    let mut num_points_to_keep = points.len();
    let mut segments_to_look_at: Vec<(usize, usize)> = vec![];

    segments_to_look_at.push((0, initial_num_points-1));

    let e = epsilon as i64;
    let e = Fraction::new(e*e, 1);

    let mut dmax_sqr: Fraction<i64> = Fraction::new(0, 1);
    let mut distance_sqr: Fraction<i64>;
    let mut index: usize;

    loop {
        let (start_idx, end_idx) = match segments_to_look_at.pop() {
            None => { break; },
            Some(x) => x,
        };

        if start_idx + 1 == end_idx {
            continue;
        }
        assert!(start_idx+1 != end_idx);
        
        index = start_idx;
        dmax_sqr = Fraction::new(0, 1);

        for (i, point) in points[start_idx+1..end_idx].iter().enumerate() {
            if points_to_keep[i+start_idx] {
                distance_sqr = point_line_distance_sqr(&point,
                                               &points[start_idx],
                                               &points[end_idx]);
                if distance_sqr > dmax_sqr {
                    index = i+start_idx+1;
                    dmax_sqr = distance_sqr;
                }
            }
        }
        if dmax_sqr > e {
            segments_to_look_at.push((start_idx, index));
            segments_to_look_at.push((index, end_idx));
        } else {
            for flag in points_to_keep[start_idx+1..end_idx].iter_mut() {
                *flag = false;
            }
        }
    }

    let new_points: Vec<_> = points.drain(..).zip(points_to_keep.into_iter()).filter_map(|(point, keep)| if keep { Some(point) } else { None }).collect();

    return new_points;
}

pub fn simplify(geom: Geometry<i32>, epsilon: i32) -> Option<Geometry<i32>> {
    match geom {
        // Can't simplify a Point. let's hope this doesn't do a copy or memory move or something
        Geometry::Point(p) => Some(Geometry::Point(p)),
        Geometry::MultiPoint(p) => Some(Geometry::MultiPoint(p)),

        Geometry::LineString(ls) => simplify_linestring(ls, epsilon, false).map(|g| g.into()),
        Geometry::MultiLineString(mls) => simplify_multilinestring(mls, epsilon).map(|g| g.into()),
        Geometry::Polygon(p) => simplify_polygon(p, epsilon).map(|g| g.into()),
        Geometry::MultiPolygon(mp) => simplify_multipolygon(mp, epsilon).map(|g| g.into()),

        Geometry::GeometryCollection(_) => unimplemented!(),
    }
}

fn simplify_linestring(geom: LineString<i32>, epsilon: i32, should_be_ring: bool) -> Option<LineString<i32>> {
    let LineString(points) = geom;
    let new_points = rdp(points, epsilon);

    if should_be_ring {
        if new_points.len() >= 4 && new_points[0] == new_points[new_points.len()-1] { 
            Some(LineString(new_points))
        } else {
            None
        }
    } else {
        if new_points.len() >= 2 {
            Some(LineString(new_points))
        } else {
            None
        }
    }
}

fn simplify_multilinestring(geom: MultiLineString<i32>, epsilon: i32) -> Option<MultiLineString<i32>> {
    Some(MultiLineString(geom.0.into_iter().filter_map(|l| simplify_linestring(l, epsilon, false)).collect()))
}

fn simplify_polygon(geom: Polygon<i32>, epsilon: i32) -> Option<Polygon<i32>> {
    let Polygon{ exterior, interiors } = geom;
    match simplify_linestring(exterior, epsilon, true) {
        None => None,
        Some(new_exterior) => {
            Some(Polygon::new(new_exterior, interiors.into_iter().filter_map(|l| simplify_linestring(l, epsilon, true)).collect()))
        }
    }

}

fn simplify_multipolygon(geom: MultiPolygon<i32>, epsilon: i32) -> Option<MultiPolygon<i32>> {
    let new_polygons: Vec<_> = geom.0.into_iter().filter_map(|p| simplify_polygon(p, epsilon)).collect();
    if new_polygons.is_empty() {
        None
    } else {
        Some(MultiPolygon(new_polygons))
    }
}

pub fn remove_unneeded_points<T: CoordinateType+Debug>(mut geom: &mut Geometry<T>) {
    remove_points_in_line(&mut geom);
    remove_duplicate_points(&mut geom);
    remove_spikes(&mut geom);
}

/// Remove unneeded points from lines.
/// 3 (or more) points in a straight line, can be simplified down to just the end points. You can
/// remove the intermediate points.
/// This is faster than a real simplification, and reduces the number of points, which makes
/// actual simplification faster.
pub fn remove_points_in_line<T: CoordinateType+Debug>(geom: &mut Geometry<T>) {
    match *geom {
        Geometry::LineString(ref mut ls) => remove_points_in_line_linestring(ls),
        Geometry::MultiLineString(ref mut mls) => {
            for mut ls in mls.0.iter_mut() {
                remove_points_in_line_linestring(&mut ls);
            }
        },
        Geometry::Polygon(ref mut p) => {
            remove_points_in_line_linestring(&mut p.exterior);
            for mut i in p.interiors.iter_mut() {
                remove_points_in_line_linestring(&mut i);
            }
            },
        Geometry::MultiPolygon(ref mut mp) => {
            for mut p in mp.0.iter_mut() {
                remove_points_in_line_linestring(&mut p.exterior);
                for mut i in p.interiors.iter_mut() {
                    remove_points_in_line_linestring(&mut i);
                }
            }
        }
        _ => {},
    }
}

fn remove_points_in_line_linestring<T: CoordinateType+Debug>(ls: &mut LineString<T>) {
    if ls.0.len() < 2 {
        return;
    }

    let slopes: Vec<_> = ls.0.windows(2).map(|points| {
        Fraction::new( points[1].x() - points[0].x(), points[1].y() - points[0].y() )
    }).collect();

    let mut keeps: Vec<bool> = vec![false; ls.0.len()];
    keeps[0] = true;
    let mut num_to_keep = 1;

    {
        // inner scope because we cause an immutable borrow with last_keep.
        let mut last_keep = &slopes[0];

        for (idx, slope) in slopes.iter().enumerate() {
            if slope != last_keep {
                keeps[idx] = true;
                last_keep = slope;
                num_to_keep += 1;
            }
        }
    }

    // we always keep the last one
    keeps[ls.0.len()-1] = true;

    if num_to_keep == ls.0.len() {
        // nothing to do, so early return
        return;
    }

    let new_points: Vec<Point<T>> = ls.0.drain(..).zip(keeps.into_iter()).filter_map(|(point, keep)| if keep { Some(point) } else { None }).collect();

    ::std::mem::replace(&mut ls.0, new_points);
}

pub fn remove_duplicate_points<T: CoordinateType+Debug>(geom: &mut Geometry<T>) {
    match *geom {
        Geometry::LineString(ref mut ls) => remove_duplicate_points_linestring(ls),
        Geometry::MultiLineString(ref mut mls) => {
            for mut ls in mls.0.iter_mut() {
                remove_duplicate_points_linestring(&mut ls);
            }
        },
        Geometry::Polygon(ref mut p) => {
            remove_duplicate_points_linestring(&mut p.exterior);
            for mut i in p.interiors.iter_mut() {
                remove_duplicate_points_linestring(&mut i);
            }
            },
        Geometry::MultiPolygon(ref mut mp) => {
            for mut p in mp.0.iter_mut() {
                remove_duplicate_points_linestring(&mut p.exterior);
                for mut i in p.interiors.iter_mut() {
                    remove_duplicate_points_linestring(&mut i);
                }
            }
        }
        _ => {},
    }
}

fn remove_duplicate_points_linestring<T: CoordinateType+Debug>(ls: &mut LineString<T>) {
    if ls.0.len() < 2 {
        return;
    }

    let mut keeps: Vec<bool> = vec![false; ls.0.len()];
    keeps[0] = true;
    let mut num_to_keep = 1;

    {
        // inner scope because we cause an immutable borrow with last_keep.
        let mut last_keep = &ls.0[0];

        for (idx, point) in ls.0.iter().enumerate().skip(1) {
            if point != last_keep {
                keeps[idx] = true;
                last_keep = point;
                num_to_keep += 1;
            }
        }
    }

    if num_to_keep == ls.0.len() {
        // nothing to do, so early return
        return;
    }

    let new_points: Vec<Point<T>> = ls.0.drain(..).zip(keeps.into_iter()).filter_map(|(point, keep)| if keep { Some(point) } else { None }).collect();

    ::std::mem::replace(&mut ls.0, new_points);


    loop {
        let len = ls.0.len();
        if len <= 2 { break; }
        if ls.0[len-1] != ls.0[len-2] { break; }
        // the last point is duplicated from the 2nd last
        ls.0.remove(len-1);
    }
}

fn remove_spikes_linestring<T: CoordinateType+Debug>(ls: &mut LineString<T>) {
    // TODO There is definitely a more effecient way to do this. There is certainly a lot of
    // memory moving.

    // we could have a few points in a line in a spike, but at each run of this loop, it'll only
    // remove the end point of the spike. So we need to run it again and again to remove all the
    // points, to entirely remove the spike.
    loop {
        let mut have_removed_points = false;

        // keep_point[i] = true means we should keep this point. false means remove it.
        let mut keep_point = vec![true; ls.0.len()];
        let mut i = 1;
        let mut points_removed = 0;

        loop {
            if i >= (ls.0.len() - 1) {
                break;
            }

            let p1 = ls.0[i-1];
            let p2 = ls.0[i];
            let p3 = ls.0[i+1];
            // This algorithm is 'twice the triangle area'. If it's 0 (i.e. both sides are equal),
            // then it's a zero area, i.e. spike.
            if (p1.x() - p3.x())*(p2.y() - p1.y()) == (p1.x() - p2.x())*(p3.y() - p1.y()) {
                keep_point[i] = false;
                points_removed += 1;
            }
            i += 1;
        }

        if points_removed > 0 {
            // Create a new vec of points but only the ones we want to keep
            let new_points: Vec<Point<T>> = ls.0.drain(..).zip(keep_point.into_iter()).filter_map(|(p, keep)| if keep { Some(p) } else { None }).collect();
            
            // and move that in for the points
            ::std::mem::replace(&mut ls.0, new_points);
            continue;
        } else {
            break;
        }
    }

}

pub fn remove_spikes<T: CoordinateType+Debug>(geom: &mut Geometry<T>) {
    match *geom {
        Geometry::LineString(ref mut ls) => remove_spikes_linestring(ls),
        Geometry::MultiLineString(ref mut mls) => {
            for mut ls in mls.0.iter_mut() {
                remove_spikes_linestring(&mut ls);
            }
        },
        Geometry::Polygon(ref mut p) => {
            remove_spikes_linestring(&mut p.exterior);
            for mut i in p.interiors.iter_mut() {
                remove_spikes_linestring(&mut i);
            }
            },
        Geometry::MultiPolygon(ref mut mp) => {
            for mut p in mp.0.iter_mut() {
                remove_spikes_linestring(&mut p.exterior);
                for mut i in p.interiors.iter_mut() {
                    remove_spikes_linestring(&mut i);
                }
            }
        }
        _ => {},
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ::validity::is_valid;

    #[test]
    fn test_remove_points_in_line1() {
        let mut geom = Geometry::LineString(LineString(vec![Point(Coordinate { x: 0, y: 0 }), Point(Coordinate { x: 10, y: 0 }), Point(Coordinate { x: 20, y: 0 })]));
        remove_points_in_line(&mut geom);
        assert_eq!(geom, Geometry::LineString(LineString(vec![Point(Coordinate { x: 0, y: 0 }), Point(Coordinate { x: 20, y: 0 })])));
    }

    #[test]
    fn remove_duplicate_points_valid() {
        let mut geom = Geometry::LineString(LineString(vec![Point(Coordinate { x: 31565, y: 20875 }), Point(Coordinate { x: 31615, y: 20887 }), Point(Coordinate { x: 31633, y: 20819 }), Point(Coordinate { x: 31593, y: 20822 }), Point(Coordinate { x: 31585, y: 20808 }), Point(Coordinate { x: 31584, y: 20850 }), Point(Coordinate { x: 31565, y: 20875 })]));
        remove_duplicate_points(&mut geom);
        assert!(is_valid(&geom), "Invalid geometry {:?}", geom);
    }


}
