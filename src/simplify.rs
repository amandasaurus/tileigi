use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{DivAssign,Rem,Mul,AddAssign};

use fraction::Fraction;

use geo::*;

/// We have a fraction a²/b², but we currently only have a & b². We want to reduce this fraction by
/// removing common multiples so that the fraction is the. It returns the new (a, b²).
/// The results of this will be used later to make the fraction when we calculate a², and we want
/// to reduce the chance of overflow
fn reduce_fraction_sqr(mut a: i64, mut b2: i64) -> (i64, i64) {
    let primes = vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
    'outer: loop {
        for prime in primes.iter() {
            if a % prime == 0 && b2 % prime*prime == 0 {
                a  /= prime;
                b2 /= prime*prime;
                continue 'outer;
            }
        }

        // Got to here => no primes used. so break out.
        break;
    }
    (a, b2)
}

fn distance_sqr(a: &Point<i32>, b: &Point<i32>) -> i64 {
    let delta_x = (a.x() as i64 - b.x() as i64).abs();
    let delta_y = (a.y() as i64 - b.y() as i64).abs();
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
        let denominator = distance_sqr(start, end);
        assert!(denominator != 0);
        if denominator < 0 {
            eprintln!("point {:?} start {:?} end {:?}", point, start, end);
            eprintln!("denominator {:?}", denominator);
            panic!();
        }

        assert!(denominator > 0);

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
        // which we turn into a+c compared with b+d to get the diff
        let ac: i64 = (end_y - start_y)*point_x + end_x*start_y;
        let bd: i64 = (end_x - start_x)*point_y + end_y*start_x;
        let numerator: i64 = if ac > bd { ac - bd } else { bd - ac };

        //println!("{} L {} numerator {} denominator {}", file!(), line!(), numerator, denominator);
        let (numerator, denominator) = reduce_fraction_sqr(numerator, denominator);
        //println!("{} L {} numerator {} denominator {}", file!(), line!(), numerator, denominator);

        assert!(numerator.checked_mul(numerator).is_some(), "{} L {}\npoint {:?} start {:?} end {:?}\nnumerator {:?}\nac {:?} bd {:?}", file!(), line!(), point, start, end, numerator, ac, bd);
        let numerator = numerator*numerator;

        Fraction::new(numerator, denominator)
    }
}

// Ramer-Douglas-Peucker line simplification algorithm
fn rdp(mut points: Vec<Point<i32>>, epsilon: i32) -> Vec<Point<i32>> {
    //println!("{} L {}", file!(), line!());
    let initial_num_points = points.len();
    if initial_num_points <= 2 {
        return points;
    }

    let mut points_to_keep: Vec<bool> = vec![true; points.len()];
    let mut segments_to_look_at: Vec<(usize, usize)> = vec![];

    segments_to_look_at.push((0, initial_num_points-1));

    let e = (epsilon as i64).pow(2);

    let mut index: usize;
    let mut wipe_segment;

    loop {
        //println!("{}:{} segments_to_look_at.len() {}, segments_to_look_at {:?}", file!(), line!(), segments_to_look_at.len(), segments_to_look_at);
        let (start_idx, end_idx) = match segments_to_look_at.pop() {
            None => { break; },
            Some(x) => x,
        };
        // There should not be a case where start/end has been 'removed'
        debug_assert!(points_to_keep[start_idx]);
        debug_assert!(points_to_keep[end_idx]);

        if start_idx + 1 == end_idx || start_idx == end_idx {
            continue;
        }
        debug_assert!(start_idx+1 != end_idx);
        debug_assert!(start_idx < end_idx);

        wipe_segment = false;

        let point1 = points[start_idx];
        let point2 = points[end_idx];

        if point1 == point2 {
            //println!("{}:{} points are the same", file!(), line!());
            // they're the same, so just look at the distance from point to all the other points
            let (max_numerator, this_index) = points[start_idx+1..end_idx].iter().enumerate()
                .zip(points_to_keep[start_idx+1..end_idx].iter())
                .filter_map(|((i, point), &keep)| if keep {
                    Some(((point.x() as i64 - point1.x() as i64).pow(2) + (point.y() as i64 - point2.y() as i64).pow(2), i))
                 } else { None } )
                .max()
                .unwrap();

            index = this_index+start_idx+1;
            debug_assert!(max_numerator > 0);

            // In this case, the numerator is the distance from the furthest point to the original
            // point, squared. i.e. numerator = distance². So we only need to compare it with e².
            wipe_segment = max_numerator < e;
            //println!("{}:{} max_numerator {} this_e {} wipe_segment {}", file!(), line!(), max_numerator, this_e, wipe_segment);
        } else {
            //println!("{}:{} points are different", file!(), line!());
            let delta_x = point2.x() as i64 - point1.x() as i64;
            let delta_y = point2.y() as i64 - point1.y() as i64;
            let end_x_start_y = point2.x() as i64*point1.y() as i64;
            let end_y_start_x = point2.y() as i64*point1.x() as i64;
        
            // the square of the distance between point1 and point2
            let mut point_distance_sqr = delta_x.pow(2) + delta_y.pow(2);
            debug_assert!(point_distance_sqr > 0);
            //println!("{}:{} point_distance_sqr {:?}", file!(), line!(), point_distance_sqr);
            
            let (max_numerator, this_index) = points[start_idx+1..end_idx].iter().enumerate()
                .zip(points_to_keep[start_idx+1..end_idx].iter())
                .filter_map(|((i, point), &keep)| if keep {
                    let ac: i64 = delta_y*(point.x() as i64) + end_x_start_y;
                    let bd: i64 = delta_x*(point.y() as i64) + end_y_start_x;
                    Some( ((ac - bd).abs(), i) )
                 } else { None } )
                .max()
                .unwrap();

            index = this_index+start_idx+1;

            debug_assert!(max_numerator > 0);
            // We want to know if numerator/distance > epsilon
            // distance = sqrt(point_distance_sqr)
            // numerator/sqrt(point_distance_sqr) > epsilon
            // numerator²/point_distance_sqr > epsilon²
            // numerator² > e*point_distance_sqr
            let this_e = e*point_distance_sqr;
            wipe_segment = max_numerator.pow(2) < this_e;
            //println!("{}:{} max_numerator {} this_e {} wipe_segment {}", file!(), line!(), max_numerator, this_e, wipe_segment);
        }


        if wipe_segment {
            if start_idx == 0 && end_idx == initial_num_points - 1 {
                // The entire line should be simplified away, keeping only the start & end points.
                // So short circuit that logic
                return vec![point1, point2];
            }

            for flag in points_to_keep[start_idx+1..end_idx].iter_mut() {
                *flag = false;
            }
        } else {
            segments_to_look_at.push((start_idx, index));
            segments_to_look_at.push((index, end_idx));
        }
    }
    //println!("{} L {}", file!(), line!());

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
    //println!("{} L {}", file!(), line!());
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

pub fn remove_unneeded_points(mut geom: Geometry<i32>) -> Option<Geometry<i32>> {
    remove_points_in_line(&mut geom);
    remove_duplicate_points(&mut geom);
    let geom = remove_spikes(geom);

    geom
}

/// Remove unneeded points from lines.
/// 3 (or more) points in a straight line, can be simplified down to just the end points. You can
/// remove the intermediate points.
/// This is faster than a real simplification, and reduces the number of points, which makes
/// actual simplification faster.
pub fn remove_points_in_line(geom: &mut Geometry<i32>) {
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

fn remove_points_in_line_linestring(ls: &mut LineString<i32>) {
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

    let new_points: Vec<Point<_>> = ls.0.drain(..).zip(keeps.into_iter()).filter_map(|(point, keep)| if keep { Some(point) } else { None }).collect();

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

pub fn remove_spikes_linestring<T: CoordinateType+Debug>(ls: LineString<T>) -> Option<LineString<T>> {
    let LineString(mut points) = ls;
    if points.len() < 2 {
        return Some(LineString(points));
    }
    //println!("{}:{} points {:?}", file!(), line!(), points);
    let mut points_to_keep = vec![true; points.len()];
    let mut have_removed_points = false;

    loop {
        have_removed_points = false;
        // Technically points will get smaller all the time, and this slice is the same size, so it
        // will be too big. But that doesn't matter
        for keep in points_to_keep.iter_mut() {
            *keep = true;
        }

        {
            //println!("{}:{} Start of loop points: {}", file!(), line!(), points.iter().map(|p| format!("({:?},{:?})", p.x(), p.y())).collect::<Vec<_>>().join(", "));
            let mut last_kept = &points[0];

            for (these_points, keep) in points[1..].windows(2).zip(points_to_keep[1..].iter_mut()) {
                let p1 = *last_kept;
                let p2 = these_points[0];
                let p3 = these_points[1];

                let zero_area = (p1.x() - p3.x())*(p2.y() - p1.y()) == (p1.x() - p2.x())*(p3.y() - p1.y());
                //println!("{}:{} these_points {:?} zero_area {} last_kept {:?}", file!(), line!(), these_points, zero_area, last_kept);
                if zero_area {
                    *keep = false;
                    have_removed_points = true;
                } else {
                    last_kept = &these_points[0];
                }
            }
        }

        //println!("{}:{} points:\n{}", file!(), line!(), points.iter().zip(points_to_keep.iter()).map(|(p, k)| format!("{:?} {:?} {}", p.x(), p.y(), k)).collect::<Vec<_>>().join("\n"));
        if have_removed_points {
            points = points.into_iter().zip(points_to_keep.iter()).filter_map(|(p, k)| if *k { Some(p) } else { None }).collect();
            continue;
        } else {
            break;
        }
    }

    Some(LineString(points))

}

pub fn remove_spikes<T: CoordinateType+Debug>(geom: Geometry<T>) -> Option<Geometry<T>> {
    match geom {
        Geometry::LineString(ls) => remove_spikes_linestring(ls).map(Geometry::LineString),
        Geometry::MultiLineString(mls) => {
            let MultiLineString( linestrings ) = mls;
            let mut new_linestrings: Vec<LineString<T>> = linestrings.into_iter().filter_map(|ls| remove_spikes_linestring(ls)).collect();

            match new_linestrings.len() {
                0 => None,
                1 => Some(Geometry::LineString(new_linestrings.remove(0))),
                _ => Some(Geometry::MultiLineString(MultiLineString(new_linestrings))),
            }
        },
        Geometry::Polygon(p) => {
            let Polygon{ exterior, interiors } = p;
            match remove_spikes_linestring(exterior) {
                None => None,
                Some(exterior) => {
                    let these_interiors = interiors.into_iter().filter_map(|ls| remove_spikes_linestring(ls)).collect();
                    Some(Geometry::Polygon(Polygon::new(exterior, these_interiors)))
                }
            }
        },
        Geometry::MultiPolygon(mp) => {
            let MultiPolygon( polygons ) = mp;
            let mut new_polygons: Vec<Polygon<T>> = polygons.into_iter().filter_map(|p| {
                let Polygon{ exterior, interiors } = p;
                match remove_spikes_linestring(exterior) {
                    None => None,
                    Some(exterior) => {
                        let these_interiors = interiors.into_iter().filter_map(|ls| remove_spikes_linestring(ls)).collect();
                        Some(Polygon::new(exterior, these_interiors))
                    }
                }
            }).collect();

            match new_polygons.len() {
                0 => None,
                1 => Some(Geometry::Polygon(new_polygons.remove(0))),
                _ => Some(Geometry::MultiPolygon(MultiPolygon(new_polygons))),
            }
        }
        x => Some(x),
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

    #[test]
    fn reduce_fraction_sqr1() {
        assert_eq!(reduce_fraction_sqr(0, 1), (0, 1));
        assert_eq!(reduce_fraction_sqr(2, 10*10), (1, 5*5));
        assert_eq!(reduce_fraction_sqr(2, 5*5), (2, 5*5));
        assert_eq!(reduce_fraction_sqr(4, 20*20), (1, 5*5));

        assert_eq!(reduce_fraction_sqr(3, 15*15), (1, 5*5));
        assert_eq!(reduce_fraction_sqr(6, 30*30), (1, 5*5));
    }

    #[test]
    fn remove_spikes_linestring1() {
        // known good simple cases
        assert_eq!(remove_spikes_linestring(Vec::<(i32, i32)>::new().into()), Some(Vec::<(i32, i32)>::new().into()));
        assert_eq!(remove_spikes_linestring(vec![(0, 0)].into()), Some(vec![(0, 0)].into()));
    }

    #[test]
    fn remove_spikes_linestring2() {
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (1, 0)].into()), Some(vec![(0, 0), (1, 0)].into()));

        // regular lines/shapes with no spikes
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (1, 0), (1, 1)].into()), Some(vec![(0, 0), (1, 0), (1, 1)].into()));
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)].into()), Some(vec![(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)].into()));
    }

    #[test]
    fn remove_spikes_linestring3() {
        // has a point that is collinear with another
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (1, 0), (2, 0)].into()), Some(vec![(0, 0), (2, 0)].into()));
        // Several points collinear
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (1, 0), (2, 0), (5, 0)].into()), Some(vec![(0, 0), (5, 0)].into()));
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (1, 0), (2, 0), (3, 0), (5, 0)].into()), Some(vec![(0, 0), (5, 0)].into()));
    }

    #[test]
    fn remove_spikes_linestring4() {
        // Has a spike
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (1, 0), (1, 1), (10, 10), (1, 1), (0, 1), (0, 0)].into()), Some(vec![(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)].into()));

        // Has a longer spike
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (1, 0), (1, 1), (10, 10), (20, 10), (10, 10), (1, 1), (0, 1), (0, 0)].into()), Some(vec![(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)].into()));

    }

    #[test]
    fn remove_spikes_linestring5() {
        // Has a 'turn'
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (0, 100), (0, 50)].into()), Some(vec![(0, 0), (0, 50)].into()));
        assert_eq!(remove_spikes_linestring(vec![(0, 0), (0, 100), (0, 50), (50, 50)].into()), Some(vec![(0, 0), (0, 50), (50, 50)].into()));
    }

}
