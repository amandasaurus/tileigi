use std::cmp::Ordering;
use std::fmt::Debug;

use geo::*;

#[derive(Debug)]
struct Fraction<T: CoordinateType> {
    numerator: T,
    denominator: T,
}

impl<T: CoordinateType> Fraction<T> {
    fn new(numerator: T, denominator: T) -> Self {
        Fraction{ numerator, denominator }
    }
}

impl<T: CoordinateType> PartialEq for Fraction<T> {
    fn eq(&self, other: &Self) -> bool {
        let a = self.numerator;
        let b = self.denominator;
        let c = other.numerator;
        let d = other.denominator;

        if b == d {
            a == c
        } else {
            a*d == b*c
        }
    }
}

impl<T: CoordinateType> PartialOrd for Fraction<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let a = self.numerator;
        let b = self.denominator;
        let c = other.numerator;
        let d = other.denominator;
        if b == d {
            a.partial_cmp(&c)
        } else {
            // In our use case, we can assume this is true
            assert!(b != T::zero());
            assert!(b > T::zero());
            assert!(d != T::zero());
            assert!(d > T::zero());

            (a*d).partial_cmp(&(b*c))
        }
    }
}

fn distance_sqr(a: &Point<i32>, b: &Point<i32>) -> i64 {
    let delta_x = if a.x() > b.x() { a.x() - b.x() } else { b.x() - a.x() };
    let delta_y = if a.y() > b.y() { a.y() - b.y() } else { b.y() - a.y() };
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
        // point = (x0, y0)
        // |(y2 - y1)x0 - (x2 - x1)y0 + x2y1 - y2x1|
        // |a - b + c - d|
        // |a+c - (b+d)|
        // which we turn into a + c compared with b+d to get the diff
        let ac = (end.y() - start.y())*point.x() + end.x()*start.y();
        let bd = (end.x() - start.x())*point.y() + end.y()*start.x();
        let ac = ac as i64;
        let bd = bd as i64;
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

// Ramer–Douglas-Peucker line simplification algorithm
fn rdp(points: &[Point<i32>], epsilon: i32) -> Vec<Point<i32>>
{
    if points.is_empty() {
        return points.to_vec();
    }
    let mut dmax_sqr: Fraction<i64> = Fraction::new(0, 1);
    let mut index: usize = 0;
    let mut distance_sqr: Fraction<i64>;

    for (i, _) in points.iter().enumerate().take(points.len() - 1).skip(1) {
        distance_sqr = point_line_distance_sqr(&points[i],
                                       &points[0],
                                       &*points.last().unwrap());
        if distance_sqr > dmax_sqr {
            index = i;
            dmax_sqr = distance_sqr;
        }
    }
    let e = epsilon as i64;
    let e = Fraction::new(e*e, 1);
    if dmax_sqr > e {
        let mut intermediate = rdp(&points[..index + 1], epsilon);
        intermediate.pop();
        intermediate.extend_from_slice(&rdp(&points[index..], epsilon));
        intermediate
    } else {
        vec![*points.first().unwrap(), *points.last().unwrap()]
    }
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
    let new_points = rdp(&geom.0, epsilon);

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

#[cfg(test)]
mod test {
}
