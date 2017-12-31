use std::cmp::Ordering;
use std::fmt::Debug;

use geo::*;

fn distance_sqr(a: &Point<i64>, b: &Point<i64>) -> i64 {
    let delta_x = if a.x() > b.x() { a.x() - b.x() } else { b.x() - a.x() };
    let delta_y = if a.y() > b.y() { a.y() - b.y() } else { b.y() - a.y() };
    assert!(delta_x >= 0);
    assert!(delta_y >= 0);

    (delta_x*delta_x + delta_y*delta_y)
}

// perpendicular distance from a point to a line
fn point_line_distance(point: &Point<i64>, start: &Point<i64>, end: &Point<i64>) -> f64
{
    if start == end {
        (distance_sqr(point, start) as f64).sqrt()
    } else {
        // point = (x0, y0)
        // |(y2 - y1)x0 - (x2 - x1)y0 + x2y1 - y2x1|
        // |a - b + c - d|
        // |a+c - (b+d)|
        // which we turn into a + c compared with b+d to get the diff
        let ac = (end.y() - start.y())*point.x() + end.x()*start.y();
        let bd = (end.x() - start.x())*point.y() + end.y()*start.x();
        let numerator = if ac > bd { ac - bd } else { bd - ac };

        let denominator: f64 = (distance_sqr(start, end) as f64).sqrt();
        assert!(denominator != 0.);
        if denominator < 0. {
            eprintln!("point {:?} start {:?} end {:?}", point, start, end);
            eprintln!("denominator {:?}", denominator);
            panic!();
        }

        assert!(denominator > 0.);

        let numerator: f64 = numerator as f64;

        numerator / denominator
    }
}

// Ramer–Douglas-Peucker line simplification algorithm
fn rdp(points: &[Point<i64>], epsilon: i64) -> Vec<Point<i64>>
{
    if points.is_empty() {
        return points.to_vec();
    }
    let mut dmax: f64 = 0.;
    let mut index: usize = 0;
    let mut distance: f64;

    for (i, _) in points.iter().enumerate().take(points.len() - 1).skip(1) {
        distance = point_line_distance(&points[i],
                                       &points[0],
                                       &*points.last().unwrap());
        if distance > dmax {
            index = i;
            dmax = distance;
        }
    }
    //let e = Fraction::new(e*e, T::one());
    let epsilon_f64: f64 = epsilon as f64;
    if dmax > epsilon_f64 {
        let mut intermediate = rdp(&points[..index + 1], epsilon);
        intermediate.pop();
        intermediate.extend_from_slice(&rdp(&points[index..], epsilon));
        intermediate
    } else {
        vec![*points.first().unwrap(), *points.last().unwrap()]
    }
}

pub fn simplify(geom: Geometry<i64>, epsilon: i64) -> Option<Geometry<i64>> {
    match geom {
        // Can't simplify a Point. let's hope this doesn't do a copy or memory move or something
        Geometry::Point(p) => Some(Geometry::Point(p)),
        Geometry::MultiPoint(p) => Some(Geometry::MultiPoint(p)),

        Geometry::LineString(ls) => simplify_linestring(ls, epsilon).map(|g| g.into()),
        Geometry::MultiLineString(mls) => simplify_multilinestring(mls, epsilon).map(|g| g.into()),
        Geometry::Polygon(p) => simplify_polygon(p, epsilon).map(|g| g.into()),
        Geometry::MultiPolygon(mp) => simplify_multipolygon(mp, epsilon).map(|g| g.into()),

        Geometry::GeometryCollection(_) => unimplemented!(),
    }
}

fn simplify_linestring(geom: LineString<i64>, epsilon: i64) -> Option<LineString<i64>> {
    Some(LineString(rdp(&geom.0, epsilon)))
}

fn simplify_multilinestring(geom: MultiLineString<i64>, epsilon: i64) -> Option<MultiLineString<i64>> {
    Some(MultiLineString(geom.0.into_iter().filter_map(|l| simplify_linestring(l, epsilon)).collect()))
}

fn simplify_polygon(geom: Polygon<i64>, epsilon: i64) -> Option<Polygon<i64>> {
    let Polygon{ exterior, interiors } = geom;
    Some(Polygon::new(simplify_linestring(exterior, epsilon).unwrap(), interiors.into_iter().filter_map(|l| simplify_linestring(l, epsilon)).collect()))
}

fn simplify_multipolygon(geom: MultiPolygon<i64>, epsilon: i64) -> Option<MultiPolygon<i64>> {
    Some(MultiPolygon(geom.0.into_iter().filter_map(|p| simplify_polygon(p, epsilon)).collect()))
}

#[cfg(test)]
mod test {
}
