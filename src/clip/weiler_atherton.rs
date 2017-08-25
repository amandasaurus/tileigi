use geo::*;

// x = (x1 - x0)t + x0
// y = (y1 - y0)t + y0

#[derive(Clone,Copy)]
enum Border<T: CoordinateType> {
    XMin(T),
    XMax(T),
    YMin(T),
    YMax(T),
}

fn is_outside<T: CoordinateType>(x: T, y: T, border: &Border<T>) -> bool {
    match *border {
        Border::XMin(xmin) => x < xmin,
        Border::XMax(xmax) => x > xmax,
        Border::YMin(ymin) => y < ymin,
        Border::YMax(ymax) => y > ymax,
    }
}


fn clip_to_border<T: CoordinateType>(ring: LineString<T>, border: &Border<T>) -> Option<LineString<T>> {
    let mut new_list = Vec::with_capacity(ring.0.len());
    let mut current_point: Option<&Point<T>> = None;
    let mut currently_outside = false;

    for point in ring.0.into_iter() {
        let x = point.x();
        let y = point.y();
        let point_is_outside = is_outside(x, y, border);
        if point_is_outside {
        } else {

            new_list.push(point);
        }
    }

    if new_list.len() == 0 {
        None
    } else {
        Some(LineString(new_list))
    }

}

pub fn clip<T: CoordinateType>(poly: Polygon<T>, bbox: &Bbox<T>) -> Option<MultiPolygon<T>> {
    None
    
}

mod test {
    use super::*;

    #[test]
    fn border_clip() {
        assert_eq!(clip_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMax(10)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );
    }
}
