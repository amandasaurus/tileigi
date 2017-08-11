use geo::*;

// x = (x1 - x0)t + x0
// y = (y1 - y0)t + y0

#[derive(Debug,Clone,Copy)]
enum Border<T: CoordinateType> {
    XMin(T),
    XMax(T),
    YMin(T),
    YMax(T),
}

fn is_inside<T: CoordinateType>(p: &Point<T>, border: &Border<T>) -> bool {
    match *border {
        Border::XMin(xmin) => p.x() >= xmin,
        Border::XMax(xmax) => p.x() <= xmax,
        Border::YMin(ymin) => p.y() >= ymin,
        Border::YMax(ymax) => p.y() <= ymax,
    }
}

fn intersection<T: CoordinateType+::std::fmt::Debug>(p1: &Point<T>, p2: &Point<T>, border: &Border<T>) -> Point<T> {
    //println!("\nstart of intersection");
    let x1 = p1.x();
    let y1 = p1.y();
    let x2 = p2.x();
    let y2 = p2.y();
    match *border {
        Border::XMin(xmin) => {
            let x = xmin;
            //println!("p1 {:?} p1 {:?} border {:?}", p1, p2, border);
            assert!(p1.x() != p2.x());
            let y = (y2-y1)*(x - x1)/(x2 - x1) + y1;
            Point::new(x, y)
        },
        Border::XMax(xmax) => {
            let x = xmax;
            assert!(p1.x() != p2.x());
            let y = (y2-y1)*(x - x1)/(x2 - x1) + y1;
            Point::new(x, y)
        },
        Border::YMin(ymin) => {
            let y = ymin;
            assert!(p1.y() != p2.y());
            let x = (x2-x1)*(y - y1)/(y2 - y1) + x1;
            Point::new(x, y)
        },
        Border::YMax(ymax) => {
            let y = ymax;
            assert!(p1.y() != p2.y());
            let x = (x2-x1)*(y - y1)/(y2 - y1) + x1;
            Point::new(x, y)
        },
    }
}


fn clip_to_border<T: CoordinateType+::std::fmt::Debug>(ring: LineString<T>, border: &Border<T>) -> Option<LineString<T>> {
    println!("\n\n\nBorder: {:?} Rings: {:?}", border, ring);
    let mut new_points = Vec::with_capacity(ring.0.len());

    // in our rings, the last point is the same as the first point, and this algorithm doesn't
    // support that.

    let mut next_idx = ring.0.len()-2;
    for i in 0..ring.0.len()-1 {
        assert_ne!(i, next_idx);
        let p1 = ring.0[i];
        let p2 = ring.0[next_idx];
        println!("\n\ni {} p1 {:?}\nnext_idx {} {:?}\nis_inside(p1) {} is_inside(p2) {}", i, p1, next_idx, p2, is_inside(&p1, border), is_inside(&p2, border));

        if is_inside(&p1, border) {
            if ! is_inside(&p2, border) {
                new_points.push(intersection(&p1, &p2, border));
            }
            new_points.push(p1);
        } else if is_inside(&p2, border) {
            new_points.push(intersection(&p1, &p2, border));
        }
        println!("all points {:?}", new_points);
        next_idx = i;
    }

    if new_points.len() == 0 {
        None
    } else {
        let first_point = new_points[0].clone();
        new_points.push(first_point);
        Some(LineString(new_points))
    }

}

pub fn clip<T: CoordinateType>(poly: Polygon<T>, bbox: &Bbox<T>) -> Option<MultiPolygon<T>> {
    None
    
}


mod test {
    use super::*;

    #[test]
    fn test_intersection() {
        assert_eq!(intersection(&Point::new(0, 0), &Point::new(10, 0), &Border::XMax(5)), Point::new(5, 0));
        assert_eq!(intersection(&Point::new(0, 0), &Point::new(10, 0), &Border::XMin(5)), Point::new(5, 0));

        assert_eq!(intersection(&Point::new(0, 0), &Point::new(0, 10), &Border::YMax(5)), Point::new(0, 5));
        assert_eq!(intersection(&Point::new(0, 0), &Point::new(0, 10), &Border::YMin(5)), Point::new(0, 5));
    }

    #[test]
    fn border_clip_simple_no_cut() {
        // Tests which should pass the polygon though directly
        assert_eq!(clip_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMax(10)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );

        assert_eq!(clip_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMin(0)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );

        assert_eq!(clip_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMin(-1)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );

        assert_eq!(clip_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMin(10)),
            None
        );
    }

    #[test]
    fn border_clip() {
        // here polygons are going to be cut
        assert_eq!(clip_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMin(1)),
            Some(vec![(1, 0), (1, 5), (5, 5), (5, 0), (1, 0)].into())
        );
        assert_eq!(clip_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::YMin(1)),
            Some(vec![(0, 1), (0, 5), (5, 5), (5, 1), (0, 1)].into())
        );

        assert_eq!(clip_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMax(2)),
            Some(vec![(0, 0), (0, 5), (2, 5), (2, 0), (0, 0)].into())
        );
    }
}
