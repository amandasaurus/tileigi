use geo::*;

use std::fmt::Debug;

// A border we want
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

fn intersection<T: CoordinateType>(p1: &Point<T>, p2: &Point<T>, border: &Border<T>) -> Point<T> {
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


fn clip_ring_to_border<T: CoordinateType+Debug>(ring: LineString<T>, border: &Border<T>) -> Option<LineString<T>> {
    //println!("\n\n\nBorder: {:?} Rings: {:?}", border, ring);
    let mut new_points = Vec::with_capacity(ring.0.len());

    // in our rings, the last point is the same as the first point, and this algorithm doesn't
    // support that.
    
    if ring.0.len() < 3 {
        //eprintln!("\n\n\nBorder: {:?} Rings: {:?}", border, ring);
        //eprintln!("{:?}", ring.0);
        // FIXME something better
        return None;
    }
    assert!(ring.0.len() >= 3);

    //if is_inside(&ring.0[0], border) {
    //    new_points.push(ring.0[0]);
    //}

    for points in ring.0.windows(2) {
        let p1 = points[0];
        let p2 = points[1];
        //println!("\n\np1 {:?}\np2 {:?}\nis_inside(p1) {} is_inside(p2) {}", p1, p2, is_inside(&p1, border), is_inside(&p2, border));
        
        if is_inside(&p1, border) {
            new_points.push(p1);
            if ! is_inside(&p2, border) {
                new_points.push(intersection(&p1, &p2, border));
            }
        } else if is_inside(&p2, border) {
            new_points.push(intersection(&p1, &p2, border));
        }
        //println!("all points {:?}", new_points);
    }

    if new_points.len() == 0 {
        None
    } else {
        let first_point = new_points[0].clone();
        new_points.push(first_point);
        Some(LineString(new_points))
    }

}

fn clip_polygon_to_border<T: CoordinateType+Debug>(poly: Polygon<T>, border: &Border<T>) -> Option<Polygon<T>> {
    let Polygon{ exterior, interiors } = poly;

    let new_exterior = clip_ring_to_border(exterior, border);
    if let Some(new_exterior) = new_exterior {
        let interiors: Vec<_> = interiors.into_iter().filter_map(|i| clip_ring_to_border(i, border)).collect();
        Some(Polygon::new(new_exterior, interiors))
    } else {
        // if the exterior ring is totally outside, then don't bother any more
        None
    }
}

fn clip_multipolygon_to_border<T: CoordinateType+Debug>(mp: MultiPolygon<T>, border: &Border<T>) -> Option<MultiPolygon<T>> {
    let polys: Vec<_> = mp.0.into_iter().filter_map(|p| clip_polygon_to_border(p, border)).collect();
    if polys.len() == 0 {
        None
    } else {
        Some(MultiPolygon(polys))
    }
}

pub fn clip_polygon_to_bbox<T: CoordinateType+Debug>(poly: Polygon<T>, bbox: &Bbox<T>) -> Option<Polygon<T>> {
    clip_polygon_to_border(poly, &Border::XMin(bbox.xmin))
           .and_then(|p| clip_polygon_to_border(p, &Border::XMax(bbox.xmax)))
           .and_then(|p| clip_polygon_to_border(p, &Border::YMin(bbox.ymin)))
           .and_then(|p| clip_polygon_to_border(p, &Border::YMax(bbox.ymax)))
}

pub fn clip_multipolygon_to_bbox<T: CoordinateType+Debug>(mp: MultiPolygon<T>, bbox: &Bbox<T>) -> Option<MultiPolygon<T>> {
    let polys: Vec<_> = mp.0.into_iter().filter_map(|p| clip_polygon_to_bbox(p, bbox)).collect();
    if polys.len() == 0 {
        None
    } else {
        Some(MultiPolygon(polys))
    }
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
        assert_eq!(clip_ring_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMax(10)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );

        assert_eq!(clip_ring_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMin(0)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );

        assert_eq!(clip_ring_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMin(-1)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );

        assert_eq!(clip_ring_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMin(10)),
            None
        );
    }

    #[test]
    fn border_clip_boxes() {
        // here polygons are going to be cut
        assert_eq!(clip_ring_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMin(1)),
            Some(vec![(1, 5), (5, 5), (5, 0), (1, 0), (1, 5)].into())
        );
        assert_eq!(clip_ring_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::YMin(1)),
            Some(vec![(0, 1), (0, 5), (5, 5), (5, 1), (0, 1)].into())
        );

        assert_eq!(clip_ring_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::XMax(2)),
            Some(vec![(0, 0), (0, 5), (2, 5), (2, 0), (0, 0)].into())
        );

        assert_eq!(clip_ring_to_border(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(), &Border::YMax(2)),
            Some(vec![(0, 0), (0, 2), (5, 2), (5, 0), (0, 0)].into())
        );
    }

    #[test]
    fn border_clip_funny_shapes(){
        // Triangle pointing up
        assert_eq!(clip_ring_to_border(vec![(0., 0.), (1., 5.), (2., 0.), (0., 0.)].into(), &Border::YMax(2.)),
            Some(vec![(0., 0.), (0.4, 2.), (1.6, 2.), (2., 0.), (0., 0.)].into())
        );
    }

    #[test]
    fn border_clip_polygon() {
        let poly = Polygon::new(
            vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(),
            vec![ vec![(1, 1), (1, 4), (4, 4), (4, 1), (1, 1)].into() ],
            );

        let new_poly = clip_polygon_to_border(poly, &Border::XMax(3));

        // yes this is a degenerate polygon
        assert_eq!(new_poly, Some(Polygon::new(
            vec![(0, 0), (0, 5), (3, 5), (3, 0), (0, 0)].into(),
            vec![ vec![(1, 1), (1, 4), (3, 4), (3, 1), (1, 1)].into() ],
            ))
        );
    }

    #[test]
    fn test_clip_polygon_to_bbox() {
        let poly = Polygon::new(
            vec![(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)].into(),
            vec![ vec![(4, 4), (4, 6), (6, 6), (6, 4), (4, 4)].into() ],
            );

        let bbox = Bbox{ xmin: 5, ymin: 5, xmax: 9, ymax: 9 };
        let new_poly = clip_polygon_to_bbox(poly, &bbox);

        // yes this is a degenerate polygon
        assert_eq!(new_poly, Some(Polygon::new(
            vec![(9, 9), (9, 5), (5, 5), (5, 9), (9, 9)].into(),
            vec![ vec![(5, 6), (6, 6), (6, 5), (5, 5), (5, 6)].into() ],
            ))
        );
    }

}
