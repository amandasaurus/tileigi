use geo::*;
use super::*;

use std::fmt::Debug;


fn clip_ring_to_border(ring: Cow<LineString<i32>>, border: &Border<i32>) -> Option<LineString<i32>> {
    //println!("\n\n\nBorder: {:?} Rings: {:?}", border, ring);

    // in our rings, the last point is the same as the first point, and this algorithm doesn't
    // support that.
    
    if ring.0.len() < 3 {
        //eprintln!("\n\n\nBorder: {:?} Rings: {:?}", border, ring);
        //eprintln!("{:?}", ring.0);
        // TODO something better
        return None;
    }
    assert!(ring.0.len() >= 3);

    // First we look if everything is all inside or all outside, and early return then, Then we
    // don't have to allocate a vec for the intersections
    let mut all_inside = true;
    let mut not_all_outside = false;

    let point_inside = is_inside(&ring.0[0], border);
    all_inside &= point_inside;
    not_all_outside |= point_inside;

    for (idx, point) in ring.0.iter().skip(1).enumerate() {
        let point_inside = is_inside(point, border);

        all_inside &= point_inside;
        not_all_outside |= point_inside;
    }

    let all_outside = ! not_all_outside;
    if all_inside {
        return Some(ring.into_owned());
    } else if all_outside {
        return None;
    }

    // No, we need to do some clipping

    let mut new_points = Vec::with_capacity(ring.0.len());

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
                new_points.push(intersection(&p1, &p2, border).into());
            }
        } else if is_inside(&p2, border) {
            new_points.push(intersection(&p1, &p2, border).into());
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

pub fn clip_polygon_to_border(poly: Cow<Polygon<i32>>, border: &Border<i32>) -> Option<Polygon<i32>> {

    match poly {
        Cow::Owned(poly) => {
            let Polygon{ exterior, interiors } = poly;
            let new_exterior = clip_ring_to_border(Cow::Owned(exterior), border);
            if let Some(new_exterior) = new_exterior {
                let interiors: Vec<_> = interiors.into_iter().filter_map(|i| clip_ring_to_border(Cow::Owned(i), border)).collect();
                Some(Polygon::new(new_exterior, interiors))
            } else {
                // if the exterior ring is totally outside, then don't bother any more
                None
            }
        },
        Cow::Borrowed(poly) => {
            let new_exterior = clip_ring_to_border(Cow::Borrowed(&poly.exterior), border);
            if let Some(new_exterior) = new_exterior {
                let interiors: Vec<_> = poly.interiors.iter().filter_map(|i| clip_ring_to_border(Cow::Borrowed(i), border)).collect();
                Some(Polygon::new(new_exterior, interiors))
            } else {
                // if the exterior ring is totally outside, then don't bother any more
                None
            }
        },
    }
}

pub fn clip_multipolygon_to_border(mp: Cow<MultiPolygon<i32>>, border: &Border<i32>) -> Option<MultiPolygon<i32>> {
    let polys: Vec<_> = match mp {
        Cow::Owned(mp) => mp.0.into_iter().filter_map(|p| clip_polygon_to_border(Cow::Owned(p), border)).collect(),
        Cow::Borrowed(mp) => mp.0.iter().filter_map(|p| clip_polygon_to_border(Cow::Borrowed(p), border)).collect(),
    };

    if polys.len() == 0 {
        None
    } else {
        Some(MultiPolygon(polys))
    }
}

pub fn clip_polygon_to_bbox(poly: Cow<Polygon<i32>>, bbox: &Bbox<i32>) -> Option<Polygon<i32>> {
    clip_polygon_to_border(poly, &Border::XMin(bbox.xmin))
           .and_then(|p| clip_polygon_to_border(Cow::Owned(p), &Border::XMax(bbox.xmax)))
           .and_then(|p| clip_polygon_to_border(Cow::Owned(p), &Border::YMin(bbox.ymin)))
           .and_then(|p| clip_polygon_to_border(Cow::Owned(p), &Border::YMax(bbox.ymax)))
}

pub fn clip_multipolygon_to_bbox(mp: Cow<MultiPolygon<i32>>, bbox: &Bbox<i32>) -> Option<MultiPolygon<i32>> {
    let polys: Vec<_> = match mp {
        Cow::Owned(mp) => mp.0.into_iter().filter_map(|p| clip_polygon_to_bbox(Cow::Owned(p), bbox)).collect(),
        Cow::Borrowed(mp) => mp.0.iter().filter_map(|p| clip_polygon_to_bbox(Cow::Borrowed(p), bbox)).collect(),
    };

    if polys.len() == 0 {
        None
    } else {
        Some(MultiPolygon(polys))
    }
}


mod test {
    use super::*;

    #[test]
    fn border_clip_simple_no_cut() {
        // Tests which should pass the polygon though directly
        assert_eq!(clip_ring_to_border(Cow::Borrowed(&vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into()), &Border::XMax(10)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );

        assert_eq!(clip_ring_to_border(Cow::Borrowed(&vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into()), &Border::XMin(0)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );

        assert_eq!(clip_ring_to_border(Cow::Borrowed(&vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into()), &Border::XMin(-1)),
            Some(vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into())
        );

        assert_eq!(clip_ring_to_border(Cow::Borrowed(&vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into()), &Border::XMin(10)),
            None
        );
    }

    #[test]
    fn border_clip_boxes() {
        // here polygons are going to be cut
        assert_eq!(clip_ring_to_border(Cow::Borrowed(&vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into()), &Border::XMin(1)),
            Some(vec![(1, 5), (5, 5), (5, 0), (1, 0), (1, 5)].into())
        );
        assert_eq!(clip_ring_to_border(Cow::Borrowed(&vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into()), &Border::YMin(1)),
            Some(vec![(0, 1), (0, 5), (5, 5), (5, 1), (0, 1)].into())
        );

        assert_eq!(clip_ring_to_border(Cow::Borrowed(&vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into()), &Border::XMax(2)),
            Some(vec![(0, 0), (0, 5), (2, 5), (2, 0), (0, 0)].into())
        );

        assert_eq!(clip_ring_to_border(Cow::Borrowed(&vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into()), &Border::YMax(2)),
            Some(vec![(0, 0), (0, 2), (5, 2), (5, 0), (0, 0)].into())
        );
    }

    #[test]
    fn border_clip_funny_shapes(){
        // Triangle pointing up
        assert_eq!(clip_ring_to_border(Cow::Borrowed(&vec![(0., 0.), (1., 5.), (2., 0.), (0., 0.)].into()), &Border::YMax(2.)),
            Some(vec![(0., 0.), (0.4, 2.), (1.6, 2.), (2., 0.), (0., 0.)].into())
        );
    }

    #[test]
    fn border_clip_polygon() {
        let poly = Polygon::new(
            vec![(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)].into(),
            vec![ vec![(1, 1), (1, 4), (4, 4), (4, 1), (1, 1)].into() ],
            );

        let new_poly = clip_polygon_to_border(Cow::Borrowed(&poly), &Border::XMax(3));

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
        let new_poly = clip_polygon_to_bbox(Cow::Borrowed(&poly), &bbox);

        // yes this is a degenerate polygon
        assert_eq!(new_poly, Some(Polygon::new(
            vec![(9, 9), (9, 5), (5, 5), (5, 9), (9, 9)].into(),
            vec![ vec![(5, 6), (6, 6), (6, 5), (5, 5), (5, 6)].into() ],
            ))
        );
    }

}
