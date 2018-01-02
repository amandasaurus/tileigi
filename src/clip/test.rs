use super::*;
use super::super::validity;

//#[test]
//fn clip_point_to_border() {
//    assert_eq!(clip_to_border(&Point::new(0, 0).into(), &Border::XMin(0)), Some(Point::new(0, 0).into()));
//    assert_eq!(clip_to_border(&Point::new(0, 0).into(), &Border::XMin(10)), None);
//}
//
//#[test]
//fn clip_point() {
//    let bbox = Bbox{ xmin: 0., ymin: 0., xmax: 4096., ymax: 4096.};
//
//    assert_eq!(clip_to_bbox(&Point::new(0., 0.).into(), &bbox), Some(Point::new(0., 0.).into()));
//    assert_eq!(clip_to_bbox(&Point::new(-1., -1.).into(), &bbox), None);
//    assert_eq!(clip_to_bbox(&Point::new(4000., 0.).into(), &bbox), Some(Point::new(4000., 0.).into()));
//    assert_eq!(clip_to_bbox(&Point::new(5000., 0.).into(), &bbox), None);
//}
//
//fn mklinestring<T: CoordinateType>(points: Vec<(T, T)>) -> LineString<T> {
//    points.into_iter().collect::<LineString<T>>()
//}
//
//#[test]
//fn test_intersections() {
//    let ls = mklinestring(vec![(0, 0), (0, 5), (5, 5)]);
//
//    // Simple all inside or all outside cases.
//    assert_eq!(calculate_intersections(&ls, &Border::XMin(0)), LineBorderIntersection::AllInside);
//    assert_eq!(calculate_intersections(&ls, &Border::YMin(0)), LineBorderIntersection::AllInside);
//    assert_eq!(calculate_intersections(&ls, &Border::XMax(10)), LineBorderIntersection::AllInside);
//    assert_eq!(calculate_intersections(&ls, &Border::YMax(10)), LineBorderIntersection::AllInside);
//
//    assert_eq!(calculate_intersections(&ls, &Border::XMin(10)), LineBorderIntersection::AllOutside);
//    assert_eq!(calculate_intersections(&ls, &Border::YMin(10)), LineBorderIntersection::AllOutside);
//    assert_eq!(calculate_intersections(&ls, &Border::XMax(-1)), LineBorderIntersection::AllOutside);
//    assert_eq!(calculate_intersections(&ls, &Border::YMax(-1)), LineBorderIntersection::AllOutside);
//
//    // actual intersections
//    assert_eq!(calculate_intersections(&ls, &Border::XMax(2)), LineBorderIntersection::Intersections(vec![IntersectionOption::Inside, IntersectionOption::Exit((2, 5)), IntersectionOption::Outside]));
//    assert_eq!(calculate_intersections(&ls, &Border::YMax(2)), LineBorderIntersection::Intersections(vec![IntersectionOption::Exit((0, 2)), IntersectionOption::Outside, IntersectionOption::Outside]));
//    assert_eq!(calculate_intersections(&ls, &Border::XMin(2)), LineBorderIntersection::Intersections(vec![IntersectionOption::Outside, IntersectionOption::Entry((2, 5)), IntersectionOption::Inside]));
//    assert_eq!(calculate_intersections(&ls, &Border::YMin(2)), LineBorderIntersection::Intersections(vec![IntersectionOption::Entry((0, 2)), IntersectionOption::Inside, IntersectionOption::Inside]));
//
//    // weird edge cases
//    assert_eq!(calculate_intersections(&ls, &Border::XMax(0)), LineBorderIntersection::Intersections(vec![IntersectionOption::Inside, IntersectionOption::Inside, IntersectionOption::Outside]));
//    assert_eq!(calculate_intersections(&ls, &Border::YMax(5)), LineBorderIntersection::AllInside);
//    assert_eq!(calculate_intersections(&ls, &Border::YMax(0)), LineBorderIntersection::Intersections(vec![IntersectionOption::Inside, IntersectionOption::Outside, IntersectionOption::Outside]));
//}
//
//#[test]
//fn clip_linestring_to_border1() {
//    let ls: Geometry<i32> = mklinestring(vec![(0, 0), (0, 5), (5, 5)]).into();
//    let ls2: Geometry<i32> = ls.clone();
//
//    // Simple all inside or all outside cases.
//    assert_eq!(clip_to_border(&ls, &Border::XMin(10)), None);
//    assert_eq!(clip_to_border(&ls, &Border::XMin(-1)), Some(ls.clone()));
//    assert_eq!(clip_to_border(&ls, &Border::YMin(0)), Some(ls.clone()));
//    assert_eq!(clip_to_border(&ls, &Border::YMax(-10)), None);
//
//    // where it intersects
//    assert_eq!(clip_to_border(&ls, &Border::XMax(1)), Some(mklinestring(vec![(0, 0), (0, 5), (1, 5)]).into()));
//    assert_eq!(clip_to_border(&ls, &Border::XMax(10)), Some(mklinestring(vec![(0, 0), (0, 5), (5, 5)]).into()));
//    assert_eq!(clip_to_border(&ls, &Border::XMax(0)), Some(mklinestring(vec![(0, 0), (0, 5)]).into()));
//
//    // where it would return just one point
//    assert_eq!(clip_to_border(&ls, &Border::YMax(0)), None);
//}
//
//#[test]
//fn clip_linestring_to_border2() {
//    let ls: Geometry<i32> = mklinestring(vec![(0, 0), (0, 5), (5, 5), (5, 0)]).into();
//    let ls2: Geometry<i32> = ls.clone();
//
//    // Simple all inside or all outside cases.
//    assert_eq!(clip_to_border(&ls, &Border::XMin(10)), None);
//    assert_eq!(clip_to_border(&ls, &Border::XMin(-1)), Some(ls.clone()));
//    assert_eq!(clip_to_border(&ls, &Border::YMin(0)), Some(ls.clone()));
//    assert_eq!(clip_to_border(&ls, &Border::YMax(-10)), None);
//
//    // where it intersects
//    assert_eq!(clip_to_border(&ls, &Border::YMax(1)), Some(MultiLineString(vec![mklinestring(vec![(0, 0), (0, 1)]), mklinestring(vec![(5, 1), (5, 0)])]).into()));
//    assert_eq!(clip_to_border(&ls, &Border::XMax(1)), Some(mklinestring(vec![(0, 0), (0, 5), (1, 5)]).into()));
//    assert_eq!(clip_to_border(&ls, &Border::XMax(10)), Some(mklinestring(vec![(0, 0), (0, 5), (5, 5), (5, 0)]).into()));
//    assert_eq!(clip_to_border(&ls, &Border::XMax(0)), Some(mklinestring(vec![(0, 0), (0, 5)]).into()));
//
//}
//
//#[test]
//fn clip_multipoint_to_border() {
//    let mp: Geometry<_> = MultiPoint(vec![Point::new(0, 0), Point::new(5, 5)]).into();
//    assert_eq!(clip_to_border(&mp, &Border::XMin(0)), Some(mp.clone()));
//    assert_eq!(clip_to_border(&mp, &Border::XMin(2)), Some(Geometry::MultiPoint(MultiPoint(vec![Point::new(5, 5)]))));
//    assert_eq!(clip_to_border(&mp, &Border::XMin(10)), None);
//}
//
////#[test]
//fn clip_linestring() {
//    let bbox = Bbox{ xmin: 0., ymin: 0., xmax: 4096., ymax: 4096.};
//
//    let p1 = Point::new(10., 10.);
//    let p2 = Point::new(20., 20.);
//    let p3 = Point::new(10., 5000.);
//    let p4 = Point::new(5000., 5000.);
//    let p5 = Point::new(10., 1000.);
//
//    assert_eq!(clip_to_bbox(&LineString(vec![p1, p2]).into(), &bbox), Some(LineString(vec![p1, p2]).into()));
//    assert_eq!(clip_to_bbox(&LineString(vec![p4, p4]).into(), &bbox), None);
//    assert_eq!(clip_to_bbox(&LineString(vec![p1, p5]).into(), &bbox), Some(LineString(vec![p1, p5]).into()));
//    assert_eq!(clip_to_bbox(&LineString(vec![p1, p3]).into(), &bbox), Some(LineString(vec![p1, Point::new(10., 4096.)]).into()));
//
//    assert_eq!(clip_to_bbox(&LineString(vec![p1, p2, p5]).into(), &bbox),
//        Some(MultiLineString(vec![
//            LineString(vec![p1, p2]),
//            LineString(vec![p2, p5]),
//                            ]).into()));
//
//    assert_eq!(clip_to_bbox(&LineString(vec![p1, p2, p5, p4]).into(), &bbox),
//        Some(MultiLineString(vec![
//            LineString(vec![p1, p2]),
//            LineString(vec![p2, p5]),
//            LineString(vec![p5, Point::new(3872.26, 4096.)]),
//                            ]).into()));
//}
//
////#[test]
//fn clip_multilinestring() {
//    let bbox = Bbox{ xmin: 0., ymin: 0., xmax: 4096., ymax: 4096.};
//
//    let p1 = Point::new(10., 10.);
//    let p2 = Point::new(20., 20.);
//    let p3 = Point::new(10., 5000.);
//    let p4 = Point::new(5000., 5000.);
//    let p5 = Point::new(10., 1000.);
//
//    assert_eq!(clip_to_bbox(
//            &MultiLineString(vec![
//                LineString(vec![p1, p2, p5]),
//                LineString(vec![p1, p2, p5, p4])
//            ]).into(), &bbox),
//        Some(MultiLineString(vec![
//            LineString(vec![p1, p2]),
//            LineString(vec![p2, p5]),
//            LineString(vec![p1, p2]),
//            LineString(vec![p2, p5]),
//            LineString(vec![p5, Point::new(3872.26, 4096.)]),
//                            ]).into()));
//}

#[test]
fn result_valid_geom() {
    let geom = Geometry::Polygon(Polygon { exterior: LineString(vec![Point(Coordinate { x: 31565, y: 20875 }), Point(Coordinate { x: 31615, y: 20887 }), Point(Coordinate { x: 31633, y: 20819 }), Point(Coordinate { x: 31593, y: 20822 }), Point(Coordinate { x: 31585, y: 20808 }), Point(Coordinate { x: 31584, y: 20850 }), Point(Coordinate { x: 31565, y: 20875 })]), interiors: vec![] });
    assert!(validity::is_valid(&geom));
    let metatile =  Metatile::new(8, 4, 0, 0).unwrap();
    let buffer = 0;
    for (t, g) in clip_geometry_to_tiles(&metatile, geom, buffer).into_iter() {
        match g {
            Some(mut g) => {
                assert!(validity::is_valid(&g), "Invalid geometry {:?}", g);
            },
            _ => {}
        }
    }
}
