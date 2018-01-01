use super::*;

#[test]
fn remap() {
//L 475 geom LineString(LineString([Point(Coordinate { x: -693741.39, y: 7049558.31 }), Point(Coordinate { x: -693886.45, y: 7049788.51 }), Point(Coordinate { x: -693905.81, y: 7049848.66 }), Point(Coordinate { x: -693923.15, y: 7049902.74 }), Point(Coordinate { x: -693956.59, y: 7050029.34 }), Point(Coordinate { x: -693985.26, y: 7050160.72 }), Point(Coordinate { x: -693997.2, y: 7050306.43 }), Point(Coordinate { x: -694009.15, y: 7050397.2 }), Point(Coordinate { x: -694022.23, y: 7050490.84 }), Point(Coordinate { x: -694037.39, y: 7050599.36 }), Point(Coordinate { x: -694166.75, y: 7051000.65 }), Point(Coordinate { x: -694400.88, y: 7051738.55 }), Point(Coordinate { x: -694427.16, y: 7051799.33 }), Point(Coordinate { x: -695009.99, y: 7052458.61 }), Point(Coordinate { x: -695055.37, y: 7052565.03 }), Point(Coordinate { x: -695093.59, y: 7052722.68 }), Point(Coordinate { x: -695103.15, y: 7053080.98 }), Point(Coordinate { x: -695072.09, y: 7054069.89 }), Point(Coordinate { x: -694990.43, y: 7054483.98 }), Point(Coordinate { x: -21474836.48, y: 20061906.38 })]))
//L 481 geom LineString(LineString([Point(Coordinate { x: 30499, y: 9711 }), Point(Coordinate { x: 30499, y: 9710 }), Point(Coordinate { x: 30498, y: 9710 }), Point(Coordinate { x: 30498, y: 9709 }), Point(Coordinate { x: 30498, y: 9708 }), Point(Coordinate { x: 30498, y: 9707 }), Point(Coordinate { x: 30497, y: 9704 }), Point(Coordinate { x: 30495, y: 9702 }), Point(Coordinate { x: 30495, y: 9701 }), Point(Coordinate { x: 30495, y: 9700 }), Point(Coordinate { x: 30495, y: 9696 }), Point(Coordinate { x: 30495, y: 9695 }), Point(Coordinate { x: -37469, y: -32848 })]))
//L 476 minx -10018754 maxx 0 miny 0 maxy 10018756 extent 32768

    let geom = Geometry::LineString(LineString(vec![Point(Coordinate { x: -693741.39, y: 7049558.31 }), Point(Coordinate { x: -693886.45, y: 7049788.51 }), Point(Coordinate { x: -693905.81, y: 7049848.66 }), Point(Coordinate { x: -693923.15, y: 7049902.74 }), Point(Coordinate { x: -693956.59, y: 7050029.34 }), Point(Coordinate { x: -693985.26, y: 7050160.72 }), Point(Coordinate { x: -693997.2, y: 7050306.43 }), Point(Coordinate { x: -694009.15, y: 7050397.2 }), Point(Coordinate { x: -694022.23, y: 7050490.84 }), Point(Coordinate { x: -694037.39, y: 7050599.36 }), Point(Coordinate { x: -694166.75, y: 7051000.65 }), Point(Coordinate { x: -694400.88, y: 7051738.55 }), Point(Coordinate { x: -694427.16, y: 7051799.33 }), Point(Coordinate { x: -695009.99, y: 7052458.61 }), Point(Coordinate { x: -695055.37, y: 7052565.03 }), Point(Coordinate { x: -695093.59, y: 7052722.68 }), Point(Coordinate { x: -695103.15, y: 7053080.98 }), Point(Coordinate { x: -695072.09, y: 7054069.89 }), Point(Coordinate { x: -694990.43, y: 7054483.98 }), Point(Coordinate { x: -21474836.48, y: 20061906.38 })]));

    let minx = -10018754.;
    let maxx = 0.;
    let miny = 0.;
    let maxy = 10018756.;
    let extent = 32768.;
    let new_geom = remap_geometry(geom, minx, maxx, miny, maxy, extent).unwrap();

    let expected = Geometry::LineString(LineString(vec![Point(Coordinate { x: 30499, y: 9711 }), Point(Coordinate { x: 30499, y: 9710 }), Point(Coordinate { x: 30498, y: 9710 }), Point(Coordinate { x: 30498, y: 9709 }), Point(Coordinate { x: 30498, y: 9708 }), Point(Coordinate { x: 30498, y: 9707 }), Point(Coordinate { x: 30497, y: 9704 }), Point(Coordinate { x: 30495, y: 9702 }), Point(Coordinate { x: 30495, y: 9701 }), Point(Coordinate { x: 30495, y: 9700 }), Point(Coordinate { x: 30495, y: 9696 }), Point(Coordinate { x: 30495, y: 9695 }), Point(Coordinate { x: -37469, y: -32848 })]));

    assert_eq!(expected, new_geom);

}
