use geo::*;

fn outcode<T: CoordinateType>(x: T, y: T, bbox: &Bbox<T>) -> u8 {
    let mut result = 0;
    if x < bbox.xmin {
        result |= 1;        // left
    } else if x > bbox.xmax {
        result |= 2;        // right
    }
    if y < bbox.ymin {
        result |= 4;        // bottom
    } else if y > bbox.ymax {
        result |= 8         // top
    }

    result
}

pub fn clip<T: CoordinateType>(p0: &Point<T>, p1: &Point<T>, bbox: &Bbox<T>) -> Option<((T, T), (T, T))> {
    let mut x0 = p0.x();
    let mut y0 = p0.y();
    let mut x1 = p1.x();
    let mut y1 = p1.y();


    // If we don't assign to it, we get a compiler warning. If we do assign, we get a warning that
    // it's never used!
    #[allow(unused_assignments)]
    let mut x = T::zero();
    #[allow(unused_assignments)]
    let mut y = T::zero();

    let mut accept = false;

    let mut outcode0 = outcode(x0, y0, bbox);
    let mut outcode1 = outcode(x1, y1, bbox);

    loop {

        if outcode0 | outcode1 == 0 {
            accept = true;
            break;
        } else if outcode0 & outcode1 != 0 {
            break;
        } else {
            let outcode_of_an_outside_point = if outcode0 != 0 { outcode0 } else { outcode1 };

            // Unlike the standard algorithm we have Y positive going down. So ymin & ymax are
            // swapped here
            if outcode_of_an_outside_point & 8 != 0 {           // top
                x = x0 + (x1 - x0) * (bbox.ymax - y0) / (y1 - y0);
				y = bbox.ymax;
			} else if outcode_of_an_outside_point & 4 != 0 {        // bottom
				x = x0 + (x1 - x0) * (bbox.ymin - y0) / (y1 - y0);
				y = bbox.ymin;
			} else if outcode_of_an_outside_point & 2 != 0 {    // right
				y = y0 + (y1 - y0) * (bbox.xmax - x0) / (x1 - x0);
				x = bbox.xmax;
			} else if outcode_of_an_outside_point & 1 != 0 {     // left
				y = y0 + (y1 - y0) * (bbox.xmin - x0) / (x1 - x0);
				x = bbox.xmin;
			} else {
                unreachable!();
            }

            if outcode_of_an_outside_point == outcode0 {
                x0 = x;
                y0 = y;
                outcode0 = outcode(x0, y0, &bbox);
            } else {
                x1 = x;
                y1 = y;
                outcode1 = outcode(x1, y1, &bbox);
            }
        }
    }

    if accept {
        Some(((x0, y0), (x1, y1)))
    } else {
        None
    }
}

mod test {
    use super::*;

    #[test]
    fn test_clip() {
        let bbox = Bbox{ xmin: 0., ymin: 0., xmax: 4096., ymax: 4096.};

        let p1 = Point::new(10., 10.);
        let p2 = Point::new(20., 20.);
        let p3 = Point::new(10., 5000.);
        let p4 = Point::new(5000., 5000.);
        let p5 = Point::new(10., 1000.);
        let p6 = Point::new(5000., 6000.);

        // entirely inside
        assert_eq!(clip(&p1, &p2, &bbox), Some(((10., 10.), (20., 20.))));

        // entirely outside
        assert_eq!(clip(&p4, &p6, &bbox), None);

        // extending over the top, left bottom and right
        assert_eq!(clip(&Point::new(10., 100.), &Point::new(10., 10000.), &bbox), Some(((10., 100.), (10., 4096.))));
        assert_eq!(clip(&Point::new(10., -100.), &Point::new(10., 100.), &bbox), Some(((10., 0.), (10., 100.))));
        assert_eq!(clip(&Point::new(-100., 10.), &Point::new(100., 10.), &bbox), Some(((0., 10.), (100., 10.))));
        assert_eq!(clip(&Point::new(100., 10.), &Point::new(10000., 10.), &bbox), Some(((100., 10.), (4096., 10.))));

        assert_eq!(clip(&Point::new(-1000., 2000.), &Point::new(2000., -1000.), &bbox), Some(((0., 1000.), (1000., 0.))));

    }

}
