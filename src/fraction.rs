use std::cmp::Ordering;
use std::ops::Neg;
use geo::*;
use std::fmt::Debug;
use std::ops::{Div,Rem,Mul,AddAssign,SubAssign};
use std::cmp::min;

fn gcd<T: CoordinateType+Rem>(mut a: T, mut b: T) -> T {
    while b != T::zero() {
        let t = b;
        b = a % b;
        a = t
    }

    a
}


fn reduce_common_factor<T: CoordinateType+Div<Output=T>+Rem>(mut a: T, mut b: T) -> (T, T) {
    let g = gcd(a, b);
    (a/g, b/g)
}


#[derive(Debug)]
pub struct Fraction<T: CoordinateType> {
    numerator: T,
    denominator: T,
}

impl<T: CoordinateType+Div<Output=T>+Rem> Fraction<T> {
    pub fn new(numerator: T, denominator: T) -> Self {
        let (numerator, denominator) = reduce_common_factor(numerator, denominator);
        Fraction{ numerator, denominator }
    }

    pub fn zero() -> Self {
        Self::new(T::zero(), T::one())
    }

    pub fn one() -> Self {
        Self::new(T::one(), T::one())
    }

}

impl PartialEq for Fraction<i64> {
    fn eq(&self, other: &Self) -> bool {
        let a = self.numerator;
        let b = self.denominator;
        let c = other.numerator;
        let d = other.denominator;

        if a == 0 || c == 0 {
            // one is zero, so this is true iff both are zero
            a == c
        } else if b == d {
            a == c
        } else {
            a.saturating_mul(d) == b.saturating_mul(c)
        }
    }
}

impl PartialEq for Fraction<i32> {
    fn eq(&self, other: &Self) -> bool {
        let a = self.numerator;
        let b = self.denominator;
        let c = other.numerator;
        let d = other.denominator;

        if a == 0 || c == 0 {
            // one is zero, so this is true iff both are zero
            a == c
        } else if b == d {
            a == c
        } else {
            a.saturating_mul(d) == b.saturating_mul(c)
        }
    }
}

impl PartialOrd for Fraction<i64> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        let mut a = self.numerator;
        let mut b = self.denominator;
        let mut c = other.numerator;
        let mut d = other.denominator;

        if b == d || a == 0 || c == 0 {
            // denominators are the same
            // or
            // one is zero, so check the numerator for which is neg
            a.partial_cmp(&c)
        } else {
            // In our use case, we can assume this is true
            // FIXME do we need this?
            //assert!(b != T::zero());
            //assert!(b > T::zero());
            //assert!(d != T::zero());
            //assert!(d > T::zero());

            //println!("{}:{} about to cmp {:?} {:?} {:?} {:?}", file!(), line!(), a, b, c, d);
            //assert!(a.checked_mul(d).is_some(), "{} L {}\nself {:?} other {:?}", file!(), line!(), self, other);
            //assert!(b.checked_mul(c).is_some(), "{} L {}\nself {:?} other {:?}", file!(), line!(), self, other);
            let ad = a.checked_mul(d);
            let bc = b.checked_mul(c);
            match (ad, bc) {
                (None, None) => Some(Ordering::Equal),
                (Some(ad), None) => Some(Ordering::Less),
                (None, Some(bc)) => Some(Ordering::Greater),
                (Some(ad), Some(bc)) => ad.partial_cmp(&bc),
            }
        }
    }
}

impl<T> Neg for Fraction<T>
    where T: CoordinateType + Neg<Output=T>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Fraction::new(-self.denominator, self.numerator)
    }
}

pub fn same_slope<T>(x: &Fraction<T>, y: &Fraction<T>) -> bool
    where T: CoordinateType + Neg<Output=T>+Rem+AddAssign+Debug
{
    let a = x.numerator;
    let b = x.denominator;
    let c = y.numerator;
    let d = y.denominator;

    if b == T::zero() && d == T::zero() {
        true
    } else if b == d {
        a == c
    } else {
        //println!("{}:{} before reduce {:?} {:?} {:?} {:?}", file!(), line!(), a, b, c, d);
        let (a, b) = reduce_common_factor(a, b);
        let (c, d) = reduce_common_factor(c, d);
        //println!("{}:{}  after reduce {:?} {:?} {:?} {:?}", file!(), line!(), a, b, c, d);
        //assert!(a.checked_mul(d).is_some(), "{} L {}\nself {:?} other {:?}", file!(), line!(), self, other);
        //assert!(b.checked_mul(c).is_some(), "{} L {}\nself {:?} other {:?}", file!(), line!(), self, other);
        a*d == b*c || -a*d == b*c
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn eq1() {
        assert_eq!(Fraction::new(0, 1), Fraction::new(0, 10));
    }

    #[test]
    fn cmp1() {
        assert_eq!(Fraction::new(0, 1).partial_cmp(&Fraction::new(0, 10)), Some(Ordering::Equal));
        assert_eq!(Fraction::new(2, 1).partial_cmp(&Fraction::new(0, 10)), Some(Ordering::Greater));
        assert_eq!(Fraction::new(2, 1).partial_cmp(&Fraction::new(10, 1)), Some(Ordering::Greater));
        assert_eq!(Fraction::new(-2, 1).partial_cmp(&Fraction::new(10, 1)), Some(Ordering::Less));
    }

    #[test]
    fn same_slope1() {
        assert!(same_slope(&Fraction::new(2, 2), &Fraction::new(1, 1)));
        assert!(same_slope(&Fraction::new(-2, 2), &Fraction::new(-1, 1)));
        assert!(same_slope(&Fraction::new(2, 2), &Fraction::new(-1, 1)));
        assert!(same_slope(&Fraction::new(-2, 2), &Fraction::new(1, 1)));
    }

    #[test]
    fn same_slope2() {
        assert!(same_slope(&Fraction::new(2, 0), &Fraction::new(1, 0)));
        assert!(same_slope(&Fraction::new(-2, 0), &Fraction::new(1, 0)));
        assert!(same_slope(&Fraction::new(2, 0), &Fraction::new(-1, 0)));
        assert!(same_slope(&Fraction::new(-2, 0), &Fraction::new(-1, 0)));
        assert!(same_slope(&Fraction::new(-1, 0), &Fraction::new(-2, 0)));
    }

    #[test]
    fn same_slope3() {
        assert!(same_slope(&Fraction::new(0, 1), &Fraction::new(0, 2)));
        assert!(!same_slope(&Fraction::new(10, 1), &Fraction::new(20, 1)));
        assert!(same_slope(&Fraction::new(10, 1), &Fraction::new(20, 2)));
    }

    #[test]
    fn gcd1() {
        assert_eq!(gcd(10, 5), 5);
        assert_eq!(gcd(5, 10), 5);
    }

    #[test]
    fn reduce_common_factor1() {
        // any minus sign seems to migrate to the denominator
        assert_eq!(reduce_common_factor(-12, 9), (4, -3));
        assert_eq!(reduce_common_factor(-12, -9), (4, 3));
        assert_eq!(reduce_common_factor(12, -9), (4, -3));
        assert_eq!(reduce_common_factor(12, 9), (4, 3));
    }

}
