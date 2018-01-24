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

impl<T: CoordinateType> PartialEq for Fraction<T> {
    fn eq(&self, other: &Self) -> bool {
        let a = self.numerator;
        let b = self.denominator;
        let c = other.numerator;
        let d = other.denominator;

        //if b == d {
        //    a == c
        //} else {
            a*d == b*c
        //}
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

impl<T> Neg for Fraction<T>
    where T: CoordinateType + Neg<Output=T>
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Fraction::new(-self.denominator, self.numerator)
    }
}

impl<T> Fraction<T>
    where T: CoordinateType + Neg<Output=T>
{
    pub fn same_slope(&self, other: &Fraction<T>) -> bool {
        self == other || self == &Fraction::new(-other.denominator, other.numerator)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn same_slope1() {
        assert!(Fraction::new(2, 2).same_slope(&Fraction::new(1, 1)));
        assert!(Fraction::new(-2, 2).same_slope(&Fraction::new(-1, 1)));
        assert!(Fraction::new(2, 2).same_slope(&Fraction::new(-1, 1)));
        assert!(Fraction::new(-2, 2).same_slope(&Fraction::new(1, 1)));
    }

    #[test]
    fn same_slope2() {
        assert!(Fraction::new(2, 0).same_slope(&Fraction::new(1, 0)));
        assert!(Fraction::new(-2, 0).same_slope(&Fraction::new(1, 0)));
        assert!(Fraction::new(2, 0).same_slope(&Fraction::new(-1, 0)));
        assert!(Fraction::new(-2, 0).same_slope(&Fraction::new(-1, 0)));

        assert!(Fraction::new(-1, 0).same_slope(&Fraction::new(-2, 0)));
    }

    #[test]
    fn same_slope3() {
        assert!(Fraction::new(0, 1).same_slope(&Fraction::new(0, 2)));
        assert!(!Fraction::new(10, 1).same_slope(&Fraction::new(20, 1)));
        assert!(Fraction::new(10, 1).same_slope(&Fraction::new(20, 2)));
    }

}
