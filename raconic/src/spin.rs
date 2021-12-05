use num_rational::Rational32;
use nalgebra as na;
use std::fmt::{Display, Formatter};
use num_traits::cast::ToPrimitive;

pub struct Spin {
    value: Rational32,
}


impl Spin {
    pub fn fractional(numer: i32, denom: i32) -> Spin {
        if denom != 1 && denom != 2 {
            panic!("Spin denominator can only be 1 or 2")
        }
        Spin { value: Rational32::new(numer, denom) }
    }

    pub fn whole(numer: i32) -> Spin {
        Spin { value: Rational32::new(numer, 1) }
    }

    pub fn numer(&self) -> i32 {
        *self.value.numer()
    }

    pub fn denom(&self) -> i32 {
        *self.value.denom()
    }

    pub fn double(&self) -> i32 {
        let double: Rational32 = Rational32::from_integer(2) * self.value;
        if !double.is_integer() {
            panic!("Double spin should be integral")
        }
        double.to_integer()
    }

    pub fn value(&self) -> Rational32 {
        self.value
    }
}


fn splus(spin: Spin) -> na::Matrix<f64, na::Dynamic, na::Dynamic, na::VecStorage<f64, na::Dynamic, na::Dynamic>> {
    let n: usize = (spin.double() + 1) as usize;
    let s = spin.value().to_f64().unwrap();
    let mut result = na::DMatrix::from_element(n, n, 0.0);
    // step from spin-1 to -spin-1 exclusive in increments of 1
    for index in 0..n-1 {
        let row = index;
        let col = row + 1;
        let m = s - 1.0 - (index as f64);
        result[(row, col)] = f64::sqrt(s*(s+1.0) - m*(m+1.0));
    }
    result
}

impl Display for Spin {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        if self.denom() == 1 {
            write!(f, "Spin {}", self.numer())
        } else {
            write!(f, "Spin {}/{}", self.numer(), self.denom())
        }
    }
}
