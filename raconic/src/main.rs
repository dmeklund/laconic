use std::cmp::Ordering;
use std::fmt::{Display, Formatter};
use std::f64::consts::PI;

use nalgebra as na;
use nalgebra::DMatrix;
use num_complex::{Complex, Complex64};
use num_rational::Rational32;
use num_traits::{Pow, cast::ToPrimitive};

const HBAR: f64 = 1.0;

pub struct Spin {
    value: Rational32,
}

trait BasisFunction {
    fn psi0(&self, x: f64) -> Complex64;
}

pub struct BoxBasis {
    n: i32,
    length: f64,
    mass: f64,
}

impl BoxBasis {
    fn new(n: i32, length: f64, mass: f64) -> BoxBasis {
        BoxBasis { n, length, mass }
    }

    fn psi(&self, x: f64, t: f64) -> num_complex::Complex64 {
        let energy = (self.n as f64*PI*HBAR/self.length).pow(2) / (2.0*self.mass);
        self.psi0(x) * (-Complex::i() * energy * t / HBAR)
    }
}

impl BasisFunction for BoxBasis {
    fn psi0(&self, x: f64) -> Complex64 {
        // psi_n(x, t=0) = sqrt(2/L) * sin(n*pi/L*(x - L/2))
        //   for -L/2 < x < L /2; 0 otherwise
        match x {
            x if (x < -self.length / 2. || x > self.length / 2.) => Complex64::from(0.0),
            _ => Complex64::from((2.0 / self.length).sqrt() * (
                (self.n as f64) * PI / self.length * (x - self.length / 2.0)
            ).sin())
        }
    }
}

pub struct LinearSpace {
    start: f64,
    end: f64,
    num: usize,
}

impl LinearSpace {
    fn len(&self) -> usize {
        self.num
    }

    fn iter(&self) -> impl Iterator<Item = f64> + '_ {
        let mut index = 0;
        std::iter::from_fn(move || {
            let result;
            if index < self.num {
                result = Some(self.start + index as f64*(self.end - self.start)/self.num as f64);
                index += 1
            } else {
                result = None
            }
            result
        })
    }
}

pub struct PositionState {
    pub xvals: LinearSpace,
    pub time: f64,
    pub vals: na::DVector<Complex64>,
    pub mass: f64,
}

// i*hbar*d/dt Psi(x, t) = [ -hbar^2/(2*m) d^2/dx^2 + V(x, t) ] Psi(x, t)
fn apply_time_step(state: PositionState, dt: f64) -> PositionState {
    let dx = 0.1;
    let d2_dx2 = eval_2nd_deriv(&state.vals, dx);
    let coeff = Complex64::from(-HBAR*HBAR / (2.0*state.mass));
    let kinetic = d2_dx2 * coeff;
    let potential = na::DVector::from_element(kinetic.len(), Complex::from(0.0));
    let d_dt = (kinetic + potential) / (Complex::i() * HBAR);
    PositionState {
        xvals: state.xvals,
        time: state.time + dt,
        vals: state.vals + d_dt * Complex::from(dt),
        mass: state.mass
    }
}

fn eval_2nd_deriv(yvals: &na::DVector<Complex64>, dx: f64) -> na::DVector<Complex64> {
    let mut result = na::DVector::from_element(yvals.len(), Complex64::new(0.0, 0.0));
    for ind in 2..result.len()-2 {
        result[ind] = 1.0/(dx*dx) * (
            - 1.0/12.0*yvals[ind-2]
            + 4.0/3.0*yvals[ind-1]
            - 5.0/2.0*yvals[ind]
            + 4.0/3.0*yvals[ind+1]
            - 1.0/12.0*yvals[ind+2]
        );
    }
    result
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

fn main() {
    let boxbasis = BoxBasis {
        n: 1,
        length: 1.0,
        mass: 1.0,
    };
    let xvals = LinearSpace { start: 0.0, end: 10.0, num: 10 };
    let psi0 = na::DVector::from_iterator(
        xvals.len(),
        xvals.iter().map(|x| { boxbasis.psi0(x) }),
    );
    let initial_state = PositionState {
        xvals,
        time: 0.0,
        vals: psi0,
        mass: 1.0
    };
    let dt = 1.0;
    // let next_psi = apply_time_step(psi0)
}
