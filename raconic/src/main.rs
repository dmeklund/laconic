#![allow(
    dead_code,
    unused_variables
)]

// mod spin;
mod plots;
mod utils;

use std::f64::consts::PI;
use std::cmp;

use itertools::Itertools;

// use nalgebra as na;
use num_complex::{Complex, Complex64};
use num_traits::{Pow};
use crate::utils::norm2;

const HBAR: f64 = 1.0;


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

    fn psi(&self, x: f64, t: f64) -> Complex64 {
        let energy = (self.n as f64*PI*HBAR/self.length).pow(2) / (2.0*self.mass);
        self.psi0(x) * Complex64::exp(-Complex::i() * energy * t / HBAR)
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

#[derive(Copy, Clone)]
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
                result = Some(self.get(index));
                index += 1
            } else {
                result = None
            }
            result
        })
    }

    fn get(&self, index: usize) -> f64{
        self.start + (index as f64)*(self.end - self.start) / (self.num as f64)
    }
}

#[derive(Clone)]
pub struct PositionState {
    pub xvals: LinearSpace,
    pub time: f64,
    pub vals: Vec<Complex64>,
    pub mass: f64,
}

// i*hbar*d/dt Psi(x, t) = [ -hbar^2/(2*m) d^2/dx^2 + V(x, t) ] Psi(x, t)
fn apply_time_step(state: PositionState, dt: f64) -> PositionState {
    let dx = 0.1;
    let d2_dx2 = eval_2nd_deriv(&state.vals, dx);
    let coeff = Complex64::from(-HBAR*HBAR / (2.0*state.mass));
    let kinetic: Vec<Complex64> = d2_dx2.iter().map(|x| x*coeff).collect();
    // let potential = na::DVector::from_element(kinetic.len(), Complex::from(0.0));
    let potential = vec![Complex::from(0.0); kinetic.len()];
    println!("Potential: {}", potential.len());
    println!("Kinetic: {}", kinetic.len());
    let d_dt: Vec<Complex64> = kinetic.iter().zip(potential.iter())
        .map(|(k, p)| k + p / (Complex::i() * HBAR))
        .collect();
    println!("d_dt: {:?}", &d_dt);
    PositionState {
        xvals: state.xvals,
        time: state.time + dt,
        vals: state.vals.iter().zip(d_dt.iter())
            .map(|(st, ddt)| st + ddt * Complex::from(dt))
            .collect(),
        mass: state.mass
    }
}

fn eval_2nd_deriv(yvals: &Vec<Complex64>, dx: f64) -> Vec<Complex64> {
    let mut result = Vec::with_capacity(yvals.len());
    println!("Evaling 2nd deriv, length {}", yvals.len());
    for ind in 0..2 {
        result.push(
            1.0 / (dx * dx) * (
                  2.0 * yvals[ind]
                - 5.0 * yvals[ind+1]
                + 4.0 * yvals[ind+2]
                - yvals[ind+3]
            )
        );
    }
    for ind in 2..yvals.len()-2 {
        result.push(
            1.0/(dx*dx) * (
                - 1.0/12.0*yvals[ind-2]
                + 4.0/3.0*yvals[ind-1]
                - 5.0/2.0*yvals[ind]
                + 4.0/3.0*yvals[ind+1]
                - 1.0/12.0*yvals[ind+2]
            )
        );
    }
    // if result.len() > 2 {
        for ind in yvals.len()-2..yvals.len() {
            result.push(1.0 / (dx * dx) * (
                2.0 * yvals[ind]
                    - 5.0 * yvals[ind - 1]
                    + 4.0 * yvals[ind - 2]
                    - yvals[ind - 3]
            ));
        }
    // }
    result
}

fn compare_timesteps(totaltime: f64, max_increments: i32) {
    for num_increments in 1..max_increments+1 {

    }
}

fn compare_single_point() {
    let boxbasis = BoxBasis::new(1, 1.0, 1.0);
    let xvals = LinearSpace { start: -0.5, end: 0.5, num: 10 };
    let psi0 = xvals
        .iter()
        .map(|x| boxbasis.psi0(x))
        .collect();
    println!("psi0: {:?}", &psi0);
    let initial_state = PositionState {
        xvals,
        time: 0.0,
        vals: psi0,
        mass: 1.0
    };

    let dt = 0.01;
    let next_psi = apply_time_step(initial_state, dt);
    let actual_psi: Vec<Complex64> = xvals
        .iter()
        .map(|x| boxbasis.psi(x, dt))
        .collect();
    println!("estimated: {:?}", &next_psi.vals);
    println!("actual: {:?}", &actual_psi);
    println!("difference: {:?}", next_psi.vals.iter().zip(actual_psi.iter())
        .map(|(est, act)| est - act));
}

fn plot_over_time() {
    let boxbasis = BoxBasis::new(1, 1.0, 1.0);
    let xvals = LinearSpace { start: -0.5, end: 0.5, num: 10 };
    let tvals = LinearSpace { start: 0.0, end: 10.0, num: 100 };
    let mut psivals = Vec::with_capacity(tvals.num);
    let mut estd = Vec::with_capacity(tvals.num);
    let mut laststate = None;
    for time_val in tvals.iter() {
        let psi: Vec<Complex64> = xvals.iter().map(|xval| { boxbasis.psi(xval, time_val) }).collect();
        println!("Psi norm: {}", norm2(psi.iter()));
        if time_val == 0.0 {
            let psi0 = psi.clone();
            laststate = Some(PositionState {
                xvals,
                time: 0.0,
                vals: psi0,
                mass: 1.0
            });
        } else {
            laststate = Some(apply_time_step(laststate.unwrap(), 0.01));
        }
        // println!("{:?}", psi);
        println!("Estd: {}", norm2(laststate.as_ref().unwrap().vals.iter()));
        psivals.push(psi);
        estd.push(laststate.as_ref().unwrap().clone());
    }

    plots::showplots(xvals, Vec::from([psivals]), 1024, 768).unwrap();

    // let root = BitMapBackend::new("images/3dtest.png", (640, 480)).into_drawing_area();
    // root.fill(&WHITE).unwrap();
    // let iter = psivals.iter().flat_map(|psi| psi.iter());
    // let (remin, remax, immin, immax) = complex_ranges(iter);
    // let mut chart = ChartBuilder::on(&root)
    //     .margin(20)
    //     .caption("Test 3d figure", ("sans-serif", 40))
    //     .build_cartesian_3d(xvals.start..xvals.end, remin..remax, immin..immax)
    //     .unwrap();
    // chart.configure_axes().draw().unwrap();
    // // println!("{:?}", psivals[0]);
    // for psi0 in psivals {
    //     chart.draw_series(LineSeries::new(
    //         (0..xvals.len()).map(|ind| (xvals.get(ind), psi0[ind].re, psi0[ind].im)),
    //         &RED
    //     )).unwrap();
    // }

    // let tvals = LinearSpace { start: 0.0, end: 2.0*PI, num: 100 };
    // let testvals: Vec<Complex64> = tvals.iter().map(|t| Complex64::new(0.0, -t).exp()).collect();
    // println!("{:?}", testvals);
    // chart.draw_series(LineSeries::new(
    //     (0..tvals.len()).map(|ind| (tvals.get(ind), testvals[ind].re, testvals[ind].im)),
    //     &RED
    // )).unwrap();
}


fn main() {
    compare_single_point();
    plot_over_time();
}
