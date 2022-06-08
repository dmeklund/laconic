use num_complex::{Complex, Complex64};

// pub fn norm2(vec: &Vec<Complex64>) -> f64 {
pub fn norm2<'a, T>(vals: T) -> f64
    where T: Iterator<Item = &'a Complex64>
{
    vals
        .into_iter()
        .map(|x| (x * x.conj()).re )
        .reduce(|accum, x| accum + x)
        .unwrap()
}
