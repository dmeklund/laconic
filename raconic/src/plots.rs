use std::error::Error;
use std::time::SystemTime;
use std::borrow::{Borrow, BorrowMut};

use minifb::{Key, Window, WindowOptions};
use num_complex::Complex64;
use plotters::prelude::*;
use plotters_bitmap::bitmap_pixel::BGRXPixel;
use plotters_bitmap::BitMapBackend;
use crate::LinearSpace;


const SAMPLE_RATE: f64 = 10_000.0;
const FRAME_RATE: f64 = 30.0;

struct BufferWrapper(Vec<u32>);

impl Borrow<[u8]> for BufferWrapper {
    fn borrow(&self) -> &[u8] {
        unsafe {
            std::slice::from_raw_parts(
                self.0.as_ptr() as *const u8,
                self.0.len() * 4
            )
        }
    }
}

impl BorrowMut<[u8]> for BufferWrapper {
    fn borrow_mut(&mut self) -> &mut [u8] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.0.as_mut_ptr() as *mut u8,
                self.0.len() * 4
            )
        }
    }
}

impl Borrow<[u32]> for BufferWrapper {
    fn borrow(&self) -> &[u32] {
        self.0.as_slice()
    }
}

impl BorrowMut<[u32]> for BufferWrapper {
    fn borrow_mut(&mut self) -> &mut [u32] {
        self.0.as_mut_slice()
    }
}

pub fn complex_ranges<'a, I>(iter: I) -> (f64, f64, f64, f64)
    where
        I: Iterator<Item = &'a Complex64>
{
    let mut remin: Option<f64> = None;
    let mut remax: Option<f64> = None;
    let mut immin: Option<f64> = None;
    let mut immax: Option<f64> = None;
    for val in iter {
        if remin == None {
            remin = Some(val.re);
            remax = Some(val.re);
            immin = Some(val.im);
            immax = Some(val.im);
        } else {
            if val.re < remin.unwrap() {
                remin = Some(val.re);
            }
            if val.re > remax.unwrap() {
                remax = Some(val.re);
            }
            if val.im < immin.unwrap() {
                immin = Some(val.im);
            }
            if val.im > immax.unwrap() {
                immax = Some(val.im);
            }
        }
    }
    (remin.unwrap(), remax.unwrap(), immin.unwrap(), immax.unwrap())
}


pub fn showplots(
        xvals: LinearSpace,
        values: Vec<Vec<Vec<Complex64>>>,
        width: usize,
        height: usize
) -> Result<(), Box<dyn Error>> {
    let mut buf = BufferWrapper(vec![0u32; width*height]);
    let mut window = Window::new(
        "My Window",
        width,
        height,
        WindowOptions::default()
    )?;
    let root = BitMapBackend::<BGRXPixel>::with_buffer_and_format(
        buf.borrow_mut(),
        (width as u32, height as u32)
    )?.into_drawing_area();
    root.fill(&BLACK)?;
    let iter = values[0].iter().flat_map(|psi| psi.iter());
    let (remin, remax, immin, immax) = complex_ranges(iter);

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .set_all_label_area_size(30)
        .build_cartesian_3d(xvals.start..xvals.end, remin..remax, immin..immax)?;
    chart
        .configure_axes()
        .label_style(("san-serif", 15).into_font().color(&GREEN))
        .axis_panel_style(&GREEN)
        .draw()?;
    let cs = chart.into_chart_state();
    drop(root);

    // let mut data = VecDeque::new();
    let start_ts = SystemTime::now();
    let mut last_flushed = 0.0;
    let mut dataind = 0;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        let epoch = SystemTime::now()
            .duration_since(start_ts)
            .unwrap()
            .as_secs_f64();
        if epoch - last_flushed > 1.0 / FRAME_RATE {
            println!("Actual frame rate: {}", 1.0 / (epoch - last_flushed));
            let root = BitMapBackend::<BGRXPixel>::with_buffer_and_format(
                buf.borrow_mut(),
                (width as u32, height as u32)
            )?.into_drawing_area();
            let mut chart = cs.clone().restore(&root);
            chart.plotting_area().fill(&BLACK)?;
            chart
                .configure_axes()
                .bold_grid_style(&GREEN.mix(0.2))
                .light_grid_style(&TRANSPARENT)
                .draw()?;
            for series in &values {
                chart.draw_series(LineSeries::new(
                    (0..xvals.len()).map(|ind| (xvals.get(ind), series[dataind][ind].re, series[dataind][ind].im)),
                    &GREEN
                ))?;
            }
            drop(root);
            drop(chart);
            window.update_with_buffer(buf.borrow(), width, height)?;
            last_flushed = epoch;
            dataind = (dataind + 1) % values[0].len();
        } else {
            let sleeptime = 1.0 / FRAME_RATE - (epoch - last_flushed);
            println!("Sleeping {}s", sleeptime);
            std::thread::sleep(std::time::Duration::from_secs_f64(sleeptime));
        }
    }
    Ok(())
}
