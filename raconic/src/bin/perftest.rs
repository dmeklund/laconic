const WIDTH: usize = 2400;
const HEIGHT: usize = 1600;

fn main() {

    let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];

    let mut window = match minifb::Window::new(
        "Test - ESC to exit",
        WIDTH,
        HEIGHT,
        minifb::WindowOptions::default()) {
        Ok(win) => win,
        Err(err) => {
            println!("Unable to create window {}", err);
            return;
        }
    };


    let mut last_time = std::time::Instant::now();
    while window.is_open() && !window.is_key_down(minifb::Key::Escape) {
        window.update_with_buffer(&buffer, WIDTH, HEIGHT);
        println!("Frame time: {:?}", last_time.elapsed());
        last_time = std::time::Instant::now();;
    }
}
