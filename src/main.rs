use nannou::prelude::*;
use glam::Vec2;
use rand::prelude::*;

struct Model {
    _window: window::Id,
    scene: Scene,
}


#[derive(Default)]
#[derive(Clone)]
#[derive(PartialEq)]
enum CellType {
    Solid,
    #[default]
    Air,
    Water,
}

#[derive(Default)]
#[derive(Clone)]
struct Cell {
    kind: CellType,
    centre: Vec2,
}

#[derive(Clone)]
struct Particle {
    position: Vec2,
    velocity: Vec2,
}

#[derive(Default)]
struct _Obstacle {
    _position: Vec2,
    _velocity: Vec2,
    _radius: f32,
}

struct Scene {
    // parameters
    gravity: Vec2,
    dt: f32,
    _flip_pic_ratio: f32,
    _over_relaxation: f32,
    initial_water_percent: f32,
    initial_particles_per_cell: u32,
    _obstacle: _Obstacle,
    cell_length: f32,

    //calculated
    grid_width: u32,
    grid_height: u32,
    particle_radius: f32,

    // grid
    cells: Vec<Vec<Cell>>,
    velocities_x: Vec<Vec<f32>>,
    velocities_y: Vec<Vec<f32>>,
    velocities_x_prev: Vec<Vec<f32>>,
    velocities_y_prev: Vec<Vec<f32>>,

    // particles
    particles: Vec<Particle>,
}

impl Default for Scene {
    fn default() -> Self {
        Scene {
            // set
            gravity: Vec2::new(0.0, -9.81),
            dt: 1.0 / 10.0,
            _flip_pic_ratio: 0.8,
            _over_relaxation: 1.9,
            initial_water_percent: 0.4,
            initial_particles_per_cell: 3,
            cell_length: 15.0,
            particle_radius: 3.0,

            // calculated
            grid_width: 0,
            grid_height: 0,
            cells: vec![vec![]],
            particles: vec![],
            _obstacle: _Obstacle::default(),
            velocities_x: vec![vec![]],
            velocities_y: vec![vec![]],
            velocities_x_prev: vec![vec![]],
            velocities_y_prev: vec![vec![]],
        }
    }
}


fn main() {
    nannou::app(model)
        .update(update)
        .run();
}

// setup
fn model(app: &App) -> Model {
    // initialise model
    let mut model = Model { 
        _window: app.new_window().view(view).build().unwrap(),
        scene: Scene::default(),
    };
    model.scene.initialise_scene(app);
    model
}


// frame by frame update
fn update(_app: &App, model: &mut Model, _update: Update) {
    model.scene.integrate_particles();
    model.scene.handle_collisions();
    // transfer velocities particles ->  grid // update cell contains particle
    // make incompressible
    model.scene.g_p_transfer_velocities();// update cell contains particle
    // compute density
    // color cell by particle density
    // give cells a color value
    // white least dense -> dark blue most dense
}

// render
fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    draw.background().color(BLACK);
    let inner_dimensions = app.main_window().inner_size_points();

    for particle in &model.scene.particles {
        draw.ellipse()
            .color(AQUAMARINE)
            .radius(model.scene.particle_radius)
            .x_y(particle.position.x - 0.5 * inner_dimensions.0, particle.position.y - 0.5 * inner_dimensions.1);
    }

    for row in &model.scene.cells {
        for cell in row {
            match cell.kind {
                CellType::Solid => {
                    draw.rect()
                        .color(SLATEGRAY)
                        .w_h(model.scene.cell_length, model.scene.cell_length)
                        .x_y(cell.centre.x - 0.5 * inner_dimensions.0, cell.centre.y - 0.5 * inner_dimensions.1);
                },
                _ => {},
            }
        }
    }

    draw.to_frame(app, &frame).unwrap();
}

impl Scene {
    fn initialise_scene(&mut self, app: &App) {
        let inner_dimensions = app.main_window().inner_size_points();

        // create initial dimensions as even number
        self.grid_width = (inner_dimensions.0 / self.cell_length).floor() as u32; 
        if self.grid_width % 2 != 0 { self.grid_width -= 1; }

        self.grid_height = (inner_dimensions.1 / self.cell_length).floor() as u32; 
        if self.grid_height % 2 != 0 { self.grid_height -= 1; }

        self.cells = vec![vec![Cell::default(); self.grid_width as usize]; self.grid_height as usize];
       
        // initialise velocity fields 
        self.velocities_x = vec![vec![0.0; self.grid_width as usize + 1]; self.grid_height as usize + 1];
        self.velocities_y = vec![vec![0.0; self.grid_width as usize + 1]; self.grid_height as usize + 1];
        self.velocities_x_prev = vec![vec![0.0; self.grid_width as usize + 1]; self.grid_height as usize + 1];
        self.velocities_y_prev = vec![vec![0.0; self.grid_width as usize + 1]; self.grid_height as usize + 1];

        // calculate cell centres for every cell
        let mut i = 0;
        for row in &mut self.cells {
            let mut j = 0;
            for cell in row {
                cell.centre = Vec2::new(
                    (j as f32 * self.cell_length) + (self.cell_length / 2.0),
                    (i as f32 * self.cell_length) + (self.cell_length / 2.0),
                );
                j += 1;
            }
            i += 1;
        }


        // initialise particles
        let initial_water_rows = (self.initial_water_percent * (self.grid_height as f32 - 3.0)).floor() as u32;
        for i in 0..initial_water_rows {
            let mut j = 0;
            for cell in &mut self.cells[(self.grid_height - 3 - i) as usize] {
                if j == 0 || j == self.grid_width - 1 { j += 1; continue; }
                if j == 1 || j == self.grid_width - 2 { j += 1; continue; }

                let half_len = self.cell_length * 0.5;

                let mut rng = rand::thread_rng();
                for _ in 0..self.initial_particles_per_cell {
                    let particle = Particle {
                        position: Vec2::new(
                            rng.gen_range(cell.centre.x - half_len..cell.centre.x + half_len),
                            rng.gen_range(cell.centre.y - half_len..cell.centre.y + half_len),
                        ),
                        velocity: Vec2::new(
                            0.0,
                            0.0,
                        ),
                    };
                    self.particles.push(particle);
                }
                j += 1;
            }
        }

        // initialise border walls 
        let mut i = 0;
        for row  in &mut self.cells {
            let mut j = 0;
            for cell in row {
                if i == 0 || i == self.grid_height - 1 { cell.kind = CellType::Solid; }
                if j == 0 || j == self.grid_width - 1 { cell.kind = CellType::Solid; }
                j += 1;
            }
            i += 1;
        }
    }

    fn integrate_particles(&mut self) {
        for particle in &mut self.particles {
            particle.velocity += self.dt * self.gravity;
            particle.position += self.dt * particle.velocity;
        }
    }

    fn handle_collisions(&mut self) {
        for particle in &mut self.particles {
            // handle collision with obstacle
            // TODO

            // handle collision with walls
            let mut x = particle.position.x;
            let mut y = particle.position.y;
            let min_x = 0.0 + self.cell_length + 0.5 * self.particle_radius;
            let min_y = 0.0 + self.cell_length + 0.5 * self.particle_radius;
            let max_x = self.grid_width as f32 * self.cell_length - 1.0 * self.cell_length;
            let max_y = self.grid_height as f32 * self.cell_length - 1.0 * self.cell_length;

            if x < min_x {
                x = min_x;
                particle.velocity.x *= 0.0;
            }
            if x > max_x {
                x = max_x;
                particle.velocity.x *= 0.0;
            }
            if y < min_y {
                y = min_y;
                particle.velocity.y *= 0.0;
            }
            if y > max_y {
                y = max_y;
                particle.velocity.y *= 0.0;
            }
            particle.position.x = x;
            particle.position.y = y;
        }

    }

    // grid -> particle
    fn g_p_transfer_velocities(&mut self) {
        for particle in &mut self.particles {
            // transfer y component
            let h = self.cell_length;
            let local_p = Vec2::new(particle.position.x, particle.position.y - (h / 2.0));
            let x_c = (local_p.x / h).floor() as usize;
            let y_c = (local_p.y / h).floor() as usize;

            let dx = local_p.x - x_c as f32 * h;
            let dy = local_p.y - y_c as f32 * h;

            let mut w1 = (1.0 - dx / h) * (1.0 - dy / h);
            let mut w2 = (dx / h)  * (1.0 - dy / h);
            let mut w3 = (dx / h) * (dy / h);
            let mut w4 = (1.0 -  dx / h) * (dy / h);

            let mut v1 = self.velocities_y[y_c][x_c];
            let mut v2 = self.velocities_y[y_c][x_c + 1];
            let mut v3 = self.velocities_y[y_c + 1][x_c + 1];
            let mut v4 = self.velocities_y[y_c + 1][x_c];

            // check if non water cells, then dont include in subsequent calcs
            if self.cells[y_c - 1][x_c - 1].kind != CellType::Water
            || self.cells[y_c - 1][x_c].kind != CellType::Water 
            || self.cells[y_c][x_c - 1].kind != CellType::Water {
                v1 = 0.0; w1 = 0.0;
            }
            if self.cells[y_c - 1][x_c].kind != CellType::Water
            || self.cells[y_c - 1][x_c + 1].kind != CellType::Water 
            || self.cells[y_c][x_c + 1].kind != CellType::Water {
                v2 = 0.0; w2 = 0.0;
            }
            if self.cells[y_c][x_c + 1].kind != CellType::Water
            || self.cells[y_c + 1][x_c + 1].kind != CellType::Water 
            || self.cells[y_c + 1][x_c].kind != CellType::Water {
                v3 = 0.0; w3 = 0.0;
            }
            if self.cells[y_c + 1][x_c].kind != CellType::Water
            || self.cells[y_c + 1][x_c - 1].kind != CellType::Water 
            || self.cells[y_c][x_c - 1].kind != CellType::Water {
                v4 = 0.0; w4 = 0.0;
            }
            let vp = (v1 * w1 + v2 * w2 + v3 * w3 + v4 * w4) / (w1 + w2 + w3 + w4);

            particle.velocity.y = vp;

            //todo set cell types
        }
    }
}
