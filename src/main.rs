use nannou::prelude::*;
use glam::{Vec2, UVec2};
use rand::prelude::*;

struct Model {
    _window: window::Id,
    scene: Scene,
}

#[derive(Default)]
#[derive(Clone, Copy)]
#[derive(PartialEq)]
enum CellType {
    Solid,
    #[default]
    Air,
    Water,
}



#[derive(Default)]
struct Grid {
    types: Vec<CellType>,
    velocities: Vec<Velocities>,
    size: UVec2,
    vel_size: UVec2,
}

#[derive(Default)]
#[derive(Clone, Copy)]
struct Velocities {
    velocity: Vec2,
    weight: Vec2,
    _prev_velocity: Vec2,
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
    particle_radius: f32,

    // grid
    grid: Grid,

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
            cell_length: 25.0,

            // calculated
            particle_radius: 0.0,
            grid: Grid::default(),
            particles: vec![],
            _obstacle: _Obstacle::default(),
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
    model.scene.p_g_transfer_velocities();
    model.scene.update_celltype();
    model.scene.g_p_transfer_velocities();
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
    let h = model.scene.cell_length;

    for y in 0..model.scene.grid.size.y {
        for x in 0..model.scene.grid.size.x {
            match model.scene.grid.get_type(x,y) {
                CellType::Solid => {
                    draw.rect()
                        .color(SLATEGRAY)
                        .w_h(h, h)
                        .x_y(((x as f32 + 0.5) * h) - 0.5 * inner_dimensions.0, ((y as f32 + 0.5) * h) - 0.5 * inner_dimensions.1);
                },
                CellType::Water => {
                    draw.rect()
                        .color(BLUE)
                        .w_h(h, h)
                        .x_y(((x as f32 + 0.5) * h) - 0.5 * inner_dimensions.0, ((y as f32 + 0.5) * h) - 0.5 * inner_dimensions.1);
                },
                _ => {},
            }
        }
    }

    for particle in &model.scene.particles {
        draw.ellipse()
            .color(AQUAMARINE)
            .radius(model.scene.particle_radius)
            .x_y(particle.position.x - 0.5 * inner_dimensions.0, particle.position.y - 0.5 * inner_dimensions.1);
    }

    draw.to_frame(app, &frame).unwrap();
}

impl Scene {
    fn initialise_scene(&mut self, app: &App) {
        let inner_dimensions = app.main_window().inner_size_points();

        self.particle_radius = self.cell_length / 8.0;

        // create initial axiss of grids, ensuring its even
        self.grid.size.x = (inner_dimensions.0 / self.cell_length).floor() as u32; 
        if self.grid.size.x % 2 != 0 { self.grid.size.x -= 1; }

        self.grid.size.y = (inner_dimensions.1 / self.cell_length).floor() as u32; 
        if self.grid.size.y % 2 != 0 { self.grid.size.y -= 1; }

        self.grid.vel_size.x = self.grid.size.x + 1;
        self.grid.vel_size.y = self.grid.size.y + 1;


        // initialise grid(s)
        self.grid.velocities = vec![
            Velocities {
                velocity: Vec2::ZERO,
                weight: Vec2::ZERO,
                _prev_velocity: Vec2::ZERO,
            }; 
            ((self.grid.size.x + 1) * (self.grid.size.y + 1)) as usize
        ];
        self.grid.types = vec![CellType::default(); (self.grid.size.x * self.grid.size.y) as usize];

        // initialise particles
        let initial_water_rows = (self.initial_water_percent * (self.grid.size.y as f32 - 3.0)).floor() as u32;
        for n in 3..initial_water_rows + 3 {

            for x in 3..(self.grid.size.x - 3) {
                let h = self.cell_length;
                let mut rng = rand::thread_rng();

                for _ in 0..self.initial_particles_per_cell {
                    let particle = Particle {
                        position: Vec2::new(
                            rng.gen_range((x as f32 ) * h..(x as f32 + 1.0) * h),
                            rng.gen_range((self.grid.size.y as f32 - n as f32) * h..(self.grid.size.y as f32 - n as f32 + 1.0) * h),
                        ),
                        velocity: Vec2::new(
                            0.0,
                            0.0,
                        ),
                    };
                    self.particles.push(particle);
                }
            }
        }

        // initialise border walls 
        for y in 0..self.grid.size.y {
            for x in 0..self.grid.size.x {
                if x == 0 || x == self.grid.size.x - 1 { *self.grid.get_type_mut(x,y) = CellType::Solid; }
                if y == 0 || y == self.grid.size.y - 1 { *self.grid.get_type_mut(x,y) = CellType::Solid; }
            }
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
            let max_x = self.grid.size.x as f32 * self.cell_length - 1.0 * self.cell_length;
            let max_y = self.grid.size.y as f32 * self.cell_length - 1.0 * self.cell_length;

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
            for axis in 0..=1 {

                let h = self.cell_length;

                let local_p = match axis { 
                    0 => Vec2::new(particle.position.x + (h / 2.0), particle.position.y),
                    _ => Vec2::new(particle.position.x, particle.position.y - (h / 2.0)),
                };

                let x_c = (local_p.x / h).floor();
                let y_c = (local_p.y / h).floor();
                let c = UVec2::new(x_c as u32, y_c as u32);

                let dx = local_p.x - c.x as f32 * h;
                let dy = local_p.y - c.y as f32 * h;

                let mut w1 = (1.0 - dx / h) * (1.0 - dy / h);
                let mut w2 = (dx / h)  * (1.0 - dy / h);
                let mut w3 = (dx / h) * (dy / h);
                let mut w4 = (1.0 -  dx / h) * (dy / h);

                let mut v1 = match axis {
                    0 => self.grid.get_val(c.x, c.y).velocity.x,
                    _ => self.grid.get_val(c.x, c.y).velocity.y,
                };
                let mut v2 = match axis {
                    0 => self.grid.get_val(c.x + 1, c.y).velocity.x,
                    _ => self.grid.get_val(c.x + 1, c.y).velocity.y,
                };
                let mut v3 = match axis {
                    0 => self.grid.get_val(c.x + 1, c.y + 1).velocity.x,
                    _ => self.grid.get_val(c.x + 1, c.y + 1).velocity.y,
                };
                let mut v4 = match axis {
                    0 => self.grid.get_val(c.x, c.y + 1).velocity.x,
                    _ => self.grid.get_val(c.x, c.y + 1).velocity.y,
                };

                // check if non water cells, then dont include in subsequent calc
                match c {
                    // handle edge cases
                    c if c.x == 0 => {
                        v1 = 0.0; v4 = 0.0;
                        w1 = 0.0; w4 = 0.0;
                    },
                    c if c.x == self.grid.size.x => {
                        v2 = 0.0; v3 = 0.0;
                        w2 = 0.0; w3 = 0.0;
                    },
                    c if c.y == 0 => {
                        v1 = 0.0; v2 = 0.0;
                        w1 = 0.0; w2 = 0.0;
                    },
                    c if c.y == self.grid.size.y => {
                        v3 = 0.0; v4 = 0.0;
                        w3 = 0.0; w4 = 0.0;
                    },
                    c if c.x == 0 && c.y == 0 => {
                        v1 = 0.0; v2 = 0.0; v4 = 0.0;
                        w1 = 0.0; w2 = 0.0; w4 = 0.0;
                    },
                    c if c.x == self.grid.size.x && c.y == self.grid.size.y => {
                        v2 = 0.0; v3 = 0.0; v4 = 0.0;
                        w2 = 0.0; w3 = 0.0; w4 = 0.0;
                    },
                    // handle general case
                    _ => {
                        // velocities stored counter clockwise from bottom left, starting at v1
                        // check for v1
                        if self.grid.get_type(c.x - 1, c.y - 1) != CellType::Water
                        || self.grid.get_type(c.x, c.y - 1) != CellType::Water
                        || self.grid.get_type(c.x - 1, c.y) != CellType::Water {
                            v1 = 0.0; w1 = 0.0;
                        }
                        // check for v2
                        if self.grid.get_type(c.x, c.y - 1) != CellType::Water
                        || self.grid.get_type(c.x + 1, c.y - 1) != CellType::Water
                        || self.grid.get_type(c.x + 1, c.y) != CellType::Water {
                            v2 = 0.0; w2 = 0.0;
                        }
                        // check for v3
                        if self.grid.get_type(c.x + 1, c.y) != CellType::Water
                        || self.grid.get_type(c.x + 1, c.y + 1) != CellType::Water
                        || self.grid.get_type(c.x, c.y + 1) != CellType::Water {
                            v3 = 0.0; w3 = 0.0;
                        }
                        // check for v4
                        if self.grid.get_type(c.x, c.y + 1) != CellType::Water
                        || self.grid.get_type(c.x - 1, c.y + 1) != CellType::Water
                        || self.grid.get_type(c.x - 1, c.y) != CellType::Water {
                            v4 = 0.0; w4 = 0.0;
                        }
                    },
                }
                let mut vp = (v1 * w1 + v2 * w2 + v3 * w3 + v4 * w4) / (w1 + w2 + w3 + w4);
                if w1 + w2 + w3 + w4 == 0.0 { vp = 0.0; }

                match axis {
                    0 => particle.velocity.x = vp,
                    _ => particle.velocity.y = vp,
                };
            }
        }
    }

    fn p_g_transfer_velocities(&mut self) {
        for y in 0..self.grid.size.y {
            for x in 0..self.grid.size.x {
                self.grid.get_val_mut(x, y).velocity = Vec2::ZERO;
                self.grid.get_val_mut(x, y).weight = Vec2::ZERO;
            }
        }
        
        for particle in &mut self.particles {
            for axis in 0..=1 {
                let h = self.cell_length;

                let local_p = match axis { 
                    0 => Vec2::new(particle.position.x + (h / 2.0), particle.position.y),
                    _ => Vec2::new(particle.position.x, particle.position.y - (h / 2.0)),
                };

                let x_c = (local_p.x / h).floor();
                let y_c = (local_p.y / h).floor();
                let c = UVec2::new(x_c as u32, y_c as u32);

                let dx = local_p.x - c.x as f32 * h;
                let dy = local_p.y - c.x as f32 * h;

                let w1 = (1.0 - dx / h) * (1.0 - dy / h);
                let w2 = (dx / h)  * (1.0 - dy / h);
                let w3 = (dx / h) * (dy / h);
                let w4 = (1.0 -  dx / h) * (dy / h);

                match axis {
                    0 => {
                        self.grid.get_val_mut(c.x, c.y).velocity.x += w1 * particle.velocity.x;
                        self.grid.get_val_mut(c.x + 1, c.y).velocity.x += w2 * particle.velocity.x;
                        self.grid.get_val_mut(c.x + 1, c.y + 1).velocity.x += w3 * particle.velocity.x;
                        self.grid.get_val_mut(c.x, c.y + 1).velocity.x += w4 * particle.velocity.x;

                        self.grid.get_val_mut(c.x, c.y).weight.x += w1;
                        self.grid.get_val_mut(c.x + 1, c.y).weight.x += w2;
                        self.grid.get_val_mut(c.x + 1, c.y + 1).weight.x += w3;
                        self.grid.get_val_mut(c.x, c.y + 1).weight.x += w4;
                    },
                    _ => {
                        self.grid.get_val_mut(c.x, c.y).velocity.y += w1 * particle.velocity.y;
                        self.grid.get_val_mut(c.x + 1, c.y).velocity.y += w2 * particle.velocity.y;
                        self.grid.get_val_mut(c.x + 1, c.y + 1).velocity.y += w3 * particle.velocity.y;
                        self.grid.get_val_mut(c.x, c.y + 1).velocity.y += w4 * particle.velocity.y;

                        self.grid.get_val_mut(c.x, c.y).weight.y += w1;
                        self.grid.get_val_mut(c.x + 1, c.y).weight.y += w2;
                        self.grid.get_val_mut(c.x + 1, c.y + 1).weight.y += w3;
                        self.grid.get_val_mut(c.x, c.y + 1).weight.y += w4;
                    },
                }
            }
        }

        for y in 0..self.grid.size.y {
            for x in 0..self.grid.size.x {
                self.grid.get_val_mut(x, y).velocity = self.grid.get_val(x, y).velocity / self.grid.get_val(x, y).weight;
            }
        }
    }

    fn solve_incompressiblity(&mut self) {
        let iterations = 1;
        for _ in 0..iterations {
            for y in 0..self.grid.size.y {
                for x in 0..self.grid.size.x {
                    for axis in 0..=1 {

                    }
                }
            }
        }
    }

    fn update_celltype(&mut self) {
        for y in 0..self.grid.size.y {
            for x in 0..self.grid.size.x {
                if self.grid.get_type(x, y) != CellType::Solid {
                    *self.grid.get_type_mut(x, y) = CellType::Air;
                }
            }
        }
        for particle in &self.particles {
            let x_c = (particle.position.x / self.cell_length).floor();
            let y_c = (particle.position.y / self.cell_length).floor();
            *self.grid.get_type_mut(x_c as u32, y_c as u32) = CellType::Water;
        }
    } 
}


impl Grid {
    fn get_type(&self, x: u32, y: u32) -> CellType {
        self.types[(x + y * self.size.x) as usize]
    }
    fn get_type_mut(&mut self, x: u32, y: u32) -> &mut CellType {
        &mut self.types[(x + y * self.size.x) as usize]
    }

    fn get_val(&self, x: u32, y: u32) -> Velocities {
        self.velocities[(x + y * self.vel_size.x) as usize]
    }

    fn get_val_mut(&mut self, x: u32, y: u32) -> &mut Velocities {
        &mut self.velocities[(x + y * self.vel_size.x) as usize]
    }
}
