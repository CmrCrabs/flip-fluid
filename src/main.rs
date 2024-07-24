use nannou::prelude::*;
use glam::{Vec2, UVec2};
use rand::prelude::*;

struct Model {
    _window: window::Id,
    scene: Scene,
}

#[derive(Default)]
#[derive(Debug)]
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
    prev_velocity: Vec2,
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
    flip_pic_ratio: f32,
    over_relaxation: f32,
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
            dt: 1.0 / 3.0,
            flip_pic_ratio: 0.8,
            over_relaxation: 1.9,
            initial_water_percent: 0.4,
            initial_particles_per_cell: 1,
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
    //model.scene.solve_incompressiblity();
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
                //CellType::Water => {
                //    draw.rect()
                //        .color(BLUE)
                //        .w_h(h, h)
                //        .x_y(((x as f32 + 0.5) * h) - 0.5 * inner_dimensions.0, ((y as f32 + 0.5) * h) - 0.5 * inner_dimensions.1);
                //},
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
        let mut inner_dimensions = app.main_window().inner_size_points();
        inner_dimensions.0 *= 1.5; inner_dimensions.1 *= 1.3;

        self.particle_radius = self.cell_length / 8.0;

        // create initial axiss of grids, ensuring its even
        self.grid.size.x = (inner_dimensions.0 / self.cell_length).floor() as u32; 
        self.grid.size.y = (inner_dimensions.1 / self.cell_length).floor() as u32; 

        self.grid.vel_size.x = self.grid.size.x + 1;
        self.grid.vel_size.y = self.grid.size.y + 1;


        // initialise grid(s)
        self.grid.velocities = vec![
            Velocities {
                velocity: Vec2::ZERO,
                weight: Vec2::ZERO,
                prev_velocity: Vec2::ZERO,
            }; 
            ((self.grid.size.x + 1) * (self.grid.size.y + 1)) as usize
        ];
        self.grid.types = vec![CellType::default(); (self.grid.size.x * self.grid.size.y) as usize];

        // initialise particles
        let initial_water_rows = (self.initial_water_percent * (self.grid.size.y as f32 - 3.0)).floor() as u32;
        for n in 5..initial_water_rows + 5 {

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
            // TODO handle collision with obstacle

            // handle collision with walls
            let mut x = particle.position.x;
            let mut y = particle.position.y;
            let min_x = self.cell_length + 0.5 * self.particle_radius;
            let min_y = self.cell_length + 0.5 * self.particle_radius;
            let max_x = self.grid.size.x as f32 * self.cell_length - 1.0 * self.cell_length - 0.5 * self.particle_radius;
            let max_y = self.grid.size.y as f32 * self.cell_length - 1.0 * self.cell_length - 0.5 * self.particle_radius;

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

    fn p_g_transfer_velocities(&mut self) {
        for y in 0..self.grid.vel_size.y {
            for x in 0..self.grid.vel_size.x {
                self.grid.get_val_mut(x, y).prev_velocity = self.grid.get_val(x, y).velocity;
                self.grid.get_val_mut(x, y).velocity = Vec2::ZERO;
                self.grid.get_val_mut(x, y).weight = Vec2::ZERO;
            }
        }
        
        for particle in &mut self.particles {
            for axis in 0..=1 {
                let h = self.cell_length;

                let local_p = match axis { 
                    0 => Vec2::new(particle.position.x, particle.position.y - (h / 2.0)),
                    _ => Vec2::new(particle.position.x - (h / 2.0), particle.position.y),
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

        for y in 0..self.grid.vel_size.y {
            for x in 0..self.grid.vel_size.x {
                self.grid.get_val_mut(x, y).velocity = self.grid.get_val(x, y).velocity / self.grid.get_val(x, y).weight;
            }
        }
    }

    // grid -> particle
    fn g_p_transfer_velocities(&mut self) {
        for particle in &mut self.particles {
            for axis in 0..=1 {

                let h = self.cell_length;

                let local_p = match axis { 
                    0 => Vec2::new(particle.position.x, particle.position.y - (h / 2.0)),
                    _ => Vec2::new(particle.position.x - (h / 2.0), particle.position.y),
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

                // check if air cell, then dont include in subsequent calc
                match c {
                    // skip edge cases
                    c if c.x == 0 => { },
                    c if c.x == self.grid.size.x - 1 => { },
                    c if c.y == 0 => { },
                    c if c.y == self.grid.size.y - 1 => { },
                    c if c.x == 0 && c.y == 0 => { },
                    c if c.x == self.grid.size.x - 1 && c.y == self.grid.size.y - 1 => { },
                    // general case
                    // velocities stored counter clockwise from bottom left, starting at v1
                    _ => {
                        if self.grid.get_type(c.x, c.y + 1) == CellType::Air {
                            v4 = 0.0; v3 = 0.0;
                            w4 = 0.0; w3 = 0.0;
                        }
                        if self.grid.get_type(c.x + 1, c.y) == CellType::Air {
                            v2 = 0.0; v3 = 0.0;
                            w2 = 0.0; w3 = 0.0;
                        }
                        if self.grid.get_type(c.x - 1, c.y) == CellType::Air {
                            v4 = 0.0; v1 = 0.0;
                            w4 = 0.0; w1 = 0.0;
                        }
                        if self.grid.get_type(c.x, c.y - 1) == CellType::Air {
                            v1 = 0.0; v2 = 0.0;
                            w1 = 0.0; w2 = 0.0;
                        }
                    },
                }

                let mut vp = (v1 * w1 + v2 * w2 + v3 * w3 + v4 * w4) / (w1 + w2 + w3 + w4);
                if w1 + w2 + w3 + w4 == 0.0 { vp = 0.0; }

                let picv = vp;
                let flipv = match axis {
                    0 => particle.velocity.x + (vp - self.grid.get_val(c.x, c.y).prev_velocity.x),
                    _ => particle.velocity.y + (vp - self.grid.get_val(c.x, c.y).prev_velocity.y),
                };
                let vel = (1.0 - self.flip_pic_ratio) * picv + self.flip_pic_ratio * flipv;

                //println!("picv: {:?}", picv);
                //println!("flipv: {:?}", flipv);
                //println!("vel: {:?}", vel);
                //println!("w1: {:?} v1: {:?}", w1, v1);
                //println!("w2: {:?} v2: {:?}", w2, v2);
                //println!("w3: {:?} v3: {:?}", w3, v3);
                //println!("w4: {:?} v4: {:?}", w4, v4);

                match axis {
                    0 => particle.velocity.x = vp,
                    _ => particle.velocity.y = vp,
                };
            }
        }
    }


    fn solve_incompressiblity(&mut self) {
        let iterations = 1;
        for _ in 0..iterations {
            for y in 0..self.grid.size.y {
                for x in 0..self.grid.size.x {
                    if self.grid.get_type(x, y) != CellType::Water { continue; }

                    let mut divergence = 
                    self.grid.get_val(x + 1, y).velocity.x 
                    - self.grid.get_val(x, y).velocity.x
                    + self.grid.get_val(x, y + 1).velocity.y 
                    - self.grid.get_val(x, y).velocity.y;

                    divergence *= self.over_relaxation;

                    let mut s = 0.0;
                    if self.grid.get_type(x - 1, y) != CellType::Solid { s += 1.0; }
                    if self.grid.get_type(x + 1, y) != CellType::Solid { s += 1.0; }
                    if self.grid.get_type(x, y - 1) != CellType::Solid { s += 1.0; }
                    if self.grid.get_type(x, y + 1) != CellType::Solid { s += 1.0; }
                    if s == 0.0 { continue; }

                    self.grid.get_val_mut(x, y).velocity.x += divergence / s;
                    self.grid.get_val_mut(x + 1, y).velocity.x -= divergence / s;
                    self.grid.get_val_mut(x, y).velocity.y += divergence / s;
                    self.grid.get_val_mut(x, y + 1).velocity.y -= divergence / s;
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
