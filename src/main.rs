use nannou::prelude::*;
use glam::Vec2;
use rand::prelude::*;

struct Model {
    _window: window::Id,
    scene: Scene,
}

type Grid = Vec<Vec<Cell>>;

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
    contains_particle: bool,
    kind: CellType,
    centre: Vec2,
    velocity: Vec2,
    q: Vec4,
    r: Vec4,
    }

#[derive(Clone)]
struct Particle {
    position: Vec2,
    velocity: Vec2,
}

struct Scene {
    // set
    gravity: Vec2,
    dt: f32,
    _flip_pic_ratio: f32,
    _over_relaxation: f32,
    initial_water_percent: f32,
    initial_particles_per_cell: u32,
    cell_length: f32,

    //calculated
    grid_width: u32,
    grid_height: u32,
    cells: Grid,
    particles: Vec<Particle>,
    particle_radius: f32,
    top_left: Vec2,
}

impl Default for Scene {
    fn default() -> Self {
        Scene {
            // set
            gravity: Vec2::new(0.0, -9.81),
            dt: 1.0 / 30.0,
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
            top_left: Vec2::default(),
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
    // simulate particles
    model.scene.integrate_particles();

    // velocity transfer particles -> grid
    // update cell contains particle

    // incompressible

    // velocity transfer grid -> particles
    // update cell contains particle

    // compute density
    
    // color cell by particle density
    // give cells a color value
    // white least dense -> dark blue most dense
    
    // update cell kind 
    model.scene.update_celltype();
}

// render
fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    draw.background().color(BLACK);

    let inner_dimensions = app.main_window().outer_size_points();

    for particle in &model.scene.particles {
        draw.ellipse()
            .color(AQUAMARINE)
            .radius(model.scene.particle_radius)
            .x_y(particle.position.x - (0.5 * inner_dimensions.0), particle.position.y - (0.5 * inner_dimensions.1));
    }

    for row in &model.scene.cells {
        for cell in row {
            match cell.kind {
                CellType::Solid => {
                    draw.rect()
                        .color(SLATEGRAY)
                        .w_h(model.scene.cell_length, model.scene.cell_length)
                        .x_y(cell.centre.x - (0.5 * inner_dimensions.0), cell.centre.y - (0.5 * inner_dimensions.1));
                },
                _ => {},
            }
        }
    }
    draw.to_frame(app, &frame).unwrap();
}

impl Scene {
    fn initialise_scene(&mut self, app: &App) {
        //get inner dimensions
        let mut inner_dimensions = app.main_window().inner_size_points();
        inner_dimensions.0 *= 1.9;
        inner_dimensions.1 *= 1.4;


        // create initial scene values
        self.grid_width = (inner_dimensions.0 / self.cell_length).floor() as u32; 
        if self.grid_width % 2 != 0 { self.grid_width -= 1; }

        self.grid_height = (inner_dimensions.1 / self.cell_length).floor() as u32; 
        if self.grid_height % 2 != 0 { self.grid_height -= 1; }

        self.cells = vec![vec![Cell::default(); self.grid_width as usize]; self.grid_height as usize];
       
        // calculate cell centres for every cell
        // really stupid and bad and will be rewritten fully
        self.top_left = Vec2::new(0.0, inner_dimensions.1);
        let mut i = 0;
        for row in &mut self.cells {
            let mut j = 0;
            for cell in row {
                // calculate cell centres for every cell
                cell.centre = Vec2::new(
                    self.top_left.x + (j as f32 * self.cell_length) + (self.cell_length / 2.0),
                    self.top_left.y - (i as f32 * self.cell_length) - (self.cell_length / 2.0),
                );
                j += 1;
            }
            i += 1;
        }

        // initialise particles
        let initial_water_rows = (self.initial_water_percent * (self.grid_height as f32 - 2.0)).floor() as u32;
        for i in 1..initial_water_rows {
            let mut j = 0;
            for cell in &mut self.cells[i as usize] {
                if j == 0 || j == self.grid_width - 1 { j += 1; continue; }
                cell.contains_particle = true;
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
                if i == 0 || i == self.grid_height -1 { cell.kind = CellType::Solid; }
                else if j == 0 || j == self.grid_width - 1 { cell.kind = CellType::Solid; }
                j += 1;
            }
            i += 1;
        }
    }

    fn integrate_particles(&mut self) {
        for particle in &mut self.particles {
            particle.velocity += self.dt * self.gravity;

            // bad obstacle detection
            //let cell_pos = ((particle.position + self.dt * particle.velocity) / self.cell_length).floor();
            //if self.cells[cell_pos.y as usize][cell_pos.x as usize].kind == CellType::Solid {
            //    particle.velocity *= 0.0;
            //}
            particle.position += self.dt * particle.velocity;
        }
    }

    // particle -> grid
    fn p_g_transfer_velocities(&mut self) {
        // set q and r to zero for all cells
        for row in &mut self.cells {
            for cell in row {
                cell.q *= 0.0;
                cell.r *= 0.0;
            }
        }

        for particle in &mut self.particles {
            let h = self.cell_length;
            let local_p = Vec2::new(particle.position.x, particle.position.y - (h / 2.0));
            let x_c = (local_p.x / h).floor();
            let y_c = (local_p.y / h).floor();
            let dx = local_p.x - x_c * h;
            let dy = local_p.y - y_c * h;

            let w = Vec4::new(
                (1.0 - dx / h) * (1.0 - dy / h),
                (dx / h)  * (1.0 - dy / h),
                (dx / h) * (dy / h),
                (1.0 -  dx / h) * (dy / h),
            );


            // deal with q maybe being undefined
            // for every defined corner 
            // sum w * q
            // divide by sum of weights
            let qp = self.cells[y_c as usize][x_c as usize].q.dot(w) / (w.x + w.y + w.z + w.w);

            self.cells[x_c as usize][y_c as usize].r += w;
            self.cells[x_c as usize][y_c as usize].q += w * qp;
        }
        for row in &mut self.cells {
            for cell in row {
                cell.q *= 1.0 / cell.r;
            }
        }
    }

    fn update_celltype(&mut self) {
        for row in &mut self.cells {
            for cell in row {
                match cell.kind {
                    CellType::Solid => continue,
                    _ => cell.kind = match cell.contains_particle {
                        true => CellType::Water,
                        false => CellType::Air,
                    },
                }
            }
        }
    }
}
