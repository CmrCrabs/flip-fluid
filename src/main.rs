use nannou::{lyon::geom::euclid::num::Floor, prelude::*};
use glam::Vec2;

struct Model {
    _window: window::Id,
    scene: Scene,
}

type Grid = Vec<Vec<Cell>>;

#[derive(Default)]
#[derive(Clone)]
enum CellType {
    _Solid,
    #[default]
    Air,
    Water,
}

#[derive(Default)]
#[derive(Clone)]
struct Cell {
    _particles: Vec<Particle>,
    kind: CellType,
    centre: Vec2,
}

#[derive(Clone)]
struct Particle {
    _position: Vec2,
    _velocity: Vec2,
}

struct Scene {
    _gravity: Vec2,
    _dt: f32,
    _flip_pic_ratio: f32,
    _over_relaxation: f32,
    grid_width: u32,
    grid_height: u32,
    initial_water_amount: u32,
    cell_length: f32,
    padding: f32,
    cells: Grid,
    particle_radius: f32,
}

impl Default for Scene {
    fn default() -> Self {
        Scene {
            // set
            _gravity: Vec2::new(0.0, -9.81),
            _dt: 1.0 / 120.0,
            _flip_pic_ratio: 0.8,
            _over_relaxation: 1.9,
            padding: 0.0,
            initial_water_amount: 530,
            cell_length: 25.0,

            // calculated
            grid_width: 0,
            grid_height: 0,
            particle_radius: 0.0,
            cells: vec![vec![]],
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
    let mut inner_dimensions = app.main_window().inner_size_points();
    inner_dimensions.0 = inner_dimensions.0 - model.scene.padding;
    inner_dimensions.1 = inner_dimensions.1 - model.scene.padding;

    // create initial scene values
    model.scene.grid_width = (inner_dimensions.0 / model.scene.cell_length).floor() as u32; 
    if model.scene.grid_width % 2 != 0 { model.scene.grid_width -= 1; }

    model.scene.grid_height = (inner_dimensions.1 / model.scene.cell_length).floor() as u32; 
    if model.scene.grid_height % 2 != 0 { model.scene.grid_height -= 1; }

    model.scene.particle_radius = model.scene.cell_length * 0.5;
    
    model.scene.cells = vec![vec![Cell::default(); model.scene.grid_width as usize]; model.scene.grid_height as usize];

    // initialise water cells
    let r = model.scene.initial_water_amount % model.scene.grid_width;
    let n = (model.scene.initial_water_amount / model.scene.grid_width).floor();
    for i in 0..n {
        for j in 0..model.scene.grid_width as usize {
            model.scene.cells[i as usize][j as usize].kind = CellType::Water;
        }
    }
    for i in 0..r {
        model.scene.cells[n as usize][i as usize].kind = CellType::Water;
    }
    
    // calculate cell centres for every cell
    let top_left = Vec2::new(-(inner_dimensions.0 / 2.0), inner_dimensions.1 / 2.0);
    println!("top left: {:?}", top_left);
    println!("inner dimensions: {:?}", inner_dimensions);
    let mut i = 0;
    for row in &mut model.scene.cells {
        let mut j = 0;
        for cell in row {
            cell.centre = Vec2::new(
                top_left.x + (j as f32 * model.scene.cell_length) + (model.scene.cell_length / 2.0),
                top_left.y - (i as f32 * model.scene.cell_length) - (model.scene.cell_length / 2.0),
            );
            //println!("centre: {:?}", cell.centre);
            j += 1;
        }
        i += 1;
    }

    println!("p1 pos: {:?}", model.scene.cells[0][0].centre);
    model
}


// frame by frame update
fn update(_app: &App, _model: &mut Model, _update: Update) {}


// render
fn view(app: &App, model: &Model, frame: Frame) {
    let draw = app.draw();
    draw.background().color(BLACK);

    for row in &model.scene.cells {
        for cell in row {
            match cell.kind {
                CellType::Water => {
                    //println!("x: {:?}\ny: {:?}", cell.centre.x, cell.centre.y);
                    draw.ellipse()
                        .color(AQUAMARINE)
                        .radius(model.scene.particle_radius)
                        .x_y(cell.centre.x, cell.centre.y);
                },
                _ => {},
            }
        }
    }
    
    draw.to_frame(app, &frame).unwrap();
}
