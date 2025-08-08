pub use wasm_bindgen_rayon::init_thread_pool;

use rayon::prelude::*;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct NeighborLinkedList {
    next: Vec<i32>, 
    first: Vec<i32>, 
    marks: Vec<i32>, 
    neighbor_indices: Vec<Vec<u32>>, 
    current_mark: i32, 
}

impl NeighborLinkedList {
    const HASH_SIZE: usize = 370111;
    const MAX_NEIGHBOR_PARTICLES: usize = 64;

    pub fn hash(gx: usize, gy: usize) -> usize {
        ((gx * 92837111) ^ (gy * 689287499)) % Self::HASH_SIZE
    }
}

#[wasm_bindgen]
impl NeighborLinkedList {
    #[wasm_bindgen(constructor)]
    pub fn new(max_particles: u32) -> Self { 
        Self { 
            next: vec![0; max_particles as usize], 
            first: vec![0; Self::HASH_SIZE as usize], 
            marks: vec![0; Self::HASH_SIZE as usize], 
            neighbor_indices: vec![vec![0; Self::MAX_NEIGHBOR_PARTICLES + 1]; max_particles as usize], 
            current_mark: 0
        }
    }

    pub fn find_neighbors(
        &mut self, 
        pos_x_ptr: *const f32, pos_y_ptr: *const f32, 
        first_neighbors_ptr: *mut u32, neighbors_ptr: *mut u32, 
        all_particles_len: usize, grid_spacing: f32,
        max_neighbors: u32
    ) {
        let pos_x = unsafe { std::slice::from_raw_parts(pos_x_ptr, all_particles_len) };
        let pos_y = unsafe { std::slice::from_raw_parts(pos_y_ptr, all_particles_len) };
        let first_neighbors = unsafe { std::slice::from_raw_parts_mut(first_neighbors_ptr, all_particles_len + 1) };
        let neighbors = unsafe { std::slice::from_raw_parts_mut(neighbors_ptr, all_particles_len * max_neighbors as usize) };

        let inv_grid_spacing = 1.0 / grid_spacing;

        self.current_mark += 1;

        pos_x.iter()
        .zip(pos_y.iter())
        .enumerate()
        .take(all_particles_len) 
        .for_each(|(i, (x, y))|{
            let gx = (x * inv_grid_spacing).floor() as usize;
            let gy = (y * inv_grid_spacing).floor() as usize;

            let h = NeighborLinkedList::hash(gx, gy) as usize;

            if self.marks[h] != self.current_mark {
                self.marks[h] = self.current_mark;
                self.first[h] = -1;
            }

            self.next[i] = self.first[h];
            self.first[h] = i as i32;
        });


        let inv_grid_spacing = 1.0 / grid_spacing;

        let h2 = grid_spacing * grid_spacing;
        pos_x.par_iter()
        .zip(pos_y.par_iter())
        .zip(self.neighbor_indices.par_iter_mut())
        .enumerate()
        .take(all_particles_len) 
        .for_each(|(i, ((x, y), neighbors_of_i))|{
            let gx = (x * inv_grid_spacing).floor() as i32;
            let gy = (y * inv_grid_spacing).floor() as i32;

            let mut neighbor_count = 0;

            for ix in gx - 1..=gx + 1 { 
                for iy in gy - 1..=gy + 1 {
                    let h = NeighborLinkedList::hash(ix as usize, iy as usize);

                    if self.marks[h] != self.current_mark {
                        continue;
                    }

                    let mut id = self.first[h];
                    while id >= 0 {
                        let dx = pos_x[id as usize] - x;
                        let dy = pos_y[id as usize] - y;
                        let r2 = dx * dx + dy * dy;
                        if r2 < h2 && id as usize != i && neighbor_count < NeighborLinkedList::MAX_NEIGHBOR_PARTICLES {
                            neighbors_of_i[neighbor_count] = id as u32;
                            neighbor_count += 1;
                        }
                        id = self.next[id as usize];
                    }
                }
            }

            // brute force
            // for id in 0..all_particles_len {
            //     let dx = pos_x[id as usize] - x;
            //     let dy = pos_y[id as usize] - y;
            //     let r2 = dx * dx + dy * dy;
            //     if r2 < h2 && id as usize != i && neighbor_count < NeighborLinkedList::MAX_NEIGHBOR_PARTICLES {
            //         neighbors_of_i[neighbor_count] = id as u32;
            //         neighbor_count += 1;
            //     }
            // }

            *neighbors_of_i.last_mut().unwrap() = neighbor_count as u32;
        });

        let mut accum: usize = 0;
        // flatten neighborhood list
        self.neighbor_indices.iter()
        .enumerate()
        .take(all_particles_len) // ここが重要。
        .for_each(|(i, neighbors_of_i)|{
            let neighbor_count = *neighbors_of_i.last().unwrap() as usize;
            first_neighbors[i] = accum as u32;
            neighbors[accum..accum + neighbor_count].copy_from_slice(&neighbors_of_i[0..neighbor_count]);
            accum += neighbor_count;
        });

        first_neighbors[all_particles_len] = accum as u32;
    }
}

#[wasm_bindgen]
pub fn alloc_buffer(n: usize) -> *mut u8 {
    let mut buf = Vec::with_capacity(n);
    buf.resize(n, 0);
    let ptr = buf.as_mut_ptr();
    std::mem::forget(buf);
    ptr
}

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

fn log_i32(num: i32) {
    log("rust:");
    log(num.to_string().as_str());
}


#[wasm_bindgen]
pub fn solve_fluid_jacobi(
    pos_x_ptr: *const f32, pos_y_ptr: *const f32, 
    lambdas_ptr: *mut f32, grads_x_ptr: *mut f32, grads_y_ptr: *mut f32, 
    fluid_deltas_x_ptr: *mut f32, fluid_deltas_y_ptr: *mut f32, 
    first_neighbors_ptr: *const u32, neighbors_ptr: *const u32, 
    all_particles_len: usize, fluid_particles_len: usize, // TODO : ここの引数がわかりにくすぎるので変える
    h: f32, rest_density: f32, eps: f32, 
    blob_inv_m: f32, max_neighbors: u32
) {
    let pos_x = unsafe { std::slice::from_raw_parts(pos_x_ptr, all_particles_len) };
    let pos_y = unsafe { std::slice::from_raw_parts(pos_y_ptr, all_particles_len) };
    let lambdas = unsafe { std::slice::from_raw_parts_mut(lambdas_ptr, all_particles_len) };
    let grads_x = unsafe { std::slice::from_raw_parts_mut(grads_x_ptr, all_particles_len) };
    let grads_y = unsafe { std::slice::from_raw_parts_mut(grads_y_ptr, all_particles_len) };
    let fluid_deltas_x = unsafe { std::slice::from_raw_parts_mut(fluid_deltas_x_ptr, all_particles_len) };
    let fluid_deltas_y = unsafe { std::slice::from_raw_parts_mut(fluid_deltas_y_ptr, all_particles_len) };
    let first_neighbors = unsafe { std::slice::from_raw_parts(first_neighbors_ptr, all_particles_len + 1) };
    let neighbors = unsafe { std::slice::from_raw_parts(neighbors_ptr, all_particles_len * max_neighbors as usize) }; // ここはちゃんと変える

    let mut inside = vec![false; all_particles_len];
    for i in 0..all_particles_len {
        if lambdas[i] > 0.5 {
            inside[i] = true;
        } else {
            inside[i] = false;
        }
    }

    let h2 = h * h;
    let pi = 3.141592;
    let poly6_scale = 4.0 / (pi * h2 * h2 * h2 * h2);
    let spiky_scale = 10.0 / (pi * h2 * h2 * h);
    let rest_density_inv = 1.0 / rest_density;

    lambdas.par_iter_mut()
    .zip(grads_x.par_iter_mut())
    .zip(grads_y.par_iter_mut())
    .enumerate()
    .for_each(|(i, ((lambda, grad_x), grad_y))|{
        let first = first_neighbors[i];
        let num = first_neighbors[i + 1] - first;
        let is_inside_for_i = i < fluid_particles_len && inside[i];

        let mut rho = poly6_scale * h2 * h2 * h2;
        let mut sum_grad_2 = 0.;
        let mut grad_ix = 0.;
        let mut grad_iy = 0.;

        if !is_inside_for_i {
            for offset in 0..num {
                let j = neighbors[(first + offset) as usize] as usize;
                let is_inside_for_j = j < fluid_particles_len && inside[j];
                if i >= fluid_particles_len && j >= fluid_particles_len { // ignore softbody & softbody interaction
                    continue;
                }

                if is_inside_for_j {
                    continue; 
                }

                let dx: f32 = pos_x[i] - pos_x[j];
                let dy: f32 = pos_y[i] - pos_y[j];
                let r2 = dx * dx + dy * dy;     
                if r2 < h2 && 0. < r2 {
                    rho += poly6_scale * (h2 - r2) * (h2 - r2) * (h2 - r2);
                    let r = r2.sqrt();
                    let w = h - r;
                    let cur_grad_x = dx / r * (spiky_scale * w * w * -3.) * rest_density_inv;
                    let cur_grad_y = dy / r * (spiky_scale * w * w * -3.) * rest_density_inv;
                    sum_grad_2 += cur_grad_x * cur_grad_x + cur_grad_y * cur_grad_y;
                    grad_ix += cur_grad_x;
                    grad_iy += cur_grad_y;
                }
            }
        }

        let inv_m = if i >= fluid_particles_len { blob_inv_m } else { 1.0 };
        sum_grad_2 += (grad_ix * grad_ix + grad_iy * grad_iy) * inv_m;
        let c = rho / rest_density - 1.;
        *lambda = -c / (sum_grad_2 + eps);
        if c < 0. {
            *lambda = 0.;
        }
        *grad_x = -grad_ix;
        *grad_y = -grad_iy;
    });        


    fluid_deltas_x.par_iter_mut()
    .zip(fluid_deltas_y.par_iter_mut())
    .enumerate()
    .for_each(|(i, (fluid_delta_x, fluid_delta_y))|{
        let first = first_neighbors[i];
        let num = first_neighbors[i + 1] - first;

        let mut delta_x = 0.;
        let mut delta_y = 0.;
        let inv_m = if i >= fluid_particles_len { blob_inv_m } else { 1.0 };

        let is_inside_for_i = i < fluid_particles_len && inside[i];

        if i >= fluid_particles_len {
            for offset in 0..num {
                let j = neighbors[(first + offset) as usize] as usize;
                let is_inside_for_j = j < fluid_particles_len && inside[j];
                if j >= fluid_particles_len || is_inside_for_j {
                    continue;
                }

                let dx: f32 = pos_x[i] - pos_x[j];
                let dy: f32 = pos_y[i] - pos_y[j];
                let r2 = dx * dx + dy * dy;     
                if r2 < h2 && 0. < r2 {
                    let r = r2.sqrt();
                    let w = h - r;
                    let cur_grad_x = dx / r * (spiky_scale * w * w * -3.) * rest_density_inv;
                    let cur_grad_y = dy / r * (spiky_scale * w * w * -3.) * rest_density_inv;
                    delta_x += cur_grad_x * (lambdas[i] + lambdas[j]);
                    delta_y += cur_grad_y * (lambdas[i] + lambdas[j]);
                }
            }
        } else if !is_inside_for_i {
            for offset in 0..num {
                let j = neighbors[(first + offset) as usize] as usize;
                let is_inside_for_j = j < fluid_particles_len && inside[j];
                let dx: f32 = pos_x[i] - pos_x[j];
                let dy: f32 = pos_y[i] - pos_y[j];
                let r2 = dx * dx + dy * dy;     

                if is_inside_for_j {
                    continue; 
                }

                if r2 < h2 && 0. < r2 {
                    let r = r2.sqrt();
                    let w = h - r;
                    let cur_grad_x = dx / r * (spiky_scale * w * w * -3.) * rest_density_inv;
                    let cur_grad_y = dy / r * (spiky_scale * w * w * -3.) * rest_density_inv;
                    if j >= fluid_particles_len {
                        delta_x += cur_grad_x * (lambdas[i] + lambdas[j]);
                        delta_y += cur_grad_y * (lambdas[i] + lambdas[j]);
                    } else {
                        let cohesion_x = dx / r * (-0.01 * 0.01);
                        let cohesion_y = dy / r * (-0.01 * 0.01);
                        let curvature_x = (grads_x[i] - grads_x[j]) * spiky_scale * w * w * w * 0.001 * 0.02;
                        let curvature_y = (grads_y[i] - grads_y[j]) * spiky_scale * w * w * w * 0.001 * 0.02;
                        delta_x += cur_grad_x * (lambdas[i] + lambdas[j]) + cohesion_x + curvature_x;
                        delta_y += cur_grad_y * (lambdas[i] + lambdas[j]) + cohesion_y + curvature_y;
                    }
                }
            }
        }

        *fluid_delta_x = delta_x * inv_m;
        *fluid_delta_y = delta_y * inv_m;
    });
}