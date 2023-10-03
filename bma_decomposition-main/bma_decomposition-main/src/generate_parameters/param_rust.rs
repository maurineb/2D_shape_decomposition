use std::fs::File;
use std::io::{self, BufRead};
use nalgebra::base::*;


// The functions below are required to extract the objects for the bma decomposition
// from the output files of the clean medial axis in 2D


// Function in order to extract the coordinates of the points and the indicies of the edges of the boundary
pub fn bnd_to_vec(filename: &str) -> (Vec<Vector2<f32>>, Vec<Vector2<i32>>) {
    let mut bnd_pts = Vec::new();
    let mut bnd_edges = Vec::new();
    let mut is_pt = false;
    let mut is_edge = false;

    let file = File::open(filename).unwrap();
    let lines = io::BufReader::new(file).lines();

    for line_ in lines {
        if let Ok(line) = line_ {
            if line == "%points" {
                is_pt = true;
                is_edge = false;
            }
            else if line == "%edges" {
                is_edge = true;
                is_pt = false;
            }
            else if is_pt {
            // This line is the coordinates of a boundary point
                let mut line_split = line.split_whitespace();
                let mut vert: Vector2<f32> = Vector2::new(0.0, 0.0);
                for i in 0..2 {
                    match line_split.next() {
                        Some(cur) => {
                            vert[i] = cur.parse::<f32>().unwrap();
                        },
                        None=> {
                        }
                    }
                }
                bnd_pts.push(vert);
            }
            else if is_edge {
            // This line is the indicies of a boundary edge
                let mut line_split = line.split_whitespace();
                let mut vert: Vector2<i32> = Vector2::new(0, 0);
                for i in 0..2 {
                    match line_split.next() {
                        Some(cur) => {
                            vert[i] = cur.parse::<i32>().unwrap();
                        },
                        None=> {
                        }
                    }
                }
                bnd_edges.push(vert);
            }
        }
    }
    bnd_pts.remove(bnd_pts.len()-1);
    
    (bnd_pts, bnd_edges)
}


// Function in order to extract the coordinates of the points and the radius of the medial axis
pub fn nod_to_vec(filename: &str) -> (Vec<Vector2<f32>>, Vec<f32>) {
    let mut ma_pts = Vec::new();
    let mut ma_rad = Vec::new();

    let file = File::open(filename).unwrap();
    let lines = io::BufReader::new(file).lines();

    for line_ in lines {
        if let Ok(line) = line_ {
            let mut line_split = line.split_whitespace();
            let mut vert: Vector2<f32> = Vector2::new(0.0, 0.0);
            for i in 0..2 {
                match line_split.next() {
                    Some(cur) => {
                        vert[i] = cur.parse::<f32>().unwrap();
                    },
                    None=> {
                    }
                }
            }
            ma_pts.push(vert);
            match line_split.next() {
                Some(rad) => {
                    ma_rad.push(rad.parse::<f32>().unwrap());
                },
                None=> {
                }
            }
        }
    }
    
    (ma_pts, ma_rad)
}


// Function in order to extract the indicies of the edges of the medial axis
pub fn edg_to_vec(filename: &str) -> Vec<Vector2<i32>> {
    let mut ma_edges = Vec::new();

    let file = File::open(filename).unwrap();
    let lines = io::BufReader::new(file).lines();

    for line_ in lines {
        if let Ok(line) = line_ {
            let mut line_split = line.split_whitespace();
            let mut vert: Vector2<i32> = Vector2::new(0, 0);
            for i in 0..2 {
                match line_split.next() {
                    Some(cur) => {
                        vert[i] = cur.parse::<i32>().unwrap();
                    },
                    None=> {
                    }
                }
            }
            ma_edges.push(vert);
        }
    }
    
    ma_edges
}


// Function in order to extract the indicies of the boundary points associated with the medial axis points
pub fn del_to_vec(filename: &str) -> Vec<Vec<i32>> {
    let mut ind_bnd_pts: Vec<Vec<i32>> = Vec::new();

    let file = File::open(filename).unwrap();
    let lines = io::BufReader::new(file).lines();
    
    for line_ in lines {
        if let Ok(line) = line_ {
            let mut line_split = line.split_whitespace();
            let mut vert: Vec<i32> = Vec::new();

            // Consuming the first two elements since we don't need them
            line_split.next();
            line_split.next();
            let line_split_len = line_split.clone();

            for _ in 0..line_split_len.count() {
                match line_split.next() {
                    Some(cur) => {
                        // Deleting the , at the end of the value
                        let mut cur_string = cur.to_string();
                        cur_string.pop();

                        vert.push(cur_string.as_str().parse::<i32>().unwrap());
                    },
                    None=> {
                    }
                }
            }
            ind_bnd_pts.push(vert);
        }
    }
    
    ind_bnd_pts
}