use nalgebra::base::*;

// Return the indices of the list containing element (Vec<Vector2>)
pub fn find_vecvector2<T: std::cmp::PartialEq>(list: &Vec<Vector2<T>>, element: T) -> (Vec<usize>, Vec<usize>) {
    let mut x: Vec<usize> = Vec::new();
    let mut y: Vec<usize> = Vec::new();

    for i in 0..list.len() {
        for j in 0..list[0].len() {
            if list[i][j] == element {
                x.push(i as usize);
                y.push(j as usize);
            }
        }
    }
    (x, y)
}

// Return the indices of the list containing element (Vec<Vec>)
pub fn find_vecvec<T: std::cmp::PartialEq>(list: &Vec<Vec<T>>, element: T) -> (Vec<usize>, Vec<usize>) {
    let mut x: Vec<usize> = Vec::new();
    let mut y: Vec<usize> = Vec::new();

    for i in 0..list.len() {
        for j in 0..list[i].len() {
            if list[i][j] == element {
                x.push(i as usize);
                y.push(j as usize);
            }
        }
    }
    (x, y)
}

// Return the indices of the list containing element (Vec)
pub fn find_vec<T: std::cmp::PartialEq>(list: &Vec<T>, element: T) -> Vec<usize> {
    let mut ind: Vec<usize> = Vec::new();

    for i in 0..list.len() {
        if list[i] == element {
            ind.push(i as usize);
        }
    }
    ind
}


// Return the unique values of a list
pub fn unique_list<T: Ord + Copy>(list: &Vec<T>) -> Vec<T> {
    let mut sorted_list = list.clone();
    sorted_list.sort();
    // let unique_list.dedup();
    let mut unique_list: Vec<T> = Vec::new();
    if sorted_list.is_empty() {
        return unique_list;
    }
    unique_list.push(sorted_list[0]);
    for i in 0..sorted_list.len() {
        if sorted_list[i] != unique_list[unique_list.len()-1] {
            unique_list.push(sorted_list[i]);
        }
    }
    unique_list
}


// Find the common values between two vectors
pub fn intersect<T: Ord + Copy>(vec_1: &Vec<T>, vec_2: &Vec<T>) -> Vec<T> {
    let mut common: Vec<T> = Vec::new();
    for i in 0..vec_1.len() {
        let ind_comm = find_vec(&vec_2, vec_1[i]);
        for j in ind_comm {
            common.push(vec_2[j]);
        }
    }
    unique_list(&common)
}


// Find the minimum value and index of this value in a list of f32
pub fn min_f32(vec: &Vec<f32>) -> (f32, Vec<usize>) {
    assert!(vec.len()>0, "The vector can't be empty.");
    let mut min = vec[0];
    for i in 1..vec.len() {
        if vec[i] < min {
            min = vec[i];
        }
    }
    let ind = find_vec(&vec, min);
    (min, ind)
}


// Find the maximum value and index of this value in a list of f32
pub fn max_f32(vec: &Vec<f32>) -> (f32, Vec<usize>) {
    assert!(vec.len()>0, "The vector can't be empty.");
    let mut max = vec[0];
    for i in 1..vec.len() {
        if vec[i] > max {
            max = vec[i];
        }
    }
    let ind = find_vec(&vec, max);
    (max, ind)
}


// Find the indices that sort a list
pub fn argsort<T: Ord>(data: &[T]) -> Vec<usize> {
    let mut indices = (0..data.len()).collect::<Vec<_>>();
    indices.sort_by_key(|&i| &data[i]);
    indices
}


// Check for counterclockwiseness
pub fn counterclockwiseness(list_points: &Vec<Vector2<f32>>) -> bool {
    let mut sum = 0.0;
    for i in 0..(list_points.len()-1) {
        sum += (list_points[i+1][0] - list_points[i][0])*(list_points[i+1][1] + list_points[i][1]);
    }
    sum += (list_points[0][0] - list_points[list_points.len()-1][0])*(list_points[0][1] + list_points[list_points.len()-1][1]);
    sum <= 0.0
}