use crate::functions::kmeans;
use rand::prelude::*;
use std::f32::consts::E;


// Returns the value of within cluster distance for k clusters
fn wcd(data_array: &Vec<Vec<f32>>, k: usize) -> f32 {
    // kmeans
    let data: Vec<Vec<f32>> = data_array.clone();
    let centroids = kmeans::initialize(k, &data);
    let clustering = kmeans::kmeans(k, &data, 200, centroids);

    // within cluster distance
    let mut wk: f32 = 0.0;
    for i in 0..k {
        let mut centroid_dist: f32 = 0.0;
        //let mut clust_nb = 0.0;
        let wedf_centroid = &clustering.centroids[i].0;
        for j in 0..data_array.len() {
            if clustering.membership[j] == i {
                let mut norm = 0.0;
                for dim in 0..data_array[0].len() {
                    norm += (data_array[j][dim] - wedf_centroid[dim] as f32)*(data_array[j][dim] - wedf_centroid[dim] as f32);
                }
                centroid_dist += norm;
                //clust_nb += 1.0;
            }
        }
        //wk += centroid_dist/clust_nb;
        wk += centroid_dist;
    }
    wk = wk.log(E);
    
    wk
}


// Returns the minimum and maximum values on each dimension of a dataset
fn bounding_box(data_array: &Vec<Vec<f32>>) -> (Vec<f32>, Vec<f32>) {
    let mut max = data_array[0].clone();
    let mut min = data_array[0].clone();

    for i in 0..data_array.len() {
        for dim in 0..data_array[0].len() {
            if data_array[i][dim] > max[dim] {
                max[dim] = data_array[i][dim];
            } else if data_array[i][dim] < min[dim] {
                min[dim] = data_array[i][dim];
            }
        }
    }

    (min, max)
}


// Returns the mean of Wk and Sk for differents simulations
fn simulate(data_array: &Vec<Vec<f32>>, iter: usize, k: usize) -> (f32, f32) {
    // Maximum  and minimum of data
    let (min_data, max_data) = bounding_box(data_array);

    let mut mean_wk: f32 = 0.0;
    let mut sk: f32 = 0.0;
    let mut wks: Vec<f32> = Vec::new();

    for _i in 0..iter {
        let mut random_data: Vec<Vec<f32>> = Vec::new();
        for _j in 0..data_array.len() {
            let mut one_random_data: Vec<f32> = Vec::new();
            for dim in 0..data_array[0].len() {
                let mut rng = rand::thread_rng();
                let y: f32 = rng.gen();
                one_random_data.push(y*(max_data[dim]-min_data[dim])+min_data[dim]);
            }
            random_data.push(one_random_data);
        }
        let wk = wcd(&random_data, k);
        wks.push(wk);
    }

    // mean
    for i in 0..iter {
        mean_wk += wks[i];
    }
    mean_wk = mean_wk/ iter as f32;

    // standard deviation
    for i in 0..iter {
        sk += (wks[i] - mean_wk)*(wks[i] - mean_wk);
    }
    sk = (sk/ iter as f32).sqrt();
    let coef: f32 = 1.0/iter as f32 + 1.0;
    sk = coef.sqrt()*sk;

    (mean_wk, sk)
}


// Gap evaluation algorithm
// Returns the optimal k for kmeans within a list of k
pub fn gap_eval(data_array: &Vec<Vec<f32>>, k_min: usize, k_max: usize) -> usize {
    let mut gaps: Vec<f32> = Vec::new();
    let mut sks: Vec<f32> = Vec::new();
    for k in k_min..k_max+1 {
        let wk = wcd(data_array, k);
        let (wk_uniform, sk) = simulate(data_array, 100, k);
        sks.push(sk);
        let gap = wk_uniform - wk;
        gaps.push(gap as f32);
    }

    // Find the optimal K for the gap statistic
    let mut k_opt = k_max;
    let ks: usize = k_max-k_min;
    for i in 0..ks {
        if gaps[i] >= gaps[i+1] - sks[i+1] {
            k_opt = k_min + i;
            break;
        }
    }
    k_opt
}