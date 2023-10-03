use std::fmt::format;

use bma_decomposition::bma::bma;
use bma_decomposition::generate_parameters::param_rust;
use bma_decomposition::functions::functions_geometry;
use bma_decomposition::functions::functions_vec;

use clap::Parser;
//use plotters::prelude::*;
//use image::GenericImageView;
//use image::{DynamicImage, GrayImage, RgbImage};
//use image::io::Reader as ImageReader;
//use nalgebra::base::*;
use bma_decomposition::functions::kmeans;
use bma_decomposition::functions::gap_evaluation;
use nalgebra::Vector2;
use rand::prelude::*;


#[derive(Parser)]
struct Cli {
    #[arg(default_value = "../../compact-skel-2d/ressources", long = "folder_skel")]
    skel_path: std::path::PathBuf,
    #[arg(default_value = "./output_files", long = "folder_results")]
    folder_results: std::path::PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    let skel_path_str = args.skel_path.to_str().unwrap_or("");
    let folder_results_path = args.folder_results.to_str().unwrap_or("");

    let bnd_in_path_str = format!("{}{}", skel_path_str, "/bnd.txt");
    let nod_in_path_str = format!("{}{}", skel_path_str, "/nod.txt");
    let edg_in_path_str = format!("{}{}", skel_path_str, "/edg.txt");
    let del_in_path_str = format!("{}{}", skel_path_str, "/del.txt");

    let (bnd_pts, bnd_edges) = param_rust::bnd_to_vec(&bnd_in_path_str);
    let (ma_pts, ma_rad) = param_rust::nod_to_vec(&nod_in_path_str);
    let ma_edges = param_rust::edg_to_vec(&edg_in_path_str);
    let ind_bnd_pts = param_rust::del_to_vec(&del_in_path_str);

    let mut bma = bma::bma::new(&bnd_pts, &bnd_edges, &ma_pts, &ma_edges, &ind_bnd_pts, &ma_rad);

    // for every folder in format for the code in matlab
    let image_skel_path_str = format!("{}{}", skel_path_str, "/skeleton.png");
    // for the rat
    //let image_skel_path_str = "../../compact-skel-2d/ressources/rat.png";
    let matlab_version = false;

    /**/
    // TO SHOW THE SHAPE AND MEDIAL AXIS
    let image_shape_path_str = format!("{}{}", folder_results_path, "/bma/shape.png");
    let image_medial_axis_path_str = format!("{}{}", folder_results_path, "/bma/medial_axis.png");
    bma.plot_shape(&image_skel_path_str, &image_shape_path_str, false);
    bma.plot_shape(&image_skel_path_str, &image_medial_axis_path_str, true);


    // TO SHOW THE BOUNDARY
    let image_bnd_outer_path_str = format!("{}{}", folder_results_path, "/bma/boundary_outer.png");
    let image_bnd_component_path_str = format!("{}{}", folder_results_path, "/bma/boundary_component.png");
    bma.plot_boundary(&image_skel_path_str, &image_bnd_outer_path_str, false);
    bma.plot_boundary(&image_skel_path_str, &image_bnd_component_path_str, true);


    // TO SHOW THE MEDIAL AXIS
    let image_bma_loop_path_str = format!("{}{}", folder_results_path, "/bma/bma_loop.png");
    let image_bma_pointtype_path_str = format!("{}{}", folder_results_path, "/bma/bma_point_type.png");
    let image_bma_branches_path_str = format!("{}{}", folder_results_path, "/bma/bma_branches.png");
    let image_bma_branch_adjacency_path_str = format!("{}{}", folder_results_path, "/bma/bma_branch_adjacency.png");

    bma.plot_bma_loops(&image_skel_path_str, &image_bma_loop_path_str);
    bma.plot_bma_pointtypes(&image_skel_path_str, &image_bma_pointtype_path_str);
    bma.plot_bma_branches(&image_skel_path_str, &image_bma_branches_path_str);
    bma.plot_bma_adjacency_branches(1, &image_skel_path_str, &image_bma_branch_adjacency_path_str);


    let image_triangles_path_str = format!("{}{}", folder_results_path, "/bma/triangles.png");
    bma.plot_triangles(&image_skel_path_str, &image_triangles_path_str);
    //bma.plot_triangle_i(&image_skel_path_str, &image_triangles_path_str, 222);
    //bma.plot_triangle_of_ma_pt(144, &image_skel_path_str, &image_triangles_path_str);
    

    // TO SHOW THE WEDF VALUES
    let image_wedf_path_str = format!("{}{}", folder_results_path, "/values/wedf.png");
    if matlab_version {
        bma.plot_wedf_matlab(&image_skel_path_str, &image_wedf_path_str);
    } else {
        bma.plot_wedf(&image_skel_path_str, &image_wedf_path_str);
    }
    

    // TO SHOW THE EDF VALUES
    let image_edf_path_str = format!("{}{}", folder_results_path, "/values/edf.png");
    if matlab_version {
        bma.plot_edf_matlab(&image_skel_path_str, &image_edf_path_str);
    } else {
        bma.plot_edf(&image_skel_path_str, &image_edf_path_str);
    }
    

    // TO SHOW THE EROSION THICKNESS VALUES
    let image_erosion_thickness_path_str = format!("{}{}", folder_results_path, "/values/erosion_thickness.png");
    if matlab_version {
        bma.plot_erosion_thickness_matlab(&image_skel_path_str, &image_erosion_thickness_path_str);
    } else {
        bma.plot_erosion_thickness(&image_skel_path_str, &image_erosion_thickness_path_str);
    }
    

    // TO SHOW THE SHAPE TUBULARITY VALUES
    let image_shape_tubularity_path_str = format!("{}{}", folder_results_path, "/values/shape_tubularity.png");
    if matlab_version {
        bma.plot_shape_tubularity_matlab(&image_skel_path_str, &image_shape_tubularity_path_str);
    } else {
        bma.plot_shape_tubularity(&image_skel_path_str, &image_shape_tubularity_path_str);
    }
    

    // TO SHOW THE WEDF HIERARCHY PARTS DECOMPOSITION
    let image_hierarchy_first_path_str = format!("{}{}", folder_results_path, "/decomposition/hierarchy_parts_before_spread.png");
    let image_hierarchy_path_str = format!("{}{}", folder_results_path, "/decomposition/hierarchy_parts.png");

    let (before_spread_hierlabel, centroids1, _centroids2, _optimal_k) = bma.hierarchy_wedf(3, 2, 8, 12);
    bma.plot_wedf_hierarchy(&image_skel_path_str, &image_hierarchy_path_str);
    bma.plot_hierarchy(&image_skel_path_str, &image_hierarchy_first_path_str, &before_spread_hierlabel);


    /*
    let image_out = format!("{}{}", folder_results_path, "/test.png");
    let mut test_bma = bma.clone();
    test_bma.medial_points.clear();
    test_bma.medial_points = vec![Vector2::new(12.0, 45.0), Vector2::new(26.0, 110.0), Vector2::new(137.0, 110.0), Vector2::new(112.0, 203.0)];
    test_bma.test(&image_skel_path_str, &image_out);
    */
    /**/
    
    // TO SHOW THE TRUNKS
    bma.shape_details_hierarchy();

    let image_trunks_path_str = format!("{}{}", folder_results_path, "/bma/trunks.png");
    //bma.plot_bma_one_trunk(2, 0, &image_skel_path_str, &image_trunks_path_str);
    //bma.plot_bma_trunks(2, &image_skel_path_str, &image_trunks_path_str);
    bma.plot_bma_all_trunks(&image_skel_path_str, &image_trunks_path_str);
    

    // TO SHOW THE PARTS DECOMPOSITION AFTER WEDF HIERARCHY BUT BEFORE COMPUTING THE TRUNKS
    let image_new_levels_path_str = format!("{}{}", folder_results_path, "/decomposition/new_levels_parts.png");

    bma.plot_new_cuts(&image_skel_path_str, &image_new_levels_path_str);


    // TO SHOW THE DETAILS PARTS DECOMPOSITION
    if !bma.trunk.is_empty() {
        let image_clustering_trunks_path_str = format!("{}{}", folder_results_path, "/decomposition/clustering_trunks.png");
        let image_cluster_1_trunks_path_str = format!("{}{}", folder_results_path, "/decomposition/cluster_1_trunks.png");
        let image_cluster_2_trunks_path_str = format!("{}{}", folder_results_path, "/decomposition/cluster_2_trunks.png");
        let image_cluster_3_trunks_path_str = format!("{}{}", folder_results_path, "/decomposition/cluster_3_trunks.png");
        let image_cluster_4_trunks_path_str = format!("{}{}", folder_results_path, "/decomposition/cluster_4_trunks.png");
        let image_cluster_5_trunks_path_str = format!("{}{}", folder_results_path, "/decomposition/cluster_5_trunks.png");
        
        let (centroids, labels) = bma.clustering_etst(1, 5, false);
        //println!("centroids : {:?}", centroids);
        //bma.plot_clusters_trunk(&image_skel_path_str, &image_clustering_trunks_path_str, &labels, 7);
        bma.plot_all_clusters_trunk(&image_skel_path_str, &image_clustering_trunks_path_str, &labels);
        bma.plot_one_cluster_trunk(&image_skel_path_str, &image_cluster_1_trunks_path_str, &labels, 0);
        let clusters = functions_vec::unique_list(&labels);
        if clusters.len() > 1 {
            bma.plot_one_cluster_trunk(&image_skel_path_str, &image_cluster_2_trunks_path_str, &labels, 1);
        }
        if clusters.len() > 2 {
            bma.plot_one_cluster_trunk(&image_skel_path_str, &image_cluster_3_trunks_path_str, &labels, 2);
        }
        if clusters.len() > 3 {
            bma.plot_one_cluster_trunk(&image_skel_path_str, &image_cluster_4_trunks_path_str, &labels, 3);
        }
        if clusters.len() > 4 {
            bma.plot_one_cluster_trunk(&image_skel_path_str, &image_cluster_5_trunks_path_str, &labels, 4);
        }
    }

    Ok(())
}