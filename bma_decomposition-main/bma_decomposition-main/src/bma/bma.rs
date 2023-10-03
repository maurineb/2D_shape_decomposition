use std::alloc::LayoutError;

//use anyhow::Result;
use nalgebra::base::*;
use image::GenericImageView;
use imageproc::point::Point;
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

use crate::functions::functions_vec;
use crate::functions::functions_geometry;
use crate::functions::kmeans;
use crate::functions::gap_evaluation;



#[derive(Clone)]
pub struct Trunk {
    level: i32,
    component: usize,
    indices: Vec<usize>,
    branchpath: Vec<usize>,
    parent_child_array: Vec<Vec<usize>>
}


#[derive(Clone)]
pub struct bma {
    // coordinates (x,y) of the boundary points
    boundary_points: Vec<Vector2<f32>>, 
    // indices of the boundary edges
    boundary_edges: Vec<Vector2<i32>>, 
    // boolean to indicate if the boundary point is on the outer border
    isouter: Vec<bool>, 
    // integer to know on which border the point is (number randomly chosen)
    component: Vec<i32>, 
    // coordinates (x,y) of the medial axis points
    medial_points: Vec<Vector2<f32>>, 
    // radius values associated with each medial axis points
    medial_radius: Vec<f32>, 
    // indices of the medial axis edges
    medial_edges: Vec<Vector2<i32>>, 
    // indicates whether medial axis point is regular, triple, single
    point_type: Vec<i32>, 
    // boolean indicating the neighbors of triple points
    istrip_nbr: Vec<bool>, 
    // boolean indicating if the boundary points associated with a ma pt are on two different components
    onloop: Vec<bool>, 
    // gives indices into boundary corr. to ma pts
    index_bndry_points: Vec<Vec<i32>>, 
    // if the points are next to each other
    adjacency_matrix: Vec<Vec<bool>>, 
    // each column gives the order of the points of each branch (0 when not in the branch)
    branch_number: Vec<Vec<i32>>, 
    // if the branches are next to each other
    branch_adjacency: Vec<Vec<bool>>,
    // edf array
    edf_array: Vec<f32>,
    // wedf array
    wedf_array: Vec<f32>,
    // gives max EDF for each branch, including distance to triple point (for storing discontinuous EDF values)
    edf_branch: Vec<f32>,
    // gives max WEDF for each branch
    wedf_branch: Vec<f32>,
    // erosion thickness array
    erosion_thickness: Vec<f32>,
    // shape tubularity
    shape_tubularity: Vec<f32>,
    // triangulation of the shape where each row is a triangle
    // each row is : [index of the corresponding ma point, index of the first point, index of the second point, index of the third point]
    triangles: Vec<Vec<i32>>,
    // TODO : add properties for shape details hierarchy
    // hierarchy label (hierlabel output from hierarchyWEDF)
    clusters_wedf: Vec<i32>,
    // pts on which the hierarchy was clustered (nubinds output from hierarchyWEDF)
    clusters_ind: Vec<usize>,
    // assigned level after cuts are made
    hierarchy_level: Vec<i32>,
    // ??? WEDF values at cuts
    level_values: Vec<i32>,
    // structure giving details of each medial fragment at each hierarchy level
    pub trunk: Vec<Trunk>,
    // # medial points x # medial fragments array indicating trunk membership of that point
    trunk_members: Vec<Vec<bool>>
}


impl bma {
    // Constructor for bma
    pub fn new(
        bndy_pts: &Vec<Vector2<f32>>,
        bndy_edges: &Vec<Vector2<i32>>,
        bma_pts: &Vec<Vector2<f32>>,
        bma_edges: &Vec<Vector2<i32>>,
        del: &Vec<Vec<i32>>,
        bma_rad: &Vec<f32>
    ) -> bma {
        let mut bma_output = bma {
            boundary_points: bndy_pts.clone(),
            boundary_edges: bndy_edges.clone(),
            isouter: Vec::new(),
            component: Vec::new(),
            medial_points: bma_pts.clone(),
            medial_radius: bma_rad.clone(),
            medial_edges: bma_edges.clone(),
            point_type: Vec::new(),
            istrip_nbr: Vec::new(),
            onloop: Vec::new(),
            index_bndry_points: del.clone(),
            adjacency_matrix: Vec::new(),
            branch_number: Vec::new(),
            branch_adjacency: Vec::new(),
            edf_array: Vec::new(),
            wedf_array: Vec::new(),
            edf_branch: Vec::new(),
            wedf_branch: Vec::new(),
            erosion_thickness: Vec::new(),
            shape_tubularity: Vec::new(),
            triangles: Vec::new(),
            clusters_wedf: Vec::new(),
            clusters_ind: Vec::new(),
            hierarchy_level: Vec::new(),
            level_values: Vec::new(),
            trunk: Vec::new(),
            trunk_members: Vec::new()
        };
        let res = bma_output.build_boundary();
        assert!(res.is_ok());
        let res = bma_output.build_bma();
        assert!(res.is_ok());
        let res = bma_output.branches_for_bma();
        assert!(res.is_ok());
        let res = bma_output.extract_triangles();
        assert!(res.is_ok());
        println!("bma created");
        let res = bma_output.calculate_wedf(true);
        println!("wedf computed");
        assert!(res.is_ok());
        let res = bma_output.calculate_edf();
        println!("edf computed");
        assert!(res.is_ok());
        let res = bma_output.calculate_et();
        println!("erosion thickness computed");
        assert!(res.is_ok());
        let res = bma_output.calculate_st();
        println!("shape tubularity computed");
        assert!(res.is_ok());

        bma_output
    }


    fn build_boundary(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Extracting the inputs from the bma object
        let bndy_pts = &self.boundary_points;
        let bndy_edges = &self.boundary_edges;
        let ind_bnd_pts = &self.index_bndry_points;
        
        let mut bndy_edges_mut = bndy_edges.clone();
        
        // Ordered list of points according to the edges
        let mut ind_order: Vec<i32> = Vec::new();
        // Current edge that we are adding to the ordered list of points
        let mut curr_ind: Vec<i32> = vec![bndy_edges_mut[0][0], bndy_edges_mut[0][1]];
        ind_order.push(curr_ind[0]);
        // Outputs
        let mut isouter: Vec<bool> = vec![false; bndy_pts.len()];
        let mut component: Vec<i32> = vec![0; bndy_pts.len()];

        let mut allinds: Vec<i32> = Vec::new();
        let mut flip = false;
        let mut comp = 1;
        bndy_edges_mut.remove(0);

        while !curr_ind.is_empty() {
            let (ind1, ind2) = functions_vec::find_vecvector2(&bndy_edges_mut, curr_ind[1]);
            let _ind = functions_vec::find_vec(&ind_order, curr_ind[1]);

            // if it's the point wasn't registered and is found in the edges
            let save_point = _ind.is_empty() && !ind1.is_empty();
            if save_point {
                ind_order.push(curr_ind[1]);
                curr_ind = vec![bndy_edges_mut[ind1[0]][0], bndy_edges_mut[ind1[0]][1]];
                if ind2[0] == 1 {
                    curr_ind.reverse();
                }
                bndy_edges_mut.remove(ind1[0]);
            } else {
                let mut temp_bndy: Vec<Vector2<f32>> = Vec::new();
                for i in 0..ind_order.len() {
                    temp_bndy.push(bndy_pts[ind_order[i] as usize]);
                }
                // Check for counterclockwiseness
                // (not supposed to go in there with the code for the skeletton)
                if !functions_vec::counterclockwiseness(&temp_bndy) {
                    temp_bndy.reverse();
                    ind_order.reverse();
                    flip = true;
                }
                // Create a list with all the points that aren't on this boundary
                let mut otherbndy_pts: Vec<Vector2<f32>> = Vec::new();
                for i in 0..bndy_pts.len() {
                    if !ind_order.contains(&(i as i32)) {
                        otherbndy_pts.push(bndy_pts[i]);
                    }
                }
                // Find if this boundary is the outer boundary
                let poly = functions_geometry::Polygon::new(&temp_bndy);
                let mut inside: i32 = 0;
                if otherbndy_pts.is_empty() {
                    inside = 10;
                } else {
                    for pt in &otherbndy_pts {
                        if poly.contains_point(pt[0], pt[1]) {
                            inside += 1;
                        }
                    }
                }

                if inside>1 {
                    for i in &ind_order {
                        isouter[*i as usize] = true;
                    }
                } else if flip {
                    ind_order.reverse();
                }
                for val in &ind_order {
                    allinds.push(*val);
                }
                for i in &ind_order {
                    component[*i as usize] = comp;
                }
                ind_order.clear();
                flip = false;
                if !bndy_edges_mut.is_empty() {
                    curr_ind = vec![bndy_edges_mut[0][0], bndy_edges_mut[0][1]];
                    bndy_edges_mut.remove(0);
                    ind_order.push(curr_ind[0]);
                } else {
                    curr_ind.clear();
                }
                comp += 1;
            }
        }
        
        let mut new_bndy_edges = bndy_edges.clone();
        for i in 0..new_bndy_edges.len() {
            let newind1 = functions_vec::find_vec(&allinds, bndy_edges[i][0]);
            let newind2 = functions_vec::find_vec(&allinds, bndy_edges[i][1]);
            new_bndy_edges[i] = Vector2::new(newind1[0] as i32, newind2[0] as i32);
        }

        let mut new_ind_bnd_pts: Vec<Vec<i32>> = Vec::new();
        for i in 0..ind_bnd_pts.len() {
            let mut newset: Vec<i32> = Vec::new();
            for j in 0..ind_bnd_pts[i].len() {
                let newind = functions_vec::find_vec(&allinds, ind_bnd_pts[i][j]);
                newset.push(newind[0] as i32);
            }
            new_ind_bnd_pts.push(newset);
        }

        let mut new_bndy_pts = bndy_pts.clone();
        for i in 0..bndy_pts.len() {
            new_bndy_pts[i] = bndy_pts[allinds[i] as usize];
        }

        let mut new_isouter = isouter.clone();
        for i in 0..isouter.len() {
            new_isouter[i] = isouter[allinds[i] as usize];
        }

        let mut new_component = component.clone();
        for i in 0..component.len() {
            new_component[i] = component[allinds[i] as usize];
        }

        self.boundary_points = new_bndy_pts;
        self.isouter = new_isouter;
        self.component = new_component;
        self.boundary_edges = new_bndy_edges;
        self.index_bndry_points = new_ind_bnd_pts;

        Ok(())
    }


    fn build_bma(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Extracting the inputs from the bma object
        let mut bma_edges = self.medial_edges.clone();

        // Adjacency matrix
        let mut adj_matrix: Vec<Vec<bool>> = vec![vec![false; self.medial_points.len()]; self.medial_points.len()];
        for i in 0..bma_edges.len() {
            if bma_edges[i][0] > bma_edges[i][1] {
                let temp = bma_edges[i][0];
                bma_edges[i][0] = bma_edges[i][1];
                bma_edges[i][1] = temp;
            }
            adj_matrix[bma_edges[i][0] as usize][bma_edges[i][1] as usize] = true;
            adj_matrix[bma_edges[i][1] as usize][bma_edges[i][0] as usize] = true;
        }

        // Types of point
        let mut point_type: Vec<i32> = Vec::new();
        // If it is a neighbour to a multiple point
        let mut multiple_neighbours: Vec<bool> = vec![false; self.medial_points.len()];
        // 0 : regular points (with two neighbours)
        // 1 : single points (with one neighbour)
        // 3 : multiple points (with at least trhee neighbours)
        for i in 0..adj_matrix.len() {
            let mut num_neighbours = 0;
            for j in 0..adj_matrix.len() {
                if adj_matrix[i][j] {
                    num_neighbours += 1;
                }
            }
            if num_neighbours == 2 {
                num_neighbours = 0;
            } else if num_neighbours > 3 {
                num_neighbours = 3;
            }
            // Update multiple_neighbours
            if num_neighbours == 3 {
                for j in 0..adj_matrix.len() {
                    if adj_matrix[i][j] {
                        multiple_neighbours[j] = true;
                    }
                }
            }
            point_type.push(num_neighbours);
        }

        // If a point is on a loop
        let mut onloop: Vec<bool> = Vec::new();
        for i in 0..self.index_bndry_points.len() {
            let mut pt_i_onloop = false;
            let oncomponent = self.component[self.index_bndry_points[i][0] as usize];
            for j in 1..self.index_bndry_points[i].len() {
                pt_i_onloop = pt_i_onloop || !(self.component[self.index_bndry_points[i][j] as usize] == oncomponent);
            }
            onloop.push(pt_i_onloop);
        }

        self.medial_edges = bma_edges;
        self.point_type = point_type;
        self.adjacency_matrix = adj_matrix;
        self.istrip_nbr = multiple_neighbours;
        self.onloop = onloop;

        Ok(())
    }


    fn branches_for_bma(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // Finding the different branches
        let mut branches: Vec<Vec<i32>> = Vec::new();

        let ind_single = functions_vec::find_vec(&self.point_type, 1);
        let ind_multiple = functions_vec::find_vec(&self.point_type, 3);

        // If the medial axis is a single loop
        if ind_single.is_empty() && ind_multiple.is_empty() {
            // Starting point
            let mut branch: Vec<usize> = Vec::new();
            branch.push(0);
            // Find the following point
            let mut last_point = 0;
            let mut nbr = functions_vec::find_vec(&self.adjacency_matrix[last_point], true);
            let mut next_point = nbr[0];

            // Looking for the neighbours of the neighbour...
            while next_point != 0 {
                branch.push(next_point);
                nbr = functions_vec::find_vec(&self.adjacency_matrix[next_point], true);
                if nbr[0] == last_point {
                    last_point = next_point;
                    next_point = nbr[1];
                } else {
                    last_point = next_point;
                    next_point = nbr[0];
                }
            }
            let mut branch_temp: Vec<i32> = vec![0; self.medial_points.len()];
            for i in 0..branch.len() {
                branch_temp[branch[i]] = i as i32 +1;
            }
            branches.push(branch_temp);
        }

        // If the medial axis is a single branch
        else if ind_multiple.is_empty() {
            let nbr_single = functions_vec::find_vec(&self.adjacency_matrix[ind_single[0]], true);
            // Starting point
            let mut branch: Vec<usize> = Vec::new();
            branch.push(ind_single[0]);
            let mut next_point = nbr_single[0];
            let mut last_point = ind_single[0];
            
            // Looking for the neighbours of the neighbour...
            while self.point_type[next_point] == 0 {
                branch.push(next_point);
                let nbr = functions_vec::find_vec(&self.adjacency_matrix[next_point], true);
                if nbr[0] == last_point {
                    last_point = next_point;
                    next_point = nbr[1];
                } else {
                    last_point = next_point;
                    next_point = nbr[0];
                }
            }
            // Add the last point of the branch
            branch.push(next_point);

            let mut branch_temp: Vec<i32> = vec![0; self.medial_points.len()];
            for i in 0..branch.len() {
                branch_temp[branch[i]] = i as i32 +1;
            }
            branches.push(branch_temp);
        }

        // If the medial axis has at least a branch point
        else {
            // Each multiple point
            for i in 0..ind_multiple.len() {
                let nbr_trip = functions_vec::find_vec(&self.adjacency_matrix[ind_multiple[i]], true);
                // Each neighbour
                for j in 0..nbr_trip.len() {
                    let mut captured = false;
                    if !branches.is_empty() {
                        for b in 0..branches.len() {
                            captured = captured || (branches[b][nbr_trip[j]] != 0 && branches[b][ind_multiple[i]] != 0);
                        }
                    }

                    if !captured {
                        let mut branch: Vec<usize> = Vec::new();
                        // Starting point
                        branch.push(ind_multiple[i]);
                        let mut next_point = nbr_trip[j];
                        let mut last_point = ind_multiple[i];

                        // Looking for the neighbours of the neighbour...
                        while self.point_type[next_point] == 0 {
                            branch.push(next_point);
                            let nbr = functions_vec::find_vec(&self.adjacency_matrix[next_point], true);
                            if nbr[0] == last_point {
                                last_point = next_point;
                                next_point = nbr[1];
                            } else {
                                last_point = next_point;
                                next_point = nbr[0];
                            }
                        }
                        // Add the last point of the branch
                        branch.push(next_point);

                        let mut branch_temp: Vec<i32> = vec![0; self.medial_points.len()];
                        for i in 0..branch.len() {
                            branch_temp[branch[i]] = i as i32 +1;
                        }
                        branches.push(branch_temp);
                    }
                }
            }
        }

        // Adjacency matrix
        let mut adj_matrix: Vec<Vec<bool>> = vec![vec![false; branches.len()]; branches.len()];
        for i in 0..branches.len() {
            // Finding the endpoints of the branch
            let end_pt1 = functions_vec::find_vec(&branches[i], 1);
            let max_val = branches[i].iter().max().unwrap();
            let end_pt2 = functions_vec::find_vec(&branches[i], *max_val);

            // Find the branches containing the endpoints
            for j in 0..branches.len() {
                if i != j && (branches[j][end_pt1[0]] != 0 || branches[j][end_pt2[0]] != 0) {
                    adj_matrix[i][j] = true;
                    adj_matrix[j][i] = true;
                }
            }
        }
        self.branch_number = branches;
        self.branch_adjacency = adj_matrix;
        self.edf_branch = vec![0.0; self.branch_number.len()];
        self.wedf_branch = vec![0.0; self.branch_number.len()];

        Ok(())
    }


    fn extract_triangles(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let mut triangles: Vec<Vec<i32>> = Vec::new();
        for i in 0..self.medial_points.len() {
            let currbounds = &self.index_bndry_points[i];
            let mut edges: Vec<Vec<usize>> = Vec::new();
            let mut new_triangle: Vec<i32> = Vec::new();

            if currbounds.len() > 3 {
                for j in 0..currbounds.len() {
                    let (ind1, ind2) = functions_vec::find_vecvector2(&self.boundary_edges, currbounds[j]);
                    for k in 0..ind1.len() {
                        let next_ind;
                        if ind2[k] == 0 {
                            next_ind = self.boundary_edges[ind1[k]][1];
                        } else {
                            next_ind = self.boundary_edges[ind1[k]][0];
                        }
                        let next_bound = functions_vec::find_vec(&currbounds, next_ind);
                        if !next_bound.is_empty() {
                            edges.push(vec![j, next_bound[0]]);
                        }
                    }
                }
                // Delaunay triangulation
                let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f32>>::new();
                // Adding the vertices for the Delaunay triangulation
                let mut points:Vec<Point2<f32>> = Vec::new();
                for k in currbounds {
                    points.push(Point2::new(self.boundary_points[*k as usize][0], self.boundary_points[*k as usize][1]));
                }
                for k in 0..points.len() {
                    cdt.insert(points[k]);
                }
                // Adding the constraints
                for k in 0..edges.len() {
                    if !cdt.intersects_constraint(points[edges[k][0]], points[edges[k][1]]) {
                        cdt.add_constraint_edge(points[edges[k][0]], points[edges[k][1]])?;
                    }
                }
                // Extracting the faces
                let inner_faces = cdt.inner_faces();
                for face in inner_faces {
                    let mut new_triangle: Vec<i32> = Vec::new();
                    new_triangle.push(i as i32);
                    let vertices = face.vertices();
                    for vertex in vertices {
                        let point = vertex.position();
                        for k in currbounds {
                            if self.boundary_points[*k as usize][0] == point.x && self.boundary_points[*k as usize][1] == point.y {
                                new_triangle.push(*k);
                            }
                        }
                    }
                    triangles.push(new_triangle);
                }
            } else {
                new_triangle.push(i as i32);
                for j in 0..3 {
                    new_triangle.push(currbounds[j]);
                }
                triangles.push(new_triangle);
            }
        }
        self.triangles = triangles;
        Ok(())
    }


    // HIERARCHY DECOMPOSITION
    pub fn hierarchy_wedf(&mut self, mainclusts: usize, mainclustnum: usize, detailclustmin: usize, detailclustmax: usize) ->
    (Vec<i32>, Vec<f32>, Vec<f32>, i32) {
        // input is bma result containing branch fields. Output is labeling of
        // medial points according to k-means clustering on WEDF values. Clustering
        // happens in two stages. The first stage, with 3 clusters, takes the top
        // mainclustnum clusters as the main shape. The second stage automatically chooses a
        // number of clusters in detailclusterange based on gap analysis, then assigns all points in the
        // main shape from the first stage to the same cluster (cluster 1). All
        // clustering is determined from WEDF values at branch points and points
        // adjacent to branch points (for loops these "adjacent" points may be in 
        // the center of the loop), then all other points are assigned to a
        // cluster based on proximity to the centroids.

        // INPUTS :
        // mainclusts: the number of clusters you want to use in determining the main shape. (has been 3)
        // mainclustnum: the number of clusters from the resulting stage 1 clustering that you want to declare as the Main Shape (has been 2)
        // detailclustmin: the lowest number of clusters you want considered for the stage 2 detail clustering (have been using 8)
        // detailclustmax: the highest number of clusters you want considered for the stage 2 detail clustering (have been using 12)

        // Output
        let mut hierlabel: Vec<i32> = vec![-1; self.medial_points.len()];

        let mut nubinds = functions_vec::find_vec(&self.point_type, 3); // triple points
        let mut nbrs: Vec<usize> = Vec::new(); // neighbors of triple points that aren't on a loop
        let mut looppts: Vec<usize> = Vec::new(); // neighbors of triple points on a loop
        for i in 0..self.istrip_nbr.len() {
            if self.istrip_nbr[i] {
                if self.onloop[i] {
                    looppts.push(i);
                } else {
                    nbrs.push(i);
                    nubinds.push(i);
                }
            }
        }
        
        // tripptnbrs on loops are different since we are burning from halfway, not trip pt
        let mut foundbranch: Vec<usize> = Vec::new();
        let mut loopnbrs: Vec<usize> = Vec::new();

        for i in 0..looppts.len() {
            let mut branchno: Vec<usize> = Vec::new();
            for b in 0..self.branch_number.len() {
                if self.branch_number[b][looppts[i]] != 0 {
                    branchno.push(b);
                }
            }

            let found = functions_vec::find_vec(&foundbranch, branchno[0]);
            if found.is_empty() {
                foundbranch.push(branchno[0]);
                let (_max_val, ind_max) = functions_vec::max_f32(&self.wedf_array);
                let (_min_val, ind_min) = functions_vec::min_f32(&self.wedf_array);
                loopnbrs.push(ind_max[0]);
                loopnbrs.push(ind_min[0]);
            }
        }
        nubinds.append(&mut loopnbrs);

        // Extract the WEDF values
        let mut data_wedf: Vec<Vec<f32>> = Vec::new();
        let mut data_wedf_normalized: Vec<Vec<f32>> = Vec::new();
        let mut mean_wedf = 0.0;
        
        for i in 0..self.wedf_array.len() {
            data_wedf.push(vec![self.wedf_array[i]]);
            mean_wedf += self.wedf_array[i];
        }
        mean_wedf = mean_wedf/ self.wedf_array.len() as f32;
        
        for i in 0..self.wedf_array.len() {
            data_wedf_normalized.push(vec![self.wedf_array[i]/mean_wedf]);
        }

        // Normalize the wedf values
        let mut nub_wedf_normalized: Vec<Vec<f32>> = Vec::new();

        if nubinds.is_empty() {
            // Only the main part if one branch or one loop
            hierlabel = vec![1; self.medial_points.len()];

            self.clusters_wedf = hierlabel.clone();
            self.clusters_ind = nubinds.clone();
            let centroids_empty: Vec<f32> = Vec::new();
            
            return (hierlabel, centroids_empty.clone(), centroids_empty.clone(), 0)
        } else {
            for i in 0..nubinds.len() {
                nub_wedf_normalized.push(vec![self.wedf_array[nubinds[i]]/mean_wedf]);
            }
        }

        // Stage 1: find main body using mainclustnum of mainclusts clusters

        // Find the seeds with the function percentile
        let centroids = kmeans::init_percentiles(mainclusts, &nub_wedf_normalized);
        let clustering = kmeans::kmeans(mainclusts, &nub_wedf_normalized, 200, centroids);

        // Sorting the centroids
        let mut centroids1: Vec<f32> = Vec::new();
        for i in 0..clustering.centroids.len() {
            centroids1.push(clustering.centroids[i].0[0] as f32);
        }
        centroids1.sort_by(|a, b| a.partial_cmp(b).unwrap());
        centroids1.reverse();

        // Finding the cluster of each data
        for i in 0..self.wedf_array.len() {
            let mut cluster = vec![0.0, (data_wedf_normalized[i][0]-centroids1[0]).abs()];
            for j in 0..centroids1.len() {
                if (data_wedf_normalized[i][0]-centroids1[j]).abs()<cluster[1] {
                    cluster[0] = j as f32;
                    cluster[1] = (data_wedf_normalized[i][0]-centroids1[j]).abs();
                }
            }
            hierlabel[i] = cluster[0] as i32 +1;
        }
        let before_spread_hierlabel = hierlabel.clone();

        // Creating the main cluster
        let mut ind_main: Vec<usize> = Vec::new();
        for i in 0..mainclustnum {
            let ind_clust_i = functions_vec::find_vec(&hierlabel, i as i32 +1);
            ind_main.extend(ind_clust_i.iter());
        }

        // LOOPS
        // Adding the whole loop if there is one of the main points inside a loop
        /* *//*
        let mut ind_loop: Vec<usize> = Vec::new();
        let mut comp_to_add: Vec<i32> = Vec::new();
        for i in &ind_main {
            if self.onloop[*i] {
                let border_pts = self.index_bndry_points[*i].clone();
                for br in &border_pts {
                    let comp = self.component[*br as usize];
                    if !self.isouter[*br as usize] && !comp_to_add.contains(&comp) {
                        comp_to_add.push(comp);
                    }
                }
            }
        }
        // Finding the indices of the loop
        for i in 0..self.medial_points.len() {
            if self.onloop[i] {
                let border_pts = self.index_bndry_points[i].clone();
                for br in &border_pts {
                    let comp = self.component[*br as usize];
                    if comp_to_add.contains(&comp) && !ind_loop.contains(&i) {
                        ind_loop.push(i);
                    }
                }
            }
        }
        // Add the indices of the loop inside the main points
        for i in &ind_loop {
            if !ind_main.contains(i) {
                ind_main.push(*i);
            }
        } /* */*/
        // END OF LOOPS


        for i in &ind_main {
            hierlabel[*i] = 1;
        }
        
        // Stage 2: use more clusters to identify hierarchy of extremities

        // Find the indices for the rest of the shape
        let mut ind_details: Vec<usize> = Vec::new();
        for i in mainclustnum..mainclusts {
            let ind_clust_i = functions_vec::find_vec(&hierlabel, i as i32 +1);
            ind_details.extend(ind_clust_i.iter());
        }

        // If there is other parts than the main one
        let mut centroids2: Vec<f32> = Vec::new();
        let mut optimal_k: usize = 0;
        if !ind_details.is_empty() {
            let mut nubinds2: Vec<usize> = Vec::new();
            for i in &ind_details {
                let ind_clust_i = functions_vec::find_vec(&nubinds, *i);
                nubinds2.extend(ind_clust_i.iter());
            }

            // Finding the WEDF values
            let mut nub_wedf_normalized: Vec<Vec<f32>> = Vec::new();
            if nubinds2.is_empty() {
                for i in &ind_details {
                    nub_wedf_normalized.push(vec![self.wedf_array[*i]/mean_wedf]);
                }
            } else {
                for i in &nubinds2 {
                    nub_wedf_normalized.push(vec![self.wedf_array[*i]/mean_wedf]);
                }
            }
            
            // Finding the optimal k for kmeans
            if detailclustmin == detailclustmax {
                optimal_k = detailclustmin;
            } else {
                if detailclustmax <= nub_wedf_normalized.len() {
                    optimal_k = gap_evaluation::gap_eval(&nub_wedf_normalized, detailclustmin, detailclustmax);
                } else {
                    optimal_k = gap_evaluation::gap_eval(&data_wedf_normalized, detailclustmin, detailclustmax);
                }
            }

            // Find the seeds with the function percentile
            let centroids = kmeans::init_percentiles(optimal_k, &nub_wedf_normalized);
            let clustering = kmeans::kmeans(optimal_k, &nub_wedf_normalized, 200, centroids);

            // Sorting the centroids
            for i in 0..clustering.centroids.len() {
                centroids2.push(clustering.centroids[i].0[0] as f32);
            }
            centroids2.sort_by(|a, b| a.partial_cmp(b).unwrap());
            centroids2.reverse();

            // Finding the cluster of each data
            for i in ind_details {
                let mut cluster = vec![0.0, (data_wedf_normalized[i][0]-centroids2[0]).abs()];
                for j in 0..centroids2.len() {
                    if (data_wedf_normalized[i][0]-centroids2[j]).abs()<cluster[1] {
                        cluster[0] = j as f32;
                        cluster[1] = (data_wedf_normalized[i][0]-centroids2[j]).abs();
                    }
                }
                hierlabel[i] = cluster[0] as i32 + 2;
            }
        }

        // OUTPUTS
        // hierlabel: a number of medial points x 1 array with integers labeling the cluster to which each medial point belongs after clustering. the lower the label value, the larger the WEDF value.
        // centroids1: the centroids of the clusters from the stage 1 clustering
        // centroids2: the centroids of the clusters from the stage 2 clustering
        // optimal_k: the number of clusters selected by the stage 2 clustering
        // nubinds
        self.clusters_wedf = hierlabel.clone();
        self.clusters_ind = nubinds.clone();

        (before_spread_hierlabel, centroids1, centroids2, optimal_k as i32)
    }


    pub fn hierarchy_wedf_v2(&mut self, mainclusts: usize, mainclustnum: usize, detailclustmin: usize, detailclustmax: usize) ->
    (Vec<i32>, Vec<f32>, Vec<f32>, i32) {
        // input is bma result containing branch fields. Output is labeling of
        // medial points according to k-means clustering on WEDF values. Clustering
        // happens in two stages. The first stage, with 3 clusters, takes the top
        // mainclustnum clusters as the main shape. The second stage automatically chooses a
        // number of clusters in detailclusterange based on gap analysis, then assigns all points in the
        // main shape from the first stage to the same cluster (cluster 1). All
        // clustering is determined from WEDF values at branch points and points
        // adjacent to branch points (for loops these "adjacent" points may be in 
        // the center of the loop), then all other points are assigned to a
        // cluster based on proximity to the centroids.

        // INPUTS :
        // mainclusts: the number of clusters you want to use in determining the main shape. (has been 3)
        // mainclustnum: the number of clusters from the resulting stage 1 clustering that you want to declare as the Main Shape (has been 2)
        // detailclustmin: the lowest number of clusters you want considered for the stage 2 detail clustering (have been using 8)
        // detailclustmax: the highest number of clusters you want considered for the stage 2 detail clustering (have been using 12)

        // Output
        let mut hierlabel: Vec<i32> = vec![-1; self.medial_points.len()];

        let mut nubinds = functions_vec::find_vec(&self.point_type, 3); // triple points
        for i in 0..self.istrip_nbr.len() {
            if self.istrip_nbr[i] {
                nubinds.push(i);
            }
        }

        // Extract the WEDF values
        let mut data_wedf: Vec<Vec<f32>> = Vec::new();
        let mut data_wedf_normalized: Vec<Vec<f32>> = Vec::new();
        let mut mean_wedf = 0.0;
        
        for i in 0..self.wedf_array.len() {
            data_wedf.push(vec![self.wedf_array[i]]);
            mean_wedf += self.wedf_array[i];
        }
        mean_wedf = mean_wedf/ self.wedf_array.len() as f32;
        
        for i in 0..self.wedf_array.len() {
            data_wedf_normalized.push(vec![self.wedf_array[i]/mean_wedf]);
        }

        // Normalize the wedf values
        let mut nub_wedf_normalized: Vec<Vec<f32>> = Vec::new();

        if nubinds.is_empty() {
            // Only the main part if one branch or one loop
            hierlabel = vec![1; self.medial_points.len()];

            self.clusters_wedf = hierlabel.clone();
            self.clusters_ind = nubinds.clone();
            let centroids_empty: Vec<f32> = Vec::new();
    
            return (hierlabel, centroids_empty.clone(), centroids_empty.clone(), 0)
        } else {
            for i in 0..nubinds.len() {
                nub_wedf_normalized.push(vec![self.wedf_array[nubinds[i]]/mean_wedf]);
            }
        }

        // Stage 1: find main body using mainclustnum of mainclusts clusters

        // Find the seeds with the function percentile
        let centroids = kmeans::init_percentiles(mainclusts, &nub_wedf_normalized);
        let clustering = kmeans::kmeans(mainclusts, &nub_wedf_normalized, 200, centroids);

        // Sorting the centroids
        let mut centroids1: Vec<f32> = Vec::new();
        for i in 0..clustering.centroids.len() {
            centroids1.push(clustering.centroids[i].0[0] as f32);
        }
        centroids1.sort_by(|a, b| a.partial_cmp(b).unwrap());
        centroids1.reverse();

        // Finding the cluster of each data
        for i in 0..nubinds.len() {
            // Find the closest centroid
            let mut cluster = vec![0.0, (data_wedf_normalized[i][0]-centroids1[0]).abs()];
            for c in 0..centroids1.len() {
                if (nub_wedf_normalized[i][0]-centroids1[c]).abs()<cluster[1] {
                    cluster[0] = c as f32;
                    cluster[1] = (nub_wedf_normalized[i][0]-centroids1[c]).abs();
                }
            }
            hierlabel[nubinds[i]] = cluster[0] as i32 +1;
        }
        
        // Find all the points of each branch
        for branch in 0..self.branch_number.len() {
            let mut ind_br: Vec<usize> = Vec::new();
            let mut wedf_nubind_br: Vec<f32> = Vec::new();
            let mut nubind_br: Vec<usize> = Vec::new();
            for j in 0..self.branch_number[branch].len() {
                if self.branch_number[branch][j] > 0 {
                    ind_br.push(j);
                    if nubinds.contains(&j) {
                        nubind_br.push(j);
                        wedf_nubind_br.push(self.wedf_array[j]);
                    }
                }
            }
            let (max_wedf, ind_max) = functions_vec::max_f32(&wedf_nubind_br);
            for j in &ind_br {
                hierlabel[*j] = hierlabel[nubind_br[ind_max[0]]] as i32;
            }
        }
        let before_spread_hierlabel = hierlabel.clone();

        // Creating the main cluster
        let mut ind_main: Vec<usize> = Vec::new();
        for i in 0..mainclustnum {
            let ind_clust_i = functions_vec::find_vec(&hierlabel, i as i32 +1);
            ind_main.extend(ind_clust_i.iter());
        }

        for i in &ind_main {
            hierlabel[*i] = 1;
        }
        
        // Stage 2: use more clusters to identify hierarchy of extremities

        // Find the indices for the rest of the shape
        let mut ind_details: Vec<usize> = Vec::new();
        for i in mainclustnum..mainclusts {
            let ind_clust_i = functions_vec::find_vec(&hierlabel, i as i32 +1);
            ind_details.extend(ind_clust_i.iter());
        }

        // If there is other parts than the main one
        let mut centroids2: Vec<f32> = Vec::new();
        let mut optimal_k: usize = 0;
        if !ind_details.is_empty() {
            let mut nubinds2: Vec<usize> = Vec::new();
            for i in &ind_details {
                let ind_clust_i = functions_vec::find_vec(&nubinds, *i);
                nubinds2.extend(ind_clust_i.iter());
            }

            // Finding the WEDF values
            let mut nub_wedf_normalized: Vec<Vec<f32>> = Vec::new();
            if nubinds2.is_empty() {
                for i in &ind_details {
                    nub_wedf_normalized.push(vec![self.wedf_array[*i]/mean_wedf]);
                }
            } else {
                for i in &nubinds2 {
                    nub_wedf_normalized.push(vec![self.wedf_array[*i]/mean_wedf]);
                }
            }
            
            // Finding the optimal k for kmeans
            if detailclustmin == detailclustmax {
                optimal_k = detailclustmin;
            } else {
                if detailclustmax <= nub_wedf_normalized.len() {
                    optimal_k = gap_evaluation::gap_eval(&nub_wedf_normalized, detailclustmin, detailclustmax);
                } else {
                    optimal_k = gap_evaluation::gap_eval(&data_wedf_normalized, detailclustmin, detailclustmax);
                }
            }

            // Find the seeds with the function percentile
            let centroids = kmeans::init_percentiles(optimal_k, &nub_wedf_normalized);
            let clustering = kmeans::kmeans(optimal_k, &nub_wedf_normalized, 200, centroids);

            // Sorting the centroids
            for i in 0..clustering.centroids.len() {
                centroids2.push(clustering.centroids[i].0[0] as f32);
            }
            centroids2.sort_by(|a, b| a.partial_cmp(b).unwrap());
            centroids2.reverse();

            // Finding the cluster of each data
            for i in ind_details {
                let mut cluster = vec![0.0, (data_wedf_normalized[i][0]-centroids2[0]).abs()];
                for j in 0..centroids2.len() {
                    if (data_wedf_normalized[i][0]-centroids2[j]).abs()<cluster[1] {
                        cluster[0] = j as f32;
                        cluster[1] = (data_wedf_normalized[i][0]-centroids2[j]).abs();
                    }
                }
                hierlabel[i] = cluster[0] as i32 + 2;
            }
        }

        // OUTPUTS
        // hierlabel: a number of medial points x 1 array with integers labeling the cluster to which each medial point belongs after clustering. the lower the label value, the larger the WEDF value.
        // centroids1: the centroids of the clusters from the stage 1 clustering
        // centroids2: the centroids of the clusters from the stage 2 clustering
        // optimal_k: the number of clusters selected by the stage 2 clustering
        // nubinds
        self.clusters_wedf = hierlabel.clone();
        self.clusters_ind = nubinds.clone();

        (before_spread_hierlabel, centroids1, centroids2, optimal_k as i32)
    }



    // TRUNKS
    fn find_ma_pts_on_loop(&self, comp: usize) -> Vec<usize> {
        let mut pts_on_loop: Vec<usize> = Vec::new();
        for i in 0..self.medial_points.len() {
            if self.onloop[i] {
                let bnd_pts = self.index_bndry_points[i].clone();
                let mut comps: Vec<usize> = Vec::new();
                for j in &bnd_pts {
                    if !comps.contains(&(self.component[*j as usize] as usize)) {
                        comps.push(self.component[*j as usize] as usize);
                    }
                }
                if comps.contains(&comp) {
                    pts_on_loop.push(i);
                }
            }
        }
        pts_on_loop
    }

    fn find_next_cut(&self, level: Vec<i32>, cut_value: Vec<i32>) -> (Vec<i32>, Vec<i32>) {
        // Finds the next level in the hierarchy and labels points as such
        let &previous_cut = level.iter().max().unwrap();
        let mut next_level = level.clone();
        let mut new_cut_value: Vec<i32> = cut_value.clone();

        let mut potential_cut_points: Vec<usize> = Vec::new();
        for i in 0..level.len() {
            let ind = functions_vec::find_vec(&self.clusters_ind, i);
            if level[i] == previous_cut && ind.is_empty() {
                potential_cut_points.push(i);
            }
        }
        let mut min_hierlabel = self.clusters_wedf[potential_cut_points[0]];
        for i in &potential_cut_points {
            if self.clusters_wedf[*i] < min_hierlabel {
                min_hierlabel = self.clusters_wedf[*i];
            }
        }
        new_cut_value.push(min_hierlabel);
        for i in 0..next_level.len() {
            if level[i] == previous_cut && self.clusters_wedf[i] >= min_hierlabel+1 {
                next_level[i] = previous_cut + 1;
            }
        }
        (next_level, new_cut_value)
    }

    fn graph_components(adjacency_matrix: &Vec<Vec<bool>>) -> (usize, Vec<usize>, Vec<Vec<usize>>) {
        // Finds the branches that belong together in each hierarchy level
        let mut members: Vec<Vec<usize>> = Vec::new();
        // Have we visited a particular node yet?
        let mut is_discovered = vec![false; adjacency_matrix.len()];
        // check every node
        for i in 0..adjacency_matrix.len() {
            if !is_discovered[i] {
                // started a new group so add it to members
                members.push(vec![i]);
                // account for discovering n
                is_discovered[i] = true;
                // set the ptr to 0
                let mut ptr = 0;
                while ptr < members[members.len()-1].len() {
                    // find neighbors
                    let nbrs = functions_vec::find_vec(&adjacency_matrix[members[members.len()-1][ptr]], true);
                    let nb_members = members.len();
                    for j in 0..nbrs.len() {
                        // here are the neighbors that are undiscovered
                        if !is_discovered[nbrs[j]] {
                            // we can now mark them as discovered
                            is_discovered[nbrs[j]] = true;
                            // add them to member list
                            members[nb_members-1].push(nbrs[j]);
                        }
                    }
                    // increment ptr so we check the next member of this component
                    ptr += 1;
                }
            }
        }
        // number of components
        let n_components = members.len();
        // compute sizes of components
        let mut sizes: Vec<usize> = Vec::new();
        for i in 0..n_components-1 {
            sizes.push(members[i].len());
        }
        (n_components, sizes, members)
    }

    fn find_next_trunk_pts_onloop(&self, trunk_pts: &Vec<usize>, branchpath: &Vec<usize>, lastpoint: usize, compmembers: &Vec<usize>, branchnums: &Vec<usize>, checkbranch: &mut Vec<bool>) -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<usize>, Vec<bool>) {
        // Find the loops of the last point
        let mut comps: Vec<usize> = Vec::new();
        let bnd_pts = &self.index_bndry_points[lastpoint];
        for i in bnd_pts {
            if !self.isouter[*i as usize] && !comps.contains(&(self.component[*i as usize] as usize)) {
                comps.push(self.component[*i as usize] as usize);
            }
        }
        
        let mut new_trunk_pts: Vec<Vec<usize>> = Vec::new();
        let mut new_branchpath: Vec<Vec<usize>> = Vec::new();
        let mut new_lastpoint: Vec<usize> = Vec::new();
        let mut ending_points: Vec<bool> = Vec::new();
        let mut num_trunk = 0;

        for comp in &comps {
            let pts_onloop = self.find_ma_pts_on_loop(*comp);
            let pts_onloop_incomp = functions_vec::intersect(&compmembers, &pts_onloop);
            let mut wedf_values: Vec<f32> = Vec::new();
            for i in &pts_onloop_incomp {
                wedf_values.push(self.wedf_array[*i]);
            }
            let (max, ind_max) = functions_vec::max_f32(&wedf_values);
            let (min, ind_min) = functions_vec::min_f32(&wedf_values);
            let ind_max = pts_onloop_incomp[ind_max[0]];
            let ind_min = pts_onloop_incomp[ind_min[0]];

            // Find neighbors outside the loop
            let mut nbrs: Vec<usize> = Vec::new();
            for i in &pts_onloop_incomp {
                let nbrs_ind = functions_vec::find_vec(&self.adjacency_matrix[*i], true);   
                for j in &nbrs_ind {
                    if compmembers.contains(j) && min > self.wedf_array[*j] && !self.onloop[*j] {
                        nbrs.push(*j);
                    }
                }
            }

            if pts_onloop.len() == pts_onloop_incomp.len() {// All the points of the loop are inside the component of the shape
                for i in 0..2 {
                    new_trunk_pts.push(trunk_pts.clone());
                    let last_index = new_trunk_pts[num_trunk].len()-1;
                    new_trunk_pts[num_trunk].remove(last_index); // Erase the point of the loop that was added since we changed it
                    new_trunk_pts[num_trunk].push(ind_max);
                    new_branchpath.push(branchpath.clone());
                    
                    // Starting point of the loop
                    let mut nbrs_max: Vec<usize> = Vec::new();
                    let nbrs_max_ind = functions_vec::find_vec(&self.adjacency_matrix[lastpoint], true);   
                    for j in &nbrs_max_ind {
                        if pts_onloop_incomp.contains(j) && self.wedf_array[ind_max] > self.wedf_array[*j] {
                            nbrs_max.push(*j);
                        }
                    }
                    new_trunk_pts[num_trunk].push(nbrs_max_ind[i]);

                    // Find the branch
                    let mut branch = 0;
                    for br in 0..self.branch_number.len() {
                        if self.branch_number[br][ind_max] > 0 && self.branch_number[br][nbrs_max_ind[i]] > 0 {
                            branch = br;
                        }
                    }
                    if !new_branchpath[num_trunk].contains(&branch) {
                        new_branchpath[num_trunk].push(branch);
                    }

                    while new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1] != ind_min {
                        let nbrs = functions_vec::find_vec(&self.adjacency_matrix[new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1]], true);
                        let mut next_pt = ind_min;  
                        for j in &nbrs {
                            if pts_onloop_incomp.contains(j) && self.wedf_array[new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1]] > self.wedf_array[*j] {
                                next_pt = *j;
                            }
                        }

                        // Find the branch
                        let mut branch = 0;
                        for br in 0..self.branch_number.len() {
                            if self.branch_number[br][new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1]] > 0 && self.branch_number[br][next_pt] > 0 {
                                branch = br;
                            }
                        }

                        new_trunk_pts[num_trunk].push(next_pt);
                        if !new_branchpath[num_trunk].contains(&branch) {
                            new_branchpath[num_trunk].push(branch);
                        }
                    }

                    if nbrs.is_empty() {
                        new_lastpoint.push(new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1]);
                        ending_points.push(true);
                        num_trunk += 1;
                    } else {
                        new_lastpoint.push(nbrs[0]);
                        ending_points.push(false);
                        num_trunk += 1;
                        for j in 1..nbrs.len() {
                            // PEUT ETRE MODIFIER ET METTRE TOUS LES POINTS DE LA BRANCHE CONCERNEE
                            new_trunk_pts.push(new_trunk_pts[num_trunk-1].clone());
                            new_branchpath.push(new_branchpath[num_trunk-1].clone());
                            new_lastpoint.push(nbrs[j]);
                            ending_points.push(false);
                            num_trunk += 1;
                        }
                    }
                }
            } else { // Only a part of the loop is inside the part
                // Find the extremities of the part of the loop
                let mut starting_pts_onloop: Vec<usize> = Vec::new();
                let pts_onloop = self.find_ma_pts_on_loop(*comp);
                let pts_onloop_incomp = functions_vec::intersect(&compmembers, &pts_onloop);
                for i in &pts_onloop_incomp {
                    let mut num_nbrs = 0;
                    for j in &pts_onloop_incomp {
                        if self.adjacency_matrix[*i][*j] {
                            num_nbrs += 1;
                        }
                    }
                    if num_nbrs == 1 {
                        starting_pts_onloop.push(*i);
                    }
                }
                assert!(starting_pts_onloop.len() == 2);
                
                for i in 0..2 {
                    new_trunk_pts.push(trunk_pts.clone());
                    let last_index = new_trunk_pts[num_trunk].len()-1;
                    new_trunk_pts[num_trunk].remove(last_index); // Erase the point of the loop that was added since we changed it
                    new_trunk_pts[num_trunk].push(starting_pts_onloop[i]);
                    new_branchpath.push(branchpath.clone());

                    while new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1] != ind_min {
                        let nbrs = functions_vec::find_vec(&self.adjacency_matrix[new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1]], true);
                        let mut next_pt = ind_min;  
                        for j in &nbrs {
                            if pts_onloop_incomp.contains(j) && self.wedf_array[new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1]] > self.wedf_array[*j] {
                                next_pt = *j;
                            }
                        }

                        // Find the branch
                        let mut branch = 0;
                        for br in 0..self.branch_number.len() {
                            if self.branch_number[br][new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1]] > 0 && self.branch_number[br][next_pt] > 0 {
                                branch = br;
                            }
                        }

                        new_trunk_pts[num_trunk].push(next_pt);
                        if !new_branchpath[num_trunk].contains(&branch) {
                            new_branchpath[num_trunk].push(branch);
                        }
                    }

                    if nbrs.is_empty() {
                        new_lastpoint.push(new_trunk_pts[num_trunk][new_trunk_pts[num_trunk].len()-1]);
                        ending_points.push(true);
                        num_trunk += 1;
                    } else {
                        new_lastpoint.push(nbrs[0]);
                        ending_points.push(false);
                        num_trunk += 1;
                        for j in 1..nbrs.len() {
                            // PEUT ETRE MODIFIER ET METTRE TOUS LES POINTS DE LA BRANCHE CONCERNEE
                            new_trunk_pts.push(new_trunk_pts[num_trunk-1].clone());
                            new_branchpath.push(new_branchpath[num_trunk-1].clone());
                            new_lastpoint.push(nbrs[j]);
                            ending_points.push(false);
                            num_trunk += 1;
                        }
                    }
                }
            }
        }

        (new_trunk_pts, new_branchpath, new_lastpoint, ending_points)
    }

    fn find_next_trunk_pts(&self, trunk_pts: &Vec<usize>, branchpath: &Vec<usize>, lastpoint: usize, compmembers: &Vec<usize>, branchnums: &Vec<usize>, checkbranch: &mut Vec<bool>) -> (Vec<Vec<usize>>, Vec<Vec<usize>>, Vec<usize>, Vec<bool>) {
        // When on a loop
        if self.onloop[lastpoint] {
            return self.find_next_trunk_pts_onloop(trunk_pts, branchpath, lastpoint, compmembers, branchnums, checkbranch);
        }
        
        // Following points of the last point on the trunk
        let mut nbrs: Vec<usize> = Vec::new();
        let nbrs_ind = functions_vec::find_vec(&self.adjacency_matrix[lastpoint], true);   
        for i in &nbrs_ind {
            if compmembers.contains(i) && self.wedf_array[lastpoint] > self.wedf_array[*i] {
                nbrs.push(*i);
            }
        }

        let mut new_trunk_pts: Vec<Vec<usize>> = Vec::new();
        let mut new_branchpath: Vec<Vec<usize>> = Vec::new();
        let mut new_lastpoint: Vec<usize> = Vec::new();
        let mut ending_points: Vec<bool> = Vec::new();
        let mut num_nbr = 0;

        let mut all_branches_checked = true;
        for i in 0..checkbranch.len() {
            all_branches_checked = all_branches_checked && checkbranch[i];
        }

        if nbrs.is_empty() || all_branches_checked {
            new_trunk_pts.push(trunk_pts.clone());
            new_branchpath.push(branchpath.clone());
            new_lastpoint.push(lastpoint);
            ending_points.push(true);
        }

        for nbr in &nbrs {
            num_nbr += 1;
            // Find the branch between the last point and its neighbors
            let mut branch = 0;
            for br in 0..self.branch_number.len() {
                if self.branch_number[br][lastpoint] > 0 && self.branch_number[br][*nbr] > 0 {
                    branch = br;
                }
            }
            new_branchpath.push(branchpath.clone());
            new_branchpath[num_nbr-1].push(branch);
            let _ind = functions_vec::find_vec(&branchnums, branch);
            checkbranch[_ind[0]] = true;

            // Find the trunk points
            new_trunk_pts.push(trunk_pts.clone());
            let mut ptinds: Vec<usize> = Vec::new();
            let mut inds: Vec<i32> = Vec::new();
            for i in 0..self.branch_number[branch].len() {
                if self.branch_number[branch][i] > 0 {
                    inds.push(self.branch_number[branch][i]);
                    ptinds.push(i);
                }
            }
            let orderinds = functions_vec::argsort(&inds);
            let mut temp_ptinds: Vec<usize> = Vec::new();
            for i in &orderinds {
                temp_ptinds.push(ptinds[*i]);
            }
            ptinds = temp_ptinds.clone();

            if self.branch_number[branch][lastpoint] < self.branch_number[branch][*nbr] {
                let ind = functions_vec::find_vec(&ptinds, lastpoint);
                for i in ind[0]+1..ptinds.len() {
                    new_trunk_pts[num_nbr-1].push(ptinds[i]);
                }
                new_lastpoint.push(ptinds[ptinds.len()-1]);
            } else {
                let ind = functions_vec::find_vec(&ptinds, lastpoint);
                for i in 0..ind[0] {
                    new_trunk_pts[num_nbr-1].push(ptinds[ind[0]-1-i]);
                }
                new_lastpoint.push(ptinds[0]);
            }

            // Find if we are on the last point of the trunk or not
            ending_points.push(false);
            if self.point_type[*nbr] == 1 {
                ending_points[num_nbr-1] = true;
            }
        }

        (new_trunk_pts, new_branchpath, new_lastpoint, ending_points)
    }

    fn collect_detail_structures(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let &max_level = self.hierarchy_level.iter().max().unwrap();
        let mut level_num = max_level;
        let mut trunknum = 0;
        let mut trunkmembers: Vec<Vec<bool>> = Vec::new();

        while level_num >= 2 {
            let mut ind_level: Vec<usize> = Vec::new();
            for e in level_num..(max_level+1) {
                let mut temp_ind = functions_vec::find_vec(&self.hierarchy_level, e);
                ind_level.append(&mut temp_ind);
            }
            let mut adj_ind_level: Vec<Vec<bool>> = Vec::new();
            for i in &ind_level {
                let mut col: Vec<bool> = Vec::new();
                for j in &ind_level {
                    col.push(self.adjacency_matrix[*i][*j]);
                }
                adj_ind_level.push(col);
            }
            let (n_components, _sizes, members) = Self::graph_components(&adj_ind_level);
            for comp_num in 0..n_components {
                let mut compmembers: Vec<usize> = Vec::new();
                for i in &members[comp_num] {
                    compmembers.push(ind_level[*i]);
                }

                // Get the branches of the component
                // EST CE QU'ON DOIT SUPPRIMER LES BRANCHES QUI SONT AUSSI INCLUES DANS DES NIVEAUX INFERIEURS ???
                let mut branchnums: Vec<usize> = Vec::new();
                let mut checkbranch: Vec<bool> = Vec::new();
                for i in 0..self.branch_number.len() {
                    for j in &compmembers {
                        if self.branch_number[i][*j] > 0 && !branchnums.contains(&i){
                            branchnums.push(i);
                            checkbranch.push(false);
                        }
                    }
                }
                
                let mut parentchild: Vec<Vec<usize>> = Vec::new();
                // Find the relation parent child between branches
                for branchno in 0..branchnums.len() {
                    let mut parents: Vec<usize> = Vec::new();
                    for br_comp in &branchnums {
                        if self.branch_adjacency[branchnums[branchno]][*br_comp] {
                            parents.push(*br_comp);
                        }
                    }
                    if !parents.is_empty() {
                        for i in &parents {
                            if self.wedf_branch[*i] > self.wedf_branch[branchnums[branchno]] && !parentchild.contains(&vec![*i, branchnums[branchno]]) {
                                parentchild.push(vec![*i, branchnums[branchno]]);
                            } else if self.wedf_branch[*i] < self.wedf_branch[branchnums[branchno]] && !parentchild.contains(&vec![branchnums[branchno], *i]) {
                                parentchild.push(vec![branchnums[branchno], *i]);
                            }
                        }
                    }
                }
                
                // find the max wedf of the component members
                let mut max_wedf_compmembers = self.wedf_array[compmembers[0]];
                let mut wedf_compmembers: Vec<f32> = Vec::new();
                for i in &compmembers {
                    wedf_compmembers.push(self.wedf_array[*i]);
                    if self.wedf_array[*i] > max_wedf_compmembers {
                        max_wedf_compmembers = self.wedf_array[*i];
                    }
                }

                // index of the first point of the trunks of this component
                let mut startptind = functions_vec::find_vec(&wedf_compmembers, max_wedf_compmembers);
                if startptind.len() > 1 {
                    let mut edf_startptind: Vec<f32> = Vec::new();
                    for i in 0..startptind.len() {
                        edf_startptind.push(self.edf_array[compmembers[startptind[i]]]);
                    }
                    let (_max, _ind) = functions_vec::max_f32(&edf_startptind);
                    startptind = vec![startptind[_ind[0]]];
                }
                let startptind = compmembers[startptind[0]];

                // find the first branch of the trunks of this component
                let mut topbranch: Vec<usize> = Vec::new();
                for br in 0..self.branch_number.len() {
                    if self.branch_number[br][startptind] > 0 {
                        if !topbranch.contains(&br) && branchnums.contains(&br) {
                            topbranch.push(br);
                        }
                    }
                }
                if topbranch.len() > 1 {
                    let mut maxval: Vec<f32> = Vec::new();
                    for br in &topbranch {
                        maxval.push(self.wedf_branch[*br]);
                    }
                    let (_maxval, topind) = functions_vec::max_f32(&maxval);
                    topbranch = vec![topbranch[topind[0]]];
                }
                let topbranch = topbranch[0];

                if self.hierarchy_level[startptind] == level_num {
                    let mut stop = true;

                    let mut all_trunk_pts: Vec<Vec<usize>> = Vec::new();
                    let mut all_branchpath: Vec<Vec<usize>> = Vec::new();
                    let mut all_lastpoint: Vec<usize> = Vec::new();
                    let mut all_endingpoint: Vec<bool> = Vec::new();

                    let mut final_trunk_pts: Vec<Vec<usize>> = Vec::new();
                    let mut final_branchpath: Vec<Vec<usize>> = Vec::new();
                    
                    let (trunk_pts, branchpath, lastpoint, endingpoint) = self.find_next_trunk_pts(&vec![startptind], &Vec::new(), startptind, &compmembers, &branchnums, &mut checkbranch);
                    
                    for i in 0..trunk_pts.len() {
                        if !endingpoint[i] {
                            all_trunk_pts.push(trunk_pts[i].clone());
                            all_branchpath.push(branchpath[i].clone());
                            all_lastpoint.push(lastpoint[i]);
                            all_endingpoint.push(endingpoint[i]);
                        } else {
                            final_trunk_pts.push(trunk_pts[i].clone());
                            final_branchpath.push(branchpath[i].clone());
                        }
                    }

                    for i in &all_endingpoint {
                        stop = stop && *i;
                    }
                    
                    while !stop {
                        let mut temp_all_trunk_pts: Vec<Vec<usize>> = Vec::new();
                        let mut temp_all_branchpath: Vec<Vec<usize>> = Vec::new();
                        let mut temp_all_lastpoint: Vec<usize> = Vec::new();
                        let mut temp_all_endingpoint: Vec<bool> = Vec::new();
                    
                        for i in 0..all_trunk_pts.len() {
                            let (trunk_pts, branchpath, lastpoint, endingpoint) = self.find_next_trunk_pts(&all_trunk_pts[i], &all_branchpath[i], all_lastpoint[i], &compmembers, &branchnums, &mut checkbranch);
                            for j in 0..trunk_pts.len() {
                                if !endingpoint[j] {
                                    temp_all_trunk_pts.push(trunk_pts[j].clone());
                                    temp_all_branchpath.push(branchpath[j].clone());
                                    temp_all_lastpoint.push(lastpoint[j]);
                                    temp_all_endingpoint.push(endingpoint[j]);
                                } else {
                                    final_trunk_pts.push(trunk_pts[j].clone());
                                    final_branchpath.push(branchpath[j].clone());
                                }
                            }
                        }
                        all_trunk_pts = temp_all_trunk_pts.clone();
                        all_branchpath = temp_all_branchpath.clone();
                        all_lastpoint = temp_all_lastpoint.clone();
                        all_endingpoint = temp_all_endingpoint.clone();

                        stop = true;
                        for i in &all_endingpoint {
                            stop = stop && *i;
                        }
                    }

                    // store all of this stuff for each trunk
                    for i in 0..final_trunk_pts.len() {
                        let new_trunk = Trunk{
                            level: level_num,
                            component: comp_num,
                            indices: final_trunk_pts[i].clone(),
                            branchpath: final_branchpath[i].clone(),
                            parent_child_array: parentchild.clone()
                        };
                        self.trunk.push(new_trunk);
                        // Take care of trunk members too
                        trunkmembers.push(vec![false; self.medial_points.len()]);
                        for i in &final_trunk_pts[i] {
                            trunkmembers[trunknum][*i] = true;
                        }
                        trunknum += 1;
                    }
                }
            }
            level_num = level_num-1;
        }
        self.trunk_members = trunkmembers.clone();
        Ok(())
    }

    pub fn shape_details_hierarchy(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        // initialize hierarchy level labels with main body = 1, the rest of the points = 2.
        let possible_values: Vec<i32> = functions_vec::unique_list(&self.clusters_wedf);
        let mut cut_value: Vec<i32> = Vec::new();
        if possible_values.len() > 1 {
            cut_value.push(possible_values[0]);
            cut_value.push(possible_values[1]);
        } else {
            cut_value = possible_values;
        }
        let mut level: Vec<i32> = vec![2; self.clusters_wedf.len()];
        let ind_main = functions_vec::find_vec(&self.clusters_wedf, 1);
        for i in ind_main {
            level[i] = 1;
        }

        // compute other levels
        let mut levelpts: Vec<usize> = Vec::new();
        for ind in &self.clusters_ind {
            if self.point_type[*ind] != 3 {
                levelpts.push(*ind);
            }
        }
        let mut finished = true;
        if !levelpts.is_empty() {
            let &max_cut_value = cut_value.iter().max().unwrap();
            let mut max_hierlabel = self.clusters_wedf[levelpts[0]];
            for i in &levelpts {
                if max_hierlabel < self.clusters_wedf[*i] {
                    max_hierlabel = self.clusters_wedf[*i];
                }
            }
            finished = max_cut_value == max_hierlabel;
        }

        while !finished {
            (level, cut_value) = self.find_next_cut(level, cut_value);
            let &max_cut_value = cut_value.iter().max().unwrap();
            let mut max_hierlabel = self.clusters_wedf[levelpts[0]];
            for i in &levelpts {
                if max_hierlabel < self.clusters_wedf[*i] {
                    max_hierlabel = self.clusters_wedf[*i];
                }
            }
            finished = max_cut_value == max_hierlabel;
        }

        self.hierarchy_level = level;
        self.level_values = cut_value;

        self.collect_detail_structures()
    }


    // COMPUTE WEDF AND EDF
    fn area_triangle(&self, ind_pts: &Vec<usize>) -> Vec<f32> {
        // Compute the area of each medial point in ind_pts
        let mut areas: Vec<f32> = vec![0.0; ind_pts.len()];
        for i in 0..self.triangles.len() {
            let in_ind_pts = functions_vec::find_vec(ind_pts, self.triangles[i][0] as usize);
            if !in_ind_pts.is_empty() {
                areas[in_ind_pts[0]] += Self::tri_area(&self.boundary_points[self.triangles[i][1] as usize], &self.boundary_points[self.triangles[i][2] as usize], &self.boundary_points[self.triangles[i][3] as usize]);
            }
        }
        areas
    }

    fn branches_on_comp(&self, si: i32) -> Vec<usize> {
        let mut boc: Vec<usize> = Vec::new();
        for i in 0..self.branch_number.len() {
            let comps = self.find_comp_of_br(i);
            for j in 0..comps.len() {
                if comps[j] == si {
                    boc.push(i);
                }
            }
        }
        boc
    }

    fn calculate_edf(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        self.calculate_wedf_edf(1, false)
    }

    fn calculate_wedf(&mut self, sum: bool) -> Result<(), Box<dyn std::error::Error>> {
        self.calculate_wedf_edf(0, sum)
    }
    
    fn calculate_wedf_edf(&mut self, wedf_edf: i32, sum: bool) -> Result<(), Box<dyn std::error::Error>> {
        // wedf_edf = 0 if wedf and = 1 if edf

        // Each WEDF or EDF is initially set to -Inf.
        if wedf_edf == 0 {
            self.wedf_array = vec![-f32::INFINITY; self.medial_points.len()];
        } else {
            self.edf_array = vec![-f32::INFINITY; self.medial_points.len()];
        }

        let mut burning_bma: bma = self.clone();

        // Finding the constrained ends
        let mut burning_ind: Vec<usize> = Vec::new();
        for i in 0..self.medial_points.len() {
            if Self::degree(i, &self.adjacency_matrix) == 1 {
                burning_ind.push(i);
            }
        }

        // Changed wedf or edf array
        let mut burning_edf: Vec<f32> = Vec::new();
        let mut burning_wedf: Vec<f32> = Vec::new();
        if wedf_edf == 0 {
            let temp_wedf = self.area_triangle(&burning_ind);
            for i in 0..burning_ind.len() {
                self.wedf_array[burning_ind[i]] = temp_wedf[i];
                burning_wedf.push(self.wedf_array[burning_ind[i]]);
            }
        } else {
            for i in 0..burning_ind.len() {
                self.edf_array[burning_ind[i]] = self.medial_radius[burning_ind[i]];
                burning_edf.push(self.edf_array[burning_ind[i]]);
            }
        }

        for br in 0..self.branch_number.len() {
            if burning_bma.is_stpl(br) {
                (burning_ind, burning_bma) = self.start_to_burn_stpl(wedf_edf, br, &burning_ind, &burning_bma);
                
                if wedf_edf == 0 {
                    burning_wedf = Vec::new();
                    for i in &burning_ind {
                        burning_wedf.push(self.wedf_array[*i]);
                    }
                } else {
                    burning_edf = Vec::new();
                    for i in &burning_ind {
                        burning_edf.push(self.edf_array[*i]);
                    }
                }
            }
        }

        let mut endloop = false;
        while !endloop {
            if burning_ind.is_empty() {
                // collect the remaining unburnt branches
                let mut ind_unburnt: Vec<usize> = Vec::new();
                for i in 0..burning_bma.branch_number.len() {
                    let mut zeros = true;
                    for j in 0..burning_bma.medial_points.len() {
                        zeros = zeros && (burning_bma.branch_number[i][j] == 0);
                    }
                    if !zeros {
                        ind_unburnt.push(i);
                    }
                }

                // keep track of handles that open loops, to update the outer boundary for next round
                let mut boundaries_toupdate: Vec<usize> = Vec::new();
                for i in &ind_unburnt {
                    if burning_bma.is_exposed_handle(*i) {
                        println!("is exposed handle");
                        (burning_ind, burning_bma) = self.start_to_burn_handle(wedf_edf, *i, &burning_ind, &burning_bma);
                        boundaries_toupdate.push(*i);
                    }
                }

                // since some indices are added, we need to update also burningWEDF/burningEDF
                if wedf_edf == 0 {
                    burning_wedf = Vec::new();
                    for i in &burning_ind {
                        burning_wedf.push(self.wedf_array[*i]);
                    }
                } else {
                    burning_edf = Vec::new();
                    for i in &burning_ind {
                        burning_edf.push(self.edf_array[*i]);
                    }
                }

                for i in &boundaries_toupdate {
                    burning_bma.update_boundary(*i);
                }
            }
            
            // burning the point of min WEDF or EDF in S
            // We find the new smallest in S=burningIndices ie x = argmin_{x in S} EDF(x)
            let smallest: f32;
            let index_of_smallest_in_burning: Vec<usize>;
            if wedf_edf == 0 {
                (smallest, index_of_smallest_in_burning) = functions_vec::min_f32(&burning_wedf);
            } else {
                (smallest, index_of_smallest_in_burning) = functions_vec::min_f32(&burning_edf);
            }
            let index_of_smallest_in_burning = index_of_smallest_in_burning[0];
            let index_of_smallest_in_points_array = burning_ind[index_of_smallest_in_burning];

            // burn the smallest (x)
            burning_bma.point_type[index_of_smallest_in_points_array] = -1;

            // Remove the smallest from S=burningIndices S = S \ {x}
            burning_ind.remove(index_of_smallest_in_burning);
            
            // Changed wedf array or edf array
            if wedf_edf == 0 {
                burning_wedf = Vec::new();
                for i in &burning_ind {
                    burning_wedf.push(self.wedf_array[*i]);
                }
            } else {
                burning_edf = Vec::new();
                for i in &burning_ind {
                    burning_edf.push(self.edf_array[*i]);
                }
            }

            // We find the index of the only point y adjacent to smallest x.
            let index_of_parent = functions_vec::find_vec(&burning_bma.adjacency_matrix[index_of_smallest_in_points_array], true);
            assert!(index_of_parent.len() == 1, "Zero or more than one parent at index.");

            // burn (x,y) and remove it from burningBma
            burning_bma.adjacency_matrix[index_of_smallest_in_points_array][index_of_parent[0]] = false;
            burning_bma.adjacency_matrix[index_of_parent[0]][index_of_smallest_in_points_array] = false;

            // update WEDF of y, that is, give it the value thru X if no other neighbor would have led to a greater WEDF
            let mut wedf_thru_x: f32 = 0.0;
            let mut edf_thru_x: f32 = 0.0;
            if wedf_edf == 0 {
                let area_parent = self.area_triangle(&index_of_parent);
                wedf_thru_x = smallest + area_parent[0];
            } else {
                edf_thru_x = smallest + ((self.medial_points[index_of_smallest_in_points_array][0]-self.medial_points[index_of_parent[0]][0]).powi(2)+(self.medial_points[index_of_smallest_in_points_array][1]-self.medial_points[index_of_parent[0]][1]).powi(2)).sqrt();
            }
            // fix for the original algo depends on that. Assign the WEDF value of y 
            if wedf_edf == 0 {
                self.update_wedf(index_of_parent[0], wedf_thru_x);
            } else {
                self.update_edf(index_of_parent[0], edf_thru_x);
            }

            // neighbors of y
            let pinds = functions_vec::find_vec(&burning_bma.adjacency_matrix[index_of_parent[0]], true);

            if pinds.len() == 1 { // if y has just 1 neighbor
                burning_ind.push(index_of_parent[0]); // add it to burn in S
                // Changed WEDF or EDF
                if wedf_edf == 0 {
                    burning_wedf = Vec::new();
                    for i in &burning_ind {
                        burning_wedf.push(self.wedf_array[*i]);
                    }
                } else {
                    burning_edf = Vec::new();
                    for i in &burning_ind {
                        burning_edf.push(self.edf_array[*i]);
                    }
                }
            }

            // if there is or was 3 branches finishing at y (y is or was a triple point)
            if self.point_type[index_of_parent[0]] == 3 {let deg_y = Self::degree(index_of_parent[0], &burning_bma.adjacency_matrix);
                if deg_y == 2 { // just drop to regular point
                    let ind_br_xy = burning_bma.find_branch(index_of_parent[0], index_of_smallest_in_points_array);
                    // remove the burnt branch from unburntBranchNumber if all interior points are burnt
                    if wedf_edf == 0 {
                        burning_bma = self.remove_branch_ifnomoreinsidepoint(wedf_edf, ind_br_xy[0], wedf_thru_x, &burning_bma);
                    } else {
                        burning_bma = self.remove_branch_ifnomoreinsidepoint(wedf_edf, ind_br_xy[0], edf_thru_x, &burning_bma);
                    }

                    // merge the two remaining branches that meet at y into a new branch in the first branch (min br #)
                    let new_br = burning_bma.merge_branches(index_of_parent[0]);

                    // and if the new branch if a STPL add it to burn
                    if new_br != 0 && burning_bma.is_stpl(new_br) {
                        (burning_ind, burning_bma) = self.start_to_burn_stpl(wedf_edf, new_br, &burning_ind, &burning_bma);
                        // Changed WEDF or EDF
                        if wedf_edf == 0 {
                            burning_wedf = Vec::new();
                            for i in &burning_ind {
                                burning_wedf.push(self.wedf_array[*i]);
                            }
                        } else {
                            burning_edf = Vec::new();
                            for i in &burning_ind {
                                burning_edf.push(self.edf_array[*i]);
                            }
                        }
                    }
                }
            }
            // We end the loop when temp.points is a single point
            let (ind1, ind2) = functions_vec::find_vecvec(&burning_bma.adjacency_matrix, true);
            if ind1.is_empty() {
                endloop = true;
            }
        }

        // Find the maximum value of each branch for wedf or edf
        if wedf_edf == 0 {
            for i in 0..self.medial_points.len() {
                for br in 0..self.branch_number.len() {
                    if self.branch_number[br][i] > 0 && self.wedf_array[i] > self.wedf_branch[br] {
                        self.wedf_branch[br] = self.wedf_array[i];
                    }
                }
            }
        } else {
            for i in 0..self.medial_points.len() {
                for br in 0..self.branch_number.len() {
                    if self.branch_number[br][i] > 0 && self.edf_array[i] > self.edf_branch[br] {
                        self.edf_branch[br] = self.edf_array[i];
                    }
                }
            }
        }

        Ok(())
    }

    fn degree(ind: usize, adj_matrix: &Vec<Vec<bool>>) -> i32 {
        let mut deg = 0;
        let mut degp = 0;
        for i in 0..adj_matrix.len() {
            if adj_matrix[ind][i] {
                deg += 1;
            }
            if adj_matrix[i][ind] {
                degp += 1;
            }
        }
        assert!(deg == degp, "unsymmetric adjacency matrix");
        deg
    }

    fn det(c1: &Vec<f32>, c2: &Vec<f32>) -> f32 {
        c1[0]*c2[1] - c1[1]*c2[0]
    }

    fn find_branch(&self, indx: usize, indy: usize) -> Vec<usize> {
        // branches containing x and y
        let mut ind_br_xy: Vec<usize> = Vec::new();
        for i in 0..self.branch_number.len() {
            if self.branch_number[i][indx] != 0 && self.branch_number[i][indy] != 0 {
                ind_br_xy.push(i);
            }
        }
        ind_br_xy
    }

    fn find_comp_of_br(&self, br: usize) -> Vec<i32> {
        let comps: Vec<i32> = Vec::new();
        let &max = self.branch_number[br].iter().max().unwrap();
        let ind_max = functions_vec::find_vec(&self.branch_number[br], max);
        let ind_max = ind_max[0];
        if max == 0 {
            // empty branch may happen when burning
            return comps;
        }
        
        let ind_first = functions_vec::find_vec(&self.branch_number[br], 1);
        let ind_first = ind_first[0];

        let mut c1: Vec<i32> = Vec::new();
        for i in 0..self.index_bndry_points[ind_first].len() {
            c1.push(self.component[self.index_bndry_points[ind_first][i] as usize]);
        }
        let mut cm: Vec<i32> = Vec::new();
        for i in 0..self.index_bndry_points[ind_max].len() {
            cm.push(self.component[self.index_bndry_points[ind_max][i] as usize]);
        }

        functions_vec::intersect(&c1, &cm)
    }

    fn flip_branch(&mut self, br: usize) -> Result<(), Box<dyn std::error::Error>> {
        let the_br = self.branch_number[br].clone();
        let max_br = the_br.iter().max().unwrap();
        for i in 0..self.branch_number[br].len() {
            if self.branch_number[br][i] != 0 {
                self.branch_number[br][i] = max_br + 1 - self.branch_number[br][i];
            }
        }
        Ok(())
    }

    fn is_exposed_handle(&self, br: usize) -> bool {
        // Output
        let is_eh: bool;

        let ind = functions_vec::find_vec(&self.branch_number[br], 0);
        let comp_br = self.find_comp_of_br(br);
        if ind.len() == self.branch_number[br].len() {
            // the branch has been burnt before
            is_eh = false;
        } else if comp_br.len() != 2 {
            // dangling branch or branch exposed both ways
            println!("pas assez de comp : {:?}", comp_br);
            is_eh = false;
        } else if self.is_outer_handle(br) {
            // We need to check that no other branch on the loop is also an ExposedHandle.
            assert!(comp_br.len()==2, "An exposed handle should have two components");
            assert!(self.is_outer_si(comp_br[0]) || self.is_outer_si(comp_br[1]), "... and one point should be outer at this point");
            let si_loop: i32;
            if self.is_outer_si(comp_br[0]) {
                si_loop = comp_br[1];
            } else {
                si_loop = comp_br[0];
            }

            // get the branches of the same components
            let boc = self.branches_on_comp(si_loop);

            for i in 0..boc.len() {
                if boc[i] != br {
                    if self.is_outer_handle(boc[i]) {
                        return false;
                    }
                }
            }
            // no other branch on this lop is outer
            is_eh = true;
        } else {
            is_eh = false;
        }

        is_eh
    }

    fn is_linked_to_outer_pt(&self, m: usize) -> bool {
        // tells if a medial point m is linked to a boundary point of the outer face
        let bnd_pts = &self.index_bndry_points[m];
        for i in 0..bnd_pts.len() {
            if self.isouter[bnd_pts[i] as usize] {
                return true;
            }
        }
        return false;
    }

    fn is_outer_handle(&self, br: usize) -> bool {
        // an handle is outer if it is a branch on a loop and its points are
        // outer, that is, they have a outer point on their boundary points.
        let mut exposed_points: Vec<bool> = Vec::new();

        // look if each point on the branch has an outer point in its corresponding boundary points
        for i in 0..self.branch_number[br].len() {
            if self.branch_number[br][i] != 0 {
                // for each point on the branch test if it is onloop and exposed
                let onloop_and_exposed = self.onloop[i] && self.is_linked_to_outer_pt(i);
                exposed_points.push(onloop_and_exposed);
            }
        }
        if exposed_points.is_empty() {
            return false;
        } else {
            let ind = functions_vec::find_vec(&exposed_points, false);
            if ind.is_empty() {
                return true;
            } else {
                return false;
            }
        }
    }

    fn is_outer_si(&self, si: i32) -> bool {
        // Return true if a point of the component si has at least one outer point
        for i in 0..self.component.len() {
            if self.component[i] == si && self.isouter[i] {
                return true;
            }
        }
        return false;
    }

    fn is_stpl(&self, br: usize) -> bool {
        let ind_first = functions_vec::find_vec(&self.branch_number[br], 1);
        assert!(ind_first.len() == 1);
        let max_val = self.branch_number[br].iter().max().unwrap();
        let ind_last = functions_vec::find_vec(&self.branch_number[br], *max_val);
        self.adjacency_matrix[ind_first[0]][ind_last[0]] && (Self::degree(ind_first[0], &self.adjacency_matrix) == 2 || Self::degree(ind_last[0], &self.adjacency_matrix) == 2)
    }

    fn merge_branches(&mut self, y: usize) -> usize {
        // Output
        let mut new_br: usize;

        assert!(Self::degree(y, &self.adjacency_matrix) == 2, "y does not have 2 neighbors");
        // maybe a branch is unburnt since another part of the branch did not burn yet, so rely on the nieghbors
        let ind_unburnt_neighbors = functions_vec::find_vec(&self.adjacency_matrix[y], true);
        assert!(ind_unburnt_neighbors.len() == 2, "for merging y should have 2 neighbors");

        let mut ind_unburntbr_y: Vec<usize> = Vec::new();
        for i in &ind_unburnt_neighbors {
            let mut vec_2 = self.find_branch(y, *i);
            ind_unburntbr_y.append(&mut vec_2);
        }
        
        let &min_br = ind_unburntbr_y.iter().min().unwrap();
        let &max_br = ind_unburntbr_y.iter().max().unwrap();

        if min_br == max_br { // last loop
            return 0;
        }

        if self.branch_number[min_br][y] == 1 { // make y the last in branch 1
            self.flip_branch(min_br);
        }
        if self.branch_number[max_br][y] != 1 { // make y the first in branch 2
            self.flip_branch(max_br);
            assert!(self.branch_number[max_br][y] == 1);
        }
        let first_in_min_br = functions_vec::find_vec(&self.branch_number[min_br], 1);
        let &last_ind_in_br_max = self.branch_number[max_br].iter().max().unwrap();
        let last_in_br_max = functions_vec::find_vec(&self.branch_number[max_br], last_ind_in_br_max);
        let in_loop = (first_in_min_br[0] == last_in_br_max[0]);

        // copy the rest of the max branch into the min branch
        let offset = self.branch_number[min_br][y] - 1; // now y should be the max index
        for i in 0..self.branch_number[min_br].len() {
            if self.branch_number[max_br][i] != 0 {
                self.branch_number[min_br][i] = offset + self.branch_number[max_br][i];
                self.branch_number[max_br][i] = 0;
            }
        }
        if in_loop {
            self.branch_number[min_br][first_in_min_br[0]] = 1;
        }

        // update the adjacency matrix for branches
        for i in 0..self.branch_adjacency[max_br].len() {
            if self.branch_adjacency[max_br][i] {
                self.branch_adjacency[min_br][i] = true;
                self.branch_adjacency[i][min_br] = true;
                self.branch_adjacency[max_br][i] = true;
                self.branch_adjacency[i][max_br] = true;
            }
        }
        self.branch_adjacency[min_br][min_br] = true;
        
        min_br
    }

    fn remove_branch_ifnomoreinsidepoint(&mut self, wedf_edf: i32, br_nb: usize, value_thru_x: f32, burn_bma: &bma) -> bma {
        // Outputs
        let mut bbma = burn_bma.clone();
        
        // if all interior point already burnt
        let mut ind_br_pa: Vec<usize> = Vec::new(); // indBrPA contains only interior points
        let mut ind_pstart: usize;
        let mut ind_pend: usize = 0;
        for i in 0..bbma.branch_number[br_nb].len() {
            if bbma.branch_number[br_nb][i] == 1 {
                ind_pstart = i;
                ind_pend = i;
            } else if bbma.branch_number[br_nb][i] == 2 {
                ind_pend = i;
            } else if bbma.branch_number[br_nb][i] != 0 {
                ind_br_pa.push(ind_pend);
                ind_pend = i;
            }
        }

        let mut br_interior_burnt = true;
        for i in &ind_br_pa {
            br_interior_burnt = br_interior_burnt && (bbma.point_type[*i] == -1);
        }

        if br_interior_burnt {
            // burn the whole branch (disappears from unburnt branches)
            bbma.branch_number[br_nb] = vec![0; bbma.medial_points.len()];
            
            for i in 0..bbma.branch_adjacency.len() {
                bbma.branch_adjacency[br_nb][i] = false;
                bbma.branch_adjacency[i][br_nb] = false;
            }
        }
        bbma
    }

    fn tri_area(c1: &Vector2<f32>, c2: &Vector2<f32>, c3: &Vector2<f32>) -> f32 {
        let p1 = vec![c3[0]-c1[0], c3[1]-c1[1]];
        let p2 = vec![c2[0]-c1[0], c2[1]-c1[1]];
        let res = Self::det(&p1, &p2)/2.0;
        res.abs()
    }

    fn start_to_burn_handle(&mut self, wedf_edf: i32, br_nb: usize, burning_ind: &Vec<usize>, burn_bma: &bma) -> (Vec<usize>, bma) {
        self.start_to_burn_loop(1, wedf_edf, br_nb, burning_ind, burn_bma)
    }

    fn start_to_burn_stpl(&mut self, wedf_edf: i32, br_nb: usize, burning_ind: &Vec<usize>, burn_bma: &bma) -> (Vec<usize>, bma) {
        self.start_to_burn_loop(0, wedf_edf, br_nb, burning_ind, burn_bma)
    }
    
    fn start_to_burn_loop(&mut self, handle_stpl: i32, wedf_edf: i32, br_nb: usize, burning_ind: &Vec<usize>, burn_bma: &bma) -> (Vec<usize>, bma) {
        // 1 for handle and 0 for STPL

        // Copy of inputs
        let mut bi = burning_ind.clone();
        let mut bbma = burn_bma.clone();
        
        let pts_a = self.medial_points.clone();

        // compute index of the midpoint z of a loop
        let mut ind1 = 1;
        let the_br = bbma.branch_number[br_nb].clone();
        let &(mut ind_max) = the_br.iter().max().unwrap();
        
        let mut len_from1 = 0.0;
        let mut len_from_max = 0.0;
        let ind1_pts_ind = functions_vec::find_vec(&the_br, 1);
        let ind1_pts_ind = ind1_pts_ind[0];
        let ind_max_pts_ind = functions_vec::find_vec(&the_br, ind_max);
        let ind_max_pts_ind = ind_max_pts_ind[0];

        if handle_stpl == 0 { // STPL
            // check the point 1 and max of theBr are neigbors
            assert!(bbma.adjacency_matrix[ind1_pts_ind][ind_max_pts_ind], "the first and last points of the stpl should be neighbors");
            len_from_max = ((pts_a[ind1_pts_ind][0]-pts_a[ind_max_pts_ind][0]).powi(2) + (pts_a[ind1_pts_ind][1]-pts_a[ind_max_pts_ind][1]).powi(2)).sqrt();
        } else { // handle
            if ind_max == 2 {
                assert!(bbma.adjacency_matrix[ind1_pts_ind][ind_max_pts_ind], "the first and last points of a single edge handle should be neighbors");
                if wedf_edf == 0 {
                    // halfArea : sum of the area of the two end triangles and divide by two
                    let areas = self.area_triangle(&vec![ind1_pts_ind, ind_max_pts_ind]);
                    let half_area = (areas[0] + areas[1])/2.0;
                    // changed WEDF
                    self.update_wedf(ind1_pts_ind, half_area);
                    self.update_wedf(ind_max_pts_ind, half_area);
                } else {
                    let half_length = ((pts_a[ind1_pts_ind][0]-pts_a[ind_max_pts_ind][0]).powi(2) + (pts_a[ind1_pts_ind][1]-pts_a[ind_max_pts_ind][1]).powi(2)).sqrt() / 2.0;
                    // changed EDF
                    self.update_edf(ind1_pts_ind, half_length);
                    self.update_edf(ind_max_pts_ind, half_length);
                }
                bbma.adjacency_matrix[ind1_pts_ind][ind_max_pts_ind] = false;
                bbma.adjacency_matrix[ind_max_pts_ind][ind1_pts_ind] = false;

                // return a value maybe
                return (bi, bbma);
            }
        }

        while ind1 != ind_max {
            if len_from1 < len_from_max {
                let first_pt = functions_vec::find_vec(&the_br, ind1);
                let second_pt = functions_vec::find_vec(&the_br, ind1+1);
                len_from1 += ((pts_a[first_pt[0]][0]-pts_a[second_pt[0]][0]).powi(2) + (pts_a[first_pt[0]][1]-pts_a[second_pt[0]][1]).powi(2)).sqrt();
                ind1 += 1;
            } else {
                let last_pt = functions_vec::find_vec(&the_br, ind_max);
                let but_last_pt = functions_vec::find_vec(&the_br, ind_max-1);
                len_from_max += ((pts_a[last_pt[0]][0]-pts_a[but_last_pt[0]][0]).powi(2) + (pts_a[last_pt[0]][1]-pts_a[but_last_pt[0]][1]).powi(2)).sqrt();
                ind_max -= 1;
            }
        }

        let z_ind_in_the_br = ind1;
        let z = functions_vec::find_vec(&the_br, z_ind_in_the_br);
        let z = z[0];
        assert!(Self::degree(z, &bbma.adjacency_matrix) == 2, "midPoint of STPL should be of degree 2");

        // find the greatest WEDF/EDF value in the branch
        let mut start_wedf = -1.0;
        let mut start_edf = -1.0;
        if wedf_edf == 0 {
            for i in 0..the_br.len() {
                if the_br[i] != 0 && self.wedf_array[i]>start_wedf {
                    start_wedf = self.wedf_array[i];
                }
            }
    
            if start_wedf < 0.0 { // all points of theBr have WEDF = -Inf 
                let temp = self.area_triangle(&vec![z]);
                start_wedf = temp[0];
            }
        } else {
            for i in 0..the_br.len() {
                if the_br[i] != 0 && self.edf_array[i]>start_edf {
                    start_edf = self.edf_array[i];
                }
            }
            if start_edf < 0.0 { // all points of theBr have WEDF = -Inf
                start_edf = self.medial_radius[z];
            }
        }

        // neighbors of z
        let z0 = functions_vec::find_vec(&the_br, z_ind_in_the_br-1);
        let mut z0 = z0[0];
        let z1 = functions_vec::find_vec(&the_br, z_ind_in_the_br+1);
        let z1 = z1[0];
        let deg_z0 = Self::degree(z0, &bbma.adjacency_matrix);
        let deg_z1 = Self::degree(z1, &bbma.adjacency_matrix);
        let mut dist0 = ((pts_a[z0][0]-pts_a[z][0]).powi(2) + (pts_a[z0][1]-pts_a[z][1]).powi(2)).sqrt();
        let dist1 = ((pts_a[z1][0]-pts_a[z][0]).powi(2) + (pts_a[z1][1]-pts_a[z][1]).powi(2)).sqrt();

        if deg_z0>2 && deg_z1>2 {
            if handle_stpl == 0 {
                assert!(true, "is you come here it means there was a STPL with just one inside point, bummer !");
            }
            // it means that we still need to erase the branch and merge the two branches attached to the neigbors of z
            bbma.adjacency_matrix[z][z0] = false;
            bbma.adjacency_matrix[z0][z] = false;
            bbma.adjacency_matrix[z][z1] = false;
            bbma.adjacency_matrix[z1][z] = false;
            bbma.point_type[z] = -1;

            if wedf_edf == 0 {
                // changed WEDF
                self.update_wedf(z, start_wedf);
                self.update_wedf(z0, start_wedf + dist0);
                self.update_wedf(z1, start_wedf + dist1);
            } else {
                // changed EDF
                self.update_edf(z, start_edf);
                self.update_edf(z0, start_edf + dist0);
                self.update_edf(z1, start_edf + dist1);
            }
        } else {
            // find the nearest of the two neighbors and assign it to z0
            if dist1<dist0 || deg_z0>2 {
                assert!(Self::degree(z1, &bbma.adjacency_matrix)==2, "beware z should be the only point on the branch and we should have caught that earlier");
                z0 = z1;
                dist0 = dist1;
            }
            bbma.adjacency_matrix[z][z0] = false;
            bbma.adjacency_matrix[z0][z] = false;

            if wedf_edf == 0 {
                // changed WEDF
                self.update_wedf(z, start_wedf);
                self.update_wedf(z0, start_wedf + dist0);
            } else {
                // changed EDF
                self.update_edf(z, start_edf);
                self.update_edf(z0, start_edf + dist0);
            }

            bi.push(z);
            bi.push(z0);
        }

        (bi, bbma)
    }

    fn update_boundary(&mut self, hnb: usize) -> Result<(), Box<dyn std::error::Error>> {
        let comps = self.find_comp_of_br(hnb);

        // for this handle one of the component should be outer
        assert!(comps.len() == 2, "comps not of length 2");
        assert!(self.is_outer_si(comps[0]) || self.is_outer_si(comps[1]), "none of Si is outer ?");
        
        let mut changing_component = comps[0];
        if self.is_outer_si(comps[0]) {
            changing_component = comps[1];
        }

        let ind_bpcg = functions_vec::find_vec(&self.component, changing_component);
        for i in ind_bpcg {
            self.isouter[i] = true; // change to outer
        }

        Ok(())
    }

    fn update_edf(&mut self, ind: usize, new_edf: f32) -> Result<(), Box<dyn std::error::Error>> {
        if self.edf_array[ind] < new_edf {
            self.edf_array[ind] = new_edf;
        }

        Ok(())
    }

    fn update_wedf(&mut self, ind: usize, new_wedf: f32) -> Result<(), Box<dyn std::error::Error>> {
        if self.wedf_array[ind] < new_wedf {
            self.wedf_array[ind] = new_wedf;
        }

        Ok(())
    }
 

    // COMPUTE ET AND ST
    fn calculate_et(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        for i in 0..self.medial_points.len() {
            let et = self.edf_array[i] - self.medial_radius[i];
            self.erosion_thickness.push(et);
        }

        // Scale invariant ?
        let (max_et, _ind) = functions_vec::max_f32(&self.erosion_thickness);
        let (min_et, _ind) = functions_vec::min_f32(&self.erosion_thickness);
        for i in 0..self.medial_points.len() {
            let mut et = (self.erosion_thickness[i]-min_et)/(max_et-min_et);
            et = et.exp();
            self.erosion_thickness.push(et);
        }
        Ok(())
    }
    
    fn calculate_st(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        for i in 0..self.medial_points.len() {
            let mut st = (self.edf_array[i] - self.medial_radius[i])/self.edf_array[i];
            st = st.exp();
            self.shape_tubularity.push(st);
        }
        /*
        // Scale invariant ?
        let (max_st, _ind) = functions_vec::max_f32(&self.shape_tubularity);
        let (min_st, _ind) = functions_vec::min_f32(&self.shape_tubularity);
        for i in 0..self.medial_points.len() {
            let mut st = (self.shape_tubularity[i]-min_st)/(max_st-min_st);
            st = st.exp();
            self.shape_tubularity.push(st);
        }*/
        Ok(())
    }


    // SIMILARITY
    fn resample_trunk(&self, nb_pts: usize, ind_trunk: Vec<usize>) -> Vec<Vector2<f32>> {
        let mut rem_nb_pts = nb_pts - 1;
        let mut new_pts: Vec<Vector2<f32>> = Vec::new();
        let mut len_trunk = 0.0;
        let mut len_edges: Vec<f32> = Vec::new();

        for i in 1..ind_trunk.len() {
            len_edges.push(((self.medial_points[ind_trunk[i-1]][0] - self.medial_points[ind_trunk[i]][0]).powi(2) + (self.medial_points[ind_trunk[i-1]][1] - self.medial_points[ind_trunk[i]][1]).powi(2)).sqrt());
            len_trunk += len_edges[len_edges.len() - 1];
        }
        
        let mut nb_pts_edges: Vec<f32> = Vec::new();
        for i in &len_edges {
            nb_pts_edges.push((*i/len_trunk*(nb_pts as f32 - 1.0)).round());
            rem_nb_pts -= nb_pts_edges[nb_pts_edges.len() - 1] as usize;
        }

        let mut ind = 0;
        while rem_nb_pts != 0 {
            ind = ind + 1;
            if ind == nb_pts_edges.len() {
                ind = 0;
            }
            nb_pts_edges[ind] += 1.0;
            rem_nb_pts -= 1;
        }

        for i in 0..nb_pts_edges.len() {
            let space_x = (self.medial_points[ind_trunk[i+1]][0] - self.medial_points[ind_trunk[i]][0])/nb_pts_edges[i];
            let space_y = (self.medial_points[ind_trunk[i+1]][1] - self.medial_points[ind_trunk[i]][1])/nb_pts_edges[i];
            for j in 0..nb_pts_edges[i] as usize {
                new_pts.push(Vector2::new(self.medial_points[i][0] + (j as f32)*space_x, self.medial_points[i][1] + (j as f32)*space_y));
            }
        }
        new_pts.push(Vector2::new(self.medial_points[nb_pts_edges.len()][0], self.medial_points[nb_pts_edges.len()][1]));

        new_pts
    } 

    fn get_st_et_of_trunk(&self, nb_pts: usize, ind_trunk: &Vec<usize>) -> (Vec<f32>, Vec<f32>) {
        let mut et_values: Vec<f32> = Vec::new();
        let mut st_values: Vec<f32> = Vec::new();

        // Trunk of 1 point
        if ind_trunk.len() == 1 {
            et_values = vec![self.erosion_thickness[ind_trunk[0]]; nb_pts];
            st_values = vec![self.shape_tubularity[ind_trunk[0]]; nb_pts];
            return (et_values, st_values)
        }

        let mut rem_nb_pts = nb_pts - 1;
        let mut len_trunk = 0.0;
        let mut len_edges: Vec<f32> = Vec::new();
        
        for i in 1..ind_trunk.len() {
            len_edges.push(((self.medial_points[ind_trunk[i-1]][0] - self.medial_points[ind_trunk[i]][0]).powi(2) + (self.medial_points[ind_trunk[i-1]][1] - self.medial_points[ind_trunk[i]][1]).powi(2)).sqrt());
            len_trunk += len_edges[len_edges.len() - 1];
        }
        
        let mut nb_pts_edges: Vec<f32> = Vec::new();
        for i in &len_edges {
            nb_pts_edges.push((*i/len_trunk*(nb_pts as f32 - 1.0)).floor());
            rem_nb_pts -= nb_pts_edges[nb_pts_edges.len() - 1] as usize;
        }

        let mut ind = 0;
        while rem_nb_pts != 0 {
            ind = ind + 1;
            if ind == nb_pts_edges.len() {
                ind = 0;
            }
            nb_pts_edges[ind] += 1.0;
            rem_nb_pts -= 1;
        }

        for i in 0..nb_pts_edges.len() {
            let et_interval = (self.erosion_thickness[ind_trunk[i+1]] - self.erosion_thickness[ind_trunk[i]])/nb_pts_edges[i];
            let st_interval = (self.shape_tubularity[ind_trunk[i+1]] - self.shape_tubularity[ind_trunk[i]])/nb_pts_edges[i];
            for j in 0..nb_pts_edges[i] as usize {
                et_values.push(self.erosion_thickness[ind_trunk[i]] + (j as f32)*et_interval);
                st_values.push(self.shape_tubularity[ind_trunk[i]] + (j as f32)*st_interval);
            }
        }
        et_values.push(self.erosion_thickness[ind_trunk[nb_pts_edges.len()]]);
        st_values.push(self.shape_tubularity[ind_trunk[nb_pts_edges.len()]]);

        (et_values, st_values)
    } 

    fn normalize_values(&self, values: &Vec<f32>, et_or_st: bool) -> Vec<f32> {
        // We want et and st to be independant of the scale
        let mut normalized_values: Vec<f32> = Vec::new();
        
        let min;
        let ind_min;
        let max;
        let ind_max;
        if et_or_st {
            (max, ind_max) = functions_vec::max_f32(&self.erosion_thickness);
            (min, ind_min) = functions_vec::min_f32(&self.erosion_thickness);
        } else {
            (max, ind_max) = functions_vec::max_f32(&self.shape_tubularity);
            (min, ind_min) = functions_vec::min_f32(&self.shape_tubularity);
        }

        for i in 0..values.len() {
            normalized_values.push((values[i]-min)/(max-min));
        }

        normalized_values
    }

    pub fn clustering_etst(&self, cluster_min: usize, cluster_max: usize, normalized: bool) ->
    (Vec<Vec<f32>>, Vec<usize>) {
        // INPUTS :
        // mainshapelabels: output of hierarchy labels from hierarchyWEDF
        // cluster_min: the lowest number of clusters you want
        // cluster_max: the highest number of clusters you want

        // Outputs
        let mut centroids: Vec<Vec<f32>> = Vec::new();

        let mut nb_pts = 0;
        for i in 0..self.trunk.len() {
            if self.trunk[i].indices.len() > nb_pts {
                nb_pts = self.trunk[i].indices.len();
            }
        }

        // Resample ET and ST values of each trunk
        let mut resampled_et_values: Vec<Vec<f32>> = Vec::new();
        let mut resampled_st_values: Vec<Vec<f32>> = Vec::new();
        for i in 0..self.trunk.len() {
            let (et, st) = self.get_st_et_of_trunk(nb_pts, &self.trunk[i].indices);
            resampled_et_values.push(et);
            resampled_st_values.push(st);
        }

        // All the values
        let mut data_normalized: Vec<Vec<f32>> = Vec::new();
        for i in 0..self.trunk.len() {
            data_normalized.push(Vec::new());
            if normalized {
                let mut et_normalized = self.normalize_values(&resampled_et_values[i], true);
                let mut st_normalized = self.normalize_values(&resampled_st_values[i], false);
                for et in &et_normalized {
                    data_normalized[i].push(*et);
                }
                for st in &st_normalized {
                    data_normalized[i].push(*st);
                }
            } else {
                for et in &resampled_et_values[i] {
                    data_normalized[i].push(*et);
                }
                for st in &resampled_st_values[i] {
                    data_normalized[i].push(*st);
                }
            }
        }
        
        // Find the number of clusters
        let optimal_k: usize;
        assert!(cluster_max >= cluster_min);
        if cluster_min != cluster_max {
            optimal_k = gap_evaluation::gap_eval(&data_normalized, cluster_min, cluster_max);
        } else {
            optimal_k = cluster_min;
        }
        //let seeds = kmeans::initialize(optimal_k, &data_normalized);
        let seeds = kmeans::init_percentiles_etst(optimal_k, &data_normalized);
        let clustering = kmeans::kmeans(optimal_k, &data_normalized, 200, seeds);
        
        // Get the centroids
        for i in 0..clustering.centroids.len() {
            centroids.push(Vec::new());
            for j in 0..clustering.centroids[i].0.len() {
                centroids[i].push(clustering.centroids[i].0[j] as f32);
            }
        }

        // OUTPUTS
        // centroids: the centroids of the clusters from the clustering
        // clustering.membership: cluster labels
        (centroids, clustering.membership)
    }


    

    // FUNCTIONS FOR PLOTTING

    // To plot the shape of an object
    pub fn plot_shape(&self, original_image:&str, new_image: &str, show_medial_axis: bool) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        if show_medial_axis {
            // Drawing the medial axis edges
            for i in 0..self.medial_edges.len() {
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.medial_points[self.medial_edges[i][0] as usize][0], self.medial_points[self.medial_edges[i][0] as usize][1]),
                    (self.medial_points[self.medial_edges[i][1] as usize][0], self.medial_points[self.medial_edges[i][1] as usize][1]),
                    image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]),
                );
            }
        }

        // Save the image
        img_rgb.save(new_image).unwrap();
    }


    // To plot the boundary of an object
    pub fn plot_boundary(&self, original_image:&str, new_image: &str, show_component: bool) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        let red = image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]);
        let blue = image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]);
        let mut colors: Vec<image::Rgba::<u8>> = Vec::new();
        if show_component {
            colors.push(red);
            colors.push(blue);
            let green = image::Rgba::<u8>([0u8, 255u8, 0u8, 255u8]);
            colors.push(green);
            let yellow = image::Rgba::<u8>([255u8, 255u8, 0u8, 255u8]);
            colors.push(yellow);
            let magenta = image::Rgba::<u8>([255u8, 0u8, 255u8, 255u8]);
            colors.push(magenta);
            let cyan = image::Rgba::<u8>([0u8, 255u8, 255u8, 255u8]);
            colors.push(cyan);
            let pink = image::Rgba::<u8>([255u8, 192u8, 203u8, 255u8]);
            colors.push(pink);
            let orange = image::Rgba::<u8>([255u8, 165u8, 0u8, 255u8]);
            colors.push(orange);
            let purple = image::Rgba::<u8>([127u8, 127u8, 255u8, 255u8]);
            colors.push(purple);
        }

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            if show_component {
                // Choose the color for the component
                let mut num_color = self.component[i];
                while num_color>=colors.len().try_into().unwrap() {
                    num_color -= colors.len() as i32; 
                }
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                    (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                    colors[num_color as usize],
                );
            } else {
                if self.isouter[i] {
                    imageproc::drawing::draw_line_segment_mut(
                        &mut img_rgb,
                        (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                        (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                        red,
                    );
                } else {
                    imageproc::drawing::draw_line_segment_mut(
                        &mut img_rgb,
                        (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                        (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                        blue,
                    );
                }
            }
        }

        // Save the image
        img_rgb.save(new_image).unwrap();
    }

    
    // To plot the medial axis of an object
    pub fn plot_bma_loops(&self, original_image:&str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        let red = image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]);
        let purple = image::Rgba::<u8>([127u8, 127u8, 255u8, 255u8]);

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        // Drawing the medial axis edges
        for i in 0..self.medial_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.medial_points[self.medial_edges[i][0] as usize][0], self.medial_points[self.medial_edges[i][0] as usize][1]),
                (self.medial_points[self.medial_edges[i][1] as usize][0], self.medial_points[self.medial_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        // Drawing the medial axis points
        for i in 0..self.medial_points.len() {
            if self.onloop[i] {
                imageproc::drawing::draw_cross_mut(
                    &mut img_rgb,
                    red,
                    self.medial_points[i][0] as i32,
                    self.medial_points[i][1] as i32,
                );
            } else {
                imageproc::drawing::draw_cross_mut(
                    &mut img_rgb,
                    purple,
                    self.medial_points[i][0] as i32,
                    self.medial_points[i][1] as i32,
                );
            }
        }

        // Save the image
        img_rgb.save(new_image).unwrap();
    }

    
    // To plot the medial axis of an object
    pub fn plot_bma_pointtypes(&self, original_image:&str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        let red = image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]);
        let green = image::Rgba::<u8>([0u8, 255u8, 0u8, 255u8]);
        let orange = image::Rgba::<u8>([255u8, 165u8, 0u8, 255u8]);
        let blue = image::Rgba::<u8>([0u8, 165u8, 255u8, 255u8]);

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        // Drawing the medial axis edges
        for i in 0..self.medial_edges.len() {
            let mut color = orange;
            if (self.point_type[self.medial_edges[i][0] as usize]==3 && self.istrip_nbr[self.medial_edges[i][1] as usize]) ||
            (self.point_type[self.medial_edges[i][1] as usize]==3 && self.istrip_nbr[self.medial_edges[i][0] as usize]) {
                color = blue;
            }
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.medial_points[self.medial_edges[i][0] as usize][0], self.medial_points[self.medial_edges[i][0] as usize][1]),
                (self.medial_points[self.medial_edges[i][1] as usize][0], self.medial_points[self.medial_edges[i][1] as usize][1]),
                color,
            );
        }
        
        // Drawing the medial axis points
        for i in 0..self.medial_points.len() {
            if self.point_type[i] == 1 {
                imageproc::drawing::draw_cross_mut(
                    &mut img_rgb,
                    red,
                    self.medial_points[i][0] as i32,
                    self.medial_points[i][1] as i32,
                );
            } else if self.point_type[i] == 3 {
                imageproc::drawing::draw_cross_mut(
                    &mut img_rgb,
                    green,
                    self.medial_points[i][0] as i32,
                    self.medial_points[i][1] as i32,
                );
            }
        }

        // Save the image
        img_rgb.save(new_image).unwrap();
    }

    
    // To plot the different branches of a medial axis
    pub fn plot_bma_branches(&self, original_image:&str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        let mut colors: Vec<image::Rgba::<u8>> = Vec::new();
        let yellow = image::Rgba::<u8>([255u8, 255u8, 0u8, 255u8]);
        colors.push(yellow);
        let green = image::Rgba::<u8>([0u8, 255u8, 0u8, 255u8]);
        colors.push(green);
        let magenta = image::Rgba::<u8>([255u8, 0u8, 255u8, 255u8]);
        colors.push(magenta);
        let cyan = image::Rgba::<u8>([0u8, 255u8, 255u8, 255u8]);
        colors.push(cyan);
        let red = image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]);
        colors.push(red);
        let pink = image::Rgba::<u8>([255u8, 150u8, 180u8, 255u8]);
        colors.push(pink);
        let blue = image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]);
        colors.push(blue);
        let orange = image::Rgba::<u8>([255u8, 165u8, 0u8, 255u8]);
        colors.push(orange);
        let purple = image::Rgba::<u8>([127u8, 127u8, 255u8, 255u8]);
        colors.push(purple);
        let dark_green = image::Rgba::<u8>([0u8, 100u8, 0u8, 255u8]);
        colors.push(dark_green);
        
        // Drawing the medial axis edges
        for num_branch in 0..self.branch_number.len() {
            let mut num_pt = 1;
            let mut ind1 = functions_vec::find_vec(&self.branch_number[num_branch], num_pt);
            let mut ind2 = functions_vec::find_vec(&self.branch_number[num_branch], num_pt+1);

            while !ind1.is_empty() && !ind2.is_empty() {
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.medial_points[ind1[0]][0], self.medial_points[ind1[0]][1]),
                    (self.medial_points[ind2[0]][0], self.medial_points[ind2[0]][1]),
                    colors[num_branch%colors.len()],
                );
                num_pt += 1;
                ind1 = functions_vec::find_vec(&self.branch_number[num_branch], num_pt);
                ind2 = functions_vec::find_vec(&self.branch_number[num_branch], num_pt+1);
            }
        }

        /*
        let startpt = 420;
        let startptnbr = 316;

        imageproc::drawing::draw_cross_mut(
            &mut img_rgb,
            image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            self.medial_points[startpt][0] as i32,
            self.medial_points[startpt][1] as i32,
        );
        imageproc::drawing::draw_cross_mut(
            &mut img_rgb,
            image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]),
            self.medial_points[startptnbr][0] as i32,
            self.medial_points[startptnbr][1] as i32,
        );*/
       
        // Save the image
        img_rgb.save(new_image).unwrap();
    }

        
    // To plot the different trunks of a medial axis
    pub fn plot_bma_trunks(&self, level: i32, original_image: &str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        // Drawing the medial axis edges
        for i in 0..self.medial_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.medial_points[self.medial_edges[i][0] as usize][0], self.medial_points[self.medial_edges[i][0] as usize][1]),
                (self.medial_points[self.medial_edges[i][1] as usize][0], self.medial_points[self.medial_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }

        // Drawing the trunks in red
        for i in 0..self.trunk.len() {
            if self.trunk[i].level == level {
                for j in 0..self.trunk[i].indices.len()-1 {
                    imageproc::drawing::draw_line_segment_mut(
                        &mut img_rgb,
                        (self.medial_points[self.trunk[i].indices[j]][0], self.medial_points[self.trunk[i].indices[j]][1]),
                        (self.medial_points[self.trunk[i].indices[j+1]][0], self.medial_points[self.trunk[i].indices[j+1]][1]),
                        image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]),
                    );
                }
            }
        }
       
        // Save the image
        img_rgb.save(new_image).unwrap();
    }


    // To plot one trunk of a medial axis
    pub fn plot_bma_one_trunk(&self, level: i32, num_trunk: i32, original_image: &str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        // Drawing the medial axis edges
        for i in 0..self.medial_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.medial_points[self.medial_edges[i][0] as usize][0], self.medial_points[self.medial_edges[i][0] as usize][1]),
                (self.medial_points[self.medial_edges[i][1] as usize][0], self.medial_points[self.medial_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }

        // Drawing the trunk in red
        let mut num = 0;
        for i in 0..self.trunk.len() {
            if self.trunk[i].level == level {
                if num == num_trunk {
                    for j in 0..self.trunk[i].indices.len()-1 {
                        imageproc::drawing::draw_line_segment_mut(
                            &mut img_rgb,
                            (self.medial_points[self.trunk[i].indices[j]][0], self.medial_points[self.trunk[i].indices[j]][1]),
                            (self.medial_points[self.trunk[i].indices[j+1]][0], self.medial_points[self.trunk[i].indices[j+1]][1]),
                            image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]),
                        );
                    }
                }
                num += 1;
            }
        }
        
        // Save the image
        img_rgb.save(new_image).unwrap();
    }
    
    
    // To plot all the trunks of a medial axis
    pub fn plot_bma_all_trunks(&self, original_image: &str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        // Drawing the medial axis edges
        for i in 0..self.medial_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.medial_points[self.medial_edges[i][0] as usize][0], self.medial_points[self.medial_edges[i][0] as usize][1]),
                (self.medial_points[self.medial_edges[i][1] as usize][0], self.medial_points[self.medial_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }

        // Extracting the different levels of the trunks
        let mut levels: Vec<i32> = Vec::new();
        for i in 0..self.trunk.len() {
            levels.push(self.trunk[i].level);
        }
        let levels = functions_vec::unique_list(&levels);
        
        // Change the interval to get the gradient of colors
        let mut new_values: Vec<f32> = Vec::new();
        for i in 0..levels.len() {
            new_values.push(levels[i] as f32/levels[levels.len() - 1] as f32);
        }

        let mut ind_level: i32 = 0;
        while ind_level <= levels.len() as i32 - 1 {
            let level = levels[ind_level as usize];
            for i in 0..self.trunk.len() {
                let value = (new_values[ind_level as usize]*1020.0) as i32;
                let mut color = image::Rgba::<u8>([0, 0, 0, 255]);
                if value <= 255 {
                    color = image::Rgba::<u8>([0, value as u8, 255, 255]);
                } else if value <= 510 {
                    color = image::Rgba::<u8>([0, 255, (510-value) as u8, 255]);
                } else if value <= 765 {
                    color = image::Rgba::<u8>([(value-510) as u8, 255, 0, 255]);
                } else if value <= 1020 {
                    color = image::Rgba::<u8>([255, (1020-value) as u8, 0, 255]);
                }
                if self.trunk[i].level == level {
                    for j in 0..self.trunk[i].indices.len()-1 {
                        imageproc::drawing::draw_line_segment_mut(
                            &mut img_rgb,
                            (self.medial_points[self.trunk[i].indices[j]][0], self.medial_points[self.trunk[i].indices[j]][1]),
                            (self.medial_points[self.trunk[i].indices[j+1]][0], self.medial_points[self.trunk[i].indices[j+1]][1]),
                            color,
                        );
                    }
                }
            }
            ind_level += 1;
        }
       
        // Save the image
        img_rgb.save(new_image).unwrap();
    }
   
    
    // To plot the neighbours a branch
    pub fn plot_bma_adjacency_branches(&self, num_branch: usize, original_image: &str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        let mut colors: Vec<image::Rgba::<u8>> = Vec::new();
        colors.push(image::Rgba::<u8>([50u8, 0u8, 0u8, 255u8]));
        colors.push(image::Rgba::<u8>([100u8, 0u8, 0u8, 255u8]));
        colors.push(image::Rgba::<u8>([150u8, 0u8, 0u8, 255u8]));
        colors.push(image::Rgba::<u8>([200u8, 0u8, 0u8, 255u8]));
        colors.push(image::Rgba::<u8>([250u8, 0u8, 0u8, 255u8]));

        let mut nb_neighbours = 0;
        // Drawing the branch in green and its neighbours in red
        for i in 0..self.branch_adjacency[num_branch].len() {
            let mut color = image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]);
            if i == num_branch {
                color = image::Rgba::<u8>([0u8, 255u8, 0u8, 255u8]);
            } else if self.branch_adjacency[num_branch][i] {
                color = colors[nb_neighbours%5];
                nb_neighbours = nb_neighbours + 1;
            }

            let mut num_pt = 1;
            let mut ind1 = functions_vec::find_vec(&self.branch_number[i], num_pt);
            let mut ind2 = functions_vec::find_vec(&self.branch_number[i], num_pt+1);

            while !ind1.is_empty() && !ind2.is_empty() {
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.medial_points[ind1[0]][0], self.medial_points[ind1[0]][1]),
                    (self.medial_points[ind2[0]][0], self.medial_points[ind2[0]][1]),
                    color,
                );
                num_pt += 1;
                ind1 = functions_vec::find_vec(&self.branch_number[i], num_pt);
                ind2 = functions_vec::find_vec(&self.branch_number[i], num_pt+1);
            }
        }
       
        // Save the image
        img_rgb.save(new_image).unwrap();
    }

    
    // To plot the triangles of one of the medial axis point
    pub fn plot_triangle_of_ma_pt(&self, num_ma: usize, original_image:&str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Drawing the medial axis edges
        for i in 0..self.medial_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.medial_points[self.medial_edges[i][0] as usize][0], self.medial_points[self.medial_edges[i][0] as usize][1]),
                (self.medial_points[self.medial_edges[i][1] as usize][0], self.medial_points[self.medial_edges[i][1] as usize][1]),
                image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]),
            );
        }

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }

        for i in 0..self.triangles.len() {
            if self.triangles[i][0] == num_ma as i32 {
                // Drawing the edges of the triangle
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.boundary_points[self.triangles[i][1] as usize][0], self.boundary_points[self.triangles[i][1] as usize][1]),
                    (self.boundary_points[self.triangles[i][2] as usize][0], self.boundary_points[self.triangles[i][2] as usize][1]),
                    image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
                );
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.boundary_points[self.triangles[i][2] as usize][0], self.boundary_points[self.triangles[i][2] as usize][1]),
                    (self.boundary_points[self.triangles[i][3] as usize][0], self.boundary_points[self.triangles[i][3] as usize][1]),
                    image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
                );
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.boundary_points[self.triangles[i][3] as usize][0], self.boundary_points[self.triangles[i][3] as usize][1]),
                    (self.boundary_points[self.triangles[i][1] as usize][0], self.boundary_points[self.triangles[i][1] as usize][1]),
                    image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
                );
            }
        }
        // Drawing the medial axis point associated
        imageproc::drawing::draw_cross_mut(
            &mut img_rgb,
            image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
            self.medial_points[num_ma][0] as i32,
            self.medial_points[num_ma][1] as i32,
        );
        // Save the image
        img_rgb.save(new_image).unwrap();
    }
    

    // To plot the triangles of an object
    pub fn plot_triangles(&self, original_image:&str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Drawing the triangles
        for i in 0..self.triangles.len() {
            // Drawing the edges of each triangle
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.triangles[i][1] as usize][0], self.boundary_points[self.triangles[i][1] as usize][1]),
                (self.boundary_points[self.triangles[i][2] as usize][0], self.boundary_points[self.triangles[i][2] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
            );
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.triangles[i][2] as usize][0], self.boundary_points[self.triangles[i][2] as usize][1]),
                (self.boundary_points[self.triangles[i][3] as usize][0], self.boundary_points[self.triangles[i][3] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
            );
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.triangles[i][3] as usize][0], self.boundary_points[self.triangles[i][3] as usize][1]),
                (self.boundary_points[self.triangles[i][1] as usize][0], self.boundary_points[self.triangles[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
            );
            // Drawing the medial axis point associated
            imageproc::drawing::draw_cross_mut(
                &mut img_rgb,
                image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]),
                self.medial_points[self.triangles[i][0] as usize][0] as i32,
                self.medial_points[self.triangles[i][0] as usize][1] as i32,
            );
        }
        // Save the image
        img_rgb.save(new_image).unwrap();
    }
    

    // To plot the triangles of an object
    pub fn plot_triangle_i(&self, original_image:&str, new_image: &str, ma_pt_i: usize) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Checking for each triangle if it is associated with the bma point
        for i in 0..self.triangles.len() {
            if self.triangles[i][0] == ma_pt_i as i32 {
                // Drawing the edges of the triangle
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.boundary_points[self.triangles[i][1] as usize][0], self.boundary_points[self.triangles[i][1] as usize][1]),
                    (self.boundary_points[self.triangles[i][2] as usize][0], self.boundary_points[self.triangles[i][2] as usize][1]),
                    image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
                );
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.boundary_points[self.triangles[i][2] as usize][0], self.boundary_points[self.triangles[i][2] as usize][1]),
                    (self.boundary_points[self.triangles[i][3] as usize][0], self.boundary_points[self.triangles[i][3] as usize][1]),
                    image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
                );
                imageproc::drawing::draw_line_segment_mut(
                    &mut img_rgb,
                    (self.boundary_points[self.triangles[i][3] as usize][0], self.boundary_points[self.triangles[i][3] as usize][1]),
                    (self.boundary_points[self.triangles[i][1] as usize][0], self.boundary_points[self.triangles[i][1] as usize][1]),
                    image::Rgba::<u8>([0u8, 0u8, 255u8, 255u8]),
                );
                // Drawing the medial axis point associated
                imageproc::drawing::draw_cross_mut(
                    &mut img_rgb,
                    image::Rgba::<u8>([255u8, 0u8, 0u8, 255u8]),
                    self.medial_points[self.triangles[i][0] as usize][0] as i32,
                    self.medial_points[self.triangles[i][0] as usize][1] as i32,
                );
            
                // Drawing the circle associated with the medial axis point and triangle
                let rad = ((self.medial_points[ma_pt_i][0]-self.boundary_points[self.triangles[i][1] as usize][0])*(self.medial_points[ma_pt_i][0]-self.boundary_points[self.triangles[i][1] as usize][0])
                         + (self.medial_points[ma_pt_i][1]-self.boundary_points[self.triangles[i][1] as usize][1])*(self.medial_points[ma_pt_i][1]-self.boundary_points[self.triangles[i][1] as usize][1])).sqrt();
                imageproc::drawing::draw_hollow_circle_mut(
                    &mut img_rgb,
                    (self.medial_points[ma_pt_i][0] as i32, self.medial_points[ma_pt_i][1] as i32),
                    rad as i32,
                    image::Rgba::<u8>([255u8, 255u8, 255u8, 255u8]),
                );
            }
        }
        
        // Drawing the boundary points associated with the medial axis point
        for i in 0..self.index_bndry_points[ma_pt_i].len() {
            imageproc::drawing::draw_cross_mut(
                &mut img_rgb,
                image::Rgba::<u8>([0u8, 255u8, 0u8, 255u8]),
                self.boundary_points[self.index_bndry_points[ma_pt_i][i] as usize][0] as i32,
                self.boundary_points[self.index_bndry_points[ma_pt_i][i] as usize][1] as i32,
            );
        }
        // Save the image
        img_rgb.save(new_image).unwrap();
    }


    // To plot values as colormap
    pub fn plot_colormap(&self, values: &Vec<f32>, original_image:&str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Change the interval of the wedf values
        let (max_value, ind_max) = functions_vec::max_f32(values);
        let (min_value, ind_min) = functions_vec::min_f32(values);
        let mut new_values = values.clone();
        for i in 0..new_values.len() {
            new_values[i] = (new_values[i]-min_value)/(max_value-min_value);
        }

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        // Drawing the medial axis points
        for i in 0..self.medial_points.len() {
            let value = (new_values[i]*1020.0) as i32;
            let mut color = image::Rgba::<u8>([255, 255, 255, 255]);
            if value <= 255 {
                color = image::Rgba::<u8>([0, value as u8, 255, 255]);
            } else if value <= 510 {
                color = image::Rgba::<u8>([0, 255, (510-value) as u8, 255]);
            } else if value <= 765 {
                color = image::Rgba::<u8>([(value-510) as u8, 255, 0, 255]);
            } else if value <= 1020 {
                color = image::Rgba::<u8>([255, (1020-value) as u8, 0, 255]);
            }
            imageproc::drawing::draw_cross_mut(
                &mut img_rgb,
                color,
                self.medial_points[i][0] as i32,
                self.medial_points[i][1] as i32,
            );
        }
        
        // Save the image
        img_rgb.save(new_image).unwrap();
    }


    // To plot values as colormap to compare with the matlab version
    pub fn plot_colormap_matlab(&self, values: &Vec<f32>, original_image:&str, new_image: &str) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Change the interval of the wedf values
        let (max_value, ind_max) = functions_vec::max_f32(values);
        let (min_value, ind_min) = functions_vec::min_f32(values);
        let mut new_values = values.clone();
        for i in 0..new_values.len() {
            new_values[i] = (new_values[i]-min_value)/(max_value-min_value);
        }
        
        // Drawing the medial axis points
        for i in 0..self.medial_points.len() {
            let value = (new_values[i]*510.0) as i32;
            let mut color = image::Rgba::<u8>([0, 0, 0, 255]);
            if value <= 255 {
                color = image::Rgba::<u8>([value as u8, 0, 255, 255]);
            } else if value <= 510 {
                color = image::Rgba::<u8>([255, 0, (510-value) as u8, 255]);
            }
            for j in 0..self.triangles.len() {
                if self.triangles[j][0] as usize == i {
                    // Points of the triangle
                    let mut poly: Vec<Point<i32>> = Vec::new();
                    for k in 1..4 {
                        let x = self.boundary_points[self.triangles[j][k] as usize][0];
                        let y = self.boundary_points[self.triangles[j][k] as usize][1];
                        poly.push(Point::new(x as i32, y as i32));
                    }
                    if !(poly[0] == poly[2]) {
                        // Drawing the polygon
                        imageproc::drawing::draw_polygon_mut(
                            &mut img_rgb,
                            &poly,
                            color,
                        );
                    }
                }
            }
        }
        
        // Save the image
        img_rgb.save(new_image).unwrap();
    }


    // To plot the WEDF values of an object
    pub fn plot_wedf(&self, original_image:&str, new_image: &str) {
        self.plot_colormap(&self.wedf_array, original_image, new_image)
    }


    // To plot the WEDF values of an object to compare with the matlab version
    pub fn plot_wedf_matlab(&self, original_image:&str, new_image: &str) {
        self.plot_colormap_matlab(&self.wedf_array, original_image, new_image)
    }


    // To plot the EDF values of an object
    pub fn plot_edf(&self, original_image:&str, new_image: &str) {
        self.plot_colormap(&self.edf_array, original_image, new_image)
    }


    // To plot the EDF values of an object to compare with the matlab version
    pub fn plot_edf_matlab(&self, original_image:&str, new_image: &str) {
        self.plot_colormap_matlab(&self.edf_array, original_image, new_image)
    }


    // To plot the erosion thickness of an object
    pub fn plot_erosion_thickness(&self, original_image:&str, new_image: &str) {
        self.plot_colormap(&self.erosion_thickness, original_image, new_image)
    }


    // To plot the erosion thickness of an object to compare with the matlab version
    pub fn plot_erosion_thickness_matlab(&self, original_image:&str, new_image: &str) {
        self.plot_colormap_matlab(&self.erosion_thickness, original_image, new_image)
    }


    // To plot the shape tubularity of an object
    pub fn plot_shape_tubularity(&self, original_image:&str, new_image: &str) {
        self.plot_colormap(&self.shape_tubularity, original_image, new_image)
    }


    // To plot the shape tubularity of an object to compare with the matlab version
    pub fn plot_shape_tubularity_matlab(&self, original_image:&str, new_image: &str) {
        self.plot_colormap_matlab(&self.shape_tubularity, original_image, new_image)
    }


    // To plot the hierarchy parts of an object
    pub fn plot_hierarchy(&self, original_image:&str, new_image: &str, hierlabel: &Vec<i32>) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);
        
        // Clusters
        let clusters = functions_vec::unique_list(&hierlabel);

        // Colors
        let mut colors: Vec<image::Rgba::<u8>> = Vec::new();
        colors.push(image::Rgba::<u8>([255, 0, 0, 255]));
        if clusters.len() > 1 {
            let step = 1020/(clusters.len()-1) as i32;
            for i in 1..clusters.len() {
                let color = step*i as i32;
                if color <= 255 {
                    colors.push(image::Rgba::<u8>([255, color as u8, 0, 255]));
                } else if color <= 510 {
                    colors.push(image::Rgba::<u8>([(510-color) as u8, 255, 0, 255]));
                } else if color <= 765 {
                    colors.push(image::Rgba::<u8>([0, 255, (color-510) as u8, 255]));
                } else if color <= 1020 {
                    colors.push(image::Rgba::<u8>([0, (1020-color) as u8, 255, 255]));
                }
            }
        }
        
        for i in 0..self.triangles.len() {
            // Points of the triangle
            let mut poly: Vec<Point<i32>> = Vec::new();
            for j in 1..4 {
                let x = self.boundary_points[self.triangles[i][j] as usize][0];
                let y = self.boundary_points[self.triangles[i][j] as usize][1];
                poly.push(Point::new(x as i32, y as i32));
            }

            // Find the number of the cluster
            let num_color = functions_vec::find_vec(&clusters, hierlabel[self.triangles[i][0] as usize]);

            if !(poly[0] == poly[2]) {
                // Drawing the polygon
                imageproc::drawing::draw_polygon_mut(
                    &mut img_rgb,
                    &poly,
                    colors[num_color[0] as usize],
                );
            }

            /*
            imageproc::drawing::draw_cross_mut(
                &mut img_rgb,
                colors[num_color[0] as usize],
                self.medial_points[self.triangles[i][0] as usize][0] as i32,
                self.medial_points[self.triangles[i][0] as usize][1] as i32,
            );
            */
        }
        
        // Save the image
        img_rgb.save(new_image).unwrap();
    }

    
    // To plot the wedf hierarchy
    pub fn plot_wedf_hierarchy(&self, original_image: &str, new_image: &str) {
        self.plot_hierarchy(original_image, new_image, &self.clusters_wedf);
    }
    

    // To plot the new levels of shape details hierarchy
    pub fn plot_new_cuts(&self, original_image: &str, new_image: &str) {
        self.plot_hierarchy(original_image, new_image, &self.hierarchy_level);
    }

    
    // To plot the clusters of trunks of an object
    pub fn plot_clusters_trunk(&self, original_image: &str, new_image: &str, labels: &Vec<usize>, level: i32) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Clusters
        let clusters = functions_vec::unique_list(&labels);

        // Colors
        let mut colors: Vec<image::Rgba::<u8>> = Vec::new();
        colors.push(image::Rgba::<u8>([255, 0, 0, 255]));
        if clusters.len() > 1 {
            let step = 1020/(clusters.len()-1) as i32;
            for i in 1..clusters.len() {
                let color = step*i as i32;
                if color <= 255 {
                    colors.push(image::Rgba::<u8>([255, color as u8, 0, 255]));
                } else if color <= 510 {
                    colors.push(image::Rgba::<u8>([(510-color) as u8, 255, 0, 255]));
                } else if color <= 765 {
                    colors.push(image::Rgba::<u8>([0, 255, (color-510) as u8, 255]));
                } else if color <= 1020 {
                    colors.push(image::Rgba::<u8>([0, (1020-color) as u8, 255, 255]));
                }
            }
        }

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        for i in 0..self.trunk.len() {
            if self.trunk[i].level == level {
                for j in 0..self.triangles.len() {
                    if self.trunk[i].indices.contains(&(self.triangles[j][0] as usize)) {
                        // Points of the triangle
                        let mut poly: Vec<Point<i32>> = Vec::new();
                        for k in 1..4 {
                            let x = self.boundary_points[self.triangles[j][k] as usize][0];
                            let y = self.boundary_points[self.triangles[j][k] as usize][1];
                            poly.push(Point::new(x as i32, y as i32));
                        }
                        if !(poly[0] == poly[2]) {
                            // Drawing the polygon
                            imageproc::drawing::draw_polygon_mut(
                                &mut img_rgb,
                                &poly,
                                colors[labels[i]],
                            );
                        }
                    }
                }
            }
        }
        
        // Save the image
        img_rgb.save(new_image).unwrap();
    }

    
    // To plot the rainbows of clusters of trunks of an object
    pub fn plot_all_clusters_trunk(&self, original_image: &str, new_image: &str, labels: &Vec<usize>) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Clusters
        let clusters = functions_vec::unique_list(&labels);

        // Colors
        let mut colors: Vec<image::Rgba::<u8>> = Vec::new();
        if clusters.len() > 5 {
            colors.push(image::Rgba::<u8>([255, 0, 0, 1]));
            let step = 1020/(clusters.len()-1) as i32;
            for i in 1..clusters.len() {
                let color = step*i as i32;
                if color <= 255 {
                    colors.push(image::Rgba::<u8>([255, color as u8, 0, 1]));
                } else if color <= 510 {
                    colors.push(image::Rgba::<u8>([(510-color) as u8, 255, 0, 1]));
                } else if color <= 765 {
                    colors.push(image::Rgba::<u8>([0, 255, (color-510) as u8, 1]));
                } else if color <= 1020 {
                    colors.push(image::Rgba::<u8>([0, (1020-color) as u8, 255, 1]));
                }
            }
        } else {
            colors.push(image::Rgba::<u8>([255, 100, 100, 1]));
            colors.push(image::Rgba::<u8>([255, 200, 0, 1]));
            colors.push(image::Rgba::<u8>([150, 255, 150, 1]));
            colors.push(image::Rgba::<u8>([200, 150, 255, 1]));
            colors.push(image::Rgba::<u8>([100, 180, 255, 1]));
        }

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        // Find the maximum level of trunk
        let mut levels: Vec<i32> = Vec::new();
        for i in 0..self.trunk.len() {
            levels.push(self.trunk[i].level);
        }
        let levels = functions_vec::unique_list(&levels);

        for level in &levels {
            for i in 0..self.trunk.len() {
                if self.trunk[i].level == *level {
                    let num_color = functions_vec::find_vec(&clusters, labels[i]);
                    for j in 0..self.triangles.len() {
                        if self.trunk[i].indices.contains(&(self.triangles[j][0] as usize)) {
                            // Points of the triangle
                            let mut poly: Vec<Point<i32>> = Vec::new();
                            for k in 1..4 {
                                let x = self.boundary_points[self.triangles[j][k] as usize][0];
                                let y = self.boundary_points[self.triangles[j][k] as usize][1];
                                poly.push(Point::new(x as i32, y as i32));
                            }
                            if !(poly[0] == poly[2]) {
                                // Drawing the polygon
                                imageproc::drawing::draw_polygon_mut(
                                    &mut img_rgb,
                                    &poly,
                                    colors[num_color[0]],
                                );
                            }
                        }
                    }
                }
            }
        }

        // Save the image
        img_rgb.save(new_image).unwrap();
    }

    
    // To plot all the trunks of one cluster
    pub fn plot_one_cluster_trunk(&self, original_image: &str, new_image: &str, labels: &Vec<usize>, label: usize) {
        // Give the same dimensions as the original image
        let img_rgb_init = image::open(original_image).unwrap();
        let mut img_rgb = image::RgbImage::new(img_rgb_init.width(), img_rgb_init.height());
        for p in img_rgb.pixels_mut() {
            p[0] = 255;
            p[1] = 255;
            p[2] = 255;
        }
        let mut img_rgb = image::DynamicImage::ImageRgb8(img_rgb);

        // Clusters
        let clusters = functions_vec::unique_list(&labels);

        // Colors
        let mut colors: Vec<image::Rgba::<u8>> = Vec::new();
        if clusters.len() > 5 {
            colors.push(image::Rgba::<u8>([255, 0, 0, 1]));
            let step = 1020/(clusters.len()-1) as i32;
            for i in 1..clusters.len() {
                let color = step*i as i32;
                if color <= 255 {
                    colors.push(image::Rgba::<u8>([255, color as u8, 0, 1]));
                } else if color <= 510 {
                    colors.push(image::Rgba::<u8>([(510-color) as u8, 255, 0, 1]));
                } else if color <= 765 {
                    colors.push(image::Rgba::<u8>([0, 255, (color-510) as u8, 1]));
                } else if color <= 1020 {
                    colors.push(image::Rgba::<u8>([0, (1020-color) as u8, 255, 1]));
                }
            }
        } else {
            colors.push(image::Rgba::<u8>([255, 100, 100, 1]));
            colors.push(image::Rgba::<u8>([255, 200, 0, 1]));
            colors.push(image::Rgba::<u8>([150, 255, 150, 1]));
            colors.push(image::Rgba::<u8>([200, 150, 255, 1]));
            colors.push(image::Rgba::<u8>([100, 180, 255, 1]));
        }

        // Drawing the boundary edges
        for i in 0..self.boundary_edges.len() {
            imageproc::drawing::draw_line_segment_mut(
                &mut img_rgb,
                (self.boundary_points[self.boundary_edges[i][0] as usize][0], self.boundary_points[self.boundary_edges[i][0] as usize][1]),
                (self.boundary_points[self.boundary_edges[i][1] as usize][0], self.boundary_points[self.boundary_edges[i][1] as usize][1]),
                image::Rgba::<u8>([0u8, 0u8, 0u8, 255u8]),
            );
        }
        
        for i in 0..self.trunk.len() {
            if labels[i] == clusters[label] {
                let num_color = functions_vec::find_vec(&clusters, labels[i]);
                for j in 0..self.triangles.len() {
                    if self.trunk[i].indices.contains(&(self.triangles[j][0] as usize)) {
                        // Points of the triangle
                        let mut poly: Vec<Point<i32>> = Vec::new();
                        for k in 1..4 {
                            let x = self.boundary_points[self.triangles[j][k] as usize][0];
                            let y = self.boundary_points[self.triangles[j][k] as usize][1];
                            poly.push(Point::new(x as i32, y as i32));
                        }
                        if !(poly[0] == poly[2]) {
                            // Drawing the polygon
                            imageproc::drawing::draw_polygon_mut(
                                &mut img_rgb,
                                &poly,
                                colors[num_color[0]],
                            );
                        }
                    }
                }
            }
        }

        // Save the image
        img_rgb.save(new_image).unwrap();
    }


} 