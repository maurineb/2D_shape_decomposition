use nalgebra::base::*;

const _EPS: f32 = 0.00001;
const _MIN: f32 = f32::MIN_POSITIVE;
const _MAX: f32 = f32::MAX;

pub struct Polygon {
    // List of coordinates of the points in order
    points: Vec<Vector2<f32>>
}

impl Polygon {

    // Constructor
    pub fn new(pts: &Vec<Vector2<f32>>) -> Polygon {
        // Error if less than 3 points
        if pts.len() < 3 {
            panic!("We need to have at least three points in order to create a new Polygon.");
        }
        Polygon{points: pts.clone()}
    }
    

    // Returns the next vertex of a vertex
    fn next(&self, num_vertex: usize) -> usize {
        let next_vertex;
        if num_vertex >= self.points.len() {
            panic!("We can't find the following vertex of this number because it doesn't exist in this polygon.");
        }
        else if num_vertex == self.points.len()-1 {
            next_vertex = 0;
        }
        else {
            next_vertex = num_vertex+1;
        }
        next_vertex
    }
    

    // Returns the previous vertex of a vertex
    fn previous(&self, num_vertex: usize) -> usize {
        let previous_vertex;
        if num_vertex >= self.points.len() {
            panic!("We can't find the following vertex of this number because it doesn't exist in this polygon.");
        }
        else if num_vertex == 0 {
            previous_vertex = self.points.len()-1;
        }
        else {
            previous_vertex = num_vertex-1;
        }
        previous_vertex
    }


    // Returns true if the point is on the edges of the polygon
    fn point_on_boundary(&self, x: f32, y: f32) -> bool {
        let x_pt = x;
        let y_pt = y;
        let mut on = false;

        for i in 0..self.points.len() {
            // Extracts the coordinates of the points of the edge
            let x1 = self.points[i][0];
            let y1 = self.points[i][1];
            let x2 = self.points[self.next(i)][0];
            let y2 = self.points[self.next(i)][1];

            // Check if it's outside the box
            if (y_pt<y1.min(y2) || y_pt>y1.max(y2)) || (x_pt>x1.max(x2) || x_pt<x1.min(x2)) {
                on = false;
            } else if (x_pt==x1 && y_pt==y1) || (x_pt==x2 && y_pt==y2) {
                on = true;
                break;
            }
            // Inside the bounding box
            else {
                // Find the coefficient of the edge
                let mut coef_segment = _MAX;
                if (x1-x2).abs() > _MIN {
                    coef_segment = (y2-y1)/(x2-x1);
                }
                // Find the coefficient of the line passing by the first point and the tested point
                let mut coef_pt = _MAX;
                if (x1-x_pt).abs() > _MIN {
                    coef_pt = (y_pt-y1)/(x_pt-x1);
                }
                // Compare the values
                if coef_pt==coef_segment || coef_pt==-coef_segment {
                    on = true;
                    break;
                } else {
                    on = false;
                }
            }
        }
        on
    }


    // Returns true if the horizontal ray starting from the point P intersects the side (segment)
    fn ray_intersects_segment(&self, x: f32, y: f32, edg_num: usize) -> bool {
        let x_pt = x;
        let mut y_pt = y;

        // Extracts the coordinates of the points of the edge
        let mut x1 = self.points[edg_num][0];
        let mut y1 = self.points[edg_num][1];
        let mut x2 = self.points[self.next(edg_num)][0];
        let mut y2 = self.points[self.next(edg_num)][1];

        // If the first point is above the second point, invert the values
        if y1>y2 {
            let temp_x = x2;
            let temp_y = y2;
            x2 = x1;
            x1 = temp_x;
            y2 = y1;
            y1 = temp_y;
        }

        if y_pt == y1 || y_pt == y2 {
            false
        }
        else if (y_pt>y2 || y_pt<y1) || x_pt>x1.max(x2) {
            false
        }
        else if x_pt<x1.min(x2) {
            true
        }
        else {
            // Find the slope of the edge
            let mut coef_segment = _MAX;
            if (x1-x2).abs() > _MIN {
                coef_segment = (y2-y1)/(x2-x1);
            }
            // Find the slope of the line passing by the first point and the tested point
            let mut coef_pt = _MAX;
            if (x1-x_pt).abs() > _MIN {
                coef_pt = (y_pt-y1)/(x_pt-x1);
            }
            // Compare the values
            coef_pt >= coef_segment
        }
    }


    // Returns true if the horizontal ray starting from the point P intersects the vertex
    fn ray_intersects_vertex(&self, x: f32, y: f32, vertex_num: usize) -> bool {
        self.points[vertex_num][1] == y && x<=self.points[vertex_num][0]
    }


    // Returns true if the point is inside (or on) the polygon
    pub fn contains_point(&self, x: f32, y: f32) -> bool {
        let mut count = 0;
        let mut vertices_on_same_line = false;

        for i in 0..self.points.len() {
            vertices_on_same_line = self.points[self.previous(i)][1]==self.points[i][1];
            if self.ray_intersects_segment(x, y, i) {
                count += 1;
                vertices_on_same_line = false;
            } else if self.ray_intersects_vertex(x, y, i) && !vertices_on_same_line {
                let mut j = self.next(i);
                while self.points[j][1]==self.points[i][1] {j = self.next(j);}
                if !((self.points[self.previous(i)][1]>self.points[i][1]) && (self.points[self.next(j)][1]>self.points[i][1])) &&
                !((self.points[self.previous(i)][1]<self.points[i][1]) && (self.points[self.next(j)][1]<self.points[i][1])) {
                    count += 1;
                }
            }
        }
        count % 2 == 1 || self.point_on_boundary(x, y)
    }
}
