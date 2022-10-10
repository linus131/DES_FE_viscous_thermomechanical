use crate::sparse_matrix::SparseMatrix;
use crate::primitives::Node;
use std::collections::HashMap;

pub(crate) struct StructuralModel{
    pub(crate) nodes: Vec<[f64;3]>,
    cells: Vec<[usize;8]>,
    to_bc_dof:Vec<usize>,
    props: [f64;3],
    loads_bcs: Vec<[f64;24]>,
    //loadsg_bcs: Vec<f64>,
    bcs: Vec<bool>,
    pub(crate) uold: Vec<f64>,
    pub(crate) uold_bcs: Vec<f64,>,
    dN:Vec<Vec<Vec<f64>>>,
    N:Vec<Vec<f64>>,
    Kloc:Vec<Vec<Vec<f64>>>,
    nbcs:usize,
    global_node_to_model_node: HashMap<usize, usize>,
    current_nbc:usize,
}

impl StructuralModel {

    pub(crate) fn new(num_nodes: usize, num_cells:usize, nbcs:usize,props:[f64;3])->StructuralModel{
        let mut stmdl = StructuralModel::new_empty(num_nodes,num_cells,nbcs,props);
        stmdl.dN = StructuralModel::generate_dN();
        stmdl.N = StructuralModel::generate_N();
//        let mut to_bc_dof = Vec::with_capacity(num_nodes*3 - nbcs);

      //  let mut count = 0;
        return stmdl
    }

    pub(crate) fn new_empty(num_nodes: usize, num_cells:usize, nbcs:usize, props:[f64;3]) ->StructuralModel{
        return StructuralModel{
            nodes: Vec::with_capacity(num_nodes),
            cells: Vec::with_capacity(num_cells),
            to_bc_dof:Vec::with_capacity(num_nodes*3 - nbcs),
            props ,
            loads_bcs: Vec::with_capacity(num_nodes*3),
            //loadsf_bcs: Vec::with_capacity(num_nodes*3 - nbcs),
            bcs: Vec::with_capacity(num_nodes*3),
            uold: Vec::with_capacity(num_nodes*3),
            uold_bcs: Vec::with_capacity(num_nodes*3 - nbcs),
            dN:Vec::with_capacity(8),
            N:Vec::with_capacity(8),
            Kloc:Vec::with_capacity(num_cells),
            nbcs,
            global_node_to_model_node: HashMap::with_capacity(num_nodes),
            current_nbc:0,
        }
    }



    /// generates derivative of shape function at 2x2x2 integration points
    fn generate_dN()->Vec<Vec<Vec<f64>>>{
        //let dN = Vec::with_capacity(8);
        let dn1 = vec![
            vec![-0.311,0.311,0.083333,-0.083333,-0.083333,0.083333,0.022329,-0.022329],
            vec![-0.311,-0.083333,0.083333,0.311,-0.083333,-0.022329,0.022329,0.083333],
            vec![-0.311,-0.083333,-0.022329,-0.083333,0.311,0.083333,0.022329,0.083333]];

        let dn2 = vec![
            vec![-0.083333,0.083333,0.022329,-0.022329,-0.311,0.311,0.083333,-0.083333],
            vec![-0.083333,-0.022329,0.022329,0.083333,-0.311,-0.083333,0.083333,0.311],
            vec![-0.311,-0.083333,-0.022329,-0.083333,0.311,0.083333,0.022329,0.083333]];

        let dn3= vec![
            vec![-0.083333,0.083333,0.311,-0.311,-0.022329,0.022329,0.083333,-0.083333],
            vec![-0.311,-0.083333,0.083333,0.311,-0.083333,-0.022329,0.022329,0.083333],
            vec![-0.083333,-0.022329,-0.083333,-0.311,0.083333,0.022329,0.083333,0.311]];

        let dn4 = vec![
            vec![-0.022329,0.022329,0.083333,-0.083333,-0.083333,0.083333,0.311,-0.311],
            vec![-0.083333,-0.022329,0.022329,0.083333,-0.311,-0.083333,0.083333,0.311],
            vec![-0.083333,-0.022329,-0.083333,-0.311,0.083333,0.022329,0.083333,0.311]];

        let dn5 = vec![
            vec![-0.311,0.311,0.083333,-0.083333,-0.083333,0.083333,0.022329,-0.022329],
            vec![-0.083333,-0.311,0.311,0.083333,-0.022329,-0.083333,0.083333,0.022329],
            vec![-0.083333,-0.311,-0.083333,-0.022329,0.083333,0.311,0.083333,0.022329]];

        let dn6 = vec![
            vec![-0.083333,0.083333,0.022329,-0.022329,-0.311,0.311,0.083333,-0.083333],
            vec![-0.022329,-0.083333,0.083333,0.022329,-0.083333,-0.311,0.311,0.083333],
            vec![-0.083333,-0.311,-0.083333,-0.022329,0.083333,0.311,0.083333,0.022329]];

        let dn7 = vec![
            vec![-0.083333,0.083333,0.311,-0.311,-0.022329,0.022329,0.083333,-0.083333],
            vec![-0.083333,-0.311,0.311,0.083333,-0.022329,-0.083333,0.083333,0.022329],
            vec![-0.022329,-0.083333,-0.311,-0.083333,0.022329,0.083333,0.311,0.083333]];

        let dn8 = vec![
            vec![-0.022329,0.022329,0.083333,-0.083333,-0.083333,0.083333,0.311,-0.311],
            vec![-0.022329,-0.083333,0.083333,0.022329,-0.083333,-0.311,0.311,0.083333],
            vec![-0.022329,-0.083333,-0.311,-0.083333,0.022329,0.083333,0.311,0.083333]];

        return vec![dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8];
    }

    /// generates shape function at 2x2x2 integration points
    fn generate_N()->Vec<Vec<f64>>{
        //let N = Vec::with_capacity(8);
        let N = vec![vec![0.49056,0.13145,0.035221,0.13145,0.13145,0.035221,0.0094374,0.035221],
                     vec![0.13145,0.035221,0.0094374,0.035221,0.49056,0.13145,0.035221,0.13145],
                     vec![0.13145,0.035221,0.13145,0.49056,0.035221,0.0094374,0.035221,0.13145],
                     vec![0.035221,0.0094374,0.035221,0.13145,0.13145,0.035221,0.13145,0.49056],
                     vec![0.13145,0.49056,0.13145,0.035221,0.035221,0.13145,0.035221,0.0094374],
                     vec![0.035221,0.13145,0.035221,0.0094374,0.13145,0.49056,0.13145,0.035221],
                     vec![0.035221,0.13145,0.49056,0.13145,0.0094374,0.035221,0.13145,0.035221],
                     vec![0.0094374,0.035221,0.13145,0.035221,0.035221,0.13145,0.49056,0.13145]];
        return N;
    }

    /// generates strain_displacement_matrix at 2x2x2 integration points
    fn generate_B(&self, dx:f64, dy:f64, dz:f64)->Vec<Vec<Vec<f64>>> {
        let mut Bmat = Vec::with_capacity(8);
        let matJac = [dx/2.0, dy/2.0, dz/2.0];
        let mut derivMat = self.dN.clone();
        for i in 0..derivMat.len() {
            for j in 0..3 {
                for k in 0..8 {
                    derivMat[i][j][k] = derivMat[i][j][k] / matJac[j];
                }
            }
        }
        for i in 0..8 {

                let mut B = SparseMatrix::create_zero_dense_matrix(6, 24);
                Bmat.push(B);

            for j in 0..8 {
                Bmat[i][0] [j * 3] = derivMat[i][0][j];
                Bmat[i][1] [j * 3 + 1] = derivMat[i][1][j];
                Bmat[i][2] [j * 3 + 2] = derivMat[i][2][j];

                Bmat[i][3][j * 3] = derivMat[i][1][j];
                Bmat[i][3] [j * 3 + 1] = derivMat[i][0] [j];

                Bmat[i][4] [j * 3 + 2] = derivMat[i][1][j];
                Bmat[i][4] [j * 3 + 1] = derivMat[i][2][j];

                Bmat[i][5] [j * 3] = derivMat[i][2][j];
                Bmat[i][5] [j * 3 + 2] = derivMat[i][0][j];
            }
            //Bmat.push(B);
        }
        return Bmat;
    }

    /// generates local stiffness matrix and stores the values in the model
    fn generate_local_stiffness_matrix(&mut self){
        let E = self.props[0];
        let nu = self.props[1];
        let density = self.props[2];

        let mut C = vec![
            vec![1_f64 - nu, nu, nu, 0., 0., 0.],
            vec![nu, 1. - nu, nu, 0., 0., 0.],
            vec![nu, nu, 1. - nu, 0., 0., 0.],
            vec![0., 0., 0., (1. - 2. * nu) / 2., 0., 0.],
            vec![0., 0., 0., 0., (1. - 2. * nu) / 2., 0.],
            vec![0., 0., 0., 0., 0., (1. - 2. * nu) / 2.],
        ];

        for i in 0..C.len(){
            for j in 0..6{
                C[i][j] = C[i][j]*E/(1.0+nu)/(1.0-2.0*nu);
            }
        }

        let a = 0.577350269189626_f64;
        let Wa = 1.0_f64;
        let quadraturevals = vec![[-a, -a, -a, Wa],
                                  [-a, -a, a, Wa],
                                  [-a, a, -a, Wa],
                                  [-a, a, a, Wa],
                                  [a, -a, -a, Wa],
                                  [a, -a, a, Wa],
                                  [a, a, -a, Wa],
                                  [a, a, a, Wa]
        ];


        // assemble matrices

        let mut nodal_coords = SparseMatrix::create_zero_dense_matrix(8, 3);
        for ii in 0..8 {
            for jj in 0..3 {
                nodal_coords[ii][jj] = self.nodes[self.cells[0][ii]][jj];
            }
        }
        let dx = f64::max(nodal_coords[7][0],nodal_coords[6][0]) - f64::min(nodal_coords[0][0],nodal_coords[1][0]);
        let dy = f64::max(nodal_coords[7][1],nodal_coords[6][1]) - f64::min(nodal_coords[0][1],nodal_coords[1][1]);
        let dz = f64::max(nodal_coords[7][2],nodal_coords[6][2]) - f64::min(nodal_coords[0][2],nodal_coords[1][2]);
        let vol = (dx*dy*dz).abs();
        let matJac = [dx,dy,dz];
        let B = self.generate_B(dx,dy,dz);


        for j in 0..self.cells.len() {
            // create local stiffness matrix
            let mut Kloc = SparseMatrix::create_zero_dense_matrix(24, 24);
            for kk in 0..8 {
                let m1 = SparseMatrix::dense_mat_mult_easy_transpose(&B[kk], &C);
                //let Bj = B[j][0][0];
                //let Bk = Bj+1.0;
                let mut m2 = SparseMatrix::dense_mat_mult_easy(&m1, &B[kk]);
                SparseMatrix::multiply_dense_by_constant(&mut m2, quadraturevals[kk][3] * vol/8.0);
                SparseMatrix::add_dense_matrix(&mut Kloc, &m2);
            }
            self.Kloc.push(Kloc);
        }


    }

    /// generates local load matrix, assembles and outputs global load matrix
    pub(crate) fn generate_load_vec_bc(&self, tempchange:Vec<f64>, dx:f64,dy:f64,dz:f64,cte:f64)->Vec<f64>{
        let mut loads= Vec::with_capacity(self.nodes.len()*3-self.nbcs);
        let mut loads_g = vec![0.0_f64;self.nodes.len()*3];
        let B = self.generate_B(dx,dy,dz);

        for i in 0..self.cells.len(){
            let load_cell = StructuralModel::generate_local_load_matrix_single_cell(&self.props,&B,tempchange[i],dx,dy,dz,1e-6);
           // println!("{:?}",load_cell);
            let map_nodes_to_dof = [
                self.cells[i][0]*3,self.cells[i][0]*3+1, self.cells[i][0]*3+2,
                self.cells[i][1]*3,self.cells[i][1]*3+1, self.cells[i][1]*3+2,
                self.cells[i][2]*3,self.cells[i][2]*3+1, self.cells[i][2]*3+2,
                self.cells[i][3]*3,self.cells[i][3]*3+1, self.cells[i][3]*3+2,
                self.cells[i][4]*3,self.cells[i][4]*3+1, self.cells[i][4]*3+2,
                self.cells[i][5]*3,self.cells[i][5]*3+1, self.cells[i][5]*3+2,
                self.cells[i][6]*3,self.cells[i][6]*3+1, self.cells[i][6]*3+2,
                self.cells[i][7]*3,self.cells[i][7]*3+1, self.cells[i][7]*3+2,
            ];

             for j in 0..24{
                 loads_g[map_nodes_to_dof[j]] = loads_g[map_nodes_to_dof[j]]+load_cell[j];
             }

        }

        //println!("{:?}",loads_g);
        for i in 0..loads_g.len(){
            if !self.bcs[i]{
                loads.push(loads_g[i]);
            }
        }
        return loads;
    }

    /// generate load vector for single cell due to temperature change
    fn generate_local_load_matrix_single_cell( props:&[f64;3], B:&Vec<Vec<Vec<f64>>>, tempchange:f64, dx:f64, dy:f64, dz:f64, cte:f64) -> [f64;24] {
        //let mut loads = Vec::with_capacity(self.cells.len());
        let E = props[0];
        let nu = props[1];
        let density = props[2];

        let mut C = vec![
            vec![1_f64 - nu, nu, nu, 0., 0., 0.],
            vec![nu, 1. - nu, nu, 0., 0., 0.],
            vec![nu, nu, 1. - nu, 0., 0., 0.],
            vec![0., 0., 0., (1. - 2. * nu) / 2., 0., 0.],
            vec![0., 0., 0., 0., (1. - 2. * nu) / 2., 0.],
            vec![0., 0., 0., 0., 0., (1. - 2. * nu) / 2.],
        ];

        for i in 0..C.len(){
            for j in 0..6{
                C[i][j] = C[i][j]*E/(1.0+nu)/(1.0-2.0*nu);
            }
        }

        let a = 0.577350269189626_f64;
        let Wa = 1.0_f64;
        let quadraturevals = vec![[-a, -a, -a, Wa],
                                  [-a, -a, a, Wa],
                                  [-a, a, -a, Wa],
                                  [-a, a, a, Wa],
                                  [a, -a, -a, Wa],
                                  [a, -a, a, Wa],
                                  [a, a, -a, Wa],
                                  [a, a, a, Wa]
        ];

        //let B = self.generate_B(dx,dy,dz);
        let vol = dx*dy*dz;
        let thermal_strains = vec![cte*tempchange, cte*tempchange, cte*tempchange, 0., 0., 0.];

        let mut loads_loc = [0.0_f64;24];
        for kk in 0..8 {
            let m1 = SparseMatrix::dense_mat_mult_easy_transpose(&B[kk], &C);
            //let z = m1.len();
            let mut m2 = SparseMatrix::dense_mat_mult_vec_easy(&m1, &thermal_strains);
            for jj in 0..m2.len(){
                loads_loc[jj] = loads_loc[jj]+ m2[jj]*vol/8.0_f64;
            }
        }

        return loads_loc;

    }


    /// generates global stiffness matrix with boundary conditions applied
    fn generate_global_stiffness_matrix_with_bc(&self,num_bcs:usize)-> SparseMatrix {
        let mut Kg_bc = SparseMatrix::with_capacity([self.nodes.len()*3-num_bcs, self.nodes.len()*3-num_bcs], (self.nodes.len()*3-num_bcs)*24);

        for i in 0..self.cells.len(){
            let map_nodes_to_dof = [
                self.cells[i][0]*3,self.cells[i][0]*3+1, self.cells[i][0]*3+2,
                self.cells[i][1]*3,self.cells[i][1]*3+1, self.cells[i][1]*3+2,
                self.cells[i][2]*3,self.cells[i][2]*3+1, self.cells[i][2]*3+2,
                self.cells[i][3]*3,self.cells[i][3]*3+1, self.cells[i][3]*3+2,
                self.cells[i][4]*3,self.cells[i][4]*3+1, self.cells[i][4]*3+2,
                self.cells[i][5]*3,self.cells[i][5]*3+1, self.cells[i][5]*3+2,
                self.cells[i][6]*3,self.cells[i][6]*3+1, self.cells[i][6]*3+2,
                self.cells[i][7]*3,self.cells[i][7]*3+1, self.cells[i][7]*3+2,
            ];
            for j in 0..24{
                if !self.bcs[map_nodes_to_dof[j]] {
                    for k in 0..24 {
                        if !self.bcs[map_nodes_to_dof[k]] {
                            // ensure symmetry
                            if k>j {
                                let nn = self.to_bc_dof[map_nodes_to_dof[j]];
                                let oo = self.to_bc_dof[map_nodes_to_dof[k]];
                                Kg_bc.insert_add(self.to_bc_dof[map_nodes_to_dof[j]], self.to_bc_dof[map_nodes_to_dof[k]], self.Kloc[i][j][k]);
                                Kg_bc.insert_add(self.to_bc_dof[map_nodes_to_dof[k]], self.to_bc_dof[map_nodes_to_dof[j]], self.Kloc[i][j][k]);
                            }
                            if k == j{
                                let val = self.Kloc[i][j][k];
                                let zz = map_nodes_to_dof[k];
                                Kg_bc.insert_add(self.to_bc_dof[map_nodes_to_dof[k]], self.to_bc_dof[map_nodes_to_dof[j]], self.Kloc[i][j][k]);
                            }
                        }
                    }
                }
            }
        }
        return Kg_bc;
    }


    /// solves Ku=F using conjugate gradient method and stores result in the model
    pub(crate) fn structural_analysis_solve(&mut self,loads_bcs:&Vec<f64>){
        self.generate_local_stiffness_matrix();
        let Kg_bc = self.generate_global_stiffness_matrix_with_bc(self.nodes.len()*3-self.current_nbc);
        //println!("inside sas");
        //println!("{:?}",loads_bcs);
        let mut uold_bcs = Vec::with_capacity(self.uold_bcs.len());
        for i in 0..self.uold_bcs.len(){
            //if !self.bcs[i] {
                uold_bcs.push( self.uold_bcs[i]);
           //% }
        }

        let soln = Kg_bc.solve_cg(loads_bcs, &uold_bcs);
        let mut count = 0;
        for i in 0.. self.uold_bcs.len(){
               self.uold_bcs[i] = soln[i] + self.uold_bcs[i];
             }


    }

    /// adds nodes and cell to the model. Needs global nodelist and connectivity of the cell based
    /// on global nodelist
    pub(crate) fn add_cell(&mut self, global_nodelist:&Vec<Node>, cell: [usize;8],  global_bcs:&mut Vec<[bool;3]>, bcs: [[bool;3];8]){
        //check if the node already exists in the model use hashtable
        //if node exists use that node number to define new cell (element)
        //if node does not exist add new node to model and define cell (element)
        let mut global_to_model_cell = Vec::with_capacity(8);
        for i in 0..8 {
            if self.global_node_to_model_node.contains_key(&cell[i]){
                let nd = self.global_node_to_model_node.get(&cell[i]).expect("no such key exists");
                global_bcs[*nd][0] = bcs[i][0];
                global_bcs[*nd][1] = bcs[i][1];
                global_bcs[*nd][2] = bcs[i][2];
                global_to_model_cell.push(*nd);
            }
            else{
                global_to_model_cell.push(self.global_node_to_model_node.len());
                self.global_node_to_model_node.insert(cell[i],self.global_node_to_model_node.len());
                global_bcs.push([bcs[i][0],bcs[i][1],bcs[i][2]]);
                let ndval = global_nodelist[cell[i]];
                self.nodes.push([ndval.pt.x, ndval.pt.y, ndval.pt.z]);
                //let pt1 = cell[i];
                for zz in 0..3 {
                    let gbci = global_bcs[cell[i]][zz];
                    self.bcs.push(gbci.clone());
                    //self.bcs.push(global_bcs[cell[i]][1]);
                    //self.bcs.push(global_bcs[cell[i]][2]);
                    if gbci {
                        self.to_bc_dof.push(0);

                    }
                    else {
                        self.to_bc_dof.push(self.current_nbc);
                        self.uold_bcs.push(0.0_f64);
                        self.current_nbc = self.current_nbc+1;
                    }
                }
            }
        }
        let mut g2m = [0;8];
        for i in 0..8{
            g2m[i] = global_to_model_cell[i];
        }
        self.cells.push(g2m);

    }

    pub(crate) fn add_cell2(&mut self, global_nodelist:&Vec<[f64;3]>, cell: [usize;8], global_bcs:&Vec<[bool;3]>){
        //check if the node already exists in the model use hashtable
        //if node exists use that node number to define new cell (element)
        //if node does not exist add new node to model and define cell (element)
        let mut global_to_model_cell = Vec::with_capacity(8);
        for i in 0..8 {
            if self.global_node_to_model_node.contains_key(&cell[i]){
                let nd = self.global_node_to_model_node.get(&cell[i]).expect("no such key exists");
                global_to_model_cell.push(*nd);
            }
            else{
                global_to_model_cell.push(self.global_node_to_model_node.len());
                self.global_node_to_model_node.insert(cell[i],self.global_node_to_model_node.len());
                let ndval = global_nodelist[cell[i]];
                self.nodes.push([ndval[0], ndval[1], ndval[2]]);
                //let mut to_dof = Vec::with_capacity(3);
                for zz in 0..3 {
                    let gbci = global_bcs[cell[i]][zz];
                    self.bcs.push(gbci.clone());
                    //self.bcs.push(global_bcs[cell[i]][1]);
                    //self.bcs.push(global_bcs[cell[i]][2]);
                    if gbci {
                        self.to_bc_dof.push(0);

                    }
                    else {
                        self.to_bc_dof.push(self.current_nbc);
                        self.uold_bcs.push(0.0_f64);
                        self.current_nbc = self.current_nbc+1;
                    }
                }
            }
        }
        let mut g2m = [0;8];
        for i in 0..8{
            g2m[i] = global_to_model_cell[i];
        }
        self.cells.push(g2m);
    }




}