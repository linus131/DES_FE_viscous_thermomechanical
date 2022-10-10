use std::borrow::{Borrow, BorrowMut};
use std::time::Instant;
use std::collections::HashMap;
use std::fmt;
use rayon::{ThreadPool, ThreadPoolBuilder};
use crate::sparse_matrix_new::SparseMatrix;
use crate::primitives::{Node, Point};
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

///shape function for eight noded hexahedral elements
const N:[[f64;8];8] = [[0.490562612162344,0.131445855765802,0.0352208109008645,0.131445855765802,0.131445855765802,0.0352208109008645,0.00943738783765592,0.0352208109008645],
    [0.131445855765802,0.0352208109008645,0.00943738783765592,0.0352208109008645,0.490562612162344,0.131445855765802,0.0352208109008645,0.131445855765802],
    [0.131445855765802,0.0352208109008645,0.131445855765802,0.490562612162344,0.0352208109008645,0.00943738783765592,0.0352208109008645,0.131445855765802],
    [0.0352208109008645,0.00943738783765592,0.0352208109008645,0.131445855765802,0.131445855765802,0.0352208109008645,0.131445855765802,0.490562612162344],
    [0.131445855765802,0.490562612162344,0.131445855765802,0.0352208109008645,0.0352208109008645,0.131445855765802,0.0352208109008645,0.00943738783765592],
    [0.0352208109008645,0.131445855765802,0.0352208109008645,0.00943738783765592,0.131445855765802,0.490562612162344,0.131445855765802,0.0352208109008645],
    [0.0352208109008645,0.131445855765802,0.490562612162344,0.131445855765802,0.00943738783765592,0.0352208109008645,0.131445855765802,0.0352208109008645],
    [0.00943738783765592,0.0352208109008645,0.131445855765802,0.0352208109008645,0.0352208109008645,0.131445855765802,0.490562612162344,0.131445855765802]];

/// derivative of shape functions for eight noded hexahedral elements
const dn1:[[f64;8];3] = [
    [ -0.311004233964073,0.311004233964073,0.0833333333333333,-0.0833333333333333,-0.0833333333333333,0.0833333333333333,0.0223290993692602,-0.0223290993692602],
    [    -0.311004233964073,-0.0833333333333333,0.0833333333333333,0.311004233964073,-0.0833333333333333,-0.0223290993692602,0.0223290993692602,0.0833333333333333],
    [   -0.311004233964073,-0.0833333333333333,-0.0223290993692602,-0.0833333333333333,0.311004233964073,0.0833333333333333,0.0223290993692602,0.0833333333333333]];

const dn2:[[f64;8];3] = [
    [ -0.0833333333333333,0.0833333333333333,0.0223290993692602,-0.0223290993692602,-0.311004233964073,0.311004233964073,0.0833333333333333,-0.0833333333333333],
    [ -0.0833333333333333,-0.0223290993692602,0.0223290993692602,0.0833333333333333,-0.311004233964073,-0.0833333333333333,0.0833333333333333,0.311004233964073],
    [-0.311004233964073,-0.0833333333333333,-0.0223290993692602,-0.0833333333333333,0.311004233964073,0.0833333333333333,0.0223290993692602,0.0833333333333333]];

const dn3:[[f64;8];3] = [
    [-0.0833333333333333,0.0833333333333333,0.311004233964073,-0.311004233964073,-0.0223290993692602,0.0223290993692602,0.0833333333333333,-0.0833333333333333],
    [    -0.311004233964073,-0.0833333333333333,0.0833333333333333,0.311004233964073,-0.0833333333333333,-0.0223290993692602,0.0223290993692602,0.0833333333333333],
    [    -0.0833333333333333,-0.0223290993692602,-0.0833333333333333,-0.311004233964073,0.0833333333333333,0.0223290993692602,0.0833333333333333,0.311004233964073]];

const dn4:[[f64;8];3] = [
    [-0.0223290993692602,0.0223290993692602,0.0833333333333333,-0.0833333333333333,-0.0833333333333333,0.0833333333333333,0.311004233964073,-0.311004233964073],
    [-0.0833333333333333,-0.0223290993692602,0.0223290993692602,0.0833333333333333,-0.311004233964073,-0.0833333333333333,0.0833333333333333,0.311004233964073],
    [-0.0833333333333333,-0.0223290993692602,-0.0833333333333333,-0.311004233964073,0.0833333333333333,0.0223290993692602,0.0833333333333333,0.311004233964073]];

const dn5:[[f64;8];3] = [
    [-0.311004233964073,0.311004233964073,0.0833333333333333,-0.0833333333333333,-0.0833333333333333,0.0833333333333333,0.0223290993692602,-0.0223290993692602],
    [-0.0833333333333333,-0.311004233964073,0.311004233964073,0.0833333333333333,-0.0223290993692602,-0.0833333333333333,0.0833333333333333,0.0223290993692602],
    [-0.0833333333333333,-0.311004233964073,-0.0833333333333333,-0.0223290993692602,0.0833333333333333,0.311004233964073,0.0833333333333333,0.0223290993692602],
];

const dn6:[[f64;8];3] = [
    [-0.0833333333333333,0.0833333333333333,0.0223290993692602,-0.0223290993692602,-0.311004233964073,0.311004233964073,0.0833333333333333,-0.0833333333333333],
    [-0.0223290993692602,-0.0833333333333333,0.0833333333333333,0.0223290993692602,-0.0833333333333333,-0.311004233964073,0.311004233964073,0.0833333333333333],
    [-0.0833333333333333,-0.311004233964073,-0.0833333333333333,-0.0223290993692602,0.0833333333333333,0.311004233964073,0.0833333333333333,0.0223290993692602]];

const dn7:[[f64;8];3] = [
    [-0.0833333333333333,0.0833333333333333,0.311004233964073,-0.311004233964073,-0.0223290993692602,0.0223290993692602,0.0833333333333333,-0.0833333333333333],
    [   -0.0833333333333333,-0.311004233964073,0.311004233964073,0.0833333333333333,-0.0223290993692602,-0.0833333333333333,0.0833333333333333,0.0223290993692602],
    [  -0.0223290993692602,-0.0833333333333333,-0.311004233964073,-0.0833333333333333,0.0223290993692602,0.0833333333333333,0.311004233964073,0.0833333333333333]];

const dn8:[[f64;8];3] = [
    [-0.0223290993692602,0.0223290993692602,0.0833333333333333,-0.0833333333333333,-0.0833333333333333,0.0833333333333333,0.311004233964073,-0.311004233964073],
    [-0.0223290993692602,-0.0833333333333333,0.0833333333333333,0.0223290993692602,-0.0833333333333333,-0.311004233964073,0.311004233964073,0.0833333333333333],
    [-0.0223290993692602,-0.0833333333333333,-0.311004233964073,-0.0833333333333333,0.0223290993692602,0.0833333333333333,0.311004233964073,0.0833333333333333],
];
const dN:[[[f64;8];3];8] = [dn1,dn2,dn3,dn4,dn5,dn6,dn7,dn8];


pub(crate) struct StructuralModel{
    pub(crate) nodes:Vec<Vec<f64>>,
    pub(crate) elements:Vec<Vec<usize>>,
    pub(crate) loads:Vec<f64>,
    pub(crate) structural_bcs: Vec<bool>,
    pub(crate) bcs:Vec<bool>,
    pub(crate) u:Vec<f64>,
    pub(crate) uel: Vec<f64>,

    pub(crate) Kmat:SparseMatrix,
    global_node_to_model_node: HashMap<usize, usize>,
    pub(crate) Bmat: Vec<[[[f64;24];6];8]>,
    pub(crate) Kloc:Vec<[[f64;24];24]>,
    //derivMat:[[[f64;8];3];8],
    pub(crate) nbcs:usize,
    pub(crate) loads_bcs:Vec<f64>,
    ubcs: Vec<f64>,
    pub(crate) logsig:bool,  //print how many iterations it took to converge
    pub(crate) isactive:Vec<bool>,
    pub(crate) elemcount:Vec<usize>,
    pub(crate) Mi: Vec<f64>,
    pub(crate) z0:Vec<f64>,
    pub(crate) r: Vec<f64>,
    pub(crate) p: Vec<f64>,
    pub(crate) tmp: Vec<f64>,
    pub(crate) ndcount: Vec<usize>,
    pub(crate) dx:f64, pub(crate) dy:f64, pub(crate) dz:f64,
    pub(crate) node_to_element: Vec<Vec<usize>>,
    pub(crate) dof_to_element: Vec<Vec<usize>>,
    pub(crate) dof_to_element_dof: Vec<Vec<usize>>,
    pub(crate) strains: Vec<[f64;6]>,
    pub(crate) has_bottom: Vec<bool>,
    pub(crate) volume: Vec<f64>,
    pub(crate) density:f64,
    pub(crate) e_pl: Vec<[f64;6]>,
    pub(crate) e_pl0: Vec<[f64;6]>,
    pub(crate) e_el: Vec<[f64;6]>,
    pub(crate) stress: Vec<[f64;6]>,

}

impl StructuralModel{
    pub fn with_capacity(capacity:usize,dx:f64,dy:f64,dz:f64, density:f64)->StructuralModel{
       // let derivMat = StructuralModel::generate_derivMat(dx,dy,dz);
       // let Bmat = StructuralModel::generate_B(&derivMat);

        StructuralModel{
            nodes: Vec::with_capacity(capacity * 4),
            elements: Vec::with_capacity(capacity),
            loads: Vec::with_capacity(capacity),
            structural_bcs: Vec::with_capacity(capacity),
            bcs: Vec::with_capacity(capacity),
            u: Vec::with_capacity(capacity*4),
            uel: Vec::with_capacity(capacity*4),
            Kmat: SparseMatrix::with_capacity(capacity*4),
            global_node_to_model_node:HashMap::with_capacity(capacity*4),
            Bmat: Vec::with_capacity(capacity),
            Kloc: Vec::with_capacity(capacity),
            //derivMat:derivMat,
            nbcs:0,
            loads_bcs:Vec::with_capacity(capacity),
            ubcs: Vec::with_capacity(capacity),
            logsig: false,
            isactive: Vec::with_capacity(capacity),
            elemcount: Vec::with_capacity(capacity),
            Mi: Vec::with_capacity(capacity),
            z0: Vec::with_capacity(capacity),
            r: Vec::with_capacity(capacity),
            p: Vec::with_capacity(capacity),
            tmp: Vec::with_capacity(capacity),
            ndcount: Vec::with_capacity(capacity),
            dx, dy, dz,
            node_to_element: Vec::with_capacity(capacity*4),
            strains: Vec::with_capacity(capacity),
            dof_to_element:Vec::with_capacity(capacity*4),
            dof_to_element_dof: Vec::with_capacity(capacity*4),
            has_bottom: Vec::with_capacity(capacity),
            volume: Vec::with_capacity(capacity),
            density,
            e_el: Vec::with_capacity(capacity),
            e_pl: Vec::with_capacity(capacity),
            e_pl0: Vec::with_capacity(capacity),
            stress: Vec::with_capacity(capacity)
        }
    }
    /// generates derivMat
    fn generate_derivMat(dx:f64,dy:f64,dz:f64)->[[[f64;8];3];8]    {
        let matJac = [dx/2.0, dy/2.0, dz/2.0];
        let mut derivMat = [[[0.0_f64;8];3];8];
        for i in 0..derivMat.len() {
            for j in 0..3 {
                for k in 0..8 {
                    derivMat[i][j][k] = dN[i][j][k] / matJac[j];
                }
            }
        }
        return derivMat
    }

    /// generates strain_displacement_matrix at 2x2x2 integration points
    fn generate_B(derivMat:&[[[f64;8];3];8])->[[[f64;24];6];8]    {
        let mut Bmat = [[[0.0;24];6];8];
        for i in 0..8 {
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
        }
        return Bmat;
    }

    /// add element to the model
    pub fn add_element(&mut self, nodes: &Vec<Node>, elem: [usize;8], bcs:&[bool;24], has_bottom:bool)  {
        let mut model_node = vec![0;8];
        self.elemcount.push(self.elemcount.len());
        self.isactive.push(true);
        self.strains.push([0.0;6]);
        self.Bmat.push([[[0.0;24];6];8]);
        self.has_bottom.push(has_bottom);
        self.volume.push(self.dx*self.dy*self.dz);
        self.e_el.push([0.0_f64;6]);
        self.e_pl.push([0.0_f64;6]);
        self.e_pl0.push([0.0_f64;6]);
        self.stress.push([0.0_f64;6]);

        for i in 0..8{
            let current_node = elem[i];
            if self.global_node_to_model_node.contains_key(&current_node){
                let convert_node = self.global_node_to_model_node.get(&current_node).expect("can't get key");
                model_node[i] = *convert_node;
                //let bc_st_pt = self.node_to_element.len()*3;
                self.node_to_element[model_node[i]].push(self.elements.len());
                self.dof_to_element[model_node[i]*3].push(self.elements.len());
                self.dof_to_element[model_node[i]*3+1].push(self.elements.len());
                self.dof_to_element[model_node[i]*3+2].push(self.elements.len());
                self.dof_to_element_dof[model_node[i]*3].push(i*3);
                self.dof_to_element_dof[model_node[i]*3+1].push(i*3+1);
                self.dof_to_element_dof[model_node[i]*3+2].push(i*3+2);

            }
            else{
                self.global_node_to_model_node.insert(current_node, self.nodes.len());
                model_node[i] = self.nodes.len();

                self.node_to_element.push(Vec::with_capacity(8)); // a maximum of 8 elements can connect to a single node in a brick model
                self.node_to_element[model_node[i]].push(self.elements.len());

                //what element?
                self.dof_to_element.push(Vec::with_capacity(24));
                self.dof_to_element.push(Vec::with_capacity(24));
                self.dof_to_element.push(Vec::with_capacity(24));
                self.dof_to_element[model_node[i]*3].push(self.elements.len());
                self.dof_to_element[model_node[i]*3+1].push(self.elements.len());
                self.dof_to_element[model_node[i]*3+2].push(self.elements.len());

                //which dof in Kloc of the element?
                self.dof_to_element_dof.push(Vec::with_capacity(24));
                self.dof_to_element_dof.push(Vec::with_capacity(24));
                self.dof_to_element_dof.push(Vec::with_capacity(24));
                self.dof_to_element_dof[model_node[i]*3].push(i*3);
                self.dof_to_element_dof[model_node[i]*3+1].push(i*3+1);
                self.dof_to_element_dof[model_node[i]*3+2].push(i*3+2);




                self.nodes.push(vec![nodes[current_node].pt.x, nodes[current_node].pt.y, nodes[current_node].pt.z]);
                for j in 0..3 {
                    self.u.push(0.0);
                    self.uel.push(0.0);
                    self.bcs.push(bcs[i*3+j]);
                    self.structural_bcs.push(bcs[i*3+j]);
                    self.Mi.push(0.0_f64);
                    self.z0.push(0.0_f64);
                    self.r.push(0.0_f64);
                    self.p.push(0.0_f64);
                    self.tmp.push(0.0_f64);
                    self.loads.push(0.0);
                    //self.ndcount.push(self.ndcount.len());
                    if bcs[i*3+j] {
                        self.nbcs += 1;
                    }
                    else{
                        // self.loads_bcs.push(0.0);
                        //self.ubcs.push(0.0)
                    }
                }
            }
        }
        self.elements.push(model_node);

        self.Kloc.push([[0.0_f64;24];24]);

        //println!("Kloc len {}", self.Kloc.len());
    }



    /// generates local stiffness matrix and stores the values in the model
    pub(crate) fn generate_local_stiffness_matrix_par( Klocs: &mut [[[f64;24];24]], E:&[f64], nu:f64, elements:&[Vec<usize>], nodes:&[Vec<f64>], Bmatrix: &mut [[[[f64;24];6];8]], C2:&[[f64;6];6],
                                                       tmp1: &mut[Vec<[f64;6]>],tmp2:&mut[Vec<[f64;24]>], Ctmp:&mut[Vec<[f64;6]>], dxdydz: &[f64;3], isactive:&[bool],volume:&mut[f64],
                                                       quadraturevals:&[[f64;4];8],elemcount: &[usize], numthreads:usize, maxthreads:usize, pool: &ThreadPool )
    {
        if numthreads < maxthreads{
            let splitpos = Klocs.len()/2;
            //if Klocs.len() < 5 {
            //   println!("Klocs len {}, elements len{}", Klocs.len(), elements.len());
            //  }
            let(k1,k2) = Klocs.split_at_mut(splitpos);
            let (E1,E2) = E.split_at(splitpos);
            let (elems1, elems2) = elements.split_at(splitpos);
            let (elc1, elc2) = elemcount.split_at(splitpos);
            let (active1, active2) = isactive.split_at(splitpos);
            let (Bmat1, Bmat2) = Bmatrix.split_at_mut(splitpos);
            let (vol1, vol2) = volume.split_at_mut(splitpos);

            let tmpsplitpos = tmp1.len()/2;
            let (tmp11, tmp12) = tmp1.split_at_mut(tmpsplitpos);
            let (tmp21, tmp22)  = tmp2.split_at_mut(tmpsplitpos);
            let (Ctmp1, Ctmp2) = Ctmp.split_at_mut(tmpsplitpos);



            pool.install(||rayon::join(
                ||StructuralModel::generate_local_stiffness_matrix_par(k1,E1, nu, elems1, nodes, Bmat1,C2,
                                                                       tmp11,tmp21, Ctmp1,
                                                                       dxdydz,active1, vol1,quadraturevals,elc1, numthreads*2, maxthreads, pool),
                ||StructuralModel::generate_local_stiffness_matrix_par(k2,E2, nu, elems2, nodes, Bmat2,C2,
                                                                       tmp12, tmp22, Ctmp2,
                                                                       dxdydz,active2, vol2,quadraturevals,elc2, numthreads*2, maxthreads,pool)
            ));
        }
        else{
            let mut nodal_coords = [[0.0;3];8];

            for j in 0..elements.len() {
                if isactive[j] {
                    //SparseMatrix::create_zero_dense_matrix(8, 3);
                    for ii in 0..8 {
                        for jj in 0..3 {
                            nodal_coords[ii][jj] = nodes[elements [j] [ii] ] [jj];
                        }
                    }
                    //let vol = 0.0;




                    // let dx =  f64::max(nodal_coords[7][0],nodal_coords[6][0]) - f64::min(nodal_coords[0][0],nodal_coords[1][0]);
                    // let dy = f64::max(nodal_coords[7][1],nodal_coords[6][1]) - f64::min(nodal_coords[0][1],nodal_coords[1][1]);
                    // let dz = f64::max(nodal_coords[7][2],nodal_coords[6][2]) - f64::min(nodal_coords[0][2],nodal_coords[1][2]);
                    // let vol = f64::abs((dx*dy*dz));
                    let mut vol = dxdydz[0] * dxdydz[1] * dxdydz[2];
                    //let matJac = dxdydz;//[dx,dy,dz];
                    //let B = Bmat;
                    // create local stiffness matrix
                    for ii in 0..Ctmp[0].len() {
                        for jj in 0..6 {
                            Ctmp[0][ii][jj] = C2[ii][jj] * E[j] / (1.0 + nu) / (1.0 - 2.0 * nu);
                        }
                    }
                    for kk in 0..24 {
                        for ll in 0..24 {
                            Klocs[j][kk][ll] = 0.0_f64;
                        }
                    }

                    let mut Bmat = [[0.0;24];6];


                    for kk in 0..8 {

                        // generate jacobian matrix for the element
                        let mj = SparseMatrix::dense_mat_mult_easy2(&dN[kk],&nodal_coords);
                       // let detmatjac = mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0];
                       // vol += detmatjac;
                        // generate derivateMatrix for the element
                        let derivMat = [  [  (dN[kk][0][0]*mj[1][1]*mj[2][2] - dN[kk][0][0]*mj[1][2]*mj[2][1] - dN[kk][1][0]*mj[0][1]*mj[2][2] + dN[kk][1][0]*mj[0][2]*mj[2][1] + dN[kk][2][0]*mj[0][1]*mj[1][2] - dN[kk][2][0]*mj[0][2]*mj[1][1])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][1]*mj[1][1]*mj[2][2] - dN[kk][0][1]*mj[1][2]*mj[2][1] - dN[kk][1][1]*mj[0][1]*mj[2][2] + dN[kk][1][1]*mj[0][2]*mj[2][1] + dN[kk][2][1]*mj[0][1]*mj[1][2] - dN[kk][2][1]*mj[0][2]*mj[1][1])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][2]*mj[1][1]*mj[2][2] - dN[kk][0][2]*mj[1][2]*mj[2][1] - dN[kk][1][2]*mj[0][1]*mj[2][2] + dN[kk][1][2]*mj[0][2]*mj[2][1] + dN[kk][2][2]*mj[0][1]*mj[1][2] - dN[kk][2][2]*mj[0][2]*mj[1][1])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][3]*mj[1][1]*mj[2][2] - dN[kk][0][3]*mj[1][2]*mj[2][1] - dN[kk][1][3]*mj[0][1]*mj[2][2] + dN[kk][1][3]*mj[0][2]*mj[2][1] + dN[kk][2][3]*mj[0][1]*mj[1][2] - dN[kk][2][3]*mj[0][2]*mj[1][1])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][4]*mj[1][1]*mj[2][2] - dN[kk][0][4]*mj[1][2]*mj[2][1] - dN[kk][1][4]*mj[0][1]*mj[2][2] + dN[kk][1][4]*mj[0][2]*mj[2][1] + dN[kk][2][4]*mj[0][1]*mj[1][2] - dN[kk][2][4]*mj[0][2]*mj[1][1])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][5]*mj[1][1]*mj[2][2] - dN[kk][0][5]*mj[1][2]*mj[2][1] - dN[kk][1][5]*mj[0][1]*mj[2][2] + dN[kk][1][5]*mj[0][2]*mj[2][1] + dN[kk][2][5]*mj[0][1]*mj[1][2] - dN[kk][2][5]*mj[0][2]*mj[1][1])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][6]*mj[1][1]*mj[2][2] - dN[kk][0][6]*mj[1][2]*mj[2][1] - dN[kk][1][6]*mj[0][1]*mj[2][2] + dN[kk][1][6]*mj[0][2]*mj[2][1] + dN[kk][2][6]*mj[0][1]*mj[1][2] - dN[kk][2][6]*mj[0][2]*mj[1][1])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][7]*mj[1][1]*mj[2][2] - dN[kk][0][7]*mj[1][2]*mj[2][1] - dN[kk][1][7]*mj[0][1]*mj[2][2] + dN[kk][1][7]*mj[0][2]*mj[2][1] + dN[kk][2][7]*mj[0][1]*mj[1][2] - dN[kk][2][7]*mj[0][2]*mj[1][1])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0])],
                            [ -(dN[kk][0][0]*mj[1][0]*mj[2][2] - dN[kk][0][0]*mj[1][2]*mj[2][0] - dN[kk][1][0]*mj[0][0]*mj[2][2] + dN[kk][1][0]*mj[0][2]*mj[2][0] + dN[kk][2][0]*mj[0][0]*mj[1][2] - dN[kk][2][0]*mj[0][2]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]), -(dN[kk][0][1]*mj[1][0]*mj[2][2] - dN[kk][0][1]*mj[1][2]*mj[2][0] - dN[kk][1][1]*mj[0][0]*mj[2][2] + dN[kk][1][1]*mj[0][2]*mj[2][0] + dN[kk][2][1]*mj[0][0]*mj[1][2] - dN[kk][2][1]*mj[0][2]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]), -(dN[kk][0][2]*mj[1][0]*mj[2][2] - dN[kk][0][2]*mj[1][2]*mj[2][0] - dN[kk][1][2]*mj[0][0]*mj[2][2] + dN[kk][1][2]*mj[0][2]*mj[2][0] + dN[kk][2][2]*mj[0][0]*mj[1][2] - dN[kk][2][2]*mj[0][2]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]), -(dN[kk][0][3]*mj[1][0]*mj[2][2] - dN[kk][0][3]*mj[1][2]*mj[2][0] - dN[kk][1][3]*mj[0][0]*mj[2][2] + dN[kk][1][3]*mj[0][2]*mj[2][0] + dN[kk][2][3]*mj[0][0]*mj[1][2] - dN[kk][2][3]*mj[0][2]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]), -(dN[kk][0][4]*mj[1][0]*mj[2][2] - dN[kk][0][4]*mj[1][2]*mj[2][0] - dN[kk][1][4]*mj[0][0]*mj[2][2] + dN[kk][1][4]*mj[0][2]*mj[2][0] + dN[kk][2][4]*mj[0][0]*mj[1][2] - dN[kk][2][4]*mj[0][2]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]), -(dN[kk][0][5]*mj[1][0]*mj[2][2] - dN[kk][0][5]*mj[1][2]*mj[2][0] - dN[kk][1][5]*mj[0][0]*mj[2][2] + dN[kk][1][5]*mj[0][2]*mj[2][0] + dN[kk][2][5]*mj[0][0]*mj[1][2] - dN[kk][2][5]*mj[0][2]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]), -(dN[kk][0][6]*mj[1][0]*mj[2][2] - dN[kk][0][6]*mj[1][2]*mj[2][0] - dN[kk][1][6]*mj[0][0]*mj[2][2] + dN[kk][1][6]*mj[0][2]*mj[2][0] + dN[kk][2][6]*mj[0][0]*mj[1][2] - dN[kk][2][6]*mj[0][2]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]), -(dN[kk][0][7]*mj[1][0]*mj[2][2] - dN[kk][0][7]*mj[1][2]*mj[2][0] - dN[kk][1][7]*mj[0][0]*mj[2][2] + dN[kk][1][7]*mj[0][2]*mj[2][0] + dN[kk][2][7]*mj[0][0]*mj[1][2] - dN[kk][2][7]*mj[0][2]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0])],
                            [  (dN[kk][0][0]*mj[1][0]*mj[2][1] - dN[kk][0][0]*mj[1][1]*mj[2][0] - dN[kk][1][0]*mj[0][0]*mj[2][1] + dN[kk][1][0]*mj[0][1]*mj[2][0] + dN[kk][2][0]*mj[0][0]*mj[1][1] - dN[kk][2][0]*mj[0][1]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][1]*mj[1][0]*mj[2][1] - dN[kk][0][1]*mj[1][1]*mj[2][0] - dN[kk][1][1]*mj[0][0]*mj[2][1] + dN[kk][1][1]*mj[0][1]*mj[2][0] + dN[kk][2][1]*mj[0][0]*mj[1][1] - dN[kk][2][1]*mj[0][1]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][2]*mj[1][0]*mj[2][1] - dN[kk][0][2]*mj[1][1]*mj[2][0] - dN[kk][1][2]*mj[0][0]*mj[2][1] + dN[kk][1][2]*mj[0][1]*mj[2][0] + dN[kk][2][2]*mj[0][0]*mj[1][1] - dN[kk][2][2]*mj[0][1]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][3]*mj[1][0]*mj[2][1] - dN[kk][0][3]*mj[1][1]*mj[2][0] - dN[kk][1][3]*mj[0][0]*mj[2][1] + dN[kk][1][3]*mj[0][1]*mj[2][0] + dN[kk][2][3]*mj[0][0]*mj[1][1] - dN[kk][2][3]*mj[0][1]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][4]*mj[1][0]*mj[2][1] - dN[kk][0][4]*mj[1][1]*mj[2][0] - dN[kk][1][4]*mj[0][0]*mj[2][1] + dN[kk][1][4]*mj[0][1]*mj[2][0] + dN[kk][2][4]*mj[0][0]*mj[1][1] - dN[kk][2][4]*mj[0][1]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][5]*mj[1][0]*mj[2][1] - dN[kk][0][5]*mj[1][1]*mj[2][0] - dN[kk][1][5]*mj[0][0]*mj[2][1] + dN[kk][1][5]*mj[0][1]*mj[2][0] + dN[kk][2][5]*mj[0][0]*mj[1][1] - dN[kk][2][5]*mj[0][1]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][6]*mj[1][0]*mj[2][1] - dN[kk][0][6]*mj[1][1]*mj[2][0] - dN[kk][1][6]*mj[0][0]*mj[2][1] + dN[kk][1][6]*mj[0][1]*mj[2][0] + dN[kk][2][6]*mj[0][0]*mj[1][1] - dN[kk][2][6]*mj[0][1]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0]),  (dN[kk][0][7]*mj[1][0]*mj[2][1] - dN[kk][0][7]*mj[1][1]*mj[2][0] - dN[kk][1][7]*mj[0][0]*mj[2][1] + dN[kk][1][7]*mj[0][1]*mj[2][0] + dN[kk][2][7]*mj[0][0]*mj[1][1] - dN[kk][2][7]*mj[0][1]*mj[1][0])/(mj[0][0]*mj[1][1]*mj[2][2] - mj[0][0]*mj[1][2]*mj[2][1] - mj[0][1]*mj[1][0]*mj[2][2] + mj[0][1]*mj[1][2]*mj[2][0] + mj[0][2]*mj[1][0]*mj[2][1] - mj[0][2]*mj[1][1]*mj[2][0])]];
                        // generate displacement to strain matrix for element

                        for jj in 0..8 {
                            Bmat[0] [jj * 3] = derivMat[0][jj];
                            Bmat[1] [jj * 3 + 1] = derivMat[1][jj];
                            Bmat[2] [jj * 3 + 2] = derivMat[2][jj];

                            Bmat[3] [jj * 3] = derivMat[1][jj];
                            Bmat[3] [jj * 3 + 1] = derivMat[0] [jj];

                            Bmat[4] [jj * 3 + 2] = derivMat[1][jj];
                            Bmat[4] [jj * 3 + 1] = derivMat[2][jj];

                            Bmat[5] [jj * 3] = derivMat[2][jj];
                            Bmat[5] [jj * 3 + 2] = derivMat[0][jj];

                            Bmatrix[j][kk] [0] [jj * 3] = derivMat[0][jj];
                            Bmatrix[j][kk] [1] [jj * 3 + 1] = derivMat[1][jj];
                            Bmatrix[j][kk] [2] [jj * 3 + 2] = derivMat[2][jj];

                            Bmatrix[j][kk] [3] [jj * 3] = derivMat[1][jj];
                            Bmatrix[j][kk] [3] [jj * 3 + 1] = derivMat[0] [jj];

                            Bmatrix[j][kk] [4] [jj * 3 + 2] = derivMat[1][jj];
                            Bmatrix[j][kk] [4] [jj * 3 + 1] = derivMat[2][jj];

                            Bmatrix[j][kk] [5] [jj * 3] = derivMat[2][jj];
                            Bmatrix[j][kk] [5] [jj * 3 + 2] = derivMat[0][jj];
                        }


                        let mut t1 = SparseMatrix::dense_mat_mult_easy_transpose(&Bmat, &Ctmp[0]);
                        let mut t2 = SparseMatrix::dense_mat_mult_easy(&t1, &Bmat);
                        SparseMatrix::multiply_dense_by_constant(&mut t2, quadraturevals[kk][3] * vol/8.0);
                        //SparseMatrix::KlocplusBtCBdetmatjac(&B[kk],&C,vol/8.0,&mut tmp1, &mut tmp2);
                        SparseMatrix::add_dense_matrix(&mut Klocs[j], &mut t2);
                    }
                    volume[j] = vol;
                }
            }

        }
    }



    pub fn reset(&mut self, numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        self.loads_bcs.clear();
        self.ubcs.clear();

        //self.Kloc.clear();
        SparseMatrix::assign_val_par(self.loads.borrow_mut(), 0.0_f64, 1, maxthreads,pool);
        //for i in 0..self.loads.len(){
        //   self.loads[i] = 0.0_f64;
        // }
    }


    pub fn calc_strain_for_all_elements(&mut self){
        for element_no in 0..self.elements.len() {
            let u = [
                self.u[self.elements[element_no][0] * 3], self.u[self.elements[element_no][0] * 3 + 1], self.u[self.elements[element_no][0] * 3 + 2],
                self.u[self.elements[element_no][1] * 3], self.u[self.elements[element_no][1] * 3 + 1], self.u[self.elements[element_no][1] * 3 + 2],
                self.u[self.elements[element_no][2] * 3], self.u[self.elements[element_no][2] * 3 + 1], self.u[self.elements[element_no][2] * 3 + 2],
                self.u[self.elements[element_no][3] * 3], self.u[self.elements[element_no][3] * 3 + 1], self.u[self.elements[element_no][3] * 3 + 2],
                self.u[self.elements[element_no][4] * 3], self.u[self.elements[element_no][4] * 3 + 1], self.u[self.elements[element_no][4] * 3 + 2],
                self.u[self.elements[element_no][5] * 3], self.u[self.elements[element_no][5] * 3 + 1], self.u[self.elements[element_no][5] * 3 + 2],
                self.u[self.elements[element_no][6] * 3], self.u[self.elements[element_no][6] * 3 + 1], self.u[self.elements[element_no][6] * 3 + 2],
                self.u[self.elements[element_no][7] * 3], self.u[self.elements[element_no][7] * 3 + 1], self.u[self.elements[element_no][7] * 3 + 2]
            ];

            for i in 0..6{
                self.strains[element_no][i] = 0.0;
            }

            for i in 0..8 {
                self.strains[element_no][0] += self.Bmat[element_no][i][0][0] * u[0] + self.Bmat[element_no][i][0][3] * u[3] + self.Bmat[element_no][i][0][6] * u[6] + self.Bmat[element_no][i][0][9] * u[9] + self.Bmat[element_no][i][0][12] * u[12] + self.Bmat[element_no][i][0][15] * u[15] + self.Bmat[element_no][i][0][18] * u[18] + self.Bmat[element_no][i][0][20] * u[20];
                self.strains[element_no][1] += self.Bmat[element_no][i][1][1] * u[1] + self.Bmat[element_no][i][1][3 + 1] * u[3 + 1] + self.Bmat[element_no][i][1][6 + 1] * u[6 + 1] + self.Bmat[element_no][i][1][9 + 1] * u[9 + 1] + self.Bmat[element_no][i][1][12 + 1] * u[12 + 1] + self.Bmat[element_no][i][1][15 + 1] * u[15 + 1] + self.Bmat[element_no][i][1][18 + 1] * u[18 + 1] + self.Bmat[element_no][i][1][20 + 1] * u[20 + 1];
                self.strains[element_no][2] += self.Bmat[element_no][i][2][2] * u[2] + self.Bmat[element_no][i][2][3 + 2] * u[3 + 2] + self.Bmat[element_no][i][2][6 + 2] * u[6 + 2] + self.Bmat[element_no][i][2][9 + 2] * u[9 + 2] + self.Bmat[element_no][i][2][12 + 2] * u[12 + 2] + self.Bmat[element_no][i][2][15 + 2] * u[15 + 2] + self.Bmat[element_no][i][2][18 + 1] * u[18 + 2] + self.Bmat[element_no][i][2][20 + 1] * u[20 + 2];
                for j in 0..24 {
                    self.strains[element_no][3] += self.Bmat[element_no][i][3][j] * u[j];
                    self.strains[element_no][4] += self.Bmat[element_no][i][4][j] * u[j];
                    self.strains[element_no][5] += self.Bmat[element_no][i][5][j] * u[j];
                }
            }
        }
    }
    pub(crate) fn calc_strain_for_all_elements_par(uel:&[f64], elements:&[Vec<usize>], Bmat:&[[[[f64;24];6];8]], strains: &mut[[f64;6]], isactive: &[bool],
                                                   numthreads: usize, maxthreads: usize, pool: &ThreadPool ){
        if numthreads < maxthreads {
            let splitpos = elements.len() / 2;
            let (el1, el2) = elements.split_at(splitpos);
            let (Bmat1, Bmat2) = Bmat.split_at(splitpos);
            let (strains1, strains2) = strains.split_at_mut(splitpos);
            let (isactive1, isactive2) = isactive.split_at(splitpos);

            pool.install(|| rayon::join(
                || StructuralModel::calc_strain_for_all_elements_par(uel, el1, Bmat1, strains1, isactive1, numthreads * 2, maxthreads, pool),
                || StructuralModel::calc_strain_for_all_elements_par(uel, el2, Bmat2, strains2, isactive2, numthreads * 2, maxthreads, pool),
            ));
        }
        else{
            for element_no in 0..elements.len() {
                if isactive[element_no] {
                    let uloc = [
                        uel[elements[element_no][0] * 3], uel[elements[element_no][0] * 3 + 1], uel[elements[element_no][0] * 3 + 2],
                        uel[elements[element_no][1] * 3], uel[elements[element_no][1] * 3 + 1], uel[elements[element_no][1] * 3 + 2],
                        uel[elements[element_no][2] * 3], uel[elements[element_no][2] * 3 + 1], uel[elements[element_no][2] * 3 + 2],
                        uel[elements[element_no][3] * 3], uel[elements[element_no][3] * 3 + 1], uel[elements[element_no][3] * 3 + 2],
                        uel[elements[element_no][4] * 3], uel[elements[element_no][4] * 3 + 1], uel[elements[element_no][4] * 3 + 2],
                        uel[elements[element_no][5] * 3], uel[elements[element_no][5] * 3 + 1], uel[elements[element_no][5] * 3 + 2],
                        uel[elements[element_no][6] * 3], uel[elements[element_no][6] * 3 + 1], uel[elements[element_no][6] * 3 + 2],
                        uel[elements[element_no][7] * 3], uel[elements[element_no][7] * 3 + 1], uel[elements[element_no][7] * 3 + 2]
                    ];
                    for i in 0..6{
                        strains[element_no][i] = 0.0;
                    }
                    for i in 0..8 {
                        for j in 0..24{
                            strains[element_no][0] += Bmat[element_no][i][0][j] * uloc[j];// + Bmat[element_no][i][0][3] * uel[3] + Bmat[element_no][i][0][6] * uel[6] + Bmat[element_no][i][0][9] * uel[9] + Bmat[element_no][i][0][12] * uel[12] + Bmat[element_no][i][0][15] * uel[15] + Bmat[element_no][i][0][18] * uel[18] + Bmat[element_no][i][0][20] * uel[20];
                        strains[element_no][1] += Bmat[element_no][i][1][j] * uloc[j];// + Bmat[element_no][i][1][3 + 1] * uel[3 + 1] + Bmat[element_no][i][1][6 + 1] * uel[6 + 1] + Bmat[element_no][i][1][9 + 1] * uel[9 + 1] + Bmat[element_no][i][1][12 + 1] * uel[12 + 1] + Bmat[element_no][i][1][15 + 1] * uel[15 + 1] + Bmat[element_no][i][1][18 + 1] * uel[18 + 1] + Bmat[element_no][i][1][20 + 1] * uel[20 + 1];
                        strains[element_no][2] += Bmat[element_no][i][2][j] * uloc[j];// + Bmat[element_no][i][2][3 + 2] * uel[3 + 2] + Bmat[element_no][i][2][6 + 2] * uel[6 + 2] + Bmat[element_no][i][2][9 + 2] * uel[9 + 2] + Bmat[element_no][i][2][12 + 2] * uel[12 + 2] + Bmat[element_no][i][2][15 + 2] * uel[15 + 2] + Bmat[element_no][i][2][18 + 1] * uel[18 + 2] + Bmat[element_no][i][2][20 + 1] * uel[20 + 2];

                            strains[element_no][3] += Bmat[element_no][i][3][j] * uloc[j];
                            strains[element_no][4] += Bmat[element_no][i][4][j] * uloc[j];
                            strains[element_no][5] += Bmat[element_no][i][5][j] * uloc[j];
                        }
                    }
                }
            }

        }

    }

    pub(crate) fn calc_stress_for_all_elements_par(uel:&[f64], elements:&[Vec<usize>], Bmat:&[[[[f64;24];6];8]], stress: &mut[[f64;6]], strains: &mut[[f64;6]], epl: &mut[[f64;6]], eel: &mut[[f64;6]], epl_0: &mut [[f64; 6]],
                                                   dt:f64,
                                                   isactive: &[bool],
                                                   G:&[f64], nu:&f64, vis:&[f64], numthreads: usize, maxthreads: usize, pool: &ThreadPool ){
        if numthreads < maxthreads {
            let splitpos = stress.len() / 2;
            let (el1, el2) = elements.split_at(splitpos);
            let (Bmat1, Bmat2) = Bmat.split_at(splitpos);
            let (stress1, stress2) = stress.split_at_mut(splitpos);
            let (strains1, strains2) = strains.split_at_mut(splitpos);
            let (eel1, eel2) = eel.split_at_mut(splitpos);
            let (epl1, epl2) = epl.split_at_mut(splitpos);
            let (epl01, epl02) = epl_0.split_at_mut(splitpos);
            let (isactive1, isactive2) = isactive.split_at(splitpos);
            let (G1, G2) = G.split_at(splitpos);
            let (vis1, vis2) = vis.split_at(splitpos);

            pool.install(|| rayon::join(
                || StructuralModel::calc_stress_for_all_elements_par(uel, el1, Bmat1, stress1,strains1, epl1,eel1, epl01,dt, isactive1, G1, nu, vis1,numthreads * 2, maxthreads, pool),
                || StructuralModel::calc_stress_for_all_elements_par(uel, el2, Bmat2, stress2, strains2, epl2,eel2, epl02,dt, isactive2, G2, nu, vis2,numthreads * 2, maxthreads, pool),
            ));
        }
        else{

            for element_no in 0..elements.len() {
                if isactive[element_no] {
                    let uloc = [
                        uel[elements[element_no][0] * 3], uel[elements[element_no][0] * 3 + 1], uel[elements[element_no][0] * 3 + 2],
                        uel[elements[element_no][1] * 3], uel[elements[element_no][1] * 3 + 1], uel[elements[element_no][1] * 3 + 2],
                        uel[elements[element_no][2] * 3], uel[elements[element_no][2] * 3 + 1], uel[elements[element_no][2] * 3 + 2],
                        uel[elements[element_no][3] * 3], uel[elements[element_no][3] * 3 + 1], uel[elements[element_no][3] * 3 + 2],
                        uel[elements[element_no][4] * 3], uel[elements[element_no][4] * 3 + 1], uel[elements[element_no][4] * 3 + 2],
                        uel[elements[element_no][5] * 3], uel[elements[element_no][5] * 3 + 1], uel[elements[element_no][5] * 3 + 2],
                        uel[elements[element_no][6] * 3], uel[elements[element_no][6] * 3 + 1], uel[elements[element_no][6] * 3 + 2],
                        uel[elements[element_no][7] * 3], uel[elements[element_no][7] * 3 + 1], uel[elements[element_no][7] * 3 + 2]
                    ];
                    for i in 0..6
                    {
                        eel[element_no][i] = 0.0;
                        strains[element_no] [i] = 0.0;
                    }

                    for i in 0..8 {
                        for j in 0..24{
                            eel[element_no][0] += Bmat[element_no][i][0][j] * uloc[j];// + Bmat[element_no][i][0][3] * uel[3] + Bmat[element_no][i][0][6] * uel[6] + Bmat[element_no][i][0][9] * uel[9] + Bmat[element_no][i][0][12] * uel[12] + Bmat[element_no][i][0][15] * uel[15] + Bmat[element_no][i][0][18] * uel[18] + Bmat[element_no][i][0][20] * uel[20];
                            eel[element_no][1] += Bmat[element_no][i][1][j] * uloc[j];// + Bmat[element_no][i][1][3 + 1] * uel[3 + 1] + Bmat[element_no][i][1][6 + 1] * uel[6 + 1] + Bmat[element_no][i][1][9 + 1] * uel[9 + 1] + Bmat[element_no][i][1][12 + 1] * uel[12 + 1] + Bmat[element_no][i][1][15 + 1] * uel[15 + 1] + Bmat[element_no][i][1][18 + 1] * uel[18 + 1] + Bmat[element_no][i][1][20 + 1] * uel[20 + 1];
                            eel[element_no][2] += Bmat[element_no][i][2][j] * uloc[j];// + Bmat[element_no][i][2][3 + 2] * uel[3 + 2] + Bmat[element_no][i][2][6 + 2] * uel[6 + 2] + Bmat[element_no][i][2][9 + 2] * uel[9 + 2] + Bmat[element_no][i][2][12 + 2] * uel[12 + 2] + Bmat[element_no][i][2][15 + 2] * uel[15 + 2] + Bmat[element_no][i][2][18 + 1] * uel[18 + 2] + Bmat[element_no][i][2][20 + 1] * uel[20 + 2];

                            eel[element_no][3] += Bmat[element_no][i][3][j] * uloc[j];
                            eel[element_no][4] += Bmat[element_no][i][4][j] * uloc[j];
                            eel[element_no][5] += Bmat[element_no][i][5][j] * uloc[j];
                        }
                    }

                    let E = G[element_no]*2.0*(1.0+nu);
                    let mult = E/(1.0-2.0*nu)/(1.0+nu);

                    let s0 = mult * (eel[element_no][0] *(1.0-nu)+eel[element_no][1]*nu+eel[element_no][2]*nu);
                    let s1 = mult * (eel[element_no][1] *(nu)+eel[element_no][1]*(1.0-nu)+eel[element_no][2]*nu);
                    let s2 = mult * (eel[element_no][2] *(nu)+eel[element_no][1]*nu+eel[element_no][2]*(1.0-nu));
                    let s3 = mult * (eel[element_no][3]) *(1.0-2.0*nu)/2.0;
                    let s4 = mult * (eel[element_no][4]) *(1.0-2.0*nu)/2.0;
                    let s5 = mult * (eel[element_no][5]) *(1.0-2.0*nu)/2.0;

                    stress[element_no][0] = s0;//mult * (strains[element_no][0] *(1.0-nu)+strains[element_no][1]*nu+strains[element_no][2]*nu);
                    stress[element_no][1] = s1;//mult * (strains[element_no][1] *(nu)+strains[element_no][1]*(1.0-nu)+strains[element_no][2]*nu);
                    stress[element_no][2] = s2;//mult * (strains[element_no][2] *(nu)+strains[element_no][1]*nu+strains[element_no][2]*(1.0-nu));
                    stress[element_no][3] = s3;//mult * (strains[element_no][3]) *(1.0-2.0*nu)/2.0;
                    stress[element_no][4] = s4;//mult * (strains[element_no][3]) *(1.0-2.0*nu)/2.0;
                    stress[element_no][5] = s5;//mult * (strains[element_no][3]) *(1.0-2.0*nu)/2.0;

                    let trace_dij:f64 = (stress[element_no][0] + stress[element_no][1] + stress[element_no][2]) / 3.0;

                    let mut dev_stress = [0.0;6];
                    dev_stress[0] = stress[element_no][0] - trace_dij;
                    dev_stress[1] = stress[element_no][1] - trace_dij;
                    dev_stress[2] = stress[element_no][2] - trace_dij;
                    dev_stress[3] = stress[element_no][3];
                    dev_stress[4] = stress[element_no][4];
                    dev_stress[5] = stress[element_no][5];

                    let mut e_dot_pl = [0.0;6];
                    for i in 0..6{
                        e_dot_pl[i] = stress[element_no][i]/2.0/vis[element_no];
                        epl[element_no][i] += e_dot_pl[i] * dt;
                        epl_0[element_no][i] = e_dot_pl[i] * dt;
                        strains[element_no][i] = eel[element_no][i]+epl[element_no][i];
                    }
                }
            }

        }

    }

    pub fn calc_strain_for_element(&self, element_no: usize)->Vec<f64>{
        let mut strainout = vec![0.0;6];
        let u = [
            self.u[self.elements[element_no][0]*3],self.u[self.elements[element_no][0]*3+1],self.u[self.elements[element_no][0]*3+2],
            self.u[self.elements[element_no][1]*3],self.u[self.elements[element_no][1]*3+1],self.u[self.elements[element_no][1]*3+2],
            self.u[self.elements[element_no][2]*3],self.u[self.elements[element_no][2]*3+1],self.u[self.elements[element_no][2]*3+2],
            self.u[self.elements[element_no][3]*3],self.u[self.elements[element_no][3]*3+1],self.u[self.elements[element_no][3]*3+2],
            self.u[self.elements[element_no][4]*3],self.u[self.elements[element_no][4]*3+1],self.u[self.elements[element_no][4]*3+2],
            self.u[self.elements[element_no][5]*3],self.u[self.elements[element_no][5]*3+1],self.u[self.elements[element_no][5]*3+2],
            self.u[self.elements[element_no][6]*3],self.u[self.elements[element_no][6]*3+1],self.u[self.elements[element_no][6]*3+2],
            self.u[self.elements[element_no][7]*3],self.u[self.elements[element_no][7]*3+1],self.u[self.elements[element_no][7]*3+2]
        ];

        for i in 0..8{
            strainout[0] += self.Bmat[element_no][i][0][0]*u[0] + self.Bmat[element_no][i][0][3]*u[3] + self.Bmat[element_no][i][0][6]*u[6] + self.Bmat[element_no][i][0][9]*u[9] + self.Bmat[element_no][i][0][12]*u[12]+ self.Bmat[element_no][i][0][15]*u[15] + self.Bmat[element_no][i][0][18]*u[18]+ self.Bmat[element_no][i][0][20]*u[20]/8.0;
            strainout[1] += self.Bmat[element_no][i][1][1]*u[1] + self.Bmat[element_no][i][1][3+1]*u[3+1] + self.Bmat[element_no][i][1][6+1]*u[6+1] + self.Bmat[element_no][i][1][9+1]*u[9+1] + self.Bmat[element_no][i][1][12+1]*u[12+1]+ self.Bmat[element_no][i][1][15+1]*u[15+1] + self.Bmat[element_no][i][1][18+1]*u[18+1]+ self.Bmat[element_no][i][1][20+1]*u[20+1]/8.0;
            strainout[2] += self.Bmat[element_no][i][2][2]*u[2] + self.Bmat[element_no][i][2][3+2]*u[3+2] + self.Bmat[element_no][i][2][6+2]*u[6+2] + self.Bmat[element_no][i][2][9+2]*u[9+2] + self.Bmat[element_no][i][2][12+2]*u[12+2]+ self.Bmat[element_no][i][2][15+2]*u[15+2] + self.Bmat[element_no][i][2][18+1]*u[18+2]+ self.Bmat[element_no][i][2][20+1]*u[20+2]/8.0;
            for j in 0..24{
                strainout[3] += self.Bmat[element_no][i][3][j]*u[j]/8.0;
                strainout[4] += self.Bmat[element_no][i][4][j]*u[j]/8.0;
                strainout[5] += self.Bmat[element_no][i][5][j]*u[j]/8.0;
            }
        }
        strainout
    }


    /// generates global stiffness matrix based on dof_to_element_dof data structure
    /// goes each row of stiffness matrix at a time
    /// looks like a good setup for lock free parallelization
    pub(crate) fn generate_global_stiffness_matrix_with_bc_rowwise(&mut self, num_bcs:usize){
        self.Kmat.reset();
        let mut dof_to_bc_dof = vec![0; self.nodes.len()*3];
        let mut dof_from_bc_dof = vec![0;self.nodes.len()*3-self.nbcs];
        let mut dofcount = 0;

        for i in 0..dof_to_bc_dof.len(){
            if !self.bcs[i]{
                dof_to_bc_dof[i] = dofcount;
                dof_from_bc_dof[dofcount] = i;
                dofcount += 1;
            }
        }

        //go through all dofs
        for i in 0..dof_from_bc_dof.len(){
            // only dofs that are not bcs
            let dofcounter = dof_from_bc_dof[i];
            //if !self.bcs[i]{
            // go through all elements that affect this dof
            for j in 0..self.dof_to_element_dof[dofcounter].len(){
                // get the row of the the element stiffness matrix that affects this dof row
                let row_Kloc = self.dof_to_element_dof[dofcounter][j];
                // get what element affects this dof row
                let elemno = self.dof_to_element[dofcounter][j];
                // get all the column dof in this row of stiffness matrix
                let cols = &self.elements[elemno];
                let Klocrow = &self.Kloc[elemno][row_Kloc];
                for k in 0..cols.len(){
                    for l in 0..3{
                        let cdof = cols[k]*3+l;
                        if !self.bcs[cdof]{
                            self.Kmat.insert(i, dof_to_bc_dof[cdof], Klocrow[k*3+l]);
                        }
                    }

                }
                // }
            }
        }

    }

    /// parallel create global stiffness matrix with bcs
    pub(crate) fn generate_global_stiffness_matrix_with_bc_rowwise_par(&mut self, num_bcs:usize, maxthreads:usize, pool:&ThreadPool){
        self.Kmat.reset();
        let mut dof_to_bc_dof = vec![0; self.nodes.len()*3];
        let mut dof_from_bc_dof = vec![0;self.nodes.len()*3-self.nbcs];
        let mut dofcount = 0;

        for i in 0..dof_to_bc_dof.len(){
            if !self.bcs[i]{
                dof_to_bc_dof[i] = dofcount;
                dof_from_bc_dof[dofcount] = i;
                dofcount += 1;
            }
        }

        SparseMatrix::insert_rows_par(self.Kmat.borrow_mut(), &dof_to_bc_dof, &dof_from_bc_dof, &self.bcs,
                                      &self.dof_to_element_dof, &self.dof_to_element, &self.elements, &self.Kloc, maxthreads, pool);

        //SparseMatrix::insert_rows_par(self.Kmat.borrow_mut(),&dof_to_bc_dof, &self.bcs,
        //&self.dof_to_element_dof, &self.dof_to_element, &self.elements,
        //&self.Kloc, maxthreads, pool);

        //go through all dofs
        /*for i in 0..self.dof_to_element_dof.len(){
            // only dofs that are not bcs
            if !self.bcs[i]{
                // go through all elements that affect this dof
                for j in 0..self.dof_to_element_dof[i].len(){
                    // get the row of the the element stiffness matrix that affects this dof row
                    let row_Kloc = self.dof_to_element_dof[i][j];
                    // get what element affects this dof row
                    let elemno = self.dof_to_element[i][j];
                    // get all the column dof in this row of stiffness matrix
                    let cols = &self.elements[elemno];
                    let Klocrow = &self.Kloc[elemno][row_Kloc];
                    for k in 0..cols.len(){
                        let cdof0 = cols[k]*3;
                        let cdof1 = cols[k]*3+1;
                        let cdof2 = cols[k]*3+2;

                        if !self.bcs[cdof0]{
                            self.Kmat.insert(dof_to_bc_dof[i], dof_to_bc_dof[cdof0], Klocrow[k*3]);
                        }
                        if !self.bcs[cdof1]{
                            self.Kmat.insert(dof_to_bc_dof[i], dof_to_bc_dof[cdof1], Klocrow[k*3+1]);
                        }
                        if !self.bcs[cdof2]{
                            self.Kmat.insert(dof_to_bc_dof[i], dof_to_bc_dof[cdof2], Klocrow[k*3+2]);
                        }

                    }
                }
            }
        }*/

    }

    /// generates global stiffness matrix with boundary conditions applied
    pub(crate) fn generate_global_stiffness_matrix_with_bc(&mut self,num_bcs:usize){//-> SparseMatrix {
        //let mut Kg_bc = SparseMatrix::with_capacity( self.nodes.len()*3-num_bcs);
        self.Kmat.reset();
        let mut dof_to_bc_dof = vec![0; self.nodes.len()*3];
        let mut dofcount = 0;

        for i in 0..dof_to_bc_dof.len(){
            if !self.bcs[i]{
                dof_to_bc_dof[i] = dofcount;
                dofcount += 1;
            }
        }
        //println!("{:?}", dof_to_bc_dof);
        // println!("{:?}", dof_to_bc_dof);
        for i in 0..self.elements.len() {
            if self.isactive[i] {
                let elem = &self.elements[i];
                let elem_dof = [
                    elem[0] * 3, elem[0] * 3 + 1, elem[0] * 3 + 2,
                    elem[1] * 3, elem[1] * 3 + 1, elem[1] * 3 + 2,
                    elem[2] * 3, elem[2] * 3 + 1, elem[2] * 3 + 2,
                    elem[3] * 3, elem[3] * 3 + 1, elem[3] * 3 + 2,
                    elem[4] * 3, elem[4] * 3 + 1, elem[4] * 3 + 2,
                    elem[5] * 3, elem[5] * 3 + 1, elem[5] * 3 + 2,
                    elem[6] * 3, elem[6] * 3 + 1, elem[6] * 3 + 2,
                    elem[7] * 3, elem[7] * 3 + 1, elem[7] * 3 + 2];

                //println!("elem dof {:?}",elem_dof);

                //map element

                for j in 0..24 {
                    if !self.bcs[elem_dof[j]] {
                        for k in 0..24 {
                            if !self.bcs[elem_dof[k]] {
                                //println!("elem_dof[j] {}", elem_dof[j]);
                                //println!("elem_dof[k] {}", elem_dof[k]);
                                // println!("dof to bc dof elem_dof[k] {}", dof_to_bc_dof[elem_dof[k]]);
                                // println!("dof to bc dof elem_dof[k] {}", dof_to_bc_dof[elem_dof[k]]);
                                self.Kmat.insert(dof_to_bc_dof[elem_dof[j]], dof_to_bc_dof[elem_dof[k]], self.Kloc[i][j][k]);
                            }
                        }
                    }
                }
            }
        }

        //return Kg_bc;
    }

    /// calculate loads from element self weight
    pub(crate) fn calc_self_wt_loads(&mut self, loads: &Vec<f64>){
        for i in 0..self.loads.len(){
            self.loads[i] = 0.0;
            self.u[i] = 0.0;
            self.uel[i] = 0.0;
        }
        for i in 0..self.elements.len() {
            if self.isactive[i]  {
                let elem = self.elements[i].clone();
                for j in 0..4 {
                    // add full loads to bottom elements? stress is a triangular distribution though... zero at top max at bottom
                    self.loads[elem[j] * 3 + 2] += self.volume[i] * -9.81 * self.density/8.0*2.0/3.0;
                }
                for j in 4..8{
                    self.loads[elem[j] * 3 + 2] += self.volume[i] * -9.81 * self.density/8.0*1.0/3.0;
                }
            }
        }
    }




    /// solves Ku=F using conjugate gradient method and stores result in the model
    pub(crate) fn structural_analysis_solve(&mut self, numthreads: usize, maxthreads: usize, pool:&ThreadPool){
        // self.generate_global_stiffness_matrix_with_bc(self.nbcs, self.bcs.borrow());
        //calculate loads_bcs
        // let mut count = 0;
        for i in 0..self.loads.len(){
            if !self.bcs[i]{
                self.loads_bcs.push(self.loads[i]);
                self.ubcs.push(0.0);
                //count += 1;
            }
        }
        //println!("self loads len {} self loadsbcs len {}", self.loads.len(), self.loads_bcs.len());

        // println!("loads len {}", self.loads.len());
        //println!("loads bcs len {}", self.loads_bcs.len());
        //println!("loads bcs");
        // for i in 0..self.loads_bcs.len()/3{
        //     println!("[{},{},{}]",self.loads_bcs[i*3],self.loads_bcs[i*3+1],self.loads_bcs[i*3+2]);
        // }


        let mut soln = self.ubcs.borrow_mut();//= vec![0.0;self.loads_bcs.len()];
        //let mut r = vec![0.0;self.loads_bcs.len()];
        // let mut p = vec![0.0;self.loads_bcs.len()];
        // let mut tmp = vec![0.0;self.loads_bcs.len()];
        self.r.resize(self.nodes.len()*3-self.nbcs,0.0_f64);
        self.p.resize(self.nodes.len()*3-self.nbcs,0.0_f64);
        self.tmp.resize(self.nodes.len()*3-self.nbcs,0.0_f64);
        self.Mi.resize(self.nodes.len()*3-self.nbcs, 0.0_f64);
        self.z0.resize(self.nodes.len()*3-self.nbcs, 0.0_f64);
        //let (act, inact) = self.elemcount.split_at(self.r.len());
        //self.Kmat.cgsolve(&self.loads_bcs, soln,self.r.borrow_mut(), self.p.borrow_mut(), self.tmp.borrow_mut(),numthreads, maxthreads, pool);
        // Kg_bc.printmatrix();
        let numsolve = self.Kmat.cgsolve_diag_scl_precon(self.loads_bcs.borrow(), soln,self.r.borrow_mut(), self.p.borrow_mut(), self.tmp.borrow_mut(),self.Mi.borrow_mut(), self.z0.borrow_mut(), numthreads, maxthreads, pool);
        //let numsolve = 400;
        if self.logsig {println!("converged in {} iterations, nbcs {}, dofs {}", numsolve, self.nbcs, self.bcs.len()-self.nbcs); self.logsig = false}
        // println!("{:?}", soln);
        let mut count = 0;
        //println!("self bcs ")
        for i in 0..self.bcs.len(){
            if self.bcs[i] {
                count += 1;
            }
            else{
                self.uel[i] = soln[i-count];                // println!("copied {} {}", i, i-count);
            }
        }
        /*for i in 0..self.nodes.len(){
            self.nodes[i][0] += self.u[i*3];
            self.nodes[i][1] += self.u[i*3+1];
            self.nodes[i][2] += self.u[i*3+2];
        }*/


    }

    /// deactivate nodes if element is inactive, compare with structural bcs, if structural bcs are true, the bc is true
    pub(crate) fn deactivate_elements(&mut self){
        //self.nbcs = 0;
        let mut activation_temp = vec![false;self.bcs.len()];
        let mut deactivation_temp = vec![false;self.bcs.len()];

        for i in 0..self.dof_to_element.len(){
            // code for deactivation only
            for j in 0..self.dof_to_element[i].len(){
                let elem = self.dof_to_element[i][j];
                if self.isactive[elem]{
                    activation_temp[i] = true;
                }
                else{
                    deactivation_temp[i] = true;
                }
            }
        }

        //go through all elements
        let mut hasbottom = vec![false; self.bcs.len() ];
        let mut hasnobottom = vec![false; self.bcs.len()];
        for i in 0..self.dof_to_element.len(){
            for j in 0..self.dof_to_element[i].len(){
            let elem = self.dof_to_element[i][j];
                if self.has_bottom[elem]{
                hasbottom[i] = true;
                }
                else{
                    hasnobottom[i] = true;
                }
            }
        }


        let mut count = 0;
        for i in 0..self.dof_to_element.len(){
            //for j in 0..3 {
                if self.structural_bcs[i] || deactivation_temp[i] || hasnobottom[i] {
                    self.bcs[i] = true;
                    self.uel[i]=0.0_f64;
                    //self.bcs[i * 3 + 1] = true;
                   // self.bcs[i * 3 + 2] = true;
                    count += 1;
                } else {
                    if activation_temp[i] {
                        self.bcs[i] = false;
                       // self.bcs[i * 3 + 1] = false;
                       // self.bcs[i * 3 + 2] = false;
                    }
              //  }
            }
        }



        self.nbcs = count;

        /*println!("bcs");
        for i in 0..self.bcs.len()/3{
            println!("{}->{},{},{}", i,  self.bcs[i*3], self.bcs[i*3+1], self.bcs[i*3+2]);
        }
        println!("structural bcs");
        for i in 0..self.bcs.len()/3{
            println!("{}->{},{},{}", i,  self.structural_bcs[i*3], self.structural_bcs[i*3+1], self.structural_bcs[i*3+2]);
        }



        //self.Mi.resize(self.nodes.len()*3-self.nbcs);
        println!("self.nbcs {}", self.nbcs);*/
    }

    pub(crate) fn write(&self, nodefile: &str, elemfile: &str, ufile: &str, loadfile: &str, bcfile: &str, strainfile:&str){
        let mut newfile1 = File::create(nodefile).expect("cant create the file");
        let mut filebuf1 = BufWriter::with_capacity(10000, newfile1);
        let mut newfile2 = File::create(elemfile).expect("cant create the file");
        let mut filebuf2 = BufWriter::with_capacity(10000, newfile2);
        let mut newfile3 = File::create(ufile).expect("cant create the file");
        let mut filebuf3 = BufWriter::with_capacity(10000, newfile3);
        let mut newfile4 = File::create(loadfile).expect("cant create the file");
        let mut filebuf4 = BufWriter::with_capacity(10000, newfile4);
        let mut newfile5 = File::create(bcfile).expect("cant create the file");
        let mut filebuf5 = BufWriter::with_capacity(10000, newfile5);
        let mut newfile6 = File::create(strainfile).expect("cant create the file");
        let mut filebuf6 = BufWriter::with_capacity(10000, newfile6);
        // writeln!(filebuf1, "nodes");
        for i in &self.nodes{
            writeln!(filebuf1, "{},{},{}", i[0],i[1],i[2]);
        }
        // writeln!(filebuf2, "elements");
        for i in &self.elements{
            writeln!(filebuf2,"{},{},{},{},{},{},{},{}",i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]);
        }
        for i in 0..self.u.len()/3{
            writeln!(filebuf3,"{},{},{}",self.u[i*3],self.u[i*3+1],self.u[i*3+2],);
        }
        //writeln!(filebuf4, "loads");
        //println!("self.loads.len ,{:?}", self.loads);
        //println!("loadfile,{:?}", loadfile);
        for i in 0..self.loads.len()/3{
            writeln!(filebuf4,"{},{},{}",self.loads[i*3],self.loads[i*3+1],self.loads[i*3+2],);
        }
        //writeln!(filebuf5, "bcs");
        for i in 0..self.bcs.len()/3{

            if self.bcs[i*3] {
                write!(filebuf5, "1,");
            }
            else{
                write!(filebuf5, "0,");
            }
            if self.bcs[i*3+1] {
                write!(filebuf5, "1,");
            }
            else{
                write!(filebuf5, "0,");
            }
            if self.bcs[i*3+2] {
                writeln!(filebuf5, "1");
            }
            else{
                writeln!(filebuf5, "0");
            }

        }

        for i in 0..self.strains.len(){
            writeln!(filebuf6,"{},{},{},{},{},{}",self.e_pl[i][0],self.e_pl[i][1],self.e_pl[i][2],self.e_pl[i][3],self.e_pl[i][4],self.e_pl[i][5]);
        }
    }
}


fn main(){
    let dx = 12.7e-3; let dy = 12.7e-3; let dz = 5.08e-3;
    let mut stm = StructuralModel::with_capacity(60, dx, dy, dz,1200.0);
    let mut nodelist = Vec::with_capacity(12);
    nodelist.push(Node{index:16, pt:Point{x:0.0,y:0.0,z:-1.0}});
    nodelist.push(Node{index:0, pt:Point{x:0.0,y:0.0,z:0.0}});
    nodelist.push(Node{index:1, pt:Point{x:1.0*dx,y:0.0,z:0.0}});
    nodelist.push(Node{index:2, pt:Point{x:1.0*dx,y:1.0*dy,z:0.0}});
    nodelist.push(Node{index:3, pt:Point{x:0.0,y:1.0*dy,z:0.0}});
    nodelist.push(Node{index:4, pt:Point{x:0.0,y:0.0,z:1.0*dz}});
    nodelist.push(Node{index:5, pt:Point{x:1.0*dx,y:0.0,z:1.0*dz}});
    nodelist.push(Node{index:6, pt:Point{x:1.0*dx,y:1.0*dy,z:1.0*dz}});
    nodelist.push(Node{index:7, pt:Point{x:0.0,y:1.0*dy,z:1.0*dz}});

    nodelist.push(Node{index:8, pt:Point{x:2.0*dx,y:0.0,z:0.0}});
    nodelist.push(Node{index:9, pt:Point{x:2.0*dx,y:1.0*dy,z:0.0}});
    nodelist.push(Node{index:10, pt:Point{x:2.0*dx,y:0.0,z:1.0*dz}});
    nodelist.push(Node{index:11, pt:Point{x:2.0*dx,y:1.0*dy,z:1.0*dz}});

    nodelist.push(Node{index:12, pt:Point{x:3.0*dx,y:0.0*dy,z:0.0*dz}});
    nodelist.push(Node{index:13, pt:Point{x:3.0*dx,y:1.0*dy,z:0.0*dz}});
    nodelist.push(Node{index:14, pt:Point{x:3.0*dx,y:0.0*dy,z:1.0*dz}});
    nodelist.push(Node{index:15, pt:Point{x:3.0*dx,y:1.0*dy,z:1.0*dz}});

    nodelist.push(Node{index:16, pt:Point{x:4.0*dx,y:0.0*dy,z:0.0*dz}});
    nodelist.push(Node{index:17, pt:Point{x:4.0*dx,y:1.0*dy,z:0.0*dz}});
    nodelist.push(Node{index:18, pt:Point{x:4.0*dx,y:0.0*dy,z:1.0*dz}});
    nodelist.push(Node{index:19, pt:Point{x:4.0*dx,y:1.0*dy,z:1.0*dz}});

    nodelist.push(Node{index:20, pt:Point{x:5.0*dx,y:0.0*dy,z:0.0*dz}});
    nodelist.push(Node{index:21, pt:Point{x:5.0*dx,y:1.0*dy,z:0.0*dz}});
    nodelist.push(Node{index:22, pt:Point{x:5.0*dx,y:0.0*dy,z:1.0*dz}});
    nodelist.push(Node{index:23, pt:Point{x:5.0*dx,y:1.0*dy,z:1.0*dz}});

    nodelist.push(Node{index:24,pt:Point{x:2.0*dx,y:0.0*dy,z:2.0*dz}});
    nodelist.push(Node{index:25,pt:Point{x:3.0*dx,y:0.0*dy,z:2.0*dz}});
    nodelist.push(Node{index:26,pt:Point{x:3.0*dx,y:1.0*dy,z:2.0*dz}});
    nodelist.push(Node{index:27,pt:Point{x:2.0*dx,y:1.0*dy,z:2.0*dz}});

    //let mut bcs = vec![[false;24];5];
    //let bcnodes = vec![0,3,4,7,20,21,22,23];


    stm.add_element(&nodelist, [1,2,3,4,5,6,7,8],&[true, true, true, false, false, false, false, false, false, true, true, true, true, true, true, false, false, false,false, false, false,true, true, true] , true);
    stm.add_element(&nodelist, [2,9,10,3,6,11,12,7], &[false;24], true);
    stm.add_element(&nodelist, [9,13,14,10,11,15,16,12], &[false;24], true);
    stm.add_element(&nodelist, [13,17,18,14,15,19,20,16], &[false;24], true);
    stm.add_element(&nodelist, [17,21,22,18,19,23,24,20], &[false, false, false, true, true,true, true, true,true, false, false, false, false, false, false, true, true, true, true, true,true, false, false, false], true);
    stm.add_element(&nodelist,[11,15,16,12,25,26,27,28], &[false;24], true);



    let pool = ThreadPoolBuilder::new().num_threads(6).build().expect("can't build the threadpool");
    let E = vec![1e8;stm.elements.len()];
    let a = 0.577350269189626_f64;
    let Wa = 1.0_f64;
    let quadraturevals = [[-a, -a, -a, Wa],
        [-a, -a, a, Wa],
        [-a, a, -a, Wa],
        [-a, a, a, Wa],
        [a, -a, -a, Wa],
        [a, -a, a, Wa],
        [a, a, -a, Wa],
        [a, a, a, Wa]
    ];
    let nu = 0.4;
    let mut C2 = [
        [1_f64 - nu, nu, nu, 0., 0., 0.],
        [nu, 1. - nu, nu, 0., 0., 0.],
        [nu, nu, 1. - nu, 0., 0., 0.],
        [0., 0., 0., (1. - 2. * nu) / 2., 0., 0.],
        [0., 0., 0., 0., (1. - 2. * nu) / 2., 0.],
        [0., 0., 0., 0., 0., (1. - 2. * nu) / 2.],
    ];

    let mut tmp2 = vec![vec![[0.0_f64;24];24];12];
    let mut tmp1 = vec![vec![[0.0_f64;6];24];12];
    let mut Ctmp = vec![vec![[0.0_f64;6];6];24];
    let dxdydz = [12.7e-3, 12.7e-3, 5.08e-3];
    stm.isactive = vec![true;stm.elements.len()];
    StructuralModel::generate_local_stiffness_matrix_par(stm.Kloc.borrow_mut(),E.borrow(),0.40,stm.elements.borrow(),stm.nodes.borrow(),
                                                         stm.Bmat.borrow_mut(),&C2,&mut tmp1, &mut tmp2,  Ctmp.borrow_mut(),&dxdydz, &stm.isactive, stm.volume.borrow_mut(),&quadraturevals,stm.elemcount.borrow(),1,6,&pool);

    stm.logsig = true; // tell me how may iterations it required to converge

    //for i in 0..stm.bcs.len()/3{
    //   println!("{},{},{}", stm.bcs[i*3], stm.bcs[i*3+1], stm.bcs[i*3+2]);
    //}
    println!("number of bcs in stm {} num dof {}", stm.nbcs, stm.nodes.len()*3);
    let selfwtlds = vec![1.0;stm.elements.len()];
    stm.calc_self_wt_loads(&selfwtlds);
    stm.generate_global_stiffness_matrix_with_bc_rowwise(stm.nbcs);
    stm.structural_analysis_solve(1,6,&pool);
    println!("bcs len {}", stm.bcs.len());



    stm.isactive[5] = false;
    stm.deactivate_elements();
    println!("bcs len {}", stm.bcs.len());
    // stm.isactive[5] = true;
    stm.deactivate_elements();
    // for i in 0..stm.bcs.len()/3{
    //   println!("{}->{},{},{}",i, stm.bcs[i*3],stm.bcs[i*3+1],stm.bcs[i*3+2]);
    // }
    //println!("Kmat data len{}",stm.Kmat.data.len());
    //stm.generate_global_stiffness_matrix_with_bc(stm.nbcs);
    //  println!("stm Kmat data{}",stm.Kmat.data.len());
    stm.reset(1,6,&pool);
    stm.calc_self_wt_loads(&selfwtlds);

    // stm.generate_global_stiffness_matrix_with_bc_rowwise(stm.nbcs);
    stm.generate_global_stiffness_matrix_with_bc_rowwise_par(stm.nbcs,6, &pool);
    stm.structural_analysis_solve(1,6,&pool);
    println!("self nbcs, {}", stm.nbcs);

    //  for i in 0..stm.bcs.len()/3{
    //     println!("{}, {}, {}",stm.bcs[i*3],stm.bcs[i*3+1],stm.bcs[i*3+2]);
    //  }

    /*  println!("strain in elements");
      for i in 0..stm.elements.len(){
          println!("element {} , {:?}", i, stm.calc_strain_for_element(i));
      }
  */
    println!("all element strains");
    stm.calc_strain_for_all_elements();
    for i in 0..stm.elements.len(){
        println!("element {} , {:?}", i, stm.strains[i]);
    }


    for i in 0..stm.u.len()/3{
        println!("[{},{},{}]",stm.u[i*3],stm.u[i*3+1],stm.u[i*3+2]);
    }
    println!("Hello world!");
    let mut a = StructuralModel::with_capacity(1000,1.0,1.0,1.0,1200.0);


}


