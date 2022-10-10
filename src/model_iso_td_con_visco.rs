use crate::primitives::Node;
use std::collections::HashMap;
use crate::interpolator::Interpolator;
use rayon::ThreadPool;
use std::io::BufWriter;
use std::fs::File;
use crate::model_updater::ModelUpdater;
use std::io::Write;
use std::borrow::{BorrowMut, Borrow};
use crate::structural_model2::StructuralModel;

const STEFAN:f64 = 5.6703e-8;
const CtoKadd:f64 = 273.15;

/// Model has list of Nodes, list of cells, list of temperature for each cell, and information about
/// six neighbors (other cell or free surface) for each cell
pub struct Model{
    nodelist:Vec<Node>,
    sp_heat_cap:Vec<f64>,
    density:Vec<f64>,
    kx:Vec<f64>,
    pub(crate) templist:Vec<f64>,
    has_neighbor:Vec<[bool;6]>,
    neighbors:Vec<[usize;6]>,
    pub(crate) temp_at_contact:Vec<[f64;6]>,
    density_vol_sp_ht_cap:Vec<f64>,
    pub(crate) cell_index_map:HashMap<isize,isize>,
    sp_heat_cap_interp:Interpolator,
    element_width:f64,
    element_height:f64,
    active_layer_first_element:usize,
    previous_active_layer:usize,
    previous_layer:usize,
    current_layer:usize,
    turn_off_layers_at:usize,
    first_element_of_each_layer:Vec<usize>,
    temp_pos_for_split:Vec<usize>,
    pub(crate) strain_time_dep:Vec<f64>,
    pub(crate) stress:Vec<f64>,
    strain_time_dep_prev:Vec<f64>,
    stress_prev:Vec<f64>,
    element_no:Vec<usize>,
    pub(crate) loads: Vec<f64>,
    pub(crate) stm:StructuralModel,
    gbcs:Vec<[bool;3]>,
    zmin:f64,
    Evec:Vec<f64>,
    tmp2 : Vec<Vec<[f64;24]>>,
    tmp1 : Vec<Vec<[f64;6]>>,
    Ctmp: Vec<Vec<[f64;6]>>,
    emissivity:f64,
}

/// methods for Model struct
impl Model{
    /// generates a new empty Model
    pub fn new(nodes:Vec<Node>,activation_cells_length:usize, sp_heat_cap_interp: Interpolator, element_width:f64, element_height:f64,
               turn_off_layers_at:usize, zmin: f64, density:f64, emissivity:f64, maxthreads: usize)->Model{
        println!("element width {} element height {}", element_width, element_height);
        let mut tmp1 = Vec::with_capacity(maxthreads*4);
        let mut tmp2 = Vec::with_capacity(maxthreads*4);
        let mut Ctmp = Vec::with_capacity(maxthreads*4);
        for i in 0..maxthreads*4{
            tmp2.push( vec![[0.0_f64;24];24]);
            tmp1.push(vec![[0.0_f64;6];24]);
            Ctmp.push(vec![[0.0;6];6]);
        }
        return Model{
            nodelist: nodes,
            sp_heat_cap:Vec::with_capacity(activation_cells_length),
            density:Vec::with_capacity(activation_cells_length),
            kx:Vec::with_capacity(activation_cells_length),
            templist: Vec::with_capacity(activation_cells_length),
            neighbors: Vec::with_capacity(activation_cells_length),
            has_neighbor:Vec::with_capacity(activation_cells_length),
            temp_at_contact: Vec::with_capacity(activation_cells_length),
            density_vol_sp_ht_cap:Vec::with_capacity(activation_cells_length),
            cell_index_map:HashMap::with_capacity(activation_cells_length),
            sp_heat_cap_interp,
            element_width,
            element_height,
            active_layer_first_element:0,
            previous_active_layer:0,
            previous_layer: 0,
            current_layer: 0,
            turn_off_layers_at,
            first_element_of_each_layer: vec![0],
            temp_pos_for_split: Vec::with_capacity(activation_cells_length),
            strain_time_dep: Vec::with_capacity(activation_cells_length),
            strain_time_dep_prev: Vec::with_capacity(activation_cells_length),
            stress: Vec::with_capacity(activation_cells_length),
            stress_prev: Vec::with_capacity(activation_cells_length),
            loads: Vec::with_capacity(activation_cells_length),
            element_no:Vec::with_capacity(activation_cells_length),
            stm: StructuralModel::with_capacity(activation_cells_length,element_width,element_width,element_height, density),
            gbcs: Vec::with_capacity(activation_cells_length),
            zmin,
            Evec: Vec::with_capacity(activation_cells_length),
            tmp1,
            tmp2,
            Ctmp,
            emissivity
        }
    }
    /// adds a cell to the model. Takes node index, specific heat capacity, density, conductivity in
    /// the x-direction, conductivity in the y and z-directions, initial temperature, and neighbor
    /// information of the cell
    pub fn addCell(&mut self,node_no:[usize;8],sp_heat_cap:f64,
                   density:f64, kx:f64, ky:f64, kz:f64, init_temp:f64,is_neighbor:[bool;6],
                   neighbor:[isize;6],cell_index:usize,temporary_templist:&mut Vec<f64>, temporary_Glist:&mut Vec<f64>,
                   temporary_nlist: &mut Vec<f64>, temporary_tdlist: &mut Vec<f64>,
                   orientation:[f64;2],layer:usize,cell_indices:&mut Vec<usize>)
    {
        let mut counter = self.sp_heat_cap.len();
        self.element_no.push(counter);
        self.cell_index_map.insert(cell_index as isize,counter as isize);
        self.Evec.push(1e8);
        let mut bcs = [false;24];
        for kk in 0..8 {

            let zval = self.nodelist[node_no[kk]].pt.z;
            if (zval-self.zmin).abs() < 1e-7{
                bcs[kk*3] = true;
                bcs[kk*3+1] = true;
                bcs[kk*3+2] = true;
            }
        }
        let elem = [node_no[0],node_no[1],node_no[3],node_no[2],node_no[4],node_no[5],node_no[7],node_no[6]];
        /*println!("zmin, {}", self.zmin);
        println!("element_coord start");
        if self.templist.len()<27{
            for i in 0..8{
                println!("{},{},{}", &self.nodelist[elem[i]].pt.x,&self.nodelist[elem[i]].pt.y,&self.nodelist[elem[i]].pt.z);
                println!("{},{},{}",&bcs[i*3],&bcs[i*3+1],&bcs[i*3+2]);
            }
        }
        println!("element_coord end");*/



        cell_indices.push(counter);
        if (layer != self.current_layer){
            self.previous_layer = self.current_layer;
            self.current_layer = layer;
            self.first_element_of_each_layer.push(counter);
            //println!("layer {}",layer);
        }
        if (layer-self.previous_active_layer)>self.turn_off_layers_at{
            self.previous_active_layer = self.previous_active_layer+1;
            self.active_layer_first_element = self.first_element_of_each_layer[self.previous_active_layer+1];

        }

        let p1 = self.nodelist[node_no[0]].clone();
        let p2 = self.nodelist[node_no[7]].clone();

        let volume = self.element_width*self.element_width*self.element_height;//((p2.pt.x-p1.pt.x)*(p2.pt.y-p1.pt.y)*(p2.pt.z-p1.pt.z)).abs();

        self.density_vol_sp_ht_cap.push(density*volume);
        temporary_templist.push(init_temp.clone());
        temporary_Glist.push(0.0);
        temporary_nlist.push(0.0);
        temporary_tdlist.push(0.0);
        self.strain_time_dep.push(0.0);
        self.strain_time_dep_prev.push(0.0);
        self.stress.push(0.0);
        self.stress_prev.push(0.0);



        // please check this one more time
        self.kx.push(kx);
        // self.ky.push(ky);

        //self.kz.push(kz);
        self.density.push(density);
        self.sp_heat_cap.push(sp_heat_cap);
        // println!("start suspicious code");

        //update neighbor index
        let mut neighbor2:[usize;6] =[0,0,0,0,0,0];
        self.temp_at_contact.push([-100.0;6]);
        for i in 0..6{
            if is_neighbor[i]{
                neighbor2[i] = self.cell_index_map.get(&(neighbor[i])).expect("can't find it").clone() as usize;
            }
        }

        self.templist.push(init_temp);
        if is_neighbor[0] {
            self.neighbors[neighbor2[0] as usize][1] = counter;
            self.has_neighbor[neighbor2[0] as usize][1] = true;
            self.temp_at_contact[neighbor2[0] as usize][1] = self.templist[neighbor2[0] as usize];
            //self.temp_at_contact[neighbor2[0] as usize][1] = self.templist[neighbor2[0] as usize][1]
            //  self.has_neighbor_f64[neighbor2[0] as usize][1] = 1.0;
            //  self.has_no_neighbor_f64[neighbor2[0] as usize][1] = 0.0;
        }
        //println!("check 1 pass");
        if is_neighbor[1] {
            self.neighbors[neighbor2[1] as usize][0]=counter;
            self.temp_at_contact[neighbor2[1] as usize][0] = self.templist[neighbor2[1] as usize];
            self.has_neighbor[neighbor2[1] as usize][0]=true;
            // self.has_neighbor_f64[neighbor2[1] as usize][0] = 1.0;
            //self.has_no_neighbor_f64[neighbor2[1] as usize][0] = 0.0;
        }
        //  println!("check 2 pass");
        if is_neighbor[2] {
            self.neighbors[neighbor2[2] as usize][3]=counter ;
            self.temp_at_contact[neighbor2[2] as usize][3] = self.templist[neighbor2[2] as usize];
            self.has_neighbor[neighbor2[2] as usize][3]=true;
            // self.has_neighbor_f64[neighbor2[2] as usize][3] = 1.0;
            // self.has_neighbor_f64[neighbor2[2] as usize][3] = 0.0;
        }
        // println!("check 3 pass");
        if is_neighbor[3] {
            self.neighbors[neighbor2[3] as usize][2]=counter ;
            self.temp_at_contact[neighbor2[3] as usize][2] = self.templist[neighbor2[3] as usize];
            self.has_neighbor[neighbor2[3] as usize][2] = true;
            //self.has_neighbor_f64[neighbor2[3] as usize][2] = 1.0;
            //self.has_neighbor_f64[neighbor2[3] as usize][2] = 0.0;
        }
        // println!("check 4 pass");
        if is_neighbor[4] {
            self.neighbors[neighbor2[4] as usize][5]=counter;
            self.temp_at_contact[neighbor2[4] as usize][5] = self.templist[neighbor2[4] as usize];
            self.has_neighbor[neighbor2[4] as usize][5]=true;
        }

        if is_neighbor[5] {
            self.neighbors[neighbor2[5] as usize][4]= counter;
            self.temp_at_contact[neighbor2[5] as usize][4] = self.templist[neighbor2[5] as usize];
            self.has_neighbor[neighbor2[5] as usize][4] = true;
        }

        self.has_neighbor.push(is_neighbor);
        self.neighbors.push(neighbor2);
        self.stm.add_element(&self.nodelist,elem, &bcs, self.has_neighbor[counter][5]);


        //go through all the active layers, add weights to cells directly below this cell being deposited
        self.loads.push(-volume * &self.density[counter] * 9.81); //9.81 is g
        //self.loads.push(0.0);
        let mut endcell = false;
        let mut cellcount = 0;
        let mut botcell:usize = 0;
       /* if !self.has_neighbor[counter][5] {
            endcell=true
        }else{
            botcell = self.neighbors[counter][5] as usize;
        }
        while !endcell && cellcount < self.turn_off_layers_at{
            self.loads[botcell] -= volume * self.density[botcell] * 9.81; //g is 9.81

            if !self.has_neighbor[botcell][5] {
                endcell=true
            }else{
                botcell = self.neighbors[botcell][5] as usize;
            }
        }
        /*println!("layer {}",self.current_layer);
        if self.loads.len() <1000 && self.current_layer != self.previous_layer {
            for i in 0..self.loads.len() {
                println!("load,{},{}", i, self.loads[i]);
            }
        }*/*/

   }


    /// calculates new cell temperatures for the model for a timestep. Takes time increment,
    /// convection coefficient and environment temperature as input. Gives out temperature of the
    /// cells as output

    pub fn find_new_cell_temp_all(&self, dt:f64, conv_coeff:f64, t_env:f64, side_area:f64, top_area:f64, sa_d_sd:f64, ta_d_td:f64,
                                  vol:f64,pool:&ThreadPool,
                                  newtemp:&mut [f64],numthreads:usize, maxthreads:usize, cell_indices_slice: &[usize]) {
        if (numthreads < maxthreads){
            let newtempsize = newtemp.len();
            let (nt0,nt1) = newtemp.split_at_mut(newtempsize/2);
            let (cellpos0, cellpos1) = cell_indices_slice.split_at(newtempsize/2);
            pool.install(
                ||rayon::join(
                    || self.find_new_cell_temp_all(dt, conv_coeff,t_env,side_area,top_area, sa_d_sd,ta_d_td,vol, pool, nt0,numthreads*2,maxthreads,cellpos0),
                    || self.find_new_cell_temp_all(dt, conv_coeff,t_env,side_area,top_area, sa_d_sd,ta_d_td,vol, pool, nt1,numthreads*2,maxthreads,cellpos1)
                ));

        }
        else{
            self.find_new_temp_all_split(dt, conv_coeff, t_env, newtemp, side_area, top_area, sa_d_sd, ta_d_td, vol, cell_indices_slice);
        }
    }



    ///finds cell temperature for values in a range
    pub fn find_new_temp_all_split(&self, dt:f64, conv_coeff:f64, t_env:f64, temp: &mut [f64], side_area:f64, top_area:f64, sa_d_sd:f64, ta_d_td:f64, vol:f64,
                                   cell_indices_slice:&[usize]){

        for cell in cell_indices_slice{
            let mut dQ = 0.0;
            let thistemp = self.templist[*cell] + CtoKadd;
            let T4 = thistemp*thistemp*thistemp*thistemp;
            let ta = t_env+CtoKadd;
            let ta4 = ta*ta*ta*ta;
            let radiation_loss_per_unit_area = STEFAN*self.emissivity*(T4-ta4);
            //front back
            //dist area and vol is  [ side area, top area, sidearea/sidedistance, toparea/topdistance, volume]
            for i in 0..2 {
                if !self.has_neighbor[*cell][i] {
                    dQ = dQ - (conv_coeff * (self.templist[*cell] - t_env)+radiation_loss_per_unit_area) * dt * side_area;// * self.has_no_neighbor_f64[*cell][i] ;
                } else {
                    let neighbor_index = self.neighbors[*cell][i];
                    let t1 = self.templist[*cell];
                    let t2 = self.templist[neighbor_index];
                    dQ = dQ - self.kx[*cell] * (t1 - t2) * sa_d_sd * dt;// *self.has_neighbor_f64[*cell][i];
                }
            }
            //right left
            for i in 2..4 {
                if !self.has_neighbor[*cell][i] {
                    dQ = dQ - (conv_coeff * (self.templist[*cell] - t_env) + radiation_loss_per_unit_area)* dt * side_area;// * self.has_no_neighbor_f64[*cell][i] ;
                } else {
                    let neighbor_index = self.neighbors[*cell][i];
                    let t1 = self.templist[*cell];
                    let t2 = self.templist[neighbor_index ];
                    dQ = dQ - self.kx[*cell] * (t1 - t2) * sa_d_sd * dt;// *self.has_neighbor_f64[*cell][i];
                }
            }

            //top down
            for i in 4..6 {
                if !self.has_neighbor[*cell][i] {
                    dQ = dQ - (conv_coeff * (self.templist[*cell] - t_env) + radiation_loss_per_unit_area) * dt * top_area;
                } else {
                    let neighbor_index = self.neighbors[*cell][i];
                    let t1 = self.templist[*cell];
                    let t2 = self.templist[neighbor_index ];
                    dQ = dQ - self.kx[*cell] * (t1 - t2) * ta_d_td * dt;
                }
            }
            let mult = &self.density_vol_sp_ht_cap[*cell] * self.sp_heat_cap[*cell];


            temp[*cell-cell_indices_slice[0]] = (mult * self.templist[*cell] + dQ)/ mult;
            //let newtemperature = temp[*cell-cell_indices_slice[0]];
            /*if newtemperature>220.0 || newtemperature<25.0{
                println!("newtemperature {} *cell-cell_indices_slice[0] {} mult {}", newtemperature, *cell-cell_indices_slice[0], mult)
            }*/
        }

    }


    /// calculates the temperature of the cells in the model during a time period with a given time
    /// step. Takes time period, time step, and a file (buffered) where output is to be written as
    /// input. Writes output to the file (buffer). Might need to change that for speed.
    pub fn run_model(&mut self, time:f64, dt:f64, conv_coeff:f64, t_env:f64, file:&mut BufWriter<File>,
                     global_time:&mut f64,pool:&ThreadPool, areas_and_dists:[f64;5], newtemp: &mut Vec<f64>, temporary_Glist: &mut Vec<f64>, temporary_nlist: &mut Vec<f64>, temporary_tdlist: &mut Vec<f64>,
                     conductivity_interp:&Interpolator, maxthreads:usize, cell_indices:&Vec<usize>, vis_interpolators: [&Interpolator;3]) {
        let count = (time / dt) as usize;

        let (inactive_temps, active_temps) = newtemp.split_at_mut(self.active_layer_first_element);
        let (old_cell_indices_slice, new_cell_indices_slice) = cell_indices.split_at(self.active_layer_first_element);
        let new_cell_slice_start = new_cell_indices_slice[0].clone();
        let ncis_len = new_cell_indices_slice.len();
        let new_cell_slice_end = new_cell_indices_slice[ncis_len - 1].clone();
        {
            let (inactive_selftemps, active_selftemps) = self.templist.split_at(self.active_layer_first_element);
            let (sp_ht_cap_inactive, sp_ht_cap_active) = self.sp_heat_cap.split_at_mut(self.active_layer_first_element);
            let (tG_inactive, tG_active) = temporary_Glist.split_at_mut(self.active_layer_first_element);
            let (tn_inactive, tn_active) = temporary_nlist.split_at_mut(self.active_layer_first_element);
            let (ttd_inactive, ttd_active) = temporary_tdlist.split_at_mut(self.active_layer_first_element);
            //G interpolator
            vis_interpolators[0].interpolate_slice(active_selftemps, tG_active, pool);
            //n interpolator
            vis_interpolators[1].interpolate_slice(active_selftemps, tn_active, pool);
            //tan delta interpolator
            vis_interpolators[2].interpolate_slice(active_selftemps, ttd_active, pool);

            self.sp_heat_cap_interp.interpolate_slice(active_selftemps, sp_ht_cap_active, pool);

            let (inactive_conductivity, active_conductivity) = self.kx.split_at_mut(self.active_layer_first_element);
            conductivity_interp.interpolate_slice(active_selftemps, active_conductivity, pool);
        }

        //println!("time {}",time);
        for i in 0..count {
            self.find_new_cell_temp_all(dt, conv_coeff, t_env, areas_and_dists[0], areas_and_dists[1], areas_and_dists[2],
                                        areas_and_dists[3], areas_and_dists[4], &pool,
                                        active_temps, 1, maxthreads, new_cell_indices_slice);

            self.templist[new_cell_slice_start..new_cell_slice_end + 1].copy_from_slice(active_temps);
            //if (self.templist.len()>0) {write!(file, "{},{},{},{}\n", self.templist[0].clone(), global_time, self.sp_heat_cap[0], self.density_vol_sp_ht_cap[0]);}
            /*           let elemno = 1250;
            if (self.templist.len()>elemno) {write!(file, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},", self.templist[elemno].clone(), global_time, self.stm.e_el[elemno][2], self.stm.e_pl[elemno][2], self.stm.stress[elemno][0],self.stm.stress[elemno][1],self.stm.stress[elemno][2],self.stm.stress[elemno][3],self.stm.stress[elemno][4],self.stm.stress[elemno][5],self.stm.loads[self.stm.elements[elemno][0]*3+2], self.stm.uel[self.stm.elements[elemno][0]*3+2], temporary_nlist[elemno], temporary_Glist[elemno]);}
            for j in 0..8{
                if (self.templist.len()>elemno) {write!(file, "{},",self.stm.uel[self.stm.elements[elemno][j]*3+0]);}
                if (self.templist.len()>elemno) {write!(file, "{},",self.stm.uel[self.stm.elements[elemno][j]*3+1]);}
                if (self.templist.len()>elemno) {write!(file, "{},",self.stm.uel[self.stm.elements[elemno][j]*3+2]);}
            }
            if (self.templist.len()>elemno) {write!(file, "\n");};*/
        }
        self.find_new_cell_temp_all(time - count as f64 * dt, conv_coeff, t_env, areas_and_dists[0], areas_and_dists[1], areas_and_dists[2],
                                    areas_and_dists[3], areas_and_dists[4], &pool,
                                    active_temps, 1, maxthreads, new_cell_indices_slice);
        let elemno = 1444;
        if (self.templist.len() > elemno) {
            write!(file, "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}", self.templist[elemno].clone(), global_time, self.stm.e_el[elemno][2], self.stm.e_pl[elemno][2], self.stm.stress[elemno][0], self.stm.stress[elemno][1], self.stm.stress[elemno][2], self.stm.stress[elemno][3], self.stm.stress[elemno][4], self.stm.stress[elemno][5], self.stm.loads[self.stm.elements[elemno][0] * 3 + 2], self.stm.uel[self.stm.elements[elemno][0] * 3 + 2], temporary_nlist[elemno], temporary_Glist[elemno], self.stm.volume[elemno]);

        for i in 0..8 {
            { write!(file, "{},", self.stm.uel[self.stm.elements[elemno][i] * 3 + 0]); }
            { write!(file, "{},", self.stm.uel[self.stm.elements[elemno][i] * 3 + 1]); }
            { write!(file, "{},", self.stm.uel[self.stm.elements[elemno][i] * 3 + 2]); }
        }
            for i in 0..8{
                write!(file,"{},",self.stm.loads[self.stm.elements[elemno][i]*3+0]);
                write!(file,"{},",self.stm.loads[self.stm.elements[elemno][i]*3+1]);
                write!(file,"{},",self.stm.loads[self.stm.elements[elemno][i]*3+2]);
            }

       // { write!(file, "{}", self.templist[self.neighbors[elemno][4]]); }//"{},{},{},{},{},{}\n", self.templist[self.neighbors[elemno][0]],self.templist[self.neighbors[elemno][1]],self.templist[self.neighbors[elemno][2]], self.templist[self.neighbors[elemno][3]], self.templist[self.neighbors[elemno][4]], self.templist[self.neighbors[elemno][5]]);};
        { write!(file, "\n"); }
    }

        let (tmplsti, tmplsta) = self.templist.split_at_mut(self.active_layer_first_element);
        let (tmpi, tmpa) = active_temps.split_at(self.active_layer_first_element);
        Model::parallel_copy(tmplsta, tmpa, 1, maxthreads, pool);

        //reset load_bcs and local stiffness matrices
        self.stm.reset(1,maxthreads,pool);

        //check if element is active -> if element has a tandelta of >1.0
       // for i in 0..self.templist.len(){
       //     if temporary_tdlist[i] < 1.0{
        //        self.stm.isactive[i] = false;
        //    }
        //    else{
        //        self.stm.isactive[i] = true;
       //     }
     //   }
        Model::parallel_active_check(&temporary_tdlist,temporary_Glist, temporary_nlist, self.Evec.borrow_mut(), self.stm.isactive.borrow_mut(), time, 2.8, 1.0,self.has_neighbor.borrow(), 1,maxthreads, pool);
       // self.stm.deactivate_elements();
       // apply bcs to element nodes, compare with structural bcs (if structural bcs true, then structural bcs take priority)



        //calculate loads
        self.stm.calc_self_wt_loads(self.loads.borrow());
       // let mut Evec = vec![0.0;temporary_Glist.len()];
        //for i in 0..self.Evec.len() {
        //    self.Evec[i] = temporary_Glist[i] * (1.0 + f64::exp(-time / (temporary_nlist[i] / temporary_Glist[i]))) * 2.8;

       // }
        self.stm.deactivate_elements();


        //
        let nu = 0.3;
        let mut C2 = [
            [1_f64 - nu, nu, nu, 0., 0., 0.],
            [nu, 1. - nu, nu, 0., 0., 0.],
            [nu, nu, 1. - nu, 0., 0., 0.],
            [0., 0., 0., (1. - 2. * nu) / 2., 0., 0.],
            [0., 0., 0., 0., (1. - 2. * nu) / 2., 0.],
            [0., 0., 0., 0., 0., (1. - 2. * nu) / 2.],
        ];
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
       // let mut tmp2 = [[0.0_f64;24];24];
       // let mut tmp1 = [[0.0_f64;6];24];


        StructuralModel::generate_local_stiffness_matrix_par(self.stm.Kloc.borrow_mut(), self.Evec.borrow(),nu,self.stm.elements.borrow(), self.stm.nodes.borrow(),self.stm.Bmat.borrow_mut(), &C2,
                                                             self.tmp1.borrow_mut(),  self.tmp2.borrow_mut(), self.Ctmp.borrow_mut(),
                                                             &[self.stm.dx, self.stm.dy, self.stm.dz], &self.stm.isactive,
                                                              self.stm.volume.borrow_mut(), &quadraturevals, self.stm.elemcount.borrow(), 1, maxthreads, pool);
        //self.stm.generate_local_stiffness_matrix(&Evec, 0.3);

        //strain calculations
        //println!("loads,{:?}",self.loads);
        self.stm.generate_global_stiffness_matrix_with_bc_rowwise_par(self.stm.nbcs,maxthreads,pool);
        self.stm.structural_analysis_solve(1,maxthreads,pool);
        //StructuralModel::calc_strain_for_all_elements_par(&self.stm.u,&self.stm.elements,self.stm.Bmat.as_slice(),self.stm.strains.borrow_mut(), &self.stm.isactive, 1, maxthreads, pool );
        StructuralModel::calc_stress_for_all_elements_par(&self.stm.uel,&self.stm.elements, &self.stm.Bmat, &mut self.stm.stress, &mut self.stm.strains, &mut self.stm.e_pl, &mut self.stm.e_el, &mut self.stm.e_pl0,
        time,&self.stm.isactive,&temporary_Glist,&nu,&temporary_nlist,1,maxthreads,pool);

        //if *global_time<10.0 {println!("time {}, global time {}", time, *global_time)}
        //println!("stmloads 1:10 {:?}", (0..10).map (|i| self.stm.loads[i]).collect::<Vec<f64>>());
        //self.stm.reset(1, maxthreads, &());
       /* println!("loads");
        let var = &self.stm.loads;
        for i in 0..var.len()/3{
            println!("[{}, {}, {}]", var[i*3],var[i*3+1],var[i*3+2]);
        }*/
       // println!("displacements");
       // let var = &self.stm.u;
       /* for i in 0..var.len()/3{
            println!("[{:.2e}, {:.2e}, {:.2e}]", var[i*3],var[i*3+1],var[i*3+2]);
        }*/
        /*if self.templist.len()<28 {
            self.stm.write("c:\\rustfiles\\feainputfile.csv");
        }*/

        //loop through elements
        //println!("no elems {}", self.templist.len());
    }

    fn parallel_active_check(tdlist:&[f64],  Glist: &[f64], nlist:&[f64], Evec: &mut[f64],isactive: &mut[bool], time: f64, two_onepnu: f64, checkvalue:f64,
                             has_neighbors: &[[bool;6]], numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        if numthreads<maxthreads{
            let splitpos = tdlist.len()/2;
            let (td1, td2) = tdlist.split_at(splitpos);
            let (G1, G2) = Glist.split_at(splitpos);
            let (n1,n2 ) = nlist.split_at(splitpos);
            let (E1, E2) = Evec.split_at_mut(splitpos);
            let (hn1, hn2) = has_neighbors.split_at(splitpos);
            let (isactive1, isactive2) = isactive.split_at_mut(splitpos);
            pool.install(||rayon::join(
                ||Model::parallel_active_check(td1,G1, n1, E1, isactive1,time, two_onepnu, checkvalue,hn1,numthreads*2, maxthreads, pool),
                ||Model::parallel_active_check(td2,G2,n2, E2,isactive2,time, two_onepnu,checkvalue,hn2,numthreads*2, maxthreads, pool)

            ));
        }
        else{
            for i in 0..tdlist.len(){
                if tdlist[i] > checkvalue && has_neighbors[i][5] {
                    isactive[i] = true;

                }
                else{
                    isactive[i] = false
                }
                Evec[i] = Glist[i] * two_onepnu;

            }
        }
    }



        fn parallel_copy(valold: &mut [f64], valnew: &[f64], numthreads: usize, maxthreads: usize, pool: &ThreadPool)
        {
            if numthreads < maxthreads {
                let splitpos = valold.len() / 2;
                let (o1, o2) = valold.split_at_mut(splitpos);
                let (n1, n2) = valnew.split_at(splitpos);

                pool.install(|| rayon::join(
                    || Model::parallel_copy(o1, n1, numthreads * 2, maxthreads, pool),
                    || Model::parallel_copy(o2, n2, numthreads * 2, maxthreads, pool)
                ));
            } else {
                valold.copy_from_slice(valnew);
            }
        }

        fn update_strains(strain_time_dep: &mut [f64], stress: &mut [f64],
                          strain_time_dep_prev: &[f64], stress_prev: &[f64], element_no: &[usize],
                          neighbors: &Vec<[usize; 6]>, has_neighbor: &Vec<[bool; 6]>,
                          dx: f64, dy: f64, dz: f64, density: f64,
                          u: &mut [f64], u_prev: &[f64], u_dot: &mut [f64], u_dot_prev: &[f64],
                          temporary_nlist: &[f64], temporary_Glist: &[f64], temporary_tdlist: &[f64],
                          loads: &[f64], time: &f64, area: &f64, num_threads: usize, max_threads: usize, pool: &ThreadPool)
        {
            if num_threads < max_threads {
                let split_loc = strain_time_dep.len() / 2;
                let (std1, std2) = strain_time_dep.split_at_mut(split_loc);
                //let (stdp1, stdp2) = strain_time_dep_prev.split_at_mut(std_len/2);
                let (strs1, strs2) = stress.split_at_mut(split_loc);
                //let (strsp1, strsp2) = stress_prev.split_at_mut(std_len/2);
                let (tnl1, tnl2) = temporary_nlist.split_at(split_loc);
                let (tGl1, tGl2) = temporary_Glist.split_at(split_loc);
                let (ttdl1, ttdl2) = temporary_tdlist.split_at(split_loc);
                let (l1, l2) = loads.split_at(split_loc);
                let (elno1, elno2) = element_no.split_at(split_loc);

                let (u1, u2) = u.split_at_mut(split_loc);
                //let (up1, up2) = u_prev.split_at(split_loc);
                let (ud1, ud2) = u_dot.split_at_mut(split_loc);
                let (udp1, udp2) = u_dot_prev.split_at(split_loc);

                pool.install(|| rayon::join(
                    || Model::update_strains(std1, strs1, strain_time_dep_prev, stress_prev, elno1, neighbors, has_neighbor, dx, dy, dz, density,
                                             u1, u_prev, ud1, udp1, tnl1, tGl1, ttdl1, l1, time, area, num_threads * 2, max_threads, pool),
                    || Model::update_strains(std2, strs2, strain_time_dep_prev, stress_prev, elno2, neighbors, has_neighbor, dx, dy, dz, density,
                                             u2, u_prev, ud2, udp2, tnl2, tGl2, ttdl2, l2, time, area, num_threads * 2, max_threads, pool)
                ));
            } else {
                for i in 0..strain_time_dep.len() {
                    let current_neighbor = neighbors[element_no[i]];
                    let current_has_neighbor = has_neighbor[element_no[i]];
                    let mut sig = [0.0; 6];// this stores disps of neighbors, if no neighbors, stores 0
                    let mut us = [0.0; 6];
                    //let sf = 0.0; let sb = 0.0; let st = 0.0; let sd = 0.0; let sl = 0.0; let sr = 0.0;
                    for j in 0..6 {
                        if current_has_neighbor[j] {
                            sig[j] = stress_prev[current_neighbor[j]];
                        }
                    }
                    for j in 4..6 {
                        if current_has_neighbor[j] {
                            us[j] = u_prev[current_neighbor[j]];
                        }
                    }

                    if temporary_tdlist[i] > 1.0 {
                        //calculate strainz from uz
                        strain_time_dep[i] = (us[5] - 2.0 * u_prev[element_no[i]] + us[4]) / (2.0 * dz);
                        let n = temporary_nlist[i];
                        let E = temporary_Glist[i];
                        let et = strain_time_dep[i];
                        let etm1 = strain_time_dep_prev[element_no[i]];
                        let s = stress_prev[element_no[i]];
                        //calculate new stress for next step
                        stress[i] = et * E;//(n*(et-etm1)+s/E)/(time/n*1.0/E);
                        //equilibrium of previous stresses
                        let rhs = (sig[0] - 2.0 * s + sig[1]) / (2.0 * dy * dz) + (sig[2] - 2.0 * s + sig[3]) / (2.0 * dx * dz) + (sig[4] - 2.0 * s + sig[5]) / (2.0 * dx * dy) + (-loads[i] / (dx * dy * dz));
                        let acc = rhs / density;
                        u_dot[i] = (u_dot_prev[i] + acc * time) * 0.001; //damping factor
                        u[i] = u_prev[element_no[i]] + u_dot[i] * time;
                        if i == 0 {
                            println!("{},{},{},{},{}", stress[i], rhs, acc, u_dot[i], u[i]);
                        }
                    }
                }
            }
        }


        /// updates model based on input. Needs to be placed in Model implementation later.
        pub fn update_model(mu: &mut ModelUpdater, mdl: &mut Model, bw: &mut BufWriter<File>, areas_and_dists: [f64; 5], time_step: f64,
                            conductivity_interp: &Interpolator, maxthreads: usize, input_data: ([f64; 3], f64, f64, f64, [f64; 2]), pool: &ThreadPool, vis_interpolators: [&Interpolator; 3]) {
            //let max_cpus = num_cpus::get_physical();

            //let pool = rayon::ThreadPoolBuilder::new().num_threads(maxthreads).build().unwrap();
            //lists for interpolation
            let mut temporary_templist = Vec::with_capacity(mu.activation_times.len());
            let mut temporary_Glist = Vec::with_capacity(mu.activation_times.len());
            let mut temporary_nlist = Vec::with_capacity(mu.activation_times.len());
            let mut temporary_tdlist = Vec::with_capacity(mu.activation_times.len());

            let mut cell_indices = Vec::with_capacity(mu.activation_times.len());
            //println!("activation times len {}" ,mu.activation_times.len());
            let mut next_cell_info = mu.get_next_cell_info();
            let mut next_cell_info2 = mu.get_next_cell_info();
            let mut global_time = 0.0;
            let time = next_cell_info2.0 - next_cell_info.0;
            let dt = time_step;
            let kx = input_data.0[0];
            let ky = input_data.0[1];
            let kz = input_data.0[2];
            let density = input_data.1;
            let init_temp = input_data.4[0];
            let t_env = input_data.4[1];
            let conv_coeff = input_data.3;
            let sp_heat_cap = input_data.2;
            mdl.addCell(next_cell_info.1, sp_heat_cap, density, kx, ky, kz, init_temp,
                        next_cell_info.3, next_cell_info.2, next_cell_info.4, &mut temporary_templist, &mut temporary_Glist, &mut temporary_nlist, &mut temporary_tdlist,
                        next_cell_info.5, next_cell_info.6, &mut cell_indices);
            //global_time = global_time + time;
            global_time = global_time + time;
            if time>1e-4 {
                mdl.run_model(time, dt, conv_coeff, t_env, bw, &mut global_time, &pool, areas_and_dists, &mut temporary_templist, &mut temporary_Glist, &mut temporary_nlist, &mut temporary_tdlist,
                              conductivity_interp, maxthreads, &cell_indices, vis_interpolators);
            }

            //println!("run model");
            let endval = mu.activation_times.len();
            for i in 1..endval - 1 {
                next_cell_info = next_cell_info2;
                next_cell_info2 = mu.get_next_cell_info();

                let time = next_cell_info2.0 - next_cell_info.0;
                // global_time = global_time + time;
                if i % 1000 == 0 { println!("{} out of {}, time {}", i, endval, global_time); mdl.stm.logsig = true}
                //println!("layer {}",next_cell_info.6);
                mdl.addCell(next_cell_info.1, sp_heat_cap, density, kx, ky, kz, init_temp,
                            next_cell_info.3, next_cell_info.2, next_cell_info.4, &mut temporary_templist,
                            &mut temporary_Glist, &mut temporary_nlist, &mut temporary_tdlist,
                            next_cell_info.5, next_cell_info.6, &mut cell_indices);
                //println!("second addcell");
                //println!("second run_model");
                global_time = global_time + time;
                if time > 1e-4 {
                    mdl.run_model(time, dt, conv_coeff, t_env, bw, &mut global_time, &pool, areas_and_dists, &mut temporary_templist,
                                  &mut temporary_Glist, &mut temporary_nlist, &mut temporary_tdlist,
                                  conductivity_interp, maxthreads, &cell_indices, vis_interpolators);
                }

                //println!("second run model");
            }
            //println!("pulupulu ");
            next_cell_info = next_cell_info2;

            mdl.addCell(next_cell_info.1, sp_heat_cap, density, kx, ky, kz, init_temp,
                        next_cell_info.3, next_cell_info.2, next_cell_info.4, &mut temporary_templist,
                        &mut temporary_Glist, &mut temporary_nlist, &mut temporary_tdlist,
                        next_cell_info.5, next_cell_info.6, &mut cell_indices);
            // global_time = global_time + time;
            let timedt = 2.0;
            /*for i in 0..(100.0 / timedt) as usize {
                global_time = global_time + timedt;
                mdl.run_model(timedt, dt, conv_coeff, t_env, bw, &mut global_time, &pool, areas_and_dists, &mut temporary_templist,
                              &mut temporary_Glist, &mut temporary_nlist, &mut temporary_tdlist,
                              conductivity_interp, maxthreads, &cell_indices, vis_interpolators);
            }*/
            println!("{} out of {}, time {}", endval, endval, global_time);
        }
    }
