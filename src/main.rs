mod deadcode;
mod gode_reader;
mod interpolator;
mod model;
mod model_generator;
mod model_iso_td_shc_td_con;
mod model_orthotropic_td_shc;
mod model_updater;
mod primitives;
mod simpleModel;
mod model_iso_td_con_visco;
mod structural_model2;
mod sparse_matrix_new;
//mod model_with_structural;

use num_cpus;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::borrow::BorrowMut;
use std::cmp::{max, Ordering};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter, Read};
use std::iter::FromIterator;
use std::ops::Index;
use std::thread::current;
use std::{fmt, fs};
//use nalgebra::{Matrix6, Vector6};
use crate::gode_reader::GCodeReader;
use crate::interpolator::*;
use crate::model::Model;
use crate::model_generator::ModelGenerator;
use crate::model_updater::ModelUpdater;
use crate::primitives::{Node, Point};
use std::str::FromStr;

use std::time::Instant;


fn main() {

    let mut maxthreads = num_cpus::get_physical()+2;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(maxthreads)
        .build()
        .expect("can't build the threadpool");
    let maxthreads = maxthreads*2;
    let gcr = GCodeReader::new("input/thorndike_seg1_solid.gcode");//_double_speed.gcode"); hollow_square_FEAtest_3ft_half.gcode//cube.gcode //fea_bowl_part_half.gcode //cylinder_hollow.gcode //thorndike_seg1_solid.gcode //mesh_conv_test.gcode
    //let gcr = GCodeReader::new("/mnt/c/rustfiles/thorndike_seg1_solid.gcode"); //cube.gcode //fea_bowl_part_half.gcode //cylinder_hollow.gcode //thorndike_seg1_solid.gcode //mesh_conv_test.gcode



    //let filename = "c:\\rustfiles\\gcodedump_med.csv";
    //gcr.write_abaqus_input(filename);
    println!("file read");

    let divs_per_bead = 3;
    let divs_per_bead_z = 3;
    let beadwidth = 12.7e-3;//12.7e-3; //8e-4;//12.7* 1e-3;
    let element_width = beadwidth / divs_per_bead as f64;
    let beadheight = 5.08e-3;//5.08e-3; //4e-4;//5.08 * 1e-3;
    let element_height = beadheight / divs_per_bead_z as f64;
    let xdiv = ((gcr.xmax - gcr.xmin + beadwidth) / element_width + 1.0).round() as usize;
    let ydiv = ((gcr.ymax - gcr.ymin + beadwidth) / element_width + 1.0).round() as usize;
    let zdiv = ((gcr.zmax - gcr.zmin + beadheight) / element_height + 1.0).round() as usize;

    let zmin = gcr.zmin.clone();
    println!(
        "xmin {} xmax {} ymin {} ymax {} zmin {} zmax {} xdiv {} ydiv {} zdiv {}",
        gcr.xmin, gcr.xmax, gcr.ymin, gcr.ymax, gcr.zmin, gcr.zmax, xdiv, ydiv, zdiv
    );
    let tic = Instant::now();
    let emissivity = 0.9;
    let convection_coefficient = 12.0;
    let input_file_speed = 46.065;
    let speed = 69.1;//input_file_speed* 60.0/input_file_speed;

    let init_temp = 200.0;
    let mut m = ModelGenerator::new(
            gcr.xmin - beadwidth / 2.0,
            gcr.xmax + beadwidth / 2.0,
            gcr.ymin - beadwidth / 2.0,
            gcr.ymax + beadwidth / 2.0,
            gcr.zmin - beadheight,
            gcr.zmax,
            xdiv,
            ydiv,
            zdiv,
            init_temp,
            element_width,
            divs_per_bead_z,
        );
    println!("time taken to generate model {:?}", tic.elapsed());

    let ndsize = m.nodelist.len();
    println!("mesh generated");

    // (to do) make a file buffer with input file data and bead parameters, calculate hash digest, compare with hash digest from previous run.
    // if hash digest is same, use activation times files, else calculate new activations


    let calc_new_act_times = true;
    let mut activation_times_all: Vec<(f64,usize,[f64;2],usize)> = Vec::with_capacity(10000000);


    let speed_mult = speed/input_file_speed;
    if calc_new_act_times {
        activation_times_all = m.generate_activation_times_all_layers(
            gcr.segment,
            gcr.is_extrusion_on,
            gcr.speed,
            beadwidth,
            beadheight,
            &pool,
            maxthreads,
        );

        for i in 0..activation_times_all.len(){
            activation_times_all[i].0 = activation_times_all[i].0/speed_mult;
        }
        let mut activation_times_file = File::create("output/activation_times_store_file.csv").expect("can't create file");
       // let mut activation_times_file = File::create("/mnt/c/rustFiles/activation_times_store_file.csv").expect("can't create file");
        let mut write_buffer = BufWriter::with_capacity(100000,activation_times_file);
        for i in activation_times_all.clone(){
            writeln!(write_buffer,"{},{},{},{},{}",i.0,i.1,i.2[0],i.2[1],i.3);
        }
    }
    else{
        let tic = Instant::now();
        let activation_times_file = File::open("output/activation_times_store_file.csv").expect("can't open file");
       // let activation_times_file = File::open("/mnt/c/rustFiles/activation_times_store_file.csv").expect("can't open file");
        let read_buffer = BufReader::with_capacity(100000, activation_times_file);
        for i in read_buffer.lines(){
            let line = i.expect("cant read line");
            let mut vals = line.split(",");
            let act_time = f64::from_str(vals.next().expect("no next")).expect("cant parse to f64");
            let el_no = usize::from_str(vals.next().expect("no next")).expect("cant parse to usize");
            let o1 = f64::from_str(vals.next().expect("no next")).expect("cant parse to f64");
            let o2 = f64::from_str(vals.next().expect("no next")).expect("cant parse to f64");
            let layer = usize::from_str(vals.next().expect("no next")).expect("cant parse to usize");
            activation_times_all.push((act_time, el_no, [o1,o2], layer));
        }
        println!("time taken to read activation times file {:?}",tic.elapsed());
    }

    let model_fraction = 1.0; //fraction of model to calculate ... for debugging purposes
    let activation_times:Vec<(f64,usize,[f64;2],usize)> = (0..(activation_times_all.len() as f32 *model_fraction) as usize ).map(|i| activation_times_all[i]).collect();

    let mut activated_elements = Vec::with_capacity(activation_times.len());
    let mut activation_times_only = Vec::with_capacity(activation_times.len());
    for i in 0..activation_times.len() {
        activated_elements.push(activation_times[i].1.clone());
        activation_times_only.push(activation_times[i].0.clone());
    }


    let cooldown_period = 60.0;
    let timestep = 0.002  ;

    let mut file1 = File::create("output/activation_times_all_layers.csv").expect("cant create the file");
    //let mut file1 = File::create("/mnt/c/rustfiles/activation_times_all_layers.csv").expect("cant create the file");
    let mut file1buf = BufWriter::with_capacity(10000, file1);
    for i in activation_times.clone() {
        write!(file1buf, "{},{},{}\n", i.0, i.1, i.3);
    }

    let mut activation_times_file = File::create("output/activation_times_with_cellno.csv").expect("cant open file");

    let mut activation_times_file_buf = BufWriter::with_capacity(10000, activation_times_file);
    //specific heat capacity table
    let sp_heat_cap_interpolation_table = Interpolator::read_data_from_file("input/sp_heat_cap_data.csv",1.0,maxthreads);
    let temps = vec![0.0, 50.0, 100.0, 150.0, 200.0, 250.0];
    let conductivity = vec![0.201; 6];
    let xstep = 50.0;
    let conductivity_interpolation_table =
        Interpolator::new(temps.clone(), conductivity.clone(), xstep, maxthreads);


    println!("specific heat capacity table read");
    for i in &activation_times {
        let mut center_point = Point {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let pts = m.celllist[i.1.clone()];
        for ii in &pts {
            center_point.x = center_point.x + m.nodelist[*ii].pt.x / 8.0;
            center_point.y = center_point.y + m.nodelist[*ii].pt.y / 8.0;
            center_point.z = center_point.z + m.nodelist[*ii].pt.z / 8.0;
        }
        write!(
            activation_times_file_buf,
            "{},{},{},{},{},{},{} \n",
            i.0.clone(),
            i.1.clone(),
            center_point.x,
            center_point.y,
            center_point.z,
            i.2.clone()[0],
            i.2.clone()[1]
        );
    }
    //let activation_times:Vec<(f64,usize)> = (0..m.celllist.len().clone()).map(|i|(i as f64,i as usize) ).collect();
    let mut mu = ModelUpdater::new(activation_times, m);
    let mut muclone = mu.clone();
    let input_data = ([0.205, 0.205, 0.205], 1210.0, 1500.0, convection_coefficient, [init_temp, 25.0]);


    let mut mdl = model_iso_td_con_visco::Model::new(
        mu.model_generator.nodelist.clone(),
        mu.activation_times.len(),
        sp_heat_cap_interpolation_table,
        element_width,
        element_height,
        12000,
        zmin-beadheight,input_data.1,emissivity,maxthreads);

    let mut file = File::create("output/unused_outuput.csv").expect("can't create file.");

    let mut bw = BufWriter::with_capacity(10000, file);

    //side area, top area, side_area/side_distance, top_area/top_distance, volume
    let areas_and_dists = [
        element_width * element_height,
        element_width * element_width,
        element_height * element_width / element_width,
        element_width * element_width / element_height,
        element_width * element_width * element_height,
    ];
    println!("model update begin");

    //input data
    //([f64;3],f64,f64,f64,[f64;2])
    //[kx,ky,kz],density,sp_heat,conv_coeff,[init_temp,tenv]



    println!("read temp vs G, viscosity, tandelta table");

    let mut viscometry_file = File::open("input/viscosity_data3.csv").expect("cant find the file");
    //let mut viscometry_file = File::open("/mnt/c/rustfiles/viscosity_data3.csv").expect("cant find the file with viscometry data");

    let mut bufreader = BufReader::with_capacity(10000, viscometry_file);
    let mut temp_vis = Vec::with_capacity(10000);
   let mut G_vis = Vec::with_capacity(10000);
   let mut n_vis = Vec::with_capacity(10000);
   let mut tandel_vis = Vec::with_capacity(10000);
    let n_units = 1.0; //MPa
    let G_units = 1.0; //MPa
    for i in bufreader.lines(){
        let dataln = i.expect("cant read data line").clone();
        let mut data = dataln.split(",");
        temp_vis.push(f64::from_str((data.next().expect("cant get next"))).expect("cant convert to f64"));
        G_vis.push(f64::from_str((data.next().expect("cant get next"))).expect("cant convert to f64") * G_units);
        n_vis.push(f64::from_str((data.next().expect("cant get next"))).expect("cant convert to f64") * n_units);
        tandel_vis.push(f64::from_str((data.next().expect("cant get next"))).expect("cant convert to f64"));
    }

    println!("viscosity data G {:?}",G_vis);
    println!("viscosity data n {:?}",n_vis);
    let tempstep = temp_vis[1]-temp_vis[0];
    let G_interp = Interpolator::new(temp_vis.clone(),G_vis,tempstep, maxthreads);
    let n_interp = Interpolator::new(temp_vis.clone(), n_vis, tempstep, maxthreads);
    let tandel_interp = Interpolator::new(temp_vis, tandel_vis, tempstep, maxthreads);




    let datastorer = model_iso_td_con_visco::Model::update_model(
        &mut mu,
        &mut mdl,
        &mut bw,
        areas_and_dists,
        0.5,
        &conductivity_interpolation_table,
        maxthreads*2,
        input_data,
        &pool,
        [&G_interp, &n_interp, &tandel_interp]
    );
    println!("simulation ended. start writing file");

    let mut strainout = File::create("input/strain_output.csv").expect("cant create strain output file");
    //let mut strainout = File::create("/mnt/c/rustfiles/strain_output.csv").expect("cant create strain output file");
    let mut stroutbuf = BufWriter::with_capacity(100000,strainout);
    for i in 0..mdl.strain_time_dep.len(){
        writeln!(stroutbuf, "{}",mdl.strain_time_dep[i]);
    }



    mdl.stm.calc_strain_for_all_elements();
    //println!("loads {:?}", mdl.stm.loads);
    mdl.stm.write("output/model_nodes.csv", "output/model_elements.csv", "output/debug_test.csv", "output/model_loads.csv", "output/model_bcs.csv", "output/model_strains.csv");
    println!("Hello, world!");


}
