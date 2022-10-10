/*let mut str_mdl = StructuralModel::new(activation_times_clone.len()*6);
let E = 2e9;
let nu = 0.3;
let G = E/(2.0*(1.0+nu));
let matprops = [E,E,E,G,G,G,nu,nu,nu];
for i in 0..activation_times_clone.len(){
    let newcell = muclone.get_next_cell_info();
    str_mdl.addCell(&mdl.nodelist,newcell.1,newcell.5,matprops);
}

str_mdl.add_bcs();
str_mdl.add_loads(&mdl.temp_at_contact);

//write nodes
let mut structure_file_nodes = File::create("C:\\rustfiles\\structure_data_nodes.csv").expect("cannot create the file");
let mut sf_nodes_buffer = BufWriter::with_capacity(100000,structure_file_nodes);
for i in str_mdl.nodes{
    write!(sf_nodes_buffer,"{},{},{}\n",i[0],i[1],i[2]);
}


//write elements
let mut structure_file_elems = File::create("c:\\rustfiles\\structure_data_elements.csv").expect("can't create file");
let mut sf_elems = BufWriter::with_capacity(100000,structure_file_elems);
for i in str_mdl.elements{
    write!(sf_elems,"{},{},{},{},{},{},{},{}\n",i[0],i[1],i[2],i[3],i[4],i[5],i[6],i[7]);
}

//write bcs
let mut structure_files_bcs = File::create("c:\\rustfiles\\structure_data_bcs.csv").expect("cant create file");
let mut sf_bcs = BufWriter::with_capacity(100000,structure_files_bcs);
for i in str_mdl.bcs{

    write!(sf_bcs,"{},{},{}\n",i[0] as usize,i[1] as usize,i[2] as usize);
}
//write loads
let mut structure_file_loads = File::create("c:\\rustfiles\\structure_data_loads.csv").expect("can't create file");
let mut sf_lds = BufWriter::with_capacity(100000,structure_file_loads);
for i in str_mdl.loads{
    write!(sf_lds,"{},{},{}\n",i[0],i[1],i[2]);
}*/


/*let mut model_file = File::create("c:\\rustfiles\\cells_and_nodes.csv").expect("can't create the file");
    let mut model_file_buffer = BufWriter::with_capacity(10000,model_file);
    for i in mdl.celllist{
    write!(model_file_buffer,"{}, {}, {}, {}, {}, {}, {} \n",i.index,mdl.nodelist[i.nodes[0]].pt.x, mdl.nodelist[i.nodes[7]].pt.x, mdl.nodelist[i.nodes[0]].pt.y, mdl.nodelist[i.nodes[7]].pt.y,
     mdl.nodelist[i.nodes[0]].pt.z, mdl.nodelist[i.nodes[7]].pt.z );
    }*/

//println!("intersects {}",is_intersect.0);
/*let mut temp_at_contact_file = File::create("C:\\rustfiles\\temp_at_contact.csv").expect("cant create the specified file");
let mut tacf_buffer = BufWriter::with_capacity(100000,temp_at_contact_file);

    for i in &mdl.temp_at_contact {
        write!(tacf_buffer, "{},{},{},{},{},{}\n", i[0], i[1], i[2], i[3], i[4], i[5]).unwrap();
    }

let mut temp_diff_at_end = File::create("c:\\rustfiles\\temp_diff_at_end.csv").expect("cant create file");
let mut tdae_buffer = BufWriter::with_capacity(10000,temp_diff_at_end);
for i in 0..mdl.temp_at_contact.len(){
    write!(tdae_buffer,"{},{},{},{},{},{}\n",mdl.temp_at_contact[i][0]-mdl.templist[i],mdl.temp_at_contact[i][1]-mdl.templist[i]
    ,mdl.temp_at_contact[i][2]-mdl.templist[i],mdl.temp_at_contact[i][3]-mdl.templist[i],mdl.temp_at_contact[i][4]-mdl.templist[i]
    ,mdl.temp_at_contact[i][5]-mdl.templist[i]).unwrap();
}

let mut cell_map_index_file = File::create("c:\\rustfiles\\cell_map_indx.csv").expect("can't create file");
let mut cmi_buffer = BufWriter::with_capacity(10000, cell_map_index_file);

for i in mdl.cell_index_map{
    write!(cmi_buffer,"{},{}\n",i.0,i.1);
}*/

/*for i in 0..gcrclone.segment.len(){
     write!(codedumpbuffer,"{},{},{},{},{},{},{},{}\n",gcrclone.segment[i][0].x,gcrclone.segment[i][0].y,gcrclone.segment[i][0].z,
     gcrclone.segment[i][1].x,gcrclone.segment[i][1].y,gcrclone.segment[i][1].z,gcrclone.is_extrusion_on[i],
     gcrclone.speed[i]);
 }*/
// println!("before center calc");
/*  for i in 0..8
  {
     // println!("{}",node_no[i]);
      let x = self.nodelist[node_no[i]].pt.x;
      let y = self.nodelist[node_no[i]].pt.y;
      let z = self.nodelist[node_no[i]].pt.z;
      centerx = x + centerx;
      centery = y + centery;
      centerz = z + centerz;
  }
 // println!("after center calc");
  let center = Point {
      x: centerx / 8.0,
      y: centery / 8.0,
      z: centerz / 8.0,
  };*/

// println!("{:?}",self.cell_index_map);

// println!("after center point");
//println!("previous_active_layer {}, active_layer_first_element {}", self.previous_active_layer, self.active_layer_first_element);

// let mut temps =  Vec::new();//    vec![ 15.0,	17.50,	20.0,	22.5000000000000	25	27.5000000000000	30	32.5000000000000	35	37.5000000000000	40	42.5000000000000	45	47.5000000000000	50	52.5000000000000	55	57.5000000000000	60	62.5000000000000	65	67.5000000000000	70	72.5000000000000	75	77.5000000000000	80	82.5000000000000	85	87.5000000000000	90	92.5000000000000	95	97.5000000000000	100	102.500000000000	105	107.500000000000	110	112.500000000000	115	117.500000000000	120	122.500000000000	125	127.500000000000	130	132.500000000000	135	137.500000000000	140	142.500000000000	145	147.500000000000	150	152.500000000000	155	157.500000000000	160	162.500000000000	165	167.500000000000	170	172.500000000000	175	177.500000000000	180	182.500000000000	185	187.500000000000	190	192.500000000000	195	197.500000000000	200	202.500000000000	205	207.500000000000	210	212.500000000000	215	217.500000000000	220	222.500000000000,	225	227.500000000000,	230.0,	232.500000000000,	235.0];
// let mut sp_ht_caps = Vec::new();//vec![ 1929.0,1929.0, 1929.0, 1929.0, 1929.0, 1929.0, 1929.0, 1929.0, 1929.0, 1929.0, 3000.0, 3000.0, 1929.0, 4500.0, 4500.0, 1929.0, 1929.0, 1929.0, 1929.0, 1929.0, 1929.0, 1929.0,1929.0];
/*let spheat_data_file = File::open("c:\\rustfiles\\sp_heat_cap_data.csv").expect("can't find file");
let spht_buf = BufReader::with_capacity(10000,spheat_data_file);
//let mut data = String::new();
for line in spht_buf.lines(){
    let line2 = line.expect("cant read the line").clone();
    let mut dat = line2.split(",");
    let sphtcap = dat.next().expect("row does not have first column");
    let temp = dat.next().expect("row does not have second column");
    temps.push(temp.parse::<f32>().expect("can't parse temp to f32"));
    sp_ht_caps.push(sphtcap.parse::<f32>().expect("can't parse sp heat cap to f32"));
}*/

// let xstep = 1.0;
/* let temps = vec![0.0,50.0,100.0,150.0,200.0,250.0];
let sp_ht_caps = vec![780.0,1040.0,1490.0,1710.0,1865.0,2020.0];
 let conductivity = vec![0.23, 0.25, 0.28, 0.29, 0.31, 0.33];
let xstep = 50.0;*/

/*  let activation_times = m.generate_activation_times(gcr.segment,gcr.is_extrusion_on,gcr.speed,beadwidth, beadheight);
   println!("activation times calculated");
   let activation_times_clone = activation_times.clone();
   /*for i in 0..activation_times.len(){
       if (activation_times[i] - 428.0).abs() < 1{
           println!("{}",i);
       }
   }*/
   */


/*
fn main() {
    let maxthreads= num_cpus::get_physical();
    let pool = rayon::ThreadPoolBuilder::new().num_threads(maxthreads).build().unwrap();
    let gcr = GCodeReader::new("c:\\rustfiles\\fea_bowl_part_half.gcode"); //cube.gcode
    let filename = "c:\\rustfiles\\gcodedump_med.csv";
    gcr.write_abaqus_input(filename);
    println!("file read");


    let beadwidth:f32 = 0.8* 1e-3;
    let beadheight = 0.4 * 1e-3;
    let xdiv = ((gcr.xmax-gcr.xmin)/beadwidth + 2.0).round() as usize ;
    let ydiv = ((gcr.ymax-gcr.ymin)/beadwidth  + 2.0).round() as usize;
    let zdiv =((gcr.zmax - gcr.zmin)/beadheight + 2.0).round() as usize;
    let m = ModelGenerator::new(gcr.xmin-beadwidth/2.0,gcr.xmax+beadwidth/2.0,gcr.ymin-beadwidth/2.0,gcr.ymax+beadwidth/2.0,gcr.zmin-beadheight,gcr.zmax,xdiv ,ydiv,zdiv,200.0);
    let ndsize = m.nodelist.len();
    println!("mesh generated");

    let activation_times = m.generate_activation_times_all_layers(gcr.segment,gcr.is_extrusion_on,gcr.speed,
    beadwidth,beadheight, &pool, maxthreads);
    let mut file1 = File::create("c:\\rustfiles\\activation_times_all_layers.csv").expect("cant create the file");
    let mut file1buf = BufWriter::with_capacity(10000,file1);
    for i in activation_times.clone(){
        write!(file1buf, "{},{},{}\n",i.0,i.1,i.3);
    }


    let mut activation_times_file = File::create("c:\\rustfiles\\activation_times_with_cellno.csv").expect("cant open file");
    let mut activation_times_file_buf = BufWriter::with_capacity(10000,activation_times_file);
    //specific heat capacity table
    let sp_heat_cap_interpolation_table = Interpolator::read_data_from_file("c:\\rustfiles\\sp_heat_cap_data.csv",1.0);
    //let conductivity_interpolation_table = Interpolator::new(temps.clone(), conductivity.clone(),xstep);

    println!("specific heat capacity table read");
    for i in &activation_times {
        let mut center_point = Point{x:0.0,y:0.0,z:0.0};
        let pts = m.celllist[i.1.clone()];
        for ii in &pts {

            center_point.x = center_point.x + m.nodelist[*ii].pt.x/8.0;
            center_point.y = center_point.y + m.nodelist[*ii].pt.y/8.0;
            center_point.z = center_point.z + m.nodelist[*ii].pt.z/8.0;
        }
        write!(activation_times_file_buf, "{},{},{},{},{},{},{} \n",i.0.clone(),i.1.clone(),
               center_point.x,center_point.y,center_point.z, i.2.clone()[0],i.2.clone()[1] );
    }
     //let activation_times:Vec<(f32,usize)> = (0..m.celllist.len().clone()).map(|i|(i as f32,i as usize) ).collect();
    let mut mu = ModelUpdater::new(activation_times,m);
    let mut muclone = mu.clone();
    let mut mdl = model::Model::new(mu.model_generator.nodelist.clone(),mu.activation_times.len(),sp_heat_cap_interpolation_table,
    beadwidth, beadheight,20);// let mut file = File::create("c:\\rustfiles\\activation_times.csv").expect("can't create file.");
    //let mut file = File::create("/mnt/c/rustfiles/final_temps.csv").expect("can't create file.");
    let mut file = File::create("c:\\rustfiles\\temp_with_time.csv").expect("can't create file.");
    let mut bw = BufWriter::with_capacity(10000,file);

    //side area, top area, side_area/side_distance, top_area/top_distance, volume
    let areas_and_dists = [beadwidth*beadheight, beadwidth*beadwidth, beadheight*beadwidth/beadwidth, beadwidth*beadwidth/beadheight, beadwidth*beadwidth*beadheight];
    println!("model update begin");


    //input data
    //([f32;3],f32,f32,f32,[f32;2])
    //[kx,ky,kz],density,sp_heat,conv_coeff,[init_temp,tenv]
    let input_data = ([0.25,0.25,0.25],1100.0,1929.0,100.0,[200.0,25.0]);

    let datastorer = model::Model::update_model( &mut mu, &mut mdl, &mut bw,areas_and_dists,1e-2, maxthreads,input_data, &pool);
    println!("simulation ended. start writing file");
  /*  for i in 0..datastorer.0.len() {
        write!(bw, "{},{},{:?}\n",datastorer.0[i],datastorer.2[i],datastorer.1[i]);
    }*/

    println!("Hello, world!");
    //println!("The system has {} physical cpus and {} logical cpus.",num_cpus::get_physical(), num_cpus::get());
}

 */

/*if has_neighbors[0]{
                self.Smat.insert(i,neighbors[i][0],0, 0.5/dx);
            }
            if has_neighbors[1]{
                self.Smat.insert(i,neighbors[i][1],1, 0.5/dx);
            }
            if has_neighbors[2]{
                self.Smat.insert(i,neighbors[i][2],2, 0.5/dy);
            }
            if has_neighbors[3]{
                self.Smat.insert(i,neighbors[i][3],3, 0.5/dy);
            }
            if has_neighbors[4]{
                self.Smat.insert(i,neighbors[i][4],4, 0.5/dz);
            }
            if has_neighbors[5]{
                self.Smat.insert(i,neighbors[i][5],5, 0.5/dz);
            }
 */