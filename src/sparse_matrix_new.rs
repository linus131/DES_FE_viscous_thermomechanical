use rayon::ThreadPool;
use std::borrow::{Borrow, BorrowMut};
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::str::FromStr;
use std::time::Instant;
use std::collections::HashMap;
use rustc_hash::{FxHashMap,FxHasher};

const ERR:f64 = 1e-6;
const MAXITER:usize = 4000;
const CONVERR:f64 = 1e-9;

pub (crate) struct SparseMatrix{
    pub (crate) data: Vec<Vec<f64>>,
    pub (crate) row_col: Vec<Vec<usize>>,
    pub (crate) current_size:usize,
    rowmaps:Vec<FxHashMap<usize,usize>>,
    //diagonal_scaler: Vec<f64>
}

impl SparseMatrix{
    pub fn with_capacity(capacity: usize)->SparseMatrix{
        let mut data = Vec::with_capacity(capacity);
        let mut row_col = Vec::with_capacity(capacity);
        let mut rowmaps = Vec::with_capacity(capacity);
        // let mut diagonal_scaler = Vec::with_capacity(capacity);
        let mut current_size = 0;
        for i in 0..capacity{
            data.push(Vec::with_capacity(40));
            row_col.push(Vec::with_capacity(40));
            rowmaps.push(FxHashMap::default());
            // diagonal_scaler.push(1.0);
        }
        SparseMatrix{
            data,
            row_col,
            current_size,
            rowmaps,
            // diagonal_scaler,
        }
    }

    pub fn insert_rows_par(Kbc: &mut SparseMatrix, dof_to_bc_dof: &[usize], dof_from_bc_dof:&[usize],
                           bcs: &[bool], dof_to_elem_dof: &[Vec<usize>], dof_to_elem:&[Vec<usize>],
                           elements: &[Vec<usize>], Kloc: &Vec<[[f64;24];24]>,
                           maxthreads:usize, pool: &ThreadPool) {
        let data = &mut Kbc.data;
        let row_col = &mut Kbc.row_col;
        let rowmaps = &mut Kbc.rowmaps;
        Kbc.current_size = dof_from_bc_dof.len();
        // let counter:Vec<usize> = (0..dof_from_elem_dof.len()).map(|i| i).collect();
        SparseMatrix::insert_rows_par_helper(
            data, row_col, rowmaps, dof_to_bc_dof, dof_from_bc_dof, bcs, dof_to_elem_dof, dof_to_elem,
            elements, Kloc, 0, 1, maxthreads, pool);
    }

    fn insert_rows_par_helper(data:&mut[Vec<f64>], row_col: &mut[Vec<usize>], rowmaps: &mut[FxHashMap<usize,usize>],
                              dof_to_bc_dof: &[usize], dof_from_bc_dof:&[usize], bcs: &[bool],
                              dof_to_element_dof: &[Vec<usize>], dof_to_element:&[Vec<usize>],
                              elements: &[Vec<usize>], Kloc: &Vec<[[f64;24];24]>, startpos:usize,
                              numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        if numthreads<maxthreads{
            let splitpos = dof_from_bc_dof.len()/2;
            let(data1, data2) = data.split_at_mut(splitpos);
            let(rc1, rc2) = row_col.split_at_mut(splitpos);
            let(rm1, rm2) = rowmaps.split_at_mut(splitpos);
            let(dfbd1, dfbd2) = dof_from_bc_dof.split_at(splitpos);
            //let(c1, c2) = counter.split_at(splitpos);
            //let(dtedof1, dtedof2) = dof_to_elem_dof.split_at(splitpos);

            pool.install(||rayon::join(
                ||SparseMatrix::insert_rows_par_helper(data1, rc1, rm1, dof_to_bc_dof,dfbd1,bcs,dof_to_element_dof, dof_to_element,
                                                       elements, Kloc, startpos,numthreads*2, maxthreads, pool),
                ||SparseMatrix::insert_rows_par_helper(data2, rc2, rm2, dof_to_bc_dof,dfbd2,bcs,dof_to_element_dof, dof_to_element,
                                                       elements, Kloc, startpos+dfbd1.len(), numthreads*2, maxthreads, pool),
            )
            );

        }
        else{
            for i in 0..dof_from_bc_dof.len(){
                let dofcounter = dof_from_bc_dof[i];
                for j in 0..dof_to_element_dof[dofcounter].len(){
                    let row_Kloc = dof_to_element_dof[dofcounter][j];
                    let elemno = dof_to_element[dofcounter][j];
                    let cols = &elements[elemno];
                    let Klocrow = &Kloc[elemno][row_Kloc];
                    for k in 0..cols.len(){
                        for l in 0..3{
                            let cdof = cols[k]*3+l;
                            if !bcs[cdof]{
                                // code for inserting into the stiffness matrix
                                let row = i;
                                let col = dof_to_bc_dof[cdof];
                                let value = Klocrow[k*3+l];

                                if rowmaps[row].contains_key(&col){
                                    let colloc = rowmaps[row].get(&col).expect("can't find key");
                                    data[row][*colloc] += value;
                                }
                                else{
                                    let newcol = data[row].len();
                                    rowmaps[row].insert(col,newcol);
                                    data[row].push(value);
                                    row_col[row].push(col);
                                }

                            }
                        }
                    }

                }
            }


        }


    }

    pub fn insert(&mut self, row:usize, col:usize, value:f64){
        if row>=self.current_size {self.current_size = row+1}

        //err if row > max size allocated
        if self.rowmaps[row].contains_key(&col){
            let colloc = self.rowmaps[row].get(&col).expect("can't find key");
            self.data[row][*colloc] += value;
        }
        else{
            let newcol = self.data[row].len();
            self.rowmaps[row].insert(col,newcol);
            self.data[row].push(value);
            self.row_col[row].push(col);
        }
    }




    pub fn reset(&mut self){
        for i in 0..self.current_size {
            self.rowmaps[i].clear();
            self.data[i].clear();
            self.row_col[i].clear();
            self.current_size = 0;
            // self.diagonal_scaler.clear();
        }
    }

    pub fn issymmetric(&self)->bool{
        for i in 0..self.current_size{
            for j in 0..self.data[i].len(){
                if (self.get_element(i,j)-self.get_element(j,i)).abs() > ERR{
                    println!("{},{},{},{}", i,j, self.get_element(i,j), self.get_element(j,i));
                    return false;
                }
            }
        }
        return true
    }


    pub fn multiply(&self, b:&Vec<f64>, c:&mut Vec<f64>){
        for i in 0..self.current_size{
            for j in 0..self.data[i].len(){
                c[i] = c[i]+self.data[i][j] * b[self.row_col[i][j]];
            }
        }
    }

    pub fn multiply2(&self, b:&[f64], c:&mut[f64]){
        for i in 0..self.current_size{
            for j in 0..self.data[i].len(){
                c[i] = c[i]+self.data[i][j] * b[self.row_col[i][j]];
            }
        }
    }

    pub fn multiply_par(&self, b:&Vec<f64>, c:&mut Vec<f64>,  numthreads:usize, maxthreads:usize, pool: &ThreadPool ){
        //println!("nel {:?}",nel);
        SparseMatrix::mult_par_inner(self.data.borrow(), self.row_col.borrow(),  b.borrow(), c.borrow_mut(),numthreads, maxthreads, pool);
    }

    pub fn multiply_par2(&self, b:&[f64], c:&mut [f64],  numthreads:usize, maxthreads:usize, pool: &ThreadPool ){
        //println!("nel {:?}",nel);
        let data_borrow_tmp:&[Vec<f64>] = (self.data.borrow());
        let (d1,d2) = data_borrow_tmp.split_at(self.current_size);
        SparseMatrix::mult_par_inner(d1, self.row_col.borrow(), b, c,numthreads, maxthreads,pool);
    }

    fn mult_par_inner(data:&[Vec<f64>], row_col:&[Vec<usize>], b:&[f64],c:&mut [f64], numthreads: usize, maxthreads:usize, pool: &ThreadPool){
        if numthreads<maxthreads{
            let splitpos = c.len()/2;
            let (d1,d2) = data.split_at(splitpos);
            let (rc1,rc2) = row_col.split_at(splitpos);
            let (c1,c2) = c.split_at_mut(splitpos);

            pool.install(||rayon::join(|| SparseMatrix::mult_par_inner(d1,rc1,b,c1,numthreads*2,maxthreads,pool),
                                       || SparseMatrix::mult_par_inner(d2,rc2,b,c2,numthreads*2,maxthreads,pool)
            ));
        }
        else{
            for i in 0..c.len(){
                for j in 0..data[i].len(){
                    c[i] = c[i] + data[i][j] * b[row_col[i][j]];
                }
            }
        }
    }


    pub fn get_element(&self, row:usize, col:usize)->f64{
        // row < numel and col<numel
        if self.rowmaps[row].contains_key(&col){
            let loc = self.rowmaps[row].get(&col).expect("key not found");
            return self.data[row][*loc];
        }
        return 0.0
    }

    ///not working properly
    pub fn printmatrix(&self)
    {
        for i in 0..self.current_size{
            print!("[");
            for j in  0..self.current_size-1 {
                print!("{}, ", self.get_element(i,j));
            }
            print!("{} ", self.get_element(i,self.current_size-1));
            println!("]");
        }
    }

    pub(crate) fn KlocplusBtCBdetmatjac(B:&[[f64;24];6], C:&[[f64;6];6], detmatjac:f64, tmp: &mut[[f64;6];24], tmp2:&mut[[f64;24];24]){
        //Bt C
        for i in 0..B.len(){
            for j in 0..C[0].len(){
                for k in 0..C.len(){
                    tmp[i][j] += B[i][k]*C[k][j];
                }
            }
        }

        // BtCB
        for i in 0..tmp.len(){
            for j in 0..B[0].len(){
                for k in 0..B.len(){
                    tmp2[i][j] += tmp[k][i] * B[k][j];
                }
            }
        }
        //BtCB matjac
        for i in 0..24{
            for j in 0..24{
                tmp2[i][j] = tmp2[i][j] * detmatjac;
            }
        }

    }

    pub(crate) fn dense_mat_mult_easy2(a: &[[f64;8];3], b: &[[f64;3];8])->[[f64;3];3]
    {
        let mut c = [[0.0_f64;3];3];

        for i in 0..a.len() {
            for j in 0..b[0].len() {
                for k in 0..b.len() {
                    c[i][j] = c[i][j] + a[i][k] * b[k][j];
                }
            }
        }
        return c

    }

    pub(crate) fn dense_mat_mult_easy(a: &[[f64;6];24], b: &[[f64;24];6])->[[f64;24];24]
    {
        let mut c = [[0.0_f64;24];24];

        for i in 0..a.len() {
            for j in 0..b[0].len() {
                for k in 0..b.len() {
                    c[i][j] = c[i][j] + a[i][k] * b[k][j];
                }
            }
        }
        return c

    }

    /// multiply dense matrix with vector
    pub(crate) fn dense_mat_mult_vec_easy(a: &Vec<Vec<f64>>, b: &Vec<f64>)->Vec<f64>
    {
        let mut vout = vec![0.0_f64;a.len()];
        for i in 0.. a.len(){
            for j in 0..b.len(){
                vout[i] = vout[i]+a[i][j]*b[j];
            }
        }
        return vout;
    }

    /// multiply transpose matrix with matrix
    pub(crate) fn dense_mat_mult_easy_transpose(a: &[[f64;24];6], b: &Vec<[f64;6]>)->[[f64;6];24]
    {
        let mut c = [[0.0;6];24];
        for i in 0..a[0].len() {
            for j in 0..b.len(){
                for k in 0..b[0].len(){
                    c[i][j] = c[i][j] + a[k][i] * b[k][j];
                }
            }
        }
        return c;
    }

    /// add two dense matrices
    pub(crate) fn add_dense_matrix(a:&mut [[f64;24];24], b:&[[f64;24];24]){
        for i in 0..a.len(){
            for j in 0..a[0].len(){
                a[i][j] = a[i][j]+b[i][j];
            }
        }
    }

    pub(crate) fn multiply_dense_by_constant(a:&mut [[f64;24];24], constant:f64){
        for i in 0..a.len(){
            for j in 0..a[0].len(){
                a[i][j] = a[i][j] * constant;
            }
        }
    }



    pub(crate) fn from_file(filename: &str, nrows: usize, ncolsmax: usize) -> SparseMatrix{
        let mut file = File::open(filename).expect("cant open file");
        let mut buffile = BufReader::with_capacity(10000,file);
        let mut out = SparseMatrix::with_capacity(nrows);
        for i in 0..nrows{
            let mut tmpstr = String::new();
            buffile.read_line(&mut tmpstr).expect("can't read the buffer");
            let tmpstrlen = tmpstr.len();
            let mut split = tmpstr[0..tmpstrlen-1].split(",");
            for j in 0..nrows{
                let tmp = split.next().expect("no next");
                // println!("tmp is {}", tmp);
                let val = f64::from_str(tmp).expect("can't parse to float");
                out.insert(i,j,val);
            }
        }
        return out
    }
    /// reads sparse matrix from file gen_rand_sym_mat_maxcol
    pub(crate) fn from_file_sparse(filename: &str, nrows: usize, ncolsmax: usize) -> SparseMatrix{
        let mut file = File::open(filename).expect("cant open file");
        let mut buffile = BufReader::with_capacity(10000,file);
        let mut out = SparseMatrix::with_capacity(nrows);
        let mut tmpstr = String::new();
        buffile.read_line(&mut tmpstr).expect("can't read buffer");
        tmpstr.pop();
        //println!("tmpstr {}", tmpstr);
        let numelems = usize::from_str(tmpstr.as_str()).expect("can't parse to usize");

        for i in 0..numelems{
            let mut tmpstr = String::new();
            buffile.read_line(&mut tmpstr).expect("can't read the buffer");
            let tmpstrlen = tmpstr.len();
            let mut split = tmpstr[0..tmpstrlen-1].split(",");
            let r = usize::from_str(split.next().expect("no next")).expect("can't parse to usize");
            let c = usize::from_str(split.next().expect("no next")).expect("can't parse to usize");
            let nbr = usize::from_str(split.next().expect("no next")).expect("can't parse to usize");
            let val = f64::from_str(split.next().expect("no next")).expect("can't parse to f64");
            out.insert(r,c,val);
        }
        return out
    }

    /// needs file with single line of values
    pub(crate) fn vec_from_file(filename: &str, nrows: usize) ->Vec<f64> {
        let mut out = Vec::with_capacity(nrows);
        let mut file = File::open(filename).expect("cant open file");
        let mut buffile = BufReader::with_capacity(10000,file);
        let mut tmpstr = String::new();
        buffile.read_line(&mut tmpstr).expect("can't read the buffer");
        let tmpstrlen = tmpstr.len();
        let mut split = tmpstr[0..tmpstrlen-1].split(",");
        for i in 0..nrows{
            out.push(f64::from_str(split.next().expect("can't get next")).expect("can't parse to f64"));
        }
        return out
    }

    fn b_minus_ax(&self, b: &[f64], x: &[f64], r: &mut[f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){

        //self.multiply2(x,r);
        //  for i in 0..r.len(){
        //      r[i] = 0.0;
        //  }
        self.multiply_par2(x, r, numthreads, maxthreads, pool);
        SparseMatrix::aminb_par(b,r, 1, maxthreads, pool);
        // for i in 0..r.len(){
        //     r[i] = b[i]- r[i];
        //  }
    }

    fn aminb_par(a:&[f64], b:&mut [f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        if numthreads< maxthreads{
            let splitpos = b.len()/2;
            let(a1,a2) = a.split_at(splitpos);
            let(b1,b2) = b.split_at_mut(splitpos);
            // let(c1, c2) =c.split_at_mut(splitpos);
            pool.install(|| rayon::join(
                ||SparseMatrix::aminb_par(a1, b1,  numthreads*2, maxthreads, pool),
                ||SparseMatrix::aminb_par(a2, b2,  numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..b.len(){
                b[i] = a[i] - b[i];
            }
        }
    }

    pub fn diagonal_scaling(&self, scaleout: &mut [f64]){
        for i in 0..self.current_size{
            let colpos = self.rowmaps[i].get(&i).expect("can't find key");
            let val = 1.0/self.data[i][*colpos];
            scaleout[i] = val;
        }
    }

    pub fn diagonal_scaling_par(scaleout:&mut [f64], data: &[Vec<f64>], rowmaps:&[FxHashMap<usize,usize>], startadd: usize, numthreads:usize, maxthreads:usize, pool: &ThreadPool){
        if numthreads<maxthreads{
            let splitpos = scaleout.len()/2;
            let (sc0, sc1) = scaleout.split_at_mut(splitpos);
            let (data1, data2) = data.split_at(splitpos);
            let (rm1, rm2) = rowmaps.split_at(splitpos);
            pool.install(|| rayon::join(
                ||SparseMatrix::diagonal_scaling_par(sc0, data1, rm1,  startadd, numthreads*2, maxthreads, pool),
                ||SparseMatrix::diagonal_scaling_par(sc1, data2, rm2,  startadd+splitpos,numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..scaleout.len(){
                //println!("startadd, {}", startadd);
                let colpos = rowmaps[i].get(&(startadd+i)).expect("cant find the key");
                let val = 1.0/data[i] [*colpos];
                scaleout[i] = val;
            }
        }
    }

    fn axb_par(a:&[f64], b:&[f64], c:&mut[f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        if numthreads< maxthreads{
            let splitpos = a.len()/2;
            let(a1,a2) = a.split_at(splitpos);
            let(b1,b2) = b.split_at(splitpos);
            let(c1, c2) =c.split_at_mut(splitpos);
            pool.install(|| rayon::join(
                ||SparseMatrix::axb_par(a1, b1, c1, numthreads*2, maxthreads, pool),
                ||SparseMatrix::axb_par(a2, b2, c2, numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..c.len(){
                c[i] = a[i] * b[i];
            }
        }
    }
    fn assign_par(a:&[f64], b:&mut[f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){

        if numthreads< maxthreads{
            let splitpos = a.len()/2;
            let(a1,a2) = a.split_at(splitpos);
            let(b1,b2) = b.split_at_mut(splitpos);
            // let(c1, c2) =c.split_at_mut(splitpos);
            pool.install(|| rayon::join(
                ||SparseMatrix::assign_par(a1, b1,  numthreads*2, maxthreads, pool),
                ||SparseMatrix::assign_par(a2, b2,  numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..b.len(){
                b[i] = a[i];
            }
        }


    }

    fn axb_accumulate_par(a:&[f64],b:&[f64], store:&mut[f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool)->f64{
        //let mut store = vec![0.0; maxthreads*2];
        for kk in 0..store.len(){
            store[kk] = 0.0_f64;
        }
        let mut c = 0.0;
        SparseMatrix::axb_accumulate_par_help(a,b,store,numthreads,maxthreads,pool);
        for i in 0..store.len(){
            c += store[i];
        }
        return c
    }

    fn axb_accumulate_par_help(a:&[f64], b:&[f64], c:&mut [f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        if numthreads< maxthreads{
            let splitpos = a.len()/2;
            let(a1,a2) = a.split_at(splitpos);
            let(b1,b2) = b.split_at(splitpos);
            let splitc = c.len()/2;
            let(c1, c2) = c.split_at_mut(splitc);
            pool.install(|| rayon::join(
                ||SparseMatrix::axb_accumulate_par_help(a1, b1, c1, numthreads*2, maxthreads, pool),
                ||SparseMatrix::axb_accumulate_par_help(a2, b2, c2, numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..a.len(){
                c[0] += a[i] * b[i];
            }
        }
    }

    pub(crate) fn assign_val_par(a:&mut[f64], b:f64, numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        if numthreads< maxthreads{
            let splitpos = a.len()/2;
            let(a1,a2) = a.split_at_mut(splitpos);
            pool.install(|| rayon::join(
                ||SparseMatrix::assign_val_par(a1, b,numthreads*2, maxthreads, pool),
                ||SparseMatrix::assign_val_par(a2, b, numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..a.len(){
                a[i] = b;
            }
        }
    }

    fn a_acc_b_times_c(a:&mut[f64], b:f64, c:&[f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        if numthreads< maxthreads{
            let splitpos = a.len()/2;
            let(a1,a2) = a.split_at_mut(splitpos);
            //let(b1,b2) = b.split_at(splitpos);
            let(c1,c2) = c.split_at(splitpos);
            pool.install(|| rayon::join(
                ||SparseMatrix::a_acc_b_times_c(a1, b, c1, numthreads*2, maxthreads, pool),
                ||SparseMatrix::a_acc_b_times_c(a2, b, c2, numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..a.len(){
                a[i] += b*c[i];
            }
        }
    }

    fn a_eq_b_plus_cxa(a:&mut[f64], b:&[f64], c:f64, numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        if numthreads< maxthreads{
            let splitpos = a.len()/2;
            let(a1,a2) = a.split_at_mut(splitpos);
            let(b1,b2) = b.split_at(splitpos);
            // let(d1,d2) = d.split_at(splitpos);
            pool.install(|| rayon::join(
                ||SparseMatrix::a_eq_b_plus_cxa(a1, b1, c, numthreads*2, maxthreads, pool),
                ||SparseMatrix::a_eq_b_plus_cxa(a2, b2, c, numthreads*2, maxthreads, pool)
            ));
        }
        else{
            for i in 0..a.len(){
                a[i] = b[i]+c*a[i];
            }
        }
    }

    // only works if the matrix is symmetric and determinant is nonzero
    pub fn cgsolve(&self, b:&[f64], x:&mut [f64], r:&mut[f64], p:&mut[f64], tmp: &mut[f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        //self.printmatrix();
        let err = CONVERR;
        let maxiter = MAXITER;
        //r = b - Ax
        self.b_minus_ax(b,x, r, numthreads, maxthreads,pool);
        //println!("b\n{:?}",b);
        //println!("x\n{:?}",x);
        //println!("r\n{:?}",r);
        //p = r
        for i in 0..r.len(){
            p[i] = r[i]
        }

        // loop start
        for j in 0..maxiter{
            //rtk * rk
            let mut rtkrk = 0.0;
            for i in 0..r.len(){
                rtkrk += r[i]*r[i];
            }

            //pkt * (A * pk)
            for i in 0..tmp.len(){
                tmp[i] = 0.0;
            }
            self.multiply_par2(p,tmp,numthreads,maxthreads,pool);
            //self.multiply2(p,tmp);
            let mut ptkApk = 0.0;
            for i in 0..p.len(){
                ptkApk += p[i]*tmp[i];
            }
            // ak = (rtk * rk) / (ptk* A* pk)
            let ak = rtkrk/ptkApk;
            //xk+1 = xk + ak * pk
            let mut sumri = 0.0;
            for i in 0..x.len(){
                x[i] = x[i] + ak * p[i];
                r[i] = r[i] - ak * tmp[i];
                sumri += r[i]*r[i];
            }
            if (sumri<err) {println!("converged in {} iterations", j);return}
            let bk = sumri/rtkrk;
            for i in 0..p.len(){
                p[i] = r[i] + bk * p[i];
            }
        }

        println!("did not converge in {} iterations", maxiter);
    }

    // preconditioning with diagonal matrix 1/sqrt(aii) Fletcher-Reeves
    pub fn cgsolve_diag_scl_precon(&self, b:&[f64], x:&mut [f64], r:&mut[f64], p:&mut[f64], tmp: &mut[f64], Mi: &mut Vec<f64>, z0: &mut Vec<f64>,  numthreads:usize, maxthreads:usize, pool:&ThreadPool) -> usize{
        let err = CONVERR;
        let maxiter = MAXITER;
        SparseMatrix::diagonal_scaling_par( Mi,self.data.borrow(), self.rowmaps.borrow(),0,1, maxthreads, pool);
        self.b_minus_ax(b,x, r, numthreads, maxthreads,pool);
        let mut store = vec![0.0;maxthreads*2];
        SparseMatrix::axb_par(Mi,r,z0, 1, maxthreads, pool);
        SparseMatrix::assign_par(z0,p, 1, maxthreads, pool);
        let mut sumri = 0.0;
        for j in 0..maxiter{
            let rtkzk = SparseMatrix::axb_accumulate_par(r, z0,  &mut store,1, maxthreads, pool);
            SparseMatrix::assign_val_par(tmp,0.0_f64, numthreads, maxthreads, pool);
            self.multiply_par2(p,tmp,numthreads,maxthreads,pool);
            let mut ptkApk = SparseMatrix::axb_accumulate_par(p,tmp,&mut store, 1, maxthreads, pool);
            let ak = rtkzk/ptkApk;
            sumri = 0.0;
            SparseMatrix::a_acc_b_times_c(x,ak,p,1,maxthreads,pool);
            SparseMatrix::a_acc_b_times_c(r,-ak,tmp,1,maxthreads,pool);
            let sumri = SparseMatrix::axb_accumulate_par(r,r,&mut store, 1, maxthreads, pool);
            if (sumri<err) {return j}
            SparseMatrix::axb_par(Mi,r,z0,1,maxthreads,pool);
            let mut rkp1tzkp1 = SparseMatrix::axb_accumulate_par(z0,r,&mut store, 1, maxthreads,pool);
            let bk = rkp1tzkp1/rtkzk;
            SparseMatrix::a_eq_b_plus_cxa(p,z0,bk, numthreads, maxthreads, pool);
        }

        println!("did not converge in {} iterations, err {}", maxiter, sumri);
        return MAXITER
    }

    // preconditioning with diagonal matrix 1/sqrt(aii) Polak_Ribiere
    pub fn cgsolve_diag_scl_precon2(&self, b:&[f64], x:&mut [f64], r:&mut[f64], p:&mut[f64], tmp: &mut[f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        //self.printmatrix();
        let err = CONVERR;
        let maxiter = MAXITER;
        let mut Mi = vec![0.0;b.len()];
        let mut z0 = vec![0.0;b.len()];
        let mut zkp1 = vec![0.0;b.len()];
        self.diagonal_scaling(&mut Mi);
        //println!("{:?}",Mi);
        //r = b - Ax
        self.b_minus_ax(b,x, r, numthreads, maxthreads,pool);
        for i in 0..Mi.len(){
            z0[i] = Mi[i] * r[i];
        }
        for i in 0..r.len(){
            p[i] = z0[i]
        }
        let mut sumri = 0.0;
        // loop start
        for j in 0..maxiter{
            //rtk * rk
            let mut rtkzk = 0.0;
            for i in 0..r.len(){
                rtkzk += r[i]*z0[i];
            }

            //pkt * (A * pk)
            for i in 0..tmp.len(){
                tmp[i] = 0.0;
            }
            self.multiply_par2(p,tmp,numthreads,maxthreads,pool);
            //self.multiply2(p,tmp);
            let mut ptkApk = 0.0;
            for i in 0..p.len(){
                ptkApk += p[i]*tmp[i];
            }
            // ak = (rtk * rk) / (ptk* A* pk)
            let ak = rtkzk/ptkApk;
            //xk+1 = xk + ak * pk
            sumri = 0.0;
            for i in 0..x.len(){
                x[i] = x[i] + ak * p[i];
                r[i] = r[i] - ak * tmp[i];
                sumri += r[i]*r[i];
            }
            if (sumri<err) {println!("converged in {} iterations", j);return}
            let mut rkp1tzkp1 = 0.0;
            for i in 0..zkp1.len(){
                zkp1[i] = Mi[i]*r[i];
                rkp1tzkp1 += (zkp1[i]-z0[i])*r[i];
            }

            let bk = rkp1tzkp1/rtkzk;
            for i in 0..p.len(){
                p[i] = zkp1[i] + bk * p[i];
            }
            z0.copy_from_slice(&zkp1);
        }

        println!("WARNING - did not converge in {} iterations, err {}", maxiter, sumri);
    }

    /// stabilized conjugate gradient with incomplete cholesky factorization preconditioner
    /// only works if determinant is nonzero (non singular matrix) and matrix is symmetric
    pub fn cgsolve_incomplete_cholesky_precon(&self, b:&[f64], x:&mut [f64], r:&mut[f64], p:&mut[f64], tmp: &mut[f64], numthreads:usize, maxthreads:usize, pool:&ThreadPool){
        todo!();
    }
}




fn main() {

    let mut mat = SparseMatrix::with_capacity(3);
    mat.insert(0,0,1.0);
    mat.insert(0,1,2.0);
    mat.insert(0,2,3.0);
    mat.insert(1,0,4.0);
    mat.insert(1,1,5.0);
    mat.insert(1,2,10.0);
    mat.insert(2,0,9.0);
    mat.insert(2,1,8.0);
    mat.insert(2,2,12.0);

    mat.printmatrix();

    let mut b = vec![-1.0,1.0,2.0];
    let mut c= vec![0.0;3];

    mat.multiply(&b, &mut c);
    println!("mm A*b -> c{:?}",c);
    let mut c= vec![0.0;3];

    let pool = rayon::ThreadPoolBuilder::new().num_threads(6).build().expect("can't build threads");
    //let nel:Vec<usize> = (0..b.len()).map(|i| i).collect();
    mat.multiply_par(&b,&mut c, 1,2, &pool);
    println!("mm A*b-> parmult{:?}",c);

    let mut A = SparseMatrix::with_capacity(3);
    A.insert(0,0,4.0);
    A.insert(0,1,1.0);
    A.insert(1,0,1.0);
    A.insert(1,1,3.0);

    let mut b = vec![1.0,2.0];



    let mut x = vec![0.0;2];
    let mut r = vec![0.0;2];
    let mut p = vec![0.0;2];
    let mut tmp = vec![0.0;2];
    A.cgsolve(b.borrow(),&mut x,&mut r,&mut p,&mut tmp,1,6,&pool);
    println!("cgsolve");

    println!("{:?}",x);

    let mut mul = vec![0.0;2];
    A.multiply_par(&x, &mut mul, 1,6,&pool);
    println!("{:?}",mul);

    let A2 = SparseMatrix::from_file("c:\\rustfiles\\testmatrix.csv",12,12);

    A2.printmatrix();

    let b = SparseMatrix::vec_from_file("c:\\rustfiles\\testvector.csv",12);
    println!("{:?}", b);

    let mut x = vec![0.0;12];
    let mut r = vec![0.0;12];
    let mut p = vec![0.0;12];
    let mut tmp = vec![0.0;12];

    A2.cgsolve(b.borrow(),x.borrow_mut(),r.borrow_mut(),p.borrow_mut(),tmp.borrow_mut(),1,6,&pool);

    println!("{:?}",x);

    println!("issymmetric A2, {}", A2.issymmetric());
    let tic = Instant::now();
    let A3 = SparseMatrix::from_file_sparse("c:\\rustfiles\\sparsetestmat.csv",1728,17);

    // A3.printmatrix();


    let b = SparseMatrix::vec_from_file("c:\\rustfiles\\testvector2.csv",1728);
    println!("sizeb, {} \n {:?}", b.len(),b);
    println!("time taken to read file {:?}", tic.elapsed());


    let mut x = vec![0.0;1728];
    let mut r = vec![0.0;1728];
    let mut p = vec![0.0;1728];
    let mut tmp = vec![0.0;1728];
    // let tic = Instant::now();
    //println!("issymmetric A3, {}", A3.issymmetric());
    //println!("time taken to check if symmetric {:?}",tic.elapsed());

    let tic = Instant::now();
    A3.cgsolve(b.borrow(),x.borrow_mut(),r.borrow_mut(),p.borrow_mut(),tmp.borrow_mut(),1,6,&pool);
    println!("time taken to solve, {:?}", tic.elapsed());

    println!("{:?}",x);



    println!("Hello, world!");
}
