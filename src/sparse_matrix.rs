use std::fmt;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use rayon::ThreadPool;

#[derive(Debug)]
pub (crate) struct SparseMatrix{
    nrows:usize,
    ncols:usize,
    hashmap: HashMap<usize,usize>,
    col_no: Vec<Vec<usize>>,
    col_no_true: Vec<Vec<usize>>,
    data: Vec<Vec<f64>>,
    nnz: usize,
}

impl SparseMatrix {
    /// create sparse matrix with 3x3 initial capacity
    fn new(size: [usize; 2]) -> SparseMatrix {
        let mut data = Vec::with_capacity(size[0]);
        let mut col_no = Vec::with_capacity(size[0]);
        let mut col_no_true = Vec::with_capacity(size[0]);
        for i in 0..size[0] {
            data.push(Vec::new());
            col_no.push(Vec::new());
            col_no_true.push(Vec::new());
        }
        return SparseMatrix {
            nrows: size[0],
            ncols: size[1],
            hashmap: HashMap::with_capacity(9),
            col_no: col_no,
            col_no_true: col_no_true,
            data: data,
            nnz: 0
        }
    }

    /// create sparse matrix with capacity
    pub(crate) fn with_capacity(size: [usize; 2], nnz: usize) -> SparseMatrix {
        let mut data = Vec::with_capacity(size[0]);
        let mut col_no = Vec::with_capacity(size[0]);
        let mut col_no_true = Vec::with_capacity(size[0]);
        for i in 0..size[0] {
            data.push(Vec::new());
            col_no.push(Vec::new());
            col_no_true.push(Vec::new());
        }
        return SparseMatrix {
            nrows: size[0],
            ncols: size[1],
            hashmap: HashMap::with_capacity(9),
            col_no: col_no,
            col_no_true: col_no_true,
            data: data,
            nnz: 0
        }
    }

    /// if row and col exists, replaces value, else inserts (replaces)
    fn insert(&mut self, row: usize, col: usize, data: f64) {
        let loc = self.nrows * row + col;
        if self.hashmap.contains_key(&loc) {
            let col_loc = self.hashmap.get(&loc).expect("cant find key");
            self.data[row][*col_loc] = data;
        } else {
            self.hashmap.insert(loc, self.col_no[row].len());
            let col_number = self.col_no[row].len();
            self.col_no[row].push(col_number);
            self.col_no_true[row].push(col);
            self.data[row].push(data);
            self.nnz = self.nnz + 1;
        }
    }

    /// if row and col exist, adds to existing value, else inserts (does not replace)
    pub(crate) fn insert_add(&mut self, row: usize, col: usize, data: f64) {
        let loc = self.nrows * row + col;
        if self.hashmap.contains_key(&loc) {
            let col_loc = self.hashmap.get(&loc).expect("cant find key");
            self.data[row][*col_loc] = self.data[row][*col_loc] + data;
        } else {

            //println!("add row {}, max row {}",row, self.nrows);
            self.hashmap.insert(loc, self.col_no[row].len());
            let col_number = self.col_no[row].len();
            self.col_no[row].push(col_number);
            self.col_no_true[row].push(col);
            self.data[row].push(data);
            self.nnz = self.nnz + 1;
        }
    }

    /// gets value for given rows and cols
    fn get(&self, row: usize, col: usize) -> f64 {
        let key = row * self.nrows + col;
        let data_out: f64;
        if self.hashmap.contains_key(&key) {
            let ncol = self.hashmap.get(&key).expect("cant find key");
            data_out = self.data[row][*ncol];
        } else {
            data_out = 0.0_f64;
        }
        return data_out;
    }

    /// multiplies with column vector (single column)
    pub(crate) fn multiply(&self, vec_in: &Vec<f64>, vec_out: &mut Vec<f64>) {
        /*for i in 0..vec_out.len(){
            vec_out[i] = 0_f64;
        }*/
        for i in 0..self.col_no.len() {
            vec_out[i] = 0_f64;
            for j in 0..self.col_no[i].len() {
                let col_no_val = self.col_no[i][j];
                let col_no_true_val = self.col_no_true[i][j];
                vec_out[i] = vec_out[i] + self.data[i][col_no_val] * vec_in[col_no_true_val];
            }
        }
    }

    pub(crate) fn multiply_par(&self, vec_in: &Vec<f64>, vec_out: &mut Vec<f64>, pool:&ThreadPool) {
        /*for i in 0..vec_out.len(){
            vec_out[i] = 0_f64;
        }*/
        for i in 0..self.col_no.len() {
            vec_out[i] = 0_f64;
            for j in 0..self.col_no[i].len() {
                let col_no_val = self.col_no[i][j];
                let col_no_true_val = self.col_no_true[i][j];
                vec_out[i] = vec_out[i] + self.data[i][col_no_val] * vec_in[col_no_true_val];
            }
        }
    }

    fn multiply_par_split(&self, vec_in: &Vec<f64>, vec_out: &mut[f64], pool: &ThreadPool, maxthreads: usize, numthreads:usize,
        indices:&[usize])
    {
        if numthreads < maxthreads {
            let (vo1,vo2) = vec_out.split_at_mut(vec_out.len()/2);
            let (id1, id2) = indices.split_at(indices.len()/2);
            pool.install(||rayon::join(
                ||SparseMatrix::multiply_par_split(&self,vec_in, vo1, pool, maxthreads, numthreads*2, id1) ,
                ||SparseMatrix::multiply_par_split(&self,vec_in, vo2, pool, maxthreads, numthreads*2, id2)
            ));
        }
        else{
            for i in 0..indices.len(){
                vec_out[indices[i]] = 0_f64;
                for j in 0..self.col_no[indices[i]].len(){
                    let col_no_val = self.col_no[indices[i]][j];
                    let col_no_true_val = self.col_no_true[indices[i]][j];
                    vec_out[indices[i]] = vec_out[indices[i]] + self.data[indices[i]][col_no_val] * vec_in[col_no_true_val];
                }
            }
        }

    }




    /// multiplies with another matrix (multi column)
    fn dense_mat_mult(a: &Vec<Vec<f64>>, b: &Vec<Vec<f64>>, c: &mut Vec<Vec<f64>>) {
        for i in 0..a.len() {
            for j in 0..b[0].len() {
                c[i][j] = 0.0_f64;
            }
        }
        for i in 0..a.len() {
            for j in 0..b[0].len() {
                for k in 0..b.len() {
                    c[i][j] = c[i][j] + a[i][k] * b[k][j];
                }
            }
        }
    }

    /// multiplies with another matrix (multi column) generates output ... reference out not supplied
    pub(crate) fn dense_mat_mult_easy(a: &Vec<Vec<f64>>, b: &Vec<Vec<f64>>) ->Vec<Vec<f64>>{
        let mut c = Vec::with_capacity(a.len());
        for i in 0..a.len() {
            let mut ov2 = Vec::with_capacity(b[0].len());
            for j in 0..b[0].len() {
                ov2.push(0.0_f64);
            }
            c.push(ov2);
        }
        for i in 0..a.len() {
            for j in 0..b[0].len() {
                for k in 0..b.len() {
                    let zz = a[i][k];
                    let zz2 = b[k][j];
                    let cij = zz+zz2;
                    let dij = cij+1.0;
                    c[i][j] = c[i][j] + a[i][k] * b[k][j];
                }
            }
        }
        return c;
    }

    /// multiply dense matrix with vector
    pub(crate) fn dense_mat_mult_vec_easy(a: &Vec<Vec<f64>>, b: &Vec<f64>)->Vec<f64>{
        let mut vout = vec![0.0_f64;a.len()];
        for i in 0.. a.len(){
            for j in 0..b.len(){
                vout[i] = vout[i]+a[i][j]*b[j];
            }
        }
        return vout;
    }

    pub(crate) fn dense_mat_mult_easy_transpose(a: &Vec<Vec<f64>>, b: &Vec<Vec<f64>>) ->Vec<Vec<f64>>{
        let mut c = Vec::with_capacity(a.len());
        for i in 0..a[0].len() {
            let mut ov2 = Vec::with_capacity(b[0].len());
            for j in 0..b.len() {
                ov2.push(0.0_f64);
            }
            c.push(ov2);
        }
        for i in 0..a[0].len() {
            for j in 0..b.len() {
                for k in 0..b[0].len() {
                    c[i][j] = c[i][j] + a[k][i] * b[k][j];
                }
            }
        }
        return c;
    }


    pub(crate) fn create_zero_dense_matrix(nrows: usize, ncols: usize) ->Vec<Vec<f64>>{
        let mut matout = Vec::with_capacity(nrows);
        for i in 0..nrows{
            matout.push(Vec::with_capacity(ncols));
            for j in 0..ncols{
                matout[i].push(0.0_f64);
            }
        }
        return matout;
    }

    pub(crate) fn add_dense_matrix(a:&mut Vec<Vec<f64>>, b:&Vec<Vec<f64>>){
        for i in 0..a.len(){
            for j in 0..a[0].len(){
                a[i][j] = a[i][j]+b[i][j];
            }
        }
    }

    pub(crate) fn multiply_dense_by_constant(a:&mut Vec<Vec<f64>>, constant:f64){
        for i in 0..a.len(){
            for j in 0..a[0].len(){
                a[i][j] = a[i][j] * constant;
            }
        }
    }

    fn to_full(&self)->Vec<Vec<f64>>{
        let mut full_mat = Vec::with_capacity(self.nrows);
        for i in 0..self.nrows{
            let mut rw = Vec::with_capacity(self.ncols);
            for j in 0..self.ncols{
                rw.push(self.get(i,j));
            }
            full_mat.push(rw);
        }
        return full_mat;
    }

    pub(crate) fn solve_cg(&self, b: &Vec<f64>, x0: &Vec<f64>) ->Vec<f64>{
        let err = 1e-12;
        let mut u = x0.clone();
        let mut temp = vec![0_f64;u.len()];
        self.multiply(&u,&mut temp);
        let mut r0:Vec<f64> = (0..u.len()).map(|i|b[i]-temp[i]).collect();
        let norm_r0:f64 = (0..r0.len()).map(|i|r0[i]*r0[i]).sum();
        if norm_r0 < err{return u}
        let p0 = r0.clone();
        let mut k =0;
        let mut rk = r0.clone();
        let mut pk = p0.clone();

        while true{
            self.multiply(&pk, &mut temp);
            let temp1:f64 = (0..rk.len()).map(|i|rk[i]*rk[i]).sum();
            let temp2:f64 = (0..rk.len()).map(|i|pk[i]*temp[i]).sum();
            let alpha_k:f64 = temp1/temp2;
            for i in 0..u.len(){
                u[i] = u[i] + alpha_k * pk[i];
            }

            let rkp1 = (0..rk.len()).map(|i|rk[i]-alpha_k*temp[i]).collect::<Vec<f64>>();

            let norm_rk = (0..r0.len()).map(|i|rk[i]*rk[i]).sum::<f64>().sqrt();
            if norm_rk < err{
                //println!("converged in {} iterations",k);
                return u

            }
            //let temp1 = (0..rkp1.len()).map(|i| rkp1[i]*rkp1[i] ).sum::<f64>();
            let mut temp3 = 0.0_f64;
            let mut temp4 = 0.0_f64;
            for i in 0..rkp1.len(){
                temp3 = temp3+rkp1[i]*rkp1[i];
                temp4 = temp4+rk[i]*rk[i];
            }

            //let temp2:f64 = (0..rkp1.len()).map(|i|rk[i]*rk[i]).sum();
            let beta_k:f64 = temp3/temp4;
            let pkp1:Vec<f64>  = (0..rkp1.len()).map(|i|rkp1[i]+beta_k*pk[i]).collect::<Vec<f64>>();
            rk = rkp1.clone();
            pk = pkp1.clone();
            let rknorm:f64 = (0..rkp1.len()).map(|i|rk[i]*rk[i]).sum();
            //println!("{}",rknorm);
            k = k+1;
            if k>10000{
                println!("did not converge in {} iterations",k);
                let mut file = File::create("c:\\rustfiles\\Kg_bc_at_failure.csv").expect("can't create the file");
                let mut fbw = BufWriter::with_capacity(self.nrows*self.ncols*2,file);
                for i in 0..self.nrows{
                  for j in 0..self.ncols{
                      write!(fbw,"{},",self.get(i,j));
                  }
                    write!(fbw,"\n");
                }

                let mut lfile = File::create("c:\\rustfiles\\load_bc_at_failure.csv").expect("cant create load file");
                let mut lfbw = BufWriter::with_capacity(self.nrows*self.ncols*2, lfile);
                for i in 0..b.len(){
                    write!(lfbw, "{}\n",b[i]);
                }

                return u;
            }
        }
        return u;


    }
}