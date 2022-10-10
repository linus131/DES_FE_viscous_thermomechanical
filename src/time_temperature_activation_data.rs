use std::fs::{File,write};
use std::io::BufWriter;
use std::io::Write;

/// TempTimeActivationData stores the temperatures in cells at different times. The temperatures are
/// updated according to updates_per_second.
pub struct TempTimeActivationData{
    pub cell_numbers: usize,
    pub temps:Vec<Vec<f64>>,
    pub updates_per_second:usize,
    pub update_times:Vec<f64>
}
impl TempTimeActivationData{
    pub fn new(cell_numbers: usize, updates_per_second: usize, model_runtime: f64)->TempTimeActivationData{
        let mut temps:Vec<Vec<f64>> = Vec::with_capacity(cell_numbers);
        let total_updates_number = (model_runtime as usize +1) * updates_per_second;
        let mut update_times:Vec<f64> = Vec::with_capacity(total_updates_number);

        for i in 0..cell_numbers{
            temps.push(Vec::new());
        }
        return TempTimeActivationData{
            cell_numbers,
            temps,
            updates_per_second,
            update_times,
        }
    }
    pub fn update_temp(&mut self, cell_no:usize, time:f64, temp:f64) {
        for i in 0..cell_no {
            self.temps[i].push(temp);
        }
        self.update_times.push(time);
    }

    pub fn write_temp_time_to_file(&self, filename:&str) {
        // first row time
        // second row to end row of temps... if cell is not activated... use NaN
        let file = File::create(filename).expect("can't open file");
        let mut file_buf = BufWriter::with_capacity(10000,file);
        for i in 0..self.update_times.len(){
            write!(file_buf, "{},", self.update_times[self.update_times.len()-1-i]);
        }
        write!(file_buf, "\n");
        for i in &self.temps{
            for j in 0..i.len() {
                write!(file_buf, "{},",i[i.len()-1-j]);
            }
            write!(file_buf,"\n");
        }
    }

}
