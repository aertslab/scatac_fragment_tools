use super::custom_errors;
use std::fmt;
use std::fmt::format;
use std::fs::File;
use std::path::Path;
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::cmp::Reverse;
use std::io::{BufRead, Write};
use pyo3::prelude::*;
use std::thread;
use noodles::{tabix, bgzf};
use noodles::csi::BinningIndex;
use noodles::core::{region::Interval, Position};

#[derive(Eq, PartialEq, Clone)]
struct GenomicRange {
    chromosome: String,
    start: usize,
    end: usize,
    cell_barcode: String,
    score: Option<usize>,
    file_index: usize,
    file_name: String
}

impl Ord for GenomicRange {
    fn cmp(&self, other: &GenomicRange) -> std::cmp::Ordering {
           self.start.cmp(&other.start)
            .then(self.end.cmp(&other.end))
    }
}

impl PartialOrd for GenomicRange {
    fn partial_cmp(&self, other: &GenomicRange) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl GenomicRange {
    fn new(line: &str, file_index: usize, filename: &str) -> Result<GenomicRange, custom_errors::InvalidFragmentFileError> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            return Err(custom_errors::InvalidFragmentFileError::new(filename));
        }
        Ok(GenomicRange{
            chromosome: fields[0].to_string(),
            start: fields[1].parse::<usize>()
                .map_err(|_| custom_errors::InvalidFragmentFileError::new(filename))?,
            end: fields[2].parse::<usize>()
                .map_err(|_| custom_errors::InvalidFragmentFileError::new(filename))?,
            cell_barcode: fields[3].to_string(),
            score: if fields.len() > 4 {
                    Some(fields[4].parse::<usize>()
                        .map_err(|_| custom_errors::InvalidFragmentFileError::new(filename))?
                    )
                } else {None},
            file_index,
            file_name: filename.to_string()
        })
    }
}

impl fmt::Display for GenomicRange {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.score {
            Some(score) => write!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                self.chromosome, self.start, self.end, self.cell_barcode, score
            ),
            None => write!(
                f,
                "{}\t{}\t{}\t{}",
                self.chromosome, self.start, self.end, self.cell_barcode
            ),
        }
    }
}


struct FragmentFileReader {
    reader: bgzf::Reader<File>,
    index: Option<tabix::Index>,
    fragment: Result<GenomicRange, custom_errors::InvalidFragmentFileError>,
    buffer: String,
    valid_cell_barcodes: Vec<String>,
    target_chrom: String,
    file_index: usize,
    file_name: String,
    at_end_of_file: bool,
}

impl FragmentFileReader {
    fn new(fragment_file_path: &str, valid_cell_barcodes: Vec<String>, target_chrom: String, file_index: usize) -> Result<FragmentFileReader, std::io::Error> {
        let reader = bgzf::Reader::new(File::open(Path::new(fragment_file_path))?);
        let mut index: Option<tabix::Index> = None;
        if Path::new(&format!("{}.tbi", fragment_file_path)).exists() {
            index = Some(tabix::fs::read(format!("{}.tbi", fragment_file_path))?);
        }
        Ok(FragmentFileReader{
            reader,
            index,
            fragment: GenomicRange::new("", file_index, fragment_file_path),
            buffer: String::new(),
            valid_cell_barcodes,
            target_chrom,
            file_index, 
            file_name: fragment_file_path.to_string(),
            at_end_of_file: false
        })
    }

    fn at_chrom(&self) -> bool {
        if let Ok(fragment) = &self.fragment{
            fragment.chromosome == self.target_chrom
        } else {
            false
        }
    }

    fn read_next(&mut self) -> Result<(), custom_errors::InvalidFragmentFileError>{
        self.buffer.clear();
        match self.reader.read_line(&mut self.buffer) {
            Ok(bytes_read) => {
                if bytes_read == 0 {
                    self.at_end_of_file = true
                }
            },
            Err(_) => {
                return Err(custom_errors::InvalidFragmentFileError::new(&self.file_name))
            }
        }
        self.buffer = self.buffer.trim().to_string();
        self.fragment = GenomicRange::new(&self.buffer, self.file_index, &self.file_name);
        Ok(())
    }

    fn skip_header(&mut self, pattern: &str) -> Result<(), custom_errors::InvalidFragmentFileError> {
        self.buffer = pattern.to_string();
        while self.buffer.starts_with(pattern) && !self.at_end_of_file {
            self.read_next()?;
        }
        Ok(())
    }

    fn skip_to_chromosome(&mut self, chromosome: &str) -> Result<(), custom_errors::InvalidFragmentFileError> {
        let mut using_index = false;
        if let Some(index) = &self.index {
            if let Some(header) = index.header() {
                let seq_names = header.reference_sequence_names();
                if let Some(chromosome_index) = seq_names.get_index_of(chromosome.as_bytes()) {
                    if let Ok(q) = index.query(
                        chromosome_index,
                        Interval::from(Position::MIN..)
                    ) {
                        let first_chunk = q[0];
                        self.reader.seek(first_chunk.start()).map_err(|_| custom_errors::InvalidFragmentFileError::new(&self.file_name))?;
                        self.read_next()?;
                        using_index = true;
                    } 
                }
            }
        } 
        if !using_index {
            // If no index, manually skip header lines.
            // otherwise this is done automatically
            self.skip_header("#")?;
        }
        while !self.at_end_of_file {
            if let Ok(fragment) = &self.fragment {
                if fragment.chromosome == chromosome {
                    return Ok(());
                }
            };
            self.read_next()?;
        }
        Ok(())
    }

    fn get_next_valid_fragment(&mut self) -> Result<Option<GenomicRange>, custom_errors::InvalidFragmentFileError> {
        while !self.at_end_of_file {
            // first check current fragment
            // it could be that after skip_header and skip_to_chromosome we are already at a 
            // valid fragment.
            if let Ok(fragment) = &self.fragment {
                if self.valid_cell_barcodes.contains(&fragment.cell_barcode) {
                    let fragment = fragment.clone();
                    if fragment.chromosome != self.target_chrom {
                        return Ok(None)
                    }
                    // read next fragment so next time we enter this function we have a new frag,et
                    self.read_next()?;
                    return Ok(Some(fragment));
                }
            }
            self.read_next()?;
        }
        Ok(None)
    }
}

fn split_fragments_by_cell_barcodes_for_chromosome(
    fragment_file_paths: &[&str],
    fragment_file_to_cell_barcode: &HashMap<String, Vec<String>>,
    chromosome: &str,
    gz_output_file: &mut bgzf::Writer<File>
) -> PyResult<()>{

    // Open fragment files which are gzipped, and pos-sorted.
    let mut readers: Vec<FragmentFileReader> = Vec::new();
    for (file_index, fragment_file_path) in fragment_file_paths.iter().enumerate() {
        if let Some(cell_barcodes) = fragment_file_to_cell_barcode.get(&fragment_file_path.to_string()) {
            readers.push(
                FragmentFileReader::new(
                    fragment_file_path,
                    cell_barcodes.to_vec(),
                    chromosome.to_string(),
                    file_index)?
            );
        }
    }

    // A binary heap will be used to write fragments in order from different files
    let mut heap = BinaryHeap::new();
    
    // go to correct chromosome and push first valid fragment to binary heap
    for reader in readers.iter_mut(){
        reader.skip_to_chromosome(chromosome)?;
        if !reader.at_end_of_file  && reader.at_chrom() {
            if let Some(fragment) = reader.get_next_valid_fragment()? {
                heap.push(Reverse(fragment));
            }
        }
    }

    // last start written
    // if we encounter a lower coordinate we know that the file is not sorted and we should crash

    let mut last_start_written: usize = 0;

    while let Some(Reverse(fragment)) = heap.pop() {
        if fragment.start < last_start_written {
            return Err(custom_errors::ValueError::new(format!("Fragment file: {} is not sorted!", fragment.file_name)).into());
        }
        gz_output_file.write_all(format!("{}\n", fragment).as_bytes())?;
        last_start_written = fragment.start;
        // read from file that currently has the smallest genomic range
        let reader = &mut readers[fragment.file_index];
        if !reader.at_end_of_file && reader.at_chrom() {
            if let Some(fragment) = reader.get_next_valid_fragment()? {
                heap.push(Reverse(fragment));
            }
        }
    }
    Ok(())
}

#[pyfunction]
pub fn split_fragment_files_by_cell_type(
    fragment_file_paths: Vec<String>,
    output_directory: &str,
    temp_directory: &str,
    cell_type_to_fragment_file_to_cell_barcode: HashMap<String, HashMap<String, Vec<String>>>,
    chromosomes: Vec<String>
) -> PyResult<()> {
    for cell_type in cell_type_to_fragment_file_to_cell_barcode.keys() {
        let mut handles: Vec<thread::JoinHandle<_>> = Vec::new();
        for chromosome in &chromosomes {
            let output_file_name = format!("{}/{}.{}.fragments.tsv.gz", temp_directory, cell_type, chromosome);
            let fragment_file_paths = fragment_file_paths.clone(); // Need to clone since threads take ownership
            let fragment_file_to_cell_barcode = cell_type_to_fragment_file_to_cell_barcode
                .get(cell_type)
                .unwrap()
                .clone();
            let file = File::create(output_file_name)?;
            let chromosome = chromosome.clone();
            let handle = thread::spawn(move || {
                let mut gz_output_file = bgzf::Writer::new(file);
                split_fragments_by_cell_barcodes_for_chromosome(
                    &fragment_file_paths.iter().map(|p| p.as_str()).collect::<Vec<_>>(),
                    &fragment_file_to_cell_barcode,
                    &chromosome,
                    &mut gz_output_file
                )
            });
            handles.push(handle);
        }
        for handle in handles {
            handle.join().expect("Thread panicked")?;
        }
        // concat all chromosomes
        let output_file_name = format!("{}/{}.fragments.tsv.gz", output_directory, cell_type);
        let output_file = File::create(&output_file_name)?;
        let mut writer = std::io::BufWriter::new(output_file);
        for chromosome in &chromosomes {
            let input_file_name = format!("{}/{}.{}.fragments.tsv.gz", temp_directory, cell_type, chromosome);
            let mut input_file = File::open(&input_file_name)?;
            std::io::copy(&mut input_file, &mut writer)?;
        }
        writer.flush()?;
    }
    Ok(())
}
