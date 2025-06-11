use bgzip::BGZFReader;
use core::fmt;
use std::io::BufRead;
use rust_htslib::bgzf::Writer;
use rust_htslib::tpool::ThreadPool;
use std::fs::File;
use std::collections::BinaryHeap;
use std::cmp::Reverse;

/// Aggregates multiple fragment files into a single file
/// This code is just a fancy implementation of the unix command `cat | sort -k1,1 -k2,2n -k3,3n | bgzip`
/// And might not be super efficient.
///
/// It would be better to make an implementation that makes use of the bgzip blocks and the fact that the files are already sorted
/// If someone wants and knows how to do that, please do!
use std::io::{Read as IoRead, Write};

fn read_fragments_file(file_name: &str, buffer: &mut String) {
    let f = File::open(file_name).unwrap_or_else(|_| panic!("Could not open file {}", file_name));
    let mut reader = BGZFReader::new(f)
        .unwrap_or_else(|_| panic!("Could not create BGZF reader for file {}", file_name));
    // Try to read file into buffer
    match reader.read_to_string(buffer) {
        Ok(_) => (),
        Err(_) => {
            println!("Could not read file {}, is it empty?", file_name);
        }
    };
}

/// Struct representing a fragment, used for sorting
///
/// # Fields
///
/// * `chrom` - Chromosome name.
/// * `start` - Start position.
/// * `end` - End position.
/// * `cell_barcode` - Cell barcode.
/// * `score` - Optional score.

#[derive(PartialEq, Eq)]
struct Fragment {
    chrom: String,
    start: usize,
    end: usize,
    cell_barcode: String,
    score: Option<usize>,
}

impl Fragment {
    /// Create a new Fragment from a string.
    ///
    /// # Arguments
    ///
    /// * `s` - String to parse.
    ///
    /// # Example
    ///
    /// ```rust
    /// let fragment = Fragment::new_from_string("chr1\t100\t200\tAACATCGATGGATG-1\t10");
    /// assert_eq!(fragment.chrom, "chr1");
    /// assert_eq!(fragment.start, 100);
    /// assert_eq!(fragment.end, 200);
    /// assert_eq!(fragment.cell_barcode, "AACATCGATGGATG-1");
    /// assert_eq!(fragment.score, Some(10));
    /// ```
    fn new_from_string(s: &str) -> Fragment {
        let fields: Vec<&str> = s.split('\t').collect();
        match fields.len() {
            4 => Fragment {
                chrom: fields[0].to_string(),
                start: fields[1].parse::<usize>().unwrap(),
                end: fields[2].parse::<usize>().unwrap(),
                cell_barcode: fields[3].to_string(),
                score: None,
            },
            5 => Fragment {
                chrom: fields[0].to_string(),
                start: fields[1].parse::<usize>().unwrap(),
                end: fields[2].parse::<usize>().unwrap(),
                cell_barcode: fields[3].to_string(),
                score: Some(fields[4].parse::<usize>().unwrap()),
            },
            _ => Fragment {
                chrom: fields[0].to_string(),
                start: fields[1].parse::<usize>().unwrap(),
                end: fields[2].parse::<usize>().unwrap(),
                cell_barcode: fields[3].to_string(),
                score: Some(fields[4].parse::<usize>().unwrap()),
            },
        }
    }
}

impl Ord for Fragment {
    fn cmp(&self, other: &Fragment) -> std::cmp::Ordering {
        let self_chrom = &self.chrom;
        let other_chrom = &other.chrom;

        let self_start = &self.start;
        let other_start = &other.start;

        let self_end = &self.end;
        let other_end = &other.end;

        let self_cell_barcode = &self.cell_barcode;
        let other_cell_barcode = &other.cell_barcode;

        if self_chrom != other_chrom {
            self_chrom.cmp(other_chrom)
        } else if self_start != other_start {
            return self_start.cmp(other_start);
        } else if self_end != other_end {
            return self_end.cmp(other_end);
        } else {
            return self_cell_barcode.cmp(other_cell_barcode);
        }
    }
}

impl PartialOrd for Fragment {
    fn partial_cmp(&self, other: &Fragment) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for Fragment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.score {
            Some(score) => write!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                self.chrom, self.start, self.end, self.cell_barcode, score
            ),
            None => write!(
                f,
                "{}\t{}\t{}\t{}",
                self.chrom, self.start, self.end, self.cell_barcode
            ),
        }
    }
}

/// Aggregates multiple fragment files into a single file.
///
/// # Arguments
/// * `path_to_fragment_files` - Paths to the fragment files.
/// * `path_to_output_file` - Path to the output file.
/// * `number_of_threads` - Number of threads to use for writing.
/// * `verbose` - Whether to print progress messages.

pub fn merge_fragment_files(
    path_to_fragment_files: &[String],
    path_to_output_file: &String,
    number_of_threads: u32,
    verbose: bool,
) {
    // initialize writer
    let tpool = ThreadPool::new(number_of_threads).unwrap_or_else(|_| {
        panic!(
            "Could not create thread pool with {} threads",
            number_of_threads
        )
    });
    let mut writer = Writer::from_path(path_to_output_file)
        .unwrap_or_else(|_| panic!("Could not open file {} for writing", path_to_output_file));
    writer
        .set_thread_pool(&tpool)
        .unwrap_or_else(|_| panic!("Could not set thread pool for file {}", path_to_output_file));

    let mut readers: Vec<(String, BGZFReader<File>)> = path_to_fragment_files
        .iter()
        .map(|path| {
            let f = File::open(path)
                .unwrap_or_else(|_| panic!("Could not open file {}", path));
            let reader = BGZFReader::new(f)
                .unwrap_or_else(|_| panic!("Could not create BGZF reader for file {}", path));
            (path.clone(), reader)
        })
        .collect();

    // Use a min-heap to keep track of the next fragment from each file
    // The heap will store (Reverse(Fragment), index_of_file)
    let mut heap = BinaryHeap::new();

    for i in 0..readers.len() {
        let (_, reader) = &mut readers[i];
        // Read the first line from each file
        let mut line_buffer = String::new();
        if let Ok(bytes_read) = reader.read_line(&mut line_buffer) {
            if bytes_read > 0 {
                let fragment = Fragment::new_from_string(line_buffer.trim());
                heap.push(Reverse((fragment, i)));
            }
        }
    }

    loop {
        // Get the smallest fragment from the heap
        let Reverse((min_fragment, file_idx)) = match heap.pop() {
            Some(f) => f,
            None => break, // All files are exhausted
        };

        // Write the smallest fragment
        writer.write_all(min_fragment.to_string().as_bytes()).unwrap();
        writer.write_all(b"\n").unwrap();

        // Read the next fragment from the file that `min_fragment` came from
        let (_, reader) = &mut readers[file_idx];
        let mut line_buffer = String::new();
        if let Ok(bytes_read) = reader.read_line(&mut line_buffer) {
            if bytes_read > 0 {
                let next_fragment = Fragment::new_from_string(line_buffer.trim());
                heap.push(Reverse((next_fragment, file_idx)));
            }
        }
    }
    writer.flush().unwrap();
    
}

fn log(message: &str, verbose: bool) {
    if verbose {
        println!("{}", message);
    }
}
