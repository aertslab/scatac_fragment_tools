use bgzip::BGZFReader;
use rust_htslib::bgzf::Writer;
use rust_htslib::tpool::ThreadPool;
use std::fs::File;
/// Aggregates multiple fragment files into a single file
/// This code is just a fancy implementation of the unix command `cat | sort -k1,1 -k2,2n -k3,3n | bgzip`
/// And might not be super efficient.
///
/// It would be better to make an implementation that makes use of the bgzip blocks and the fact that the files are already sorted
/// If someone wants and knows how to do that, please do!
use std::io::{Read as IoRead, Write};

fn read_fragments_file(file_name: &String, buffer: &mut String) {
    let f = File::open(file_name).expect(&format!("Could not open file {}", file_name));
    let mut reader = BGZFReader::new(f).expect(&format!(
        "Could not create BGZF reader for file {}",
        file_name
    ));
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

#[derive(PartialEq, Eq, PartialOrd)]
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
        let fields: Vec<&str> = s.split("\t").collect();
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
            _ => panic!("Invalid number of fields in fragment file!"),
        }
    }

    /// Convert a Fragment to a string.
    ///
    /// # Example
    ///
    /// ```rust
    /// let fragment = Fragment::new_from_string("chr1\t100\t200\tAACATCGATGGATG-1\t10");
    /// assert_eq!(fragment.to_string(), "chr1\t100\t200\tAACATCGATGGATG-1\t10");
    /// ```
    fn to_string(&self) -> String {
        match self.score {
            Some(score) => format!(
                "{}\t{}\t{}\t{}\t{}",
                self.chrom, self.start, self.end, self.cell_barcode, score
            ),
            None => format!(
                "{}\t{}\t{}\t{}",
                self.chrom, self.start, self.end, self.cell_barcode
            ),
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

        if self_chrom != other_chrom {
            return self_chrom.cmp(other_chrom);
        } else if self_start != other_start {
            return self_start.cmp(other_start);
        } else {
            return self_end.cmp(other_end);
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
    path_to_fragment_files: &Vec<String>,
    path_to_output_file: &String,
    number_of_threads: u32,
    verbose: bool,
) {
    // initialize writer
    let tpool = ThreadPool::new(number_of_threads).expect(&format!(
        "Could not create thread pool with {} threads",
        number_of_threads
    ));
    let mut writer = Writer::from_path(&path_to_output_file).expect(&format!(
        "Could not open file {} for writing",
        path_to_output_file
    ));
    writer.set_thread_pool(&tpool).expect(&format!(
        "Could not set thread pool for file {}",
        path_to_output_file
    ));

    // initialize buffer
    let mut buffer = String::new();

    // read all files into buffer
    for path_to_fragment_file in path_to_fragment_files.iter() {
        log(&format!("Reading file {}", path_to_fragment_file), verbose);
        read_fragments_file(path_to_fragment_file, &mut buffer);
    }

    // split buffer and remove empty lines
    let mut fragments: Vec<Fragment> = buffer
        .split("\n")
        .filter(|s| !s.is_empty())
        .map(|s| Fragment::new_from_string(s))
        .collect();

    // sort fragments
    log("Sorting fragments", verbose);
    fragments.sort();

    // write fragments
    log("Writing fragments", verbose);
    for fragment in fragments {
        writer.write(fragment.to_string().as_bytes()).unwrap();
        writer.write(b"\n").unwrap();
    }
    writer.flush().unwrap();
}

fn log(message: &str, verbose: bool) {
    if verbose {
        println!("{}", message);
    }
}
