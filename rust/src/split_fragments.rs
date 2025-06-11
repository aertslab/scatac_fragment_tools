use itertools::Itertools;
use rust_htslib::bgzf::Writer;
use rust_htslib::tbx::{self, Read as TbxRead};
use rust_htslib::tpool::ThreadPool;
use std::collections::HashMap;
/// Splits a tabix-index fragment file into multiple files based on cell type.
use std::io::Write;

/// A lazy BGZF writer that only opens the file when the first write is called.
///
/// # Fields
///
/// * `writer` - The BGZF writer.
/// * `path` - The path to the file.
/// * `tpool` - The thread pool to use for writing.
/// * `written` - Whether the file has been written to yet.
///
/// # Methods
///
/// * `new` - Creates a new LazyBgzfWriter.
/// * `write` - Opens the file, if it has not been opened yet, and writes the given bytes to it.

struct LazyBgzfWriter<'a> {
    writer: Option<Writer>,
    path: String,
    tpool: &'a ThreadPool,
    written: bool,
}

impl LazyBgzfWriter<'_> {
    /// Creates a new LazyBgzfWriter.
    ///
    /// # Arguments
    ///
    /// * `path` - The path to the file.
    /// * `tpool` - The thread pool to use for writing.

    fn new(path: String, tpool: &ThreadPool) -> LazyBgzfWriter {
        LazyBgzfWriter {
            writer: None,
            path,
            tpool,
            written: false,
        }
    }

    /// Opens the file, if it has not been opened yet, and writes the given bytes to it.
    ///
    /// # Arguments
    ///
    /// * `bytes` - The bytes to write.
    fn write_all(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        self.written = true;
        if self.writer.is_none() {
            let mut writer = Writer::from_path(&self.path)
                .unwrap_or_else(|_| panic!("Could not open file \"{}\" for writing", self.path));
            writer
                .set_thread_pool(self.tpool)
                .unwrap_or_else(|_| panic!("Could not set thread pool for \"{}\"", self.path));
            self.writer = Some(writer);
        }
        self.writer.as_mut().unwrap().write_all(bytes)
    }
}

fn sanitize_string_for_filename(s: String) -> String {
    s.replace([' ', '/'], "_")
}

/// Splits a tabix-index fragment file into multiple files based on cell type.
///
/// # Arguments
///
/// * `path_to_fragments` - Path to the fragments file.
/// * `path_to_output_folder` - Path to the output folder,
///     one file per cell type will be written here and the cell type name will be used as the filename.
///     If there are no fragments for a cell type, no file will be written for that cell type.
/// * `cell_barcode_to_cell_type` - A HashMap mapping cell barcodes to cell types.
/// * `chromsizes` - A HashMap mapping contig names to contig sizes.
/// * `number_of_threads` - Number of threads to use for writing.
/// * `verbose` - Whether to print progress messages.

pub fn split_fragments_by_cell_barcode(
    path_to_fragments: &String,
    path_to_output_folder: &String,
    cell_barcode_to_cell_type: HashMap<String, Vec<String>>,
    chromsizes: HashMap<String, u64>,
    number_of_threads: u32,
    verbose: bool,
) {
    // Initialize reader
    let mut tbx_reader = tbx::Reader::from_path(path_to_fragments)
        .unwrap_or_else(|_| panic!("Could not open file \"{}\"", path_to_fragments));

    // Initialize writers
    // Use lazy writer to avoid generating empty files
    let writer_tpool = ThreadPool::new(number_of_threads).unwrap_or_else(|_| {
        panic!(
            "Could not create thread pool with {} threads",
            number_of_threads
        )
    });
    let mut cell_type_to_writer: HashMap<&String, LazyBgzfWriter> = HashMap::new();
    let unique_cell_types: Vec<&String> = cell_barcode_to_cell_type
        .values()
        .flatten()
        .unique()
        .collect();
    for cell_type in unique_cell_types {
        let cell_type_name = sanitize_string_for_filename(cell_type.clone().to_string());
        let path_to_output = format!(
            "{}/{}.fragments.tsv.gz",
            path_to_output_folder, cell_type_name
        );
        let lazy_writer = LazyBgzfWriter::new(path_to_output, &writer_tpool);
        cell_type_to_writer.insert(cell_type, lazy_writer);
    }

    // initialize variables to store read data
    let mut read: Vec<u8> = Vec::new();

    let contigs_in_fragments_file = tbx_reader.seqnames();

    for contig in chromsizes.keys().sorted() {
        if !contigs_in_fragments_file.contains(contig) {
            log(
                &format!(
                    "Skipping contig \"{}\" because it is not in the fragments file",
                    contig
                ),
                verbose,
            );
            continue;
        }
        log(&format!("Processing contig \"{}\"", contig), verbose);
        // get contig id and size and fetch whole contig
        let contig_id = tbx_reader
            .tid(contig)
            .unwrap_or_else(|_| panic!("Could not get contig id for contig \"{}\"", contig));
        let contig_size = chromsizes.get(contig).unwrap();
        tbx_reader
            .fetch(contig_id, 0, *contig_size)
            .unwrap_or_else(|_| {
                panic!("Could not fetch contig \"{}\" from fragments file", contig)
            });

        // read first read of contig
        let mut not_at_end = tbx_reader
            .read(&mut read)
            .unwrap_or_else(|_| panic!("Could not read from fragments file"));
        let mut read_as_str = String::from_utf8(read.clone()).unwrap();

        // loop over reads
        while not_at_end {
            let read_cb = read_as_str.split('\t').nth(3).unwrap().to_string();
            if let Some(cell_types) = cell_barcode_to_cell_type.get(&read_cb) {
                for cell_type in cell_types {
                    let writer = cell_type_to_writer.get_mut(cell_type).unwrap();
                    writer.write_all(&read).unwrap_or_else(|_| {
                        panic!(
                            "Could not write contig \"{}\" to \"{}\" fragments file",
                            contig, &writer.path
                        )
                    });
                    writer.write_all(b"\n").unwrap_or_else(|_| {
                        panic!(
                            "Could not write contig \"{}\" to \"{}\" fragments file",
                            contig, &writer.path
                        )
                    });
                }
            }
            read.clear();
            not_at_end = tbx_reader.read(&mut read).unwrap();
            read_as_str = String::from_utf8(read.clone()).unwrap();
        }

        // flush buffers
        for writer in cell_type_to_writer.values_mut() {
            if writer.written {
                log(
                    &format!(
                        "Flush reads for contig \"{}\" to \"{}\" fragments file",
                        contig, writer.path
                    ),
                    verbose,
                );
                writer.writer.as_mut().unwrap().flush().unwrap_or_else(|_| {
                    panic!(
                        "Could not flush reads for contig \"{}\" to \"{}\" fragments file",
                        contig, &writer.path
                    )
                });
            }
        }
    }
}

fn log(message: &str, verbose: bool) {
    if verbose {
        println!("{}", message);
    }
}
