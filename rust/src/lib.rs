mod aggregate_fragments;
mod split_fragments;

use pyo3::prelude::*;
use std::collections::HashMap;

/// Split fragments by cell barcode.
///
/// # Arguments
///
/// * `path_to_fragments` - Path to the fragments file.
/// * `path_to_output_folder` - Path to the output folder,
///    one file per cell type will be written here and the cell type name will be used as the filename.
///    If there are no fragments for a cell type, no file will be written for that cell type.
/// * `cell_type_to_cell_barcodes` - A HashMap mapping cell types to cell barcodes.
/// * `chromsizes` - A HashMap mapping chromosome names to chromosome sizes.
/// * `verbose` - Whether to print progress messages.
///
/// # Example
///
/// ```python
/// from scatac_fragment_tools import _rust_scatac_fragment_tools
/// rust_scatac_fragment_tools.split_fragments_by_cell_barcode(
///     path_to_fragments="fragments.tsv.gz",
///     path_to_output_folder="fragments_by_cell_type",
///     cell_type_to_cell_barcodes={
///         "cell_type_1": ["AACATCGATGGATG-1", "AACATCGATGGTTG-1"],
///         "cell_type_2": ["TTGATCGATGGATG-1", "AACATCGCTAGATG-1"]
///     },
///     chromsizes={
///         "chr1": 248956422,
///         "chr2": 242193529
///     },
///     verbose=True
/// )
/// ```

#[pyfunction]
fn split_fragments_by_cell_barcode(
    path_to_fragments: String,
    path_to_output_folder: String,
    cell_type_to_cell_barcodes: HashMap<String, Vec<String>>,
    chromsizes: HashMap<String, u64>,
    verbose: bool,
) -> PyResult<()> {
    // Invert cell_type_to_cell_barcodes
    let mut cell_barcode_to_cell_type: HashMap<String, Vec<String>> = HashMap::new();
    for (cell_type, cell_barcodes) in cell_type_to_cell_barcodes.iter() {
        for cell_barcode in cell_barcodes.iter() {
            cell_barcode_to_cell_type
                .entry(cell_barcode.to_string())
                .or_insert(Vec::new())
                .push(cell_type.to_string());
        }
    }
    split_fragments::split_fragments_by_cell_barcode(
        &path_to_fragments,
        &path_to_output_folder,
        cell_barcode_to_cell_type,
        chromsizes,
        5,
        verbose,
    );
    Ok(())
}

/// Merge fragment files.
///
/// # Arguments
///
/// * `path_to_fragment_files` - Paths to the fragment files.
/// * `path_to_output_file` - Path to the output file.
/// * `number_of_threads` - Number of threads to use for writing.
/// * `verbose` - Whether to print progress messages.
///
/// # Example
///
/// ```python
/// from scatac_fragment_tools import _rust_scatac_fragment_tools
/// rust_scatac_fragment_tools.merge_fragment_files(
///     path_to_fragment_files=[
///         "fragments_by_cell_type/cell_type_1.fragments.tsv.gz",
///         "fragments_by_cell_type/cell_type_2.fragments.tsv.gz"
///     ],
///     path_to_output_file="merged_fragments.tsv.gz",
///     number_of_threads=5,
///     verbose=True
/// )
/// ```

#[pyfunction]
fn merge_fragment_files(
    path_to_fragment_files: Vec<String>,
    path_to_output_file: String,
    number_of_threads: u32,
    verbose: bool,
) -> PyResult<()> {
    aggregate_fragments::merge_fragment_files(
        &path_to_fragment_files,
        &path_to_output_file,
        number_of_threads,
        verbose,
    );
    Ok(())
}

#[pymodule]
fn _rust_scatac_fragment_tools(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // set version dunder
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    // add functions
    m.add_function(wrap_pyfunction!(split_fragments_by_cell_barcode, m)?)?;
    m.add_function(wrap_pyfunction!(merge_fragment_files, m)?)?;
    Ok(())
}
