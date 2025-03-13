mod pseudobulk;
mod custom_errors;

use pyo3::prelude::*;


/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn _rust_scatac_fragment_tools(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // set version dunder
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    // add functions
    m.add_function(wrap_pyfunction!(pseudobulk::split_fragment_files_by_cell_type, m)?)?;
    Ok(())
}
