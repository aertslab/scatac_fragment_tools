use pyo3::prelude::*;

#[pymodule]
fn _rust_scatac_fragment_tools(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // set version dunder
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}