use pyo3::prelude::*;
use std::os::raw::c_char;
use std::str;

enum AllowedSpecies {
    human,
    mouse,
}



/// Accepts a query sequence and a species from the Python caller.
/// Returns a string containing the IMG numbering.
#[pyfunction]
fn align_imgt(query_seq_ptr: &str, species_ptr: &str) -> PyResult<String> {
    let query_bytes = unsafe { CStr::from_ptr(query_seq_ptr).to_bytes() };
    let query_seq = str::from_utf8(query_bytes).unwrap();
    let species_bytes = unsafe { CStr::from_ptr(species_ptr).to_bytes() };
    let species = str::from_utf8(species_bytes).unwrap();

    match species {
        human => ,
        mouse => ,
        _ => Err(exceptions::PyValueError::new_err("Unrecognized species supplied.")),
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn antpack(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(align_query_sequence, m)?)?;
    Ok(())
}
