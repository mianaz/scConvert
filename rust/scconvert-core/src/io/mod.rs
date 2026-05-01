//! Thin helpers over `hdf5-metno`.
//!
//! Purpose: centralise the HDF5 patterns we use repeatedly — streaming a
//! 1-D array from one dataset to another, copying attributes, creating
//! a destination group with matching compression and chunking — so the
//! format-specific converters in `routes::` can stay readable.
//!
//! Nothing here knows about h5ad or h5seurat; format semantics live in
//! `routes::`.

pub mod sparse;

use crate::{ConvertError, Result};
use hdf5_metno::Group;

/// Read an HDF5 attribute as a `String`, handling both scalar and
/// 1-element-array shapes. AnnData writes scalar UTF-8; some h5seurat
/// variants wrap the same data in a length-1 array.
pub fn read_string_attr(location: &hdf5_metno::Location, name: &str) -> Result<Option<String>> {
    let attr = match location.attr(name) {
        Ok(a) => a,
        Err(_) => return Ok(None),
    };
    let dtype = attr.dtype().map_err(hdf5_err)?;
    let td = dtype.to_descriptor().map_err(hdf5_err)?;
    match td {
        hdf5_metno::types::TypeDescriptor::VarLenUnicode
        | hdf5_metno::types::TypeDescriptor::VarLenAscii
        | hdf5_metno::types::TypeDescriptor::FixedUnicode(_)
        | hdf5_metno::types::TypeDescriptor::FixedAscii(_) => {
            // `read_scalar::<hdf5_metno::types::VarLenUnicode>()` is the
            // canonical path; for fixed-length strings we fall back to
            // `read_raw`.
            match attr.read_scalar::<hdf5_metno::types::VarLenUnicode>() {
                Ok(v) => Ok(Some(v.to_string())),
                Err(_) => {
                    let buf: Vec<hdf5_metno::types::VarLenUnicode> =
                        attr.read_raw().map_err(hdf5_err)?;
                    Ok(buf.into_iter().next().map(|s| s.to_string()))
                }
            }
        }
        _ => Ok(None),
    }
}

/// Write a scalar UTF-8 string attribute on a location (group or dataset).
pub fn write_string_attr(location: &hdf5_metno::Location, name: &str, value: &str) -> Result<()> {
    use std::str::FromStr;
    let v = hdf5_metno::types::VarLenUnicode::from_str(value)
        .map_err(|e| ConvertError::Other(format!("attr {name}: {e}")))?;
    let attr = location
        .new_attr::<hdf5_metno::types::VarLenUnicode>()
        .shape(())
        .create(name)
        .map_err(hdf5_err)?;
    attr.write_scalar(&v).map_err(hdf5_err)?;
    Ok(())
}

/// Write a fixed-length int64 array attribute.
pub fn write_i64_array_attr(
    location: &hdf5_metno::Location,
    name: &str,
    values: &[i64],
) -> Result<()> {
    let attr = location
        .new_attr::<i64>()
        .shape([values.len()])
        .create(name)
        .map_err(hdf5_err)?;
    attr.write(values).map_err(hdf5_err)?;
    Ok(())
}

/// Create an empty group at `path` if it does not exist.
pub fn ensure_group(parent: &Group, path: &str) -> Result<Group> {
    match parent.group(path) {
        Ok(g) => Ok(g),
        Err(_) => parent.create_group(path).map_err(hdf5_err),
    }
}

pub(crate) fn hdf5_err(e: hdf5_metno::Error) -> ConvertError {
    ConvertError::Other(format!("hdf5: {e}"))
}
