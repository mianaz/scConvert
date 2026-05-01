//! Sparse matrix streaming between HDF5 files.
//!
//! A scanpy/anndata CSR matrix has identical on-disk bytes to a
//! SeuratDisk/h5seurat CSC of the transposed shape — the three arrays
//! (`data`, `indices`, `indptr`) encode both interpretations. So copying
//! one into the other is a straight stream copy, with only the
//! shape-or-dims attribute flipped. No values touched, no transpose,
//! no sort.
//!
//! This module handles *only* the array copy. Format-specific dressings
//! (group attributes, layer paths, assay names, …) live in
//! `routes::h5ad_to_h5seurat` etc.

use crate::{io::hdf5_err, ConvertError, Result};
use hdf5_metno::{Dataset, Group};

/// Shape of a CSR/CSC group: `[n_rows, n_cols]`.
#[derive(Debug, Clone, Copy)]
pub struct MatrixShape {
    pub n_rows: i64,
    pub n_cols: i64,
    pub nnz: usize,
}

/// Peek at `data`/`indices`/`indptr` inside `src` and return the matrix
/// shape without reading the payload. Prefers the CSR/CSC group `shape`
/// attribute; falls back to `dims` (h5seurat convention); else derives
/// from `indptr.len() - 1` and `indices.max() + 1` (last-resort path,
/// not free — avoids scanning when the attr is present).
pub fn read_shape(src: &Group) -> Result<MatrixShape> {
    let data = src.dataset("data").map_err(hdf5_err)?;
    let indptr = src.dataset("indptr").map_err(hdf5_err)?;
    let nnz = data_size(&data)?;
    let indptr_len = data_size(&indptr)?;

    // shape attribute (h5ad convention).
    if let Ok(attr) = src.attr("shape") {
        let dims: Vec<i64> = attr.read_raw().map_err(hdf5_err)?;
        if dims.len() == 2 {
            return Ok(MatrixShape {
                n_rows: dims[0],
                n_cols: dims[1],
                nnz,
            });
        }
    }
    // dims attribute (h5seurat convention).
    if let Ok(attr) = src.attr("dims") {
        let dims: Vec<i64> = attr.read_raw().map_err(hdf5_err)?;
        if dims.len() == 2 {
            return Ok(MatrixShape {
                n_rows: dims[0],
                n_cols: dims[1],
                nnz,
            });
        }
    }
    // Inferred — expensive; acceptable because this only fires on
    // unusual files where both attrs are missing.
    let n_rows = indptr_len.saturating_sub(1) as i64;
    let indices = src.dataset("indices").map_err(hdf5_err)?;
    let n_cols = indices_max_plus_one(&indices)?;
    Ok(MatrixShape {
        n_rows,
        n_cols,
        nnz,
    })
}

/// Stream a CSR-style sparse group `src` (with `data`, `indices`,
/// `indptr` datasets) into a newly-created group `dst_name` under
/// `dst_parent`. Bytes are copied in fixed-size chunks so peak memory is
/// bounded regardless of `nnz`.
///
/// `compression` is a 0–9 gzip level; `None` disables compression.
///
/// This copy is agnostic to how the caller labels the result (CSR vs CSC,
/// h5ad vs h5seurat). The caller is responsible for writing the shape /
/// dims attribute with the right orientation.
pub fn stream_csr_group(
    src: &Group,
    dst_parent: &Group,
    dst_name: &str,
    compression: Option<u8>,
    stream_chunk_bytes: usize,
) -> Result<MatrixShape> {
    let src_data = src.dataset("data").map_err(hdf5_err)?;
    let src_indices = src.dataset("indices").map_err(hdf5_err)?;
    let src_indptr = src.dataset("indptr").map_err(hdf5_err)?;

    let shape = read_shape(src)?;
    let nnz = shape.nnz;
    let indptr_len = data_size(&src_indptr)?;

    let dst_group = dst_parent.create_group(dst_name).map_err(hdf5_err)?;

    // On-disk chunking: match the source when possible, else use 64 Ki.
    let disk_chunk = src_data
        .chunk()
        .and_then(|v| v.into_iter().next())
        .unwrap_or(65_536)
        .min(nnz.max(1));

    let dst_data = build_1d::<f64>(&dst_group, nnz, disk_chunk, compression)
        .create("data")
        .map_err(hdf5_err)?;
    let dst_indices = build_1d::<i32>(&dst_group, nnz, disk_chunk, compression)
        .create("indices")
        .map_err(hdf5_err)?;

    let chunk_elems_f64 = (stream_chunk_bytes / 8).max(disk_chunk);
    let chunk_elems_i32 = (stream_chunk_bytes / 4).max(disk_chunk);
    stream_1d::<f64>(&src_data, &dst_data, nnz, chunk_elems_f64)?;
    stream_1d::<i32>(&src_indices, &dst_indices, nnz, chunk_elems_i32)?;

    // indptr is (n_rows + 1) int32; small enough to read in full.
    let indptr_vals: Vec<i32> = src_indptr.read_raw().map_err(hdf5_err)?;
    if indptr_vals.len() != indptr_len {
        return Err(ConvertError::MalformedInput {
            reason: format!(
                "indptr length mismatch: metadata says {indptr_len}, \
                 read {} values",
                indptr_vals.len()
            ),
        });
    }
    let indptr_chunk = disk_chunk.min(indptr_len.max(1));
    let dst_indptr = build_1d::<i32>(&dst_group, indptr_len, indptr_chunk, compression)
        .create("indptr")
        .map_err(hdf5_err)?;
    dst_indptr.write(&indptr_vals).map_err(hdf5_err)?;

    Ok(shape)
}

// ---- internals ----

fn data_size(ds: &Dataset) -> Result<usize> {
    let shape = ds.shape();
    if shape.len() != 1 {
        return Err(ConvertError::MalformedInput {
            reason: format!("expected 1-D dataset, got shape {shape:?}"),
        });
    }
    Ok(shape[0])
}

fn indices_max_plus_one(ds: &Dataset) -> Result<i64> {
    let v: Vec<i64> = ds.read_raw().map_err(hdf5_err)?;
    let max = v.into_iter().max().unwrap_or(-1);
    Ok(max + 1)
}

fn build_1d<T: hdf5_metno::H5Type>(
    parent: &Group,
    dims: usize,
    chunk: usize,
    compression: Option<u8>,
) -> hdf5_metno::DatasetBuilderEmptyShape {
    let mut builder = parent.new_dataset::<T>().shape([dims]);
    if dims > 0 {
        builder = builder.chunk([chunk]);
        if let Some(level) = compression {
            builder = builder.deflate(level);
        }
    }
    builder
}

fn stream_1d<T>(src: &Dataset, dst: &Dataset, nnz: usize, chunk_elems: usize) -> Result<()>
where
    T: hdf5_metno::H5Type + Copy,
{
    if nnz == 0 {
        return Ok(());
    }
    let chunk = chunk_elems.max(1);
    let mut pos = 0usize;
    while pos < nnz {
        let end = (pos + chunk).min(nnz);
        let slice: ndarray::Array1<T> = src.read_slice_1d(pos..end).map_err(hdf5_err)?;
        dst.write_slice(&slice, pos..end).map_err(hdf5_err)?;
        pos = end;
    }
    Ok(())
}
