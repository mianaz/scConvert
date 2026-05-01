//! scconvert-core — native core for scConvert.
//!
//! B1 scope: data model, error taxonomy, fidelity report schema, and route
//! planner schema. No HDF5/Zarr IO yet. [`convert`] is a failing stub that
//! higher layers (CLI, C ABI, R/Python wrappers) can already link against.

mod error;
mod format;
pub mod io;
mod model;
mod report;
mod route;
mod routes;

pub use error::{ConvertError, Result};
pub use format::{detect_format, Format};
pub use model::{CategoricalMeta, DType, SparseLayout};
pub use report::{ArrayDigest, FidelityReport, FieldFate, FieldRecord, ReportFormat};
pub use route::{Route, RouteEstimate, RouteKind, RuntimeClass};

use std::path::Path;

/// Convert a file at `src` into `dst`, returning a fidelity report.
///
/// B2: h5ad → h5seurat is native (sparse X only, for now). Every other
/// format pair still returns [`ConvertError::NotYetImplemented`] and the
/// caller should fall back to the R/C CLI path.
pub fn convert(src: &Path, dst: &Path, options: &ConvertOptions) -> Result<FidelityReport> {
    let src_fmt = detect_format(src);
    let dst_fmt = detect_format(dst);
    match (src_fmt, dst_fmt) {
        (Format::H5ad, Format::H5seurat) => routes::h5ad_to_h5seurat::convert(src, dst, options),
        _ => Err(ConvertError::NotYetImplemented {
            kind: format!("convert({:?} -> {:?})", src_fmt, dst_fmt),
        }),
    }
}

/// Plan a conversion without executing it.
///
/// B1 stub — returns a rejected Route so callers can exercise the schema
/// end to end even before the planner is wired up.
pub fn plan(src: &Path, dst: &Path, _options: &ConvertOptions) -> Result<Route> {
    let src_fmt = detect_format(src);
    let dst_fmt = detect_format(dst);
    Ok(Route {
        src: src_fmt,
        dst: dst_fmt,
        kind: RouteKind::Rejected {
            reason: "route planner not implemented (B5)".into(),
        },
        estimate: RouteEstimate::default(),
    })
}

/// User-tunable conversion options, hydrated from CLI flags or a JSON blob
/// at the C ABI boundary.
#[derive(Debug, Default, Clone, serde::Serialize, serde::Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct ConvertOptions {
    /// Advisory memory cap in bytes. 0 = unset.
    pub max_memory_bytes: u64,
    /// Number of worker threads. 0 = auto.
    pub threads: u32,
    /// Override the destination format (useful for URI-style paths where
    /// the extension is ambiguous).
    pub dst_format: Option<Format>,
    /// If true, overwrite `dst` when it already exists.
    pub overwrite: bool,
    /// Compression level 0–9 (gzip/zstd depending on format).
    pub compression_level: Option<u8>,
    /// Populate the fidelity report's per-array digests. Off by default
    /// because digests require a second pass over the data.
    pub compute_digests: bool,
}
