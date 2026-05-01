//! h5ad → h5seurat conversion.
//!
//! B2 first slice: copy the primary sparse `X` matrix as Seurat v5
//! `assays/{assay}/layers/{counts,data}`, write the minimal skeleton of
//! required top-level groups, and emit a best-effort fidelity report.
//!
//! Not yet wired: obs metadata, var features, obsm embeddings, obsp
//! graphs, raw counts, uns tree, dense X. Each lands in a follow-up
//! commit on this file.

use crate::io::{self, sparse};
use crate::report::{FidelityReport, FieldFate, FieldRecord};
use crate::{ConvertError, ConvertOptions, Format, Result};
use hdf5_metno::File;
use std::path::Path;
use std::time::Instant;

/// Default assay name — matches the R package and the C CLI.
const DEFAULT_ASSAY: &str = "RNA";

/// Top-level groups an h5seurat file must contain for `LoadH5Seurat`
/// to succeed, even when empty.
const REQUIRED_EMPTY_GROUPS: &[&str] = &[
    "tools",
    "commands",
    "images",
    "neighbors",
    "misc",
    "graphs",
    "reductions",
];

pub fn convert(src: &Path, dst: &Path, options: &ConvertOptions) -> Result<FidelityReport> {
    if !src.exists() {
        return Err(ConvertError::SourceMissing {
            path: src.display().to_string(),
        });
    }
    if dst.exists() && !options.overwrite {
        return Err(ConvertError::DestinationExists {
            path: dst.display().to_string(),
        });
    }

    let started = Instant::now();
    let mut report = FidelityReport::new(Format::H5ad, Format::H5seurat);

    // Atomic write: stream into a sibling temp directory, rename on success.
    // We use a TempDir + fixed filename rather than NamedTempFile because on
    // Windows, dropping a NamedTempFile and immediately re-creating at the
    // same path races against the OS's pending-delete state and trips
    // ERROR_SHARING_VIOLATION.
    let parent = dst
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    std::fs::create_dir_all(parent)?;
    let tmp_dir = tempfile::Builder::new()
        .prefix("scconvert-")
        .tempdir_in(parent)?;
    let tmp_path = tmp_dir.path().join("out.h5seurat");

    let src_file = File::open(src).map_err(io::hdf5_err)?;
    let dst_file = File::create(&tmp_path).map_err(io::hdf5_err)?;

    // -- 1. top-level h5seurat attributes --
    io::write_string_attr(&dst_file, "active.assay", DEFAULT_ASSAY)?;

    // -- 2. sparse X → assays/{assay}/layers/{counts, data} --
    let x_group = src_file
        .group("X")
        .map_err(|_| ConvertError::MalformedInput {
            reason: "h5ad source has no /X group".into(),
        })?;
    let has_raw = src_file.group("raw").is_ok();

    let assays = io::ensure_group(&dst_file, "assays")?;
    let assay_grp = io::ensure_group(&assays, DEFAULT_ASSAY)?;

    // Required assay-level attributes (matches src/sc_h5ad.c:50-63).
    let key = format!("{}_", DEFAULT_ASSAY.to_lowercase());
    io::write_string_attr(&assay_grp, "key", &key)?;
    io::write_string_attr(&assay_grp, "s4class", "SeuratObject::Assay5")?;

    let layers_grp = io::ensure_group(&assay_grp, "layers")?;

    // Stream-read budget per iteration. Clamped to [1 MiB, 64 MiB]: too
    // small and we do more HDF5 round-trips than necessary; too large and
    // the R fallback's memory guarantees no longer apply.
    let chunk_bytes = options
        .max_memory_bytes
        .clamp(1024 * 1024, 64 * 1024 * 1024) as usize;

    // Shape in h5ad is [n_obs, n_vars] (rows = cells). Seurat expects
    // [n_vars, n_obs] (features × cells). The on-disk CSR bytes are
    // already a valid CSC of the transposed shape — so we copy as-is
    // and just swap the shape when writing the `dims` attribute.
    let layer_name = if has_raw { "data" } else { "counts" };
    let shape = sparse::stream_csr_group(
        &x_group,
        &layers_grp,
        layer_name,
        options.compression_level,
        chunk_bytes,
    )?;
    let dims_xyswapped = [shape.n_cols, shape.n_rows];
    let layer = layers_grp.group(layer_name).map_err(io::hdf5_err)?;
    io::write_i64_array_attr(&layer, "dims", &dims_xyswapped)?;

    report.fields.push(FieldRecord {
        path: format!("/assays/{DEFAULT_ASSAY}/layers/{layer_name}"),
        fate: FieldFate::Transformed,
        note: Some(format!(
            "h5ad CSR {}x{} -> h5seurat CSC {}x{} (zero-copy, shape swap)",
            shape.n_rows, shape.n_cols, shape.n_cols, shape.n_rows
        )),
        digest: None,
    });

    // -- 3. If the source has raw/X, also emit counts for it. --
    if has_raw {
        if let Ok(raw_grp) = src_file.group("raw") {
            if let Ok(raw_x) = raw_grp.group("X") {
                let raw_shape = sparse::stream_csr_group(
                    &raw_x,
                    &layers_grp,
                    "counts",
                    options.compression_level,
                    chunk_bytes,
                )?;
                let raw_dims = [raw_shape.n_cols, raw_shape.n_rows];
                let counts = layers_grp.group("counts").map_err(io::hdf5_err)?;
                io::write_i64_array_attr(&counts, "dims", &raw_dims)?;
                report.fields.push(FieldRecord {
                    path: format!("/assays/{DEFAULT_ASSAY}/layers/counts"),
                    fate: FieldFate::Transformed,
                    note: Some("h5ad raw/X CSR -> h5seurat counts CSC (zero-copy)".into()),
                    digest: None,
                });
            }
        }
    }

    // -- 4. Empty required groups so LoadH5Seurat does not bail. --
    for name in REQUIRED_EMPTY_GROUPS {
        if dst_file.group(name).is_err() {
            dst_file.create_group(name).map_err(io::hdf5_err)?;
        }
    }

    // Warnings for fields not yet wired in this slice.
    for missing in [
        "obs", "var", "obsm", "obsp", "varm", "varp", "uns", "layers",
    ] {
        if src_file.group(missing).is_ok() {
            report.warnings.push(format!(
                "/{missing} present in source but not yet migrated by the \
                 Rust backend (falls back to the C CLI in the R package)"
            ));
        }
    }

    // Flush + close before the rename so readers see a complete file.
    drop(assay_grp);
    drop(assays);
    drop(layers_grp);
    drop(x_group);
    drop(src_file);
    dst_file.flush().map_err(io::hdf5_err)?;
    drop(dst_file);

    // Atomic rename over any existing dst (only reachable if overwrite
    // was set; we checked earlier).
    if dst.exists() {
        std::fs::remove_file(dst)?;
    }
    std::fs::rename(&tmp_path, dst)?;

    report.bytes_written = std::fs::metadata(dst).ok().map(|m| m.len());
    report.duration_ms = Some(started.elapsed().as_millis() as u64);
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use hdf5_metno::File;
    use ndarray::Array1;
    use tempfile::TempDir;

    /// Build a minimal synthetic h5ad with just a sparse CSR `/X`. The
    /// values are such that the roundtrip is easy to verify.
    fn synthetic_h5ad(path: &Path) {
        let f = File::create(path).unwrap();
        let x = f.create_group("X").unwrap();
        let data: Array1<f64> = Array1::from(vec![1.0, 2.0, 3.0, 4.0]);
        let indices: Array1<i32> = Array1::from(vec![0, 2, 1, 2]);
        let indptr: Array1<i32> = Array1::from(vec![0, 2, 2, 4]);
        x.new_dataset::<f64>()
            .shape([4])
            .chunk([4])
            .create("data")
            .unwrap()
            .write(&data)
            .unwrap();
        x.new_dataset::<i32>()
            .shape([4])
            .chunk([4])
            .create("indices")
            .unwrap()
            .write(&indices)
            .unwrap();
        x.new_dataset::<i32>()
            .shape([4])
            .chunk([4])
            .create("indptr")
            .unwrap()
            .write(&indptr)
            .unwrap();
        // shape attr: [n_obs=3 rows, n_vars=3 cols].
        let attr = x.new_attr::<i64>().shape([2]).create("shape").unwrap();
        attr.write(&[3i64, 3]).unwrap();
    }

    #[test]
    fn minimum_x_roundtrips() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("in.h5ad");
        let dst = dir.path().join("out.h5seurat");
        synthetic_h5ad(&src);

        let opts = ConvertOptions {
            overwrite: true,
            ..Default::default()
        };
        let report = convert(&src, &dst, &opts).unwrap();

        // Output file should exist and contain the counts layer.
        assert!(dst.exists());
        let f = File::open(&dst).unwrap();
        let counts = f.group("assays/RNA/layers/counts").unwrap();

        let data: Vec<f64> = counts.dataset("data").unwrap().read_raw().unwrap();
        assert_eq!(data, vec![1.0, 2.0, 3.0, 4.0]);

        let indices: Vec<i32> = counts.dataset("indices").unwrap().read_raw().unwrap();
        assert_eq!(indices, vec![0, 2, 1, 2]);

        let indptr: Vec<i32> = counts.dataset("indptr").unwrap().read_raw().unwrap();
        assert_eq!(indptr, vec![0, 2, 2, 4]);

        // dims must be swapped: h5ad was [3 obs, 3 vars] so Seurat
        // sees [3 features, 3 cells].
        let dims: Vec<i64> = counts.attr("dims").unwrap().read_raw().unwrap();
        assert_eq!(dims, vec![3, 3]);

        // key / s4class / active.assay present.
        let active: hdf5_metno::types::VarLenUnicode =
            f.attr("active.assay").unwrap().read_scalar().unwrap();
        assert_eq!(active.to_string(), "RNA");
        let assay = f.group("assays/RNA").unwrap();
        let key: hdf5_metno::types::VarLenUnicode =
            assay.attr("key").unwrap().read_scalar().unwrap();
        assert_eq!(key.to_string(), "rna_");

        // The layer moniker: no raw/X so we emitted `counts`.
        assert_eq!(
            report
                .fields
                .iter()
                .find(|r| r.path.ends_with("/counts"))
                .map(|r| r.fate),
            Some(FieldFate::Transformed)
        );
        assert!(report.duration_ms.is_some());
        assert!(report.bytes_written.is_some());
    }

    #[test]
    fn refuses_to_overwrite_without_flag() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("in.h5ad");
        let dst = dir.path().join("out.h5seurat");
        synthetic_h5ad(&src);
        std::fs::write(&dst, b"existing content").unwrap();

        let opts = ConvertOptions::default();
        let err = convert(&src, &dst, &opts).unwrap_err();
        assert!(matches!(err, ConvertError::DestinationExists { .. }));
    }
}
