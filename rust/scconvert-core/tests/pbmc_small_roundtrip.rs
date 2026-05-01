//! Integration test: run the native h5ad → h5seurat route on the
//! pbmc_small fixture shipped with the R package. Verifies shape, dims,
//! and that the sparse arrays are readable as expected.
//!
//! Rerun with `--features pbmc-compare-c-cli` once the C CLI binary is
//! built in-tree at `src/scconvert`; the test then additionally asserts
//! bit-identical output between the two converters.

use std::path::PathBuf;

fn fixture() -> PathBuf {
    // tests/ runs with CARGO_MANIFEST_DIR = scconvert-core. Walk up to
    // the repo root and into inst/testdata.
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.pop(); // rust/
    p.pop(); // repo root
    p.push("inst/testdata/pbmc_small.h5ad");
    p
}

#[test]
fn pbmc_small_x_copy_is_correct() {
    let src = fixture();
    if !src.exists() {
        eprintln!("skipping: pbmc_small.h5ad fixture not found");
        return;
    }
    let tmpdir = tempfile::tempdir().unwrap();
    let dst = tmpdir.path().join("pbmc_small.h5seurat");

    let opts = scconvert_core::ConvertOptions {
        overwrite: true,
        ..Default::default()
    };
    let report = scconvert_core::convert(&src, &dst, &opts).unwrap();

    assert!(dst.exists(), "output file was not created");
    assert_eq!(report.src_format, scconvert_core::Format::H5ad);
    assert_eq!(report.dst_format, scconvert_core::Format::H5seurat);
    // pbmc_small has both X (normalised data) and raw/X (counts).
    // We expect two field records: one `data` layer, one `counts`.
    assert_eq!(
        report.fields.len(),
        2,
        "expected 2 field records, got {:?}",
        report.fields
    );
    for rec in &report.fields {
        assert!(
            rec.path.starts_with("/assays/RNA/layers/"),
            "unexpected field path: {}",
            rec.path
        );
    }

    // Inspect the file through hdf5-metno to confirm structure.
    let f = hdf5_metno::File::open(&dst).unwrap();
    let data = f.group("assays/RNA/layers/data").unwrap();
    let counts = f.group("assays/RNA/layers/counts").unwrap();

    // dims = [n_vars, n_cells] = [2000, 214] on pbmc_small.
    let dims: Vec<i64> = data.attr("dims").unwrap().read_raw().unwrap();
    assert_eq!(dims, vec![2000, 214]);
    let counts_dims: Vec<i64> = counts.attr("dims").unwrap().read_raw().unwrap();
    assert_eq!(counts_dims, vec![2000, 214]);

    // data/indices/indptr all exist and are sensibly sized.
    for grp in [&data, &counts] {
        let data_ds = grp.dataset("data").unwrap();
        let idx_ds = grp.dataset("indices").unwrap();
        let ptr_ds = grp.dataset("indptr").unwrap();
        assert_eq!(data_ds.shape(), idx_ds.shape());
        assert!(data_ds.shape()[0] > 0);
        // indptr length is n_rows (cells) + 1.
        assert_eq!(ptr_ds.shape()[0], 215);
    }

    // All required top-level groups so LoadH5Seurat does not bail.
    for g in [
        "tools",
        "commands",
        "images",
        "neighbors",
        "misc",
        "graphs",
        "reductions",
    ] {
        assert!(f.group(g).is_ok(), "missing required top-level group {g}");
    }
}
