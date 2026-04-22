# scConvert — Developer Notes

This file captures developer-facing context that is **not** meant for
users or CRAN reviewers but is important for future maintainers. Normal
user-facing changes belong in `NEWS.md`; the manuscript-specific
analysis lives in the companion repository under
`scConvert-manuscript/`.

Last updated: 2026-04-07

------------------------------------------------------------------------

## 1. Recent bug fixes (and their regression tests)

All three bugs below were discovered while preparing the Nature Methods
/ GB submission benchmark suite. Each has a matching test in
`tests/testthat/test-regression-fixes.R`; if one of those tests fails, a
fix has regressed.

### 1.1 Unsorted CSR column indices (`.sort_dgc_indices`)

**Symptom.** `readH5AD("pbmc3k_processed.h5ad")` (from
`scanpy.datasets.pbmc3k_processed()`) used to fail or return a
structurally invalid `dgCMatrix`: column indices within rows were not
monotonically increasing, violating the `dgCMatrix` contract (`@i` must
be sorted within each column defined by `@p`).

**Root cause.** scipy’s CSR/CSC matrices are **not** required to carry
sorted indices. `scanpy` sometimes ships files in that state. When we
reinterpret a CSR (cells x genes) as a CSC (genes x cells) during
[`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md),
the unsorted “column indices within a row” become unsorted “row indices
within a column” in the dgCMatrix view, which is a validity bug.

**Fix.** `LoadH5AD.R:95-148` defines `.sort_dgc_indices()` which (a)
detects whether any column has unsorted row indices via a fast scan, and
(b) sorts within each column if so. Only runs the sort if needed. Mirror
fix in `LoadH5MU.R:94-114`.

**Regression test.** `test-regression-fixes.R` constructs a 5x4 CSR with
deliberately shuffled indices, loads it, and verifies (1) successful
load, (2) correct numerical reconstruction, (3) all columns of the
resulting `dgCMatrix` have sorted `@i`.

**Status.** Confirmed fixed; test passes.

### 1.2 `H5SeuratToZarr` 1D dense dataset crash (“Chain D”)

**Symptom.** The direct h5Seurat -\> Zarr converter (invoked via
`scConvert("foo.h5seurat", "foo.zarr")` or the standalone
[`H5SeuratToZarr()`](https://mianaz.github.io/scConvert/reference/H5SeuratToZarr.md))
used to crash with a dimensionality mismatch error when encountering a
1D dense HDF5 dataset inside the h5Seurat structure.

**Root cause.** A 1D dense dataset cannot be reliably interpreted as a
(rows, cols) matrix without external shape metadata; the converter was
attempting a blind reshape.

**Fix.** `R/Convert.R:6489-6493` – warns and skips the offending dataset
rather than crashing. The caller loses that metadata field but gets a
valid Zarr store.

**Regression test.** Tiny Seurat object -\> writeH5Seurat -\>
H5SeuratToZarr -\> readZarr; asserts no crash and `ncol`/`nrow` match.

**Status.** Confirmed fixed; test passes. Referenced in the manuscript
Limitations as “the h5Seurat -\> Zarr direct converter had a
dimensionality mismatch bug at the time of initial benchmarking; this
has since been fixed and a regression test is included in the test
suite.”

### 1.3 `writeZarr` scale.data shape mismatch

**Symptom.** When writing a Seurat object that had been processed with
`ScaleData()`, the resulting Zarr store could not be read back by
`anndata.read_zarr()`: it raised
`ValueError: Value passed for key 'scale.data' is of incorrect shape`.

**Root cause.**
[`Seurat::ScaleData()`](https://satijalab.org/seurat/reference/ScaleData.html)
produces a matrix of shape `(n_hvg x n_cells)`, where `n_hvg` is usually
much smaller than `n_genes` (the default is 2000). AnnData requires
every entry in `layers/*` to match the shape of `X`, so writing
`scale.data` as a layer violates the AnnData contract.

**Fix.** `R/SaveZarr.R:131-133` skips any layer whose dimensions do not
match `X`. Emits an informational message.

**Regression test.** Creates a 50-gene x 30-cell Seurat object, runs
`NormalizeData` + `FindVariableFeatures(nfeatures=20)` + `ScaleData`,
confirms `dim(scale.data) == (20, 30)`, calls
[`writeZarr()`](https://mianaz.github.io/scConvert/reference/writeZarr.md),
and asserts (a) the store is created, (b) no `layers/scale.data`
directory exists, and (c) the store round-trips through
[`readZarr()`](https://mianaz.github.io/scConvert/reference/readZarr.md).

**Status.** Confirmed fixed; test passes.

------------------------------------------------------------------------

## 2. Dataset-pipeline investigation: Stereo-seq and CosMx (2026-04-07)

Both appeared as “failed” rows in `task5_spatial.csv` with the error
message `HDF5-API Errors: file signature not found`. This is **not** a
scConvert bug; it is a dataset download/manifest bug in the companion
benchmark repo.

**What we found.** The files
`scConvert-manuscript/data/tier2/spatial/stereoseq_mosta.h5ad` and
`scConvert-manuscript/data/tier2/spatial/cosmx_nsclc.h5ad` are **HTML
documents**, not HDF5 files. Their first eight bytes are
`3c 21 44 4f 43 54 59 50` (“\<!DOCTYP”) instead of the HDF5 magic number
`89 48 44 46 0d 0a 1a 0a`.

The cause: the original manifest entries pointed the HTTP downloader at
**landing pages** (`https://db.cngb.org/stomics/mosta/download/` and
`https://nanostring.com/.../ffpe-dataset/`), not direct file URLs. The
downloader saved the HTML response bodies under `.h5ad` extensions
without verifying the content type. scConvert was correct to reject
them.

**Upstream reality.**

- **Stereo-seq** native format is `.gef` (HDF5-based but with a
  Stereo-seq schema that is *not* AnnData-compatible) or `.gem` (plain
  text). The scverse-compatible path is to convert via
  `stereopy.io.stereo_to_anndata()`. CNGB hosts the MOSTA `.gef` files
  behind a free account.
- **CosMx SMI** native format is a flat-file bundle: cell metadata CSV,
  expression matrix CSV, FOV positions CSV, and tiff stacks. The
  canonical conversion to AnnData is `squidpy.read.nanostring()`.
  NanoString distributes the bundles behind a gated vendor portal.

**Fixes applied to the companion repo (not this package).**

1.  `benchmark/80_download_datasets.R` gains a `.is_hdf5_file()` helper
    that validates the HDF5 magic number on any downloaded file whose
    declared format is HDF5-based. The downloader now refuses to save
    non-HDF5 payloads under an `.h5ad` extension and emits a clear
    diagnostic (“looks like an HTML landing page; update the URL in the
    manifest”).
2.  `benchmark/00_config.R::check_dataset_availability()` also applies
    the signature check, so any pre-existing corrupt files are filtered
    out of downstream task runs even if they are not re-downloaded.
3.  Two new download methods `nanostring_convert` and
    `stereoseq_convert` implement the `squidpy.read.nanostring()` and
    `stereopy` conversion paths, respectively. Both expect a
    manually-staged raw bundle under `data/raw/{cosmx,stereoseq}/`
    because the upstream sources are gated.
4.  The dataset manifest entries for `stereoseq_mosta` and `cosmx_nsclc`
    have been rewritten to document all of the above and to point at the
    new conversion methods.
5.  The two corrupt HTML-as-h5ad files have been deleted.

**Manuscript consequence.** The spatial Results section now honestly
notes: “Two additional technologies (Stereo-seq MOSTA and CosMx SMI)
could not be processed because the public exemplar files are not stored
in HDF5; scConvert therefore cannot be evaluated on them without an
upstream format converter.”

**Take-away for maintainers.** If someone adds a new spatial or
multi-vendor benchmark dataset, verify **before** relying on it that (a)
the download URL is a direct file link, not a landing page, and (b) the
file passes `.is_hdf5_file()` or the equivalent signature test for its
declared format.

------------------------------------------------------------------------

## 3. Known limitations that are not bugs

Several items from the original limitation list have been addressed in
the development branch; what remains below is truly format-inherent or
deferred.

**Closed during the 0.1.0.9001 integration pass (2026-04-07):**

- **CLI varp preservation.** `sc_stream_varp()` now handles both sparse
  (CSR/CSC group) and dense (2D array) varp entries through h5ad \<-\>
  h5seurat streaming. See `src/sc_groups.c` and the CLI block in
  `tests/testthat/test-varp-roundtrip.R`.
- **IMC multi-image double-roundtrip.** `R/SpatialConversion.R:245` now
  sorts `Images()` and iterates every image instead of processing only
  `Images()[1]`. Synthetic regression test in
  `tests/testthat/test-imc-images-ordering.R`.
- **FOV segmentation/molecules silent loss.**
  [`WriteFOVToH5AD()`](https://mianaz.github.io/scConvert/reference/WriteFOVToH5AD.md)
  and
  [`ReadFOVFromH5AD()`](https://mianaz.github.io/scConvert/reference/ReadFOVFromH5AD.md)
  in `R/SpatialFOV.R` serialise FOV boundaries and molecules through
  `uns/spatial/{library}/segmentation/` and `molecules/` subgroups.
  Regression test in `tests/testthat/test-fov-roundtrip.R`.
- **Vendor raw formats (Stereo-seq GEF, CosMx SMI bundles).** Read-only
  support lives in `R/LoadStereoSeq.R` and `R/LoadCosMx.R`; the C CLI
  auto-detects them and delegates to the R backend via
  `execvp("Rscript")`.
- **Dense embedding overflow risk in CLI.** `sc_xmalloc()` +
  `sc_check_mul_size()` guard the zarr/loom transpose allocations in
  `src/sc_zarr.c:1043, 1668, 1674` and `src/sc_loom.c:138`.
- **Loom factor-level / h5mu per-modality `uns`.** `R/SaveLoom.R` and
  `R/LoadLoom.R` round-trip factor levels via
  `scConvert_extensions/col_factor_levels/{name}`.
  `R/Convert.R::writeH5MU` and `R/LoadH5MU.R` mirror per-modality uns
  via `misc[["__h5mu_uns_per_mod__"]]`.

**Still open / truly inherent:**

- **Non-Visium spatial fidelity remains \<100% for some fields.**
  Command logs and proprietary metadata that have no Seurat slot still
  vanish. The FOV fix above closes the silent-loss gap for CosMx/Xenium,
  but IMC panel metadata and 4i channel ordering remain format-level
  problems.
- **C CLI has no formal memory-safety audit.** `sc_xmalloc()` and
  overflow guards were added in the dominant dense-alloc sites, but a
  fuzzer pass (AFL++/libFuzzer on malformed h5ad) is still planned for
  pre-submission.
- **HDF5 2.0 forward compatibility untested.** The CLI depends on
  `H5Oget_info3` which is present since HDF5 1.10. HDF5 2.0 has not been
  released as of 2026-04-07 so no test is possible.
- **No BPCells streaming integration in the CLI.**
  `readH5AD(use.bpcells=TRUE)` works at the R level but the CLI does not
  produce BPCells-compatible output.

------------------------------------------------------------------------

## 4. Manuscript ↔︎ package consistency

The manuscript (`scConvert-manuscript/scConvert.tex`) contains several
specific numerical claims that must stay in sync with package reality:

- **Test counts**: `Availability` section claims 140 `test_that` blocks,
  465 assertions, 19 cross-validation tests, 3 regression tests. The
  development branch (0.1.0.9001) now has **153 test_that blocks** and
  **621 assertions** after the vendor-format + limitation sweep; update
  the manuscript before the next submission. Re-run
  `Rscript benchmark/94_audit_manuscript.R` in the companion repo.
- **Direct converter count**: `Methods` claims “36 direct `AToB()`
  converter functions”. If you add or remove exports of the form
  `AToB()`, update the manuscript.
- **CLI format list**: “h5ad, h5Seurat, h5mu, Loom, Zarr”. Should stay
  in sync with `src/main.c`.
- **Seurat v5 path adaptation**: manuscript describes reading from
  `assays/{assay}/data` (v4) and `assays/{assay}/layers/data` (v5). The
  logic lives in `R/h5Seurat.R:406-416`.

The audit script `scConvert-manuscript/benchmark/94_audit_manuscript.R`
extracts numbers from the benchmark CSVs and greps the manuscript for
stale markers. It exits non-zero if any known-bad phrase is detected.
Run `make audit` in the companion repo before every submission.
