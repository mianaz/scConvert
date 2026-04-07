# scConvert 0.1.0.9001 (development)

> **Status:** unreleased; rolling development between 0.1.0 and the next tag.

## New features

- **Native Stereo-seq GEF reader (`LoadStereoSeqGef()`).** Pure-R reader for
  BGI `.gef` and `.cellbin.gef` files using `hdf5r` only. Handles both the
  square-bin and cell-bin schemas documented by STOmics. No Python,
  stereopy, or reticulate dependency. Spot coordinates are stored in
  `meta.data$spatial_x/y` and `misc$spatial_technology = "StereoSeq"`.
- **Native CosMx SMI reader (`LoadCosMx()`).** Thin R wrapper around
  `Seurat::LoadNanostring()` that validates the canonical flat-file bundle
  (`*exprMat*.csv`, `*metadata*.csv`, `*fov_positions*.csv`, `*tx_file*.csv`)
  and tags the result with `misc$spatial_technology = "CosMx"`. No squidpy
  or reticulate dependency.
- **CLI auto-delegation for vendor raw formats.** The `scconvert` C binary
  now auto-detects `.gef`, `.cellbin.gef`, and CosMx bundle directories and
  transparently delegates to the R backend via `Rscript` + `execvp()`. Users
  can write `scconvert mosta.gef mosta.h5ad` directly. `Rscript` lookup
  happens up-front; paths are absolutised; no shell-parsed command strings.
- **FOV round-trip through h5ad.** `writeH5AD()` now serializes
  `FOV@boundaries` and `FOV@molecules` into a stable
  `uns/spatial/{library}/segmentation/` and `uns/spatial/{library}/molecules/`
  contract. `readH5AD()` automatically rebuilds any FOV library it finds on
  load. Backward-compatible with squidpy and scanpy (they ignore unknown
  `uns/spatial/{lib}/` children).
- **CLI `varp` preservation.** The new `sc_stream_varp()` in `src/sc_groups.c`
  mirrors `sc_stream_obsp()` and maps `/varp/` (h5ad) to `/misc/__varp__/`
  (h5seurat). Handles both sparse (CSR/CSC group) and dense (array dataset)
  varp entries. Closes the manuscript limitation "CLI does not preserve
  varp".
- **Loom factor-level preservation.** `writeLoom()` now stores factor levels
  and `ordered` flag under `/scConvert_extensions/col_factor_levels/{name}`,
  outside `col_attrs` so loompy/scanpy continue to read the file without
  errors. `readLoom()` restores the factors on load.
- **h5mu per-modality `uns` round-trip.** `writeH5MU()` and `readH5MU()` now
  mirror each modality's `uns` group under
  `obj@misc[["__h5mu_uns_per_mod__"]][[modality]]` so per-modality uns
  entries survive h5mu round-trips instead of being flattened.
- **IMC multi-image support.** `SeuratSpatialToH5AD()` now iterates over
  every image in deterministic sorted order instead of processing only
  `Images()[1]`. Fixes the IMC 14/15 -> 11/13 double-roundtrip degradation
  documented in `NOTES.md` section 3.

## Robustness

- **C CLI memory-safety helpers.** New `sc_xmalloc()`, `sc_xcalloc()`,
  `sc_xrealloc()`, and `sc_check_mul_size()` in `src/sc_util.c` replace raw
  `malloc()` calls with overflow-checked allocations at the dense-embedding
  transpose sites in `src/sc_zarr.c:1043, 1668, 1674` and the column-buffer
  allocation in `src/sc_loom.c:138`. Prevents SIZE_MAX overflow on
  chip-scale embeddings.

## Bug fixes

- **`readH5AD()`: handle unsorted CSR column indices.** scanpy-processed files
  such as `scanpy.datasets.pbmc3k_processed()` ship with CSR matrices whose
  column indices are not sorted within each row. Previously these produced an
  invalid `dgCMatrix` on read. `readH5AD()` and `readH5MU()` now detect this
  condition via `.sort_dgc_indices()` in `R/LoadH5AD.R` (and the sibling helper
  in `R/LoadH5MU.R`) and sort indices column-wise before constructing the
  sparse matrix. A regression test has been added in
  `tests/testthat/test-regression-fixes.R`.
- **`H5SeuratToZarr()`: do not crash on 1D dense datasets.** The direct
  h5Seurat -> Zarr converter previously attempted to infer (rows, cols) from
  a 1D dense HDF5 dataset and crashed with a dimensionality error. It now
  emits a warning and skips the offending dataset. Tracked as "Chain D" in
  the benchmark manuscript; regression test included.
- **`writeZarr()`: skip scale.data layer when its shape differs from X.**
  `Seurat::ScaleData()` produces an `(n_hvg x n_cells)` matrix whose row
  count does not match `X`'s `n_genes`. AnnData requires all `layers/*` to
  match `X`'s shape, so `writeZarr()` now skips any layer whose dimensions
  differ from the default assay's data matrix. A regression test covers
  this case. (See `R/SaveZarr.R:131`.)

## Testing

- Added `tests/testthat/test-regression-fixes.R` with three regression tests
  pinning the bugs listed above. The test suite is now **140 `test_that`
  blocks** with **465 assertions** (up from 137 / 448).

## Documentation

- Added `NOTES.md` describing the state of recent bug fixes, benchmark
  findings that are specific to the manuscript (not to users), and
  investigations into dataset-pipeline issues (Stereo-seq MOSTA and CosMx
  SMI -- upstream format is not HDF5; see the notes file for details).

# scConvert 0.1.0

> **Release Date:** 2026-03-10

## Highlights

Initial public release of scConvert — a universal single-cell format converter for R.

### Universal Format Conversion
- Support for 7 formats: h5ad, h5Seurat, h5mu, Loom, Zarr, RDS, and SingleCellExperiment
- Hub architecture with 30+ conversion paths via `scConvert()`
- Direct HDF5 paths for h5ad/h5Seurat without intermediate loading

### Direct h5ad Loading
- `readH5AD()` for native h5ad-to-Seurat conversion without intermediate files
- Sparse (CSR/CSC) and dense matrix support
- Categorical metadata, dimensional reductions, neighbor graphs, and spatial data

### MuData (h5mu) Multimodal Support
- `readH5MU()` / `writeH5MU()` for multimodal single-cell data
- Automatic modality-to-assay name mapping (rna->RNA, prot->ADT, atac->ATAC)
- No MuDataSeurat or Python dependency required

### Zarr AnnData Support
- `readZarr()` and `writeZarr()` for Zarr-based AnnData stores (v2 format)
- Sparse CSR/CSC and dense matrix support
- Categorical metadata and dimensional reduction preservation

### Spatial Data (Visium)
- Bidirectional Visium spatial data conversion with image reconstruction
- Proper coordinate handling and scale factor preservation
- Compatible with scanpy/squidpy spatial analysis workflows

### C CLI Binary
- Standalone `scconvert` binary for h5ad/h5Seurat/h5mu conversions
- Streaming on-disk conversion without R or Python runtime
- Options: `--assay`, `--gzip`, `--overwrite`, `--quiet`

### BPCells On-Disk Loading
- `readH5AD(..., use.bpcells = TRUE)` for memory-efficient atlas-scale analysis
- Compatible with all Seurat analysis functions

### Seurat v5 Compatibility
- Full support for Seurat v5 Assay5 objects
- Proper handling of V5 layered data structure in conversions
