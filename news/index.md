# Changelog

## scConvert 0.2.1 (development)

### New features

- **C CLI binary now ships on Windows.** Previous releases built and
  bundled `scconvert` (Linux + macOS) only, and Windows users silently
  fell through to the slower R streaming path. The package’s CI now
  compiles a Windows `scconvert.exe` against MSYS2’s ucrt64 toolchain
  and bundles the four required runtime DLLs (`hdf5.dll`, `zlib1.dll`,
  `libgcc_s_seh-1.dll`, `libwinpthread-1.dll`) alongside the exe in
  `inst/bin/`. Windows resolves implicit-load DLLs from the exe’s own
  directory first, so users do not need MSYS2 or any HDF5 installation.
  Adds roughly 7 MB to the installed package on Windows.

### Bug fixes

- **[`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
  no longer aborts on `close_all()` errors.** The Windows CRAN binary of
  hdf5r 1.3.12 bundles HDF5 1.12.1, which refuses to close a file when
  any leaf-object ID is still held by the reference counter. The error
  propagated out of [`on.exit()`](https://rdrr.io/r/base/on.exit.html),
  surfaced as `scConvert_cli(h5ad -> rds)` returning `FALSE`, and broke
  a handful of Windows test expectations. The three `close_all()` sites
  in `R/LoadH5AD.R` are now wrapped in
  `tryCatch(..., error = function(e) NULL)`, matching the pattern
  already used in `R/WriteH5AD.R`. The read has fully completed by the
  time `close_all()` runs; cleanup failure is best-effort.

- **C reader: respect source string encoding (HDF5 \>= 2.0 compat).**
  `sc_get_str_attr`, `sc_get_str_array_attr`, `sc_copy_group_attrs`, the
  `copy_attr_cb` in `sc_dataframe.c`, the h5seurat-factor and
  h5ad-categorical readers, and the group-copy attribute paths in
  `sc_groups.c` all forced `sc_create_vlen_str_type()` (UTF-8 vlen) as
  the H5Aread/H5Dread memtype. hdf5r writes vlen strings with
  CSET=ASCII; HDF5 \>= 2.0 refuses the implicit ASCII\<-\>UTF-8
  conversion and the read silently fails. Symptom: on HDF5 2.x,
  [`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
  dropped every obs metadata column past `nFeature_RNA` because the
  `column-order` attribute read returned all empty strings. Fix: use the
  source attribute / dataset’s own datatype as the memtype and check
  H5Aread / H5Dread return values explicitly. CI did not catch this
  because the CI matrix uses HDF5 1.14, which permits the implicit
  conversion.

- **writeH5AD: normalize obsm embedding orientation.** Both
  `.writeH5AD_c` and
  [`DirectSeuratToH5AD()`](https://mianaz.github.io/scConvert/reference/DirectSeuratToH5AD.md)
  now defensively transpose reduction embeddings whose first dim doesn’t
  match `n_cells`, and skip with a warning if neither dim matches.
  Prevents the C writer from emitting `obsm` datasets with HDF5 shape
  `(n_dims, n_cells)`, which fails `anndata.read_h5ad()`. Test:
  `tests/testthat/test-regression-fixes.R` asserts
  `f['obsm/X_pca'].shape == (n_cells, n_dims)` via h5py for both writer
  paths.

### Documentation

- **pkgdown site no longer publishes internal developer notes.** Moved
  `NOTES.md` into `dev/NOTES.md` so pkgdown’s home-page scanner stops
  rendering `NOTES.html` onto the public site. Workflow deploy now uses
  `clean: true` to purge stale internal pages (`findings.html`,
  `NOTES.html`) from gh-pages.

## scConvert 0.2.0

> **Release Date:** 2026-05-04

### Breaking changes

- **HDF5 1.14+ required.** `SystemRequirements` bumped from
  `HDF5 (>= 1.10.0)` to `HDF5 (>= 1.14.0)`. The hdf5r close path
  segfaults on libhdf5 1.10.x when closing files with many open child
  IDs (groups/datasets opened via `[[]]` subsetting); this is upstream
  and not catchable in R. Modern Linux distros (Ubuntu 24.04+, Debian
  12+, Fedora 38+), recent macOS packages, and Bioconductor’s `Rhdf5lib`
  all ship 1.14+, so this matches real deployments.

### New features

#### Existing native readers (carried from development)

- **Native Stereo-seq GEF reader
  ([`LoadStereoSeqGef()`](https://mianaz.github.io/scConvert/reference/LoadStereoSeqGef.md)).**
  Pure-R reader for BGI `.gef` and `.cellbin.gef` files using `hdf5r`
  only. Handles both the square-bin and cell-bin schemas documented by
  STOmics. No Python, stereopy, or reticulate dependency. Spot
  coordinates are stored in `meta.data$spatial_x/y` and
  `misc$spatial_technology = "StereoSeq"`.
- **Native CosMx SMI reader
  ([`LoadCosMx()`](https://mianaz.github.io/scConvert/reference/LoadCosMx.md)).**
  Thin R wrapper around
  [`Seurat::LoadNanostring()`](https://satijalab.org/seurat/reference/ReadNanostring.html)
  that validates the canonical flat-file bundle (`*exprMat*.csv`,
  `*metadata*.csv`, `*fov_positions*.csv`, `*tx_file*.csv`) and tags the
  result with `misc$spatial_technology = "CosMx"`. No squidpy or
  reticulate dependency.
- **CLI auto-delegation for vendor raw formats.** The `scconvert` C
  binary now auto-detects `.gef`, `.cellbin.gef`, and CosMx bundle
  directories and transparently delegates to the R backend via
  `Rscript` + `execvp()`. Users can write
  `scconvert mosta.gef mosta.h5ad` directly. `Rscript` lookup happens
  up-front; paths are absolutised; no shell-parsed command strings.
- **FOV round-trip through h5ad.**
  [`writeH5AD()`](https://mianaz.github.io/scConvert/reference/writeH5AD.md)
  now serializes `FOV@boundaries` and `FOV@molecules` into a stable
  `uns/spatial/{library}/segmentation/` and
  `uns/spatial/{library}/molecules/` contract.
  [`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
  automatically rebuilds any FOV library it finds on load.
  Backward-compatible with squidpy and scanpy (they ignore unknown
  `uns/spatial/{lib}/` children).
- **CLI `varp` preservation.** The new `sc_stream_varp()` in
  `src/sc_groups.c` mirrors `sc_stream_obsp()` and maps `/varp/` (h5ad)
  to `/misc/__varp__/` (h5seurat). Handles both sparse (CSR/CSC group)
  and dense (array dataset) varp entries. Closes the manuscript
  limitation “CLI does not preserve varp”.
- **Loom factor-level preservation.**
  [`writeLoom()`](https://mianaz.github.io/scConvert/reference/writeLoom.md)
  now stores factor levels and `ordered` flag under
  `/scConvert_extensions/col_factor_levels/{name}`, outside `col_attrs`
  so loompy/scanpy continue to read the file without errors.
  [`readLoom()`](https://mianaz.github.io/scConvert/reference/readLoom.md)
  restores the factors on load.
- **h5mu per-modality `uns` round-trip.**
  [`writeH5MU()`](https://mianaz.github.io/scConvert/reference/writeH5MU.md)
  and
  [`readH5MU()`](https://mianaz.github.io/scConvert/reference/readH5MU.md)
  now mirror each modality’s `uns` group under
  `obj@misc[["__h5mu_uns_per_mod__"]][[modality]]` so per-modality uns
  entries survive h5mu round-trips instead of being flattened.
- **IMC multi-image support.**
  [`SeuratSpatialToH5AD()`](https://mianaz.github.io/scConvert/reference/SeuratSpatialToH5AD.md)
  now iterates over every image in deterministic sorted order instead of
  processing only `Images()[1]`. Fixes the IMC 14/15 -\> 11/13
  double-roundtrip degradation documented in `dev/NOTES.md` section 3.

### P0 robustness fixes (Codex review response, Part A)

- **SOMA / SpatialData generic dispatch.** Lambdas registered for the
  `scConvert` generic now accept `filename =` so
  [`scConvert.character()`](https://mianaz.github.io/scConvert/reference/scConvert.md)
  reaches the right method. Adds
  `tests/testthat/test-generic-dispatch.R`.
- **CLI build hygiene.** CLI `.o` files are isolated to `src/cli_obj/`;
  `make -f Makefile.cli install-bin` copies the binary to
  `inst/bin/scconvert`;
  [`sc_find_cli()`](https://mianaz.github.io/scConvert/reference/sc_find_cli.md)
  prefers it over the source-tree `src/scconvert`.
- **Bounded-memory R sparse streaming.**
  [`.stream_sparse_group()`](https://mianaz.github.io/scConvert/reference/dot-stream_sparse_group.md)
  now reads in 64 MiB chunks (tunable via
  `options("scConvert.stream_chunk_bytes")`) instead of materialising
  whole sparse matrices. Adds `tests/testthat/test-stream-memory.R`.
- **Canonical h5mu layout on write.**
  [`writeH5MU()`](https://mianaz.github.io/scConvert/reference/writeH5MU.md)
  now writes the top-level `/var` (concat of modalities),
  `/obsmap/{mod}` (always `0..n-1` for Seurat sources), and
  `/varmap/{mod}` (block-diagonal with `-1` sentinels) in the muon
  convention.
- **Atomic SOMA / SpatialData writes.** Both writers now build under a
  sibling temp name and rename on success, so a mid-write crash leaves
  the user’s existing path untouched.
  ([`writeSpatialData()`](https://mianaz.github.io/scConvert/reference/writeSpatialData.md)
  initially shipped with an `on.exit` that deleted its own
  freshly-renamed output; fixed in 2026-05-01 via a `write_succeeded`
  disarm flag.)
- **Python-validation tests in CI.** `tests/scverse-env.yml` +
  `setup-micromamba` make `tests/testthat/test-python-validation.R`
  runnable on the GitHub runner. Test no longer hardcodes the macOS
  conda path.

### Robustness

- **C CLI memory-safety helpers.** New `sc_xmalloc()`, `sc_xcalloc()`,
  `sc_xrealloc()`, and `sc_check_mul_size()` in `src/sc_util.c` replace
  raw `malloc()` calls with overflow-checked allocations at the
  dense-embedding transpose sites in `src/sc_zarr.c:1043, 1668, 1674`
  and the column-buffer allocation in `src/sc_loom.c:138`. Prevents
  SIZE_MAX overflow on chip-scale embeddings.
- **Defensive `close_all` wrap on direct-path conversion.** R/Convert.R
  wraps `hfile$close_all()` in `tryCatch` for HDF5 1.10.x graceful
  degradation. (1.14+ is required and tested; the wrap is no-op there.)

### Bug fixes

- **Native Stereo-seq GEF reader
  ([`LoadStereoSeqGef()`](https://mianaz.github.io/scConvert/reference/LoadStereoSeqGef.md)).**
  Pure-R reader for BGI `.gef` and `.cellbin.gef` files using `hdf5r`
  only. Handles both the square-bin and cell-bin schemas documented by
  STOmics. No Python, stereopy, or reticulate dependency. Spot
  coordinates are stored in `meta.data$spatial_x/y` and
  `misc$spatial_technology = "StereoSeq"`.
- **Native CosMx SMI reader
  ([`LoadCosMx()`](https://mianaz.github.io/scConvert/reference/LoadCosMx.md)).**
  Thin R wrapper around
  [`Seurat::LoadNanostring()`](https://satijalab.org/seurat/reference/ReadNanostring.html)
  that validates the canonical flat-file bundle (`*exprMat*.csv`,
  `*metadata*.csv`, `*fov_positions*.csv`, `*tx_file*.csv`) and tags the
  result with `misc$spatial_technology = "CosMx"`. No squidpy or
  reticulate dependency.
- **CLI auto-delegation for vendor raw formats.** The `scconvert` C
  binary now auto-detects `.gef`, `.cellbin.gef`, and CosMx bundle
  directories and transparently delegates to the R backend via
  `Rscript` + `execvp()`. Users can write
  `scconvert mosta.gef mosta.h5ad` directly. `Rscript` lookup happens
  up-front; paths are absolutised; no shell-parsed command strings.
- **FOV round-trip through h5ad.**
  [`writeH5AD()`](https://mianaz.github.io/scConvert/reference/writeH5AD.md)
  now serializes `FOV@boundaries` and `FOV@molecules` into a stable
  `uns/spatial/{library}/segmentation/` and
  `uns/spatial/{library}/molecules/` contract.
  [`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
  automatically rebuilds any FOV library it finds on load.
  Backward-compatible with squidpy and scanpy (they ignore unknown
  `uns/spatial/{lib}/` children).
- **CLI `varp` preservation.** The new `sc_stream_varp()` in
  `src/sc_groups.c` mirrors `sc_stream_obsp()` and maps `/varp/` (h5ad)
  to `/misc/__varp__/` (h5seurat). Handles both sparse (CSR/CSC group)
  and dense (array dataset) varp entries. Closes the manuscript
  limitation “CLI does not preserve varp”.
- **Loom factor-level preservation.**
  [`writeLoom()`](https://mianaz.github.io/scConvert/reference/writeLoom.md)
  now stores factor levels and `ordered` flag under
  `/scConvert_extensions/col_factor_levels/{name}`, outside `col_attrs`
  so loompy/scanpy continue to read the file without errors.
  [`readLoom()`](https://mianaz.github.io/scConvert/reference/readLoom.md)
  restores the factors on load.
- **h5mu per-modality `uns` round-trip.**
  [`writeH5MU()`](https://mianaz.github.io/scConvert/reference/writeH5MU.md)
  and
  [`readH5MU()`](https://mianaz.github.io/scConvert/reference/readH5MU.md)
  now mirror each modality’s `uns` group under
  `obj@misc[["__h5mu_uns_per_mod__"]][[modality]]` so per-modality uns
  entries survive h5mu round-trips instead of being flattened.
- **IMC multi-image support.**
  [`SeuratSpatialToH5AD()`](https://mianaz.github.io/scConvert/reference/SeuratSpatialToH5AD.md)
  now iterates over every image in deterministic sorted order instead of
  processing only `Images()[1]`. Fixes the IMC 14/15 -\> 11/13
  double-roundtrip degradation documented in `dev/NOTES.md` section 3.

### Robustness

- **C CLI memory-safety helpers.** New `sc_xmalloc()`, `sc_xcalloc()`,
  `sc_xrealloc()`, and `sc_check_mul_size()` in `src/sc_util.c` replace
  raw `malloc()` calls with overflow-checked allocations at the
  dense-embedding transpose sites in `src/sc_zarr.c:1043, 1668, 1674`
  and the column-buffer allocation in `src/sc_loom.c:138`. Prevents
  SIZE_MAX overflow on chip-scale embeddings.

### Bug fixes

- **[`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md):
  handle unsorted CSR column indices.** scanpy-processed files such as
  `scanpy.datasets.pbmc3k_processed()` ship with CSR matrices whose
  column indices are not sorted within each row. Previously these
  produced an invalid `dgCMatrix` on read.
  [`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
  and
  [`readH5MU()`](https://mianaz.github.io/scConvert/reference/readH5MU.md)
  now detect this condition via `.sort_dgc_indices()` in `R/LoadH5AD.R`
  (and the sibling helper in `R/LoadH5MU.R`) and sort indices
  column-wise before constructing the sparse matrix. A regression test
  has been added in `tests/testthat/test-regression-fixes.R`.
- **[`H5SeuratToZarr()`](https://mianaz.github.io/scConvert/reference/H5SeuratToZarr.md):
  do not crash on 1D dense datasets.** The direct h5Seurat -\> Zarr
  converter previously attempted to infer (rows, cols) from a 1D dense
  HDF5 dataset and crashed with a dimensionality error. It now emits a
  warning and skips the offending dataset. Tracked as “Chain D” in the
  benchmark manuscript; regression test included.
- **[`writeZarr()`](https://mianaz.github.io/scConvert/reference/writeZarr.md):
  skip scale.data layer when its shape differs from X.**
  [`Seurat::ScaleData()`](https://satijalab.org/seurat/reference/ScaleData.html)
  produces an `(n_hvg x n_cells)` matrix whose row count does not match
  `X`’s `n_genes`. AnnData requires all `layers/*` to match `X`’s shape,
  so
  [`writeZarr()`](https://mianaz.github.io/scConvert/reference/writeZarr.md)
  now skips any layer whose dimensions differ from the default assay’s
  data matrix. A regression test covers this case. (See
  `R/SaveZarr.R:131`.)

### Testing

- Added `tests/testthat/test-regression-fixes.R` with three regression
  tests pinning the bugs listed above.
- Added `tests/testthat/test-generic-dispatch.R` (SOMA/SpatialData
  lambda signature pinning).
- Added `tests/testthat/test-stream-memory.R` (bounded-memory
  verification on tiny chunk budgets).
- Appended a canonical-h5mu-layout block to
  `tests/testthat/test-h5mu-multimodal.R`.
- Test suite is now **166 `test_that` blocks** with **773 assertions**
  on macOS / Ubuntu / Windows (was 137 / 448 at 0.1.0).
- 19 vignettes build cleanly under
  [`tools::buildVignettes()`](https://rdrr.io/r/tools/buildVignettes.html)
  including the 6 with live Python interop chunks via reticulate.

### CI

- **Python validation in CI** via `setup-micromamba` and
  `tests/scverse-env.yml`. anndata 0.12, scanpy 1.11, squidpy 1.8,
  mudata 0.3, loompy 3.0.

### Documentation

- Added `dev/NOTES.md` describing the state of recent bug fixes,
  benchmark findings that are specific to the manuscript (not to users),
  and investigations into dataset-pipeline issues (Stereo-seq MOSTA and
  CosMx SMI – upstream format is not HDF5; see the notes file for
  details).

## scConvert 0.1.0

> **Release Date:** 2026-03-10

### Highlights

Initial public release of scConvert — a universal single-cell format
converter for R.

#### Universal Format Conversion

- Support for 7 formats: h5ad, h5Seurat, h5mu, Loom, Zarr, RDS, and
  SingleCellExperiment
- Hub architecture with 30+ conversion paths via
  [`scConvert()`](https://mianaz.github.io/scConvert/reference/scConvert.md)
- Direct HDF5 paths for h5ad/h5Seurat without intermediate loading

#### Direct h5ad Loading

- [`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
  for native h5ad-to-Seurat conversion without intermediate files
- Sparse (CSR/CSC) and dense matrix support
- Categorical metadata, dimensional reductions, neighbor graphs, and
  spatial data

#### MuData (h5mu) Multimodal Support

- [`readH5MU()`](https://mianaz.github.io/scConvert/reference/readH5MU.md)
  /
  [`writeH5MU()`](https://mianaz.github.io/scConvert/reference/writeH5MU.md)
  for multimodal single-cell data
- Automatic modality-to-assay name mapping (rna-\>RNA, prot-\>ADT,
  atac-\>ATAC)
- No MuDataSeurat or Python dependency required

#### Zarr AnnData Support

- [`readZarr()`](https://mianaz.github.io/scConvert/reference/readZarr.md)
  and
  [`writeZarr()`](https://mianaz.github.io/scConvert/reference/writeZarr.md)
  for Zarr-based AnnData stores (v2 format)
- Sparse CSR/CSC and dense matrix support
- Categorical metadata and dimensional reduction preservation

#### Spatial Data (Visium)

- Bidirectional Visium spatial data conversion with image reconstruction
- Proper coordinate handling and scale factor preservation
- Compatible with scanpy/squidpy spatial analysis workflows

#### C CLI Binary

- Standalone `scconvert` binary for h5ad/h5Seurat/h5mu conversions
- Streaming on-disk conversion without R or Python runtime
- Options: `--assay`, `--gzip`, `--overwrite`, `--quiet`

#### BPCells On-Disk Loading

- `readH5AD(..., use.bpcells = TRUE)` for memory-efficient atlas-scale
  analysis
- Compatible with all Seurat analysis functions

#### Seurat v5 Compatibility

- Full support for Seurat v5 Assay5 objects
- Proper handling of V5 layered data structure in conversions
