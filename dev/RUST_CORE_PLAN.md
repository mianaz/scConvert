# Rust Core Backend Plan

## Goal

Keep the current C CLI as the production native backend while Rust matures as a
parallel backend. Users should eventually be able to choose `c`, `rust`, `r`, or
`auto` based on route support, performance, and preference.

## Backend Contract

Both native backends must accept the same core options:

- `input`
- `output`
- `assay`
- `compression` / `gzip`
- `overwrite`
- `verbose`

Both backends must also share these semantics:

- Exit code `0` means success.
- Nonzero exit codes map to stable error classes.
- Successful output is validated before the caller returns success.
- Failure cleanup never deletes outside the requested output path.
- Format support is explicit per route; unsupported same-format pairs must not
  default to another conversion.

## User-Facing Backend Selection

R API:

```r
options(scConvert.native_backend = "auto") # auto | c | rust | r
options(scConvert.cli_path_c = "/path/to/scconvert-c")
options(scConvert.cli_path_rust = "/path/to/scconvert-rust")

scConvert_cli(
  input,
  output,
  backend = getOption("scConvert.native_backend", "auto")
)
```

CLI:

```bash
scconvert --backend auto in.h5ad out.h5seurat
scconvert --backend c in.h5ad out.h5seurat
scconvert --backend rust in.h5ad out.h5seurat
```

Recommended binary layout:

- `scconvert-c`: current C CLI.
- `scconvert-rust`: Rust CLI.
- `scconvert`: future dispatcher once both backends have stable route coverage.

## Backend Selection Policy

- `backend = "c"`: force C; error if unsupported.
- `backend = "rust"`: force Rust; error if unsupported. Experimental routes
  should require an explicit opt-in until promoted.
- `backend = "r"`: bypass native paths.
- `backend = "auto"`: choose the best stable backend for the route. Fallback is
  allowed only after a failed backend validates that no partial output survived.

Rust becomes the default for a route only after:

- Output loads in R and, where applicable, Python/scverse.
- Parity tests pass against C/R for the same fixtures.
- Benchmarks show equal or better performance or materially better diagnostics.
- Errors and fallback behavior are stable.

## Capability Registry

Add one route registry consumed by both R dispatch and CLI dispatch:

```r
native_backends <- data.frame(
  backend = c("c", "rust"),
  source = c("h5ad", "h5ad"),
  dest = c("h5seurat", "h5seurat"),
  status = c("stable", "experimental"),
  supports_dense = c(TRUE, FALSE),
  supports_zarr = c(TRUE, FALSE),
  supports_fidelity_report = c(FALSE, TRUE)
)
```

The registry should be data, not scattered `if` logic, so support can be
promoted route by route.

## Milestones

### B1: Stabilize Current Rust Slice

- Add required h5Seurat root attributes: `version`, `project`,
  `active.assay`.
- Add an R-load validation test for Rust-produced h5Seurat.
- Update Rust status documentation. The branch is no longer contract-only; it
  has a partial `h5ad -> h5seurat` implementation.
- Rename the Rust binary to avoid collision with the C CLI.

### B2: Complete `h5ad -> h5seurat`

- Sparse `X` and `raw/X`.
- Dense `X`.
- Common data/index dtypes: float32, float64, int32, int64.
- `obs -> meta.data`.
- `var -> features` and retained feature metadata.
- `obsm -> reductions`.
- `obsp -> graphs`.
- `varp`.
- `uns`.
- Required empty h5Seurat groups.
- Output passes `readH5Seurat()`.
- Output passes structural comparison against C CLI on fixtures.

### B3: Add `h5seurat -> h5ad`

- Establish the first Rust roundtrip pair.
- Validate:
  - `h5ad -> h5seurat -> h5ad`
  - `h5seurat -> h5ad -> h5seurat`
- Compare matrix values, cell names, feature names, metadata, reductions,
  graphs, and pairwise annotations.

### B4: Add h5mu and Loom After h5ad/h5seurat Are Stable

Suggested order:

- `h5mu -> h5seurat`
- `h5seurat -> h5mu`
- `loom -> h5seurat`
- `h5seurat -> loom`

C remains default for these routes until Rust parity is proven.

### B5: Zarr

Use a real Rust Zarr implementation, preferably `zarrs`, rather than expanding
hand-written Zarr parsing.

Candidate routes:

- `h5ad <-> zarr`
- `h5seurat <-> zarr`
- `zarr <-> h5mu`
- `zarr <-> loom`

### B6: Benchmark Harness

Compare:

- R streaming
- C CLI
- Rust CLI

Metrics:

- Wall time
- Peak RSS
- Output size
- Backend selected
- Validation result
- Fidelity report

Fixtures:

- Tiny synthetic sparse
- `pbmc_small`
- 100K sparse synthetic
- Dense small
- Multimodal h5mu
- Zarr fixture
- Malformed fixture

### B7: Public Backend Switching

Expose the stable backend controls after at least one Rust route is
production-ready. Keep startup quiet unless a forced backend is unavailable or
no native backend exists.

## Immediate Priority Steps

1. Patch Rust `h5ad -> h5seurat` to write the missing root attributes.
2. Add a Rust integration test that fails if the root attributes are missing.
3. Add an R validation test or script that calls `readH5Seurat()` on
   Rust-produced output.
4. Update Rust README status from B1 contract-only to B2 partial implementation.
5. Decide the Rust binary name before wiring it into R dispatch.
