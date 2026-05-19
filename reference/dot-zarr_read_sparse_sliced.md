# Sparse-matrix slice via indptr pushdown

AnnData stores sparse X as rows = cells / cols = features in CSR, or the
reverse in CSC. For CSR we use indptr to locate the data/indices ranges
for the requested rows (cells = obs_idx) and only fetch those chunks of
`data`/`indices`; columns (var_idx) are then filtered in memory. For CSC
the same trick works on var_idx; obs_idx falls back to a memory subset.

## Usage

``` r
.zarr_read_sparse_sliced(
  store_path,
  rel_path,
  encoding_type,
  shape,
  transpose,
  obs_idx,
  var_idx
)
```
