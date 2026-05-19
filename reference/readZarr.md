# Load an AnnData Zarr store as a Seurat object

Direct conversion from AnnData Zarr format (.zarr directory) to Seurat
object. Supports zarr v2 stores as produced by Python anndata. This
enables reading of cloud-native AnnData stores from CELLxGENE, Human
Cell Atlas, and SpatialData workflows.

## Usage

``` r
readZarr(
  file,
  assay.name = "RNA",
  verbose = TRUE,
  cache = TRUE,
  obs_idx = NULL,
  var_idx = NULL,
  layers = NULL,
  obsm = NULL,
  obsp = NULL,
  varm = NULL,
  varp = NULL,
  uns = NULL,
  include_x = TRUE,
  ...
)
```

## Arguments

- file:

  Path to a local .zarr directory, or a remote `s3://`/`gs://` URL.

- assay.name:

  Name for the primary assay (default: "RNA")

- verbose:

  Show progress messages

- cache:

  For remote URLs: `TRUE` to persist the download under
  `tools::R_user_dir("scConvert", "cache")`; `FALSE` to use a tempdir
  that is discarded with the R session. Ignored for local paths.

- obs_idx, var_idx:

  Integer indices (1-based) selecting cells / features to keep. `NULL`
  (default) reads all. Index slicing pushes down to chunk fetches: only
  chunks of the cell / feature index, obs columns, obsm embeddings, var
  columns, and sparse `indptr`-resolved data blocks containing the
  requested indices are pulled.

- layers:

  Character vector of `/layers/*` groups to read, or `NULL` (default)
  for all, or `character(0)` to skip all layers entirely. Skipping
  unused layers avoids their chunk fetches on remote stores.

- obsm, obsp, varm, varp, uns:

  Same semantics as `layers` for the corresponding AnnData groups:
  `NULL` keeps everything,
  [`character()`](https://rdrr.io/r/base/character.html) drops the
  entire group, otherwise reads only the named items.

- include_x:

  If `FALSE`, skip reading the main expression matrix entirely. Useful
  for metadata-only reads from a remote store. Default `TRUE`. (If
  `FALSE`, a zero-row placeholder assay is constructed so the returned
  Seurat object still has the right number of cells.)

- ...:

  Additional arguments (currently unused)

## Value

A `Seurat` object

## Details

`file` may also be a remote URL: `s3://bucket/key.zarr`,
`gs://bucket/key.zarr`. Anonymous (public) buckets only; private buckets
requiring SigV4 signing are not supported. When `cache` is `TRUE`, the
downloaded store is kept under `tools::R_user_dir("scConvert", "cache")`
and reused on subsequent calls with the same URL.
