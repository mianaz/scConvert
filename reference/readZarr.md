# Load an AnnData Zarr store as a Seurat object

Direct conversion from AnnData Zarr format (.zarr directory) to Seurat
object. Supports zarr v2 stores as produced by Python anndata. This
enables reading of cloud-native AnnData stores from CELLxGENE, Human
Cell Atlas, and SpatialData workflows.

## Usage

``` r
readZarr(file, assay.name = "RNA", verbose = TRUE, cache = TRUE, ...)
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
