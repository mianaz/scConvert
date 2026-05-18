# Save a Seurat object as an AnnData Zarr store

Writes a Seurat object to AnnData Zarr v2 format (.zarr directory). The
output is readable by Python `anndata.read_zarr()`, scanpy, and any tool
supporting the AnnData on-disk specification.

## Usage

``` r
writeZarr(
  object,
  filename,
  assay = DefaultAssay(object),
  overwrite = FALSE,
  compressor = NULL,
  compression_level = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A Seurat object

- filename:

  Path for the output .zarr directory

- assay:

  Name of assay to export (default: `DefaultAssay(object)`)

- overwrite:

  If TRUE, overwrite existing zarr store

- compressor:

  Compression codec. `NULL` (default) auto-selects: Zstd when a Zstd
  codec package is on the search path, otherwise zlib. As of 2026 there
  is no maintained CRAN Zstd-bytes wrapper for R, so the default falls
  through to zlib for most users; install
  [zstdlite](https://github.com/coolbutuseless/zstdlite) from source to
  enable Zstd. Other accepted values: `"zstd"`, `"zlib"`/`"gzip"`,
  `"blosc"` (requires the blosc package), `"none"`, or a list
  `list(id = ..., level = ...)` matching the Zarr v2 compressor schema.
  Zstd at level 3 typically gives 1.5–2x smaller files than zlib level 6
  and decompresses 3–5x faster.

- compression_level:

  Override the per-codec default level. `NULL` uses zstd 3, zlib
  [`GetCompressionLevel()`](https://mianaz.github.io/scConvert/reference/GetCompressionLevel.md),
  blosc 5.

- verbose:

  Show progress messages

- ...:

  Additional arguments (currently unused)

## Value

Invisibly returns `filename`
