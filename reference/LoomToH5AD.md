# Convert a loom file to h5ad format

Reads the loom into a Seurat object and writes it as h5ad.

## Usage

``` r
LoomToH5AD(
  source,
  dest,
  assay = "RNA",
  overwrite = FALSE,
  gzip = 4L,
  verbose = TRUE
)
```

## Arguments

- source:

  Path to input .loom file

- dest:

  Path for output .h5ad file

- assay:

  Assay name (default: "RNA")

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`
