# Convert a loom file to h5seurat format

Converts a loom file to h5seurat format. By default uses streaming
conversion. The loom dense matrix is converted to sparse CSC for
h5seurat.

## Usage

``` r
LoomToH5Seurat(
  source,
  dest,
  assay = "RNA",
  overwrite = FALSE,
  gzip = 4L,
  verbose = TRUE,
  stream = TRUE
)
```

## Arguments

- source:

  Path to input .loom file

- dest:

  Path for output .h5seurat file

- assay:

  Assay name to use (default: "RNA")

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

- stream:

  If TRUE (default), stream fields directly. If FALSE, load as Seurat
  then save as h5seurat.

## Value

Invisibly returns `dest`
