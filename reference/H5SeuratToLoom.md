# Convert an h5seurat file to loom format

Converts an h5seurat file to loom format. By default uses streaming
conversion. Sparse matrices are densified for loom output.

## Usage

``` r
H5SeuratToLoom(
  source,
  dest,
  overwrite = FALSE,
  gzip = 4L,
  verbose = TRUE,
  stream = TRUE
)
```

## Arguments

- source:

  Path to input .h5seurat file

- dest:

  Path for output .loom file

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

- stream:

  If TRUE (default), stream fields directly. If FALSE, load as Seurat
  then save as loom.

## Value

Invisibly returns `dest`
