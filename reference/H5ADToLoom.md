# Convert an h5ad file to loom format

Converts via streaming: h5ad -\> h5seurat (temp) -\> loom.

## Usage

``` r
H5ADToLoom(
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

  Path to input .h5ad file

- dest:

  Path for output .loom file

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

- stream:

  If TRUE (default), use streaming. If FALSE, load into R.

## Value

Invisibly returns `dest`
