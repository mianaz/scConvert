# Convert a loom file to zarr format

Converts via streaming: loom -\> h5seurat (temp) -\> zarr.

## Usage

``` r
LoomToZarr(
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

  Path to input .loom file

- dest:

  Path for output .zarr directory

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
