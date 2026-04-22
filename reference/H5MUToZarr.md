# Convert an h5mu file to zarr format

Converts via streaming: h5mu -\> h5seurat (temp) -\> zarr.

## Usage

``` r
H5MUToZarr(
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

  Path to input .h5mu file

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
