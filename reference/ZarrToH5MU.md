# Convert a zarr store to h5mu format

Converts via streaming: zarr -\> h5seurat (temp) -\> h5mu.

## Usage

``` r
ZarrToH5MU(
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

  Path to input .zarr directory

- dest:

  Path for output .h5mu file

- assay:

  Assay name (default: "RNA")

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
