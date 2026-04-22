# Convert a zarr store to loom format

Converts via streaming: zarr -\> h5seurat (temp) -\> loom.

## Usage

``` r
ZarrToLoom(
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

  Path for output .loom file

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
