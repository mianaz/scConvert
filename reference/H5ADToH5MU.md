# Convert an h5ad file to h5mu format

Converts via streaming: h5ad -\> h5seurat (temp) -\> h5mu. Creates a
single-modality h5mu file.

## Usage

``` r
H5ADToH5MU(
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

  Path to input .h5ad file

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
