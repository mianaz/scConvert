# Convert an h5mu file to h5ad format

Converts via streaming: h5mu -\> h5seurat (temp) -\> h5ad. Only the
first modality (or `assay`) is exported to h5ad.

## Usage

``` r
H5MUToH5AD(
  source,
  dest,
  assay = NULL,
  overwrite = FALSE,
  gzip = 4L,
  verbose = TRUE
)
```

## Arguments

- source:

  Path to input .h5mu file

- dest:

  Path for output .h5ad file

- assay:

  Assay name to export (default: first modality)

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`
