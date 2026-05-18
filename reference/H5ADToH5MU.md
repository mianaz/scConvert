# Convert an h5ad file to h5mu format

Reads the h5ad into a Seurat object and writes it as a single-modality
h5mu file.

## Usage

``` r
H5ADToH5MU(
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

## Value

Invisibly returns `dest`
