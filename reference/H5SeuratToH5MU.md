# Convert an h5seurat file to h5mu format

Converts an h5seurat file to MuData h5mu format. By default uses
streaming conversion. Each Seurat assay becomes an h5mu modality.

## Usage

``` r
H5SeuratToH5MU(
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

  Path for output .h5mu file

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

- stream:

  If TRUE (default), stream fields directly. If FALSE, load as Seurat
  then save as h5mu.

## Value

Invisibly returns `dest`
