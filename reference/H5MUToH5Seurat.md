# Convert an h5mu file to h5seurat format

Converts a MuData h5mu file to h5seurat format. By default uses
streaming conversion that copies fields directly between HDF5 files
without loading into R memory. Each h5mu modality becomes a Seurat
assay.

## Usage

``` r
H5MUToH5Seurat(
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

  Path for output .h5seurat file

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

- stream:

  If TRUE (default), stream fields directly. If FALSE, load as Seurat
  then save as h5seurat.

## Value

Invisibly returns `dest`
