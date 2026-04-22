# Convert a zarr store to h5ad format

Converts an AnnData zarr store to h5ad (HDF5) format. By default uses
streaming conversion that copies fields directly without loading into R
memory.

## Usage

``` r
ZarrToH5AD(
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

  Path to input .zarr directory

- dest:

  Path for output .h5ad file

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

- stream:

  If TRUE (default), stream fields directly without creating a Seurat
  intermediate. If FALSE, load as Seurat then save as h5ad.

## Value

Invisibly returns `dest`
