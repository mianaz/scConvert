# Convert an h5ad file to zarr format

Converts an AnnData h5ad file to zarr format. By default uses streaming
conversion that copies fields directly without loading into R memory.
Both formats follow the same AnnData on-disk specification.

## Usage

``` r
H5ADToZarr(
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

  Path to input .h5ad file

- dest:

  Path for output .zarr directory

- overwrite:

  If TRUE, overwrite existing output

- gzip:

  Gzip compression level (0-9)

- verbose:

  Show progress messages

- stream:

  If TRUE (default), stream fields directly without creating a Seurat
  intermediate. If FALSE, load as Seurat then save as zarr.

## Value

Invisibly returns `dest`
