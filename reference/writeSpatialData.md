# Write a Seurat object to SpatialData zarr format

Writes a Seurat object to the scverse SpatialData zarr specification.
The expression data is written as an anndata table to `tables/table/`,
spatial coordinates are written to `shapes/`, and tissue images are
written to `images/` in OME-NGFF format.

## Usage

``` r
writeSpatialData(
  object,
  path,
  table = "table",
  region = NULL,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A
  [`Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
  object

- path:

  Path for the output .zarr directory

- table:

  Name for the table (default: "table")

- region:

  Name for the spatial region (default: auto-detect from image names, or
  "spots")

- overwrite:

  If TRUE, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `path`
