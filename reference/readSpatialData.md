# Read a SpatialData zarr store into a Seurat object

Reads a SpatialData zarr store following the scverse SpatialData
specification. The expression data is read from the `tables/`
subdirectory (standard anndata zarr), spatial coordinates from `shapes/`
or `points/`, and tissue images from `images/` in OME-NGFF format.

## Usage

``` r
readSpatialData(path, table = "table", images = TRUE, verbose = TRUE)
```

## Arguments

- path:

  Path to the SpatialData .zarr directory

- table:

  Name of the table within `tables/` to read (default: "table")

- images:

  Logical; if TRUE (default), attempt to read OME-NGFF images from the
  `images/` directory

- verbose:

  Show progress messages

## Value

A
[`Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
object with spatial coordinates and optionally tissue images
