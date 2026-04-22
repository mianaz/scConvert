# Load an AnnData Zarr store as a Seurat object

Direct conversion from AnnData Zarr format (.zarr directory) to Seurat
object. Supports zarr v2 stores as produced by Python anndata. This
enables reading of cloud-native AnnData stores from CELLxGENE, Human
Cell Atlas, and SpatialData workflows.

## Usage

``` r
readZarr(file, assay.name = "RNA", verbose = TRUE, ...)
```

## Arguments

- file:

  Path to .zarr directory

- assay.name:

  Name for the primary assay (default: "RNA")

- verbose:

  Show progress messages

- ...:

  Additional arguments (currently unused)

## Value

A `Seurat` object
