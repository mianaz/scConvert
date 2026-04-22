# Convert Seurat spatial data to h5ad format

Convert Seurat spatial data to h5ad format

## Usage

``` r
SeuratSpatialToH5AD(
  seurat_obj,
  h5ad_file,
  library_id = "library_1",
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial data

- h5ad_file:

  H5AD file handle or path to write to

- library_id:

  Library ID for spatial data (default: "library_1")

- verbose:

  Print progress messages
