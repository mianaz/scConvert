# Load a 10x Xenium output bundle into a Seurat object

Thin R-only wrapper around
[`Seurat::LoadXenium()`](https://satijalab.org/seurat/reference/ReadXenium.html)
that validates the bundle layout, sets the scConvert spatial-technology
tag, and returns a Seurat object ready for conversion to h5ad / h5Seurat
/ Zarr via
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md).

## Usage

``` r
LoadXenium(data.dir, fov = "fov", assay = "Xenium", verbose = TRUE)
```

## Arguments

- data.dir:

  Path to the Xenium output directory

- fov:

  FOV name attached to the resulting Seurat object (default "fov")

- assay:

  Seurat assay name (default "Xenium")

- verbose:

  Print progress messages

## Value

A Seurat object with an FOV slot containing cell centroids, segmentation
boundaries, and transcript molecules

## Details

The bundle directory must contain the canonical 10x Xenium output files
(`cell_feature_matrix.h5` or `cell_feature_matrix/`, `cells.csv.gz` or
`cells.parquet`, and optional `transcripts.parquet` /
`cell_boundaries.parquet` / `nucleus_boundaries.parquet`).

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- LoadXenium("/path/to/xenium_lung/")
scConvert(seurat_obj, "xenium_lung.h5ad")
} # }
```
