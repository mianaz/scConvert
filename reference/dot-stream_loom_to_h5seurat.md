# Stream loom to h5seurat (direct HDF5-to-HDF5)

Reads loom dense matrix, converts to sparse, and writes h5seurat layout.
Loom stores genes x cells (dense), h5seurat stores genes x cells (CSC
sparse).

## Usage

``` r
.stream_loom_to_h5seurat(
  loom_path,
  h5seurat_path,
  assay = "RNA",
  gzip = 4L,
  verbose = TRUE
)
```

## Arguments

- loom_path:

  Path to input loom file

- h5seurat_path:

  Path to output h5seurat file

- assay:

  Character assay name (default: "RNA")

- gzip:

  Integer gzip compression level

- verbose:

  Logical
