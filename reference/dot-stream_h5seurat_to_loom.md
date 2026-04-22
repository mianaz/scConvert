# Stream h5seurat to loom (direct HDF5-to-HDF5)

Reads h5seurat sparse matrix (CSC), densifies, transposes to genes x
cells, and writes to loom format.

## Usage

``` r
.stream_h5seurat_to_loom(h5seurat_path, loom_path, gzip = 4L, verbose = TRUE)
```

## Arguments

- h5seurat_path:

  Path to input h5seurat file

- loom_path:

  Path to output loom file

- gzip:

  Integer gzip compression level

- verbose:

  Logical
