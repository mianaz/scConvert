# Stream a sparse matrix group from one HDF5 file to another

Copies data/indices/indptr arrays with optional CSR\<-\>CSC
reinterpretation. CSR of \[A x B\] has identical on-disk layout to CSC
of \[B x A\], so conversion between h5ad CSR (cells x genes) and
h5seurat CSC (genes x cells) requires only swapping the shape/dims
attribute.

## Usage

``` r
.stream_sparse_group(
  src_group,
  dst_parent,
  dst_name,
  gzip,
  src_format = "h5ad",
  dst_format = "h5seurat",
  assay = NULL
)
```

## Arguments

- src_group:

  H5Group containing data/indices/indptr

- dst_parent:

  H5Group to write into

- dst_name:

  Name for new group in dst_parent

- gzip:

  Integer gzip compression level

- src_format:

  Character, "h5ad" or "h5seurat" or "h5mu"

- dst_format:

  Character, "h5ad" or "h5seurat"

- assay:

  Character assay name (for h5seurat graphs)
