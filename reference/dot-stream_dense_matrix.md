# Stream a dense matrix from one HDF5 file to another

Reads a dense dataset and writes it with gzip compression.

## Usage

``` r
.stream_dense_matrix(
  src_dataset,
  dst_parent,
  dst_name,
  gzip,
  transpose = FALSE
)
```

## Arguments

- src_dataset:

  H5D source dataset

- dst_parent:

  H5Group to write into

- dst_name:

  Name for new dataset

- gzip:

  Integer gzip compression level

- transpose:

  Logical, whether to transpose
