# Read HDF5 Files

Read data from HDF5 files

## Usage

``` r
# S3 method for class 'H5D'
as.array(x, ...)

# S3 method for class 'H5D'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'H5Group'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'H5D'
as.list(x, ...)

# S3 method for class 'H5Group'
as.list(x, recursive = TRUE, ...)

# S3 method for class 'H5D'
as.matrix(x, transpose = FALSE, ...)

# S3 method for class 'H5Group'
as.matrix(x, ...)
```

## Arguments

- x:

  An HDF5 dataset (H5D) object

- ...:

  Arguments passed to other methods

- row.names:

  `NULL` or a character vector giving the row names for the data frame.
  Missing values are not allowed.

- optional:

  logical. If `TRUE`, setting row names and converting column names (to
  syntactic names: see
  [`make.names`](https://rdrr.io/r/base/make.names.html)) is optional.
  Note that all of R's base package
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  use `optional` only for column names treatment, basically with the
  meaning of
  [`data.frame`](https://rdrr.io/r/base/data.frame.html)`(*, check.names = !optional)`.
  See also the `make.names` argument of the `matrix` method.

- recursive:

  Whether to recursively convert subgroups to nested lists

- transpose:

  scTranspose the matrix before returning

## Value

Varies depending on the method being called

`as.array`: an array with the data from the HDF5 dataset

`as.data.frame`: returns a
[`data.frame`](https://rdrr.io/r/base/data.frame.html) with the data
from the HDF5 dataset

`as.list`: a list with the data from the HDF5 dataset

`as.matrix`: a matrix with the data from the HDF5 dataset

## See also

[`as.sparse`](https://satijalab.github.io/seurat-object/reference/as.sparse.html)
