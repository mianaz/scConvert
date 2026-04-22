# scTranspose a matrix

scTranspose a matrix

## Usage

``` r
scTranspose(x, ...)

# S3 method for class 'dgCMatrix'
scTranspose(x, ...)

# S3 method for class 'H5D'
scTranspose(
  x,
  dest = GetParent(x = x),
  dname = paste0("t_", basename(path = x$get_obj_name())),
  overwrite = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'H5Group'
scTranspose(
  x,
  dest = GetParent(x = x),
  dname = paste0("t_", basename(path = x$get_obj_name())),
  overwrite = FALSE,
  ...
)
```

## Arguments

- x:

  A matrix to transpose

- ...:

  Arguments passed to other methods

- dest:

  ...

- dname:

  ...

- overwrite:

  ...

- verbose:

  Show progress updates

## Value

[`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
method: returns a
[`dgCMatrix`](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html) with
the data of `x` transposed

[`H5D`](http://hhoeflin.github.io/hdf5r/reference/H5D-class.md) and
[`H5Group`](http://hhoeflin.github.io/hdf5r/reference/H5Group-class.md)
methods: Invisibly returns `NULL`
