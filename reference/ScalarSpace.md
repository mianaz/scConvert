# Get cached scalar HDF5 dataspace

Returns a singleton H5S scalar dataspace, avoiding repeated allocation.

## Usage

``` r
ScalarSpace()
```

## Value

An [`H5S`](http://hhoeflin.github.io/hdf5r/reference/H5S-class.md)
object of type "scalar"
