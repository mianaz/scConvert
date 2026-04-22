# Resolve a nested HDF5 path to its object

Resolve a nested HDF5 path to its object

## Usage

``` r
ResolveNestedH5Path(base_group, path)
```

## Arguments

- base_group:

  The starting H5Group

- path:

  A path string, possibly containing "/" for nested paths

## Value

The H5 object at the resolved path, or NULL if not found
