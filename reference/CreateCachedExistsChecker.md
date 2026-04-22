# Create a cached HDF5 existence checker

Creates a closure that checks for HDF5 path existence with result
caching to improve performance for repeated checks.

## Usage

``` r
CreateCachedExistsChecker()
```

## Value

A function with signature `function(group, path)` that checks if `path`
exists in `group`, using a cache to avoid repeated HDF5 calls.

## Examples

``` r
if (FALSE) { # \dontrun{
safe_exists <- CreateCachedExistsChecker()
if (safe_exists(h5_group, "data/matrix")) {
  # path exists
}
} # }
```
