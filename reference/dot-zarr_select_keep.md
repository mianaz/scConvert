# Apply a user-supplied selection filter to a vector of available names

- `NULL` returns `available` unchanged (keep all).

- `character(0)` returns `character(0)` (drop all).

- Otherwise returns the intersection.

## Usage

``` r
.zarr_select_keep(filter, available)
```
