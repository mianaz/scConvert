# Read a single column from a zarr DataFrame

Dispatches to numeric or string reading based on dtype.

## Usage

``` r
.sd_read_column(store_path, col_name)
```

## Arguments

- store_path:

  Path to the zarr DataFrame group

- col_name:

  Column name

## Value

Vector of values, or NULL if column not found
