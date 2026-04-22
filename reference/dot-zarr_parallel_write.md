# Write multiple zarr arrays in parallel

Uses mclapply on Unix for parallel compression, falls back to lapply on
Windows.

## Usage

``` r
.zarr_parallel_write(write_specs, compressor)
```

## Arguments

- write_specs:

  List of list(dir, data, dtype) specs

- compressor:

  Compressor spec
