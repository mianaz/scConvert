# Copy HDF5 matrix data handling nested paths

Copy HDF5 matrix data handling nested paths

## Usage

``` r
CopyH5MatrixData(src_group, src_path, dst_loc, dst_name, verbose = FALSE)
```

## Arguments

- src_group:

  The source H5Group containing the data

- src_path:

  Path to the source data (can be nested)

- dst_loc:

  Destination H5File or H5Group

- dst_name:

  Name for the copied data at destination

- verbose:

  Show progress messages

## Value

Invisible NULL
