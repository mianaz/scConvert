# Convert an h5ad file to SpatialData zarr format

Reads an h5ad file and writes it as a SpatialData zarr store. If the
h5ad contains spatial information (`obsm/spatial`, Visium-style
`uns/spatial`), it is placed in the appropriate SpatialData elements.

## Usage

``` r
H5ADToSpatialData(source, dest, overwrite = FALSE, verbose = TRUE)
```

## Arguments

- source:

  Path to input .h5ad file

- dest:

  Path for output SpatialData .zarr directory

- overwrite:

  If TRUE, overwrite existing output

- verbose:

  Show progress messages

## Value

Invisibly returns `dest`
