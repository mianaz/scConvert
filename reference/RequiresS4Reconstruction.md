# Check if assay requires S4 slot reconstruction

Determines if an assay object requires S4 slot reconstruction (V4 and
below) or should be left intact (V5).

## Usage

``` r
RequiresS4Reconstruction(object, h5_group = NULL, verbose = FALSE)
```

## Arguments

- object:

  Assay or Assay5 object

- h5_group:

  Optional HDF5 group to check for s4class attribute

- verbose:

  Print diagnostic messages

## Value

Logical indicating if S4 reconstruction is needed
