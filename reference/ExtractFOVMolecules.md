# Extract molecule coordinates from an FOV object

Returns a named list where each entry is a data frame of `(x, y, gene)`
transcript positions. The top-level name mirrors
`names(fov_obj@molecules)` (typically "molecules" for CosMx/Xenium).

## Usage

``` r
ExtractFOVMolecules(fov_obj)
```

## Arguments

- fov_obj:

  FOV object

## Value

named list or NULL
