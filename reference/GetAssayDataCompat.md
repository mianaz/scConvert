# Get assay data with V4/V5 compatibility

Retrieves data from an Assay or Assay5 object, automatically using the
appropriate parameter name (layer for V5, slot for V4).

## Usage

``` r
GetAssayDataCompat(object, layer_or_slot, verbose = FALSE)
```

## Arguments

- object:

  Assay or Assay5 object

- layer_or_slot:

  Name of layer (V5) or slot (V4) to retrieve

- verbose:

  Print diagnostic messages

## Value

Matrix or sparse matrix of assay data

## See also

V5Compatibility.R for the original implementation
