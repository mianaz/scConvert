# Set assay data with V4/V5 compatibility

Sets data in an Assay or Assay5 object, automatically using the
appropriate parameter name (layer for V5, slot for V4).

## Usage

``` r
SetAssayDataCompat(object, layer_or_slot, new.data, verbose = FALSE)
```

## Arguments

- object:

  Assay or Assay5 object

- layer_or_slot:

  Name of layer (V5) or slot (V4) to set

- new.data:

  New data to set

- verbose:

  Print diagnostic messages

## Value

Modified Assay or Assay5 object
