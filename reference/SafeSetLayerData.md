# Safely set layer data in a V5 Assay

Handles the case where LayerData\<- might not be exported from
SeuratObject

## Usage

``` r
SafeSetLayerData(object, layer, value)
```

## Arguments

- object:

  An Assay or Assay5 object

- layer:

  Layer name

- value:

  Data to set

## Value

Modified assay object
