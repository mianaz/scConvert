# Extract segmentation boundaries from an FOV object

Walks `fov_obj@boundaries`; returns a list keyed by boundary name. Each
element is a list with fields `type` ("Segmentation" or "Centroids"),
`cell_ids`, `coords` (n_vertices x 2 numeric) and `polygon_offsets`
(CSR-style row pointers, length n_cells + 1; present only for
Segmentation).

## Usage

``` r
ExtractFOVSegmentation(fov_obj)
```

## Arguments

- fov_obj:

  FOV object

## Value

named list or NULL
