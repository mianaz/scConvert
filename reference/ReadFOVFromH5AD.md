# Read FOV structure back from an h5ad library group

Inverse of
[`WriteFOVToH5AD`](https://mianaz.github.io/scConvert/reference/WriteFOVToH5AD.md):
reads segmentation/ and molecules/ subgroups from an already-open
`uns/spatial/{library}/` H5Group, and reconstitutes a `FOV` object.
Returns NULL when neither subgroup is present (non-FOV libraries such as
classic Visium).

## Usage

``` r
ReadFOVFromH5AD(lib_group, key = "fov", assay = "Spatial")
```

## Arguments

- lib_group:

  hdf5r H5Group handle for the library

- key:

  FOV key (trailing underscore added if missing)

- assay:

  Assay name to attach to the FOV (default "Spatial")

## Value

an FOV object or NULL
