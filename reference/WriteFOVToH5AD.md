# Write FOV segmentation and molecules into an open h5ad library group

Writes `@boundaries` (as segmentation/) and `@molecules` (as molecules/)
into `lib_group`, which must be an already-open H5Group handle at
`/uns/spatial/{library_id}/`. Used by
[`SeuratSpatialToH5AD`](https://mianaz.github.io/scConvert/reference/SeuratSpatialToH5AD.md)
during Seurat -\> h5ad conversion.

## Usage

``` r
WriteFOVToH5AD(fov_obj, lib_group, library_id = "fov", verbose = TRUE)
```

## Arguments

- fov_obj:

  FOV object

- lib_group:

  hdf5r H5Group handle for the library

- library_id:

  Identifier string (for error messages)

- verbose:

  Print progress messages

## Value

invisible NULL
