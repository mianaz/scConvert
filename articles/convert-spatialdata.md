# Working with SpatialData (.spatialdata.zarr)

## Overview

[SpatialData](https://spatialdata.scverse.org/) is the scverse standard
for representing spatial omics data. It uses a Zarr-based directory
layout that combines OME-NGFF images, geometric shapes, transcript-level
point annotations, segmentation labels, and anndata expression tables
into a single hierarchical store. Technologies such as Visium, MERFISH,
Xenium, Slide-seq, CODEX, and Stereo-seq all have SpatialData
representations, making it a unifying format across the spatial omics
ecosystem.

scConvert provides native R support for reading and writing SpatialData
stores with no Python dependency, and includes direct pair converters
for moving data between SpatialData and other formats (h5ad, h5Seurat,
Zarr).

``` r

library(Seurat)
library(scConvert)
```

## SpatialData store layout

A `.spatialdata.zarr` directory follows this structure:

    sample.spatialdata.zarr/
    +-- .zattrs              # {"spatialdata_attrs": {"version": "0.2.0"}}
    +-- .zgroup              # {"zarr_format": 2}
    +-- images/              # OME-NGFF multiscale tissue images
    +-- labels/              # Segmentation masks (cell/nucleus outlines)
    +-- points/              # Point annotations (e.g., transcript coordinates)
    +-- shapes/              # Geometric shapes (spot positions with radii)
    +-- tables/              # AnnData expression tables (one or more)
        +-- table/           # Default table (standard anndata zarr)

The `tables/table/` subdirectory is a standard anndata zarr store, so
scConvert’s existing Zarr infrastructure handles the expression data,
metadata, embeddings, and graphs. The SpatialData-specific elements
(shapes, points, images, labels) are handled by dedicated readers and
writers.

## Reading SpatialData

Use
[`readSpatialData()`](https://mianaz.github.io/scConvert/reference/readSpatialData.md)
to load a SpatialData zarr store into a Seurat object:

``` r

obj <- readSpatialData("visium.spatialdata.zarr")
obj
```

The reader extracts the following components:

| SpatialData Element | Seurat Destination | Notes |
|----|----|----|
| `tables/table/` | Expression assay, metadata, reductions | Via [`readZarr()`](https://mianaz.github.io/scConvert/reference/readZarr.md) |
| `shapes/{region}/` | Spatial coordinates (FOV/image) | x, y, radius columns |
| `points/{region}/` | `misc$points_{name}` | Transcript-level data |
| `images/{name}/` | Image object (if coords present) | OME-NGFF, normalized to \[0,1\] |
| `labels/` | `misc$__spatialdata_labels__` | Names only (see Limitations) |

SpatialData attributes (`region`, `region_key`, `instance_key`) are
stored in `misc$__spatialdata_attrs__` and used during write-back to
ensure round-trip fidelity.

### Specifying a table

SpatialData stores can contain multiple tables. Use the `table` argument
to select which one to read:

``` r

# Read the default table
obj <- readSpatialData("multi_table.spatialdata.zarr", table = "table")

# Read a specific table
obj <- readSpatialData("multi_table.spatialdata.zarr", table = "rna_counts")
```

### Skipping images

For large stores where only the expression data and coordinates are
needed, set `images = FALSE` to skip OME-NGFF image reading:

``` r

obj <- readSpatialData("large_sample.spatialdata.zarr", images = FALSE)
```

## Writing SpatialData

Use
[`writeSpatialData()`](https://mianaz.github.io/scConvert/reference/writeSpatialData.md)
to write a Seurat object as a SpatialData zarr store:

``` r

writeSpatialData(obj, "output.spatialdata.zarr", overwrite = TRUE)
```

The writer creates the full SpatialData directory layout:

| Seurat Source | SpatialData Destination | Notes |
|----|----|----|
| Expression, metadata, reductions | `tables/table/` | Via [`writeZarr()`](https://mianaz.github.io/scConvert/reference/writeZarr.md) |
| Tissue coordinates | `shapes/{region}/` | DataFrame with x, y, radius |
| `misc$points_*` | `points/{name}/` | Preserved from read |
| Image arrays | `images/{name}/` | OME-NGFF with multiscales metadata |

The region name is auto-detected from the Seurat image name, or you can
set it explicitly:

``` r

writeSpatialData(obj, "output.spatialdata.zarr",
                 region = "spots", overwrite = TRUE)
```

If the object has metadata columns `spatial_x` and `spatial_y` but no
Seurat image object,
[`writeSpatialData()`](https://mianaz.github.io/scConvert/reference/writeSpatialData.md)
uses those as the shape coordinates.

## Direct pair converters

scConvert provides six direct conversion functions that combine reading
and writing in a single call, without requiring you to manage
intermediate Seurat objects:

``` r

# SpatialData <-> h5ad
SpatialDataToH5AD("sample.spatialdata.zarr", "sample.h5ad")
H5ADToSpatialData("sample.h5ad", "sample.spatialdata.zarr")

# SpatialData <-> h5Seurat
SpatialDataToH5Seurat("sample.spatialdata.zarr", "sample.h5seurat")
H5SeuratToSpatialData("sample.h5seurat", "sample.spatialdata.zarr")

# SpatialData <-> standard Zarr (anndata zarr without SpatialData wrapper)
SpatialDataToZarr("sample.spatialdata.zarr", "sample.zarr")
ZarrToSpatialData("sample.zarr", "sample.spatialdata.zarr")
```

All six functions accept `overwrite = TRUE` and `verbose = TRUE/FALSE`
arguments. The `SpatialDataToH5AD`, `SpatialDataToH5Seurat`, and
`SpatialDataToZarr` functions also accept a `table` argument to specify
which table to extract from the store.

These converters also work through the `scConvert()` dispatcher:

``` r

scConvert("sample.spatialdata.zarr", dest = "sample.h5ad", overwrite = TRUE)
scConvert("sample.h5ad", dest = "output.spatialdata.zarr", overwrite = TRUE)
```

## Python interoperability

SpatialData stores written by scConvert are readable by the Python
[spatialdata](https://spatialdata.scverse.org/) library and its
ecosystem (squidpy, napari-spatialdata):

``` python
import spatialdata as sd

sdata = sd.read_zarr("output.spatialdata.zarr")
print(sdata)
print(sdata.tables)    # Expression tables
print(sdata.shapes)    # Spot/cell geometries

# Visualize with spatialdata-plot
sdata.pl.render_shapes().pl.show()
```

Conversely, SpatialData stores created by the Python `spatialdata`
library (e.g., from `squidpy.datasets` or `spatialdata-io` readers) can
be loaded directly with
[`readSpatialData()`](https://mianaz.github.io/scConvert/reference/readSpatialData.md).

## Known limitations

- **Labels/segmentation masks**: scConvert records the names of label
  layers found in `labels/` but does not read the full mask arrays.
  These are typically large (cell/nucleus segmentation at pixel
  resolution) and have no natural representation in Seurat. Label names
  are preserved in `misc$__spatialdata_labels__` for metadata purposes.

- **Multiple tables**: Only one table is read per call. If your
  SpatialData store contains multiple tables (e.g., separate RNA and
  protein assays), call
  [`readSpatialData()`](https://mianaz.github.io/scConvert/reference/readSpatialData.md)
  once per table and merge the resulting Seurat objects.

- **Large images**: OME-NGFF images larger than 100 million pixels are
  skipped to avoid memory issues. Use `images = FALSE` and work with
  coordinates only, or downsample the image externally before reading.

- **Coordinate transforms**: SpatialData supports affine and sequence
  coordinate transformations between elements. These are not currently
  applied during reading; coordinates are used as-is from the zarr
  arrays.

## See also

- [Spatial Transcriptomics:
  Visium](https://mianaz.github.io/scConvert/articles/spatial-visium.md)
  – Visium-specific workflows with tissue images and squidpy validation
- [Spatial
  Technologies](https://mianaz.github.io/scConvert/articles/spatial-technologies.md)
  – support for MERFISH, Slide-seq, Xenium, CODEX, and other platforms
- [Conversions: Seurat and
  Zarr](https://mianaz.github.io/scConvert/articles/convert-zarr.md) –
  general Zarr format details, streaming conversion, and Python interop

## Session Info

``` r

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.3
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/Indiana/Indianapolis
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.56         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.30    lifecycle_1.0.5   cli_3.6.5         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4   systemfonts_1.3.1
#> [17] compiler_4.5.2    tools_4.5.2       ragg_1.5.0        bslib_0.10.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.1.7       fs_1.6.7          htmlwidgets_1.6.4
```
