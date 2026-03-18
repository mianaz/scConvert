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

## Create a demo SpatialData store

We create a Seurat object with spatial coordinates and write it as a
SpatialData store using
[`writeSpatialData()`](https://mianaz.github.io/scConvert/reference/writeSpatialData.md):

``` r

library(Seurat)
set.seed(42)

# Synthetic spatial expression data: 100 genes × 200 cells
counts <- Matrix::rsparsematrix(100, 200, density = 0.05, rand.x = function(n) rpois(n, 5) + 1)
dimnames(counts) <- list(paste0("Gene", seq_len(100)), paste0("Cell", seq_len(200)))
demo <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)

# Add spatial coordinates and metadata
demo$spatial_x <- runif(200, 0, 1000)
demo$spatial_y <- runif(200, 0, 1000)
demo$cell_type <- factor(sample(c("Epithelial", "Stromal", "Immune"), 200, replace = TRUE))
demo$region <- "tissue_section"

# Write as SpatialData
writeSpatialData(demo, "demo.spatialdata.zarr", overwrite = TRUE, verbose = TRUE)
```

## Reading SpatialData

Use
[`readSpatialData()`](https://mianaz.github.io/scConvert/reference/readSpatialData.md)
to load the store into a Seurat object:

``` r

obj <- readSpatialData("demo.spatialdata.zarr", verbose = TRUE)
obj
#> An object of class Seurat 
#> 100 features across 200 samples within 1 assay 
#> Active assay: RNA (100 features, 0 variable features)
#>  1 layer present: counts
#>  1 spatial field of view present: spots
```

### What gets preserved

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

### Inspect the Seurat object

``` r

cat("Cells:    ", ncol(obj), "\n")
#> Cells:     200
cat("Features: ", nrow(obj), "\n")
#> Features:  100
cat("Layers:   ", paste(Layers(obj), collapse = ", "), "\n")
#> Layers:    counts
cat("Images:   ", paste(Images(obj), collapse = ", "), "\n")
#> Images:    spots

# SpatialData metadata preserved for roundtrip
sd_attrs <- obj@misc[["__spatialdata_attrs__"]]
if (!is.null(sd_attrs)) {
  cat("Region:       ", paste(sd_attrs$region, collapse = ", "), "\n")
  cat("Region key:   ", sd_attrs$region_key, "\n")
  cat("Instance key: ", sd_attrs$instance_key, "\n")
}
#> Region:        spots 
#> Region key:    region 
#> Instance key:  instance_id

# Points stored in misc
point_keys <- grep("^points_", names(obj@misc), value = TRUE)
if (length(point_keys) > 0) {
  for (pk in point_keys) {
    cat("Points (", pk, "): ", nrow(obj@misc[[pk]]), " entries\n", sep = "")
  }
}

# Label names
label_names <- obj@misc[["__spatialdata_labels__"]]
if (!is.null(label_names)) {
  cat("Labels:   ", paste(label_names, collapse = ", "), "\n")
}
```

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
to write the Seurat object back as a SpatialData zarr store, completing
the round trip:

``` r

writeSpatialData(obj, "demo_roundtrip.spatialdata.zarr", overwrite = TRUE)
cat("Written to: demo_roundtrip.spatialdata.zarr\n")
#> Written to: demo_roundtrip.spatialdata.zarr
cat("Store size: ",
    round(sum(file.info(
      list.files("demo_roundtrip.spatialdata.zarr",
                 recursive = TRUE, full.names = TRUE)
    )$size) / 1024^2, 1), " MB\n")
#> Store size:  0  MB
```

The writer creates the full SpatialData directory layout:

| Seurat Source | SpatialData Destination | Notes |
|----|----|----|
| Expression, metadata, reductions | `tables/table/` | Via [`writeZarr()`](https://mianaz.github.io/scConvert/reference/writeZarr.md) |
| Tissue coordinates | `shapes/{region}/` | DataFrame with x, y, radius |
| `misc$points_*` | `points/{name}/` | Preserved from read |
| Image arrays | `images/{name}/` | OME-NGFF with multiscales metadata |

### Verify the roundtrip

``` r

obj_rt <- readSpatialData("demo_roundtrip.spatialdata.zarr", verbose = FALSE)

cat("Original cells:   ", ncol(obj), "\n")
#> Original cells:    200
cat("Roundtrip cells:  ", ncol(obj_rt), "\n")
#> Roundtrip cells:   200
cat("Original features:", nrow(obj), "\n")
#> Original features: 100
cat("Roundtrip features:", nrow(obj_rt), "\n")
#> Roundtrip features: 100
cat("Dimensions match: ", ncol(obj) == ncol(obj_rt) && nrow(obj) == nrow(obj_rt), "\n")
#> Dimensions match:  TRUE
```

## Direct pair converters

scConvert provides six direct conversion functions that combine reading
and writing in a single call, without requiring you to manage
intermediate Seurat objects.

### SpatialData to h5ad

``` r

SpatialDataToH5AD("demo.spatialdata.zarr", "demo_spatial.h5ad", overwrite = TRUE)
cat("h5ad size:", round(file.size("demo_spatial.h5ad") / 1024^2, 1), "MB\n")
#> h5ad size: 0.1 MB
```

### All six converters

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
ecosystem (squidpy, napari-spatialdata).

### Read the h5ad output in Python

The h5ad produced by
[`SpatialDataToH5AD()`](https://mianaz.github.io/scConvert/reference/SpatialDataToH5AD.md)
is a standard AnnData file readable by scanpy:

``` python
import scanpy as sc

adata = sc.read_h5ad("demo_spatial.h5ad")
print(adata)
#> AnnData object with n_obs × n_vars = 200 × 100
#>     obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'spatial_x', 'spatial_y', 'cell_type', 'region'
#>     uns: '__spatialdata_version__', 'spatial'
#>     obsm: 'spatial'
print(f"\nobs columns: {list(adata.obs.columns)}")
#> 
#> obs columns: ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'spatial_x', 'spatial_y', 'cell_type', 'region']
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
#> Cells: 200, Genes: 100
```

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
#> other attached packages:
#> [1] scConvert_0.1.0    Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-1          
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.4        
#>   [4] spatstat.utils_3.2-2   farver_2.1.2           rmarkdown_2.30        
#>   [7] fs_1.6.7               ragg_1.5.0             vctrs_0.7.1           
#>  [10] ROCR_1.0-12            spatstat.explore_3.7-0 htmltools_0.5.9       
#>  [13] sass_0.4.10            sctransform_0.4.3      parallelly_1.46.1     
#>  [16] KernSmooth_2.23-26     bslib_0.10.0           htmlwidgets_1.6.4     
#>  [19] desc_1.4.3             ica_1.0-3              plyr_1.8.9            
#>  [22] plotly_4.12.0          zoo_1.8-15             cachem_1.1.0          
#>  [25] igraph_2.2.2           mime_0.13              lifecycle_1.0.5       
#>  [28] pkgconfig_2.0.3        Matrix_1.7-4           R6_2.6.1              
#>  [31] fastmap_1.2.0          MatrixGenerics_1.22.0  fitdistrplus_1.2-6    
#>  [34] future_1.69.0          shiny_1.13.0           digest_0.6.39         
#>  [37] S4Vectors_0.48.0       patchwork_1.3.2        tensor_1.5.1          
#>  [40] RSpectra_0.16-2        irlba_2.3.7            GenomicRanges_1.62.1  
#>  [43] textshaping_1.0.4      progressr_0.18.0       spatstat.sparse_3.1-0 
#>  [46] httr_1.4.8             polyclip_1.10-7        abind_1.4-8           
#>  [49] compiler_4.5.2         bit64_4.6.0-1          withr_3.0.2           
#>  [52] S7_0.2.1               fastDummies_1.7.5      MASS_7.3-65           
#>  [55] tools_4.5.2            lmtest_0.9-40          otel_0.2.0            
#>  [58] httpuv_1.6.16          future.apply_1.20.2    goftest_1.2-3         
#>  [61] glue_1.8.0             nlme_3.1-168           promises_1.5.0        
#>  [64] grid_4.5.2             Rtsne_0.17             cluster_2.1.8.2       
#>  [67] reshape2_1.4.5         generics_0.1.4         hdf5r_1.3.12          
#>  [70] gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2           
#>  [73] data.table_1.18.2.1    XVector_0.50.0         BiocGenerics_0.56.0   
#>  [76] BPCells_0.2.0          spatstat.geom_3.7-0    RcppAnnoy_0.0.23      
#>  [79] ggrepel_0.9.7          RANN_2.6.2             pillar_1.11.1         
#>  [82] stringr_1.6.0          spam_2.11-3            RcppHNSW_0.6.0        
#>  [85] later_1.4.8            splines_4.5.2          dplyr_1.2.0           
#>  [88] lattice_0.22-9         survival_3.8-6         bit_4.6.0             
#>  [91] deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2          
#>  [94] pbapply_1.7-4          knitr_1.51             gridExtra_2.3         
#>  [97] Seqinfo_1.0.0          IRanges_2.44.0         scattermore_1.2       
#> [100] stats4_4.5.2           xfun_0.56              matrixStats_1.5.0     
#> [103] UCSC.utils_1.6.1       stringi_1.8.7          lazyeval_0.2.2        
#> [106] yaml_2.3.12            evaluate_1.0.5         codetools_0.2-20      
#> [109] tibble_3.3.1           cli_3.6.5              uwot_0.2.4            
#> [112] xtable_1.8-8           reticulate_1.45.0      systemfonts_1.3.1     
#> [115] jquerylib_0.1.4        GenomeInfoDb_1.46.2    dichromat_2.0-0.1     
#> [118] Rcpp_1.1.1             globals_0.19.1         spatstat.random_3.4-4 
#> [121] png_0.1-8              spatstat.univar_3.1-6  parallel_4.5.2        
#> [124] pkgdown_2.2.0          ggplot2_4.0.2          dotCall64_1.2         
#> [127] listenv_0.10.1         viridisLite_0.4.3      scales_1.4.0          
#> [130] ggridges_0.5.7         purrr_1.2.1            crayon_1.5.3          
#> [133] rlang_1.1.7            cowplot_1.2.0
```
