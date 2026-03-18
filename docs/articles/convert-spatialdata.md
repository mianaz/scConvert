# SpatialData Zarr Stores

## Introduction

[SpatialData](https://spatialdata.scverse.org/) is the scverse standard
for spatial omics data. A `.spatialdata.zarr` store combines expression
tables, spatial coordinates, shapes, and tissue images in a single Zarr
directory. Technologies like Visium, MERFISH, Xenium, Slide-seq, and
CODEX all have SpatialData representations, making it a unifying format
for the spatial transcriptomics ecosystem.

scConvert provides
[`readSpatialData()`](https://mianaz.github.io/scConvert/reference/readSpatialData.md)
and
[`writeSpatialData()`](https://mianaz.github.io/scConvert/reference/writeSpatialData.md)
to move between SpatialData stores and Seurat objects, with no Python
dependency.

## Read spatial data from h5ad

We start by reading a 400-spot Visium mouse brain dataset shipped with
scConvert as an h5ad file. This dataset includes spatial coordinates, a
tissue image, PCA/UMAP embeddings, and 15 clusters across 1500 genes.

``` r

h5ad_path <- system.file("extdata", "spatial_demo.h5ad", package = "scConvert")
brain <- readH5AD(h5ad_path, verbose = FALSE)
brain
#> An object of class Seurat 
#> 1500 features across 400 samples within 1 assay 
#> Active assay: RNA (1500 features, 1500 variable features)
#>  2 layers present: counts, data
#>  2 dimensional reductions calculated: pca, umap
#>  1 spatial field of view present: anterior1
```

``` r

SpatialFeaturePlot(brain, features = "Ttr") +
  ggtitle("Ttr expression -- Visium mouse brain (from h5ad)")
```

![](convert-spatialdata_files/figure-html/spatial-plot-original-1.png)

Ttr (transthyretin) marks the choroid plexus. The spatial expression
pattern confirms the data was loaded correctly from the h5ad file.

``` r

DimPlot(brain, reduction = "umap", label = TRUE, pt.size = 1) +
  ggtitle("Cluster assignments (from h5ad)") + NoLegend()
```

![](convert-spatialdata_files/figure-html/dimplot-original-1.png)

## Write to SpatialData and read back

SpatialData stores have a complex directory structure with separate
groups for tables, shapes, images, and coordinate systems. We write the
Seurat object to a SpatialData Zarr store and read it back to verify the
round-trip.

``` r

sd_path <- file.path(tempdir(), "brain.spatialdata.zarr")
writeSpatialData(brain, sd_path, overwrite = TRUE, verbose = FALSE)
cat("SpatialData store written to:", sd_path, "\n")
#> SpatialData store written to: /var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T//RtmpD6izHU/brain.spatialdata.zarr
```

``` r

brain_rt <- readSpatialData(sd_path, verbose = FALSE)
cat("Cells:", ncol(brain_rt), "| Genes:", nrow(brain_rt), "\n")
#> Cells: 400 | Genes: 1500
cat("Reductions:", paste(Reductions(brain_rt), collapse = ", "), "\n")
#> Reductions: pca, umap
cat("Images:", paste(Images(brain_rt), collapse = ", "), "\n")
#> Images: anterior1
```

## Compare original and round-trip

``` r

library(patchwork)

p1 <- FeaturePlot(brain, features = "Ttr", pt.size = 1) +
  ggtitle("Ttr -- Original (h5ad)")
p2 <- FeaturePlot(brain_rt, features = "Ttr", pt.size = 1) +
  ggtitle("Ttr -- After SpatialData round-trip")
p1 + p2
```

![](convert-spatialdata_files/figure-html/compare-feature-1.png)

### Fidelity check

``` r

stopifnot(ncol(brain_rt) == ncol(brain))
stopifnot(nrow(brain_rt) == nrow(brain))
stopifnot(identical(sort(colnames(brain_rt)), sort(colnames(brain))))
cat("All checks passed:", ncol(brain_rt), "spots x", nrow(brain_rt), "genes.\n")
#> All checks passed: 400 spots x 1500 genes.
```

## Direct pair converters

Convert between SpatialData and other formats in a single call:

``` r

# SpatialData <-> h5ad
SpatialDataToH5AD("sample.spatialdata.zarr", "sample.h5ad")
H5ADToSpatialData("sample.h5ad", "sample.spatialdata.zarr")

# SpatialData <-> h5Seurat
SpatialDataToH5Seurat("sample.spatialdata.zarr", "sample.h5seurat")
H5SeuratToSpatialData("sample.h5seurat", "sample.spatialdata.zarr")

# Or use the universal dispatcher
scConvert("sample.spatialdata.zarr", dest = "sample.h5ad", overwrite = TRUE)
```

## Limitations

- **Images**: Large OME-NGFF images (\>100M pixels) are skipped to avoid
  memory issues. Use `images = FALSE` to read coordinates only.
- **Coordinate transforms**: SpatialData supports affine transformations
  between elements. These are not currently applied during reading;
  coordinates are used as-is.
- **Labels/segmentation**: Cell and nucleus segmentation masks are
  recorded by name but not loaded (they have no natural Seurat
  representation).
- **Multiple tables**: Only one table is read per call. Use the `table`
  argument to select which one.

## Python interoperability

SpatialData stores written by scConvert are readable by the Python
[spatialdata](https://spatialdata.scverse.org/) library. Requires Python
with spatialdata installed.

``` python
import spatialdata as sd

sdata = sd.read_zarr("brain.spatialdata.zarr")
print(sdata)
sdata.pl.render_shapes().pl.show()
```

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
#> [1] patchwork_1.3.2    ggplot2_4.0.2      Seurat_5.4.0       SeuratObject_5.3.0
#> [5] sp_2.2-1           scConvert_0.1.0   
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
#>  [31] fastmap_1.2.0          fitdistrplus_1.2-6     future_1.69.0         
#>  [34] shiny_1.13.0           digest_0.6.39          tensor_1.5.1          
#>  [37] RSpectra_0.16-2        irlba_2.3.7            textshaping_1.0.4     
#>  [40] labeling_0.4.3         progressr_0.18.0       spatstat.sparse_3.1-0 
#>  [43] httr_1.4.8             polyclip_1.10-7        abind_1.4-8           
#>  [46] compiler_4.5.2         bit64_4.6.0-1          withr_3.0.2           
#>  [49] S7_0.2.1               fastDummies_1.7.5      MASS_7.3-65           
#>  [52] tools_4.5.2            lmtest_0.9-40          otel_0.2.0            
#>  [55] httpuv_1.6.16          future.apply_1.20.2    goftest_1.2-3         
#>  [58] glue_1.8.0             nlme_3.1-168           promises_1.5.0        
#>  [61] grid_4.5.2             Rtsne_0.17             cluster_2.1.8.2       
#>  [64] reshape2_1.4.5         generics_0.1.4         hdf5r_1.3.12          
#>  [67] gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2           
#>  [70] data.table_1.18.2.1    spatstat.geom_3.7-0    RcppAnnoy_0.0.23      
#>  [73] ggrepel_0.9.7          RANN_2.6.2             pillar_1.11.1         
#>  [76] stringr_1.6.0          spam_2.11-3            RcppHNSW_0.6.0        
#>  [79] later_1.4.8            splines_4.5.2          dplyr_1.2.0           
#>  [82] lattice_0.22-9         survival_3.8-6         bit_4.6.0             
#>  [85] deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2          
#>  [88] pbapply_1.7-4          knitr_1.51             gridExtra_2.3         
#>  [91] scattermore_1.2        xfun_0.56              matrixStats_1.5.0     
#>  [94] stringi_1.8.7          lazyeval_0.2.2         yaml_2.3.12           
#>  [97] evaluate_1.0.5         codetools_0.2-20       tibble_3.3.1          
#> [100] cli_3.6.5              uwot_0.2.4             xtable_1.8-8          
#> [103] reticulate_1.45.0      systemfonts_1.3.1      jquerylib_0.1.4       
#> [106] dichromat_2.0-0.1      Rcpp_1.1.1             globals_0.19.1        
#> [109] spatstat.random_3.4-4  png_0.1-8              spatstat.univar_3.1-6 
#> [112] parallel_4.5.2         pkgdown_2.2.0          dotCall64_1.2         
#> [115] listenv_0.10.1         viridisLite_0.4.3      scales_1.4.0          
#> [118] ggridges_0.5.7         purrr_1.2.1            crayon_1.5.3          
#> [121] rlang_1.1.7            cowplot_1.2.0
```
