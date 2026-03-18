# CITE-seq Analysis Workflow

## Introduction

CITE-seq (Cellular Indexing of Transcriptomes and Epitopes by
Sequencing) measures RNA transcription and surface protein abundance
from the same cells. This produces two assays per cell: RNA gene
expression and ADT (antibody-derived tag) protein counts. The h5mu
format is ideal for storing both modalities in a single file. scConvert
handles both assays during format conversion, and provides two export
paths:

- **h5ad**: Exports a single assay (RNA or ADT) per file, for scanpy
  workflows
- **h5mu**: Exports all assays in one file, for muon/MuData workflows

## Read CITE-seq data from h5mu

The shipped `citeseq_demo.h5mu` contains 500 cells with RNA (2,000
genes) and ADT (10 antibodies: CD3, CD4, CD8, CD45RA, CD56, CD16, CD11c,
CD14, CD19, CD34), pre-processed with PCA, UMAP, and clustering.

``` r

h5mu_file <- system.file("extdata", "citeseq_demo.h5mu", package = "scConvert")
obj <- readH5MU(h5mu_file)

cat("Cells:", ncol(obj), "\n")
#> Cells: 500
cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
#> Assays: ADT, RNA
cat("RNA features:", nrow(obj[["RNA"]]), "\n")
#> RNA features: 2000
cat("ADT features:", nrow(obj[["ADT"]]), "\n")
#> ADT features: 10
cat("ADT markers:", paste(rownames(obj[["ADT"]]), collapse = ", "), "\n")
#> ADT markers: CD3, CD4, CD8, CD45RA, CD56, CD16, CD11c, CD14, CD19, CD34
```

## Visualize RNA clusters

``` r

DimPlot(obj, group.by = "seurat_clusters", label = TRUE, pt.size = 0.8) +
  ggtitle("RNA-based clusters")
```

![](multimodal-citeseq_files/figure-html/dimplot-rna-1.png)

## Visualize protein expression

ADT markers reveal cell-surface protein levels that complement the
transcriptomic clusters. After reading h5mu, the ADT assay contains only
raw counts, so we normalize with CLR before plotting. CD3 marks T cells;
CD14 marks monocytes.

``` r

DefaultAssay(obj) <- "ADT"
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)
```

``` r

library(patchwork)

p1 <- FeaturePlot(obj, features = "CD3", pt.size = 0.8) + ggtitle("CD3 (T cells)")
p2 <- FeaturePlot(obj, features = "CD14", pt.size = 0.8) + ggtitle("CD14 (Monocytes)")
p1 + p2
```

![](multimodal-citeseq_files/figure-html/featureplot-adt-1.png)

### ADT expression across clusters

``` r

VlnPlot(obj, features = "CD3", group.by = "seurat_clusters", pt.size = 0.1) +
  ggtitle("CD3 protein across clusters") +
  NoLegend()
```

![](multimodal-citeseq_files/figure-html/vlnplot-adt-1.png)

## Export: h5ad vs h5mu

### h5ad – single assay

[`writeH5AD()`](https://mianaz.github.io/scConvert/reference/writeH5AD.md)
writes one assay at a time. By default it writes the active assay (RNA).
This is the right choice for standard scanpy workflows.

``` r

h5ad_path <- file.path(tempdir(), "citeseq_rna.h5ad")
writeH5AD(obj, h5ad_path, overwrite = TRUE)
cat("Wrote RNA to:", basename(h5ad_path), "\n")
#> Wrote RNA to: citeseq_rna.h5ad
cat("File size:", round(file.size(h5ad_path) / 1024^2, 1), "MB\n")
#> File size: 0.7 MB
```

### h5mu – all assays

[`writeH5MU()`](https://mianaz.github.io/scConvert/reference/writeH5MU.md)
writes every assay as a separate modality, keeping RNA and ADT together
in a single file.

``` r

h5mu_path <- file.path(tempdir(), "citeseq_roundtrip.h5mu")
writeH5MU(obj, h5mu_path, overwrite = TRUE)
cat("Wrote RNA + ADT to:", basename(h5mu_path), "\n")
#> Wrote RNA + ADT to: citeseq_roundtrip.h5mu
cat("File size:", round(file.size(h5mu_path) / 1024^2, 1), "MB\n")
#> File size: 0.7 MB
```

## Reload and verify

### From h5ad (RNA only)

``` r

loaded_h5ad <- readH5AD(h5ad_path)
cat("h5ad loaded:", ncol(loaded_h5ad), "cells x", nrow(loaded_h5ad), "features\n")
#> h5ad loaded: 500 cells x 2000 features
cat("Assays:", paste(Assays(loaded_h5ad), collapse = ", "), "\n")
#> Assays: RNA
```

### From h5mu (RNA + ADT)

``` r

loaded_h5mu <- readH5MU(h5mu_path)
cat("h5mu loaded:", ncol(loaded_h5mu), "cells\n")
#> h5mu loaded: 500 cells
cat("Assays:", paste(Assays(loaded_h5mu), collapse = ", "), "\n")
#> Assays: ADT, RNA
cat("RNA features:", nrow(loaded_h5mu[["RNA"]]), "\n")
#> RNA features: 2000
cat("ADT features:", nrow(loaded_h5mu[["ADT"]]), "\n")
#> ADT features: 10
```

### Verify expression preservation

``` r

# Compare against the original h5mu-loaded object
orig <- readH5MU(h5mu_file)
common_cells <- intersect(colnames(orig), colnames(loaded_h5mu))

common_genes <- intersect(rownames(orig[["RNA"]]), rownames(loaded_h5mu[["RNA"]]))
orig_rna <- as.numeric(GetAssayData(orig, assay = "RNA", layer = "counts")[
  head(common_genes, 100), head(common_cells, 100)])
rt_rna <- as.numeric(GetAssayData(loaded_h5mu, assay = "RNA", layer = "counts")[
  head(common_genes, 100), head(common_cells, 100)])
cat("RNA counts identical:", identical(orig_rna, rt_rna), "\n")
#> RNA counts identical: TRUE

common_adt <- intersect(rownames(orig[["ADT"]]), rownames(loaded_h5mu[["ADT"]]))
orig_adt <- as.numeric(GetAssayData(orig, assay = "ADT", layer = "counts")[
  common_adt, head(common_cells, 100)])
rt_adt <- as.numeric(GetAssayData(loaded_h5mu, assay = "ADT", layer = "counts")[
  common_adt, head(common_cells, 100)])
cat("ADT counts identical:", identical(orig_adt, rt_adt), "\n")
#> ADT counts identical: TRUE
```

## When to use h5ad vs h5mu

| Scenario | Recommended format | Why |
|----|----|----|
| scanpy RNA analysis | h5ad | scanpy expects single-AnnData input |
| Share with Python collaborator (multi-assay) | h5mu | Keeps RNA + ADT together in one file |
| Archive a CITE-seq experiment | h5mu | Preserves all modalities |
| Feed into cellxgene | h5ad | cellxgene reads h5ad, not h5mu |
| Downstream muon/MOFA+ analysis | h5mu | muon operates on MuData objects |

## Python interoperability

Both exported files are directly readable in Python.

``` python
# Requires Python: pip install scanpy
import scanpy as sc

adata = sc.read_h5ad("citeseq_rna.h5ad")
print(adata)
sc.pl.umap(adata, color="seurat_clusters")
```

``` python
# Requires Python: pip install mudata
import mudata as md

mdata = md.read_h5mu("citeseq.h5mu")
print(mdata)
for name, mod in mdata.mod.items():
    print(f"  {name}: {mod.n_obs} cells x {mod.n_vars} features")
```

## Clean up

``` r

unlink(c(h5ad_path, h5mu_path))
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
#>   [4] ggbeeswarm_0.7.3       spatstat.utils_3.2-2   farver_2.1.2          
#>   [7] rmarkdown_2.30         fs_1.6.7               ragg_1.5.0            
#>  [10] vctrs_0.7.1            ROCR_1.0-12            spatstat.explore_3.7-0
#>  [13] htmltools_0.5.9        sass_0.4.10            sctransform_0.4.3     
#>  [16] parallelly_1.46.1      KernSmooth_2.23-26     bslib_0.10.0          
#>  [19] htmlwidgets_1.6.4      desc_1.4.3             ica_1.0-3             
#>  [22] plyr_1.8.9             plotly_4.12.0          zoo_1.8-15            
#>  [25] cachem_1.1.0           igraph_2.2.2           mime_0.13             
#>  [28] lifecycle_1.0.5        pkgconfig_2.0.3        Matrix_1.7-4          
#>  [31] R6_2.6.1               fastmap_1.2.0          MatrixGenerics_1.22.0 
#>  [34] fitdistrplus_1.2-6     future_1.69.0          shiny_1.13.0          
#>  [37] digest_0.6.39          S4Vectors_0.48.0       tensor_1.5.1          
#>  [40] RSpectra_0.16-2        irlba_2.3.7            GenomicRanges_1.62.1  
#>  [43] textshaping_1.0.4      labeling_0.4.3         progressr_0.18.0      
#>  [46] spatstat.sparse_3.1-0  httr_1.4.8             polyclip_1.10-7       
#>  [49] abind_1.4-8            compiler_4.5.2         bit64_4.6.0-1         
#>  [52] withr_3.0.2            S7_0.2.1               fastDummies_1.7.5     
#>  [55] MASS_7.3-65            tools_4.5.2            vipor_0.4.7           
#>  [58] lmtest_0.9-40          otel_0.2.0             beeswarm_0.4.0        
#>  [61] httpuv_1.6.16          future.apply_1.20.2    goftest_1.2-3         
#>  [64] glue_1.8.0             nlme_3.1-168           promises_1.5.0        
#>  [67] grid_4.5.2             Rtsne_0.17             cluster_2.1.8.2       
#>  [70] reshape2_1.4.5         generics_0.1.4         hdf5r_1.3.12          
#>  [73] gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2           
#>  [76] data.table_1.18.2.1    XVector_0.50.0         BiocGenerics_0.56.0   
#>  [79] BPCells_0.2.0          spatstat.geom_3.7-0    RcppAnnoy_0.0.23      
#>  [82] ggrepel_0.9.7          RANN_2.6.2             pillar_1.11.1         
#>  [85] stringr_1.6.0          spam_2.11-3            RcppHNSW_0.6.0        
#>  [88] later_1.4.8            splines_4.5.2          dplyr_1.2.0           
#>  [91] lattice_0.22-9         survival_3.8-6         bit_4.6.0             
#>  [94] deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2          
#>  [97] pbapply_1.7-4          knitr_1.51             gridExtra_2.3         
#> [100] Seqinfo_1.0.0          IRanges_2.44.0         scattermore_1.2       
#> [103] stats4_4.5.2           xfun_0.56              matrixStats_1.5.0     
#> [106] UCSC.utils_1.6.1       stringi_1.8.7          lazyeval_0.2.2        
#> [109] yaml_2.3.12            evaluate_1.0.5         codetools_0.2-20      
#> [112] tibble_3.3.1           cli_3.6.5              uwot_0.2.4            
#> [115] xtable_1.8-8           reticulate_1.45.0      systemfonts_1.3.1     
#> [118] jquerylib_0.1.4        GenomeInfoDb_1.46.2    dichromat_2.0-0.1     
#> [121] Rcpp_1.1.1             globals_0.19.1         spatstat.random_3.4-4 
#> [124] png_0.1-8              ggrastr_1.0.2          spatstat.univar_3.1-6 
#> [127] parallel_4.5.2         pkgdown_2.2.0          dotCall64_1.2         
#> [130] listenv_0.10.1         viridisLite_0.4.3      scales_1.4.0          
#> [133] ggridges_0.5.7         purrr_1.2.1            crayon_1.5.3          
#> [136] rlang_1.1.7            cowplot_1.2.0
```
