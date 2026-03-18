# Converting Between Seurat and TileDB-SOMA

## What is TileDB-SOMA?

[TileDB-SOMA](https://github.com/single-cell-data/TileDB-SOMA) is a
cloud-native, array-based format for storing and querying single-cell
data at scale. It is the storage backend behind [CELLxGENE
Census](https://chanzuckerberg.github.io/cellxgene-census/), the largest
public atlas of single-cell data (61M+ cells across 900+ datasets).

Key properties of SOMA:

- **Columnar storage** backed by TileDB arrays – efficient slicing by
  cells or features without loading the full dataset
- **Cloud-native**: local paths, S3, GCS, and TileDB Cloud URIs
- **Query filtering**: select cells by metadata predicates at the
  storage level
- **Multi-measurement**: RNA, ATAC, and protein assays in one experiment

scConvert provides
[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md)
and
[`writeSOMA()`](https://mianaz.github.io/scConvert/reference/writeSOMA.md)
for Seurat interoperability, plus pair converters for all supported
formats. All SOMA functions require the
[tiledbsoma](https://github.com/single-cell-data/TileDB-SOMA) R package.

``` r

library(scConvert)
library(Seurat)
```

## Writing a Seurat object to SOMA

We start with the shipped `pbmc_small` dataset, write it to a SOMA
experiment, then read it back to verify the roundtrip.

``` r

# Load the shipped demo data
pbmc <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
pbmc
#> An object of class Seurat 
#> 2000 features across 214 samples within 1 assay 
#> Active assay: RNA (2000 features, 2000 variable features)
#>  3 layers present: counts, data, scale.data
#>  2 dimensional reductions calculated: pca, umap
```

Use
[`writeSOMA()`](https://mianaz.github.io/scConvert/reference/writeSOMA.md)
to save any Seurat object as a SOMA experiment:

``` r

soma_uri <- file.path(tempdir(), "pbmc_small.soma")

# Ensure graphs have a valid DefaultAssay (required by tiledbsoma)
for (gn in names(pbmc@graphs)) {
  if (length(DefaultAssay(pbmc@graphs[[gn]])) == 0) {
    DefaultAssay(pbmc@graphs[[gn]]) <- DefaultAssay(pbmc)
  }
}

writeSOMA(pbmc, uri = soma_uri, overwrite = TRUE)
cat("SOMA experiment written to:", soma_uri, "\n")
#> SOMA experiment written to: /var/folders/9l/bl67cpdj3rzgkx2pfk0flmhc0000gn/T//RtmpCi1jNj/pbmc_small.soma
```

[`writeSOMA()`](https://mianaz.github.io/scConvert/reference/writeSOMA.md)
delegates to
[`tiledbsoma::write_soma()`](https://single-cell-data.github.io/TileDB-SOMA/reference/write_soma.html),
which converts all assays, metadata, embeddings, and graphs. Each Seurat
assay becomes a separate SOMA measurement.

## Reading from SOMA

Use
[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md)
to load a SOMA experiment as a Seurat object:

``` r

obj <- readSOMA(soma_uri)
obj
#> An object of class Seurat 
#> 2000 features across 214 samples within 1 assay 
#> Active assay: RNA (2000 features, 0 variable features)
#>  1 layer present: counts
#>  2 dimensional reductions calculated: pca, umap
```

### Verify roundtrip fidelity

``` r

# Cell and gene counts match
cat("Original cells:", ncol(pbmc), " Roundtrip cells:", ncol(obj), "\n")
#> Original cells: 214  Roundtrip cells: 214
cat("Original genes:", nrow(pbmc), " Roundtrip genes:", nrow(obj), "\n")
#> Original genes: 2000  Roundtrip genes: 2000
stopifnot(ncol(obj) == ncol(pbmc))
stopifnot(nrow(obj) == nrow(pbmc))

# Metadata columns are preserved
orig_meta_cols <- sort(names(pbmc@meta.data))
rt_meta_cols <- sort(names(obj@meta.data))
shared_cols <- intersect(orig_meta_cols, rt_meta_cols)
cat("Shared metadata columns:", length(shared_cols), "/", length(orig_meta_cols), "\n")
#> Shared metadata columns: 7 / 7
cat("Columns:", paste(shared_cols, collapse = ", "), "\n")
#> Columns: nCount_RNA, nFeature_RNA, orig.ident, percent.mt, RNA_snn_res.0.5, seurat_annotations, seurat_clusters
```

### Filtering cells and features

[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md)
supports server-side filtering via value filter queries, which is
particularly useful for large atlases:

``` r

# Filter by cell type and tissue
obj <- readSOMA(
  "path/to/experiment",
  obs_query = "cell_type == 'T cell' & tissue == 'blood' & disease == 'normal'",
  obs_column_names = c("cell_type", "tissue", "donor_id")
)

# Read a specific measurement (default is "RNA")
obj <- readSOMA("path/to/experiment", measurement = "ATAC")
```

### What gets read

[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md)
maps the following SOMA components to Seurat:

| SOMA Source | Seurat Destination | Description |
|----|----|----|
| `obs` | `meta.data` | Cell metadata (factors, numerics, strings) |
| `ms[measurement]/var` | Feature metadata | Gene-level annotations |
| `ms[measurement]/X/data` | Default assay counts | Expression matrix (sparse) |
| `ms[measurement]/obsm/*` | `reductions` | PCA, UMAP, etc. |
| `ms[measurement]/obsp/*` | `graphs` | Neighbor graphs (SNN, distances) |
| Other measurements | Additional assays | Multi-modal support |

The `soma_joinid` column (internal to TileDB-SOMA) is automatically
removed from cell and feature metadata.

## Direct pair converters

scConvert provides direct conversion functions between SOMA and all
other supported formats. All route through a Seurat intermediate.

### SOMA and h5ad

``` r

h5ad_path <- file.path(tempdir(), "demo.h5ad")
soma_from_h5ad <- file.path(tempdir(), "from_h5ad.soma")
h5ad_from_soma <- file.path(tempdir(), "from_soma.h5ad")

# Write pbmc_small to h5ad first
writeH5AD(pbmc, filename = h5ad_path, overwrite = TRUE)

# h5ad -> SOMA
H5ADToSOMA(h5ad_path, soma_from_h5ad, overwrite = TRUE)

# SOMA -> h5ad
SOMAToH5AD(soma_from_h5ad, h5ad_from_soma, overwrite = TRUE)
```

Verify the h5ad roundtrip:

``` r

rt <- readH5AD(h5ad_from_soma)
cat("h5ad -> SOMA -> h5ad roundtrip:\n")
#> h5ad -> SOMA -> h5ad roundtrip:
cat("  Cells:", ncol(rt), " Genes:", nrow(rt), "\n")
#>   Cells: 214  Genes: 2000
stopifnot(ncol(rt) == ncol(pbmc))
stopifnot(nrow(rt) == nrow(pbmc))
```

### SOMA and h5Seurat

``` r

h5seurat_path <- file.path(tempdir(), "demo.h5seurat")
soma_from_h5s <- file.path(tempdir(), "from_h5seurat.soma")
h5s_from_soma <- file.path(tempdir(), "from_soma.h5seurat")

# Write to h5Seurat, then convert through SOMA
writeH5Seurat(pbmc, filename = h5seurat_path, overwrite = TRUE)
H5SeuratToSOMA(h5seurat_path, soma_from_h5s, overwrite = TRUE)
SOMAToH5Seurat(soma_from_h5s, h5s_from_soma, overwrite = TRUE)

rt2 <- readH5Seurat(h5s_from_soma)
cat("h5Seurat -> SOMA -> h5Seurat roundtrip:\n")
#> h5Seurat -> SOMA -> h5Seurat roundtrip:
cat("  Cells:", ncol(rt2), " Genes:", nrow(rt2), "\n")
#>   Cells: 214  Genes: 2000
stopifnot(ncol(rt2) == ncol(pbmc))
```

### Other format pairs

Converters also exist for Loom, Zarr, and h5mu – each follows the same
`{Format}ToSOMA()` / `SOMATo{Format}()` pattern:

``` r

LoomToSOMA("data.loom", "data.soma")    # Loom -> SOMA
ZarrToSOMA("data.zarr", "data.soma")    # Zarr -> SOMA
H5MUToSOMA("data.h5mu", "data.soma")    # h5mu -> SOMA
SOMAToLoom("data.soma", "output.loom")   # SOMA -> Loom
SOMAToZarr("data.soma", "output.zarr")   # SOMA -> Zarr
SOMAToH5MU("data.soma", "output.h5mu")  # SOMA -> h5mu
```

## CLI and dispatcher usage

The C CLI binary handles HDF5-to-HDF5 conversions only, so SOMA is not
supported in the compiled binary. However, the R-level
[`scConvert_cli()`](https://mianaz.github.io/scConvert/reference/scConvert_cli.md)
and the universal `scConvert()` dispatcher auto-detect SOMA from the
`.soma` extension and route through the R hub:

``` r

# scConvert_cli auto-dispatches
scConvert_cli("data.h5ad", "data.soma")
scConvert_cli("data.soma", "output.h5ad")

# Universal dispatcher
scConvert("data.h5ad", dest = "data.soma", overwrite = TRUE)
scConvert("data.soma", dest = "output.h5seurat", overwrite = TRUE)
```

## CELLxGENE Census access

[CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/)
hosts 61M+ cells as a single SOMA experiment. Query it with
[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md),
then convert the result to any format:

``` r

library(cellxgene.census)

census <- open_soma(census_version = "stable")
human_uri <- census$get("census_data")$get("homo_sapiens")$uri

# Read a subset directly as a Seurat object
tcells <- readSOMA(
  uri = human_uri,
  measurement = "RNA",
  obs_query = "cell_type == 'T cell' & tissue_general == 'blood'",
  obs_column_names = c("cell_type", "tissue", "donor_id")
)

# Save locally in any format
writeH5AD(tcells, "census_tcells.h5ad")
writeSOMA(tcells, "census_tcells.soma")
```

See the [cellxgene-census
documentation](https://chanzuckerberg.github.io/cellxgene-census/) for
more details on the Census API.

## Data mapping reference

| Seurat | SOMA | Direction | Notes |
|----|----|----|----|
| Default assay counts | `ms[RNA]/X/data` | Read/Write | Sparse matrix |
| `meta.data` | `obs` | Read/Write | Factors, numerics, strings |
| Feature metadata | `ms[RNA]/var` | Read/Write | Gene annotations |
| `Embeddings(, "pca")` | `ms[RNA]/obsm/X_pca` | Read/Write | Dense array |
| `Embeddings(, "umap")` | `ms[RNA]/obsm/X_umap` | Read/Write | Dense array |
| `Graphs(, "RNA_snn")` | `ms[RNA]/obsp/RNA_snn` | Read/Write | Sparse matrix |
| Additional assays | Other measurements | Read/Write | Multi-modal support |

## Limitations

- **Requires tiledbsoma**: install with
  `install.packages("tiledbsoma")`.
  [`writeSOMA()`](https://mianaz.github.io/scConvert/reference/writeSOMA.md)
  further requires
  [`tiledbsoma::write_soma()`](https://single-cell-data.github.io/TileDB-SOMA/reference/write_soma.html)
  to be available.
- **No streaming conversion**: unlike h5ad/zarr converters, SOMA always
  routes through a Seurat intermediate (TileDB arrays differ
  fundamentally from HDF5).
- **No C binary path**: the compiled C CLI handles HDF5-to-HDF5 only.
- **scale.data** is not written to SOMA – re-run
  [`ScaleData()`](https://satijalab.org/seurat/reference/ScaleData.html)
  after loading.
- **Active identity**
  ([`Idents()`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md))
  is not stored – set manually after reading.

## See also

- [`vignette("convert-anndata")`](https://mianaz.github.io/scConvert/articles/convert-anndata.md)
  – Converting between Seurat and h5ad format
- [`vignette("cli-usage")`](https://mianaz.github.io/scConvert/articles/cli-usage.md)
  – Command-line interface for batch conversions
- [`vignette("convert-zarr")`](https://mianaz.github.io/scConvert/articles/convert-zarr.md)
  – Converting between Seurat and Zarr format

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
#> [1] Seurat_5.4.0       SeuratObject_5.3.0 sp_2.2-1           scConvert_0.1.0   
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
#>  [37] tiledb_0.34.0          S4Vectors_0.48.0       patchwork_1.3.2       
#>  [40] tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.7           
#>  [43] GenomicRanges_1.62.1   textshaping_1.0.4      progressr_0.18.0      
#>  [46] spatstat.sparse_3.1-0  httr_1.4.8             polyclip_1.10-7       
#>  [49] abind_1.4-8            compiler_4.5.2         bit64_4.6.0-1         
#>  [52] S7_0.2.1               fastDummies_1.7.5      MASS_7.3-65           
#>  [55] tiledbsoma_2.3.0       tools_4.5.2            lmtest_0.9-40         
#>  [58] otel_0.2.0             httpuv_1.6.16          future.apply_1.20.2   
#>  [61] goftest_1.2-3          glue_1.8.0             nlme_3.1-168          
#>  [64] promises_1.5.0         grid_4.5.2             Rtsne_0.17            
#>  [67] cluster_2.1.8.2        reshape2_1.4.5         generics_0.1.4        
#>  [70] hdf5r_1.3.12           gtable_0.3.6           spatstat.data_3.1-9   
#>  [73] tidyr_1.3.2            data.table_1.18.2.1    XVector_0.50.0        
#>  [76] BiocGenerics_0.56.0    BPCells_0.2.0          spatstat.geom_3.7-0   
#>  [79] RcppAnnoy_0.0.23       ggrepel_0.9.7          RANN_2.6.2            
#>  [82] pillar_1.11.1          stringr_1.6.0          spam_2.11-3           
#>  [85] nanoarrow_0.8.0        RcppHNSW_0.6.0         later_1.4.8           
#>  [88] splines_4.5.2          dplyr_1.2.0            lattice_0.22-9        
#>  [91] survival_3.8-6         bit_4.6.0              deldir_2.0-4          
#>  [94] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-4         
#>  [97] knitr_1.51             gridExtra_2.3          Seqinfo_1.0.0         
#> [100] IRanges_2.44.0         RcppCCTZ_0.2.14        scattermore_1.2       
#> [103] stats4_4.5.2           xfun_0.56              matrixStats_1.5.0     
#> [106] UCSC.utils_1.6.1       stringi_1.8.7          lazyeval_0.2.2        
#> [109] yaml_2.3.12            evaluate_1.0.5         codetools_0.2-20      
#> [112] tibble_3.3.1           cli_3.6.5              uwot_0.2.4            
#> [115] arrow_23.0.1.1         xtable_1.8-8           reticulate_1.45.0     
#> [118] systemfonts_1.3.1      jquerylib_0.1.4        GenomeInfoDb_1.46.2   
#> [121] dichromat_2.0-0.1      Rcpp_1.1.1             globals_0.19.1        
#> [124] spatstat.random_3.4-4  RcppSpdlog_0.0.27      png_0.1-8             
#> [127] spatstat.univar_3.1-6  parallel_4.5.2         pkgdown_2.2.0         
#> [130] ggplot2_4.0.2          assertthat_0.2.1       dotCall64_1.2         
#> [133] listenv_0.10.1         spdl_0.0.5             viridisLite_0.4.3     
#> [136] scales_1.4.0           ggridges_0.5.7         crayon_1.5.3          
#> [139] purrr_1.2.1            rlang_1.1.7            cowplot_1.2.0         
#> [142] nanotime_0.3.13
```
