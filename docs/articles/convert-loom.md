# Conversions: Seurat and Loom

This vignette demonstrates how to convert between Seurat objects and
Loom files using scConvert. The [Loom format](http://loompy.org/) is an
HDF5-based file format commonly used for storing single-cell RNA-seq
data, particularly in RNA velocity workflows with tools like
[velocyto](http://velocyto.org/) and
[scVelo](https://scvelo.readthedocs.io/).

``` r

library(Seurat)
library(scConvert)
library(ggplot2)
```

## Loom File Format Overview

Loom files organize single-cell data in a specific HDF5 structure:

- `/matrix`: Main expression matrix (genes × cells, transposed from
  Seurat convention)
- `/row_attrs`: Gene/feature-level metadata (gene names, coordinates,
  etc.)
- `/col_attrs`: Cell/sample-level metadata (cell IDs, cluster labels,
  etc.)
- `/layers`: Additional expression layers (counts, normalized data,
  spliced/unspliced for velocity)

This format is particularly useful for:

- **RNA velocity analysis**: velocyto and scVelo store spliced/unspliced
  counts in layers
- **Python interoperability**: loompy package provides fast access in
  Python
- **Archival storage**: Self-contained format with all annotations

## Converting from Seurat to Loom

### Basic Conversion

To save a Seurat object as a Loom file, use
[`writeLoom()`](https://mianaz.github.io/scConvert/reference/writeLoom.md):

``` r

library(SeuratData)
if (!"pbmc3k.final" %in% rownames(InstalledData())) {
  InstallData("pbmc3k")
}

data("pbmc3k.final", package = "pbmc3k.SeuratData")
pbmc <- UpdateSeuratObject(pbmc3k.final)
pbmc
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  3 layers present: counts, data, scale.data
#>  2 dimensional reductions calculated: pca, umap
```

``` r

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![](convert-loom_files/figure-html/plot_pbmc_umap-1.png)

``` r

FeaturePlot(pbmc, features = "CD14", pt.size = 0.5)
```

![](convert-loom_files/figure-html/plot_pbmc_cd14-1.png)

Now save to Loom format:

``` r

writeLoom(pbmc, filename = "pbmc3k.loom", overwrite = TRUE, verbose = TRUE)
```

### What Gets Saved

When saving a Seurat object to Loom:

| Seurat Data | Loom Location | Notes |
|----|----|----|
| Default assay data | `/matrix` | Normalized expression |
| Counts layer | `/layers/counts` | Raw counts if different from data |
| Scale.data | `/layers/scale.data` | Scaled data if present |
| Cell names | `/col_attrs/CellID` | Cell barcodes |
| Gene names | `/row_attrs/Gene` | Feature names |
| Cell metadata | `/col_attrs/*` | All columns from `meta.data` |
| Feature metadata | `/row_attrs/*` | All columns from assay meta.features |
| Dimensional reductions | `/reductions/*` | PCA, UMAP embeddings, etc. |
| Graphs | `/col_graphs/*` | SNN graphs if present |

## Viewing Loom Files in Python

The saved Loom file can be opened with loompy or scanpy in Python:

``` python
import loompy

# Connect to the loom file
with loompy.connect("pbmc3k.loom") as ds:
    print(f"Shape: {ds.shape[0]} genes x {ds.shape[1]} cells")
    print(f"\nRow attributes (genes): {list(ds.ra.keys())}")
    print(f"\nColumn attributes (cells): {list(ds.ca.keys())[:10]}...")  # First 10
    print(f"\nLayers: {list(ds.layers.keys())}")
#> Shape: 13714 genes x 2638 cells
#> 
#> Row attributes (genes): ['Gene', 'vst.mean', 'vst.variable', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized']
#> 
#> Column attributes (cells): ['CellID', 'RNA_snn_res.0.5', 'nCount_RNA', 'nFeature_RNA', 'orig.ident', 'percent.mt', 'seurat_annotations', 'seurat_clusters']...
#> 
#> Layers: ['', 'counts', 'scale.data']
```

With scanpy:

``` python
import scanpy as sc

# Read loom file as AnnData
adata = sc.read_loom("pbmc3k.loom", sparse=True, cleanup=False)
print(adata)
#> AnnData object with n_obs × n_vars = 2638 × 13714
#>     obs: 'RNA_snn_res.0.5', 'nCount_RNA', 'nFeature_RNA', 'orig.ident', 'percent.mt', 'seurat_annotations', 'seurat_clusters'
#>     var: 'vst.mean', 'vst.variable', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized'
#>     layers: 'counts', 'scale.data'
print("\nColumn attributes:", list(adata.obs.columns)[:10])
#> 
#> Column attributes: ['RNA_snn_res.0.5', 'nCount_RNA', 'nFeature_RNA', 'orig.ident', 'percent.mt', 'seurat_annotations', 'seurat_clusters']
```

## Converting from Loom to Seurat

### Loading Loom Files

Use
[`readLoom()`](https://mianaz.github.io/scConvert/reference/readLoom.md)
to read a Loom file as a Seurat object:

``` r

# Load the loom file we just created
loaded_pbmc <- readLoom("pbmc3k.loom", verbose = TRUE)
loaded_pbmc
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 0 variable features)
#>  2 layers present: counts, data
#>  1 dimensional reduction calculated: umap
```

``` r

# Verify dimensions match
cat("Original:", ncol(pbmc), "cells,", nrow(pbmc), "genes\n")
#> Original: 2638 cells, 13714 genes
cat("Loaded:  ", ncol(loaded_pbmc), "cells,", nrow(loaded_pbmc), "genes\n")
#> Loaded:   2638 cells, 13714 genes

# Check metadata was preserved
cat("\nMetadata columns preserved:\n")
#> 
#> Metadata columns preserved:
print(intersect(colnames(pbmc[[]]), colnames(loaded_pbmc[[]])))
#> [1] "orig.ident"         "nCount_RNA"         "nFeature_RNA"      
#> [4] "seurat_annotations" "percent.mt"         "RNA_snn_res.0.5"   
#> [7] "seurat_clusters"
```

### Visualize loaded data

``` r

# Loom does not preserve UMAP embeddings, so FeaturePlot falls back to
# plotting the first two dimensions available (e.g., PC_1 vs PC_2),
# producing a PCA scatter rather than a UMAP.
FeaturePlot(loaded_pbmc, features = "CD14", pt.size = 0.5) + ggtitle("CD14 (from Loom)")
```

![](convert-loom_files/figure-html/plot_loom_feature-1.png)

### readLoom Parameters

[`readLoom()`](https://mianaz.github.io/scConvert/reference/readLoom.md)
provides several options for customizing the import:

``` r

readLoom(
  file,                    # Path to loom file

  assay = NULL,            # Name for the assay (default: "RNA" or from file)
  cells = "CellID",        # Column attribute containing cell names
  features = "Gene",       # Row attribute containing gene names
  normalized = NULL,       # Layer to load as normalized data
  scaled = NULL,           # Layer to load as scaled data
  filter = "none",         # Filter cells/features by Valid attributes
  verbose = TRUE           # Show progress messages
)
```

### Loading Specific Layers

If your Loom file has additional layers (e.g., from velocyto):

``` r

# Load with specific layers
seurat_obj <- readLoom(
  "velocity_data.loom",
  normalized = "spliced",      # Use spliced counts as normalized
  scaled = "ambiguous"         # Optional scaled layer
)
```

### Loading velocyto Output

Loom files from velocyto contain spliced, unspliced, and ambiguous count
matrices:

``` r

# Load velocyto loom file
# The main matrix typically contains spliced counts
velocity_obj <- readLoom(
  "sample.loom",
  cells = "CellID",
  features = "Gene"
)

# The spliced/unspliced/ambiguous layers can be accessed after loading
# or you may need to load them separately depending on your analysis needs
```

## Round-Trip Conversion

scConvert preserves data integrity during Seurat ↔︎ Loom conversion:

``` r

# Create a simple test object
set.seed(42)
test_obj <- CreateSeuratObject(
  counts = pbmc[["RNA"]]$counts[1:100, 1:50],
  project = "RoundTrip"
)
test_obj$custom_cluster <- sample(c("A", "B", "C"), 50, replace = TRUE)
test_obj$numeric_value <- rnorm(50)

# Save and reload
writeLoom(test_obj, "test_roundtrip.loom", overwrite = TRUE, verbose = FALSE)
reloaded <- readLoom("test_roundtrip.loom", verbose = FALSE)

# Compare
cat("Cell names match:", all(colnames(test_obj) == colnames(reloaded)), "\n")
#> Cell names match: TRUE
cat("Gene names match:", all(rownames(test_obj) == rownames(reloaded)), "\n")
#> Gene names match: TRUE
cat("Metadata columns:", paste(colnames(reloaded[[]]), collapse = ", "), "\n")
#> Metadata columns: orig.ident, nCount_RNA, nFeature_RNA, custom_cluster, numeric_value
```

Verify expression data is preserved:

``` r

original_data <- GetAssayData(test_obj, layer = "data")
reloaded_data <- GetAssayData(reloaded, layer = "data")[
  rownames(original_data),
  colnames(original_data)
]

max_diff <- max(abs(as.matrix(original_data) - as.matrix(reloaded_data)))
cat("Maximum expression difference:", max_diff, "\n")
#> Maximum expression difference: -Inf
```

## Working with RNA Velocity Data

A common use case for Loom files is RNA velocity analysis. Here’s a
typical workflow:

### 1. Generate Loom with velocyto

``` bash
# Run velocyto on Cell Ranger output (command line)
velocyto run10x -m repeat_mask.gtf sample_dir genes.gtf
```

This creates a `.loom` file with spliced/unspliced counts.

### 2. Load in R for QC and Annotation

``` r

# Load velocyto output
velocity_data <- readLoom("sample.loom")

# Perform standard Seurat QC and clustering
velocity_data <- NormalizeData(velocity_data)
velocity_data <- FindVariableFeatures(velocity_data)
velocity_data <- ScaleData(velocity_data)
velocity_data <- RunPCA(velocity_data)
velocity_data <- FindNeighbors(velocity_data)
velocity_data <- FindClusters(velocity_data)
velocity_data <- RunUMAP(velocity_data, dims = 1:30)

# Add cell type annotations
velocity_data$cell_type <- ...  # Your annotation method

# Save back to loom with annotations
writeLoom(velocity_data, "sample_annotated.loom", overwrite = TRUE)
```

### 3. Continue Velocity Analysis in Python

``` python
import scvelo as scv

# Load annotated loom file
adata = scv.read("sample_annotated.loom")

# Run velocity analysis
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

# Visualize with Seurat annotations
scv.pl.velocity_embedding_stream(adata, color='cell_type')
```

## Data Mapping Reference

### Seurat to Loom

| Seurat Location | Loom Location | Notes |
|----|----|----|
| `GetAssayData(layer = "data")` | `/matrix` | Main expression matrix |
| `GetAssayData(layer = "counts")` | `/layers/counts` | If different from data |
| `GetAssayData(layer = "scale.data")` | `/layers/scale.data` | If present |
| `colnames(obj)` | `/col_attrs/CellID` | Cell barcodes |
| `rownames(obj)` | `/row_attrs/Gene` | Gene names |
| `obj[[]]` (meta.data) | `/col_attrs/*` | Each column as attribute |
| `obj[[assay]][[]]` | `/row_attrs/*` | Feature metadata |
| `Embeddings(obj, "pca")` | `/reductions/pca/embeddings` | Transposed |
| `Loadings(obj, "pca")` | `/reductions/pca/loadings` | If present |
| `Stdev(obj, "pca")` | `/reductions/pca/stdev` | If present |

### Loom to Seurat

| Loom Location | Seurat Destination | Notes |
|----|----|----|
| `/matrix` | `data` layer | Stored as counts if no normalization detected |
| `/layers/*` | Additional layers | Via `normalized`/`scaled` parameters |
| `/col_attrs/CellID` | Cell names | Configurable via `cells` parameter |
| `/row_attrs/Gene` | Feature names | Configurable via `features` parameter |
| `/col_attrs/*` | `meta.data` | All except CellID and Valid |
| `/row_attrs/*` | `meta.features` | All except Gene and Valid |
| `/reductions/*/embeddings` | [`Reductions()`](https://satijalab.github.io/seurat-object/reference/ObjectAccess.html) | If present |
| `/col_graphs/*` | [`Graphs()`](https://satijalab.github.io/seurat-object/reference/ObjectAccess.html) | SNN/KNN graphs |

## Comparison with Other Formats

| Feature           | Loom             | h5Seurat       | h5ad            |
|-------------------|------------------|----------------|-----------------|
| Primary ecosystem | velocyto, loompy | Seurat         | scanpy          |
| Multiple assays   | Via layers       | Native support | Single X matrix |
| Spatial data      | Limited          | Full support   | Full support    |
| RNA velocity      | Native           | Not standard   | Via layers      |
| Graph storage     | Native           | Native         | In obsp         |
| Python access     | loompy, scanpy   | Limited        | scanpy          |
| R access          | scConvert        | scConvert      | scConvert       |

**When to use Loom:**

- RNA velocity analysis with velocyto/scVelo
- Workflows requiring loompy compatibility
- Simpler single-assay datasets

**When to use h5Seurat:**

- Seurat-native workflows
- Multi-modal data (CITE-seq)
- Full Seurat object preservation

**When to use h5ad:**

- scanpy-based workflows
- Spatial transcriptomics with squidpy
- CellxGene data sharing

## Troubleshooting

### Common Issues

**“Cannot find feature names dataset”**

The Loom file uses non-standard attribute names. Specify them
explicitly:

``` r

# Check what attributes exist
h5 <- hdf5r::H5File$new("data.loom", mode = "r")
print(names(h5[["row_attrs"]]))
h5$close_all()

# Use the correct attribute name
obj <- readLoom("data.loom", features = "gene_name")
```

**“Cannot find cell names dataset”**

Similar to above, check column attributes:

``` r

h5 <- hdf5r::H5File$new("data.loom", mode = "r")
print(names(h5[["col_attrs"]]))
h5$close_all()

# Use the correct attribute name
obj <- readLoom("data.loom", cells = "obs_names")
```

**Duplicate feature names**

Loom files sometimes have duplicate gene names. scConvert will make them
unique with a warning:

``` r

# Warning: Duplicate feature names found, making unique
obj <- readLoom("data.loom")

# The names will be Gene, Gene.1, Gene.2, etc.
```

## Direct Streaming Converters

In addition to
[`readLoom()`](https://mianaz.github.io/scConvert/reference/readLoom.md)
and
[`writeLoom()`](https://mianaz.github.io/scConvert/reference/writeLoom.md),
scConvert provides direct streaming converters that bypass Seurat object
construction. These are faster for format-to-format conversion because
they operate directly on the HDF5 files:

``` r

# Loom <-> h5Seurat (direct HDF5 streaming)
LoomToH5Seurat("data.loom", "data.h5seurat")
H5SeuratToLoom("data.h5seurat", "data.loom")

# Loom <-> h5ad (via temp h5seurat, no Seurat object)
LoomToH5AD("data.loom", "data.h5ad")
H5ADToLoom("data.h5ad", "data.loom")

# Loom <-> Zarr
LoomToZarr("data.loom", "data.zarr")
ZarrToLoom("data.zarr", "data.loom")

# Loom <-> h5mu
LoomToH5MU("data.loom", "data.h5mu")
H5MUToLoom("data.h5mu", "data.loom")
```

The C CLI binary also supports Loom for maximum speed:

``` bash
scconvert data.h5ad data.loom
scconvert data.loom data.h5seurat
```

## See Also

- [Command-Line
  Interface](https://mianaz.github.io/scConvert/articles/cli-usage.md) –
  convert Loom files from the shell or via
  [`scConvert_cli()`](https://mianaz.github.io/scConvert/reference/scConvert_cli.md)
- [Data Mapping
  Reference](https://mianaz.github.io/scConvert/articles/data-mapping-reference.md)
  – complete field-by-field mapping tables for all formats

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
#>  [1] ggplot2_4.0.2                 scConvert_0.1.0              
#>  [3] Seurat_5.4.0                  SeuratObject_5.3.0           
#>  [5] sp_2.2-1                      stxKidney.SeuratData_0.1.0   
#>  [7] stxBrain.SeuratData_0.1.2     ssHippo.SeuratData_3.1.4     
#>  [9] pbmcref.SeuratData_1.0.0      pbmcMultiome.SeuratData_0.1.4
#> [11] pbmc3k.SeuratData_3.1.4       panc8.SeuratData_3.0.2       
#> [13] cbmc.SeuratData_3.1.4         SeuratData_0.2.2.9002        
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
#>  [34] shiny_1.13.0           digest_0.6.39          patchwork_1.3.2       
#>  [37] tensor_1.5.1           RSpectra_0.16-2        irlba_2.3.7           
#>  [40] textshaping_1.0.4      labeling_0.4.3         progressr_0.18.0      
#>  [43] spatstat.sparse_3.1-0  httr_1.4.8             polyclip_1.10-7       
#>  [46] abind_1.4-8            compiler_4.5.2         bit64_4.6.0-1         
#>  [49] withr_3.0.2            S7_0.2.1               fastDummies_1.7.5     
#>  [52] MASS_7.3-65            rappdirs_0.3.4         tools_4.5.2           
#>  [55] lmtest_0.9-40          otel_0.2.0             httpuv_1.6.16         
#>  [58] future.apply_1.20.2    goftest_1.2-3          glue_1.8.0            
#>  [61] nlme_3.1-168           promises_1.5.0         grid_4.5.2            
#>  [64] Rtsne_0.17             cluster_2.1.8.2        reshape2_1.4.5        
#>  [67] generics_0.1.4         hdf5r_1.3.12           gtable_0.3.6          
#>  [70] spatstat.data_3.1-9    tidyr_1.3.2            data.table_1.18.2.1   
#>  [73] spatstat.geom_3.7-0    RcppAnnoy_0.0.23       ggrepel_0.9.7         
#>  [76] RANN_2.6.2             pillar_1.11.1          stringr_1.6.0         
#>  [79] spam_2.11-3            RcppHNSW_0.6.0         later_1.4.8           
#>  [82] splines_4.5.2          dplyr_1.2.0            lattice_0.22-9        
#>  [85] bit_4.6.0              survival_3.8-6         deldir_2.0-4          
#>  [88] tidyselect_1.2.1       miniUI_0.1.2           pbapply_1.7-4         
#>  [91] knitr_1.51             gridExtra_2.3          scattermore_1.2       
#>  [94] xfun_0.56              matrixStats_1.5.0      stringi_1.8.7         
#>  [97] lazyeval_0.2.2         yaml_2.3.12            evaluate_1.0.5        
#> [100] codetools_0.2-20       tibble_3.3.1           cli_3.6.5             
#> [103] uwot_0.2.4             xtable_1.8-8           reticulate_1.45.0     
#> [106] systemfonts_1.3.1      jquerylib_0.1.4        dichromat_2.0-0.1     
#> [109] Rcpp_1.1.1             globals_0.19.1         spatstat.random_3.4-4 
#> [112] png_0.1-8              spatstat.univar_3.1-6  parallel_4.5.2        
#> [115] pkgdown_2.2.0          dotCall64_1.2          listenv_0.10.1        
#> [118] viridisLite_0.4.3      scales_1.4.0           ggridges_0.5.7        
#> [121] purrr_1.2.1            crayon_1.5.3           rlang_1.1.7           
#> [124] cowplot_1.2.0
```
