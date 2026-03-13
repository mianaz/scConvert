# Data Mapping Reference

## Overview

This reference documents exactly how scConvert maps data structures
between Seurat objects and each supported file format. Understanding
these mappings is essential for verifying conversion fidelity and
debugging unexpected behavior.

## h5ad (AnnData) Mapping

### Seurat -\> h5ad (via h5Seurat intermediate)

| Seurat Component | h5ad Location | Notes |
|----|----|----|
| `GetAssayData(layer = "data")` | `X` | Log-normalized expression (primary matrix) |
| `GetAssayData(layer = "counts")` | `raw/X` | Raw UMI counts (if available) |
| `GetAssayData(layer = "scale.data")` | *(skipped)* | Only variable features; recompute with `sc.pp.scale()` |
| [`VariableFeatures()`](https://satijalab.github.io/seurat-object/reference/VariableFeatures.html) | `var['highly_variable']` | Boolean column in var |
| `meta.data` | `obs` | Cell-level metadata |
| `meta.features` | `var` | Gene-level metadata |
| `Embeddings(, "pca")` | `obsm['X_pca']` | PCA coordinates |
| `Embeddings(, "umap")` | `obsm['X_umap']` | UMAP coordinates |
| `Embeddings(, "tsne")` | `obsm['X_tsne']` | t-SNE coordinates |
| `Graphs(, "RNA_snn")` | `obsp['connectivities']` | SNN graph |
| `Graphs(, "RNA_nn")` | `obsp['distances']` | KNN distances |
| [`GetTissueCoordinates()`](https://satijalab.github.io/seurat-object/reference/GetTissueCoordinates.html) | `obsm['spatial']` | Spatial spot coordinates |
| [`Images()`](https://satijalab.github.io/seurat-object/reference/Images.html) tissue image | `uns['spatial'][lib]['images']` | H&E images |
| Scale factors | `uns['spatial'][lib]['scalefactors']` | Visium pixel scaling |
| `misc` | `uns` | Unstructured metadata |

**Important notes:**

- `data` (all genes) is prioritized over `scale.data` (variable features
  only) for `X`
- Variable features are stored as boolean `highly_variable` in `var`,
  not as a separate list
- Scaled data is NOT stored in standard h5ad files; recompute with
  [`ScaleData()`](https://satijalab.org/seurat/reference/ScaleData.html)
  after loading

### h5ad -\> Seurat (readH5AD)

| h5ad Location | Seurat Destination | Condition |
|----|----|----|
| `X` | Default assay `data` layer | Always (log-normalized expected) |
| `raw/X` | `counts` layer | If `raw` group exists |
| `X` | `counts` layer | Fallback if no `raw` group |
| `obs` | `meta.data` | Categorical columns -\> R factors |
| `var` | Feature metadata (`meta.features`) | All columns preserved |
| `var['highly_variable']` | [`VariableFeatures()`](https://satijalab.github.io/seurat-object/reference/VariableFeatures.html) | Boolean TRUE -\> variable |
| `obsm/X_umap` | `reductions$umap` | Auto-detected by prefix |
| `obsm/X_pca` | `reductions$pca` | Auto-detected by prefix |
| `obsm/X_tsne` | `reductions$tsne` | Auto-detected by prefix |
| `obsm/X_draw_graph_*` | `reductions$draw_graph_*` | Force-directed layouts |
| `obsm/spatial` | Spatial coordinates | Via [`H5ADSpatialToSeurat()`](https://mianaz.github.io/scConvert/reference/H5ADSpatialToSeurat.md) |
| `obsp/connectivities` | `graphs$RNA_snn` | SNN graph |
| `obsp/distances` | `graphs$RNA_nn` | Distance graph |
| `uns/*` | `misc` | Unstructured annotations |

### Categorical encoding

AnnData stores categorical columns using HDF5 enum types:

    obs/cell_type/
      categories: ["B cell", "T cell", "Monocyte"]
      codes: [0, 2, 1, 0, ...]

scConvert decodes these to R factors with correct level mapping.
SeuratDisk handles this incorrectly (level corruption).

## h5Seurat Mapping

h5Seurat is a native Seurat format. The mapping is 1:1:

| Seurat Slot                    | h5Seurat Group/Dataset              |
|--------------------------------|-------------------------------------|
| `assays$RNA$counts`            | `/assays/RNA/counts` (sparse group) |
| `assays$RNA$data`              | `/assays/RNA/data` (sparse group)   |
| `assays$RNA$scale.data`        | `/assays/RNA/scale.data` (dense)    |
| `assays$RNA$features`          | `/assays/RNA/features`              |
| `assays$RNA$variable.features` | `/assays/RNA/variable.features`     |
| `assays$RNA$meta.features`     | `/assays/RNA/meta.features`         |
| `meta.data`                    | `/meta.data`                        |
| `reductions$pca`               | `/reductions/pca/cell.embeddings`   |
| `reductions$umap`              | `/reductions/umap/cell.embeddings`  |
| `graphs$RNA_snn`               | `/graphs/RNA_snn` (sparse)          |
| `images$anterior1`             | `/images/anterior1` (S4 object)     |
| `misc`                         | `/misc`                             |
| `commands`                     | `/commands`                         |

## h5mu (MuData) Mapping

### Seurat -\> h5mu (writeH5MU)

| Seurat Component      | h5mu Location                                 |
|-----------------------|-----------------------------------------------|
| Each assay            | `/mod/{modality}/` (one AnnData per modality) |
| Assay counts          | `/mod/{modality}/X`                           |
| Assay features        | `/mod/{modality}/var`                         |
| Cell metadata         | `/obs` (global)                               |
| Per-modality metadata | `/mod/{modality}/obs`                         |

**Modality name mapping:**

| Seurat Assay | h5mu Modality  |
|--------------|----------------|
| RNA          | rna            |
| ADT          | prot           |
| ATAC         | atac           |
| Spatial      | spatial        |
| SCT          | sct            |
| Other        | lowercase name |

### h5mu -\> Seurat (readH5MU)

Reverse mapping with automatic name conversion:

| h5mu Modality | Seurat Assay |
|---------------|--------------|
| rna           | RNA          |
| prot          | ADT          |
| atac          | ATAC         |
| spatial       | Spatial      |
| Other         | preserved    |

## Loom Mapping

### Seurat -\> Loom (writeLoom)

| Seurat Component    | Loom Location                         |
|---------------------|---------------------------------------|
| Default assay data  | `/matrix` (genes x cells, transposed) |
| Counts layer        | `/layers/counts`                      |
| Scale.data          | `/layers/scale.data`                  |
| Cell barcodes       | `/col_attrs/CellID`                   |
| Gene names          | `/row_attrs/Gene`                     |
| `meta.data` columns | `/col_attrs/*`                        |
| `meta.features`     | `/row_attrs/*`                        |
| Embeddings          | `/reductions/*/embeddings`            |
| Graphs              | `/col_graphs/*`                       |

### Loom -\> Seurat (readLoom)

| Loom Location       | Seurat Destination                                 |
|---------------------|----------------------------------------------------|
| `/matrix`           | `data` layer (or `counts`)                         |
| `/layers/*`         | Additional layers via `normalized`/`scaled` params |
| `/col_attrs/CellID` | Cell names (configurable)                          |
| `/row_attrs/Gene`   | Feature names (configurable)                       |
| `/col_attrs/*`      | `meta.data`                                        |
| `/row_attrs/*`      | `meta.features`                                    |

**Known limitations:**

- Loom stores the main matrix as normalized data; raw counts may not
  roundtrip perfectly
- Variable features are not natively supported
- PCA standard deviations are not preserved

## RDS Mapping

RDS is R’s native serialization format. The mapping is exact:

``` r

# Save
saveRDS(seurat_obj, "data.rds")

# Load
obj <- readRDS("data.rds")
```

All Seurat slots are preserved exactly, including BPCells on-disk
references (which may need path adjustment).

## SingleCellExperiment (SCE) Mapping

| Seurat Component  | SCE Component             |
|-------------------|---------------------------|
| `counts` layer    | `counts(sce)`             |
| `data` layer      | `logcounts(sce)`          |
| `meta.data`       | `colData(sce)`            |
| `meta.features`   | `rowData(sce)`            |
| `reductions$pca`  | `reducedDim(sce, "PCA")`  |
| `reductions$umap` | `reducedDim(sce, "UMAP")` |

## Metadata column name conventions

Common column names differ between ecosystems:

| Seurat (`meta.data`) | scanpy (`obs`) | Description |
|----|----|----|
| `seurat_clusters` | `leiden` / `louvain` | Cluster assignments |
| `orig.ident` | `batch` / `sample` | Sample ID |
| `nCount_RNA` | `n_counts` / `total_counts` | Total UMI per cell |
| `nFeature_RNA` | `n_genes` / `n_genes_by_counts` | Genes detected |
| `percent.mt` | `pct_counts_mt` | Mitochondrial fraction |
| `cell_type` | `cell_type` / `celltype` | Cell type annotation |
| `Phase` | `phase` / `cell_cycle_phase` | Cell cycle |

Use `standardize = TRUE` in
[`scConvert()`](https://mianaz.github.io/scConvert/reference/scConvert-package.html)
for automatic name translation.

## Indexing conventions

| Convention            | Python (h5ad)     | R (Seurat)          |
|-----------------------|-------------------|---------------------|
| Array indices         | 0-based           | 1-based             |
| Categorical codes     | 0-indexed         | 1-indexed factors   |
| Cluster labels        | Preserved as-is   | Preserved as-is     |
| Sparse matrix indices | 0-based (CSR/CSC) | 0-based (dgCMatrix) |

scConvert handles all index conversions automatically.

## Verification workflow

After any conversion, verify data integrity:

``` r

library(scConvert)

# Load original and roundtripped
orig <- readH5AD("original.h5ad")
writeH5Seurat(orig, "temp.h5Seurat", overwrite = TRUE)
scConvert("temp.h5Seurat", dest = "h5ad", overwrite = TRUE)
rt <- readH5AD("temp.h5ad")

# Check dimensions
stopifnot(ncol(orig) == ncol(rt))
stopifnot(nrow(orig) == nrow(rt))

# Check expression correlation
common_c <- intersect(colnames(orig), colnames(rt))
common_g <- intersect(rownames(orig), rownames(rt))
set.seed(42)
sc <- sample(common_c, 100)
sg <- sample(common_g, 100)
corr <- cor(
  as.numeric(GetAssayData(orig, layer = "data")[sg, sc]),
  as.numeric(GetAssayData(rt, layer = "data")[sg, sc])
)
cat("Correlation:", corr, "\n")  # Should be 1.0

# Check metadata
stopifnot(all(colnames(orig[[]]) %in% colnames(rt[[]])))

# Check reductions
stopifnot(all(names(orig@reductions) %in% names(rt@reductions)))
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
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.56         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.30    lifecycle_1.0.5   cli_3.6.5         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4   systemfonts_1.3.1
#> [17] compiler_4.5.2    tools_4.5.2       ragg_1.5.0        bslib_0.10.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.1.7       fs_1.6.6          htmlwidgets_1.6.4
```
