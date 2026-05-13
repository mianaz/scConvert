# Convert Between Seurat and SingleCellExperiment

Many Bioconductor workflows – scran normalization, scater quality
control, slingshot trajectory inference – operate on the
SingleCellExperiment (SCE) container. Seurat provides built-in
converters
([`as.SingleCellExperiment()`](https://satijalab.org/seurat/reference/as.SingleCellExperiment.html)
and
[`as.Seurat()`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md))
that make it straightforward to move data between the two ecosystems.
scConvert extends this by letting you route SCE objects to any file
format (h5ad, h5Seurat, loom, zarr, rds) through its hub architecture.

## Load demo data

scConvert ships a small PBMC dataset. We start by reading it into a
Seurat object.

``` r

pbmc <- readH5AD(system.file("extdata", "pbmc_demo.h5ad", package = "scConvert"))
pbmc
#> An object of class Seurat 
#> 2000 features across 500 samples within 1 assay 
#> Active assay: RNA (2000 features, 2000 variable features)
#>  2 layers present: counts, data
#>  2 dimensional reductions calculated: pca, umap
```

## Seurat to SingleCellExperiment

Convert the Seurat object to SCE. Counts, normalized data, dimensional
reductions, and cell metadata are all transferred.

``` r

library(SingleCellExperiment)

sce <- as.SingleCellExperiment(pbmc)
sce
#> class: SingleCellExperiment 
#> dim: 2000 500 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(2000): PPBP LYZ ... CLIC2 HEMGN
#> rowData names(0):
#> colnames(500): AAACATTGAGCTAC AAACATTGATCAGC ... TTTCGAACTCTCAT
#>   TTTGCATGCCTCAC
#> colData names(8): orig.ident nCount_RNA ... seurat_clusters ident
#> reducedDimNames(2): PCA UMAP
#> mainExpName: RNA
#> altExpNames(0):
```

Inspect what was carried over:

``` r

cat("Assays:", paste(assayNames(sce), collapse = ", "), "\n")
#> Assays: counts, logcounts
cat("Reduced dims:", paste(reducedDimNames(sce), collapse = ", "), "\n")
#> Reduced dims: PCA, UMAP
cat("colData columns:", paste(colnames(colData(sce)), collapse = ", "), "\n")
#> colData columns: orig.ident, nCount_RNA, nFeature_RNA, seurat_annotations, percent.mt, RNA_snn_res.0.5, seurat_clusters, ident
```

## SingleCellExperiment to Seurat

Convert the SCE back to a Seurat object with
[`as.Seurat()`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md).
The converter detects which assays exist (counts, logcounts) and maps
them to the appropriate Seurat layers.

``` r

seurat_rt <- as.Seurat(sce, counts = "counts", data = "logcounts")
seurat_rt
#> An object of class Seurat 
#> 2000 features across 500 samples within 1 assay 
#> Active assay: RNA (2000 features, 0 variable features)
#>  2 layers present: counts, data
#>  2 dimensional reductions calculated: PCA, UMAP
```

## Verify the round-trip

Confirm that expression values, cell barcodes, and embeddings survive
the Seurat -\> SCE -\> Seurat round-trip.

``` r

# Cell barcodes
cat("Barcodes match:", identical(Cells(pbmc), Cells(seurat_rt)), "\n")
#> Barcodes match: TRUE

# Gene names
cat("Features match:", identical(rownames(pbmc), rownames(seurat_rt)), "\n")
#> Features match: TRUE

# Count matrix values
orig_counts <- GetAssayData(pbmc, layer = "counts")
rt_counts <- GetAssayData(seurat_rt, layer = "counts")
cat("Counts identical:", identical(orig_counts, rt_counts), "\n")
#> Counts identical: TRUE

# Dimensionality reductions (SCE stores these in reducedDims)
cat("Reductions in original:", paste(names(pbmc@reductions), collapse = ", "), "\n")
#> Reductions in original: pca, umap
cat("Reductions after round-trip:", paste(names(seurat_rt@reductions), collapse = ", "), "\n")
#> Reductions after round-trip: PCA, UMAP
```

## SCE to h5ad via scConvert

To share an SCE dataset with Python collaborators, pass it directly to
[`scConvert()`](https://mianaz.github.io/scConvert/reference/scConvert.md).
The SCE is converted to Seurat internally, then written to the target
format – no manual intermediate step needed.

``` r

h5ad_path <- tempfile(fileext = ".h5ad")
scConvert(sce, dest = h5ad_path, overwrite = TRUE)
cat("h5ad file size:", round(file.size(h5ad_path) / 1e6, 1), "MB\n")
#> h5ad file size: 0.9 MB
```

Read it back to confirm the data survived the full SCE -\> Seurat -\>
h5ad -\> Seurat chain:

``` r

pbmc_h5ad <- readH5AD(h5ad_path)
cat("Cells:", ncol(pbmc_h5ad), "\n")
#> Cells: 500
cat("Genes:", nrow(pbmc_h5ad), "\n")
#> Genes: 2000
cat("Reductions:", paste(names(pbmc_h5ad@reductions), collapse = ", "), "\n")
#> Reductions: PCA, UMAP
```

## Seurat to SCE via scConvert

The reverse also works. Use `dest = "sce"` to get a SingleCellExperiment
directly from any Seurat object:

``` r

sce2 <- scConvert(pbmc, dest = "sce")
cat("Class:", class(sce2), "\n")
#> Class: SingleCellExperiment
cat("Dimensions:", nrow(sce2), "x", ncol(sce2), "\n")
#> Dimensions: 2000 x 500
```

## Practical example: Bioconductor QC on SCE, then back to Seurat

A common workflow is to use Bioconductor tools for quality control or
normalization, then bring the result back into Seurat for downstream
analysis. Here we compute log-library-size and detected-feature counts
using base R (the same quantities that `scater::addPerCellQCMetrics`
would compute).

``` r

# Compute QC metrics on the SCE
sce_qc <- as.SingleCellExperiment(pbmc)
counts_mat <- assay(sce_qc, "counts")
colData(sce_qc)$log_library_size <- log1p(Matrix::colSums(counts_mat))
colData(sce_qc)$detected_features <- Matrix::colSums(counts_mat > 0)

cat("Added QC columns:",
    paste(c("log_library_size", "detected_features"), collapse = ", "), "\n")
#> Added QC columns: log_library_size, detected_features

# Convert back to Seurat -- the new colData columns become metadata
seurat_qc <- as.Seurat(sce_qc, counts = "counts", data = "logcounts")
cat("Seurat metadata columns:",
    paste(colnames(seurat_qc[[]]), collapse = ", "), "\n")
#> Seurat metadata columns: orig.ident, nCount_RNA, nFeature_RNA, seurat_annotations, percent.mt, RNA_snn_res.0.5, seurat_clusters, ident, log_library_size, detected_features
```

The new QC columns (`log_library_size`, `detected_features`) are now
available as Seurat metadata for filtering or visualization.

## Clean up

``` r

unlink(h5ad_path)
```
