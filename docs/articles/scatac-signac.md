# scATAC-seq Data Conversion with Signac

## Introduction

Single-cell ATAC-seq (scATAC-seq) measures chromatin accessibility at
single-cell resolution. The primary data structure is a **peak-by-cell
count matrix**, where each row represents a genomic region (peak) and
each value indicates how many Tn5 insertion events were observed in that
region for a given cell.

### scATAC-seq data components

| Component | Description | Format |
|----|----|----|
| Peak matrix | Peaks x cells count matrix | Sparse matrix (HDF5, MEX) |
| Fragment file | Per-read genomic coordinates | `fragments.tsv.gz` + `.tbi` index |
| Cell metadata | Per-cell QC and annotations | Data frame |
| Gene annotations | Genomic coordinates for genes | GRanges (GTF/GFF) |
| Motif matrix | TF binding motif enrichment | Matrix in ChromatinAssay |

In R, the [Signac](https://stuartlab.org/signac/) package extends Seurat
with the `ChromatinAssay` class, which wraps the peak matrix alongside
fragment file paths, gene annotations, and motif information. In Python,
scATAC-seq data is typically stored as an AnnData object (h5ad) with
peaks as `var` and cells as `obs`, following the same conventions as
scRNA-seq.

### What scConvert handles

scConvert converts the **peak count matrix** and associated **cell
metadata** and **embeddings** between R and Python formats.
Signac-specific annotations (fragment files, motif matrices, gene
annotations) are not part of the h5ad specification and require separate
handling.

``` r

library(Seurat)
library(Signac)
library(scConvert)
```

## Download a real scATAC-seq dataset

We use the 10x Genomics PBMC 10K scATAC-seq dataset (Chromium
Controller, ATAC v2 chemistry). The filtered peak-barcode matrix is
distributed as an HDF5 file that can be read directly by
[`Seurat::Read10X_h5()`](https://satijalab.org/seurat/reference/Read10X_h5.html).

``` r

atac_url <- paste0(
 "https://cf.10xgenomics.com/samples/cell-atac/2.1.0/",
 "10k_pbmc_ATACv2_nextgem_Chromium_Controller/",
 "10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5"
)

if (!file.exists("pbmc_atac_peaks.h5")) {
  download.file(atac_url, destfile = "pbmc_atac_peaks.h5", mode = "wb")
}

cat("File size:", round(file.size("pbmc_atac_peaks.h5") / 1024^2, 1), "MB\n")
```

## Load into Signac ChromatinAssay

[`Read10X_h5()`](https://satijalab.org/seurat/reference/Read10X_h5.html)
reads the filtered peak-barcode matrix. We then create a
`ChromatinAssay` and wrap it in a Seurat object.

``` r

# Read the 10x HDF5 peak matrix
counts <- Read10X_h5("pbmc_atac_peaks.h5")

# Create a ChromatinAssay (without fragment file or annotations for portability)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  min.cells = 10,
  min.features = 200
)

# Create Seurat object
pbmc_atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC"
)

cat("Cells:", ncol(pbmc_atac), "\n")
cat("Peaks:", nrow(pbmc_atac), "\n")
cat("Assay class:", class(pbmc_atac[["ATAC"]])[1], "\n")

# Inspect peak names (genomic coordinates)
head(rownames(pbmc_atac), 5)
#> [1] "chr1:565107-565550"  "chr1:569174-569639" ...
```

### Basic QC and processing

``` r

# QC metrics
pbmc_atac$peak_region_fragments <- colSums(GetAssayData(pbmc_atac, layer = "counts"))
pbmc_atac$pct_reads_in_peaks <- pbmc_atac$peak_region_fragments / pbmc_atac$nCount_ATAC

# Normalization (TF-IDF is standard for scATAC)
pbmc_atac <- RunTFIDF(pbmc_atac)
pbmc_atac <- FindTopFeatures(pbmc_atac, min.cutoff = "q0")
pbmc_atac <- RunSVD(pbmc_atac)

# UMAP (skip first LSI component -- it correlates with sequencing depth)
pbmc_atac <- RunUMAP(pbmc_atac, reduction = "lsi", dims = 2:30)
pbmc_atac <- FindNeighbors(pbmc_atac, reduction = "lsi", dims = 2:30)
pbmc_atac <- FindClusters(pbmc_atac, resolution = 0.5, algorithm = 3)

DimPlot(pbmc_atac, label = TRUE, pt.size = 0.3) + NoLegend()
```

## Convert Signac to h5ad

The
[`writeH5AD()`](https://mianaz.github.io/scConvert/reference/writeH5AD.md)
function writes the peak matrix, cell metadata, and embeddings directly
to h5ad format. The `ChromatinAssay` peak matrix is extracted as a
standard sparse matrix during conversion.

``` r

writeH5AD(pbmc_atac, filename = "pbmc_atac.h5ad", overwrite = TRUE)
cat("h5ad file size:", round(file.size("pbmc_atac.h5ad") / 1024^2, 1), "MB\n")
```

The resulting h5ad file can be loaded in Python with scanpy or
episcanpy:

``` python
import scanpy as sc

adata = sc.read_h5ad("pbmc_atac.h5ad")
print(adata)
print(f"\nPeak names (first 5): {list(adata.var_names[:5])}")
print(f"obsm keys: {list(adata.obsm.keys())}")
print(f"obs columns: {list(adata.obs.columns)[:10]}")

if "X_umap" in adata.obsm:
    sc.pl.umap(adata, color="seurat_clusters", title="scATAC UMAP — scanpy (from scConvert h5ad)")
```

### CLI equivalent

``` r

scConvert_cli("pbmc_atac.h5ad", "pbmc_atac_cli.h5ad")
```

Or from the command line:

``` bash
scconvert pbmc_atac.h5ad pbmc_atac_cli.h5ad
```

## Convert h5ad back to Seurat with Signac

When loading an h5ad file that contains a peak matrix,
[`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
returns a standard Seurat object. You can then upgrade the assay to a
Signac `ChromatinAssay` if needed.

``` r

# Load the h5ad as a standard Seurat object
atac_loaded <- readH5AD("pbmc_atac.h5ad")

cat("Cells:", ncol(atac_loaded), "\n")
cat("Peaks:", nrow(atac_loaded), "\n")
cat("Assay class:", class(atac_loaded[["ATAC"]])[1], "\n")

# Inspect: peak names should be chr:start-end format
head(rownames(atac_loaded), 5)
```

### Upgrade to ChromatinAssay

If you need Signac functionality (e.g., motif analysis, coverage plots),
upgrade the loaded assay to a `ChromatinAssay`. This requires
re-attaching any Signac-specific annotations.

``` r

# Extract the counts matrix from the loaded object
peak_counts <- GetAssayData(atac_loaded, layer = "counts")

# Create a ChromatinAssay from the loaded data
chrom_assay_rt <- CreateChromatinAssay(
  counts = peak_counts,
  sep = c(":", "-"),
  min.cells = 0,
  min.features = 0
)

# Replace the assay in the loaded object
atac_loaded[["ATAC"]] <- chrom_assay_rt
cat("Upgraded assay class:", class(atac_loaded[["ATAC"]])[1], "\n")

# Optional: re-attach fragment file if available locally
# Fragments(atac_loaded) <- CreateFragmentObject("fragments.tsv.gz")

# Optional: add gene annotations
# library(EnsDb.Hsapiens.v86)
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# Annotation(atac_loaded) <- annotations
```

### CLI equivalent

``` r

scConvert_cli("pbmc_atac.h5ad", "pbmc_atac_roundtrip.rds")
```

## Round-trip fidelity verification

After converting Signac -\> h5ad -\> Seurat, we verify that key data
structures are preserved.

``` r

# Reload from h5ad
atac_rt <- readH5AD("pbmc_atac.h5ad")

# --- Dimension check ---
cat("=== Dimensions ===\n")
cat("Original:  ", ncol(pbmc_atac), "cells x", nrow(pbmc_atac), "peaks\n")
cat("Roundtrip: ", ncol(atac_rt), "cells x", nrow(atac_rt), "peaks\n")
stopifnot(ncol(pbmc_atac) == ncol(atac_rt))
stopifnot(nrow(pbmc_atac) == nrow(atac_rt))
cat("Dimension check PASSED\n\n")

# --- Peak name preservation ---
cat("=== Peak Names ===\n")
orig_peaks <- sort(rownames(pbmc_atac))
rt_peaks <- sort(rownames(atac_rt))
cat("Peaks match:", identical(orig_peaks, rt_peaks), "\n")
stopifnot(identical(orig_peaks, rt_peaks))
cat("Peak name check PASSED\n\n")

# --- Cell barcode preservation ---
cat("=== Cell Barcodes ===\n")
orig_cells <- sort(colnames(pbmc_atac))
rt_cells <- sort(colnames(atac_rt))
cat("Barcodes match:", identical(orig_cells, rt_cells), "\n")
stopifnot(identical(orig_cells, rt_cells))
cat("Barcode check PASSED\n\n")

# --- Expression correlation ---
cat("=== Count Matrix Fidelity ===\n")
set.seed(42)
sc <- sample(colnames(pbmc_atac), min(500, ncol(pbmc_atac)))
sp <- sample(rownames(pbmc_atac), min(200, nrow(pbmc_atac)))
corr <- cor(
  as.numeric(GetAssayData(pbmc_atac, layer = "counts")[sp, sc]),
  as.numeric(GetAssayData(atac_rt, layer = "counts")[sp, sc])
)
cat("Count matrix correlation:", corr, "\n")
stopifnot(corr > 0.999)
cat("Count matrix check PASSED\n\n")

# --- Metadata preservation ---
cat("=== Cell Metadata ===\n")
orig_meta_cols <- sort(colnames(pbmc_atac[[]]))
rt_meta_cols <- sort(colnames(atac_rt[[]]))
shared_cols <- intersect(orig_meta_cols, rt_meta_cols)
cat("Original metadata columns:", length(orig_meta_cols), "\n")
cat("Roundtrip metadata columns:", length(rt_meta_cols), "\n")
cat("Shared columns:", length(shared_cols), "\n")

# --- Embedding preservation ---
cat("\n=== Embeddings ===\n")
if ("umap" %in% Reductions(atac_rt)) {
  orig_umap <- Embeddings(pbmc_atac, "umap")
  rt_umap <- Embeddings(atac_rt, "umap")
  umap_rmse <- sqrt(mean((orig_umap[sc, ] - rt_umap[sc, ])^2))
  cat("UMAP RMSE:", umap_rmse, "\n")
  stopifnot(umap_rmse < 1e-6)
  cat("UMAP check PASSED\n")
} else {
  cat("UMAP not found in roundtrip object (may need recomputation)\n")
}

if ("lsi" %in% Reductions(atac_rt)) {
  cat("LSI reduction preserved\n")
} else {
  cat("LSI not found in roundtrip object (stored as obsm['X_lsi'] in h5ad)\n")
}
```

## Multiome: RNA + ATAC to h5mu

For 10x Multiome experiments (simultaneous RNA + ATAC from the same
cells), the h5mu format stores both modalities in a single file. Each
assay becomes a separate modality.

``` r

# Example: create a combined RNA + ATAC Seurat object
# In practice, you would load from Cell Ranger ARC output:
#   rna_counts <- Read10X_h5("filtered_feature_bc_matrix.h5", ...)
# Here we simulate a small RNA assay for demonstration:

# Add a dummy RNA assay (in real workflows, this comes from Multiome data)
set.seed(42)
n_cells <- ncol(pbmc_atac)
dummy_rna <- matrix(
  rpois(500 * n_cells, lambda = 2),
  nrow = 500, ncol = n_cells,
  dimnames = list(
    paste0("Gene", seq_len(500)),
    colnames(pbmc_atac)
  )
)
pbmc_atac[["RNA"]] <- CreateAssay5Object(counts = as(dummy_rna, "dgCMatrix"))

cat("Assays:", paste(Assays(pbmc_atac), collapse = ", "), "\n")
cat("ATAC peaks:", nrow(pbmc_atac[["ATAC"]]), "\n")
cat("RNA genes:", nrow(pbmc_atac[["RNA"]]), "\n")

# Export both modalities to h5mu
writeH5MU(pbmc_atac, filename = "pbmc_atac_multiome.h5mu", overwrite = TRUE)
cat("h5mu file size:", round(file.size("pbmc_atac_multiome.h5mu") / 1024^2, 1), "MB\n")
```

### Load h5mu back

``` r

multiome_rt <- readH5MU("pbmc_atac_multiome.h5mu")

cat("Assays:", paste(Assays(multiome_rt), collapse = ", "), "\n")
cat("ATAC peaks:", nrow(multiome_rt[["ATAC"]]), "\n")
cat("RNA genes:", nrow(multiome_rt[["RNA"]]), "\n")
cat("Cells:", ncol(multiome_rt), "\n")
```

### CLI equivalent

``` r

scConvert_cli("multiome.rds", "multiome.h5mu")
```

## What is and is not preserved

### Preserved during conversion

| Component | Seurat/Signac | h5ad/h5mu | Notes |
|----|----|----|----|
| Peak count matrix | `counts` layer | `X` | Sparse matrix, lossless |
| Normalized matrix | `data` layer | `X` or `layers` | TF-IDF values |
| Cell metadata | `meta.data` | `obs` | All columns preserved |
| Peak metadata | `meta.features` | `var` | Feature-level annotations |
| LSI embeddings | `reductions$lsi` | `obsm['X_lsi']` | Latent semantic indexing |
| UMAP coordinates | `reductions$umap` | `obsm['X_umap']` | Visualization coords |
| Variable features | [`VariableFeatures()`](https://satijalab.github.io/seurat-object/reference/VariableFeatures.html) | `var['highly_variable']` | Top accessible peaks |
| Cluster assignments | `meta.data$seurat_clusters` | `obs['seurat_clusters']` | Cell labels |

### NOT preserved (Signac-specific)

| Component | Why | Workaround |
|----|----|----|
| Fragment file paths | Local filesystem reference, not data | Copy `fragments.tsv.gz` + `.tbi` separately; re-attach with [`CreateFragmentObject()`](https://stuartlab.org/signac/reference/CreateFragmentObject.html) |
| Gene annotations | GRanges object, no h5ad equivalent | Re-add with [`Annotation()`](https://stuartlab.org/signac/reference/Annotation.html) after loading (e.g., from EnsDb.Hsapiens.v86) |
| Motif matrices | Signac-specific slot | Re-compute with [`AddMotifs()`](https://stuartlab.org/signac/reference/AddMotifs.html) using a PWM database (e.g., JASPAR) |
| ChromatinAssay class | R-specific S4 class | Upgrade with [`CreateChromatinAssay()`](https://stuartlab.org/signac/reference/CreateChromatinAssay.html) after loading (see above) |
| Links (peak-gene) | Signac-specific slot | Re-compute with [`LinkPeaks()`](https://stuartlab.org/signac/reference/LinkPeaks.html) after re-attaching annotations |

### Practical recommendations

1.  **Always transfer fragment files separately.** They are large (often
    \>10 GB) and are not embedded in h5ad. Copy both `fragments.tsv.gz`
    and `fragments.tsv.gz.tbi`.

2.  **Re-attach annotations after loading.** Gene annotations and motif
    databases are organism-specific and best loaded from canonical
    sources (Ensembl, JASPAR) rather than serialized.

3.  **Peak names encode genomic coordinates.** The `chr:start-end`
    format is preserved as `var_names` in h5ad, so peak identity
    survives round-trip without loss.

4.  **For Multiome data, prefer h5mu.** It keeps RNA and ATAC in
    separate modalities within a single file, matching the muon/MuData
    convention in Python.

## See Also

- [Multimodal: CITE-seq and
  Multiome](https://mianaz.github.io/scConvert/articles/multimodal-citeseq.md)
  – combining RNA + ATAC or RNA + ADT in a single object

## Session Info

``` r

sessionInfo()
```
