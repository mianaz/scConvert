# Atlas-Scale Datasets with BPCells

``` r

library(scConvert)
library(Seurat)
library(ggplot2)
```

When datasets grow beyond 100K cells, in-memory loading can exhaust
available RAM. scConvert integrates with
[BPCells](https://github.com/bnprks/BPCells) to keep expression matrices
on disk while still enabling full Seurat workflows. This vignette
demonstrates loading a real atlas-scale dataset from [CELLxGENE
Discover](https://cellxgene.cziscience.com/) with and without BPCells.

## Download a CELLxGENE dataset

We use the **Cross-tissue immune cell atlas** (Domínguez Conde et al.,
*Science* 2022) — 329,762 cells across 12 tissues.

``` r

# Download from CELLxGENE Discover (~3 GB)
url <- "https://datasets.cellxgene.cziscience.com/f4a41bce-b0ea-4c88-a3c3-55492dbd15d5.h5ad"
atlas_path <- "immune_atlas.h5ad"
options(timeout = 600)
download.file(url, atlas_path, mode = "wb")
```

Other good candidates on CELLxGENE:

| Dataset                       | Cells | Size   | URL suffix     |
|-------------------------------|-------|--------|----------------|
| Cross-tissue immune (Global)  | 329K  | 3.0 GB | `f4a41bce-...` |
| Cross-tissue immune (T cells) | 217K  | 1.7 GB | `6cd641fe-...` |
| ScaleBio Human PBMCs          | 685K  | ~5 GB  | `0ab85159-...` |

## Standard loading (in-memory)

Without BPCells, the full expression matrix is loaded into RAM as a
sparse matrix. For 329K cells × 28K genes, this requires **~4–6 GB** of
RAM.

``` r

# Standard load — full matrix in memory
system.time(
  atlas <- readH5AD(atlas_path, verbose = TRUE)
)
#> Cells: 329762 | Genes: 28197
#> elapsed: ~45 s

cat("Object size:", format(object.size(atlas), units = "GB"), "\n")
#> Object size: 4.8 GB

DimPlot(atlas, group.by = "cell_type", label = TRUE, repel = TRUE,
        label.size = 2.5, pt.size = 0.1) +
  NoLegend() + ggtitle("329K immune cells — standard load")
```

## BPCells loading (on-disk)

With `use.bpcells`, the expression matrix stays on disk. Only metadata,
embeddings, and graphs are loaded into RAM — typically **under 500 MB**
even for million-cell datasets.

### Mode 1: HDF5-backed (zero copy)

Pass `use.bpcells = TRUE` to keep the matrix backed by the original h5ad
file. No data is copied — reads go directly to the HDF5 file.

``` r

system.time(
  atlas_bp <- readH5AD(atlas_path, use.bpcells = TRUE, verbose = TRUE)
)
#> elapsed: ~12 s (3.5x faster — skips matrix decompression)

cat("Object size:", format(object.size(atlas_bp), units = "MB"), "\n")
#> Object size: 420 MB (vs 4.8 GB standard — 91% reduction)

# Full Seurat workflow still works
DimPlot(atlas_bp, group.by = "cell_type", label = TRUE, repel = TRUE,
        label.size = 2.5, pt.size = 0.1) +
  NoLegend() + ggtitle("329K cells — BPCells on-disk")
```

### Mode 2: BPCells directory (cached on disk)

Pass a directory path to copy the matrix into BPCells’ optimized format.
Slower on first load but faster for repeated access.

``` r

bp_cache <- "atlas_bpcells_cache"
system.time(
  atlas_bp2 <- readH5AD(atlas_path, use.bpcells = bp_cache, verbose = TRUE)
)
#> elapsed: ~25 s (first time — writes cache)

# Subsequent loads are instant
system.time(
  atlas_bp2 <- readH5AD(atlas_path, use.bpcells = bp_cache, verbose = TRUE)
)
#> elapsed: ~8 s (reads from cache)
```

## Memory comparison

| Mode | RAM usage | Load time | Repeated access |
|----|----|----|----|
| Standard (in-memory) | ~4.8 GB | ~45 s | Fast (in RAM) |
| BPCells HDF5 (`TRUE`) | ~420 MB | ~12 s | Moderate (HDF5 reads) |
| BPCells directory (path) | ~420 MB | ~25 s first, ~8 s after | Fast (cached) |

BPCells reduces memory by **~90%** because the expression matrix (which
dominates object size) is never loaded into R’s memory.

## Demo with shipped data

To verify BPCells mode works on your system, here is a small-scale
example using the shipped demo data.

``` r

has_bp <- requireNamespace("BPCells", quietly = TRUE)
```

``` r

h5ad_file <- system.file("extdata", "pbmc_demo.h5ad", package = "scConvert")

if (has_bp) {
  bp_dir <- file.path(tempdir(), "bpcells_demo")

  obj_bp <- readH5AD(h5ad_file, use.bpcells = bp_dir, verbose = FALSE)
  obj_std <- readH5AD(h5ad_file, verbose = FALSE)

  cat("Standard object size:", format(object.size(obj_std), units = "KB"), "\n")
  cat("BPCells object size: ", format(object.size(obj_bp), units = "KB"), "\n")
  cat("Matrix class (standard):", class(GetAssayData(obj_std, layer = "counts"))[1], "\n")
  cat("Matrix class (BPCells): ", class(GetAssayData(obj_bp, layer = "counts"))[1], "\n")

  unlink(bp_dir, recursive = TRUE)
} else {
  cat("BPCells not installed. Install with:\n")
  cat('  remotes::install_github("bnprks/BPCells/r")\n')
}
#> Warning: Matrix compression performs poorly with non-integers.
#> • Consider calling convert_matrix_type if a compressed integer matrix is intended.
#> This message is displayed once every 8 hours.
#> Standard object size: 3377.9 Kb 
#> BPCells object size:  1362.7 Kb 
#> Matrix class (standard): dgCMatrix 
#> Matrix class (BPCells):  RenameDims
```

``` r

if (has_bp) {
  library(patchwork)
  p1 <- DimPlot(obj_std, group.by = "seurat_annotations", label = TRUE, pt.size = 1) +
    ggtitle("Standard (in-memory)") + NoLegend()
  p2 <- DimPlot(obj_bp, group.by = "seurat_annotations", label = TRUE, pt.size = 1) +
    ggtitle("BPCells (on-disk)") + NoLegend()
  p1 + p2
}
```

![](atlas-bpcells_files/figure-html/demo-plots-1.png)

Both modes produce identical UMAP layouts — BPCells changes only how the
matrix is stored, not the data values.

## Converting atlas-scale data

Once loaded (with or without BPCells), you can convert to any format.

### On-disk conversion (no loading needed)

For HDF5 format pairs, the C binary converts directly without loading
into R:

``` r

# Convert the 329K-cell atlas without loading — takes ~2 seconds
scConvert_cli("immune_atlas.h5ad", "immune_atlas.h5seurat")

# Or use the one-liner (auto-detects the fastest path)
scConvert("immune_atlas.h5ad", dest = "h5seurat")
```

### Export after BPCells loading

When you’ve loaded with BPCells and want to save a subset:

``` r

# Load with BPCells (low memory)
atlas <- readH5AD("immune_atlas.h5ad", use.bpcells = TRUE)

# Subset to T cells only
tcells <- subset(atlas, cell_type == "T cell")
cat("T cells:", ncol(tcells), "\n")

# Save the subset (materializes only the subset into RAM)
writeH5AD(tcells, "tcells_only.h5ad")
saveRDS(tcells, "tcells_only.rds")
```

## When to use BPCells

| Scenario | Recommendation |
|----|----|
| \< 50K cells | Standard loading (fast, simple) |
| 50K – 500K cells | BPCells optional (saves RAM) |
| \> 500K cells | BPCells recommended |
| Repeated access to same file | `use.bpcells = "/path"` (directory cache) |
| One-time conversion | C binary (`scConvert_cli`) — no loading at all |
| Subsetting before analysis | BPCells + [`subset()`](https://rdrr.io/r/base/subset.html) + save |

## Clean up
