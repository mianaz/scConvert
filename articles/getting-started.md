# Getting Started with scConvert

``` r

suppressPackageStartupMessages({
  library(scConvert)
  library(Seurat)
})
```

scConvert converts between single-cell data formats entirely in R, with
no Python dependency. This vignette walks through installation, basic
conversion, and file I/O in under five minutes.

## Installation

Install scConvert from GitHub:

``` r

remotes::install_github("mianaz/scConvert")
```

## Quick conversion with scConvert()

The
[`scConvert()`](https://mianaz.github.io/scConvert/reference/scConvert.md)
function is a universal dispatcher: give it a source and a destination,
and it picks the fastest conversion path automatically. File formats are
detected from extensions.

``` r

h5ad_file <- system.file("testdata", "pbmc_small.h5ad", package = "scConvert")

h5seurat_out <- file.path(tempdir(), "pbmc.h5seurat")
scConvert(h5ad_file, dest = h5seurat_out, overwrite = TRUE)

zarr_out <- file.path(tempdir(), "pbmc.zarr")
scConvert(h5ad_file, dest = zarr_out, overwrite = TRUE)
```

[`scConvert()`](https://mianaz.github.io/scConvert/reference/scConvert.md)
works for any supported format pair – h5ad, h5Seurat, Zarr, Loom, RDS,
and more.

## Loading files into Seurat

When you need to work with the data in R, use the format-specific
readers. Each returns a standard Seurat object.

### From h5ad

``` r

obj <- readH5AD(h5ad_file)
obj
#> An object of class Seurat 
#> 2000 features across 214 samples within 1 assay 
#> Active assay: RNA (2000 features, 2000 variable features)
#>  2 layers present: counts, data
#>  2 dimensional reductions calculated: pca, umap
```

``` r

dim(obj)
#> [1] 2000  214
head(obj[[]], n = 4)
#>                orig.ident nCount_RNA nFeature_RNA seurat_annotations percent.mt
#> AACCAGTGATACCG     pbmc3k       1539          335       FCGR3A+ Mono   2.786458
#> AAGATTACCGCCTT     pbmc3k        898          284                 DC   1.763553
#> AAGCAAGAGCTTAG     pbmc3k        668          201                 NK   2.214022
#> AAGCCATGAACTGC     pbmc3k       2329          435                 DC   1.401472
#>                RNA_snn_res.0.5 seurat_clusters
#> AACCAGTGATACCG               5               5
#> AAGATTACCGCCTT               7               7
#> AAGCAAGAGCTTAG               6               6
#> AAGCCATGAACTGC               7               7
```

### From h5Seurat

``` r

obj2 <- readH5Seurat(h5seurat_out)
obj2
#> An object of class Seurat 
#> 2000 features across 214 samples within 1 assay 
#> Active assay: RNA (2000 features, 2000 variable features)
#>  2 layers present: counts, data
#>  2 dimensional reductions calculated: pca, umap
```

### From Zarr

``` r

obj3 <- readZarr(zarr_out)
obj3
#> An object of class Seurat 
#> 2000 features across 214 samples within 1 assay 
#> Active assay: RNA (2000 features, 2000 variable features)
#>  2 layers present: counts, data
#>  2 dimensional reductions calculated: pca, umap
```

## Writing files from Seurat

Starting from any Seurat object, write to the format your collaborators
or downstream tools need.

``` r

h5ad_out <- file.path(tempdir(), "output.h5ad")
writeH5AD(obj, h5ad_out, verbose = FALSE)

h5s_out <- file.path(tempdir(), "output.h5seurat")
writeH5Seurat(obj, h5s_out, overwrite = TRUE, verbose = FALSE)

zarr_out2 <- file.path(tempdir(), "output.zarr")
writeZarr(obj, zarr_out2, verbose = FALSE)
```

``` r

sizes <- data.frame(
  Format = c("h5ad", "h5Seurat", "Zarr"),
  Size_MB = round(c(
    file.size(h5ad_out),
    file.size(h5s_out),
    sum(file.info(list.files(zarr_out2, recursive = TRUE, full.names = TRUE))$size)
  ) / 1024^2, 2)
)
knitr::kable(sizes, col.names = c("Format", "Size (MB)"))
```

| Format   | Size (MB) |
|:---------|----------:|
| h5ad     |      0.54 |
| h5Seurat |      0.59 |
| Zarr     |      0.34 |

## Supported formats

| Format               | Extension   | Ecosystem         | Read | Write |
|----------------------|-------------|-------------------|:----:|:-----:|
| AnnData              | `.h5ad`     | scanpy, CELLxGENE | yes  |  yes  |
| h5Seurat             | `.h5seurat` | Seurat            | yes  |  yes  |
| MuData               | `.h5mu`     | muon (multimodal) | yes  |  yes  |
| Loom                 | `.loom`     | loompy, HCA       | yes  |  yes  |
| Zarr                 | `.zarr`     | cloud AnnData     | yes  |  yes  |
| TileDB-SOMA          | `soma://`   | CELLxGENE Census  | yes  |  yes  |
| SpatialData          | `.zarr`     | scverse spatial   | yes  |  yes  |
| RDS                  | `.rds`      | R native          | yes  |  yes  |
| SingleCellExperiment | in-memory   | Bioconductor      | yes  |   –   |

**When to use which format:**

- **h5ad** – sharing with Python users or submitting to CELLxGENE.
- **h5Seurat** – archiving a Seurat object with selective loading
  support.
- **Zarr** – cloud-friendly, chunk-based access (e.g., S3 / GCS).
- **Loom** – interoperability with loompy, velocyto, or legacy
  pipelines.
- **RDS** – quick local save/load within R (no HDF5 dependency).

## Verifying a conversion

After converting, a quick sanity check confirms that dimensions and cell
identifiers are preserved.

``` r

original <- readH5AD(h5ad_file)
converted <- readH5Seurat(h5seurat_out)

stopifnot(identical(dim(original), dim(converted)))
stopifnot(identical(sort(colnames(original)), sort(colnames(converted))))

cat("Dimensions match:", paste(dim(original), collapse = " x "), "\n")
#> Dimensions match: 2000 x 214
cat("Cell names match:", length(colnames(original)), "cells verified\n")
#> Cell names match: 214 cells verified
```

## Next steps

- **[Convert Between Seurat and
  AnnData](https://mianaz.github.io/scConvert/articles/convert-anndata.md)**
  – layer mapping, Python interop, and round-trip verification.
- **[In-Memory vs On-Disk
  Conversion](https://mianaz.github.io/scConvert/articles/conversion-modes.md)**
  – hub path, streaming, and the C binary for large datasets.
- **[Direct H5AD
  Loading](https://mianaz.github.io/scConvert/articles/h5Seurat-load.md)**
  – selective loading and BPCells for atlas-scale data.
- **[Multimodal
  H5MU](https://mianaz.github.io/scConvert/articles/multimodal-h5mu.md)**
  – CITE-seq and ATAC+RNA via MuData.
- **[Zarr
  Format](https://mianaz.github.io/scConvert/articles/convert-zarr.md)**
  – cloud-native storage and streaming converters.
- **[Spatial
  Technologies](https://mianaz.github.io/scConvert/articles/spatial-technologies.md)**
  – Visium, MERFISH, and SpatialData support.
- **[CLI
  Usage](https://mianaz.github.io/scConvert/articles/cli-usage.md)** –
  the standalone C binary for batch conversion without R.

## Clean up
