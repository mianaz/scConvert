
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scConvert <img src="man/figures/graphical_abstract.png" align="right" width="45%" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/scConvert)](https://CRAN.R-project.org/package=scConvert)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://github.com/mianaz/scConvert)
[![R-CMD-check](https://github.com/mianaz/scConvert/actions/workflows/r-cmd-check.yaml/badge.svg)](https://github.com/mianaz/scConvert/actions/workflows/r-cmd-check.yaml)
[![Codecov](https://codecov.io/gh/mianaz/scConvert/graph/badge.svg)](https://codecov.io/gh/mianaz/scConvert)
<!-- badges: end -->

**Universal single-cell format conversion for R. No Python required.**

scConvert converts between **9 single-cell data formats** through a hub
architecture using Seurat as a universal intermediate, providing **50+
conversion paths** from a single `scConvert()` interface.

## Supported Formats

| Format               | Extension   | Ecosystem          | Read | Write |
|----------------------|-------------|--------------------|:----:|:-----:|
| AnnData              | `.h5ad`     | scanpy / CELLxGENE | yes  |  yes  |
| h5Seurat             | `.h5Seurat` | Seurat             | yes  |  yes  |
| MuData               | `.h5mu`     | muon / multimodal  | yes  |  yes  |
| Loom                 | `.loom`     | loompy / HCA       | yes  |  yes  |
| Zarr                 | `.zarr`     | cloud AnnData      | yes  |  yes  |
| TileDB-SOMA          | `soma://`   | CELLxGENE Census   | yes  |  yes  |
| SpatialData          | `.zarr`     | scverse spatial    | yes  |  yes  |
| RDS                  | `.rds`      | R native           | yes  |  yes  |
| SingleCellExperiment | in-memory   | Bioconductor       | yes  |   –   |

## Installation

``` r
# Install from GitHub
devtools::install_github("mianaz/scConvert")
```

## Quick Start

``` r
library(scConvert)

# One-line conversion between any format pair
scConvert("data.h5ad", dest = "h5seurat")
scConvert("data.h5mu", dest = "rds")
scConvert(seurat_obj, dest = "output.h5ad")

# Direct h5ad loading with full metadata preservation
obj <- readH5AD("data.h5ad")

# Memory-efficient atlas-scale loading via BPCells
obj <- readH5AD("atlas_1M_cells.h5ad", use.bpcells = TRUE)

# Multimodal (CITE-seq, ATAC+RNA) support
obj <- readH5MU("citeseq.h5mu")   # auto-maps rna->RNA, prot->ADT
writeH5MU(obj, "output.h5mu")

# TileDB-SOMA (CELLxGENE Census)
obj <- readSOMA("soma://collection/measurement")
writeSOMA(obj, "output.soma")

# SpatialData zarr stores
obj <- readSpatialData("experiment.spatialdata.zarr")
writeSpatialData(obj, "output.spatialdata.zarr")
```

## Hub Architecture

Three conversion tiers provide optimal speed for each format pair:

                      ┌─── h5ad ───┐
                      │      ╎     │
                      ├─── h5Seurat│
                      │      ╎     ├── C CLI (fastest, streaming)
       ┌────────┐     ├─── h5mu ───┤
       │ Seurat │─────┤      ╎     │
       └────────┘     ├─── Loom ───┘
                      │
                      ├─── Zarr ········ R streaming (no Seurat)
                      │
                      ├─── TileDB-SOMA
                      │
                      ├─── SpatialData
                      │
                      ├─── RDS
                      │
                      └─── SCE (Bioconductor)

       ─── hub path (via Seurat)
       ╎╎╎ direct HDF5 streaming (no intermediate)
       ··· R streaming via temp h5seurat

- **C CLI**: All HDF5 pairs (h5ad, h5Seurat, h5mu, Loom) – streaming,
  constant memory
- **R streaming**: Direct format-to-format without materializing a
  Seurat object
- **Hub path**: Load → Seurat → Save for RDS, SCE, and cross-tier pairs

## Key Features

| Feature | Description |
|----|----|
| **No Python** | Pure R via hdf5r. No reticulate, no conda environments. |
| **Seurat v5** | Full Assay5 support with layered counts/data/scale.data. |
| **Spatial** | Visium roundtrip with image reconstruction and scale factors. |
| **SpatialData** | Read/write scverse SpatialData zarr stores with OME-NGFF images. |
| **Multimodal** | Native h5mu read/write for CITE-seq, ATAC+RNA, etc. |
| **TileDB-SOMA** | Read/write SOMA for CELLxGENE Census interoperability. |
| **BPCells** | On-disk matrix loading – 87% memory reduction at atlas scale. |
| **C CLI** | Standalone binary for streaming on-disk HDF5 conversion. |

## Performance

Benchmarked against 9 competing tools on datasets from 100 to 1,000,000
cells (Apple M4 Max, 48 GB RAM):

| Operation              | 100K cells | 500K cells | Comparison                    |
|------------------------|------------|------------|-------------------------------|
| Read h5ad              | 2.84 s     | 9.69 s     | 3–9x faster than alternatives |
| Write h5ad             | 21.3 s     | –          | Comparable to anndataR        |
| BPCells load           | 2.55 s     | 6.24 s     | 87% less memory               |
| CLI (h5ad to h5seurat) | 0.12 s     | 0.59 s     | Streaming, constant memory    |

Benchmarked on Apple M4 Max, 48 GB RAM.

## Command-Line Interface

A standalone C binary for HDF5-based conversions (h5ad, h5Seurat, h5mu,
Loom) without R or Python:

``` bash
# Build (requires HDF5 headers)
cd src && make

# Convert between any HDF5 format pair
scconvert data.h5ad data.h5seurat
scconvert data.h5seurat data.h5ad --assay RNA
scconvert multimodal.h5mu multimodal.h5seurat
scconvert data.h5ad data.loom
scconvert data.loom data.h5seurat
```

Options: `--assay <name>`, `--gzip <0-9>`, `--overwrite`, `--quiet`,
`--version`, `--help`

## Documentation

- [Conversions: h5Seurat and
  AnnData](https://mianaz.github.io/scConvert/articles/convert-anndata.html)
- [Direct H5AD
  Loading](https://mianaz.github.io/scConvert/articles/direct-h5ad-loading.html)
- [Multimodal
  H5MU](https://mianaz.github.io/scConvert/articles/multimodal-h5mu.html)
- [Loom
  Format](https://mianaz.github.io/scConvert/articles/convert-loom.html)
- [Zarr
  Format](https://mianaz.github.io/scConvert/articles/convert-zarr.html)
- [Spatial
  Technologies](https://mianaz.github.io/scConvert/articles/spatial-technologies.html)
- [CLI
  Usage](https://mianaz.github.io/scConvert/articles/cli-usage.html)
- [h5Seurat
  Specification](https://mianaz.github.io/scConvert/articles/h5Seurat-spec.html)

## Citation

If you use scConvert in your research, please cite:

> Zeng Z. scConvert: a pure R, universal converter for single-cell data
> formats. *bioRxiv* (2026).

## License

GPL-3
