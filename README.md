
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scConvert <img src="man/figures/graphical_abstract.png" align="right" width="45%" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/scConvert)](https://CRAN.R-project.org/package=scConvert)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://github.com/mianaz/scConvert)
[![R-CMD-check](https://github.com/mianaz/scConvert/actions/workflows/r-cmd-check.yaml/badge.svg)](https://github.com/mianaz/scConvert/actions/workflows/r-cmd-check.yaml)
[![Codecov](https://codecov.io/gh/mianaz/scConvert/graph/badge.svg)](https://codecov.io/gh/mianaz/scConvert)
<!-- badges: end -->

**Universal single-cell format conversion for R. No Python required.**

scConvert converts between **7 single-cell data formats** through a hub
architecture using Seurat as a universal intermediate, providing **30+
conversion paths** from a single `scConvert()` interface. Direct streaming
converters bypass the hub for h5ad, h5Seurat, and Zarr pairs.

## Supported Formats

| Format | Extension | Ecosystem | Read | Write |
|--------|-----------|-----------|:----:|:-----:|
| AnnData | `.h5ad` | scanpy / CELLxGENE | yes | yes |
| h5Seurat | `.h5Seurat` | Seurat | yes | yes |
| MuData | `.h5mu` | muon / multimodal | yes | yes |
| Loom | `.loom` | loompy / HCA | yes | yes |
| Zarr | `.zarr` | cloud AnnData | yes | yes |
| RDS | `.rds` | R native | yes | yes |
| SingleCellExperiment | in-memory | Bioconductor | yes | -- |

## Installation

```r
# Install from GitHub
devtools::install_github("mianaz/scConvert")
```

## Quick Start

```r
library(scConvert)

# One-line conversion between any format pair
scConvert("data.h5ad", dest = "h5seurat")
scConvert("data.h5mu", dest = "rds")
scConvert(seurat_obj, dest = "output.h5ad")

# Direct h5ad loading with full metadata preservation
obj <- LoadH5AD("data.h5ad")

# Memory-efficient atlas-scale loading via BPCells
obj <- LoadH5AD("atlas_1M_cells.h5ad", use.bpcells = TRUE)

# Multimodal (CITE-seq, ATAC+RNA) support
obj <- LoadH5MU("citeseq.h5mu")   # auto-maps rna->RNA, prot->ADT
SaveH5MU(obj, "output.h5mu")

# Direct format converters (streaming, no Seurat intermediate)
H5ADToZarr("data.h5ad", "data.zarr")
ZarrToH5AD("data.zarr", "output.h5ad")
H5SeuratToZarr("data.h5seurat", "data.zarr")
ZarrToH5Seurat("data.zarr", "data.h5seurat")
```

## Architecture

All conversions route through Seurat as a universal intermediate.
Direct paths bypass the hub for performance-critical format pairs.

```
             ┌─── h5ad (AnnData) ───┐
             │        ╎             │
             ├─── h5Seurat          ├── Zarr
             │        ╎             │    ┊
  ┌────────┐ ├─── h5mu (MuData)    │
  │ Seurat │─┤                     │
  └────────┘ ├─── Loom             │
             │                     │
             ├─── RDS              │
             │                     │
             └─── SCE (Bioc)       │

  ─── hub path (via Seurat)
  ╎╎╎ direct HDF5 path (C binary)
  ┊┊┊ streaming path (field-by-field)
```

Three conversion tiers:

| Tier | Formats | Method |
|------|---------|--------|
| **C binary** | h5ad / h5Seurat / h5mu | On-disk HDF5 remapping, constant memory |
| **Streaming** | h5ad / h5Seurat ↔ Zarr | Field-by-field copy, no Seurat intermediate |
| **Hub** | All other pairs | Load → Seurat → Save |

## Key Features

| Feature | Description |
|---------|-------------|
| **No Python** | Pure R via hdf5r. No reticulate, no conda environments. |
| **Seurat v5** | Full Assay5 support with layered counts/data/scale.data. |
| **Spatial** | Visium roundtrip with image reconstruction and scale factors. |
| **Multimodal** | Native h5mu read/write for CITE-seq, ATAC+RNA, etc. |
| **BPCells** | On-disk matrix loading -- 87% memory reduction at atlas scale. |
| **Streaming** | Direct h5ad/h5Seurat/Zarr conversion without Seurat intermediate. |
| **C CLI** | Standalone binary for on-disk HDF5 conversion. |

## Direct Converters

Named converters provide direct format-to-format conversion with streaming
enabled by default. These avoid the Seurat hub entirely, preserving all data
bit-for-bit:

```r
# h5ad <-> Zarr (same AnnData spec, different storage backend)
H5ADToZarr("data.h5ad", "data.zarr")
ZarrToH5AD("data.zarr", "data.h5ad")

# h5Seurat <-> Zarr (structural translation: factor encoding, matrix layout)
H5SeuratToZarr("data.h5seurat", "data.zarr")
ZarrToH5Seurat("data.zarr", "data.h5seurat")

# h5ad <-> h5Seurat (direct HDF5 path)
H5ADToH5Seurat("data.h5ad", "data.h5seurat")
H5SeuratToH5AD("data.h5seurat", "data.h5ad")

# Fall back to Seurat hub if needed
H5ADToZarr("data.h5ad", "data.zarr", stream = FALSE)
```

## Performance

Benchmarked against 9 competing tools on datasets from 100 to 1,000,000
cells (Apple M4 Max, 48 GB RAM):

| Operation | 100K cells | 500K cells | Comparison |
|-----------|-----------|-----------|------------|
| Read h5ad | 2.84 s | 9.69 s | 3--9x faster than alternatives |
| Write h5ad | 21.3 s | -- | Comparable to anndataR |
| BPCells load | 2.55 s | 6.24 s | 87% less memory |
| CLI (h5ad to h5seurat) | 0.12 s | 0.59 s | Streaming, constant memory |

## Command-Line Interface

A standalone C binary for HDF5-based conversions without R or Python:

```bash
# Build (requires HDF5 headers)
cd src && make

# Convert between HDF5 formats
scconvert data.h5ad data.h5seurat
scconvert data.h5seurat data.h5ad --assay RNA
scconvert multimodal.h5mu multimodal.h5seurat
scconvert data.h5ad data.h5seurat --gzip 6 --overwrite
```

The R-level `scConvert_cli()` extends this to all format pairs, automatically
selecting the fastest available path (C binary → streaming → hub):

```r
scConvert_cli("data.h5ad", "data.zarr")       # streaming
scConvert_cli("data.h5seurat", "data.zarr")   # streaming
scConvert_cli("data.h5ad", "data.h5seurat")   # C binary
scConvert_cli("data.rds", "data.zarr")        # hub
```

Options: `--assay <name>`, `--gzip <0-9>`, `--overwrite`, `--quiet`, `--version`, `--help`

## Documentation

- [Conversions: h5Seurat and AnnData](https://mianaz.github.io/scConvert/articles/convert-anndata.html)
- [Conversions: Zarr Format](https://mianaz.github.io/scConvert/articles/convert-zarr.html)
- [Multimodal H5MU](https://mianaz.github.io/scConvert/articles/multimodal-h5mu.html)
- [Loom Format](https://mianaz.github.io/scConvert/articles/convert-loom.html)
- [CLI Usage](https://mianaz.github.io/scConvert/articles/cli-usage.html)
- [h5Seurat Specification](https://mianaz.github.io/scConvert/articles/h5Seurat-spec.html)
- [Data Mapping Reference](https://mianaz.github.io/scConvert/articles/data-mapping-reference.html)

## Citation

If you use scConvert in your research, please cite:

> Zeng Z. scConvert: a pure R, universal converter for single-cell data formats. *bioRxiv* (2026).

## License

GPL-3
