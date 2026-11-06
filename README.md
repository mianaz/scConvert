# scConvert <img src="man/figures/graphical_abstract.png" alt="scConvert graphical abstract showing format conversion paths" align="right" width="45%" />

<!-- badges: start -->

[![GitHub
release](https://img.shields.io/github/v/release/mianaz/scConvert?label=release)](https://github.com/mianaz/scConvert/releases)
[![Lifecycle:
experimental](https://lifecycle.r-lib.org/articles/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/mianaz/scConvert/actions/workflows/r-cmd-check.yaml/badge.svg)](https://github.com/mianaz/scConvert/actions/workflows/r-cmd-check.yaml)
[![pkgdown](https://github.com/mianaz/scConvert/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/mianaz/scConvert/actions/workflows/pkgdown.yaml)
[![Codecov](https://codecov.io/gh/mianaz/scConvert/graph/badge.svg)](https://codecov.io/gh/mianaz/scConvert)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
<!-- badges: end -->

**Load any single-cell format as a Seurat object. Write a Seurat object as any format. Convert any two formats in your terminal.**

scConvert bridges **11 single-cell data formats** through a hub architecture
using Seurat as a universal intermediate, providing **50+ conversion paths**
from a single `scConvert()` call in R or one `scconvert` command in the shell.

## Supported Formats

| Format               | Extension        | Ecosystem          | Read | Write |
|----------------------|------------------|--------------------|:----:|:-----:|
| AnnData              | `.h5ad`          | scanpy / CELLxGENE | yes  |  yes  |
| Spatial AnnData      | `.h5ad`          | Visium / squidpy   | yes  |  yes  |
| h5Seurat             | `.h5Seurat`      | Seurat             | yes  |  yes  |
| MuData               | `.h5mu`          | muon / multimodal  | yes  |  yes  |
| Loom                 | `.loom`          | loompy / HCA       | yes  |  yes  |
| Zarr                 | `.zarr`          | cloud AnnData      | yes  |  yes  |
| TileDB-SOMA          | `soma://`        | CELLxGENE Census   | yes  |  yes  |
| SpatialData          | `.zarr`          | scverse spatial    | yes  |  yes  |
| RDS                  | `.rds`           | R native           | yes  |  yes  |
| SingleCellExperiment | in-memory        | Bioconductor       | yes  |   –   |
| Stereo-seq GEF       | `.gef` / `.cellbin.gef` | STOmics   | yes  |   –   |

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
| **Large data** | Streaming conversion between h5ad, h5Seurat, h5mu and Loom that never loads the expression matrices (`stream = TRUE`), and parallel batch conversion (`scConvertBatch()`). |
| **BPCells** | On-disk matrix loading – 87% memory reduction at atlas scale. |
| **C CLI** | Standalone binary for fast HDF5 conversion (opt-in; see below). |
| **Version compatibility** | Reads every anndata layout from 0.7 to 0.13 (nullable strings, `null`, zarr v3 sharding + zstd); `upgradeH5AD()`/`downgradeH5AD()` rewrite h5ad files between layouts and `upgradeSeurat()`/`downgradeSeurat()` convert objects between Seurat v3/v4 and v5 losslessly. |

## Performance

Synthetic sparse h5ad (20K genes, 5% density), median of 3 runs:

| Operation              | 100K cells | 500K cells | Comparison                    |
|------------------------|------------|------------|-------------------------------|
| Read h5ad              | 2.56 s     | 12.68 s    | Native C reader + Seurat      |
| Write h5ad (gzip=0)    | 0.61 s     | 3.29 s     | C writer, 17× faster          |
| Write h5ad (gzip=1)    | 6.03 s     | 30.66 s    | C writer, compressed          |
| Write h5ad (old R)     | 10.49 s    | 52.32 s    | R/hdf5r, gzip=4              |
| CLI (h5ad ↔ h5seurat)  | 0.19 s     | 0.63 s     | Streaming, constant memory    |

The C writer writes h5ad directly via `.Call()`, using the zero-copy
CSC↔CSR reinterpretation (no sparse transpose needed). Use `gzip = 0` for
local/temporary files (6× faster), `gzip = 1` for a good speed/size balance
(252 MB vs 1.2 GB at 100K cells). BPCells on-disk loading (`use.bpcells`)
reduces memory by ~87% at atlas scale. Benchmarked on Apple M4 Max, 128 GB
RAM.

## Large Datasets

For files that do not fit comfortably in memory, convert without loading the
object:

``` r
# any pair of h5ad, h5Seurat, h5mu and Loom; expression matrices copied on disk
scConvert("atlas.h5ad", dest = "atlas.h5seurat", stream = TRUE)
scConvert("atlas.h5seurat", dest = "atlas.h5ad", stream = TRUE)
scConvert("citeseq.h5mu", dest = "citeseq.h5seurat", stream = TRUE)

# many files at once, one worker process per file
files <- list.files("samples", pattern = "\\.h5ad$", full.names = TRUE)
scConvertBatch(files, dest = "h5seurat", workers = 4, stream = TRUE)
```

Streaming reads only the small parts of the file (metadata, embeddings,
graphs, images, `uns`) and converts them with the regular reader and writer;
the expression matrices (h5ad `X`, `layers` and `raw`, h5Seurat layers, h5mu
modalities, the Loom matrix and layers) are copied between the files by HDF5,
or in 64 MB blocks (`options(scConvert.stream_chunk_bytes)`) when they must be
transposed or densified, or features subset or reordered. The output is identical to an
in-memory conversion. `stream = "auto"` (the default) streams files larger
than 2 GB (`options(scConvert.stream_threshold)`).

Synthetic h5ad, 100,000 cells x 5,000 genes, 100 million non-zeros (Linux,
4 cores):

| h5ad → h5Seurat | Peak memory | Time |
|-----------------|------------:|-----:|
| in memory       | 6.5 GB      | 101 s |
| `stream = TRUE` | 1.2 GB      | 34 s |

Streaming memory depends on the number of cells (metadata), not on the size
of the matrices: 60 million non-zeros over 200,000 cells also peak at 1.2 GB,
of which 0.4 GB is R with Seurat loaded. Converting four 50,000-cell files
with `scConvertBatch(workers = 4)` took 17 s instead of 47 s.

## Command-Line Interface

A standalone C binary for HDF5-based conversions (h5ad, h5Seurat, h5mu,
Loom) without R or Python:

``` bash
# Build (requires HDF5 headers)
cd src && make -f Makefile.cli

# Optional: install into the R package tree so `library(scConvert)` picks it up
make -f Makefile.cli install-bin   # copies to ../inst/bin/scconvert

# Convert between any HDF5 format pair
./scconvert data.h5ad data.h5seurat
./scconvert data.h5seurat data.h5ad --assay RNA
./scconvert multimodal.h5mu multimodal.h5seurat
./scconvert data.h5ad data.loom
./scconvert data.loom data.h5seurat
```

Options: `--assay <name>`, `--gzip <0-9>`, `--overwrite`, `--quiet`,
`--version`, `--help`

The CLI is fast but does not carry everything the R paths do (feature
metadata, variable features and spatial images are dropped), so `scConvert()`
only uses it with `options(scConvert.use_cli = TRUE)`. For large files with
full fidelity, use `stream = TRUE`.

## Documentation

- [Conversions: h5Seurat and
  AnnData](https://mianaz.github.io/scConvert/articles/convert-anndata.html)
- [Getting
  Started](https://mianaz.github.io/scConvert/articles/getting-started.html)
- [Multimodal
  H5MU](https://mianaz.github.io/scConvert/articles/multimodal-h5mu.html)
- [Loom
  Format](https://mianaz.github.io/scConvert/articles/convert-loom.html)
- [Zarr
  Format](https://mianaz.github.io/scConvert/articles/convert-zarr.html)
- [Spatial: Visium](https://mianaz.github.io/scConvert/articles/spatial-visium.html)
  and [CosMx](https://mianaz.github.io/scConvert/articles/spatial-cosmx.html)
- [CLI
  Usage](https://mianaz.github.io/scConvert/articles/cli-usage.html)
- [Atlas-scale
  conversion](https://mianaz.github.io/scConvert/articles/atlas-scale-cli.html)

## Citation

If you use scConvert in your research, please cite:

> Zeng Z. (2026). scConvert: a pure R, universal converter for single-cell data
> formats. GitHub. https://github.com/mianaz/scConvert

## License

scConvert is released under the [MIT License](LICENSE.md), © 2026 Ziyu Zeng.
Versions before 0.3.0 were distributed under GPL-3 because they included code
from SeuratDisk; that code has been replaced by an independent
implementation. scConvert reads and writes the h5Seurat file layout used by
SeuratDisk.
