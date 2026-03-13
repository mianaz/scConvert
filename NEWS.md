# scConvert 0.1.0

> **Release Date:** 2026-03-10

## Highlights

Initial public release of scConvert — a universal single-cell format converter for R.

### Universal Format Conversion
- Support for 7 formats: h5ad, h5Seurat, h5mu, Loom, Zarr, RDS, and SingleCellExperiment
- Hub architecture with 30+ conversion paths via `scConvert()`
- Direct HDF5 paths for h5ad/h5Seurat without intermediate loading

### Direct h5ad Loading
- `readH5AD()` for native h5ad-to-Seurat conversion without intermediate files
- Sparse (CSR/CSC) and dense matrix support
- Categorical metadata, dimensional reductions, neighbor graphs, and spatial data

### MuData (h5mu) Multimodal Support
- `readH5MU()` / `writeH5MU()` for multimodal single-cell data
- Automatic modality-to-assay name mapping (rna->RNA, prot->ADT, atac->ATAC)
- No MuDataSeurat or Python dependency required

### Zarr AnnData Support
- `readZarr()` and `writeZarr()` for Zarr-based AnnData stores (v2 format)
- Sparse CSR/CSC and dense matrix support
- Categorical metadata and dimensional reduction preservation

### Spatial Data (Visium)
- Bidirectional Visium spatial data conversion with image reconstruction
- Proper coordinate handling and scale factor preservation
- Compatible with scanpy/squidpy spatial analysis workflows

### C CLI Binary
- Standalone `scconvert` binary for h5ad/h5Seurat/h5mu conversions
- Streaming on-disk conversion without R or Python runtime
- Options: `--assay`, `--gzip`, `--overwrite`, `--quiet`

### BPCells On-Disk Loading
- `readH5AD(..., use.bpcells = TRUE)` for memory-efficient atlas-scale analysis
- Compatible with all Seurat analysis functions

### Seurat v5 Compatibility
- Full support for Seurat v5 Assay5 objects
- Proper handling of V5 layered data structure in conversions
