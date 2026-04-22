# Fast Seurat object constructor (bypass CreateSeuratObject validation)

Constructs a Seurat object directly from a pre-validated sparse matrix,
skipping the expensive validation in CreateSeuratObject (cell/feature
filtering, name checks, etc.). Use only when data is known-good (e.g.,
freshly read from a validated HDF5 file).

## Usage

``` r
.fast_create_seurat(counts, assay.name = "RNA")
```

## Arguments

- counts:

  dgCMatrix with rownames (features) and colnames (cells)

- assay.name:

  Name for the assay

## Value

Seurat object
