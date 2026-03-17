# Converting Between Seurat and TileDB-SOMA

## What is TileDB-SOMA?

[TileDB-SOMA](https://github.com/single-cell-data/TileDB-SOMA) is a
cloud-native, array-based format for storing and querying single-cell
data at scale. It is the storage backend behind [CELLxGENE
Census](https://chanzuckerberg.github.io/cellxgene-census/), the largest
public atlas of single-cell data (61M+ cells across 900+ datasets).

Key properties of SOMA:

- **Columnar storage** backed by TileDB arrays – efficient slicing by
  cells or features without loading the full dataset
- **Cloud-native**: local paths, S3, GCS, and TileDB Cloud URIs
- **Query filtering**: select cells by metadata predicates at the
  storage level
- **Multi-measurement**: RNA, ATAC, and protein assays in one experiment

scConvert provides
[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md)
and
[`writeSOMA()`](https://mianaz.github.io/scConvert/reference/writeSOMA.md)
for Seurat interoperability, plus pair converters for all supported
formats. All SOMA functions require the
[tiledbsoma](https://github.com/single-cell-data/TileDB-SOMA) R package.

``` r

# Install tiledbsoma (required for SOMA support)
install.packages("tiledbsoma")

library(scConvert)
library(Seurat)
```

## Reading from SOMA

Use
[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md)
to load a SOMA experiment as a Seurat object:

``` r

# Read a local SOMA experiment
obj <- readSOMA("path/to/experiment")

# Read a specific measurement (default is "RNA")
obj <- readSOMA("path/to/experiment", measurement = "ATAC")
```

### Filtering cells and features

[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md)
supports server-side filtering via value filter queries, which is
particularly useful for large atlases:

``` r

# Filter by cell type and tissue
obj <- readSOMA(
  "path/to/experiment",
  obs_query = "cell_type == 'T cell' & tissue == 'blood' & disease == 'normal'",
  obs_column_names = c("cell_type", "tissue", "donor_id")
)
```

### What gets read

[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md)
maps the following SOMA components to Seurat:

| SOMA Source | Seurat Destination | Description |
|----|----|----|
| `obs` | `meta.data` | Cell metadata (factors, numerics, strings) |
| `ms[measurement]/var` | Feature metadata | Gene-level annotations |
| `ms[measurement]/X/data` | Default assay counts | Expression matrix (sparse) |
| `ms[measurement]/obsm/*` | `reductions` | PCA, UMAP, etc. |
| `ms[measurement]/obsp/*` | `graphs` | Neighbor graphs (SNN, distances) |
| Other measurements | Additional assays | Multi-modal data (ATAC, protein) |

The `soma_joinid` column (internal to TileDB-SOMA) is automatically
removed from cell and feature metadata.

## Writing to SOMA

Use
[`writeSOMA()`](https://mianaz.github.io/scConvert/reference/writeSOMA.md)
to save a Seurat object as a SOMA experiment:

``` r

# Write to a local SOMA experiment
writeSOMA(seurat_obj, uri = "output.soma")

# Overwrite an existing experiment
writeSOMA(seurat_obj, uri = "output.soma", overwrite = TRUE)
```

[`writeSOMA()`](https://mianaz.github.io/scConvert/reference/writeSOMA.md)
delegates to `tiledbsoma::write_soma()`, which converts all assays,
metadata, embeddings, and graphs. Each Seurat assay becomes a separate
SOMA measurement.

## Direct pair converters

scConvert provides direct conversion functions between SOMA and all
other supported formats. All route through a Seurat intermediate.

### SOMA and h5ad

``` r

# h5ad -> SOMA
H5ADToSOMA("data.h5ad", "data.soma")

# SOMA -> h5ad
SOMAToH5AD("data.soma", "output.h5ad")
```

### SOMA and h5Seurat

``` r

# h5Seurat -> SOMA
H5SeuratToSOMA("data.h5seurat", "data.soma")

# SOMA -> h5Seurat
SOMAToH5Seurat("data.soma", "output.h5seurat")
```

### Other format pairs

Converters also exist for Loom, Zarr, and h5mu – each follows the same
`{Format}ToSOMA()` / `SOMATo{Format}()` pattern:

``` r

LoomToSOMA("data.loom", "data.soma")    # Loom -> SOMA
ZarrToSOMA("data.zarr", "data.soma")    # Zarr -> SOMA
H5MUToSOMA("data.h5mu", "data.soma")    # h5mu -> SOMA
SOMAToLoom("data.soma", "output.loom")   # SOMA -> Loom
SOMAToZarr("data.soma", "output.zarr")   # SOMA -> Zarr
SOMAToH5MU("data.soma", "output.h5mu")  # SOMA -> h5mu
```

## CLI and dispatcher usage

Both
[`scConvert_cli()`](https://mianaz.github.io/scConvert/reference/scConvert_cli.md)
and the universal `scConvert()` dispatcher auto-detect SOMA from the
`.soma` extension:

``` r

# scConvert_cli auto-dispatches
scConvert_cli("data.h5ad", "data.soma")
scConvert_cli("data.soma", "output.h5ad")

# Universal dispatcher
scConvert("data.h5ad", dest = "data.soma", overwrite = TRUE)
scConvert("data.soma", dest = "output.h5seurat", overwrite = TRUE)
```

## CELLxGENE Census access

[CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/)
hosts 61M+ cells as a single SOMA experiment. Query it with
[`readSOMA()`](https://mianaz.github.io/scConvert/reference/readSOMA.md),
then convert the result to any format:

``` r

library(cellxgene.census)

census <- open_soma(census_version = "stable")
human_uri <- census$get("census_data")$get("homo_sapiens")$uri

# Read a subset directly as a Seurat object
tcells <- readSOMA(
  uri = human_uri,
  measurement = "RNA",
  obs_query = "cell_type == 'T cell' & tissue_general == 'blood'",
  obs_column_names = c("cell_type", "tissue", "donor_id")
)

# Save locally in any format
writeH5AD(tcells, "census_tcells.h5ad")
writeSOMA(tcells, "census_tcells.soma")
```

See the [cellxgene-census
documentation](https://chanzuckerberg.github.io/cellxgene-census/) for
more details on the Census API.

## Data mapping reference

| Seurat | SOMA | Direction | Notes |
|----|----|----|----|
| Default assay counts | `ms[RNA]/X/data` | Read/Write | Sparse matrix |
| `meta.data` | `obs` | Read/Write | Factors, numerics, strings |
| Feature metadata | `ms[RNA]/var` | Read/Write | Gene annotations |
| `Embeddings(, "pca")` | `ms[RNA]/obsm/X_pca` | Read/Write | Dense array |
| `Embeddings(, "umap")` | `ms[RNA]/obsm/X_umap` | Read/Write | Dense array |
| `Graphs(, "RNA_snn")` | `ms[RNA]/obsp/RNA_snn` | Read/Write | Sparse matrix |
| Additional assays | Other measurements | Read/Write | Multi-modal support |

## Limitations

- **Requires tiledbsoma**: install with
  `install.packages("tiledbsoma")`.
  [`writeSOMA()`](https://mianaz.github.io/scConvert/reference/writeSOMA.md)
  further requires `tiledbsoma::write_soma()` to be available.
- **No streaming conversion**: unlike h5ad/zarr converters, SOMA always
  routes through a Seurat intermediate (TileDB arrays differ
  fundamentally from HDF5).
- **No C binary path**: the compiled C CLI handles HDF5-to-HDF5 only.
- **scale.data** is not written to SOMA – re-run
  [`ScaleData()`](https://satijalab.org/seurat/reference/ScaleData.html)
  after loading.
- **Active identity**
  ([`Idents()`](https://mianaz.github.io/scConvert/reference/h5Seurat-bindings.md))
  is not stored – set manually after reading.

## See also

- [`vignette("convert-anndata")`](https://mianaz.github.io/scConvert/articles/convert-anndata.md)
  – Converting between Seurat and h5ad format
- [`vignette("cli-usage")`](https://mianaz.github.io/scConvert/articles/cli-usage.md)
  – Command-line interface for batch conversions
- [`vignette("convert-zarr")`](https://mianaz.github.io/scConvert/articles/convert-zarr.md)
  – Converting between Seurat and Zarr format

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
#> [25] rlang_1.1.7       fs_1.6.7          htmlwidgets_1.6.4
```
