# Conversions: h5Seurat and AnnData

This vignette showcases how to convert between `Seurat` objects and
AnnData files via h5Seurat files. This allows interoperability between
Seurat and [Scanpy](https://scanpy.readthedocs.io/).

``` r

library(Seurat)
library(scConvert)
library(ggplot2)
library(patchwork)
```

## Converting from Seurat to AnnData via h5Seurat

To demonstrate conversion from a `Seurat` object to an AnnData file,
we’ll use the `pbmc3k.final` dataset from SeuratData - a processed PBMC
dataset with clustering and UMAP. If SeuratData is not available, you
can download an h5ad directly from [CZ CellxGene
Discover](https://cellxgene.cziscience.com/collections):

``` r

# Alternative: download from CellxGene (Tabula Microcebus lung, ~1.7K cells)
cxg_url <- "https://datasets.cellxgene.cziscience.com/582ffa61-1c54-4492-8763-2ecfa9efa070.h5ad"
download.file(cxg_url, destfile = "cxg_demo.h5ad", mode = "wb")
pbmc <- readH5AD("cxg_demo.h5ad")
```

``` r

library(SeuratData)
if (!"pbmc3k.final" %in% rownames(InstalledData())) {
  InstallData("pbmc3k")
}

data("pbmc3k.final", package = "pbmc3k.SeuratData")
pbmc <- UpdateSeuratObject(pbmc3k.final)
pbmc
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 2000 variable features)
#>  3 layers present: counts, data, scale.data
#>  2 dimensional reductions calculated: pca, umap
```

This is a fully processed Seurat object with clustering and dimensional
reductions:

``` r

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

![](convert-anndata_files/figure-html/plot_pbmc-1.png)

``` r

FeaturePlot(pbmc, features = "CD14", pt.size = 0.5)
```

![](convert-anndata_files/figure-html/plot_pbmc_cd14-1.png)

Converting the `Seurat` object to an AnnData file is a two-step process:

1.  Save the `Seurat` object as an h5Seurat file using
    [`writeH5Seurat()`](https://mianaz.github.io/scConvert/reference/writeH5Seurat.md)
2.  Convert to AnnData using `scConvert()`

``` r

cat("Seurat layers:", paste(Layers(pbmc), collapse = ", "), "\n")
#> Seurat layers: counts, data, scale.data
writeH5Seurat(pbmc, filename = "pbmc3k.h5Seurat", overwrite = TRUE)
scConvert("pbmc3k.h5Seurat", dest = "h5ad", overwrite = TRUE)
```

We can view the AnnData file in Scanpy:

``` python
import scanpy as sc
adata = sc.read_h5ad("pbmc3k.h5ad")
print(adata)
#> AnnData object with n_obs × n_vars = 2638 × 13714
#>     obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'seurat_annotations', 'percent.mt', 'RNA_snn_res.0.5', 'seurat_clusters'
#>     var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable'
#>     uns: 'n_variable_features', 'neighbors', 'seurat'
#>     obsm: 'X_pca', 'X_umap'
#>     varm: 'PCs'
#>     obsp: 'connectivities', 'distances'
```

And visualize with cluster annotations:

``` python
sc.pl.umap(adata, color='seurat_annotations', legend_loc='on data', legend_fontsize=8)
```

![](convert-anndata_files/figure-html/plot_adata-1.png)

``` python
sc.pl.umap(adata, color='CD14', use_raw=False)
```

![](convert-anndata_files/figure-html/plot_adata_cd14-3.png)

The conversion preserves expression patterns - CD14 shows consistent
distribution in both tools.

## Direct Loading with readH5AD

In addition to the two-step `scConvert()` +
[`readH5Seurat()`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md)
workflow shown above, scConvert provides
[`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
for loading h5ad files directly into Seurat objects without creating an
intermediate h5Seurat file.

``` r

# One-step: load h5ad directly into a Seurat object
pbmc_direct <- readH5AD("pbmc3k.h5ad", verbose = TRUE)

# Compare with the two-step approach
cat("Direct cells:", ncol(pbmc_direct), "| Two-step cells:", ncol(pbmc), "\n")
#> Direct cells: 2638 | Two-step cells: 2638
cat("Direct reductions:", paste(names(pbmc_direct@reductions), collapse = ", "), "\n")
#> Direct reductions: pca, umap
```

[`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md)
reads the following from h5ad files:

| h5ad Location | Seurat Destination | Description |
|----|----|----|
| `X` | Default assay `data` layer | Expression matrix (sparse or dense) |
| `raw/X` | `counts` layer | Raw counts if present |
| `layers/*` | Additional layers | Named layers mapped to Seurat slots |
| `obs` | `meta.data` | Cell metadata (categorical preserved as factors) |
| `var` | Feature metadata | Gene-level annotations |
| `var['highly_variable']` | [`VariableFeatures()`](https://satijalab.github.io/seurat-object/reference/VariableFeatures.html) | Variable feature selection |
| `obsm/X_umap` | `reductions$umap` | UMAP coordinates |
| `obsm/X_pca` | `reductions$pca` | PCA embeddings |
| `obsm/X_tsne` | `reductions$tsne` | tSNE coordinates |
| `obsm/spatial` | Spatial coordinates | Via [`H5ADSpatialToSeurat()`](https://mianaz.github.io/scConvert/reference/H5ADSpatialToSeurat.md) |
| `obsp/connectivities` | `graphs$RNA_snn` | SNN graph |
| `obsp/distances` | `graphs$RNA_nn` | Distance graph |
| `uns/*` | `misc` | Unstructured annotations |

**When to use each approach:**

| Scenario | Recommended |
|----|----|
| Quick exploration of an h5ad file | [`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md) |
| Round-trip editing (load, modify, re-export) | `scConvert()` + [`readH5Seurat()`](https://mianaz.github.io/scConvert/reference/readH5Seurat.md) |
| Need h5Seurat for other tools | `scConvert()` |
| Loading scanpy-processed data for Seurat analysis | [`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md) |
| Working with spatial h5ad from CellxGene | [`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md) |

## Converting from AnnData to Seurat via h5Seurat

To demonstrate conversion from AnnData to Seurat, we’ll use a colorectal
cancer sample from [CellxGene](https://cellxgene.cziscience.com).

``` r

h5ad_path <- system.file("testdata", "crc_sample.h5ad", package = "scConvert")
if (file.exists(h5ad_path)) {
  file.copy(h5ad_path, "crc_sample.h5ad", overwrite = TRUE)
} else {
  download.file(
    "https://datasets.cellxgene.cziscience.com/91cf9a95-0b9a-4ece-b8eb-7b9e3409a0d3.h5ad",
    "crc_sample.h5ad", mode = "wb"
  )
}
```

View the h5ad file in Scanpy:

``` python
import scanpy as sc
adata_crc = sc.read_h5ad("crc_sample.h5ad")
print(adata_crc)
#> AnnData object with n_obs × n_vars = 935 × 25344
#>     obs: 'total_counts', 'log1p_total_counts', 'Sample ID', 'PhenoGraph_clusters', 'Patient', 'Primary Site', 'Sample Type', 'Site', 'DC 1', 'DC 2', 'DC 3', 'DC 4', 'Module Absorptive Intestine Score', 'Module EMT Score', 'Module Injury Repair Score', 'Module Squamous Score', 'Module Neuroendocrine Score', 'Module Endoderm Development Score', 'Module Tumor ISC-like Score', 'Module Secretory Intestine Score', 'Module Intestine Score', 'palantir_pseudotime', 'palantir_neuroendocrine_branch_probability', 'palantir_squamous_branch_probability', 'Fetal, Conserved', 'Module Osteoblast Score', 'Treatment', 'donor_id', 'development_stage_ontology_term_id', 'sex_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'disease_ontology_term_id', 'tissue_type', 'tissue_ontology_term_id', 'cell_type_ontology_term_id', 'assay_ontology_term_id', 'suspension_type', 'Tumor Status', 'is_primary_data', 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
#>     var: 'total_counts', 'highly_variable', 'gene', 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
#>     uns: 'citation', 'default_embedding', 'neighbors', 'organism', 'organism_ontology_term_id', 'schema_reference', 'schema_version', 'title'
#>     obsm: 'X_umap'
```

CellxGene datasets use Ensembl IDs as `var_names` by default. The
`feature_name` column contains gene symbols. During conversion,
scConvert automatically uses gene symbols when `feature_name` is
available.

> **Note on data layers**: During conversion, AnnData’s `X` matrix
> (typically log-normalized data indicated by `log1p_total_counts` in
> obs) is stored in Seurat’s `data` layer. If `raw/X` exists, it becomes
> the `counts` layer.

Visualize with scanpy before conversion:

``` python
sc.pl.umap(adata_crc, color='tissue')
```

![](convert-anndata_files/figure-html/plot_crc_tissue-1.png)

``` python
import pandas as pd

sc.pp.normalize_total(adata_crc, target_sum=1e4)
sc.pp.log1p(adata_crc)

# Set gene symbols as var_names (force string conversion from Categorical)
adata_crc.var_names = pd.Index(adata_crc.var['feature_name'].astype(str).values)
adata_crc.var_names_make_unique()

adata_crc.write_h5ad("crc_normalized.h5ad")

sc.pl.umap(adata_crc, color='EPCAM', use_raw=False, title='EPCAM')
```

![](convert-anndata_files/figure-html/normalize_crc-3.png)

Convert to Seurat:

``` r

scConvert("crc_normalized.h5ad", dest = "h5seurat", overwrite = TRUE)
crc <- readH5Seurat("crc_normalized.h5seurat")

# Verify layer mapping: X -> data (log-normalized), raw/X -> counts (if exists)
cat("Layers:", paste(Layers(crc), collapse = ", "), "\n")
#> Layers: counts, data
cat("Data layer range:", round(range(GetAssayData(crc, layer = "data")[1:100, 1:10]), 2), "\n")
#> Data layer range: 0 2.77
crc
#> An object of class Seurat 
#> 25344 features across 935 samples within 1 assay 
#> Active assay: RNA (25344 features, 3466 variable features)
#>  2 layers present: counts, data
#>  1 dimensional reduction calculated: umap
```

### Visualize Converted Data

The UMAP coordinates and normalized expression from scanpy are
preserved:

``` r

DimPlot(crc, reduction = "umap", group.by = "tissue", pt.size = 0.5)
```

![](convert-anndata_files/figure-html/plot_crc-1.png)

``` r

FeaturePlot(crc, features = "EPCAM", pt.size = 0.5)
```

![](convert-anndata_files/figure-html/plot_crc_epcam-1.png)

The conversion preserves expression patterns - EPCAM shows consistent
distribution in both tools.

## Visium Spatial Data Conversion

For spatial transcriptomics data, we use the stxBrain dataset from
SeuratData (Visium v2 format):

``` r

library(SeuratData)
if (!"stxBrain" %in% rownames(InstalledData())) {
  InstallData("stxBrain")
}

brain <- UpdateSeuratObject(LoadData("stxBrain", type = "anterior1"))
brain <- NormalizeData(brain)
cat("Layers:", paste(Layers(brain), collapse = ", "), "\n")
#> Layers: counts, data

SpatialFeaturePlot(brain, features = "Hpca")
```

![](convert-anndata_files/figure-html/get_stxBrain-1.png)

Convert to h5ad using the direct pipeline:

``` r

writeH5AD(brain, "stxBrain.h5ad", overwrite = TRUE)
```

View in Python with Squidpy:

``` python
import squidpy as sq
import scanpy as sc

adata_spatial = sc.read_h5ad("stxBrain.h5ad")
print(adata_spatial)
#> AnnData object with n_obs × n_vars = 2696 × 31053
#>     obs: 'orig.ident', 'nCount_Spatial', 'nFeature_Spatial', 'slice', 'region'
#>     uns: 'spatial'
#>     obsm: 'spatial'

sq.pl.spatial_scatter(adata_spatial, color="Hpca", library_id="anterior1",
                      img_res_key="lowres", size=1, alpha=0.5, use_raw=False)
```

![](convert-anndata_files/figure-html/load_stxBrain-1.png)

## Multi-assay Conversion (CITE-seq)

For multi-modal data like CITE-seq, each assay must be converted
separately since h5ad format only supports a single matrix per file.

**Conversion behavior:**

- **No assay specified**: Only the default assay is converted
- **Single assay specified**: That specific assay is converted
- **Multiple assays**: Must call the function multiple times, once per
  assay

This example uses the `cbmc` dataset from SeuratData:

``` r

library(SeuratData)
if (!"cbmc" %in% rownames(InstalledData())) {
  InstallData("cbmc")
}
data("cbmc", package = "cbmc.SeuratData")
cbmc <- UpdateSeuratObject(cbmc)

cat("Assays:", paste(Assays(cbmc), collapse = ", "), "\n")
#> Assays: RNA, ADT

# Process RNA for visualization
DefaultAssay(cbmc) <- "RNA"
cbmc <- NormalizeData(cbmc, verbose = FALSE)
cbmc <- FindVariableFeatures(cbmc, verbose = FALSE)
cbmc <- ScaleData(cbmc, verbose = FALSE)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30, verbose = FALSE)

# Normalize ADT with CLR
cbmc <- NormalizeData(cbmc, assay = "ADT", normalization.method = "CLR", verbose = FALSE)
```

Visualize RNA and ADT markers before conversion:

``` r

library(patchwork)
p1 <- FeaturePlot(cbmc, features = "CD3D", pt.size = 0.3) + ggtitle("CD3D (RNA)")
DefaultAssay(cbmc) <- "ADT"
p2 <- FeaturePlot(cbmc, features = "adt_CD3", pt.size = 0.3) + ggtitle("CD3 protein (ADT)")
p1 + p2
```

![](convert-anndata_files/figure-html/cite_seq_featureplot-1.png)

Convert RNA assay to h5ad:

``` r

DefaultAssay(cbmc) <- "RNA"
writeH5Seurat(cbmc, "cbmc_rna.h5seurat", overwrite = TRUE)
scConvert("cbmc_rna.h5seurat", dest = "h5ad", overwrite = TRUE)
```

Verify and visualize the converted file in scanpy:

``` python
import scanpy as sc
adata_rna = sc.read_h5ad("cbmc_rna.h5ad")
print(adata_rna)
#> AnnData object with n_obs × n_vars = 8617 × 20501
#>     obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'nCount_ADT', 'nFeature_ADT', 'rna_annotations', 'protein_annotations'
#>     var: 'vst.mean', 'vst.variance', 'vst.variance.expected', 'vst.variance.standardized', 'vst.variable'
#>     uns: 'n_variable_features', 'seurat'
#>     obsm: 'X_pca', 'X_umap'
#>     varm: 'PCs'

sc.pl.umap(adata_rna, color='CD3D', use_raw=False, title='CD3D (RNA, from h5ad)')
```

![](convert-anndata_files/figure-html/verify_cbmc_h5ad-1.png)

The h5ad file preserves the full structure: expression matrix (X), cell
metadata (obs), gene info (var), and UMAP coordinates.

Convert the ADT assay separately:

``` r

# Convert ADT assay
DefaultAssay(cbmc) <- "ADT"
writeH5Seurat(cbmc, "cbmc_adt.h5seurat", overwrite = TRUE)
scConvert("cbmc_adt.h5seurat", dest = "h5ad", overwrite = TRUE)
```

``` python
import scanpy as sc
adata_adt = sc.read_h5ad("cbmc_adt.h5ad")
print(adata_adt)
#> AnnData object with n_obs × n_vars = 8617 × 10
#>     obs: 'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'nCount_ADT', 'nFeature_ADT', 'rna_annotations', 'protein_annotations'
#>     var: 'features'
#>     uns: 'seurat'
#>     obsm: 'X_umap'
print("ADT features:", list(adata_adt.var_names)[:10])
#> ADT features: ['CD3', 'CD4', 'CD8', 'CD45RA', 'CD56', 'CD16', 'CD11c', 'CD14', 'CD19', 'CD34']
```

Each assay produces a separate h5ad file with the same cell metadata but
different features.

> **Note on multi-modal formats**: For true multi-modal
> interoperability, consider the [MuData/h5mu
> format](https://muon.scverse.org/) from the scverse ecosystem.
> scConvert provides native
> [`writeH5MU()`](https://mianaz.github.io/scConvert/reference/writeH5MU.md)
> and
> [`readH5MU()`](https://mianaz.github.io/scConvert/reference/readH5MU.md)
> functions for reading and writing h5mu files with no external
> dependencies. See
> [`vignette("multimodal-h5mu")`](https://mianaz.github.io/scConvert/articles/multimodal-h5mu.md)
> for details.

## Spatial h5ad to Seurat

Converting native spatial h5ad files to Seurat is fully supported. We
use a Visium colon sample from
[CellxGene](https://cellxgene.cziscience.com) that was processed with
scanpy/squidpy standard workflows:

``` r

cache_dir <- tools::R_user_dir("scConvert", which = "cache")
cache_path <- file.path(cache_dir, "visium_colon.h5ad")

if (!file.exists(cache_path)) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  message("Downloading Visium colon dataset (~1.7GB)...")
  download.file(
    "https://datasets.cellxgene.cziscience.com/ab9f4860-c0e3-444b-a982-38c13f0be6f5.h5ad",
    cache_path, mode = "wb"
  )
}
file.copy(cache_path, "visium_colon.h5ad", overwrite = TRUE)
#> [1] TRUE
```

View the native spatial h5ad in Python:

``` python
import scanpy as sc
adata_spatial = sc.read_h5ad("visium_colon.h5ad")
print(adata_spatial)
#> AnnData object with n_obs × n_vars = 4992 × 32397
#>     obs: 'in_tissue', 'array_row', 'array_col', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'sangerID', 'region', 'donor_type', 'age', 'facility', 'flushed', 'annotation_final', 'Adip1', 'Adip2', 'Adip3', 'B', 'B_plasma', 'CD14+Mo', 'CD16+Mo', 'CD4+T_act', 'CD4+T_naive', 'CD8+T_cytox', 'CD8+T_em', 'CD8+T_te', 'CD8+T_trans', 'DC', 'EC10_CMC-like', 'EC1_cap', 'EC2_cap', 'EC3_cap', 'EC4_immune', 'EC5_art', 'EC6_ven', 'EC7_endocardial', 'EC8_ln', 'FB1', 'FB2', 'FB3', 'FB4_activated', 'FB5', 'FB6', 'ILC', 'LYVE1+IGF1+MP', 'LYVE1+MP_cycling', 'LYVE1+TIMD4+MP', 'MAIT-like', 'Mast', 'Meso', 'MoMP', 'NC1_glial', 'NC2_glial_NGF+', 'NK_CD16hi', 'NK_CD56hi', 'Neut', 'PC1_vent', 'PC2_atria', 'PC3_str', 'SAN_P_cell', 'SMC1_basic', 'SMC2_art', 'T/NK_cycling', 'aCM1', 'aCM2', 'aCM3', 'aCM4', 'AVN_bundle_cell', 'PC4_CMC-like', 'vCM1', 'vCM2', 'vCM3_stressed', 'vCM4', 'vCM5', 'AVN_P_cell', 'CD4+T_Tfh', 'CD4+T_Th1', 'CD4+T_Th2', 'CD4+T_reg', 'NC5_glial', 'aCM5', 'Adip4', 'NC3_glial', 'NC6_schwann', 'EC9_FB-like', 'gdT', 'Adip1_abundance', 'Adip2_abundance', 'Adip3_abundance', 'B_abundance', 'B_plasma_abundance', 'CD14+Mo_abundance', 'CD16+Mo_abundance', 'CD4+T_act_abundance', 'CD4+T_naive_abundance', 'CD8+T_cytox_abundance', 'CD8+T_em_abundance', 'CD8+T_te_abundance', 'CD8+T_trans_abundance', 'DC_abundance', 'EC10_CMC-like_abundance', 'EC1_cap_abundance', 'EC2_cap_abundance', 'EC3_cap_abundance', 'EC4_immune_abundance', 'EC5_art_abundance', 'EC6_ven_abundance', 'EC7_endocardial_abundance', 'EC8_ln_abundance', 'FB1_abundance', 'FB2_abundance', 'FB3_abundance', 'FB4_activated_abundance', 'FB5_abundance', 'FB6_abundance', 'ILC_abundance', 'LYVE1+IGF1+MP_abundance', 'LYVE1+MP_cycling_abundance', 'LYVE1+TIMD4+MP_abundance', 'MAIT-like_abundance', 'Mast_abundance', 'Meso_abundance', 'MoMP_abundance', 'NC1_glial_abundance', 'NC2_glial_NGF+_abundance', 'NK_CD16hi_abundance', 'NK_CD56hi_abundance', 'Neut_abundance', 'PC1_vent_abundance', 'PC2_atria_abundance', 'PC3_str_abundance', 'SAN_P_cell_abundance', 'SMC1_basic_abundance', 'SMC2_art_abundance', 'T/NK_cycling_abundance', 'aCM1_abundance', 'aCM2_abundance', 'aCM3_abundance', 'aCM4_abundance', 'AVN_bundle_cell_abundance', 'PC4_CMC-like_abundance', 'vCM1_abundance', 'vCM2_abundance', 'vCM3_stressed_abundance', 'vCM4_abundance', 'vCM5_abundance', 'AVN_P_cell_abundance', 'CD4+T_Tfh_abundance', 'CD4+T_Th1_abundance', 'CD4+T_Th2_abundance', 'CD4+T_reg_abundance', 'NC5_glial_abundance', 'aCM5_abundance', 'Adip4_abundance', 'NC3_glial_abundance', 'NC6_schwann_abundance', 'EC9_FB-like_abundance', 'gdT_abundance', 'assay_ontology_term_id', 'cell_type_ontology_term_id', 'donor_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'is_primary_data', 'self_reported_ethnicity_ontology_term_id', 'sex_ontology_term_id', 'suspension_type', 'tissue_ontology_term_id', 'tissue_type', 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
#>     var: 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'feature_length', 'feature_type'
#>     uns: 'cell_type_ontology_term_id_colors', 'citation', 'fullres_xml_metadata', 'organism', 'organism_ontology_term_id', 'schema_reference', 'schema_version', 'spatial', 'spatial_metadata', 'title'
#>     obsm: 'X_means_cell_abundance_w_sf', 'X_prop', 'X_q05_cell_abundance_w_sf', 'X_q95_cell_abundance_w_sf', 'X_stds_cell_abundance_w_sf', 'spatial'

lib_id = list(adata_spatial.uns['spatial'].keys())[0]
spatial_data = adata_spatial.uns['spatial'][lib_id]
print("\nSpatial library:", lib_id)
#> 
#> Spatial library: HCAHeartST13228106
print("Image keys:", list(spatial_data.get('images', {}).keys()))
#> Image keys: ['fullres', 'hires']
print("in_tissue column:", 'in_tissue' in adata_spatial.obs.columns)
#> in_tissue column: True
```

Normalize, set gene symbols, and save for conversion:

``` python
import squidpy as sq
import scanpy as sc
import pandas as pd

adata_spatial = sc.read_h5ad("visium_colon.h5ad")

# Filter to tissue spots only
adata_spatial = adata_spatial[adata_spatial.obs['in_tissue'] == 1].copy()

sc.pp.normalize_total(adata_spatial, target_sum=1e4)
sc.pp.log1p(adata_spatial)

adata_spatial.var_names = pd.Index(adata_spatial.var['feature_name'].astype(str).values)
adata_spatial.var_names_make_unique()

adata_spatial.write_h5ad("visium_normalized.h5ad")

lib_id = list(adata_spatial.uns['spatial'].keys())[0]
sq.pl.spatial_scatter(adata_spatial, color='ACTC1', library_id=lib_id,
                      size=1, alpha=0.5, use_raw=False, title='ACTC1')
```

![](convert-anndata_files/figure-html/squidpy_spatial_ACTC1-1.png)

Convert to Seurat:

``` r

scConvert("visium_normalized.h5ad", dest = "h5seurat", overwrite = TRUE)
visium <- readH5Seurat("visium_normalized.h5seurat")

# Verify layer mapping: X -> data (log-normalized)
cat("Layers:", paste(Layers(visium), collapse = ", "), "\n")
#> Layers: counts, data
visium
#> An object of class Seurat 
#> 32397 features across 3452 samples within 1 assay 
#> Active assay: RNA (32397 features, 0 variable features)
#>  2 layers present: counts, data
#>  5 dimensional reductions calculated: means_cell_abundance_w_sf, prop, q05_cell_abundance_w_sf, q95_cell_abundance_w_sf, stds_cell_abundance_w_sf
#>  1 spatial field of view present: HCAHeartST13228106
```

Verify spatial data was preserved:

``` r

SpatialFeaturePlot(visium, features = "ACTC1", pt.size.factor = 3, alpha = 1)
```

![](convert-anndata_files/figure-html/plot_spatial_native-1.png)

> **Note**: Spatial images and coordinates are preserved during
> conversion. Some scanpy-specific structures (like neighbor graphs in
> `obsp`) may need to be recomputed in Seurat using
> [`FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html).

## Data Mapping Reference

For complete data mapping tables covering all supported formats (h5ad,
h5Seurat, h5mu, Loom, Zarr), including layer mapping, metadata
conventions, indexing, and structure verification, see
[`vignette("data-mapping-reference")`](https://mianaz.github.io/scConvert/articles/data-mapping-reference.md).

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
#> other attached packages:
#>  [1] patchwork_1.3.2               ggplot2_4.0.2                
#>  [3] scConvert_0.1.0               Seurat_5.4.0                 
#>  [5] SeuratObject_5.3.0            sp_2.2-1                     
#>  [7] stxKidney.SeuratData_0.1.0    stxBrain.SeuratData_0.1.2    
#>  [9] ssHippo.SeuratData_3.1.4      pbmcref.SeuratData_1.0.0     
#> [11] pbmcMultiome.SeuratData_0.1.4 pbmc3k.SeuratData_3.1.4      
#> [13] panc8.SeuratData_3.0.2        cbmc.SeuratData_3.1.4        
#> [15] SeuratData_0.2.2.9002        
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.4        
#>   [4] spatstat.utils_3.2-2   farver_2.1.2           rmarkdown_2.30        
#>   [7] fs_1.6.7               ragg_1.5.0             vctrs_0.7.1           
#>  [10] ROCR_1.0-12            spatstat.explore_3.7-0 htmltools_0.5.9       
#>  [13] sass_0.4.10            sctransform_0.4.3      parallelly_1.46.1     
#>  [16] KernSmooth_2.23-26     bslib_0.10.0           htmlwidgets_1.6.4     
#>  [19] desc_1.4.3             ica_1.0-3              plyr_1.8.9            
#>  [22] plotly_4.12.0          zoo_1.8-15             cachem_1.1.0          
#>  [25] igraph_2.2.2           mime_0.13              lifecycle_1.0.5       
#>  [28] pkgconfig_2.0.3        Matrix_1.7-4           R6_2.6.1              
#>  [31] fastmap_1.2.0          MatrixGenerics_1.22.0  fitdistrplus_1.2-6    
#>  [34] future_1.69.0          shiny_1.13.0           digest_0.6.39         
#>  [37] S4Vectors_0.48.0       tensor_1.5.1           RSpectra_0.16-2       
#>  [40] irlba_2.3.7            GenomicRanges_1.62.1   textshaping_1.0.4     
#>  [43] labeling_0.4.3         progressr_0.18.0       spatstat.sparse_3.1-0 
#>  [46] httr_1.4.8             polyclip_1.10-7        abind_1.4-8           
#>  [49] compiler_4.5.2         bit64_4.6.0-1          withr_3.0.2           
#>  [52] S7_0.2.1               fastDummies_1.7.5      MASS_7.3-65           
#>  [55] rappdirs_0.3.4         tools_4.5.2            lmtest_0.9-40         
#>  [58] otel_0.2.0             httpuv_1.6.16          future.apply_1.20.2   
#>  [61] goftest_1.2-3          glue_1.8.0             nlme_3.1-168          
#>  [64] promises_1.5.0         grid_4.5.2             Rtsne_0.17            
#>  [67] cluster_2.1.8.2        reshape2_1.4.5         generics_0.1.4        
#>  [70] hdf5r_1.3.12           gtable_0.3.6           spatstat.data_3.1-9   
#>  [73] tidyr_1.3.2            data.table_1.18.2.1    XVector_0.50.0        
#>  [76] BiocGenerics_0.56.0    BPCells_0.2.0          spatstat.geom_3.7-0   
#>  [79] RcppAnnoy_0.0.23       ggrepel_0.9.7          RANN_2.6.2            
#>  [82] pillar_1.11.1          stringr_1.6.0          spam_2.11-3           
#>  [85] RcppHNSW_0.6.0         later_1.4.8            splines_4.5.2         
#>  [88] dplyr_1.2.0            lattice_0.22-9         bit_4.6.0             
#>  [91] survival_3.8-6         deldir_2.0-4           tidyselect_1.2.1      
#>  [94] miniUI_0.1.2           pbapply_1.7-4          knitr_1.51            
#>  [97] gridExtra_2.3          Seqinfo_1.0.0          IRanges_2.44.0        
#> [100] scattermore_1.2        stats4_4.5.2           xfun_0.56             
#> [103] matrixStats_1.5.0      UCSC.utils_1.6.1       stringi_1.8.7         
#> [106] lazyeval_0.2.2         yaml_2.3.12            evaluate_1.0.5        
#> [109] codetools_0.2-20       tibble_3.3.1           cli_3.6.5             
#> [112] uwot_0.2.4             xtable_1.8-8           reticulate_1.45.0     
#> [115] systemfonts_1.3.1      jquerylib_0.1.4        GenomeInfoDb_1.46.2   
#> [118] dichromat_2.0-0.1      Rcpp_1.1.1             globals_0.19.1        
#> [121] spatstat.random_3.4-4  png_0.1-8              spatstat.univar_3.1-6 
#> [124] parallel_4.5.2         pkgdown_2.2.0          dotCall64_1.2         
#> [127] listenv_0.10.1         viridisLite_0.4.3      scales_1.4.0          
#> [130] ggridges_0.5.7         purrr_1.2.1            crayon_1.5.3          
#> [133] rlang_1.1.7            cowplot_1.2.0
```
