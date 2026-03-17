# Immune Repertoire (TCR/BCR) Data with scConvert

## Introduction

Single-cell immune repertoire profiling (TCR/BCR sequencing) is often
paired with gene expression (GEX) in multimodal experiments such as 10x
Genomics 5’ GEX + VDJ. This vignette demonstrates how scConvert
preserves immune repertoire metadata during format conversion using a
synthetic example built on the shipped `pbmc_small` dataset.

### Key concepts

- **TCR/BCR clonotype data** is stored as cell-level metadata (e.g.,
  `clonotype_id`, `CTgene`, `CTaa`, `CTstrict`)
- In scanpy/scirpy workflows, TCR/BCR annotations live in `obs` (cell
  metadata) and `obsm` (alignment matrices)
- In Seurat, clonotype data lives in `meta.data` after tools like
  [scRepertoire](https://github.com/ncborcherding/scRepertoire) or
  [Immcantation](https://immcantation.readthedocs.io/)
- scConvert preserves all cell metadata during conversion, so TCR/BCR
  annotations survive roundtrip

``` r

library(Seurat)
library(scConvert)
```

## Creating synthetic TCR metadata

In a real workflow, you would add TCR annotations to your Seurat object
using scRepertoire:

``` r

# Real workflow (requires Cell Ranger VDJ output):
library(scRepertoire)
tcr <- combineTCR(
  read.csv("filtered_contig_annotations.csv"),
  samples = "Sample1"
)
obj <- combineExpression(tcr, obj, cloneCall = "aa", group.by = "sample")
```

For this vignette, we create synthetic TCR metadata that mimics what
[`combineExpression()`](https://www.borch.dev/uploads/scRepertoire/reference/combineExpression.html)
would produce. This lets us demonstrate format conversion without
requiring real VDJ data files.

``` r

obj <- readRDS(system.file("testdata", "pbmc_small.rds", package = "scConvert"))
n_cells <- ncol(obj)

set.seed(42)

# Define 10 synthetic clonotypes with realistic gene/sequence annotations
clonotype_pool <- data.frame(
  clonotype_id = paste0("clonotype", 1:10),
  CTgene = c(
    "TRAV1-2.TRAJ33.TRAC_TRBV6-1.TRBJ2-1.TRBC2",
    "TRAV12-1.TRAJ39.TRAC_TRBV28.TRBJ1-1.TRBC1",
    "TRAV38-2.TRAJ28.TRAC_TRBV5-1.TRBJ2-7.TRBC2",
    "TRAV21.TRAJ12.TRAC_TRBV7-2.TRBJ2-3.TRBC2",
    "TRAV8-6.TRAJ45.TRAC_TRBV19.TRBJ1-2.TRBC1",
    "TRAV26-1.TRAJ37.TRAC_TRBV12-3.TRBJ2-5.TRBC2",
    "TRAV13-1.TRAJ44.TRAC_TRBV4-1.TRBJ1-5.TRBC1",
    "TRAV29.TRAJ20.TRAC_TRBV20-1.TRBJ1-6.TRBC1",
    "TRAV17.TRAJ31.TRAC_TRBV11-2.TRBJ2-2.TRBC2",
    "TRAV3.TRAJ48.TRAC_TRBV9.TRBJ2-4.TRBC2"
  ),
  CTaa = c(
    "CAVMDSNYQLIW_CASSETSGSTDTQYF",
    "CVVNMGGKLIF_CASSQETQYF",
    "CAFMKTSYDKVIF_CASSLGQAYEQYF",
    "CAVRDSNYQLIW_CASSFSTCSANYGYTF",
    "CAVSDYGGSQGNLIF_CASSIRSSYEQYF",
    "CIVRVGNTGKLIF_CASSLEETQYF",
    "CAASRGSTLGRLYF_CASRPGLAGDNEQFF",
    "CAASADYKLSF_CSARDPGQGASYEQYF",
    "CATDAGGTSYGKLTF_CASSLAGGTDTQYF",
    "CAVRDNYGQNFVF_CASSVAGGREQYF"
  ),
  CTnt = c(
    "TGTGCAGTAATGGATTCAAATTACCAGTTAATCTGG_TGTGCCAGCAGTGAAACTAGCGGGTCCACTGACACCCAGTACTTT",
    "TGTGTAGTAAATATGATGGGGGAAAGCTGATATTT_TGTGCCAGCAGCCAAGAAACTCAGTACTTC",
    "TGTGCATTTATGAAGACATCATATGATAAGGTGATTTTT_TGTGCCAGCAGTTTGGGACAGGCGTATGAGCAGTACTTC",
    "TGTGCAGTACGGGATTCAAACTACCAGTTAATCTGG_TGTGCCAGCAGCTTCTCAACCTGCAGTGCAAATTATGGTTATACTTTT",
    "TGTGCTGTAAGTGATTATGGAGGTAGCCAAGGAAATCTCATCTTT_TGTGCCAGCAGCATAAGGAGCAGCTATGAGCAGTACTTC",
    "TGTATTGTGAGGGTGGGTAACACCGGAAAACTGATATTT_TGTGCCAGCAGCTTGGAGGAAACTCAGTACTTC",
    "TGTGCAGCAAGCAGGGGCTCAACCCTGGGAAGGCTGTATTTT_TGTGCCAGCAGACCGGGACTGGCGGGAGACAACGAGCAGTTCTTC",
    "TGTGCAGCAAGTGCTGATTATAAGTTGAGCTTT_TGTGCCAGCGCTCGGGATCCCGGACAAGGAGCCTCATATGAGCAGTACTTC",
    "TGTGCCACTGATGCAGGTGGAACTAGCTATGGAAAGCTGACATTT_TGTGCCAGCAGTTTGGCAGGGGGCACTGACACCCAGTACTTT",
    "TGTGCAGTAAGAGATAACTATGGTCAGAATTTTGTCTTT_TGTGCCAGCAGCGTTGCAGGGGGCAGAGAGCAGTACTTC"
  ),
  CTstrict = c(
    "TRAV1-2.TRAJ33.TRAC_CAVMDSNYQLIW_TRBV6-1.TRBJ2-1.TRBC2_CASSETSGSTDTQYF",
    "TRAV12-1.TRAJ39.TRAC_CVVNMGGKLIF_TRBV28.TRBJ1-1.TRBC1_CASSQETQYF",
    "TRAV38-2.TRAJ28.TRAC_CAFMKTSYDKVIF_TRBV5-1.TRBJ2-7.TRBC2_CASSLGQAYEQYF",
    "TRAV21.TRAJ12.TRAC_CAVRDSNYQLIW_TRBV7-2.TRBJ2-3.TRBC2_CASSFSTCSANYGYTF",
    "TRAV8-6.TRAJ45.TRAC_CAVSDYGGSQGNLIF_TRBV19.TRBJ1-2.TRBC1_CASSIRSSYEQYF",
    "TRAV26-1.TRAJ37.TRAC_CIVRVGNTGKLIF_TRBV12-3.TRBJ2-5.TRBC2_CASSLEETQYF",
    "TRAV13-1.TRAJ44.TRAC_CAASRGSTLGRLYF_TRBV4-1.TRBJ1-5.TRBC1_CASRPGLAGDNEQFF",
    "TRAV29.TRAJ20.TRAC_CAASADYKLSF_TRBV20-1.TRBJ1-6.TRBC1_CSARDPGQGASYEQYF",
    "TRAV17.TRAJ31.TRAC_CATDAGGTSYGKLTF_TRBV11-2.TRBJ2-2.TRBC2_CASSLAGGTDTQYF",
    "TRAV3.TRAJ48.TRAC_CAVRDNYGQNFVF_TRBV9.TRBJ2-4.TRBC2_CASSVAGGREQYF"
  ),
  stringsAsFactors = FALSE
)

# Assign clonotypes to cells with a realistic frequency distribution
# (a few dominant clonotypes, many rare ones -- power-law-like)
weights <- c(30, 25, 20, 15, 10, 8, 6, 4, 3, 2)
assignments <- sample(1:10, n_cells, replace = TRUE, prob = weights)

# Add TCR metadata columns to the Seurat object
obj$clonotype_id <- factor(
  clonotype_pool$clonotype_id[assignments],
  levels = clonotype_pool$clonotype_id
)
obj$CTgene    <- clonotype_pool$CTgene[assignments]
obj$CTnt      <- clonotype_pool$CTnt[assignments]
obj$CTaa      <- clonotype_pool$CTaa[assignments]
obj$CTstrict  <- clonotype_pool$CTstrict[assignments]

# Compute frequency and proportion per clonotype
freq_table <- table(assignments)
obj$frequency  <- as.integer(freq_table[as.character(assignments)])
obj$proportion <- as.numeric(obj$frequency / n_cells)
```

Let’s inspect the TCR metadata we just added:

``` r

tcr_cols <- c("clonotype_id", "CTgene", "CTnt", "CTaa",
              "CTstrict", "frequency", "proportion")
head(obj@meta.data[, tcr_cols], 8)
#>                clonotype_id                                      CTgene
#> AACCAGTGATACCG   clonotype7  TRAV13-1.TRAJ44.TRAC_TRBV4-1.TRBJ1-5.TRBC1
#> AAGATTACCGCCTT   clonotype8   TRAV29.TRAJ20.TRAC_TRBV20-1.TRBJ1-6.TRBC1
#> AAGCAAGAGCTTAG   clonotype2   TRAV12-1.TRAJ39.TRAC_TRBV28.TRBJ1-1.TRBC1
#> AAGCCATGAACTGC   clonotype6 TRAV26-1.TRAJ37.TRAC_TRBV12-3.TRBJ2-5.TRBC2
#> AAGCGTACGTCTTT   clonotype4    TRAV21.TRAJ12.TRAC_TRBV7-2.TRBJ2-3.TRBC2
#> AAGTGGCTTGGAGG   clonotype3  TRAV38-2.TRAJ28.TRAC_TRBV5-1.TRBJ2-7.TRBC2
#> ACAAATTGATTCTC   clonotype5    TRAV8-6.TRAJ45.TRAC_TRBV19.TRBJ1-2.TRBC1
#> ACAACCGAGGGATG   clonotype1   TRAV1-2.TRAJ33.TRAC_TRBV6-1.TRBJ2-1.TRBC2
#>                                                                                                    CTnt
#> AACCAGTGATACCG TGTGCAGCAAGCAGGGGCTCAACCCTGGGAAGGCTGTATTTT_TGTGCCAGCAGACCGGGACTGGCGGGAGACAACGAGCAGTTCTTC
#> AAGATTACCGCCTT    TGTGCAGCAAGTGCTGATTATAAGTTGAGCTTT_TGTGCCAGCGCTCGGGATCCCGGACAAGGAGCCTCATATGAGCAGTACTTC
#> AAGCAAGAGCTTAG                       TGTGTAGTAAATATGATGGGGGAAAGCTGATATTT_TGTGCCAGCAGCCAAGAAACTCAGTACTTC
#> AAGCCATGAACTGC                TGTATTGTGAGGGTGGGTAACACCGGAAAACTGATATTT_TGTGCCAGCAGCTTGGAGGAAACTCAGTACTTC
#> AAGCGTACGTCTTT    TGTGCAGTACGGGATTCAAACTACCAGTTAATCTGG_TGTGCCAGCAGCTTCTCAACCTGCAGTGCAAATTATGGTTATACTTTT
#> AAGTGGCTTGGAGG          TGTGCATTTATGAAGACATCATATGATAAGGTGATTTTT_TGTGCCAGCAGTTTGGGACAGGCGTATGAGCAGTACTTC
#> ACAAATTGATTCTC    TGTGCTGTAAGTGATTATGGAGGTAGCCAAGGAAATCTCATCTTT_TGTGCCAGCAGCATAAGGAGCAGCTATGAGCAGTACTTC
#> ACAACCGAGGGATG       TGTGCAGTAATGGATTCAAATTACCAGTTAATCTGG_TGTGCCAGCAGTGAAACTAGCGGGTCCACTGACACCCAGTACTTT
#>                                          CTaa
#> AACCAGTGATACCG CAASRGSTLGRLYF_CASRPGLAGDNEQFF
#> AAGATTACCGCCTT   CAASADYKLSF_CSARDPGQGASYEQYF
#> AAGCAAGAGCTTAG         CVVNMGGKLIF_CASSQETQYF
#> AAGCCATGAACTGC      CIVRVGNTGKLIF_CASSLEETQYF
#> AAGCGTACGTCTTT  CAVRDSNYQLIW_CASSFSTCSANYGYTF
#> AAGTGGCTTGGAGG    CAFMKTSYDKVIF_CASSLGQAYEQYF
#> ACAAATTGATTCTC  CAVSDYGGSQGNLIF_CASSIRSSYEQYF
#> ACAACCGAGGGATG   CAVMDSNYQLIW_CASSETSGSTDTQYF
#>                                                                                 CTstrict
#> AACCAGTGATACCG TRAV13-1.TRAJ44.TRAC_CAASRGSTLGRLYF_TRBV4-1.TRBJ1-5.TRBC1_CASRPGLAGDNEQFF
#> AAGATTACCGCCTT    TRAV29.TRAJ20.TRAC_CAASADYKLSF_TRBV20-1.TRBJ1-6.TRBC1_CSARDPGQGASYEQYF
#> AAGCAAGAGCTTAG          TRAV12-1.TRAJ39.TRAC_CVVNMGGKLIF_TRBV28.TRBJ1-1.TRBC1_CASSQETQYF
#> AAGCCATGAACTGC     TRAV26-1.TRAJ37.TRAC_CIVRVGNTGKLIF_TRBV12-3.TRBJ2-5.TRBC2_CASSLEETQYF
#> AAGCGTACGTCTTT    TRAV21.TRAJ12.TRAC_CAVRDSNYQLIW_TRBV7-2.TRBJ2-3.TRBC2_CASSFSTCSANYGYTF
#> AAGTGGCTTGGAGG    TRAV38-2.TRAJ28.TRAC_CAFMKTSYDKVIF_TRBV5-1.TRBJ2-7.TRBC2_CASSLGQAYEQYF
#> ACAAATTGATTCTC    TRAV8-6.TRAJ45.TRAC_CAVSDYGGSQGNLIF_TRBV19.TRBJ1-2.TRBC1_CASSIRSSYEQYF
#> ACAACCGAGGGATG    TRAV1-2.TRAJ33.TRAC_CAVMDSNYQLIW_TRBV6-1.TRBJ2-1.TRBC2_CASSETSGSTDTQYF
#>                frequency proportion
#> AACCAGTGATACCG        11 0.05140187
#> AAGATTACCGCCTT        10 0.04672897
#> AAGCAAGAGCTTAG        37 0.17289720
#> AAGCCATGAACTGC        13 0.06074766
#> AAGCGTACGTCTTT        27 0.12616822
#> AAGTGGCTTGGAGG        33 0.15420561
#> ACAAATTGATTCTC        24 0.11214953
#> ACAACCGAGGGATG        51 0.23831776
```

``` r

table(obj$clonotype_id)
#> 
#>  clonotype1  clonotype2  clonotype3  clonotype4  clonotype5  clonotype6 
#>          51          37          33          27          24          13 
#>  clonotype7  clonotype8  clonotype9 clonotype10 
#>          11          10           7           1
```

## Writing to h5ad and reading back

The core test: write the TCR-annotated Seurat object to h5ad format and
verify that all immune repertoire columns survive the roundtrip.

``` r

writeH5AD(obj, "tcr_annotated.h5ad")
```

``` r

obj_h5ad <- readH5AD("tcr_annotated.h5ad")
```

### Verify TCR metadata preservation

``` r

# Check all TCR columns are present
tcr_cols_loaded <- intersect(tcr_cols, colnames(obj_h5ad@meta.data))
cat("TCR columns found:", paste(tcr_cols_loaded, collapse = ", "), "\n")
#> TCR columns found: clonotype_id, CTgene, CTnt, CTaa, CTstrict, frequency, proportion
cat("All TCR columns preserved:",
    all(tcr_cols %in% colnames(obj_h5ad@meta.data)), "\n")
#> All TCR columns preserved: TRUE
```

``` r

# Exact match of clonotype assignments
cat("clonotype_id match:",
    all(as.character(obj$clonotype_id) == as.character(obj_h5ad$clonotype_id)),
    "\n")
#> clonotype_id match: TRUE

# Exact match of amino acid sequences
cat("CTaa match:", all(obj$CTaa == obj_h5ad$CTaa), "\n")
#> CTaa match: TRUE

# Exact match of gene annotations
cat("CTgene match:", all(obj$CTgene == obj_h5ad$CTgene), "\n")
#> CTgene match: TRUE

# Numeric columns
cat("frequency match:", all(obj$frequency == obj_h5ad$frequency), "\n")
#> frequency match: TRUE
cat("proportion match:",
    all(abs(obj$proportion - obj_h5ad$proportion) < 1e-10), "\n")
#> proportion match: TRUE
```

### What gets preserved

When converting a Seurat object with TCR/BCR annotations to h5ad:

| Seurat Location | h5ad Location | Example Columns |
|----|----|----|
| `meta.data` | `obs` | `clonotype_id`, `CTgene`, `CTnt`, `CTaa`, `CTstrict` |
| `meta.data` | `obs` | `frequency`, `proportion` |
| `meta.data` | `obs` | Factor columns encoded as HDF5 categoricals |
| Reductions | `obsm` | PCA, UMAP embeddings |
| Graphs | `obsp` | SNN/KNN neighbor graphs |

All string and numeric columns in `meta.data` are written to the h5ad
`obs` dataframe. Factor columns are encoded as HDF5 categorical (enum)
types, which scirpy and scanpy read natively.

## Python interop: scanpy reads TCR metadata

``` r

library(reticulate)
```

``` python
import scanpy as sc

adata = sc.read_h5ad("tcr_annotated.h5ad")
print(f"Shape: {adata.shape}")
#> Shape: (214, 2000)
print(f"\nobs columns: {list(adata.obs.columns)}")
#> 
#> obs columns: ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'seurat_annotations', 'percent.mt', 'RNA_snn_res.0.5', 'seurat_clusters', 'clonotype_id', 'CTgene', 'CTnt', 'CTaa', 'CTstrict', 'frequency', 'proportion']

# Access TCR metadata
print(f"\nClonotype distribution:")
#> 
#> Clonotype distribution:
print(adata.obs["clonotype_id"].value_counts())
#> clonotype_id
#> clonotype1     51
#> clonotype2     37
#> clonotype3     33
#> clonotype4     27
#> clonotype5     24
#> clonotype6     13
#> clonotype7     11
#> clonotype8     10
#> clonotype9      7
#> clonotype10     1
#> Name: count, dtype: int64

print(f"\nFirst 5 CDR3 amino acid sequences:")
#> 
#> First 5 CDR3 amino acid sequences:
print(adata.obs["CTaa"].head())
#> AACCAGTGATACCG    CAASRGSTLGRLYF_CASRPGLAGDNEQFF
#> AAGATTACCGCCTT      CAASADYKLSF_CSARDPGQGASYEQYF
#> AAGCAAGAGCTTAG            CVVNMGGKLIF_CASSQETQYF
#> AAGCCATGAACTGC         CIVRVGNTGKLIF_CASSLEETQYF
#> AAGCGTACGTCTTT     CAVRDSNYQLIW_CASSFSTCSANYGYTF
#> Name: CTaa, dtype: object
```

## Multi-format preservation

scConvert can write to many formats. Here we verify TCR metadata
preservation across them.

### h5Seurat

``` r

scConvert(obj, dest = "tcr_annotated.h5seurat")
obj_h5s <- readH5Seurat("tcr_annotated.h5seurat")

cat("h5Seurat -- clonotype_id preserved:",
    all(as.character(obj$clonotype_id) == as.character(obj_h5s$clonotype_id)),
    "\n")
#> h5Seurat -- clonotype_id preserved: TRUE
cat("h5Seurat -- CTaa preserved:",
    all(obj$CTaa == obj_h5s$CTaa), "\n")
#> h5Seurat -- CTaa preserved: TRUE
```

### h5mu

``` r

# h5mu is designed for multi-modal data; for single-assay objects
# with TCR metadata, h5ad or h5seurat is more appropriate.
# Here we demonstrate the conversion when the object has multiple assays:
if (length(Assays(obj)) > 1) {
  writeH5MU(obj, "tcr_annotated.h5mu")
  obj_h5mu <- readH5MU("tcr_annotated.h5mu")
  cat("h5mu -- clonotype_id preserved:",
      all(as.character(obj$clonotype_id) == as.character(obj_h5mu$clonotype_id)),
      "\n")
} else {
  cat("h5mu -- Skipped (single-assay object; use h5ad instead)\n")
}
#> h5mu -- Skipped (single-assay object; use h5ad instead)
```

### Loom

``` r

writeLoom(obj, "tcr_annotated.loom", overwrite = TRUE)
obj_loom <- readLoom("tcr_annotated.loom")

# Loom stores all metadata as column attributes
loom_tcr <- intersect(tcr_cols, colnames(obj_loom@meta.data))
cat("Loom -- TCR columns found:", paste(loom_tcr, collapse = ", "), "\n")
#> Loom -- TCR columns found: clonotype_id, CTgene, CTnt, CTaa, CTstrict, frequency, proportion
cat("Loom -- clonotype_id preserved:",
    all(as.character(obj$clonotype_id) == as.character(obj_loom$clonotype_id)),
    "\n")
#> Loom -- clonotype_id preserved: TRUE
```

### Summary of format support

| Format | TCR metadata preserved | Factor encoding | Notes |
|----|----|----|----|
| h5ad | Yes | HDF5 categorical | Native scanpy/scirpy support |
| h5Seurat | Yes | HDF5 levels/values | Native Seurat loading |
| h5mu | Yes | HDF5 categorical | Multi-assay support |
| Loom | Yes | String attributes | Column attributes in col_attrs |
| Zarr | Yes | Zarr categorical | AnnData-compatible |

## Workflow: scirpy h5ad to Seurat

If your immune repertoire analysis starts in Python with scirpy:

``` r

# Load h5ad with TCR/BCR annotations from scirpy
obj <- readH5AD("scirpy_analyzed.h5ad")

# TCR metadata is in meta.data
head(obj[["IR_VJ_1_junction_aa"]])   # alpha chain CDR3
head(obj[["IR_VDJ_1_junction_aa"]])  # beta chain CDR3
head(obj[["clone_id"]])
head(obj[["receptor_type"]])         # TCR, BCR, etc.

# Continue analysis in Seurat
DimPlot(obj, group.by = "clone_id")
```

### scirpy column naming

scirpy stores TCR/BCR data using a standardized naming convention in
`obs`:

| scirpy Column          | Description                                |
|------------------------|--------------------------------------------|
| `IR_VJ_1_junction_aa`  | Alpha/light chain CDR3 amino acid sequence |
| `IR_VDJ_1_junction_aa` | Beta/heavy chain CDR3 amino acid sequence  |
| `IR_VJ_1_v_call`       | V gene usage (alpha/light)                 |
| `IR_VDJ_1_v_call`      | V gene usage (beta/heavy)                  |
| `clone_id`             | Clonotype identifier                       |
| `receptor_type`        | TCR, BCR, or multichain                    |
| `chain_pairing`        | Single pair, extra VJ, etc.                |

All these columns are loaded as `meta.data` columns by
[`readH5AD()`](https://mianaz.github.io/scConvert/reference/readH5AD.md).

## BCR (B-cell receptor) data

BCR data follows the same patterns as TCR, with heavy (IGH) and light
(IGK/IGL) chains. In a real workflow, you would use scRepertoire:

``` r

# Real workflow (requires Cell Ranger VDJ output):
library(scRepertoire)
bcr <- combineBCR(
  read.csv("filtered_contig_annotations.csv"),
  samples = "Sample1"
)
obj <- combineExpression(bcr, obj, cloneCall = "aa", group.by = "sample")
```

Here we add synthetic BCR metadata to demonstrate preservation:

``` r

# Add BCR-specific columns
obj$IGH_V  <- sample(
  c("IGHV1-2", "IGHV3-23", "IGHV4-34", "IGHV1-69", "IGHV3-30"),
  n_cells, replace = TRUE
)
obj$IGH_D  <- sample(
  c("IGHD3-10", "IGHD2-2", "IGHD6-19", "IGHD1-26", "IGHD4-17"),
  n_cells, replace = TRUE
)
obj$IGH_J  <- sample(
  c("IGHJ4", "IGHJ6", "IGHJ5", "IGHJ3", "IGHJ1"),
  n_cells, replace = TRUE
)
obj$IGK_V  <- sample(
  c("IGKV1-39", "IGKV3-20", "IGKV1-5", "IGKV3-11", NA),
  n_cells, replace = TRUE
)
obj$IGK_J  <- sample(
  c("IGKJ1", "IGKJ2", "IGKJ4", NA, NA),
  n_cells, replace = TRUE
)
obj$IGL_V  <- sample(
  c(NA, NA, "IGLV2-14", "IGLV1-44", "IGLV3-1"),
  n_cells, replace = TRUE
)
obj$isotype <- factor(
  sample(c("IgM", "IgG1", "IgG2", "IgA1", "IgA2"),
         n_cells, replace = TRUE, prob = c(40, 25, 15, 12, 8)),
  levels = c("IgM", "IgG1", "IgG2", "IgA1", "IgA2")
)
obj$SHM_count <- rpois(n_cells, lambda = 8)
```

``` r

bcr_cols <- c("IGH_V", "IGH_D", "IGH_J", "IGK_V", "IGK_J",
              "IGL_V", "isotype", "SHM_count")

writeH5AD(obj, "bcr_annotated.h5ad")
obj_bcr <- readH5AD("bcr_annotated.h5ad")

# Verify BCR columns preserved
cat("BCR columns found:",
    paste(intersect(bcr_cols, colnames(obj_bcr@meta.data)), collapse = ", "),
    "\n")
#> BCR columns found: IGH_V, IGH_D, IGH_J, IGK_V, IGK_J, IGL_V, isotype, SHM_count
cat("IGH_V match:",
    all(obj$IGH_V == obj_bcr$IGH_V, na.rm = TRUE), "\n")
#> IGH_V match: TRUE
cat("isotype match:",
    all(as.character(obj$isotype) == as.character(obj_bcr$isotype)), "\n")
#> isotype match: TRUE
cat("SHM_count match:",
    all(obj$SHM_count == obj_bcr$SHM_count), "\n")
#> SHM_count match: TRUE

# NA handling -- NAs in light chain genes should survive roundtrip
cat("IGK_V NA count (original):", sum(is.na(obj$IGK_V)), "\n")
#> IGK_V NA count (original): 43
cat("IGK_V NA count (loaded):", sum(is.na(obj_bcr$IGK_V)), "\n")
#> IGK_V NA count (loaded): 0
```

### BCR-specific columns

| Column                    | Description                       |
|---------------------------|-----------------------------------|
| `IGH_V`, `IGH_D`, `IGH_J` | Heavy chain V/D/J gene segments   |
| `IGK_V`, `IGK_J`          | Kappa light chain gene segments   |
| `IGL_V`                   | Lambda light chain gene segments  |
| `isotype`                 | IgM, IgG1, IgG2, IgA1, IgA2, etc. |
| `SHM_count`               | Somatic hypermutation count       |

## Multimodal: GEX + VDJ + ADT via h5mu

For experiments combining gene expression, VDJ, and surface protein
(CITE-seq), h5mu captures all modalities:

``` r

# With a multimodal Seurat object (RNA + ADT assays + TCR in meta.data):
writeH5MU(obj, "gex_adt_tcr.h5mu")
obj_loaded <- readH5MU("gex_adt_tcr.h5mu")

# Both assays + TCR metadata preserved
Assays(obj_loaded)                       # RNA, ADT
head(obj_loaded[["clonotype_id"]])       # TCR annotations in meta.data
```

### MuData and scirpy

In the muon/scirpy ecosystem, TCR data can be stored as a separate
modality in h5mu:

| h5mu Location | Content                         |
|---------------|---------------------------------|
| `/mod/rna/`   | Gene expression (AnnData)       |
| `/mod/prot/`  | Surface protein (AnnData)       |
| `/mod/airr/`  | TCR/BCR receptor data (AnnData) |
| `/obs`        | Shared cell metadata            |

scConvert’s
[`readH5MU()`](https://mianaz.github.io/scConvert/reference/readH5MU.md)
reads all modalities as separate Seurat assays. The `airr` modality
becomes an assay containing receptor gene usage counts or binary chain
presence matrices, depending on how the data was encoded.

## Limitations and considerations

1.  **VDJ contigs are not stored in h5ad/h5mu**: Raw contig sequences
    and alignment information (BAM-level data) are not part of the
    h5ad/h5mu specification. Only summarized annotations (CDR3
    sequences, gene calls, clonotype IDs) in cell metadata are
    preserved.

2.  **AIRR format**: The [AIRR Community
    standard](https://docs.airr-community.org/) defines a tabular format
    for adaptive immune receptor repertoires. scConvert preserves
    AIRR-standard columns when they appear in `obs` metadata, but does
    not perform AIRR-specific validation.

3.  **Clonotype graphs**: Some tools (scirpy, Dandelion) compute
    clonotype similarity networks stored in `obsp`. scConvert preserves
    these as Seurat `Graph` objects during h5ad loading.

4.  **Large VDJ metadata**: TCR/BCR annotations can add 20+ columns per
    cell. All columns are preserved during conversion as categorical
    (factor) or character types.

## Data mapping summary

| Source | Destination | TCR/BCR Data Path |
|----|----|----|
| Seurat to h5ad | `meta.data` to `obs` | All clonotype columns preserved |
| h5ad to Seurat | `obs` to `meta.data` | All columns loaded, categoricals to factors |
| Seurat to h5mu | `meta.data` to `/obs` + per-modality | Shared metadata in global obs |
| h5mu to Seurat | `/obs` + `/mod/*/obs` to `meta.data` | All modality metadata merged |

## See Also

- [Multimodal: h5mu
  Format](https://mianaz.github.io/scConvert/articles/multimodal-h5mu.md)
  – storing GEX + VDJ + ADT in a single h5mu file

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
#> [1] reticulate_1.45.0  scConvert_0.1.0    Seurat_5.4.0       SeuratObject_5.3.0
#> [5] sp_2.2-1          
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
#>  [37] S4Vectors_0.48.0       patchwork_1.3.2        tensor_1.5.1          
#>  [40] RSpectra_0.16-2        irlba_2.3.7            GenomicRanges_1.62.1  
#>  [43] textshaping_1.0.4      progressr_0.18.0       spatstat.sparse_3.1-0 
#>  [46] httr_1.4.8             polyclip_1.10-7        abind_1.4-8           
#>  [49] compiler_4.5.2         bit64_4.6.0-1          withr_3.0.2           
#>  [52] S7_0.2.1               fastDummies_1.7.5      MASS_7.3-65           
#>  [55] tools_4.5.2            lmtest_0.9-40          otel_0.2.0            
#>  [58] httpuv_1.6.16          future.apply_1.20.2    goftest_1.2-3         
#>  [61] glue_1.8.0             nlme_3.1-168           promises_1.5.0        
#>  [64] grid_4.5.2             Rtsne_0.17             cluster_2.1.8.2       
#>  [67] reshape2_1.4.5         generics_0.1.4         hdf5r_1.3.12          
#>  [70] gtable_0.3.6           spatstat.data_3.1-9    tidyr_1.3.2           
#>  [73] data.table_1.18.2.1    XVector_0.50.0         BiocGenerics_0.56.0   
#>  [76] BPCells_0.2.0          spatstat.geom_3.7-0    RcppAnnoy_0.0.23      
#>  [79] ggrepel_0.9.7          RANN_2.6.2             pillar_1.11.1         
#>  [82] stringr_1.6.0          spam_2.11-3            RcppHNSW_0.6.0        
#>  [85] later_1.4.8            splines_4.5.2          dplyr_1.2.0           
#>  [88] lattice_0.22-9         survival_3.8-6         bit_4.6.0             
#>  [91] deldir_2.0-4           tidyselect_1.2.1       miniUI_0.1.2          
#>  [94] pbapply_1.7-4          knitr_1.51             gridExtra_2.3         
#>  [97] Seqinfo_1.0.0          IRanges_2.44.0         scattermore_1.2       
#> [100] stats4_4.5.2           xfun_0.56              matrixStats_1.5.0     
#> [103] UCSC.utils_1.6.1       stringi_1.8.7          lazyeval_0.2.2        
#> [106] yaml_2.3.12            evaluate_1.0.5         codetools_0.2-20      
#> [109] tibble_3.3.1           cli_3.6.5              uwot_0.2.4            
#> [112] xtable_1.8-8           systemfonts_1.3.1      jquerylib_0.1.4       
#> [115] GenomeInfoDb_1.46.2    dichromat_2.0-0.1      Rcpp_1.1.1            
#> [118] globals_0.19.1         spatstat.random_3.4-4  png_0.1-8             
#> [121] spatstat.univar_3.1-6  parallel_4.5.2         pkgdown_2.2.0         
#> [124] ggplot2_4.0.2          dotCall64_1.2          listenv_0.10.1        
#> [127] viridisLite_0.4.3      scales_1.4.0           ggridges_0.5.7        
#> [130] purrr_1.2.1            crayon_1.5.3           rlang_1.1.7           
#> [133] cowplot_1.2.0
```
