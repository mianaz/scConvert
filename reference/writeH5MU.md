# Write a Seurat object to h5mu format

Writes a multi-assay Seurat object to an h5mu file, with each assay
becoming a separate modality under `/mod/{modality}/`. Works with Seurat
V5 Assay5 objects. No external dependencies required.

## Usage

``` r
writeH5MU(
  object,
  filename = NULL,
  assays = NULL,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object with one or more assays

- filename:

  Output h5mu file path. If NULL, derived from project name.

- assays:

  Character vector of assay names to export (default: all)

- overwrite:

  Overwrite existing file (default: FALSE)

- verbose:

  Print progress messages (default: TRUE)

## Value

Invisibly returns the output filename

## Details

The h5mu format is designed for multimodal data and stores each Seurat
assay as a separate modality under `/mod/{modality_name}`. This
function:

- Extracts each specified assay from the Seurat object

- Converts assays to modality structure

- Writes counts and data layers for each modality

- Preserves cell metadata in global /obs and per-modality obs

- Maintains dimensional reductions and graphs

## Assay Mapping

By default, Seurat assay names are mapped to standard MuData modality
names:

- RNA -\> rna

- ADT -\> prot

- ATAC -\> atac

- Spatial -\> spatial

- Other names are converted to lowercase

## See also

[`readH5MU`](https://mianaz.github.io/scConvert/reference/readH5MU.md),
[`scConvert`](https://mianaz.github.io/scConvert/reference/scConvert.md),
[`as.h5mu`](https://mianaz.github.io/scConvert/reference/as.h5mu.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Write a multi-assay Seurat object (e.g. CITE-seq with RNA + ADT)
writeH5MU(seurat_obj, filename = "multimodal.h5mu")

# Write specific assays only
writeH5MU(seurat_obj, filename = "rna_only.h5mu", assays = "RNA")
} # }
```
