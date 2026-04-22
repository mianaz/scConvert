# Load a MuData H5MU file as a Seurat object

Read multimodal MuData (.h5mu) files and convert to a Seurat object with
multiple assays. Uses native HDF5 reading with no external dependencies
(no MuDataSeurat or Python required).

## Usage

``` r
readH5MU(
  file,
  modalities = NULL,
  assay.names = NULL,
  restore.spatial = TRUE,
  verbose = TRUE
)
```

## Arguments

- file:

  Path to .h5mu file

- modalities:

  Character vector of modality names to load. If NULL (default), loads
  all modalities found in the file

- assay.names:

  Named vector mapping modality names to desired Seurat assay names. If
  NULL, uses standard mapping (rna-\>RNA, prot-\>ADT, atac-\>ATAC, etc.)

- restore.spatial:

  Logical; if TRUE, attempts to restore spatial data from the h5mu file
  structure (coordinates, images, scalefactors)

- verbose:

  Show progress messages

## Value

A `Seurat` object with multiple assays corresponding to each modality

## Details

The h5mu format stores multimodal data where each modality is stored as
a separate AnnData-like structure under `/mod/{modality_name}`. This
function:

- Reads each modality's expression matrix (X) and layers

- Creates Seurat assays with appropriate names

- Preserves global cell metadata from `/obs`

- Reads per-modality embeddings (obsm) and graphs (obsp)

- Restores spatial data if present

## Modality Mapping

By default, modality names are mapped to standard Seurat assay names:

- rna -\> RNA

- prot -\> ADT

- atac -\> ATAC

- spatial -\> Spatial

- Other names are preserved as-is

## Examples

``` r
if (FALSE) { # \dontrun{
# Load all modalities from an h5mu file
seurat_obj <- readH5MU("multimodal_data.h5mu")

# Load specific modalities only
seurat_obj <- readH5MU("multimodal_data.h5mu", modalities = c("rna", "prot"))

# Custom assay name mapping
seurat_obj <- readH5MU(
  "data.h5mu",
  assay.names = c(rna = "RNA", prot = "Protein")
)
} # }
```
