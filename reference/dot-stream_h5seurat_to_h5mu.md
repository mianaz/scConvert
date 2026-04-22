# Stream h5seurat to h5mu (direct HDF5-to-HDF5)

For each assay, creates a modality under /mod/ with the full AnnData
layout. Uses the zero-copy CSC\<-\>CSR reinterpretation trick.

## Usage

``` r
.stream_h5seurat_to_h5mu(h5seurat_path, h5mu_path, gzip = 4L, verbose = TRUE)
```

## Arguments

- h5seurat_path:

  Path to input h5seurat file

- h5mu_path:

  Path to output h5mu file

- gzip:

  Integer gzip compression level

- verbose:

  Logical
