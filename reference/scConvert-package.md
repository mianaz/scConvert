# scConvert: Universal Single-Cell Format Conversion for R

Universal single-cell format conversion for R. Converts between h5ad
(AnnData), h5Seurat, h5mu (MuData), Loom, Zarr, RDS, and
SingleCellExperiment formats through a hub architecture using Seurat as
an intermediate. No Python dependency required.

## Package options

scConvert uses the following options to control behavior, users can
configure these with [`options`](https://rdrr.io/r/base/options.html):

- `scConvert.dtypes.logical_to_int`:

  When writing [logical](https://rdrr.io/r/base/logical.html) vectors,
  coerce to integer types to ensure compatibility across languages (see
  [`BoolToInt`](https://mianaz.github.io/scConvert/reference/BoolToInt.md)
  for more details)

- `scConvert.dtypes.dataframe_as_group`:

  When writing [data.frame](https://rdrr.io/r/base/data.frame.html)s,
  always write out as a group regardless of factor presence

- `scConvert.chunking.MARGIN`:

  Default direction for chunking datasets; choose from:

  largest

  :   Chunk along the largest dimension of a dataset

  smallest

  :   Chunk along the smallest dimension

  first

  :   Chunk along the first dimension

  last

  :   Chunk along the last dimension

- `scConvert.dimreducs.allglobal`:

  Treat all DimReducs as global, regardless of actual global status

- `scConvert.compression.level`:

  Gzip compression level for HDF5 dataset writes (0-9). Level 0 means no
  compression, level 4 (default) is a good balance of speed vs size.
  Only applies to newly created datasets, not to data copied via
  `obj_copy_from`.

## See also

Useful links:

- <https://mianaz.github.io/scConvert/>

- <https://github.com/mianaz/scConvert>

- Report bugs at <https://github.com/mianaz/scConvert/issues>

## Author

**Maintainer**: Ziyu Zeng <zzeng4@nd.edu>
([ORCID](https://orcid.org/0009-0005-4055-8492)) \[copyright holder\]
