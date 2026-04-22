# Add AnnData encoding attributes to an HDF5 dataset or group

Helper function to add standard AnnData encoding-type and
encoding-version attributes to an HDF5 object. This is used for string
arrays, categorical variables, and dataframes to ensure compatibility
with AnnData/scanpy.

## Usage

``` r
AddAnndataEncoding(
  h5obj,
  encoding_type = "string-array",
  encoding_version = "0.2.0"
)
```

## Arguments

- h5obj:

  An hdf5r H5D or H5Group object

- encoding_type:

  The encoding type (e.g., 'string-array', 'categorical', 'dataframe')

- encoding_version:

  The encoding version (default: '0.2.0')

## Value

NULL (modifies h5obj in place)
