# Read var from a SOMA measurement

Read var from a SOMA measurement

## Usage

``` r
.soma_read_var(ms, var_query = NULL, column_names = NULL)
```

## Arguments

- ms:

  A SOMAMeasurement object

- var_query:

  Optional value filter string

- column_names:

  Optional column names to read

## Value

A data.frame with feature metadata, rownames set to feature names
