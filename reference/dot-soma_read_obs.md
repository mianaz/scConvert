# Read obs from a SOMA experiment

Read obs from a SOMA experiment

## Usage

``` r
.soma_read_obs(exp, obs_query = NULL, column_names = NULL)
```

## Arguments

- exp:

  A SOMAExperiment object

- obs_query:

  Optional value filter string

- column_names:

  Optional column names to read

## Value

A data.frame with cell metadata, rownames set to cell identifiers
