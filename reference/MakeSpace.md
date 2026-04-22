# Make a space

Generate a blank space `n` characters long; useful for aligning text to
be printed to console

## Usage

``` r
MakeSpace(n)
```

## Arguments

- n:

  Length space should be

## Value

A space (' ') of length `n`

## Examples

``` r
# \donttest{
scConvert:::MakeSpace(n = 10)
#> [1] "          "
cat('hello', scConvert:::MakeSpace(n = 10), 'world\n', sep = '')
#> hello          world
# }
```
