# Cached version of GuessDType for known string constants

Checks an internal cache before calling
[`GuessDType`](https://mianaz.github.io/scConvert/reference/GuessDType.md).
This avoids repeated H5T allocation for the same AnnData encoding
strings.

## Usage

``` r
CachedGuessDType(x, stype = "utf8")
```

## Arguments

- x:

  An R object to guess the HDF5 dtype for

- stype:

  String encoding type (default 'utf8')

## Value

An object of class
[`H5T`](http://hhoeflin.github.io/hdf5r/reference/H5T-class.md)
