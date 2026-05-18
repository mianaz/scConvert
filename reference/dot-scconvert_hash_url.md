# Stable hash of a URL string for cache directory naming

Uses base R \`digest\`-free MD5-like via \`tools::md5sum\` on a
tempfile. Falls back to a CRC32-ish hex of the byte sum if md5sum is
unavailable.

## Usage

``` r
.scconvert_hash_url(url)
```
