# Translate cloud URL scheme to HTTPS base URL

s3://bucket/key -\> https://bucket.s3.amazonaws.com/key (virtual-hosted
style, anonymous public buckets only). gs://bucket/key -\>
https://storage.googleapis.com/bucket/key. http(s):// is returned as-is.

## Usage

``` r
.zarr_translate_url(url)
```

## Arguments

- url:

  A character URL

## Value

A list with \`base\` (https URL, trailing slash stripped) and
\`list_url\` (for ListObjectsV2 enumeration on S3/GCS) and \`kind\`
("s3", "gs", "http").
