#' @include AnnDataEncoding.R
NULL

#' Detect remote zarr URL
#'
#' @param x A character string (path or URL)
#'
#' @return TRUE if x looks like a remote zarr URL (s3://, gs://, http(s)://)
#'
#' @keywords internal
#'
.is_remote_zarr_url <- function(x) {
  if (!is.character(x) || length(x) != 1L) return(FALSE)
  grepl("^(s3|gs|https?)://", x, ignore.case = TRUE)
}

#' Translate cloud URL scheme to HTTPS base URL
#'
#' s3://bucket/key -> https://bucket.s3.amazonaws.com/key (virtual-hosted style,
#' anonymous public buckets only).
#' gs://bucket/key -> https://storage.googleapis.com/bucket/key.
#' http(s):// is returned as-is.
#'
#' @param url A character URL
#'
#' @return A list with `base` (https URL, trailing slash stripped) and
#'   `list_url` (for ListObjectsV2 enumeration on S3/GCS) and `kind`
#'   ("s3", "gs", "http").
#'
#' @keywords internal
#'
.zarr_translate_url <- function(url) {
  url <- sub("/+$", "", url)
  if (grepl("^s3://", url, ignore.case = TRUE)) {
    rest <- sub("^s3://", "", url, ignore.case = TRUE)
    parts <- strsplit(rest, "/", fixed = TRUE)[[1]]
    bucket <- parts[1L]
    key <- if (length(parts) > 1L) paste(parts[-1L], collapse = "/") else ""
    base <- paste0("https://", bucket, ".s3.amazonaws.com")
    full <- if (nzchar(key)) paste0(base, "/", key) else base
    return(list(kind = "s3", base = full,
                list_url = base, prefix = key, bucket = bucket))
  }
  if (grepl("^gs://", url, ignore.case = TRUE)) {
    rest <- sub("^gs://", "", url, ignore.case = TRUE)
    parts <- strsplit(rest, "/", fixed = TRUE)[[1]]
    bucket <- parts[1L]
    key <- if (length(parts) > 1L) paste(parts[-1L], collapse = "/") else ""
    base <- paste0("https://storage.googleapis.com/", bucket)
    full <- if (nzchar(key)) paste0(base, "/", key) else base
    return(list(kind = "gs", base = full,
                list_url = base, prefix = key, bucket = bucket))
  }
  list(kind = "http", base = url, list_url = NA_character_,
       prefix = NA_character_, bucket = NA_character_)
}

#' List all object keys under a prefix in an S3/GCS bucket
#'
#' Uses ListObjectsV2 (S3) or list_v1 (GCS) over anonymous HTTP. Paginates
#' via continuation tokens. Returns keys relative to the prefix.
#'
#' @param info Output of `.zarr_translate_url()`.
#'
#' @return Character vector of keys relative to the prefix (no leading
#'   slash). Empty character if listing failed.
#'
#' @keywords internal
#'
.zarr_list_remote_keys <- function(info) {
  if (info$kind == "http") {
    stop("Remote zarr listing requires s3:// or gs:// URL; got a raw HTTP URL. ",
         "scConvert does not enumerate arbitrary HTTP servers.", call. = FALSE)
  }
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Reading remote zarr requires the 'httr' package. ",
         "Install with: install.packages('httr')", call. = FALSE)
  }
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Reading remote zarr requires the 'xml2' package. ",
         "Install with: install.packages('xml2')", call. = FALSE)
  }

  keys <- character(0)
  prefix <- if (nzchar(info$prefix)) paste0(info$prefix, "/") else ""
  token <- NULL

  repeat {
    query <- list(`list-type` = "2", prefix = prefix)
    if (!is.null(token)) query[["continuation-token"]] <- token
    resp <- httr::GET(info$list_url, query = query, httr::timeout(60))
    if (httr::http_error(resp)) {
      stop("Failed to list bucket ", info$bucket, " (HTTP ",
           httr::status_code(resp), "). ",
           "For private buckets you must sign requests yourself; ",
           "scConvert supports anonymous public buckets only.", call. = FALSE)
    }
    doc <- xml2::read_xml(httr::content(resp, "raw"))
    ns <- xml2::xml_ns(doc)
    key_nodes <- xml2::xml_find_all(doc, ".//d1:Contents/d1:Key", ns)
    page_keys <- xml2::xml_text(key_nodes)
    if (nzchar(prefix)) page_keys <- sub(paste0("^", prefix), "", page_keys, fixed = FALSE)
    keys <- c(keys, page_keys)
    truncated <- xml2::xml_find_first(doc, ".//d1:IsTruncated", ns)
    if (length(truncated) == 0 || !identical(tolower(xml2::xml_text(truncated)), "true")) break
    next_token <- xml2::xml_find_first(doc, ".//d1:NextContinuationToken", ns)
    if (length(next_token) == 0) break
    token <- xml2::xml_text(next_token)
  }
  keys[nzchar(keys)]
}

#' Mirror a remote zarr store to a local directory
#'
#' Fetches every object under the prefix in parallel-friendly serial GETs,
#' preserving relative paths. Uses an on-disk cache keyed by URL when
#' `cache_dir` is non-NULL: a subsequent call against the same URL returns
#' the cached directory without re-downloading.
#'
#' @param url Remote URL (s3://, gs://).
#' @param cache_dir Cache directory, or NULL for tempdir (no persistence).
#' @param verbose Print progress.
#'
#' @return Path to local directory containing the mirrored store.
#'
#' @keywords internal
#'
.zarr_fetch_remote <- function(url, cache_dir = NULL, verbose = TRUE) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Reading remote zarr requires the 'httr' package. ",
         "Install with: install.packages('httr')", call. = FALSE)
  }
  info <- .zarr_translate_url(url)

  hash <- substr(.scconvert_hash_url(url), 1L, 16L)
  if (is.null(cache_dir)) {
    dest <- tempfile(pattern = paste0("zarr_", hash, "_"))
    dir.create(dest, recursive = TRUE)
    use_cache <- FALSE
  } else {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    dest <- file.path(cache_dir, hash)
    use_cache <- TRUE
  }

  marker <- file.path(dest, ".scconvert_zarr_complete")
  if (use_cache && file.exists(marker)) {
    if (verbose) message("Using cached zarr store: ", dest)
    return(dest)
  }
  if (use_cache) {
    unlink(dest, recursive = TRUE, force = TRUE)
    dir.create(dest, recursive = TRUE)
  }

  if (verbose) message("Listing remote zarr keys: ", url)
  keys <- .zarr_list_remote_keys(info)
  if (length(keys) == 0L) {
    stop("Remote zarr store is empty or unreadable: ", url, call. = FALSE)
  }

  if (verbose) message("Downloading ", length(keys), " objects to ", dest)
  for (i in seq_along(keys)) {
    key <- keys[i]
    rel <- key
    local_path <- file.path(dest, rel)
    parent <- dirname(local_path)
    if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
    src_url <- paste0(info$base, "/", key)
    resp <- httr::GET(src_url, httr::write_disk(local_path, overwrite = TRUE),
                      httr::timeout(300))
    if (httr::http_error(resp)) {
      stop("Failed to fetch ", src_url, " (HTTP ",
           httr::status_code(resp), ")", call. = FALSE)
    }
    if (verbose && (i %% 50L == 0L || i == length(keys))) {
      message("  fetched ", i, "/", length(keys))
    }
  }
  if (use_cache) file.create(marker)
  dest
}

#' Stable hash of a URL string for cache directory naming
#'
#' Uses base R `digest`-free MD5-like via `tools::md5sum` on a tempfile.
#' Falls back to a CRC32-ish hex of the byte sum if md5sum is unavailable.
#'
#' @keywords internal
.scconvert_hash_url <- function(url) {
  tf <- tempfile()
  writeLines(url, tf)
  on.exit(unlink(tf))
  md5 <- tools::md5sum(tf)
  if (!is.na(md5)) return(unname(md5))
  format(as.hexmode(sum(as.integer(charToRaw(url)))), width = 8L)
}

#' User cache directory for scConvert
#'
#' @keywords internal
.scconvert_cache_dir <- function() {
  base <- tryCatch(tools::R_user_dir("scConvert", which = "cache"),
                   error = function(e) NULL)
  if (is.null(base)) base <- file.path(tempdir(), "scConvert-cache")
  base
}
