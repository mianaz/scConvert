#' @include AnnDataEncoding.R
#' @include ZarrRemote.R
NULL

#' Storage abstraction for zarr v2 stores
#'
#' Returns a list of method closures that hide the difference between a
#' local-filesystem store and a remote HTTP-backed store. Every \code{.zarr_*}
#' read helper that previously took a \code{store_path} string now also
#' accepts the list returned here.
#'
#' Methods:
#' \describe{
#'   \item{\code{kind}}{\code{"local"} or \code{"http"}.}
#'   \item{\code{root}}{The store root URI (filesystem path or URL).}
#'   \item{\code{exists(rel)}}{Logical: does \code{rel} exist in the store?}
#'   \item{\code{read_bytes(rel)}}{Raw vector with the contents of \code{rel}.
#'     Errors on missing keys; callers should \code{store$exists(rel)} first
#'     if a missing key is expected.}
#'   \item{\code{read_json(rel)}}{Parsed JSON list, or empty list if missing.}
#'   \item{\code{list_dirs(rel)}}{Character vector of immediate subgroup
#'     names under \code{rel} (zarr-meaningful only; filters \code{c/},
#'     dotfiles, and \code{__*}).}
#' }
#'
#' @param x Either a filesystem path, an \code{s3://}/\code{gs://}/HTTP URL,
#'   or an already-constructed store list.
#' @param cache For HTTP stores: \code{TRUE} to persist fetched objects
#'   under \code{tools::R_user_dir("scConvert", "cache")}, \code{FALSE} for
#'   a tempdir that is discarded with the session, \code{NA} for default
#'   (TRUE).
#'
#' @return A store list.
#'
#' @keywords internal
#'
.zarr_make_store <- function(x, cache = NA) {
  if (is.list(x) && !is.null(x$kind)) return(x)
  if (!is.character(x) || length(x) != 1L) {
    stop("zarr store input must be a path, URL, or store list", call. = FALSE)
  }
  if (.is_remote_zarr_url(x)) {
    cache_dir <- if (isFALSE(cache)) NULL else .scconvert_cache_dir()
    return(.zarr_store_http(x, cache_dir = cache_dir))
  }
  .zarr_store_local(x)
}

#' Local filesystem zarr store
#'
#' @keywords internal
.zarr_store_local <- function(root) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  abs <- function(rel) {
    if (!nzchar(rel)) return(root)
    file.path(root, rel)
  }
  list(
    kind = "local",
    root = root,
    exists = function(rel) {
      p <- abs(rel)
      file.exists(p) || dir.exists(p)
    },
    read_bytes = function(rel) {
      p <- abs(rel)
      if (!file.exists(p)) {
        stop("zarr: missing key '", rel, "' in ", root, call. = FALSE)
      }
      readBin(p, "raw", file.info(p)$size)
    },
    read_json = function(rel) {
      p <- abs(rel)
      if (!file.exists(p)) return(list())
      jsonlite::fromJSON(p, simplifyVector = TRUE)
    },
    list_dirs = function(rel) {
      p <- abs(rel)
      if (!dir.exists(p)) return(character(0))
      entries <- list.dirs(p, full.names = FALSE, recursive = FALSE)
      entries[!grepl("^\\.|^c$|^__", entries)]
    }
  )
}

#' HTTP-backed zarr store with per-chunk lazy fetch
#'
#' Lazy: each \code{read_bytes()} issues a fresh HTTP GET (cached on disk
#' if \code{cache_dir} is non-NULL). \code{list_dirs()} uses S3
#' ListObjectsV2 with the slash delimiter. Once a store is opened, the
#' full object listing is fetched up front and memoised; subsequent
#' \code{exists()} / \code{list_dirs()} checks hit the memoised manifest.
#'
#' @keywords internal
.zarr_store_http <- function(url, cache_dir = NULL) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("HTTP zarr store requires the 'httr' package. ",
         "Install with: install.packages('httr')", call. = FALSE)
  }
  info <- .zarr_translate_url(url)
  if (info$kind == "http") {
    stop("Lazy-fetch HTTP zarr store requires s3:// or gs:// URL ",
         "(needed for ListObjectsV2). Got a raw HTTP URL: ", url,
         call. = FALSE)
  }

  # Fetch the full key manifest once up front (single S3 ListObjects walk).
  # This is O(1) HTTP roundtrips per ~1000 keys; tiny vs chunk fetches.
  keys <- .zarr_list_remote_keys(info)
  manifest <- setNames(rep(TRUE, length(keys)), keys)

  # Cache directory: keyed by URL hash so multiple URLs coexist.
  cache_root <- NULL
  if (!is.null(cache_dir)) {
    hash <- substr(.scconvert_hash_url(url), 1L, 16L)
    cache_root <- file.path(cache_dir, hash)
    dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  }

  abs_url <- function(rel) {
    if (!nzchar(rel)) return(info$base)
    paste0(info$base, "/", rel)
  }
  cache_path <- function(rel) {
    if (is.null(cache_root)) return(NULL)
    file.path(cache_root, rel)
  }

  fetch <- function(rel) {
    cp <- cache_path(rel)
    if (!is.null(cp) && file.exists(cp)) {
      return(readBin(cp, "raw", file.info(cp)$size))
    }
    resp <- httr::GET(abs_url(rel), httr::timeout(120))
    if (httr::http_error(resp)) {
      stop("zarr: HTTP ", httr::status_code(resp), " fetching ",
           abs_url(rel), call. = FALSE)
    }
    bytes <- httr::content(resp, as = "raw")
    if (!is.null(cp)) {
      parent <- dirname(cp)
      if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
      writeBin(bytes, cp)
    }
    bytes
  }

  list(
    kind = "http",
    root = info$base,
    exists = function(rel) {
      # Membership check via manifest; cheap, no roundtrip.
      if (!nzchar(rel)) return(TRUE)
      if (isTRUE(manifest[rel])) return(TRUE)
      # Also TRUE if any key sits under rel (i.e., rel is a "directory").
      prefix <- paste0(rel, "/")
      any(startsWith(names(manifest), prefix))
    },
    read_bytes = function(rel) {
      if (!nzchar(rel)) stop("zarr: cannot read root", call. = FALSE)
      fetch(rel)
    },
    read_json = function(rel) {
      if (!isTRUE(manifest[rel])) return(list())
      bytes <- fetch(rel)
      jsonlite::fromJSON(rawToChar(bytes), simplifyVector = TRUE)
    },
    list_dirs = function(rel) {
      prefix <- if (nzchar(rel)) paste0(rel, "/") else ""
      under <- names(manifest)[startsWith(names(manifest), prefix)]
      if (length(under) == 0L) return(character(0))
      under <- substr(under, nchar(prefix) + 1L, .Machine$integer.max)
      # Immediate-child segments only.
      segs <- vapply(strsplit(under, "/", fixed = TRUE), `[[`, character(1L), 1L)
      segs <- unique(segs)
      segs[!grepl("^\\.|^c$|^__", segs)]
    }
  )
}
