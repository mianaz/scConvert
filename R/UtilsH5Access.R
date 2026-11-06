#' HDF5 Access Utility Functions
#'
#' This file contains utility functions for safely accessing HDF5 files,
#' including existence checking with caching, dataset reading with V5
#' compatibility, and feature retrieval.
#'
#' @name UtilsH5Access
#' @keywords internal
NULL

#' Create a cached HDF5 existence checker
#'
#' Creates a closure that checks for HDF5 path existence with result caching
#' to improve performance for repeated checks.
#'
#' @return A function with signature \code{function(group, path)} that checks
#'   if \code{path} exists in \code{group}, using a cache to avoid repeated
#'   HDF5 calls.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' safe_exists <- CreateCachedExistsChecker()
#' if (safe_exists(h5_group, "data/matrix")) {
#'   # path exists
#' }
#' }
#'
CreateCachedExistsChecker <- function() {
  cache <- new.env(hash = TRUE, parent = emptyenv())

  function(group, path) {
    cache_key <- paste0(group$get_obj_name(), "::", path)

    if (exists(cache_key, envir = cache, inherits = FALSE)) {
      return(get(cache_key, envir = cache))
    }

    result <- tryCatch(
      expr = {
        group$exists(name = path)
      },
      error = function(e) {
        path %in% names(group)
      }
    )

    assign(cache_key, result, envir = cache)
    return(result)
  }
}

#' Get feature names with V5 fallback
#'
#' Retrieves feature names from an HDF5 group, trying V5 location first
#' (meta.data/_index), then falling back to direct features dataset.
#'
#' @param h5_group HDF5 group containing feature information
#' @param verbose Print diagnostic messages
#'
#' @return Character vector of feature names, or NULL if not found
#' @keywords internal
#'
GetFeaturesV5Safe <- function(h5_group, verbose = FALSE) {
  v5_index_exists <- tryCatch(
    h5_group$exists("meta.data/_index"),
    error = function(e) FALSE
  )

  if (v5_index_exists) {
    if (verbose) {
      message("Using V5 feature index from meta.data/_index")
    }
    return(as.character(h5_group[["meta.data/_index"]][]))
  }

  features_exists <- tryCatch(
    h5_group$exists("features"),
    error = function(e) FALSE
  )

  if (features_exists) {
    if (verbose) {
      message("Using direct features dataset")
    }
    result <- tryCatch(
      expr = as.character(h5_group[["features"]][]),
      error = function(e) {
        if (length(h5_group[["features"]]$dims) > 1) {
          return(as.character(h5_group[["features"]][, 1]))
        }
        return(NULL)
      }
    )
    return(result)
  }

  if (verbose) {
    message("No features found in h5_group")
  }
  return(NULL)
}


