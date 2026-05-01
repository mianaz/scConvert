#' @include Convert.R
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CLI-accelerated conversion helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Find the scconvert CLI binary
#'
#' Searches for the compiled scconvert binary in the package source tree
#' and on the system PATH.
#'
#' @return Path to the scconvert binary, or NULL if not found
#'
#' @keywords internal
#'
sc_find_cli <- function() {
  # 0. Explicit user override
  override <- getOption("scConvert.cli_path")
  if (!is.null(override) && nzchar(override) && file.exists(override)) {
    return(override)
  }

  # 1. Installed package inst/bin (becomes bin/ in installed tree)
  pkg_bin <- system.file("bin", "scconvert", package = "scConvert")
  if (nzchar(pkg_bin) && file.exists(pkg_bin)) {
    return(pkg_bin)
  }

  # 2. Source tree src/ (devtools::load_all during development)
  pkg_src <- system.file("src", "scconvert", package = "scConvert")
  if (nzchar(pkg_src) && file.exists(pkg_src)) {
    return(pkg_src)
  }

  # 3. Relative to package root (another load_all shape)
  pkg_dir <- find.package("scConvert", quiet = TRUE)
  if (length(pkg_dir) > 0) {
    for (rel in c(file.path("src", "scconvert"),
                  file.path("bin", "scconvert"))) {
      p <- file.path(pkg_dir, rel)
      if (file.exists(p)) return(p)
    }
  }

  # 4. System PATH
  on_path <- Sys.which("scconvert")
  if (nzchar(on_path)) {
    return(unname(on_path))
  }

  return(NULL)
}

#' Run the scconvert CLI for file-to-file conversion
#'
#' Performs file-to-file format conversion with three tiers:
#' \enumerate{
#'   \item \strong{C binary}: HDF5 format pairs (h5ad, h5seurat, h5mu)
#'   \item \strong{Streaming}: zarr pairs (h5ad/h5seurat \eqn{\leftrightarrow}{<->} zarr) -- copies
#'     fields directly without creating a Seurat intermediate
#'   \item \strong{R hub}: all other format pairs via \code{scConvert()}
#' }
#'
#' Supported formats: h5ad, h5seurat, h5mu, loom, rds, zarr.
#'
#' @param input Input file path
#' @param output Output file path
#' @param assay Assay name (passed as --assay)
#' @param gzip Gzip compression level (0-9, passed as --gzip)
#' @param overwrite If TRUE, remove output file before conversion
#' @param verbose If TRUE, print CLI output
#'
#' @return TRUE on success, FALSE on failure
#'
#' @export
#'
scConvert_cli <- function(
  input,
  output,
  assay = "RNA",
  gzip = 4L,
  overwrite = FALSE,
  verbose = TRUE
) {
  is_uri <- grepl("^[a-zA-Z][a-zA-Z0-9+.-]*://", input)
  if (!is_uri && !file.exists(input) && !dir.exists(input)) {
    stop("Input file not found: ", input, call. = FALSE)
  }

  stype <- FileType(file = input)
  dtype <- FileType(file = output)

  # Tier 1: Try C binary first for all supported HDF5/loom/zarr format pairs
  cli_formats <- c("h5ad", "h5seurat", "h5mu", "loom", "zarr")
  use_c_cli <- stype %in% cli_formats && dtype %in% cli_formats

  if (use_c_cli) {
    cli_bin <- sc_find_cli()
    if (!is.null(cli_bin)) {
      if (file.exists(output) || dir.exists(output)) {
        if (overwrite) {
          unlink(output, recursive = TRUE)
        } else {
          stop("Output exists: ", output, ". Use overwrite = TRUE.", call. = FALSE)
        }
      }

      args <- c(input, output, "--assay", assay, "--gzip", as.character(gzip))
      if (verbose) message("Using scconvert CLI for fast on-disk conversion")

      result <- tryCatch({
        ret <- system2(cli_bin, args = args, stdout = TRUE, stderr = TRUE)
        status <- attr(ret, "status")
        if (!is.null(status) && status != 0L) {
          if (verbose) {
            message("CLI conversion failed (exit ", status, "): ",
                    paste(ret, collapse = "\n"))
          }
          FALSE
        } else {
          if (verbose && length(ret) > 0) message(paste(ret, collapse = "\n"))
          file.exists(output)
        }
      }, error = function(e) {
        if (verbose) message("CLI binary error: ", e$message)
        FALSE
      })

      if (!result && (file.exists(output) || dir.exists(output))) {
        try(unlink(output, recursive = TRUE), silent = TRUE)
      }
      if (result) return(TRUE)
      if (verbose) message("C binary failed, falling back to R streaming")
    }
  }

  # Tier 2: R streaming converters (zarr pairs + h5mu/loom pairs)
  streaming_converters <- list(
    # zarr pairs
    "h5ad|zarr"      = function() H5ADToZarr(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "zarr|h5ad"      = function() ZarrToH5AD(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5seurat|zarr"  = function() H5SeuratToZarr(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "zarr|h5seurat"  = function() ZarrToH5Seurat(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # h5mu <-> zarr
    "h5mu|zarr"      = function() H5MUToZarr(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "zarr|h5mu"      = function() ZarrToH5MU(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # loom <-> zarr
    "loom|zarr"      = function() LoomToZarr(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "zarr|loom"      = function() ZarrToLoom(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # h5mu <-> h5seurat (R streaming)
    "h5mu|h5seurat"  = function() H5MUToH5Seurat(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5seurat|h5mu"  = function() H5SeuratToH5MU(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # loom <-> h5seurat (R streaming)
    "loom|h5seurat"  = function() LoomToH5Seurat(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5seurat|loom"  = function() H5SeuratToLoom(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # h5mu <-> h5ad (R streaming via h5seurat)
    "h5mu|h5ad"      = function() H5MUToH5AD(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5ad|h5mu"      = function() H5ADToH5MU(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # loom <-> h5ad (R streaming via h5seurat)
    "loom|h5ad"      = function() LoomToH5AD(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5ad|loom"      = function() H5ADToLoom(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose),
    # loom <-> h5mu (R streaming via h5seurat)
    "loom|h5mu"      = function() LoomToH5MU(input, output, assay = assay, overwrite = overwrite, gzip = gzip, verbose = verbose),
    "h5mu|loom"      = function() H5MUToLoom(input, output, overwrite = overwrite, gzip = gzip, verbose = verbose)
  )

  pair_key <- paste0(stype, "|", dtype)
  if (pair_key %in% names(streaming_converters)) {
    result <- tryCatch({
      streaming_converters[[pair_key]]()
      TRUE
    }, error = function(e) {
      if (verbose) message("R streaming failed: ", e$message,
                           " -- falling back to R hub")
      FALSE
    })
    if (result) return(TRUE)
  }

  # R-level conversion: supports all format pairs
  if (verbose) message("Converting ", stype, " -> ", dtype, " via R")
  tryCatch({
    scConvert(source = input, dest = output, assay = assay,
            overwrite = overwrite, verbose = verbose)
    TRUE
  }, error = function(e) {
    if (verbose) message("Conversion failed: ", e$message)
    FALSE
  })
}
