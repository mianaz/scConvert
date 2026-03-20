#' @include zzz.R
#' @include AnnDataEncoding.R
#' @importFrom Matrix t
#' @importFrom Seurat GetAssayData DefaultAssay VariableFeatures Embeddings
#' @importFrom SeuratObject Cells Layers
NULL

#' Save a Seurat object as an AnnData Zarr store
#'
#' Writes a Seurat object to AnnData Zarr v2 format (.zarr directory).
#' The output is readable by Python \code{anndata.read_zarr()}, scanpy,
#' and any tool supporting the AnnData on-disk specification.
#'
#' @param object A Seurat object
#' @param filename Path for the output .zarr directory
#' @param assay Name of assay to export (default: \code{DefaultAssay(object)})
#' @param overwrite If TRUE, overwrite existing zarr store
#' @param verbose Show progress messages
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns \code{filename}
#'
#' @export
#'
writeZarr <- function(object, filename, assay = DefaultAssay(object),
                     overwrite = FALSE, verbose = TRUE, ...) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package required for writeZarr. ",
         "Install with: install.packages('jsonlite')", call. = FALSE)
  }

  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }

  if (dir.exists(filename) || file.exists(filename)) {
    if (!overwrite) {
      stop("Destination exists: ", filename,
           ". Use overwrite = TRUE to replace.", call. = FALSE)
    }
  }

  # Atomic write: write to temp dir, rename on success
  use_atomic <- dir.exists(filename) || file.exists(filename)
  final_dest <- filename
  if (use_atomic) {
    filename <- tempfile(tmpdir = dirname(final_dest), fileext = ".zarr")
  }

  if (verbose) message("Saving Seurat object to Zarr: ", final_dest)

  cell.names <- Cells(object)
  assay_obj <- object[[assay]]
  feature.names <- rownames(assay_obj)

  compressor <- list(id = "zlib", level = GetCompressionLevel())

  # 1. Create root group
  .zarr_create_group(filename, attrs = list(
    `encoding-type` = "anndata",
    `encoding-version` = "0.1.0"
  ))

  # 2. Write obs (cell metadata)
  if (verbose) message("Writing cell metadata (obs)...")
  .zarr_write_anndata_dataframe(
    store_path = filename,
    group_path = "obs",
    df = object@meta.data,
    index = cell.names,
    compressor = compressor,
    verbose = verbose
  )

  # 3. Write var (feature metadata)
  if (verbose) message("Writing feature metadata (var)...")
  var_df <- tryCatch(
    assay_obj[[]],
    error = function(e) data.frame(row.names = feature.names)
  )
  if (is.null(var_df) || !is.data.frame(var_df)) {
    var_df <- data.frame(row.names = feature.names)
  }
  .zarr_write_anndata_dataframe(
    store_path = filename,
    group_path = "var",
    df = var_df,
    index = feature.names,
    compressor = compressor,
    verbose = verbose
  )

  # 4. Write X (primary expression matrix)
  if (verbose) message("Writing expression matrix (X)...")
  x_data <- tryCatch(
    GetAssayData(object, assay = assay, layer = "counts"),
    error = function(e) {
      tryCatch(
        GetAssayData(object, assay = assay, layer = "data"),
        error = function(e2) NULL
      )
    }
  )
  # Convert BPCells if needed
  x_data <- ConvertBPCellsMatrix(x_data, verbose = verbose)

  if (!is.null(x_data) && !IsMatrixEmpty(x_data)) {
    .zarr_write_anndata_matrix(filename, "X", x_data, compressor, verbose)
  }

  # 5. Write layers
  layer_names <- tryCatch(Layers(assay_obj), error = function(e) character(0))
  # Remove the layer we already wrote as X
  x_layer <- if ("counts" %in% layer_names) "counts" else "data"
  extra_layers <- setdiff(layer_names, x_layer)

  if (length(extra_layers) > 0) {
    if (verbose) message("Writing layers...")
    .zarr_create_group(file.path(filename, "layers"))

    # Get X dimensions for layer validation
    x_dims <- dim(x_data)  # (genes, cells) in R; anndata stores as (cells, genes)

    for (ln in extra_layers) {
      if (verbose) message("  Writing layer: ", ln)
      tryCatch({
        layer_data <- GetAssayData(object, assay = assay, layer = ln)
        layer_data <- ConvertBPCellsMatrix(layer_data, verbose = verbose)
        if (!is.null(layer_data) && !IsMatrixEmpty(layer_data)) {
          # AnnData requires all layers to have the same shape as X.
          # scale.data only covers HVGs and will be smaller — skip it.
          if (!identical(dim(layer_data), x_dims)) {
            if (verbose) message("    Skipping layer '", ln,
                                 "': shape ", paste(dim(layer_data), collapse = "x"),
                                 " != X shape ", paste(x_dims, collapse = "x"))
            next
          }
          anndata_name <- SeuratLayerToAnnData(ln)
          .zarr_write_anndata_matrix(
            filename, file.path("layers", anndata_name),
            layer_data, compressor, verbose
          )
        }
      }, error = function(e) {
        if (verbose) warning("Could not write layer ", ln, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  } else {
    # Write empty layers group
    .zarr_create_group(file.path(filename, "layers"))
  }

  # 6. Write obsm (dimensional reductions)
  if (length(object@reductions) > 0) {
    if (verbose) message("Writing dimensional reductions (obsm)...")
    .zarr_create_group(file.path(filename, "obsm"))

    for (reduc_name in names(object@reductions)) {
      if (verbose) message("  Writing reduction: ", reduc_name)
      tryCatch({
        embeddings <- Embeddings(object, reduction = reduc_name)
        anndata_name <- paste0("X_", reduc_name)
        .zarr_write_numeric(
          dir = file.path(filename, "obsm", anndata_name),
          data = embeddings,
          dtype = "<f8",
          compressor = compressor
        )
      }, error = function(e) {
        if (verbose) warning("Could not write reduction ", reduc_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  } else {
    .zarr_create_group(file.path(filename, "obsm"))
  }

  # 7. Write obsp (neighbor graphs)
  if (length(object@graphs) > 0) {
    if (verbose) message("Writing neighbor graphs (obsp)...")
    .zarr_create_group(file.path(filename, "obsp"))

    for (graph_name in names(object@graphs)) {
      if (verbose) message("  Writing graph: ", graph_name)
      tryCatch({
        graph_mat <- object@graphs[[graph_name]]
        # Map Seurat graph names to AnnData conventions
        anndata_name <- graph_name
        if (grepl("_snn$", graph_name)) anndata_name <- "connectivities"
        if (grepl("_nn$", graph_name)) anndata_name <- "distances"

        .zarr_write_anndata_matrix(
          filename, file.path("obsp", anndata_name),
          as(graph_mat, "dgCMatrix"), compressor, verbose,
          transpose = FALSE
        )
      }, error = function(e) {
        if (verbose) warning("Could not write graph ", graph_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  } else {
    .zarr_create_group(file.path(filename, "obsp"))
  }

  # 8. Write uns (unstructured metadata)
  if (length(object@misc) > 0) {
    if (verbose) message("Writing unstructured data (uns)...")
    .zarr_create_group(file.path(filename, "uns"))

    for (item_name in names(object@misc)) {
      tryCatch({
        item <- object@misc[[item_name]]
        item_path <- file.path("uns", item_name)

        if (is.character(item)) {
          .zarr_write_strings(
            dir = file.path(filename, item_path),
            strings = item,
            compressor = compressor
          )
        } else if (is.numeric(item) || is.logical(item)) {
          .zarr_write_numeric(
            dir = file.path(filename, item_path),
            data = item,
            compressor = compressor
          )
        }
        # Skip complex objects silently
      }, error = function(e) {
        if (verbose) warning("Could not write uns item ", item_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  } else {
    .zarr_create_group(file.path(filename, "uns"))
  }

  # Atomic rename: replace original only after successful write
  if (use_atomic) {
    unlink(final_dest, recursive = TRUE)
    file.rename(filename, final_dest)
    filename <- final_dest
  }

  if (verbose) {
    message("\nSuccessfully saved Zarr store")
    message("  Cells: ", length(cell.names))
    message("  Features: ", length(feature.names))
    message("  Path: ", filename)
  }

  invisible(filename)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal AnnData Zarr Write Helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write an AnnData DataFrame (obs or var) to zarr store
#'
#' @param store_path Path to zarr store
#' @param group_path Relative path for the group (e.g., "obs", "var")
#' @param df Data frame to write
#' @param index Character vector of index values (cell/gene names)
#' @param compressor Compressor spec
#' @param verbose Show progress
#'
#' @keywords internal
#'
.zarr_write_anndata_dataframe <- function(store_path, group_path, df, index,
                                           compressor, verbose = TRUE) {
  col_names <- colnames(df)

  # Create group with DataFrame encoding
  .zarr_create_group(
    file.path(store_path, group_path),
    attrs = list(
      `encoding-type` = "dataframe",
      `encoding-version` = "0.2.0",
      `_index` = "_index",
      `column-order` = col_names
    )
  )

  # Write _index (string array of row names)
  .zarr_write_strings(
    dir = file.path(store_path, group_path, "_index"),
    strings = as.character(index),
    compressor = compressor,
    attrs = list(
      `encoding-type` = "string-array",
      `encoding-version` = "0.2.0"
    )
  )

  # Write each column
  for (col in col_names) {
    col_data <- df[[col]]
    col_dir <- file.path(store_path, group_path, col)

    if (is.factor(col_data)) {
      # Write as categorical
      encoded <- EncodeCategorical(col_data)
      .zarr_create_group(col_dir, attrs = list(
        `encoding-type` = "categorical",
        `encoding-version` = "0.2.0",
        ordered = FALSE
      ))
      .zarr_write_numeric(
        dir = file.path(col_dir, "codes"),
        data = as.integer(encoded$codes),
        dtype = "|i1",
        compressor = compressor
      )
      .zarr_write_strings(
        dir = file.path(col_dir, "categories"),
        strings = encoded$categories,
        compressor = compressor,
        attrs = list(
          `encoding-type` = "string-array",
          `encoding-version` = "0.2.0"
        )
      )
    } else if (is.character(col_data)) {
      # String array
      .zarr_write_strings(
        dir = col_dir,
        strings = col_data,
        compressor = compressor,
        attrs = list(
          `encoding-type` = "string-array",
          `encoding-version` = "0.2.0"
        )
      )
    } else if (is.logical(col_data)) {
      # Boolean as int8
      .zarr_write_numeric(
        dir = col_dir,
        data = as.integer(col_data),
        dtype = "|b1",
        compressor = compressor
      )
    } else if (is.integer(col_data)) {
      .zarr_write_numeric(
        dir = col_dir,
        data = col_data,
        dtype = "<i4",
        compressor = compressor
      )
    } else if (is.numeric(col_data)) {
      .zarr_write_numeric(
        dir = col_dir,
        data = col_data,
        dtype = "<f8",
        compressor = compressor
      )
    }
  }
}

#' Write an AnnData matrix (sparse or dense) to zarr store
#'
#' @param store_path Path to zarr store
#' @param rel_path Relative path within store
#' @param mat Matrix (dgCMatrix or dense)
#' @param compressor Compressor spec
#' @param verbose Show progress
#' @param transpose If TRUE, transpose genes x cells to cells x genes for AnnData
#'
#' @keywords internal
#' @noRd
NULL

#' Write multiple zarr arrays in parallel
#'
#' Uses mclapply on Unix for parallel compression, falls back to lapply on Windows.
#'
#' @param write_specs List of list(dir, data, dtype) specs
#' @param compressor Compressor spec
#'
#' @keywords internal
#'
.zarr_parallel_write <- function(write_specs, compressor) {
  par_fn <- if (.Platform$OS.type == "unix" && length(write_specs) > 1) {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = min(length(X), 2L))
  } else {
    lapply
  }
  par_fn(write_specs, function(spec) {
    .zarr_write_numeric(
      dir = spec$dir,
      data = spec$data,
      dtype = spec$dtype,
      compressor = compressor
    )
  })
  invisible(NULL)
}

.zarr_write_anndata_matrix <- function(store_path, rel_path, mat, compressor,
                                        verbose = TRUE, transpose = TRUE) {
  # Determine if sparse
  is_sparse <- inherits(mat, c("dgCMatrix", "dgTMatrix", "dgRMatrix"))

  if (is_sparse) {
    # Write as CSR sparse matrix
    mat <- as(mat, "dgCMatrix")
    csr <- DeconstructSparseCSR(mat)

    if (!transpose) {
      # For obsp graphs: already in correct orientation, just extract CSC as-is
      csr <- list(
        data = mat@x,
        indices = mat@i,     # 0-based row indices for CSC
        indptr = mat@p,
        shape = as.integer(dim(mat))
      )
    }

    .zarr_create_group(
      file.path(store_path, rel_path),
      attrs = list(
        `encoding-type` = if (transpose) "csr_matrix" else "csc_matrix",
        `encoding-version` = "0.1.0",
        shape = as.list(as.integer(csr$shape))
      )
    )
    .zarr_parallel_write(
      write_specs = list(
        list(dir = file.path(store_path, rel_path, "data"),
             data = csr$data, dtype = "<f8"),
        list(dir = file.path(store_path, rel_path, "indices"),
             data = as.integer(csr$indices), dtype = "<i4"),
        list(dir = file.path(store_path, rel_path, "indptr"),
             data = as.integer(csr$indptr), dtype = "<i4")
      ),
      compressor = compressor
    )
  } else {
    # Dense matrix
    if (transpose) {
      # genes x cells -> cells x genes
      mat <- t(mat)
    }
    .zarr_write_numeric(
      dir = file.path(store_path, rel_path),
      data = mat,
      dtype = "<f8",
      compressor = compressor,
      attrs = list(
        `encoding-type` = "array",
        `encoding-version` = "0.2.0"
      )
    )
  }
}
