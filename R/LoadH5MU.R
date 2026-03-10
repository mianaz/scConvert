#' Load a MuData H5MU file as a Seurat object
#'
#' Read multimodal MuData (.h5mu) files and convert to a Seurat object with
#' multiple assays. Uses native HDF5 reading with no external dependencies
#' (no MuDataSeurat or Python required).
#'
#' @param file Path to .h5mu file
#' @param modalities Character vector of modality names to load. If NULL (default),
#'   loads all modalities found in the file
#' @param assay.names Named vector mapping modality names to desired Seurat assay
#'   names. If NULL, uses standard mapping (rna->RNA, prot->ADT, atac->ATAC, etc.)
#' @param restore.spatial Logical; if TRUE, attempts to restore spatial data from
#'   the h5mu file structure (coordinates, images, scalefactors)
#' @param verbose Show progress messages
#'
#' @return A \code{Seurat} object with multiple assays corresponding to each modality
#'
#' @details
#' The h5mu format stores multimodal data where each modality is stored as a
#' separate AnnData-like structure under \code{/mod/{modality_name}}. This function:
#' \itemize{
#'   \item Reads each modality's expression matrix (X) and layers
#'   \item Creates Seurat assays with appropriate names
#'   \item Preserves global cell metadata from \code{/obs}
#'   \item Reads per-modality embeddings (obsm) and graphs (obsp)
#'   \item Restores spatial data if present
#' }
#'
#' @section Modality Mapping:
#' By default, modality names are mapped to standard Seurat assay names:
#' \itemize{
#'   \item rna -> RNA
#'   \item prot -> ADT
#'   \item atac -> ATAC
#'   \item spatial -> Spatial
#'   \item Other names are preserved as-is
#' }
#'
#' @examples
#' \dontrun{
#' # Load all modalities from an h5mu file
#' seurat_obj <- LoadH5MU("multimodal_data.h5mu")
#'
#' # Load specific modalities only
#' seurat_obj <- LoadH5MU("multimodal_data.h5mu", modalities = c("rna", "prot"))
#'
#' # Custom assay name mapping
#' seurat_obj <- LoadH5MU(
#'   "data.h5mu",
#'   assay.names = c(rna = "RNA", prot = "Protein")
#' )
#' }
#'
#' @importFrom hdf5r H5File h5attr h5attr_names
#' @importFrom Seurat CreateSeuratObject CreateAssayObject Assays DefaultAssay<-
#'   CreateDimReducObject SetAssayData
#' @importFrom SeuratObject Cells AddMetaData
#' @importFrom Matrix sparseMatrix
#'
#' @export
#'
LoadH5MU <- function(file,
                     modalities = NULL,
                     assay.names = NULL,
                     restore.spatial = TRUE,
                     verbose = TRUE) {

  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }

  h5mu <- H5File$new(file, mode = "r")
  on.exit(tryCatch(h5mu$close_all(), error = function(e) NULL))

  if (verbose) message("Loading H5MU file: ", file)

  # --- Helper: read sparse or dense matrix from an HDF5 group/dataset ---
  ReadMatrix <- function(h5_obj, transpose = TRUE) {
    if (inherits(h5_obj, "H5Group")) {
      if (h5_obj$exists("data") && h5_obj$exists("indices") && h5_obj$exists("indptr")) {
        data_vals <- h5_obj[["data"]][]
        indices <- h5_obj[["indices"]][]
        indptr <- h5_obj[["indptr"]][]

        if (h5_obj$attr_exists("shape")) {
          shape <- h5attr(h5_obj, "shape")
          n_rows <- shape[1]
          n_cols <- shape[2]
        } else {
          n_rows <- length(indptr) - 1
          n_cols <- if (length(indices) > 0) max(indices) + 1L else 0L
        }

        if (transpose) {
          # Zero-copy CSR->CSC: CSR of (n_rows x n_cols) == CSC of (n_cols x n_rows)
          mat <- new("dgCMatrix",
            i = as.integer(indices),
            p = as.integer(indptr),
            x = as.numeric(data_vals),
            Dim = c(as.integer(n_cols), as.integer(n_rows))
          )
        } else {
          indices_1based <- indices + 1L
          row_indices <- rep(seq_len(n_rows), diff(indptr))
          mat <- sparseMatrix(
            i = row_indices, j = indices_1based, x = data_vals,
            dims = c(n_rows, n_cols), index1 = TRUE
          )
        }
        return(mat)
      }
    } else if (inherits(h5_obj, "H5D")) {
      # Dense matrix
      mat <- h5_obj[, ]
      if (transpose) mat <- t(mat)
      return(mat)
    }
    stop("Unknown matrix format", call. = FALSE)
  }

  # --- Helper: read index (cell/feature names) from a df group ---
  ReadIndex <- function(grp) {
    if (grp$exists("_index")) {
      return(as.character(grp[["_index"]][]))
    } else if (grp$exists("index")) {
      return(as.character(grp[["index"]][]))
    }
    NULL
  }

  # --- Helper: read dataframe columns from an obs/var group ---
  ReadDFColumns <- function(grp, exclude = c("_index", "index", "__categories")) {
    col_names <- setdiff(names(grp), exclude)
    if (length(col_names) == 0) return(NULL)

    # Pre-classify columns using ls()
    grp_ls <- grp$ls()
    grp_types <- setNames(as.character(grp_ls$obj_type), grp_ls$name)

    # Legacy categorical support
    has_legacy_cats <- grp$exists("__categories")
    legacy_cats <- if (has_legacy_cats) grp[["__categories"]] else NULL
    legacy_cat_names <- if (has_legacy_cats) names(legacy_cats) else character(0)

    result <- list()
    for (col in col_names) {
      if (!col %in% names(grp_types)) next
      tryCatch({
        if (grp_types[col] == "H5I_GROUP") {
          # Modern categorical
          col_obj <- grp[[col]]
          enc_type <- tryCatch(h5attr(col_obj, "encoding-type"), error = function(e) "")
          if (enc_type == "categorical" && col_obj$exists("codes") && col_obj$exists("categories")) {
            codes <- col_obj[["codes"]]$read()
            categories <- as.character(col_obj[["categories"]]$read())
            result[[col]] <- DecodeCategorical(codes, categories)
          }
        } else if (grp_types[col] == "H5I_DATASET") {
          if (has_legacy_cats && col %in% legacy_cat_names) {
            codes <- grp[[col]]$read()
            categories <- as.character(legacy_cats[[col]]$read())
            result[[col]] <- DecodeCategorical(codes, categories)
          } else {
            vals <- grp[[col]]$read()
            if (is.character(vals)) vals <- as.character(vals)
            result[[col]] <- vals
          }
        }
      }, error = function(e) {
        # skip unreadable columns
      })
    }
    if (length(result) == 0) return(NULL)
    result
  }

  # ========== Get available modalities ==========
  if (!h5mu$exists("mod")) {
    stop("No modalities found in h5mu file (missing /mod group)", call. = FALSE)
  }
  available_modalities <- names(h5mu[["mod"]])
  if (length(available_modalities) == 0) {
    stop("No modalities found in h5mu file", call. = FALSE)
  }

  if (verbose) {
    message("Found ", length(available_modalities), " modalities: ",
            paste(available_modalities, collapse = ", "))
  }

  # ========== Determine which modalities to load ==========
  if (is.null(modalities)) {
    modalities_to_load <- available_modalities
  } else {
    modalities_to_load <- modalities
    missing <- setdiff(modalities_to_load, available_modalities)
    if (length(missing) > 0) {
      warning("Requested modalities not found in file: ",
              paste(missing, collapse = ", "), immediate. = TRUE)
      modalities_to_load <- intersect(modalities_to_load, available_modalities)
    }
  }
  if (length(modalities_to_load) == 0) {
    stop("No valid modalities to load", call. = FALSE)
  }

  # ========== Set up assay name mapping ==========
  if (is.null(assay.names)) {
    assay.names <- GetDefaultModalityMapping(modalities_to_load)
  } else {
    default_mapping <- GetDefaultModalityMapping(modalities_to_load)
    for (mod in modalities_to_load) {
      if (!mod %in% names(assay.names)) {
        assay.names[mod] <- default_mapping[mod]
      }
    }
  }

  # ========== Read global cell names ==========
  cell_names <- NULL
  if (h5mu$exists("obs")) {
    cell_names <- ReadIndex(h5mu[["obs"]])
  }
  # Fallback: read from first modality
  if (is.null(cell_names)) {
    first_mod <- modalities_to_load[1]
    mod_path <- paste0("mod/", first_mod, "/obs")
    if (h5mu$exists(mod_path)) {
      cell_names <- ReadIndex(h5mu[[mod_path]])
    }
  }

  # ========== Load each modality ==========
  seurat_obj <- NULL
  for (mod_name in modalities_to_load) {
    assay_name <- assay.names[mod_name]
    mod_path <- paste0("mod/", mod_name)
    mod_grp <- h5mu[[mod_path]]

    if (verbose) message("  Loading modality '", mod_name, "' -> assay '", assay_name, "'")

    # --- Read feature names from var ---
    gene_names <- NULL
    if (mod_grp$exists("var")) {
      gene_names <- ReadIndex(mod_grp[["var"]])
    }

    # --- Read X (expression matrix) ---
    if (!mod_grp$exists("X")) {
      if (verbose) warning("Modality '", mod_name, "' has no X matrix, skipping",
                           immediate. = TRUE)
      next
    }

    mat <- ReadMatrix(mod_grp[["X"]], transpose = TRUE)

    # Read cell names for this modality
    mod_cell_names <- NULL
    if (mod_grp$exists("obs")) {
      mod_cell_names <- ReadIndex(mod_grp[["obs"]])
    }
    n_genes_expected <- if (!is.null(gene_names)) length(gene_names) else nrow(mat)
    n_cells_expected <- if (!is.null(mod_cell_names)) length(mod_cell_names)
                        else if (!is.null(cell_names)) length(cell_names)
                        else ncol(mat)

    # Smart dimension detection: transpose if matrix is in (cells x genes) orientation
    if (nrow(mat) == n_cells_expected && ncol(mat) == n_genes_expected &&
        nrow(mat) != n_genes_expected) {
      mat <- t(mat)
    }

    # Assign names
    if (!is.null(gene_names)) {
      if (length(gene_names) == nrow(mat)) {
        rownames(mat) <- gene_names
      } else {
        rownames(mat) <- gene_names[seq_len(nrow(mat))]
      }
    } else {
      rownames(mat) <- paste0(mod_name, "_Gene", seq_len(nrow(mat)))
    }

    if (!is.null(mod_cell_names) && length(mod_cell_names) == ncol(mat)) {
      colnames(mat) <- mod_cell_names
    } else if (!is.null(cell_names) && length(cell_names) == ncol(mat)) {
      colnames(mat) <- cell_names
    } else {
      colnames(mat) <- paste0("Cell", seq_len(ncol(mat)))
    }

    # --- Create or add assay ---
    if (is.null(seurat_obj)) {
      seurat_obj <- CreateSeuratObject(
        counts = mat,
        project = "H5MU",
        assay = assay_name,
        min.cells = 0,
        min.features = 0
      )
      # Update cell_names from the Seurat object
      cell_names <- Cells(seurat_obj)
    } else {
      # Add as additional assay
      new_assay <- CreateAssayObject(counts = mat)
      seurat_obj[[assay_name]] <- new_assay
    }

    # --- Read additional layers ---
    if (mod_grp$exists("layers")) {
      layer_names <- names(mod_grp[["layers"]])
      for (lname in layer_names) {
        tryCatch({
          lmat <- ReadMatrix(mod_grp[["layers"]][[lname]], transpose = TRUE)
          if (nrow(lmat) == nrow(mat) && ncol(lmat) == ncol(mat)) {
            rownames(lmat) <- rownames(mat)
            colnames(lmat) <- colnames(mat)
            seurat_slot <- switch(lname,
              "counts" = "counts",
              "data" = "data",
              "log_normalized" = "data",
              "scale.data" = "scale.data",
              "scaled" = "scale.data",
              lname
            )
            seurat_obj[[assay_name]] <- SetAssayData(
              object = seurat_obj[[assay_name]],
              layer = seurat_slot,
              new.data = lmat
            )
          }
        }, error = function(e) {
          if (verbose) warning("Could not read layer '", lname,
                               "' from modality '", mod_name, "'", immediate. = TRUE)
        })
      }
    }

    # --- Read embeddings (obsm) ---
    if (mod_grp$exists("obsm")) {
      obsm_grp <- mod_grp[["obsm"]]
      for (reduc_name in names(obsm_grp)) {
        # Skip empty dict groups
        if (!inherits(obsm_grp[[reduc_name]], "H5D") &&
            !(inherits(obsm_grp[[reduc_name]], "H5Group") &&
              obsm_grp[[reduc_name]]$exists("data"))) {
          # Try reading as dense dataset
          emb_ok <- tryCatch({
            emb <- obsm_grp[[reduc_name]][, ]
            !is.null(emb) && length(emb) > 0
          }, error = function(e) FALSE)
          if (!emb_ok) next
        }

        tryCatch({
          emb <- obsm_grp[[reduc_name]][, ]
          if (nrow(emb) != length(cell_names) && ncol(emb) == length(cell_names)) {
            emb <- t(emb)
          }
          if (nrow(emb) != length(cell_names)) {
            next
          }

          rownames(emb) <- cell_names
          clean_name <- gsub("^X_", "", reduc_name)
          key <- paste0(clean_name, "_")

          reduc_obj <- CreateDimReducObject(
            embeddings = emb,
            key = key,
            assay = assay_name
          )
          seurat_obj[[clean_name]] <- reduc_obj
          if (verbose) message("    Added reduction: ", clean_name)
        }, error = function(e) {
          # skip unreadable embeddings
        })
      }
    }

    # --- Read graphs (obsp) ---
    if (mod_grp$exists("obsp")) {
      obsp_grp <- mod_grp[["obsp"]]
      for (graph_name in names(obsp_grp)) {
        tryCatch({
          graph_mat <- ReadMatrix(obsp_grp[[graph_name]], transpose = FALSE)
          if (nrow(graph_mat) == length(cell_names) && ncol(graph_mat) == length(cell_names)) {
            rownames(graph_mat) <- cell_names
            colnames(graph_mat) <- cell_names
            full_name <- paste0(assay_name, "_", graph_name)
            graph_obj <- as.Graph(graph_mat)
            DefaultAssay(graph_obj) <- assay_name
            seurat_obj[[full_name]] <- graph_obj
            if (verbose) message("    Added graph: ", full_name)
          }
        }, error = function(e) {
          # skip unreadable graphs
        })
      }
    }

    # --- Read varp (pairwise variable annotations) ---
    if (mod_grp$exists("varp")) {
      varp_list <- list()
      for (varp_name in names(mod_grp[["varp"]])) {
        tryCatch({
          varp_mat <- ReadMatrix(mod_grp[["varp"]][[varp_name]], transpose = FALSE)
          feat_names <- rownames(mat)
          if (nrow(varp_mat) == length(feat_names) && ncol(varp_mat) == length(feat_names)) {
            rownames(varp_mat) <- feat_names
            colnames(varp_mat) <- feat_names
          }
          varp_list[[varp_name]] <- varp_mat
          if (verbose) message("    Added varp: ", varp_name)
        }, error = function(e) {
          # skip unreadable varp items
        })
      }
      if (length(varp_list) > 0) {
        seurat_obj@misc[[paste0("__varp__.", assay_name)]] <- varp_list
      }
    }

    # --- Read var metadata (feature metadata) ---
    if (mod_grp$exists("var")) {
      var_cols <- ReadDFColumns(mod_grp[["var"]])
      if (!is.null(var_cols) && length(var_cols) > 0) {
        var_df <- as.data.frame(var_cols, stringsAsFactors = FALSE)
        if (nrow(var_df) == nrow(mat)) {
          rownames(var_df) <- rownames(mat)
          tryCatch({
            seurat_obj[[assay_name]][[names(var_cols)]] <- var_df
          }, error = function(e) {
            # skip if can't set feature metadata
          })
        }
      }
    }
  }

  if (is.null(seurat_obj)) {
    stop("No modalities could be loaded", call. = FALSE)
  }

  # ========== Read global obs metadata ==========
  if (h5mu$exists("obs")) {
    if (verbose) message("  Reading global metadata...")
    global_cols <- ReadDFColumns(h5mu[["obs"]])
    if (!is.null(global_cols)) {
      for (col_name in names(global_cols)) {
        vals <- global_cols[[col_name]]
        if (length(vals) == ncol(seurat_obj)) {
          seurat_obj[[col_name]] <- vals
        }
      }
      if (verbose) message("  Added ", length(global_cols), " global metadata columns")
    }
  }

  # ========== Restore spatial data if requested ==========
  if (restore.spatial) {
    seurat_obj <- RestoreSpatialFromH5MU(
      h5mu_file = h5mu,
      seurat_obj = seurat_obj,
      modalities = modalities_to_load,
      assay.names = assay.names,
      verbose = verbose
    )
  }

  # ========== Set default assay ==========
  final_assays <- Assays(seurat_obj)
  if ("RNA" %in% final_assays) {
    DefaultAssay(seurat_obj) <- "RNA"
  } else if (length(final_assays) > 0) {
    DefaultAssay(seurat_obj) <- final_assays[1]
  }

  if (verbose) {
    message("\nSuccessfully loaded H5MU file")
    message("  Cells: ", ncol(seurat_obj))
    message("  Assays: ", paste(Assays(seurat_obj), collapse = ", "))
    if (length(seurat_obj@reductions) > 0) {
      message("  Reductions: ", paste(names(seurat_obj@reductions), collapse = ", "))
    }
    if (length(seurat_obj@graphs) > 0) {
      message("  Graphs: ", paste(names(seurat_obj@graphs), collapse = ", "))
    }
    message("  Metadata columns: ", ncol(seurat_obj@meta.data))
  }

  return(seurat_obj)
}


#' Get default modality to assay name mapping
#'
#' @param modalities Character vector of modality names
#' @return Named character vector mapping modality names to assay names
#' @keywords internal
#'
GetDefaultModalityMapping <- function(modalities) {
  # Standard mapping based on muon/MuData conventions
  default_map <- c(
    rna = "RNA",
    sct = "SCT",      # SCTransform assay
    prot = "ADT",
    atac = "ATAC",
    spatial = "Spatial",
    hto = "HTO",
    adt = "ADT",
    cite = "ADT",
    peaks = "ATAC"
  )

  mapping <- setNames(modalities, modalities)
  for (mod in modalities) {
    if (tolower(mod) %in% names(default_map)) {
      mapping[mod] <- default_map[tolower(mod)]
    }
  }

  return(mapping)
}


#' Restore spatial data from h5mu file to Seurat object
#'
#' @param h5mu_file H5File connection to h5mu file
#' @param seurat_obj Seurat object to add spatial data to
#' @param modalities Character vector of modalities that were loaded
#' @param assay.names Named vector mapping modality to assay names
#' @param verbose Logical; print messages
#' @return Modified Seurat object with spatial data
#' @keywords internal
#'
RestoreSpatialFromH5MU <- function(h5mu_file, seurat_obj, modalities,
                                   assay.names, verbose = FALSE) {

  # Check each modality for spatial data
  for (mod in modalities) {
    mod_path <- paste0("mod/", mod)
    if (!h5mu_file$exists(mod_path)) next

    mod_group <- h5mu_file[[mod_path]]

    # Check for spatial coordinates in obsm
    if (mod_group$exists("obsm") && "spatial" %in% names(mod_group[["obsm"]])) {
      assay_name <- assay.names[mod]

      if (verbose) message("  Found spatial data in modality '", mod, "'")

      seurat_obj <- ConvertH5MUSpatialToSeurat(
        h5mu_file = h5mu_file,
        seurat_obj = seurat_obj,
        modality = mod,
        assay_name = assay_name,
        verbose = verbose
      )
    }
  }

  return(seurat_obj)
}


#' Convert h5mu spatial data to Seurat format for a specific modality
#'
#' @param h5mu_file H5File connection to h5mu file
#' @param seurat_obj Seurat object
#' @param modality Modality name in h5mu file
#' @param assay_name Corresponding assay name in Seurat
#' @param verbose Logical
#' @return Modified Seurat object
#' @keywords internal
#'
#' @importFrom hdf5r H5File h5attr
#'
ConvertH5MUSpatialToSeurat <- function(h5mu_file, seurat_obj, modality,
                                       assay_name, verbose = FALSE) {

  mod_path <- paste0("mod/", modality)
  mod_group <- h5mu_file[[mod_path]]

  if (!mod_group$exists("obsm") || !"spatial" %in% names(mod_group[["obsm"]])) {
    return(seurat_obj)
  }

  spatial_coords <- mod_group[["obsm"]][["spatial"]][, ]
  spatial_coords <- as.matrix(spatial_coords)

  cell_names <- Cells(seurat_obj)

  if (nrow(spatial_coords) != length(cell_names)) {
    if (verbose) {
      warning("Spatial coordinate count doesn't match cell count in modality '",
              modality, "'", immediate. = TRUE)
    }
    return(seurat_obj)
  }

  rownames(spatial_coords) <- cell_names

  # Check for Visium-style data in uns
  visium_data <- NULL
  if (mod_group$exists("uns") && "spatial" %in% names(mod_group[["uns"]])) {
    spatial_uns <- mod_group[["uns/spatial"]]
    library_ids <- names(spatial_uns)

    if (length(library_ids) > 0) {
      visium_data <- list()
      for (lib_id in library_ids) {
        lib_group <- spatial_uns[[lib_id]]
        lib_data <- list()

        if ("scalefactors" %in% names(lib_group)) {
          sf_group <- lib_group[["scalefactors"]]
          lib_data$scalefactors <- list()
          for (sf_name in names(sf_group)) {
            lib_data$scalefactors[[sf_name]] <- sf_group[[sf_name]][]
          }
        }

        if ("images" %in% names(lib_group)) {
          lib_data$has_images <- TRUE
          lib_data$image_names <- names(lib_group[["images"]])
        }

        visium_data[[lib_id]] <- lib_data
      }
    }
  }

  # Add coordinates to metadata
  seurat_obj@meta.data$spatial_x <- spatial_coords[, 1]
  seurat_obj@meta.data$spatial_y <- spatial_coords[, 2]

  if (!is.null(visium_data)) {
    seurat_obj@misc[[paste0(modality, "_spatial_metadata")]] <- visium_data
    seurat_obj@misc[[paste0(modality, "_spatial_technology")]] <- "Visium"

    if (verbose) {
      message("    Found Visium spatial metadata for ",
              length(visium_data), " library(ies)")
    }
  }

  return(seurat_obj)
}
