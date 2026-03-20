#' @include zzz.R
#' @include Convert.R
#' @include AnnDataEncoding.R
#' @importFrom hdf5r H5File h5attr H5S GuessDType
#' @importFrom Matrix sparseMatrix t
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Modality / Assay Name Mapping Helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Map h5mu modality name to Seurat assay name
#'
#' @param mod Character modality name from h5mu (e.g., "rna", "prot")
#'
#' @return Character assay name (e.g., "RNA", "ADT")
#'
#' @keywords internal
#'
.modality_to_assay <- function(mod) {
  map <- c(rna = "RNA", prot = "ADT", atac = "ATAC", spatial = "Spatial",
           hto = "HTO", sct = "SCT")
  match <- map[tolower(mod)]
  if (!is.na(match)) unname(match) else mod
}

#' Map Seurat assay name to h5mu modality name
#'
#' @param assay Character assay name (e.g., "RNA", "ADT")
#'
#' @return Character modality name (e.g., "rna", "prot")
#'
#' @keywords internal
#'
.assay_to_modality <- function(assay) {
  map <- c(RNA = "rna", ADT = "prot", ATAC = "atac", Spatial = "spatial",
           HTO = "hto", SCT = "sct")
  match <- map[assay]
  if (!is.na(match)) unname(match) else tolower(assay)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Shared streaming helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Stream a sparse matrix group from one HDF5 file to another
#'
#' Copies data/indices/indptr arrays with optional CSR<->CSC reinterpretation.
#' CSR of [A x B] has identical on-disk layout to CSC of [B x A], so conversion
#' between h5ad CSR (cells x genes) and h5seurat CSC (genes x cells) requires
#' only swapping the shape/dims attribute.
#'
#' @param src_group H5Group containing data/indices/indptr
#' @param dst_parent H5Group to write into
#' @param dst_name Name for new group in dst_parent
#' @param gzip Integer gzip compression level
#' @param src_format Character, "h5ad" or "h5seurat" or "h5mu"
#' @param dst_format Character, "h5ad" or "h5seurat"
#' @param assay Character assay name (for h5seurat graphs)
#'
#' @keywords internal
#'
.stream_sparse_group <- function(src_group, dst_parent, dst_name, gzip,
                                 src_format = "h5ad", dst_format = "h5seurat",
                                 assay = NULL) {
  if (!src_group$exists("data") || !src_group$exists("indices") ||
      !src_group$exists("indptr")) {
    return(invisible(NULL))
  }

  data_vals <- src_group[["data"]]$read()
  indices   <- src_group[["indices"]]$read()
  indptr    <- src_group[["indptr"]]$read()

  # Read shape from source
  shape <- NULL
  if (src_group$attr_exists("shape")) {
    shape <- hdf5r::h5attr(src_group, "shape")
  } else if (src_group$attr_exists("dims")) {
    shape <- hdf5r::h5attr(src_group, "dims")
  }

  chunk_size <- 65536L
  grp <- dst_parent$create_group(dst_name)
  grp$create_dataset("data", robj = as.numeric(data_vals),
                     chunk_dims = min(length(data_vals), chunk_size),
                     gzip_level = gzip)
  grp$create_dataset("indices", robj = as.integer(indices),
                     chunk_dims = min(length(indices), chunk_size),
                     gzip_level = gzip)
  grp$create_dataset("indptr", robj = as.integer(indptr),
                     chunk_dims = min(length(indptr), chunk_size),
                     gzip_level = gzip)

  # Write shape/dims attribute
  # CSR [cells x genes] <-> CSC [genes x cells]: zero-copy, just swap shape
  if (dst_format == "h5seurat") {
    dims_val <- if (!is.null(shape)) as.integer(rev(shape)) else NULL
    if (!is.null(dims_val)) {
      grp$create_attr(attr_name = "dims", robj = dims_val,
                      dtype = GuessDType(dims_val))
    }
    if (!is.null(assay)) {
      grp$create_attr(attr_name = "assay.used", robj = assay,
                      dtype = CachedGuessDType(assay), space = ScalarSpace())
    }
  } else {
    # h5ad / h5mu target
    shape_val <- if (!is.null(shape)) as.integer(rev(shape)) else NULL
    if (!is.null(shape_val)) {
      grp$create_attr(attr_name = "shape", robj = shape_val,
                      dtype = GuessDType(shape_val))
    }
    grp$create_attr(attr_name = "encoding-type", robj = "csr_matrix",
                    dtype = CachedGuessDType("csr_matrix"), space = ScalarSpace())
    grp$create_attr(attr_name = "encoding-version", robj = "0.1.0",
                    dtype = CachedGuessDType("0.1.0"), space = ScalarSpace())
  }

  invisible(grp)
}

#' Stream a dense matrix from one HDF5 file to another
#'
#' Reads a dense dataset and writes it with gzip compression.
#'
#' @param src_dataset H5D source dataset
#' @param dst_parent H5Group to write into
#' @param dst_name Name for new dataset
#' @param gzip Integer gzip compression level
#' @param transpose Logical, whether to transpose
#'
#' @keywords internal
#'
.stream_dense_matrix <- function(src_dataset, dst_parent, dst_name, gzip,
                                 transpose = FALSE) {
  mat <- src_dataset$read()
  if (transpose && is.matrix(mat)) mat <- t(mat)
  dims <- if (is.matrix(mat)) dim(mat) else length(mat)
  dst_parent$create_dataset(dst_name, robj = mat,
                            chunk_dims = dims, gzip_level = gzip)
  invisible(NULL)
}

#' Stream an h5ad-style DataFrame group to h5seurat meta.data
#'
#' Reads obs columns (categorical or plain) from an AnnData DataFrame group
#' and writes them to h5seurat meta.data format (factors use levels/values).
#'
#' @param src_group H5Group for the source DataFrame (obs)
#' @param dst_parent H5Group or H5File to write meta.data into
#' @param gzip Integer compression level
#' @param verbose Logical
#'
#' @return Character vector of column names written
#'
#' @keywords internal
#'
.stream_h5ad_df_to_h5seurat_md <- function(src_group, dst_parent, gzip,
                                           verbose = TRUE) {
  md_grp <- dst_parent$create_group("meta.data")
  col_order <- character(0)

  if (src_group$attr_exists("column-order")) {
    col_order <- hdf5r::h5attr(src_group, "column-order")
  } else {
    col_order <- setdiff(names(src_group), c("_index", "index", "__categories"))
  }

  # Legacy categorical support
  has_cats <- src_group$exists("__categories")
  cats_grp <- if (has_cats) src_group[["__categories"]] else NULL
  cat_names <- if (has_cats) names(cats_grp) else character(0)

  col_written <- character(0)
  for (col in col_order) {
    if (!src_group$exists(col)) next
    tryCatch({
      col_obj <- src_group[[col]]

      if (inherits(col_obj, "H5Group")) {
        # Modern categorical -> h5seurat factor
        enc_type <- tryCatch(hdf5r::h5attr(col_obj, "encoding-type"),
                             error = function(e) "")
        if (enc_type == "categorical" &&
            col_obj$exists("categories") && col_obj$exists("codes")) {
          cats <- as.character(col_obj[["categories"]]$read())
          codes <- col_obj[["codes"]]$read()
          codes_int <- as.integer(codes)
          values_1based <- codes_int + 1L
          values_1based[codes_int == -1L] <- NA_integer_

          fac_grp <- md_grp$create_group(col)
          fac_grp$create_dataset("levels", robj = cats,
                                 dtype = CachedUtf8Type(),
                                 chunk_dims = length(cats), gzip_level = gzip)
          fac_grp$create_dataset("values", robj = values_1based,
                                 chunk_dims = length(values_1based),
                                 gzip_level = gzip)
          col_written <- c(col_written, col)
        }
      } else if (inherits(col_obj, "H5D")) {
        if (has_cats && col %in% cat_names) {
          # Legacy categorical
          codes <- col_obj$read()
          cats <- as.character(cats_grp[[col]]$read())
          codes_int <- as.integer(codes)
          values_1based <- codes_int + 1L
          values_1based[codes_int == -1L] <- NA_integer_

          fac_grp <- md_grp$create_group(col)
          fac_grp$create_dataset("levels", robj = cats,
                                 dtype = CachedUtf8Type(),
                                 chunk_dims = length(cats), gzip_level = gzip)
          fac_grp$create_dataset("values", robj = values_1based,
                                 chunk_dims = length(values_1based),
                                 gzip_level = gzip)
          col_written <- c(col_written, col)
        } else {
          vals <- col_obj$read()
          if (is.character(vals)) {
            md_grp$create_dataset(col, robj = vals, dtype = CachedUtf8Type(),
                                  chunk_dims = length(vals), gzip_level = gzip)
          } else if (is.logical(vals)) {
            md_grp$create_dataset(col, robj = as.integer(vals),
                                  chunk_dims = length(vals), gzip_level = gzip)
          } else {
            md_grp$create_dataset(col, robj = vals,
                                  chunk_dims = length(vals), gzip_level = gzip)
          }
          col_written <- c(col_written, col)
        }
      }
    }, error = function(e) {
      if (verbose) warning("Could not stream obs column ", col, ": ",
                           e$message, immediate. = TRUE)
    })
  }

  if (length(col_written) > 0) {
    md_grp$create_attr(attr_name = "colnames", robj = col_written,
                       dtype = GuessDType(col_written))
  }

  invisible(col_written)
}

#' Stream h5seurat meta.data to an h5ad-style DataFrame group
#'
#' Reads meta.data columns (factors or plain datasets) and writes them
#' as AnnData DataFrame format (categoricals use categories/codes).
#'
#' @param md_group H5Group for meta.data
#' @param dst_parent H5Group or H5File to write into
#' @param dst_name Name for the DataFrame group (e.g., "obs")
#' @param cell_names Character vector of cell barcodes for _index
#' @param gzip Integer compression level
#' @param verbose Logical
#'
#' @return Invisible NULL
#'
#' @keywords internal
#'
.stream_h5seurat_md_to_h5ad_df <- function(md_group, dst_parent, dst_name,
                                           cell_names, gzip, verbose = TRUE) {
  md_cols <- character(0)
  if (md_group$attr_exists("colnames")) {
    md_cols <- hdf5r::h5attr(md_group, "colnames")
  } else {
    md_cols <- names(md_group)
  }

  grp <- dst_parent$create_group(dst_name)
  grp$create_attr(attr_name = "encoding-type", robj = "dataframe",
                  dtype = CachedGuessDType("dataframe"), space = ScalarSpace())
  grp$create_attr(attr_name = "encoding-version", robj = "0.2.0",
                  dtype = CachedGuessDType("0.2.0"), space = ScalarSpace())
  grp$create_attr(attr_name = "_index", robj = "_index",
                  dtype = CachedGuessDType("_index"), space = ScalarSpace())

  # Write _index
  grp$create_dataset("_index", robj = cell_names,
                     dtype = CachedUtf8Type(),
                     chunk_dims = length(cell_names), gzip_level = gzip)
  AddAnndataEncoding(grp[["_index"]], encoding_type = "string-array")

  col_written <- character(0)
  for (col in md_cols) {
    if (!md_group$exists(col)) next
    tryCatch({
      col_obj <- md_group[[col]]

      if (inherits(col_obj, "H5Group")) {
        # h5seurat factor -> h5ad categorical
        if (col_obj$exists("levels") && col_obj$exists("values")) {
          levs <- as.character(col_obj[["levels"]]$read())
          vals <- col_obj[["values"]]$read()
          codes <- as.integer(vals) - 1L  # 1-based -> 0-based

          cat_grp <- grp$create_group(col)
          cat_grp$create_attr(attr_name = "encoding-type", robj = "categorical",
                              dtype = CachedGuessDType("categorical"),
                              space = ScalarSpace())
          cat_grp$create_attr(attr_name = "encoding-version", robj = "0.2.0",
                              dtype = CachedGuessDType("0.2.0"),
                              space = ScalarSpace())
          cat_grp$create_attr(attr_name = "ordered", robj = FALSE,
                              dtype = GuessDType(TRUE), space = ScalarSpace())
          cat_grp$create_dataset("codes", robj = codes,
                                 chunk_dims = length(codes), gzip_level = gzip)
          cat_grp$create_dataset("categories", robj = levs,
                                 dtype = CachedUtf8Type(),
                                 chunk_dims = length(levs), gzip_level = gzip)
          AddAnndataEncoding(cat_grp[["categories"]],
                             encoding_type = "string-array")
          col_written <- c(col_written, col)
        }
      } else if (inherits(col_obj, "H5D")) {
        vals <- col_obj$read()
        if (is.character(vals)) {
          grp$create_dataset(col, robj = vals, dtype = CachedUtf8Type(),
                             chunk_dims = length(vals), gzip_level = gzip)
          AddAnndataEncoding(grp[[col]], encoding_type = "string-array")
        } else if (is.logical(vals)) {
          grp$create_dataset(col, robj = BoolToInt(vals),
                             chunk_dims = length(vals), gzip_level = gzip)
        } else {
          grp$create_dataset(col, robj = vals,
                             chunk_dims = length(vals), gzip_level = gzip)
        }
        col_written <- c(col_written, col)
      }
    }, error = function(e) {
      if (verbose) warning("Could not stream meta.data column ", col, ": ",
                           e$message, immediate. = TRUE)
    })
  }

  grp$create_attr(attr_name = "column-order", robj = col_written,
                  dtype = CachedUtf8Type())
  invisible(NULL)
}

#' Write required h5seurat root attributes and empty groups
#'
#' @param h5 H5File handle (mode = "w")
#' @param assay Character active assay name
#' @param n_cells Integer number of cells (for active.ident)
#'
#' @keywords internal
#'
.h5seurat_write_scaffold <- function(h5, assay, n_cells) {
  # Root attributes
  h5$create_attr(attr_name = "active.assay", robj = assay,
                 dtype = CachedGuessDType(assay), space = ScalarSpace())
  h5$create_attr(attr_name = "project", robj = "scConvert",
                 dtype = CachedGuessDType("scConvert"), space = ScalarSpace())
  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("SeuratObject")),
    error = function(e) "5.0.0"
  )
  h5$create_attr(attr_name = "version", robj = pkg_ver,
                 dtype = CachedGuessDType(pkg_ver), space = ScalarSpace())

  # Required empty groups
  for (g in c("tools", "commands", "images", "neighbors")) {
    if (!h5$exists(g)) h5$create_group(g)
  }
  if (!h5$exists("misc")) h5$create_group("misc")

  # active.ident (default: all cells = single identity)
  if (n_cells > 0) {
    ident_grp <- h5$create_group("active.ident")
    ident_grp$create_dataset("levels", robj = "scConvert",
                             dtype = CachedUtf8Type(),
                             chunk_dims = 1L, gzip_level = 4L)
    ident_grp$create_dataset("values", robj = rep(1L, n_cells),
                             chunk_dims = n_cells, gzip_level = 4L)
  }

  invisible(NULL)
}

#' Write an AnnData dict-style empty group
#'
#' @param parent H5Group or H5File
#' @param name Group name
#'
#' @keywords internal
#'
.h5ad_write_dict_group <- function(parent, name) {
  if (parent$exists(name)) return(invisible(parent[[name]]))
  grp <- parent$create_group(name)
  grp$create_attr(attr_name = "encoding-type", robj = "dict",
                  dtype = CachedGuessDType("dict"), space = ScalarSpace())
  grp$create_attr(attr_name = "encoding-version", robj = "0.1.0",
                  dtype = CachedGuessDType("0.1.0"), space = ScalarSpace())
  invisible(grp)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Core streaming converter: h5mu -> h5seurat
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Stream h5mu to h5seurat (direct HDF5-to-HDF5)
#'
#' For each modality in /mod/, converts X, obs, var, obsm, obsp to h5seurat
#' layout. Uses the zero-copy CSR<->CSC reinterpretation trick.
#'
#' @param h5mu_path Path to input h5mu file
#' @param h5seurat_path Path to output h5seurat file
#' @param gzip Integer gzip compression level
#' @param verbose Logical
#'
#' @keywords internal
#'
.stream_h5mu_to_h5seurat <- function(h5mu_path, h5seurat_path,
                                     gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for h5mu -> h5seurat streaming", call. = FALSE)

  src <- hdf5r::H5File$new(h5mu_path, mode = "r")
  on.exit(src$close_all(), add = TRUE)
  dst <- hdf5r::H5File$new(h5seurat_path, mode = "w")
  on.exit(dst$close_all(), add = TRUE)

  # Discover modalities
  if (!src$exists("mod"))
    stop("No /mod/ group found in h5mu file", call. = FALSE)
  mod_root <- src[["mod"]]
  mod_names <- if (mod_root$attr_exists("mod-order")) {
    hdf5r::h5attr(mod_root, "mod-order")
  } else {
    names(mod_root)
  }

  if (length(mod_names) == 0)
    stop("No modalities found in h5mu file", call. = FALSE)

  # First modality determines active assay
  first_assay <- .modality_to_assay(mod_names[1])

  # Read global cell names from /obs
  cell_names <- NULL
  if (src$exists("obs")) {
    obs_grp <- src[["obs"]]
    if (obs_grp$exists("_index")) {
      cell_names <- as.character(obs_grp[["_index"]]$read())
    } else if (obs_grp$exists("index")) {
      cell_names <- as.character(obs_grp[["index"]]$read())
    }
  }
  # Fall back to first modality's obs
  if (is.null(cell_names) && src$exists(paste0("mod/", mod_names[1], "/obs"))) {
    mod_obs <- src[[paste0("mod/", mod_names[1], "/obs")]]
    if (mod_obs$exists("_index")) {
      cell_names <- as.character(mod_obs[["_index"]]$read())
    }
  }
  if (is.null(cell_names)) cell_names <- character(0)
  n_cells <- length(cell_names)

  # Write cell.names
  if (n_cells > 0) {
    dst$create_dataset("cell.names", robj = cell_names,
                       dtype = CachedUtf8Type(),
                       chunk_dims = n_cells, gzip_level = gzip)
  }

  # Create assays root
  assays_grp <- dst$create_group("assays")
  reductions_grp <- dst$create_group("reductions")
  graphs_grp <- dst$create_group("graphs")

  for (mod_name in mod_names) {
    assay <- .modality_to_assay(mod_name)
    mod_path <- paste0("mod/", mod_name)
    if (!src$exists(mod_path)) next

    mod_grp <- src[[mod_path]]
    if (verbose) message("  Streaming modality '", mod_name, "' -> assay '",
                         assay, "'...")

    # Create assay group
    a_grp <- assays_grp$create_group(assay)
    a_grp$create_attr(attr_name = "key", robj = paste0(tolower(assay), "_"),
                      dtype = CachedGuessDType(paste0(tolower(assay), "_")),
                      space = ScalarSpace())
    a_grp$create_attr(attr_name = "s4class", robj = "SeuratObject::Assay5",
                      dtype = CachedGuessDType("SeuratObject::Assay5"),
                      space = ScalarSpace())

    # --- features (var/_index) ---
    feature_names <- character(0)
    if (mod_grp$exists("var")) {
      var_grp <- mod_grp[["var"]]
      if (var_grp$exists("_index")) {
        feature_names <- as.character(var_grp[["_index"]]$read())
      } else if (var_grp$exists("index")) {
        feature_names <- as.character(var_grp[["index"]]$read())
      }
    }
    if (length(feature_names) > 0) {
      a_grp$create_dataset("features", robj = feature_names,
                           dtype = CachedUtf8Type(),
                           chunk_dims = length(feature_names), gzip_level = gzip)
    }

    # --- X -> layers/data ---
    layers_grp <- a_grp$create_group("layers")
    if (mod_grp$exists("X")) {
      x_obj <- mod_grp[["X"]]
      if (inherits(x_obj, "H5Group")) {
        # Sparse: CSR [cells x genes] -> CSC [genes x cells] (zero-copy)
        .stream_sparse_group(x_obj, layers_grp, "data", gzip,
                             src_format = "h5ad", dst_format = "h5seurat")
      } else if (inherits(x_obj, "H5D")) {
        # Dense: read, transpose, write
        .stream_dense_matrix(x_obj, layers_grp, "data", gzip, transpose = TRUE)
      }
    }

    # --- layers/counts -> layers/counts ---
    raw_x_path <- NULL
    if (mod_grp$exists("raw") && mod_grp[["raw"]]$exists("X")) {
      raw_x_path <- "raw/X"
    } else if (mod_grp$exists("layers") && mod_grp[["layers"]]$exists("counts")) {
      raw_x_path <- "layers/counts"
    }
    if (!is.null(raw_x_path)) {
      raw_obj <- mod_grp[[raw_x_path]]
      if (inherits(raw_obj, "H5Group")) {
        .stream_sparse_group(raw_obj, layers_grp, "counts", gzip,
                             src_format = "h5ad", dst_format = "h5seurat")
      } else if (inherits(raw_obj, "H5D")) {
        .stream_dense_matrix(raw_obj, layers_grp, "counts", gzip,
                             transpose = TRUE)
      }
    }

    # --- obsm -> reductions (strip X_ prefix) ---
    if (mod_grp$exists("obsm")) {
      obsm_grp <- mod_grp[["obsm"]]
      for (rn in names(obsm_grp)) {
        clean_name <- gsub("^X_", "", rn)
        tryCatch({
          r_grp <- reductions_grp$create_group(clean_name)
          r_grp$create_attr(attr_name = "active.assay", robj = assay,
                            dtype = CachedGuessDType(assay),
                            space = ScalarSpace())
          key <- AnnDataReductionKey(clean_name)
          r_grp$create_attr(attr_name = "key", robj = key,
                            dtype = CachedGuessDType(key), space = ScalarSpace())
          r_grp$create_attr(attr_name = "global", robj = 0L,
                            dtype = GuessDType(0L), space = ScalarSpace())
          r_grp$create_group("misc")

          emb_obj <- obsm_grp[[rn]]
          if (inherits(emb_obj, "H5D")) {
            emb <- emb_obj$read()
            if (is.matrix(emb) && nrow(emb) != n_cells && ncol(emb) == n_cells) {
              emb <- t(emb)
            }
            r_grp$create_dataset("cell.embeddings", robj = emb,
                                 chunk_dims = dim(emb), gzip_level = gzip)
          }
        }, error = function(e) {
          if (verbose) warning("Could not stream reduction ", clean_name,
                               ": ", e$message, immediate. = TRUE)
        })
      }
    }

    # --- obsp -> graphs (prefix with assay_) ---
    if (mod_grp$exists("obsp")) {
      obsp_grp <- mod_grp[["obsp"]]
      for (gn in names(obsp_grp)) {
        graph_name <- paste0(assay, "_", gn)
        tryCatch({
          g_obj <- obsp_grp[[gn]]
          if (inherits(g_obj, "H5Group")) {
            .stream_sparse_group(g_obj, graphs_grp, graph_name, gzip,
                                 src_format = "h5ad", dst_format = "h5seurat",
                                 assay = assay)
          }
        }, error = function(e) {
          if (verbose) warning("Could not stream graph ", graph_name,
                               ": ", e$message, immediate. = TRUE)
        })
      }
    }
  }

  # --- Global obs -> meta.data ---
  if (verbose) message("  Streaming global obs -> meta.data...")
  if (src$exists("obs")) {
    .stream_h5ad_df_to_h5seurat_md(src[["obs"]], dst, gzip, verbose = verbose)
  } else {
    # Create minimal meta.data
    md_grp <- dst$create_group("meta.data")
    md_grp$create_attr(attr_name = "colnames", robj = character(0),
                       dtype = CachedUtf8Type())
  }

  # Root scaffold (attrs, empty groups, active.ident)
  .h5seurat_write_scaffold(dst, first_assay, n_cells)

  dst$flush()
  if (verbose) message("Streaming h5mu -> h5seurat complete: ", h5seurat_path)
  invisible(h5seurat_path)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Core streaming converter: h5seurat -> h5mu
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Stream h5seurat to h5mu (direct HDF5-to-HDF5)
#'
#' For each assay, creates a modality under /mod/ with the full AnnData layout.
#' Uses the zero-copy CSC<->CSR reinterpretation trick.
#'
#' @param h5seurat_path Path to input h5seurat file
#' @param h5mu_path Path to output h5mu file
#' @param gzip Integer gzip compression level
#' @param verbose Logical
#'
#' @keywords internal
#'
.stream_h5seurat_to_h5mu <- function(h5seurat_path, h5mu_path,
                                     gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for h5seurat -> h5mu streaming", call. = FALSE)

  src <- hdf5r::H5File$new(h5seurat_path, mode = "r")
  on.exit(src$close_all(), add = TRUE)
  dst <- hdf5r::H5File$new(h5mu_path, mode = "w")
  on.exit(dst$close_all(), add = TRUE)

  # MuData root attrs
  dst$create_attr(attr_name = "encoding-type", robj = "MuData",
                  dtype = CachedGuessDType("MuData"), space = ScalarSpace())
  dst$create_attr(attr_name = "encoding-version", robj = "0.1.0",
                  dtype = CachedGuessDType("0.1.0"), space = ScalarSpace())

  # Read cell names
  cell_names <- if (src$exists("cell.names")) {
    as.character(src[["cell.names"]]$read())
  } else {
    character(0)
  }
  n_cells <- length(cell_names)

  # Discover assays
  assay_names <- character(0)
  if (src$exists("assays")) {
    assay_names <- names(src[["assays"]])
  }
  if (length(assay_names) == 0)
    stop("No assays found in h5seurat file", call. = FALSE)

  # Modality mapping
  mod_order <- vapply(assay_names, .assay_to_modality, character(1),
                      USE.NAMES = FALSE)

  # Create /mod/ group
  mod_root <- dst$create_group("mod")
  mod_root$create_attr(attr_name = "mod-order", robj = mod_order,
                       dtype = CachedUtf8Type())

  # Read reduction info for assay assignment
  reduc_assay_map <- list()
  if (src$exists("reductions")) {
    red_root <- src[["reductions"]]
    for (rn in names(red_root)) {
      tryCatch({
        r_grp <- red_root[[rn]]
        r_assay <- if (r_grp$attr_exists("active.assay")) {
          hdf5r::h5attr(r_grp, "active.assay")
        } else {
          assay_names[1]
        }
        reduc_assay_map[[rn]] <- r_assay
      }, error = function(e) NULL)
    }
  }

  # Read graph info for assay assignment
  graph_assay_map <- list()
  if (src$exists("graphs")) {
    g_root <- src[["graphs"]]
    for (gn in names(g_root)) {
      tryCatch({
        gg <- g_root[[gn]]
        g_assay <- if (gg$attr_exists("assay.used")) {
          hdf5r::h5attr(gg, "assay.used")
        } else {
          assay_names[1]
        }
        graph_assay_map[[gn]] <- g_assay
      }, error = function(e) NULL)
    }
  }

  for (i in seq_along(assay_names)) {
    assay <- assay_names[i]
    modality <- mod_order[i]
    assay_path <- paste0("assays/", assay)

    if (!src$exists(assay_path)) next
    a_grp <- src[[assay_path]]

    if (verbose) message("  Streaming assay '", assay, "' -> modality '",
                         modality, "'...")

    mod_grp <- mod_root$create_group(modality)

    # --- features -> var/_index ---
    feature_names <- character(0)
    if (a_grp$exists("features")) {
      feature_names <- as.character(a_grp[["features"]]$read())
    }

    # var group (empty DataFrame with just _index)
    var_grp <- mod_grp$create_group("var")
    var_grp$create_attr(attr_name = "encoding-type", robj = "dataframe",
                        dtype = CachedGuessDType("dataframe"),
                        space = ScalarSpace())
    var_grp$create_attr(attr_name = "encoding-version", robj = "0.2.0",
                        dtype = CachedGuessDType("0.2.0"), space = ScalarSpace())
    var_grp$create_attr(attr_name = "_index", robj = "_index",
                        dtype = CachedGuessDType("_index"), space = ScalarSpace())
    var_grp$create_attr(attr_name = "column-order", robj = character(0),
                        dtype = CachedUtf8Type())
    if (length(feature_names) > 0) {
      var_grp$create_dataset("_index", robj = feature_names,
                             dtype = CachedUtf8Type(),
                             chunk_dims = length(feature_names),
                             gzip_level = gzip)
      AddAnndataEncoding(var_grp[["_index"]], encoding_type = "string-array")
    }

    # --- X (layers/data -> CSR, or scale.data or data) ---
    data_src <- NULL
    # Seurat v5 paths
    if (a_grp$exists("layers") && a_grp[["layers"]]$exists("data")) {
      data_src <- a_grp[["layers"]][["data"]]
    } else if (a_grp$exists("data")) {
      data_src <- a_grp[["data"]]
    }

    if (!is.null(data_src)) {
      if (inherits(data_src, "H5Group")) {
        # CSC [genes x cells] -> CSR [cells x genes] (zero-copy)
        .stream_sparse_group(data_src, mod_grp, "X", gzip,
                             src_format = "h5seurat", dst_format = "h5ad")
      } else if (inherits(data_src, "H5D")) {
        .stream_dense_matrix(data_src, mod_grp, "X", gzip, transpose = TRUE)
      }
    }

    # --- layers/counts -> layers/counts ---
    counts_src <- NULL
    if (a_grp$exists("layers") && a_grp[["layers"]]$exists("counts")) {
      counts_src <- a_grp[["layers"]][["counts"]]
    } else if (a_grp$exists("counts")) {
      counts_src <- a_grp[["counts"]]
    }

    if (!is.null(counts_src)) {
      layers_grp <- .h5ad_write_dict_group(mod_grp, "layers")
      if (inherits(counts_src, "H5Group")) {
        .stream_sparse_group(counts_src, layers_grp, "counts", gzip,
                             src_format = "h5seurat", dst_format = "h5ad")
      } else if (inherits(counts_src, "H5D")) {
        .stream_dense_matrix(counts_src, layers_grp, "counts", gzip,
                             transpose = TRUE)
      }
    }

    # --- obs (per-modality: just cell names) ---
    obs_grp <- mod_grp$create_group("obs")
    obs_grp$create_attr(attr_name = "encoding-type", robj = "dataframe",
                        dtype = CachedGuessDType("dataframe"),
                        space = ScalarSpace())
    obs_grp$create_attr(attr_name = "encoding-version", robj = "0.2.0",
                        dtype = CachedGuessDType("0.2.0"), space = ScalarSpace())
    obs_grp$create_attr(attr_name = "_index", robj = "_index",
                        dtype = CachedGuessDType("_index"), space = ScalarSpace())
    obs_grp$create_attr(attr_name = "column-order", robj = character(0),
                        dtype = CachedUtf8Type())
    if (n_cells > 0) {
      obs_grp$create_dataset("_index", robj = cell_names,
                             dtype = CachedUtf8Type(),
                             chunk_dims = n_cells, gzip_level = gzip)
      AddAnndataEncoding(obs_grp[["_index"]], encoding_type = "string-array")
    }

    # --- reductions -> obsm ---
    assay_reducs <- names(reduc_assay_map)[vapply(reduc_assay_map,
                                                   identical, logical(1), assay)]
    if (length(assay_reducs) > 0) {
      obsm_grp <- .h5ad_write_dict_group(mod_grp, "obsm")
      for (rn in assay_reducs) {
        red_path <- paste0("reductions/", rn)
        if (!src$exists(red_path)) next
        tryCatch({
          red_grp <- src[[red_path]]
          if (red_grp$exists("cell.embeddings")) {
            emb <- red_grp[["cell.embeddings"]]$read()
            key_name <- paste0("X_", rn)
            # h5seurat stores [n_cells, n_comp] in R. HDF5 on disk is
            # Fortran order [n_comp, n_cells]. AnnData expects C order
            # [n_cells, n_comp]. The transposition matches scanpy convention.
            emb_t <- if (is.matrix(emb)) t(emb) else emb
            obsm_grp$create_dataset(key_name, robj = emb_t,
                                    chunk_dims = dim(emb_t),
                                    gzip_level = gzip)
            obsm_grp[[key_name]]$create_attr(
              attr_name = "encoding-type", robj = "array",
              dtype = CachedGuessDType("array"), space = ScalarSpace())
            obsm_grp[[key_name]]$create_attr(
              attr_name = "encoding-version", robj = "0.2.0",
              dtype = CachedGuessDType("0.2.0"), space = ScalarSpace())
          }
        }, error = function(e) {
          if (verbose) warning("Could not stream reduction ", rn,
                               ": ", e$message, immediate. = TRUE)
        })
      }
    }

    # --- graphs -> obsp (strip assay_ prefix) ---
    assay_graphs <- names(graph_assay_map)[vapply(graph_assay_map,
                                                   identical, logical(1), assay)]
    if (length(assay_graphs) > 0) {
      obsp_grp <- .h5ad_write_dict_group(mod_grp, "obsp")
      for (gn in assay_graphs) {
        g_path <- paste0("graphs/", gn)
        if (!src$exists(g_path)) next
        # Strip assay prefix for h5mu name
        clean_name <- sub(paste0("^", assay, "_"), "", gn)
        tryCatch({
          g_src <- src[[g_path]]
          if (inherits(g_src, "H5Group")) {
            .stream_sparse_group(g_src, obsp_grp, clean_name, gzip,
                                 src_format = "h5seurat", dst_format = "h5ad")
          }
        }, error = function(e) {
          if (verbose) warning("Could not stream graph ", gn,
                               ": ", e$message, immediate. = TRUE)
        })
      }
    }

    # Ensure empty dict groups for required h5mu/anndata groups
    for (empty_name in c("obsm", "obsp", "varm", "varp", "layers", "uns")) {
      if (!mod_grp$exists(empty_name)) {
        .h5ad_write_dict_group(mod_grp, empty_name)
      }
    }

    # AnnData encoding on modality group
    mod_grp$create_attr(attr_name = "encoding-type", robj = "anndata",
                        dtype = CachedGuessDType("anndata"),
                        space = ScalarSpace())
    mod_grp$create_attr(attr_name = "encoding-version", robj = "0.1.0",
                        dtype = CachedGuessDType("0.1.0"), space = ScalarSpace())
  }

  # --- Global obs (meta.data -> /obs) ---
  if (verbose) message("  Streaming meta.data -> global obs...")
  if (src$exists("meta.data")) {
    .stream_h5seurat_md_to_h5ad_df(src[["meta.data"]], dst, "obs",
                                   cell_names, gzip, verbose = verbose)
  } else {
    # Write empty obs DataFrame
    WriteDFGroup(dst, "obs", data.frame(row.names = cell_names),
                 cell_names, gzip = gzip)
  }

  # Global obsm (empty dict)
  .h5ad_write_dict_group(dst, "obsm")

  dst$flush()
  if (verbose) message("Streaming h5seurat -> h5mu complete: ", h5mu_path)
  invisible(h5mu_path)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Core streaming converter: loom -> h5seurat
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Stream loom to h5seurat (direct HDF5-to-HDF5)
#'
#' Reads loom dense matrix, converts to sparse, and writes h5seurat layout.
#' Loom stores genes x cells (dense), h5seurat stores genes x cells (CSC sparse).
#'
#' @param loom_path Path to input loom file
#' @param h5seurat_path Path to output h5seurat file
#' @param assay Character assay name (default: "RNA")
#' @param gzip Integer gzip compression level
#' @param verbose Logical
#'
#' @keywords internal
#'
.stream_loom_to_h5seurat <- function(loom_path, h5seurat_path, assay = "RNA",
                                     gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for loom -> h5seurat streaming", call. = FALSE)

  src <- hdf5r::H5File$new(loom_path, mode = "r")
  on.exit(src$close_all(), add = TRUE)
  dst <- hdf5r::H5File$new(h5seurat_path, mode = "w")
  on.exit(dst$close_all(), add = TRUE)

  # Detect assay name from loom metadata if present
  if (src$attr_exists("SEURAT_ASSAY")) {
    assay <- hdf5r::h5attr(src, "SEURAT_ASSAY")
  }

  # --- Read cell and feature names ---
  cell_names <- NULL
  if (src$exists("col_attrs") && src[["col_attrs"]]$exists("CellID")) {
    cell_names <- as.character(src[["col_attrs"]][["CellID"]]$read())
  } else if (src$exists("col_attrs")) {
    # Try first available string array in col_attrs
    ca <- src[["col_attrs"]]
    for (n in names(ca)) {
      vals <- tryCatch(as.character(ca[[n]]$read()), error = function(e) NULL)
      if (!is.null(vals) && is.character(vals)) {
        cell_names <- vals
        break
      }
    }
  }

  feature_names <- NULL
  if (src$exists("row_attrs") && src[["row_attrs"]]$exists("Gene")) {
    feature_names <- as.character(src[["row_attrs"]][["Gene"]]$read())
  } else if (src$exists("row_attrs")) {
    ra <- src[["row_attrs"]]
    for (n in names(ra)) {
      vals <- tryCatch(as.character(ra[[n]]$read()), error = function(e) NULL)
      if (!is.null(vals) && is.character(vals)) {
        feature_names <- vals
        break
      }
    }
  }

  if (is.null(cell_names)) cell_names <- character(0)
  if (is.null(feature_names)) feature_names <- character(0)
  n_cells <- length(cell_names)
  n_features <- length(feature_names)

  # Write cell.names
  if (n_cells > 0) {
    dst$create_dataset("cell.names", robj = cell_names,
                       dtype = CachedUtf8Type(),
                       chunk_dims = n_cells, gzip_level = gzip)
  }

  # Create assay structure
  assays_grp <- dst$create_group("assays")
  a_grp <- assays_grp$create_group(assay)
  a_grp$create_attr(attr_name = "key", robj = paste0(tolower(assay), "_"),
                    dtype = CachedGuessDType(paste0(tolower(assay), "_")),
                    space = ScalarSpace())
  a_grp$create_attr(attr_name = "s4class", robj = "SeuratObject::Assay5",
                    dtype = CachedGuessDType("SeuratObject::Assay5"),
                    space = ScalarSpace())

  # features
  if (n_features > 0) {
    a_grp$create_dataset("features", robj = feature_names,
                         dtype = CachedUtf8Type(),
                         chunk_dims = n_features, gzip_level = gzip)
  }

  layers_grp <- a_grp$create_group("layers")
  chunk_size <- 65536L

  # --- /matrix -> layers/data ---
  # Loom /matrix is genes x cells (dense). Convert to sparse CSC for h5seurat.
  if (src$exists("matrix")) {
    if (verbose) message("  Reading loom /matrix (genes x cells dense)...")
    mat <- src[["matrix"]]$read()
    # hdf5r reverses HDF5 C-order [genes, cells] → R [cells, genes].
    # Transpose to get R [genes, cells] for h5seurat CSC format.
    mat <- t(mat)

    if (verbose) message("  Converting to sparse CSC...")
    sp <- as(mat, "dgCMatrix")
    rm(mat)

    data_grp <- layers_grp$create_group("data")
    data_grp$create_dataset("data", robj = sp@x,
                            chunk_dims = min(length(sp@x), chunk_size),
                            gzip_level = gzip)
    data_grp$create_dataset("indices", robj = sp@i,
                            chunk_dims = min(length(sp@i), chunk_size),
                            gzip_level = gzip)
    data_grp$create_dataset("indptr", robj = sp@p,
                            chunk_dims = min(length(sp@p), chunk_size),
                            gzip_level = gzip)
    data_grp$create_attr(attr_name = "dims", robj = dim(sp),
                         dtype = GuessDType(dim(sp)))
    rm(sp)
  }

  # --- /layers -> additional layers ---
  if (src$exists("layers")) {
    loom_layers <- src[["layers"]]
    for (layer_name in names(loom_layers)) {
      tryCatch({
        if (verbose) message("  Streaming layer: ", layer_name)
        layer_mat <- t(loom_layers[[layer_name]]$read())  # transpose: hdf5r reverses dims
        sp <- as(layer_mat, "dgCMatrix")
        rm(layer_mat)

        # Map common loom layer names to h5seurat
        h5s_layer_name <- switch(tolower(layer_name),
          "spliced" = layer_name,
          "unspliced" = layer_name,
          layer_name
        )

        l_grp <- layers_grp$create_group(h5s_layer_name)
        l_grp$create_dataset("data", robj = sp@x,
                             chunk_dims = min(length(sp@x), chunk_size),
                             gzip_level = gzip)
        l_grp$create_dataset("indices", robj = sp@i,
                             chunk_dims = min(length(sp@i), chunk_size),
                             gzip_level = gzip)
        l_grp$create_dataset("indptr", robj = sp@p,
                             chunk_dims = min(length(sp@p), chunk_size),
                             gzip_level = gzip)
        l_grp$create_attr(attr_name = "dims", robj = dim(sp),
                          dtype = GuessDType(dim(sp)))
        rm(sp)
      }, error = function(e) {
        if (verbose) warning("Could not stream layer ", layer_name, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- col_attrs -> meta.data ---
  if (verbose) message("  Streaming col_attrs -> meta.data...")
  md_grp <- dst$create_group("meta.data")
  col_written <- character(0)
  if (src$exists("col_attrs")) {
    ca <- src[["col_attrs"]]
    for (col in names(ca)) {
      if (col == "CellID") next  # already used as cell.names
      tryCatch({
        col_obj <- ca[[col]]
        if (!inherits(col_obj, "H5D")) next
        vals <- col_obj$read()
        # Skip multi-dimensional arrays (e.g. embeddings stored in col_attrs)
        if (!is.null(dim(vals)) && length(dim(vals)) > 1) next
        if (is.character(vals)) {
          md_grp$create_dataset(col, robj = vals, dtype = CachedUtf8Type(),
                                chunk_dims = length(vals), gzip_level = gzip)
        } else if (is.logical(vals)) {
          md_grp$create_dataset(col, robj = as.integer(vals),
                                chunk_dims = length(vals), gzip_level = gzip)
        } else {
          md_grp$create_dataset(col, robj = vals,
                                chunk_dims = length(vals), gzip_level = gzip)
        }
        col_written <- c(col_written, col)
      }, error = function(e) {
        if (verbose) warning("Could not stream col_attr ", col, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }
  if (length(col_written) > 0) {
    md_grp$create_attr(attr_name = "colnames", robj = col_written,
                       dtype = GuessDType(col_written))
  }

  # --- col_graphs -> graphs ---
  if (verbose) message("  Streaming col_graphs -> graphs...")
  graphs_grp <- dst$create_group("graphs")
  if (src$exists("col_graphs")) {
    cg <- src[["col_graphs"]]
    for (gn in names(cg)) {
      tryCatch({
        g_grp <- cg[[gn]]
        if (!inherits(g_grp, "H5Group")) next
        if (!g_grp$exists("a") || !g_grp$exists("b") || !g_grp$exists("w")) next

        # COO format (0-based a, b indices + weights)
        a_idx <- g_grp[["a"]]$read()
        b_idx <- g_grp[["b"]]$read()
        w_vals <- g_grp[["w"]]$read()

        # Convert COO to CSC (dgCMatrix)
        mat <- Matrix::sparseMatrix(
          i = as.integer(a_idx) + 1L,
          j = as.integer(b_idx) + 1L,
          x = as.numeric(w_vals),
          dims = c(n_cells, n_cells),
          repr = "C"  # CSC
        )

        graph_dst <- graphs_grp$create_group(gn)
        graph_dst$create_dataset("data", robj = mat@x,
                                 chunk_dims = min(length(mat@x), chunk_size),
                                 gzip_level = gzip)
        graph_dst$create_dataset("indices", robj = mat@i,
                                 chunk_dims = min(length(mat@i), chunk_size),
                                 gzip_level = gzip)
        graph_dst$create_dataset("indptr", robj = mat@p,
                                 chunk_dims = min(length(mat@p), chunk_size),
                                 gzip_level = gzip)
        graph_dst$create_attr(attr_name = "dims", robj = dim(mat),
                              dtype = GuessDType(dim(mat)))
        graph_dst$create_attr(attr_name = "assay.used", robj = assay,
                              dtype = CachedGuessDType(assay),
                              space = ScalarSpace())
      }, error = function(e) {
        if (verbose) warning("Could not stream graph ", gn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- reductions (loom scConvert extension or col_attrs embeddings) ---
  reductions_grp <- dst$create_group("reductions")
  if (src$exists("reductions")) {
    red_root <- src[["reductions"]]
    for (rn in names(red_root)) {
      tryCatch({
        red_src <- red_root[[rn]]
        if (!inherits(red_src, "H5Group")) next

        r_grp <- reductions_grp$create_group(rn)
        r_grp$create_attr(attr_name = "active.assay", robj = assay,
                          dtype = CachedGuessDType(assay),
                          space = ScalarSpace())
        key <- AnnDataReductionKey(rn)
        r_grp$create_attr(attr_name = "key", robj = key,
                          dtype = CachedGuessDType(key), space = ScalarSpace())
        r_grp$create_attr(attr_name = "global", robj = 0L,
                          dtype = GuessDType(0L), space = ScalarSpace())
        r_grp$create_group("misc")

        if (red_src$exists("embeddings")) {
          emb <- red_src[["embeddings"]]$read()
          # Loom stores embeddings transposed (writeLoom convention).
          # hdf5r reverses: loom HDF5 C [n_cells, n_comp] → R [n_comp, n_cells].
          # h5seurat needs R [n_cells, n_comp], so transpose.
          emb <- t(emb)
          r_grp$create_dataset("cell.embeddings", robj = emb,
                               chunk_dims = dim(emb), gzip_level = gzip)
        }
      }, error = function(e) {
        if (verbose) warning("Could not stream reduction ", rn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # Root scaffold (attrs, empty groups, active.ident)
  .h5seurat_write_scaffold(dst, assay, n_cells)

  dst$flush()
  if (verbose) message("Streaming loom -> h5seurat complete: ", h5seurat_path)
  invisible(h5seurat_path)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Core streaming converter: h5seurat -> loom
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Stream h5seurat to loom (direct HDF5-to-HDF5)
#'
#' Reads h5seurat sparse matrix (CSC), densifies, transposes to genes x cells,
#' and writes to loom format.
#'
#' @param h5seurat_path Path to input h5seurat file
#' @param loom_path Path to output loom file
#' @param gzip Integer gzip compression level
#' @param verbose Logical
#'
#' @keywords internal
#'
.stream_h5seurat_to_loom <- function(h5seurat_path, loom_path,
                                     gzip = 4L, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE))
    stop("hdf5r required for h5seurat -> loom streaming", call. = FALSE)

  src <- hdf5r::H5File$new(h5seurat_path, mode = "r")
  on.exit(src$close_all(), add = TRUE)
  dst <- hdf5r::H5File$new(loom_path, mode = "w")
  on.exit(dst$close_all(), add = TRUE)

  # Detect active assay
  assay <- "RNA"
  if (src$attr_exists("active.assay")) {
    assay <- hdf5r::h5attr(src, "active.assay")
  }
  assay_path <- paste0("assays/", assay)

  # Read cell names
  cell_names <- if (src$exists("cell.names")) {
    as.character(src[["cell.names"]]$read())
  } else {
    character(0)
  }
  n_cells <- length(cell_names)

  # Read feature names
  feature_names <- character(0)
  if (src$exists(assay_path) && src[[assay_path]]$exists("features")) {
    feature_names <- as.character(src[[assay_path]][["features"]]$read())
  }
  n_features <- length(feature_names)

  # Helper to read a sparse group from h5seurat and return as dense
  .read_h5seurat_sparse_as_dense <- function(grp) {
    if (!grp$exists("data") || !grp$exists("indices") || !grp$exists("indptr"))
      return(NULL)
    data_vals <- grp[["data"]]$read()
    indices <- grp[["indices"]]$read()
    indptr <- grp[["indptr"]]$read()
    dims <- NULL
    if (grp$attr_exists("dims")) dims <- hdf5r::h5attr(grp, "dims")
    if (is.null(dims)) dims <- c(n_features, n_cells)

    # CSC format: indices are row indices, indptr are column pointers
    mat <- new("dgCMatrix",
      i = as.integer(indices),
      p = as.integer(indptr),
      x = as.numeric(data_vals),
      Dim = as.integer(dims)
    )
    as.matrix(mat)
  }

  chunk_size <- 65536L

  # --- layers/data -> /matrix (dense, genes x cells) ---
  if (verbose) message("  Reading h5seurat expression data...")
  dense_mat <- NULL

  if (src$exists(assay_path)) {
    a_grp <- src[[assay_path]]
    # Try v5 path first, then v4
    data_path <- NULL
    if (a_grp$exists("layers") && a_grp[["layers"]]$exists("data")) {
      data_path <- "layers/data"
    } else if (a_grp$exists("data")) {
      data_path <- "data"
    }

    if (!is.null(data_path)) {
      data_obj <- a_grp[[data_path]]
      if (inherits(data_obj, "H5Group")) {
        dense_mat <- .read_h5seurat_sparse_as_dense(data_obj)
      } else if (inherits(data_obj, "H5D")) {
        dense_mat <- data_obj$read()
      }
    }
  }

  if (!is.null(dense_mat)) {
    if (verbose) message("  Writing /matrix (", nrow(dense_mat), " x ",
                         ncol(dense_mat), " dense)...")
    # Loom /matrix is [genes x cells] in HDF5 C-order.
    # hdf5r reverses R dims when writing, so we transpose: R [cells x genes]
    # → HDF5 C-order [genes x cells].
    loom_mat <- t(dense_mat)
    dst$create_dataset("matrix", robj = loom_mat,
                       chunk_dims = dim(loom_mat), gzip_level = gzip)
  }

  # --- Additional layers ---
  dst_layers <- dst$create_group("layers")
  if (src$exists(assay_path)) {
    a_grp <- src[[assay_path]]
    # Check for counts
    counts_path <- NULL
    if (a_grp$exists("layers") && a_grp[["layers"]]$exists("counts")) {
      counts_path <- "layers/counts"
    } else if (a_grp$exists("counts")) {
      counts_path <- "counts"
    }

    if (!is.null(counts_path)) {
      tryCatch({
        counts_obj <- a_grp[[counts_path]]
        counts_mat <- if (inherits(counts_obj, "H5Group")) {
          .read_h5seurat_sparse_as_dense(counts_obj)
        } else {
          counts_obj$read()
        }
        if (!is.null(counts_mat)) {
          if (verbose) message("  Writing /layers/counts...")
          # Transpose for loom convention (same as /matrix)
          loom_counts <- t(counts_mat)
          dst_layers$create_dataset("counts", robj = loom_counts,
                                    chunk_dims = dim(loom_counts),
                                    gzip_level = gzip)
        }
      }, error = function(e) {
        if (verbose) warning("Could not write counts layer: ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- col_attrs (cell metadata + CellID) ---
  if (verbose) message("  Writing col_attrs...")
  ca_grp <- dst$create_group("col_attrs")
  if (n_cells > 0) {
    ca_grp$create_dataset("CellID", robj = cell_names,
                          dtype = CachedUtf8Type(),
                          chunk_dims = n_cells, gzip_level = gzip)
  }

  # meta.data columns -> col_attrs
  if (src$exists("meta.data")) {
    md <- src[["meta.data"]]
    md_cols <- character(0)
    if (md$attr_exists("colnames")) {
      md_cols <- hdf5r::h5attr(md, "colnames")
    } else {
      md_cols <- names(md)
    }

    for (col in md_cols) {
      if (!md$exists(col)) next
      tryCatch({
        col_obj <- md[[col]]
        if (inherits(col_obj, "H5Group")) {
          # Factor: reconstruct as string array for loom
          if (col_obj$exists("levels") && col_obj$exists("values")) {
            levs <- as.character(col_obj[["levels"]]$read())
            vals <- col_obj[["values"]]$read()
            str_vals <- levs[vals]  # 1-based indexing
            ca_grp$create_dataset(col, robj = str_vals,
                                  dtype = CachedUtf8Type(),
                                  chunk_dims = length(str_vals),
                                  gzip_level = gzip)
          }
        } else if (inherits(col_obj, "H5D")) {
          vals <- col_obj$read()
          if (!is.null(dim(vals)) && length(dim(vals)) > 1) next
          if (is.character(vals)) {
            ca_grp$create_dataset(col, robj = vals, dtype = CachedUtf8Type(),
                                  chunk_dims = length(vals), gzip_level = gzip)
          } else {
            ca_grp$create_dataset(col, robj = vals,
                                  chunk_dims = length(vals), gzip_level = gzip)
          }
        }
      }, error = function(e) {
        if (verbose) warning("Could not write col_attr ", col, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- row_attrs (feature metadata + Gene) ---
  if (verbose) message("  Writing row_attrs...")
  ra_grp <- dst$create_group("row_attrs")
  if (n_features > 0) {
    ra_grp$create_dataset("Gene", robj = feature_names,
                          dtype = CachedUtf8Type(),
                          chunk_dims = n_features, gzip_level = gzip)
  }

  # --- col_graphs (graphs -> COO format) ---
  if (verbose) message("  Writing col_graphs...")
  cg_grp <- dst$create_group("col_graphs")
  dst$create_group("row_graphs")  # Required loom structure
  if (src$exists("graphs")) {
    g_root <- src[["graphs"]]
    for (gn in names(g_root)) {
      tryCatch({
        g_src <- g_root[[gn]]
        if (!inherits(g_src, "H5Group")) next
        if (!g_src$exists("data") || !g_src$exists("indices") ||
            !g_src$exists("indptr")) next

        data_vals <- g_src[["data"]]$read()
        indices <- g_src[["indices"]]$read()
        indptr <- g_src[["indptr"]]$read()
        dims <- if (g_src$attr_exists("dims")) {
          hdf5r::h5attr(g_src, "dims")
        } else {
          c(n_cells, n_cells)
        }

        # CSC to COO: expand column pointers to get (row, col, val) triplets
        n_cols_g <- length(indptr) - 1L
        col_idx <- rep(seq_len(n_cols_g) - 1L, diff(as.integer(indptr)))

        g_dst <- cg_grp$create_group(gn)
        g_dst$create_dataset("a", robj = as.integer(indices),  # row indices (0-based)
                             chunk_dims = min(length(indices), chunk_size),
                             gzip_level = gzip)
        g_dst$create_dataset("b", robj = as.integer(col_idx),  # col indices (0-based)
                             chunk_dims = min(length(col_idx), chunk_size),
                             gzip_level = gzip)
        g_dst$create_dataset("w", robj = as.numeric(data_vals),
                             chunk_dims = min(length(data_vals), chunk_size),
                             gzip_level = gzip)
      }, error = function(e) {
        if (verbose) warning("Could not stream graph ", gn, ": ",
                             e$message, immediate. = TRUE)
      })
    }
  }

  # --- reductions (scConvert extension) ---
  if (src$exists("reductions")) {
    red_root <- src[["reductions"]]
    red_names <- names(red_root)
    if (length(red_names) > 0) {
      dst_red <- dst$create_group("reductions")
      for (rn in red_names) {
        tryCatch({
          r_grp <- red_root[[rn]]
          if (!inherits(r_grp, "H5Group")) next
          if (!r_grp$exists("cell.embeddings")) next

          emb <- r_grp[["cell.embeddings"]]$read()
          # Transpose to match loom/writeLoom convention: store as
          # R [n_comp, n_cells] so hdf5r writes HDF5 C [n_cells, n_comp],
          # and readLoom's t() restores [n_cells, n_comp].
          emb <- t(emb)
          dr_grp <- dst_red$create_group(rn)
          dr_grp$create_dataset("embeddings", robj = emb,
                                chunk_dims = dim(emb), gzip_level = gzip)
        }, error = function(e) {
          if (verbose) warning("Could not stream reduction ", rn, ": ",
                               e$message, immediate. = TRUE)
        })
      }
    }
  }

  # Loom version and assay attributes
  dst$create_attr(attr_name = "LOOM_SPEC_VERSION", robj = "3.0.0",
                  dtype = CachedGuessDType("3.0.0"), space = ScalarSpace())
  dst$create_attr(attr_name = "SEURAT_ASSAY", robj = assay,
                  dtype = CachedGuessDType(assay), space = ScalarSpace())

  dst$flush()
  if (verbose) message("Streaming h5seurat -> loom complete: ", loom_path)
  invisible(loom_path)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Composite streaming helper
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Run two-step streaming conversion via a temp h5seurat file
#'
#' Creates a temp h5seurat file, streams source -> h5seurat, then
#' h5seurat -> dest. The temp file is deleted on exit.
#'
#' @param source Input file path
#' @param dest Output file path
#' @param step1_fn Function to convert source -> h5seurat (temp)
#' @param step2_fn Function to convert h5seurat (temp) -> dest
#' @param step1_args Named list of additional args for step1_fn
#' @param step2_args Named list of additional args for step2_fn
#' @param gzip Integer gzip compression level
#' @param verbose Logical
#'
#' @keywords internal
#'
.stream_composite_via_h5seurat <- function(source, dest,
                                           step1_fn, step2_fn,
                                           step1_args = list(),
                                           step2_args = list(),
                                           gzip = 4L, verbose = TRUE) {
  tmp <- tempfile(fileext = ".h5seurat")
  on.exit(unlink(tmp), add = TRUE)

  # Step 1: source -> temp h5seurat
  step1_call <- c(list(source, tmp, gzip = gzip, verbose = verbose), step1_args)
  do.call(step1_fn, step1_call)

  # Step 2: temp h5seurat -> dest
  step2_call <- c(list(tmp, dest, gzip = gzip, verbose = verbose), step2_args)
  do.call(step2_fn, step2_call)

  invisible(dest)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported converters: Core (h5seurat hub)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert an h5mu file to h5seurat format
#'
#' Converts a MuData h5mu file to h5seurat format. By default uses streaming
#' conversion that copies fields directly between HDF5 files without loading
#' into R memory. Each h5mu modality becomes a Seurat assay.
#'
#' @param source Path to input .h5mu file
#' @param dest Path for output .h5seurat file
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly. If FALSE,
#'   load as Seurat then save as h5seurat.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5MUToH5Seurat <- function(source, dest, overwrite = FALSE, gzip = 4L,
                            verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    .stream_h5mu_to_h5seurat(source, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readH5MU(source, verbose = verbose)
    writeH5Seurat(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert an h5seurat file to h5mu format
#'
#' Converts an h5seurat file to MuData h5mu format. By default uses streaming
#' conversion. Each Seurat assay becomes an h5mu modality.
#'
#' @param source Path to input .h5seurat file
#' @param dest Path for output .h5mu file
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly. If FALSE,
#'   load as Seurat then save as h5mu.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5SeuratToH5MU <- function(source, dest, overwrite = FALSE, gzip = 4L,
                            verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    .stream_h5seurat_to_h5mu(source, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readH5Seurat(source, verbose = verbose)
    writeH5MU(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert a loom file to h5seurat format
#'
#' Converts a loom file to h5seurat format. By default uses streaming
#' conversion. The loom dense matrix is converted to sparse CSC for h5seurat.
#'
#' @param source Path to input .loom file
#' @param dest Path for output .h5seurat file
#' @param assay Assay name to use (default: "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly. If FALSE,
#'   load as Seurat then save as h5seurat.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
LoomToH5Seurat <- function(source, dest, assay = "RNA", overwrite = FALSE,
                            gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    .stream_loom_to_h5seurat(source, dest, assay = assay, gzip = gzip,
                             verbose = verbose)
  } else {
    obj <- readLoom(source, verbose = verbose)
    writeH5Seurat(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert an h5seurat file to loom format
#'
#' Converts an h5seurat file to loom format. By default uses streaming
#' conversion. Sparse matrices are densified for loom output.
#'
#' @param source Path to input .h5seurat file
#' @param dest Path for output .loom file
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), stream fields directly. If FALSE,
#'   load as Seurat then save as loom.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5SeuratToLoom <- function(source, dest, overwrite = FALSE, gzip = 4L,
                            verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    .stream_h5seurat_to_loom(source, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readH5Seurat(source, verbose = verbose)
    writeLoom(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported converters: Composite (via temp h5seurat)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert an h5mu file to h5ad format
#'
#' Converts via streaming: h5mu -> h5seurat (temp) -> h5ad.
#' Only the first modality (or \code{assay}) is exported to h5ad.
#'
#' @param source Path to input .h5mu file
#' @param dest Path for output .h5ad file
#' @param assay Assay name to export (default: first modality)
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5MUToH5AD <- function(source, dest, assay = NULL, overwrite = FALSE,
                        gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    .stream_h5mu_to_h5seurat(source, tmp, gzip = gzip, verbose = verbose)
    # Detect assay from the temp h5seurat if not specified
    if (is.null(assay)) {
      h5tmp <- hdf5r::H5File$new(tmp, mode = "r")
      assay <- tryCatch(hdf5r::h5attr(h5tmp, "active.assay"),
                        error = function(e) "RNA")
      h5tmp$close_all()
    }
    # h5seurat -> h5ad via existing infrastructure
    obj <- readH5Seurat(tmp, verbose = verbose)
    writeH5AD(obj, dest, assay = assay, overwrite = FALSE, verbose = verbose)
  } else {
    obj <- readH5MU(source, verbose = verbose)
    if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
    writeH5AD(obj, dest, assay = assay, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert an h5ad file to h5mu format
#'
#' Converts via streaming: h5ad -> h5seurat (temp) -> h5mu.
#' Creates a single-modality h5mu file.
#'
#' @param source Path to input .h5ad file
#' @param dest Path for output .h5mu file
#' @param assay Assay name (default: "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5ADToH5MU <- function(source, dest, assay = "RNA", overwrite = FALSE,
                        gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    # h5ad -> h5seurat via existing infrastructure, then stream to h5mu
    obj <- readH5AD(source, verbose = verbose)
    writeH5Seurat(obj, tmp, overwrite = FALSE, verbose = verbose)
    .stream_h5seurat_to_h5mu(tmp, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readH5AD(source, verbose = verbose)
    writeH5MU(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert a loom file to h5ad format
#'
#' Converts via streaming: loom -> h5seurat (temp) -> h5ad.
#'
#' @param source Path to input .loom file
#' @param dest Path for output .h5ad file
#' @param assay Assay name (default: "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
LoomToH5AD <- function(source, dest, assay = "RNA", overwrite = FALSE,
                        gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    .stream_loom_to_h5seurat(source, tmp, assay = assay, gzip = gzip,
                             verbose = verbose)
    # h5seurat -> h5ad via existing infrastructure
    obj <- readH5Seurat(tmp, verbose = verbose)
    writeH5AD(obj, dest, assay = assay, overwrite = FALSE, verbose = verbose)
  } else {
    obj <- readLoom(source, verbose = verbose)
    writeH5AD(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert an h5ad file to loom format
#'
#' Converts via streaming: h5ad -> h5seurat (temp) -> loom.
#'
#' @param source Path to input .h5ad file
#' @param dest Path for output .loom file
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5ADToLoom <- function(source, dest, overwrite = FALSE, gzip = 4L,
                        verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    # h5ad -> h5seurat via existing infrastructure, then stream to loom
    obj <- readH5AD(source, verbose = verbose)
    writeH5Seurat(obj, tmp, overwrite = FALSE, verbose = verbose)
    .stream_h5seurat_to_loom(tmp, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readH5AD(source, verbose = verbose)
    writeLoom(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert a loom file to h5mu format
#'
#' Converts via streaming: loom -> h5seurat (temp) -> h5mu.
#' Creates a single-modality h5mu file.
#'
#' @param source Path to input .loom file
#' @param dest Path for output .h5mu file
#' @param assay Assay name (default: "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
LoomToH5MU <- function(source, dest, assay = "RNA", overwrite = FALSE,
                        gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    .stream_loom_to_h5seurat(source, tmp, assay = assay, gzip = gzip,
                             verbose = verbose)
    .stream_h5seurat_to_h5mu(tmp, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readLoom(source, verbose = verbose)
    writeH5MU(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert an h5mu file to loom format
#'
#' Converts via streaming: h5mu -> h5seurat (temp) -> loom.
#' Exports the primary modality.
#'
#' @param source Path to input .h5mu file
#' @param dest Path for output .loom file
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5MUToLoom <- function(source, dest, overwrite = FALSE, gzip = 4L,
                        verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    .stream_h5mu_to_h5seurat(source, tmp, gzip = gzip, verbose = verbose)
    .stream_h5seurat_to_loom(tmp, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readH5MU(source, verbose = verbose)
    writeLoom(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported converters: Zarr composites
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert an h5mu file to zarr format
#'
#' Converts via streaming: h5mu -> h5seurat (temp) -> zarr.
#'
#' @param source Path to input .h5mu file
#' @param dest Path for output .zarr directory
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
H5MUToZarr <- function(source, dest, overwrite = FALSE, gzip = 4L,
                        verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    .stream_h5mu_to_h5seurat(source, tmp, gzip = gzip, verbose = verbose)
    # Detect assay from temp
    h5tmp <- hdf5r::H5File$new(tmp, mode = "r")
    tmp_assay <- tryCatch(hdf5r::h5attr(h5tmp, "active.assay"),
                          error = function(e) "RNA")
    h5tmp$close_all()
    .stream_h5seurat_to_zarr(tmp, dest, assay = tmp_assay, gzip = gzip,
                             verbose = verbose)
  } else {
    obj <- readH5MU(source, verbose = verbose)
    writeZarr(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert a zarr store to h5mu format
#'
#' Converts via streaming: zarr -> h5seurat (temp) -> h5mu.
#'
#' @param source Path to input .zarr directory
#' @param dest Path for output .h5mu file
#' @param assay Assay name (default: "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
ZarrToH5MU <- function(source, dest, assay = "RNA", overwrite = FALSE,
                        gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!dir.exists(source))
    stop("Zarr store not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    .stream_zarr_to_h5seurat(source, tmp, assay = assay, gzip = gzip,
                             verbose = verbose)
    .stream_h5seurat_to_h5mu(tmp, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readZarr(source, verbose = verbose)
    writeH5MU(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert a loom file to zarr format
#'
#' Converts via streaming: loom -> h5seurat (temp) -> zarr.
#'
#' @param source Path to input .loom file
#' @param dest Path for output .zarr directory
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
LoomToZarr <- function(source, dest, overwrite = FALSE, gzip = 4L,
                        verbose = TRUE, stream = TRUE) {
  if (!file.exists(source))
    stop("Input file not found: ", source, call. = FALSE)
  if (dir.exists(dest) || file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest, recursive = TRUE)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    .stream_loom_to_h5seurat(source, tmp, assay = "RNA", gzip = gzip,
                             verbose = verbose)
    .stream_h5seurat_to_zarr(tmp, dest, assay = "RNA", gzip = gzip,
                             verbose = verbose)
  } else {
    obj <- readLoom(source, verbose = verbose)
    writeZarr(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}

#' Convert a zarr store to loom format
#'
#' Converts via streaming: zarr -> h5seurat (temp) -> loom.
#'
#' @param source Path to input .zarr directory
#' @param dest Path for output .loom file
#' @param assay Assay name (default: "RNA")
#' @param overwrite If TRUE, overwrite existing output
#' @param gzip Gzip compression level (0-9)
#' @param verbose Show progress messages
#' @param stream If TRUE (default), use streaming. If FALSE, load into R.
#'
#' @return Invisibly returns \code{dest}
#'
#' @export
#'
ZarrToLoom <- function(source, dest, assay = "RNA", overwrite = FALSE,
                        gzip = 4L, verbose = TRUE, stream = TRUE) {
  if (!dir.exists(source))
    stop("Zarr store not found: ", source, call. = FALSE)
  if (file.exists(dest)) {
    if (!overwrite)
      stop("Output exists: ", dest, ". Use overwrite = TRUE.", call. = FALSE)
    unlink(dest)
  }
  if (stream) {
    tmp <- tempfile(fileext = ".h5seurat")
    on.exit(unlink(tmp), add = TRUE)
    .stream_zarr_to_h5seurat(source, tmp, assay = assay, gzip = gzip,
                             verbose = verbose)
    .stream_h5seurat_to_loom(tmp, dest, gzip = gzip, verbose = verbose)
  } else {
    obj <- readZarr(source, verbose = verbose)
    writeLoom(obj, dest, overwrite = FALSE, verbose = verbose)
  }
  invisible(dest)
}
