#' @include zzz.R
#' @importFrom hdf5r h5types H5File
#' @importFrom Seurat GetTissueCoordinates
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Field-of-View (FOV) spatial object helpers
#
# FOV objects (Seurat 5+) store imaging-based spatial transcriptomics data for
# Xenium, CosMx, MERFISH, Vizgen. Their boundaries slot holds Segmentation or
# Centroids objects; their molecules slot holds per-gene SpatialPoints.
#
# scConvert serializes FOVs into the h5ad container under the convention:
#   /obsm/spatial                               (n_cells, 2) universal centroids
#   /uns/spatial/{library_id}/segmentation/     dict of boundary sub-groups
#   /uns/spatial/{library_id}/molecules/        dict of molecule sub-groups
# This layout is backward-compatible with squidpy and scanpy readers (which
# ignore unknown uns children) and lets scConvert round-trip FOVs without
# silent data loss.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Write a scalar string attribute on an H5 location using the idiom already
# used throughout SpatialConversion.R.
.sc_h5_str_attr <- function(h5_obj, name, value) {
  if (name %in% hdf5r::list.attributes(h5_obj)) {
    h5_obj$attr_delete(attr_name = name)
  }
  h5_obj$create_attr(
    attr_name = name, robj = value,
    dtype = CachedGuessDType(x = value), space = ScalarSpace()
  )
  invisible(NULL)
}

#' Extract centroid coordinates from an FOV object
#'
#' @param fov_obj FOV (or nested Centroids/Segmentation) object
#' @return numeric matrix (n_cells x 2) with colnames c("x","y"), or NULL
#' @keywords internal
ExtractFOVCentroids <- function(fov_obj) {
  if (!inherits(fov_obj, "FOV")) return(NULL)
  boundaries <- tryCatch(fov_obj@boundaries, error = function(e) NULL)
  if (is.null(boundaries) || length(boundaries) == 0) return(NULL)

  # Prefer a Centroids boundary; fall back to the flattened polygon centroid
  # per cell when only a Segmentation is present.
  for (bn in names(boundaries)) {
    b <- boundaries[[bn]]
    if (inherits(b, "Centroids")) {
      coords <- tryCatch(Seurat::GetTissueCoordinates(b),
                          error = function(e) NULL)
      if (is.null(coords)) next
      if (all(c("x", "y") %in% colnames(coords))) {
        m <- as.matrix(coords[, c("x", "y"), drop = FALSE])
        rownames(m) <- coords$cell %||% rownames(coords)
        return(m)
      }
    }
  }
  # Fall-back: mean vertex per cell from the first Segmentation
  for (bn in names(boundaries)) {
    b <- boundaries[[bn]]
    if (inherits(b, "Segmentation")) {
      coords <- tryCatch(Seurat::GetTissueCoordinates(b),
                          error = function(e) NULL)
      if (is.null(coords)) next
      if (!all(c("x", "y", "cell") %in% colnames(coords))) next
      agg <- stats::aggregate(coords[, c("x", "y")],
                              by = list(cell = coords$cell), FUN = mean)
      m <- as.matrix(agg[, c("x", "y")])
      rownames(m) <- as.character(agg$cell)
      return(m)
    }
  }
  NULL
}

#' Extract segmentation boundaries from an FOV object
#'
#' Walks \code{fov_obj@boundaries}; returns a list keyed by boundary name.
#' Each element is a list with fields \code{type} ("Segmentation" or
#' "Centroids"), \code{cell_ids}, \code{coords} (n_vertices x 2 numeric) and
#' \code{polygon_offsets} (CSR-style row pointers, length n_cells + 1; present
#' only for Segmentation).
#'
#' @param fov_obj FOV object
#' @return named list or NULL
#' @keywords internal
ExtractFOVSegmentation <- function(fov_obj) {
  if (!inherits(fov_obj, "FOV")) return(NULL)
  boundaries <- tryCatch(fov_obj@boundaries, error = function(e) NULL)
  if (is.null(boundaries) || length(boundaries) == 0) return(NULL)

  out <- list()
  for (bn in names(boundaries)) {
    b <- boundaries[[bn]]
    coords_df <- tryCatch(Seurat::GetTissueCoordinates(b),
                           error = function(e) NULL)
    if (is.null(coords_df) || nrow(coords_df) == 0) next
    if (!all(c("x", "y", "cell") %in% colnames(coords_df))) next

    if (inherits(b, "Segmentation")) {
      # Preserve input cell order via rle on the 'cell' column — polygons
      # already arrive grouped (one contiguous run per cell).
      r <- rle(as.character(coords_df$cell))
      cell_ids <- r$values
      lens <- as.integer(r$lengths)
      offsets <- c(0L, cumsum(lens))
      out[[bn]] <- list(
        type            = "Segmentation",
        cell_ids        = cell_ids,
        coords          = as.matrix(coords_df[, c("x", "y"), drop = FALSE]),
        polygon_offsets = as.integer(offsets)
      )
    } else if (inherits(b, "Centroids")) {
      out[[bn]] <- list(
        type     = "Centroids",
        cell_ids = as.character(coords_df$cell),
        coords   = as.matrix(coords_df[, c("x", "y"), drop = FALSE])
      )
    }
  }
  if (length(out) == 0) NULL else out
}

#' Extract molecule coordinates from an FOV object
#'
#' Returns a named list where each entry is a data frame of
#' \code{(x, y, gene)} transcript positions. The top-level name mirrors
#' \code{names(fov_obj@molecules)} (typically "molecules" for CosMx/Xenium).
#'
#' @param fov_obj FOV object
#' @return named list or NULL
#' @keywords internal
ExtractFOVMolecules <- function(fov_obj) {
  if (!inherits(fov_obj, "FOV")) return(NULL)
  mols_slot <- tryCatch(fov_obj@molecules, error = function(e) NULL)
  if (is.null(mols_slot) || length(mols_slot) == 0) return(NULL)

  out <- list()
  for (mn in names(mols_slot)) {
    m <- mols_slot[[mn]]
    if (inherits(m, "Molecules")) {
      df <- tryCatch(Seurat::GetTissueCoordinates(m),
                      error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) next
      gene_col <- if ("molecule" %in% colnames(df)) "molecule" else
                  if ("gene"     %in% colnames(df)) "gene"     else NA_character_
      if (is.na(gene_col)) next
      out[[mn]] <- data.frame(
        x    = as.numeric(df$x),
        y    = as.numeric(df$y),
        gene = as.character(df[[gene_col]]),
        stringsAsFactors = FALSE
      )
    } else if (is.data.frame(m) && all(c("x", "y") %in% colnames(m))) {
      gene_col <- if ("gene" %in% colnames(m)) m$gene else NA_character_
      out[[mn]] <- data.frame(
        x    = as.numeric(m$x),
        y    = as.numeric(m$y),
        gene = as.character(gene_col),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(out) == 0) NULL else out
}

#' Convert FOV object to an h5ad-compatible list structure
#'
#' Produced structure is consumed by \code{\link{WriteFOVToH5AD}}. Users do
#' not normally call this directly.
#'
#' @param fov_obj FOV object
#' @param verbose Print progress messages
#' @return named list with \code{centroids}, \code{segmentation},
#'   \code{molecules}, \code{technology}, \code{fov_name}
#' @keywords internal
ConvertFOVToH5AD <- function(fov_obj, verbose = TRUE) {
  if (!inherits(fov_obj, "FOV")) {
    stop("Object is not an FOV object", call. = FALSE)
  }
  result <- list()

  cents <- ExtractFOVCentroids(fov_obj)
  if (!is.null(cents)) {
    result$centroids <- cents
    if (verbose) message("    Found ", nrow(cents), " centroids")
  }

  seg <- ExtractFOVSegmentation(fov_obj)
  if (!is.null(seg)) {
    result$segmentation <- seg
    if (verbose) {
      n_total <- sum(vapply(seg, function(s) length(s$cell_ids), integer(1)))
      message("    Found ", length(seg), " boundary set(s), ",
              n_total, " cell(s) total")
    }
  }

  mols <- ExtractFOVMolecules(fov_obj)
  if (!is.null(mols)) {
    result$molecules <- mols
    if (verbose) {
      n_total <- sum(vapply(mols, nrow, integer(1)))
      message("    Found ", length(mols), " molecule set(s), ",
              n_total, " transcript(s) total")
    }
  }

  result$technology <- "imaging"
  result$fov_name   <- tryCatch(fov_obj@key, error = function(e) "fov")
  result
}

#' Write FOV segmentation and molecules into an open h5ad library group
#'
#' Writes \code{@boundaries} (as segmentation/) and \code{@molecules} (as
#' molecules/) into \code{lib_group}, which must be an already-open H5Group
#' handle at \code{/uns/spatial/{library_id}/}. Used by
#' \code{\link{SeuratSpatialToH5AD}} during Seurat -> h5ad conversion.
#'
#' @param fov_obj FOV object
#' @param lib_group hdf5r H5Group handle for the library
#' @param library_id Identifier string (for error messages)
#' @param verbose Print progress messages
#' @return invisible NULL
#' @export
WriteFOVToH5AD <- function(fov_obj, lib_group, library_id = "fov",
                            verbose = TRUE) {
  if (!inherits(fov_obj, "FOV")) return(invisible(NULL))
  fov_data <- ConvertFOVToH5AD(fov_obj, verbose = verbose)

  # Segmentation / Centroids boundaries
  if (!is.null(fov_data$segmentation)) {
    if (lib_group$exists("segmentation")) lib_group$link_delete("segmentation")
    seg_grp <- lib_group$create_group("segmentation")
    .sc_h5_str_attr(seg_grp, "encoding-type", "dict")
    .sc_h5_str_attr(seg_grp, "encoding-version", "0.1.0")

    for (bn in names(fov_data$segmentation)) {
      seg <- fov_data$segmentation[[bn]]
      bgrp <- seg_grp$create_group(bn)
      .sc_h5_str_attr(bgrp, "type", seg$type)

      bgrp$create_dataset(
        name = "cell_ids",
        robj = as.character(seg$cell_ids)
      )
      # coords (n_vtx, 2): transpose for hdf5r so Python reads (n_vtx, 2).
      bgrp$create_dataset(
        name = "coords",
        robj = t(as.matrix(seg$coords)),
        dtype = h5types$H5T_NATIVE_DOUBLE
      )
      if (!is.null(seg$polygon_offsets)) {
        bgrp$create_dataset(
          name = "polygon_offsets",
          robj = as.integer(seg$polygon_offsets),
          dtype = h5types$H5T_NATIVE_INT32
        )
      }
    }
    if (verbose) {
      message("Wrote ", length(fov_data$segmentation),
              " boundary set(s) for library '", library_id, "'")
    }
  }

  # Transcript molecules
  if (!is.null(fov_data$molecules)) {
    if (lib_group$exists("molecules")) lib_group$link_delete("molecules")
    mol_grp <- lib_group$create_group("molecules")
    .sc_h5_str_attr(mol_grp, "encoding-type", "dict")
    .sc_h5_str_attr(mol_grp, "encoding-version", "0.1.0")

    for (mn in names(fov_data$molecules)) {
      df <- fov_data$molecules[[mn]]
      sub <- mol_grp$create_group(mn)
      sub$create_dataset(name = "x",    robj = as.numeric(df$x),
                         dtype = h5types$H5T_NATIVE_DOUBLE)
      sub$create_dataset(name = "y",    robj = as.numeric(df$y),
                         dtype = h5types$H5T_NATIVE_DOUBLE)
      sub$create_dataset(name = "gene", robj = as.character(df$gene))
    }
    if (verbose) {
      message("Wrote ", length(fov_data$molecules),
              " molecule set(s) for library '", library_id, "'")
    }
  }

  invisible(NULL)
}

#' Read FOV structure back from an h5ad library group
#'
#' Inverse of \code{\link{WriteFOVToH5AD}}: reads segmentation/ and molecules/
#' subgroups from an already-open \code{uns/spatial/{library}/} H5Group, and
#' reconstitutes a \code{FOV} object. Returns NULL when neither subgroup is
#' present (non-FOV libraries such as classic Visium).
#'
#' @param lib_group hdf5r H5Group handle for the library
#' @param key FOV key (trailing underscore added if missing)
#' @param assay Assay name to attach to the FOV (default "Spatial")
#' @return an FOV object or NULL
#' @export
ReadFOVFromH5AD <- function(lib_group, key = "fov", assay = "Spatial") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) return(NULL)
  if (!exists("CreateFOV", where = asNamespace("SeuratObject"),
              inherits = FALSE)) return(NULL)

  has_seg <- tryCatch(lib_group$exists("segmentation"), error = function(e) FALSE)
  has_mol <- tryCatch(lib_group$exists("molecules"),    error = function(e) FALSE)
  if (!has_seg && !has_mol) return(NULL)

  boundaries_list <- list()
  if (isTRUE(has_seg)) {
    seg_grp <- lib_group[["segmentation"]]
    for (bn in names(seg_grp)) {
      bgrp <- seg_grp[[bn]]
      seg_type <- tryCatch(h5attr(bgrp, "type"),
                            error = function(e) "Segmentation")
      cell_ids <- as.character(bgrp[["cell_ids"]]$read())
      # coords on disk is (n_vtx, 2) after the hdf5r transpose; reading it
      # back into R gives a (2, n_vtx) matrix (column-major) which we
      # transpose to the familiar (n_vtx, 2) layout.
      coords_raw <- bgrp[["coords"]]$read()
      if (is.matrix(coords_raw) && nrow(coords_raw) == 2L) {
        coords <- t(coords_raw)
      } else {
        coords <- as.matrix(coords_raw)
      }
      colnames(coords) <- c("x", "y")

      if (seg_type == "Segmentation" && bgrp$exists("polygon_offsets")) {
        offsets <- as.integer(bgrp[["polygon_offsets"]]$read())
        # Expand per-polygon offsets back to a flat (x, y, cell) data.frame.
        n_polys <- length(cell_ids)
        lens <- diff(offsets)
        cell_col <- rep.int(cell_ids, times = lens)
        df <- data.frame(
          x    = as.numeric(coords[, 1]),
          y    = as.numeric(coords[, 2]),
          cell = cell_col,
          stringsAsFactors = FALSE
        )
        boundaries_list[[bn]] <- tryCatch(
          SeuratObject::CreateSegmentation(coords = df),
          error = function(e) NULL
        )
      } else {
        df <- data.frame(
          x    = as.numeric(coords[, 1]),
          y    = as.numeric(coords[, 2]),
          cell = as.character(cell_ids),
          stringsAsFactors = FALSE
        )
        boundaries_list[[bn]] <- tryCatch(
          SeuratObject::CreateCentroids(coords = df),
          error = function(e) NULL
        )
      }
    }
    boundaries_list <- Filter(Negate(is.null), boundaries_list)
  }

  molecules_list <- list()
  if (isTRUE(has_mol)) {
    mol_grp <- lib_group[["molecules"]]
    for (mn in names(mol_grp)) {
      sub <- mol_grp[[mn]]
      if (!sub$exists("x") || !sub$exists("y") || !sub$exists("gene")) next
      df <- data.frame(
        x    = as.numeric(sub[["x"]]$read()),
        y    = as.numeric(sub[["y"]]$read()),
        gene = as.character(sub[["gene"]]$read()),
        stringsAsFactors = FALSE
      )
      molecules_list[[mn]] <- tryCatch(
        SeuratObject::CreateMolecules(coords = df),
        error = function(e) NULL
      )
    }
    molecules_list <- Filter(Negate(is.null), molecules_list)
  }

  # An FOV requires at least one boundary (Centroids or Segmentation).
  if (length(boundaries_list) == 0) return(NULL)

  fov_key <- if (grepl("_$", key)) key else paste0(key, "_")
  fov <- tryCatch(
    SeuratObject::CreateFOV(
      coords    = boundaries_list,
      type      = if (inherits(boundaries_list[[1]], "Segmentation"))
                    "segmentation" else "centroids",
      molecules = if (length(molecules_list) > 0) molecules_list[[1]] else NULL,
      assay     = assay,
      key       = fov_key
    ),
    error = function(e) NULL
  )
  fov
}

# Historical alias kept for backward compatibility with any downstream code
# that imported scConvert:::H5ADToFOV before the FOV serialization contract
# was added. New code should use ReadFOVFromH5AD().
#
#' @keywords internal
H5ADToFOV <- function(spatial_data, key = "fov") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) return(NULL)
  centroids <- NULL
  if (!is.null(spatial_data$centroids) && is.matrix(spatial_data$centroids)) {
    coords <- spatial_data$centroids
    if (ncol(coords) >= 2L) {
      cell_names <- rownames(coords) %||% paste0("cell_", seq_len(nrow(coords)))
      centroids <- tryCatch(
        SeuratObject::CreateCentroids(coords = data.frame(
          x = coords[, 1], y = coords[, 2], cell = cell_names,
          stringsAsFactors = FALSE)),
        error = function(e) NULL
      )
    }
  }
  if (is.null(centroids)) return(NULL)
  tryCatch(
    SeuratObject::CreateFOV(coords = centroids, type = "centroids", key = key),
    error = function(e) NULL
  )
}
