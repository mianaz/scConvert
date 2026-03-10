#' Save a Seurat object to an H5MU file
#'
#' Export multimodal Seurat objects to MuData (.h5mu) format. Uses scConvert's
#' native writer by default, which works with Seurat V5 Assay5 objects and
#' requires no external dependencies. Falls back to MuDataSeurat::WriteH5MU
#' if \code{use.mudataseurat = TRUE} is set explicitly.
#'
#' @param object A \code{Seurat} object with one or more assays
#' @param filename Path to output .h5mu file
#' @param assays Character vector of assay names to export. If NULL (default),
#'   exports all assays
#' @param overwrite Overwrite existing file
#' @param verbose Show progress messages
#' @param use.mudataseurat Logical; if TRUE, use MuDataSeurat::WriteH5MU instead
#'   of the native writer. Note: MuDataSeurat may not work with Seurat V5 Assay5.
#'   Default: FALSE.
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the filename
#'
#' @details
#' The h5mu format is designed for multimodal data and stores each Seurat assay
#' as a separate modality under \code{/mod/{modality_name}}. This function:
#' \itemize{
#'   \item Extracts each specified assay from the Seurat object
#'   \item Converts assays to modality structure
#'   \item Writes counts and data layers for each modality
#'   \item Preserves cell metadata in global /obs and per-modality obs
#'   \item Maintains dimensional reductions and graphs
#' }
#'
#' @section Assay Mapping:
#' By default, Seurat assay names are mapped to standard MuData modality names:
#' \itemize{
#'   \item RNA -> rna
#'   \item ADT -> prot
#'   \item ATAC -> atac
#'   \item Spatial -> spatial
#'   \item Other names are converted to lowercase
#' }
#'
#' @examples
#' \dontrun{
#' # Save a multimodal Seurat object (CITE-seq)
#' SaveH5MU(seurat_obj, "multimodal_data.h5mu")
#'
#' # Save specific assays only
#' SaveH5MU(seurat_obj, "rna_and_protein.h5mu", assays = c("RNA", "ADT"))
#'
#' # Include spatial data
#' SaveH5MU(visium_obj, "spatial_multimodal.h5mu")
#' }
#'
#' @importFrom Seurat Assays Images DefaultAssay Project
#' @importFrom SeuratObject Cells
#'
#' @export
#'
SaveH5MU <- function(object,
                     filename,
                     assays = NULL,
                     overwrite = FALSE,
                     verbose = TRUE,
                     use.mudataseurat = FALSE,
                     ...) {

  if (!inherits(object, "Seurat")) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }

  if (use.mudataseurat) {
    # Legacy path: delegate to MuDataSeurat
    if (!requireNamespace("MuDataSeurat", quietly = TRUE)) {
      stop(
        "Package 'MuDataSeurat' is required when use.mudataseurat = TRUE.\n",
        "Install it with: remotes::install_github('PMBio/MuDataSeurat')\n",
        "Or use the default native writer (use.mudataseurat = FALSE).",
        call. = FALSE
      )
    }
    if (file.exists(filename) && !overwrite) {
      stop("File already exists: ", filename,
           "\nSet overwrite = TRUE to replace it.", call. = FALSE)
    }
    if (file.exists(filename) && overwrite) {
      file.remove(filename)
    }
    get("WriteH5MU", envir = asNamespace("MuDataSeurat"))(object, filename, ...)
    return(invisible(filename))
  }

  # Native path: use SeuratToH5MU (no external dependencies, V5 compatible)
  SeuratToH5MU(
    object = object,
    filename = filename,
    assays = assays,
    overwrite = overwrite,
    verbose = verbose
  )
}


#' Get default assay to modality name mapping
#'
#' @param assays Character vector of assay names
#' @return Named character vector mapping assay names to modality names
#' @keywords internal
#'
GetDefaultAssayToModalityMapping <- function(assays) {
  # Reverse mapping from modality to assay
  default_map <- c(
    RNA = "rna",
    ADT = "prot",
    ATAC = "atac",
    Spatial = "spatial",
    HTO = "hto",
    SCT = "rna",  # SCTransform maps to rna by default
    Protein = "prot",
    Peaks = "atac"
  )

  mapping <- setNames(tolower(assays), assays)
  for (assay in assays) {
    if (assay %in% names(default_map)) {
      mapping[assay] <- default_map[assay]
    }
  }

  # Detect and resolve conflicts (e.g., both RNA and SCT mapping to "rna")
  modality_counts <- table(mapping)
  conflicts <- names(modality_counts)[modality_counts > 1]

  if (length(conflicts) > 0) {
    for (conflict in conflicts) {
      # Find all assays mapping to this conflicting modality
      conflicting_assays <- names(mapping)[mapping == conflict]

      # Special handling for RNA/SCT conflict
      if (conflict == "rna" && "SCT" %in% conflicting_assays && "RNA" %in% conflicting_assays) {
        # Keep RNA as "rna", rename SCT to "sct"
        mapping["SCT"] <- "sct"
        message("Note: Both RNA and SCT assays detected. Mapping SCT -> 'sct' modality to avoid conflict.")
      } else if (conflict == "prot" && length(conflicting_assays) > 1) {
        # Handle protein assay conflicts (ADT, Protein, etc.)
        for (i in seq_along(conflicting_assays)) {
          assay <- conflicting_assays[i]
          if (i > 1) {  # Keep first as "prot", rename others
            mapping[assay] <- paste0("prot_", tolower(assay))
            message("Note: Multiple protein assays detected. Mapping ", assay, " -> '", mapping[assay], "' to avoid conflict.")
          }
        }
      } else {
        # Generic conflict resolution: append assay name
        for (i in seq_along(conflicting_assays)) {
          assay <- conflicting_assays[i]
          if (i > 1) {  # Keep first, rename others
            mapping[assay] <- tolower(assay)
            message("Note: Modality conflict detected. Mapping ", assay, " -> '", mapping[assay], "'.")
          }
        }
      }
    }
  }

  return(mapping)
}


#' Validate Seurat object for multimodal export
#'
#' @param object Seurat object
#' @param assays Character vector of assays to validate
#' @param verbose Logical
#' @return List with 'valid' (logical) and 'message' (character) elements
#' @keywords internal
#'
#' @importFrom Seurat Assays
#' @importFrom SeuratObject Cells
#'
ValidateMultimodalObject <- function(object, assays, verbose = FALSE) {

  # Check that all assays have the same cells
  cell_names <- Cells(object)
  n_cells <- length(cell_names)

  for (assay in assays) {
    assay_obj <- object[[assay]]
    assay_cells <- colnames(assay_obj)

    if (length(assay_cells) != n_cells) {
      return(list(
        valid = FALSE,
        message = paste0(
          "Assay '", assay, "' has ", length(assay_cells),
          " cells but object has ", n_cells, " cells. ",
          "All assays must have the same cells for h5mu export."
        )
      ))
    }

    if (!identical(assay_cells, cell_names)) {
      return(list(
        valid = FALSE,
        message = paste0(
          "Assay '", assay, "' has different cell names than the main object. ",
          "All assays must have matching cell names for h5mu export."
        )
      ))
    }
  }

  # Check for unique feature names across modalities (important for h5mu)
  all_features <- character()
  for (assay in assays) {
    features <- rownames(object[[assay]])
    duplicates <- intersect(all_features, features)
    if (length(duplicates) > 0 && verbose) {
      message("  Warning: ", length(duplicates),
              " features in '", assay, "' overlap with other assays. ",
              "MuData handles this, but be aware when analyzing.")
    }
    all_features <- c(all_features, features)
  }

  return(list(valid = TRUE, message = ""))
}


#' Prepare spatial data for h5mu export
#'
#' @param object Seurat object
#' @param assays Assay names
#' @param modality.names Modality name mapping
#' @param verbose Logical
#' @return Modified Seurat object with spatial metadata prepared
#' @keywords internal
#'
#' @importFrom Seurat Images
#'
PrepareSpatialForH5MU <- function(object, assays, modality.names, verbose = FALSE) {

  images <- Images(object)

  if (length(images) == 0) {
    return(object)
  }

  # Store spatial data in misc for MuDataSeurat to pick up
  for (img_name in images) {
    img_obj <- object[[img_name]]

    # Determine which modality this spatial data belongs to
    target_assay <- if ("Spatial" %in% assays) {
      "Spatial"
    } else {
      assays[1]
    }

    target_modality <- modality.names[target_assay]

    if (verbose) {
      message("  Associating spatial data '", img_name,
              "' with modality '", target_modality, "'")
    }

    # Store spatial metadata for later processing
    object@misc[[paste0("h5mu_spatial_", target_modality)]] <- list(
      image_name = img_name,
      image_obj = img_obj
    )
  }

  return(object)
}


#' Save a Seurat object to h5mu format (alias for SaveH5MU)
#'
#' @param object A Seurat object
#' @param filename Output filename
#' @param ... Additional arguments passed to SaveH5MU
#'
#' @return Invisibly returns the filename
#'
#' @rdname SaveH5MU
#' @export
#'
as.h5mu <- function(object, filename, ...) {
  SaveH5MU(object = object, filename = filename, ...)
}
