#' @include Convert.R
NULL


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

  # Store spatial data in misc for later processing
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


#' Save a Seurat object to h5mu format
#'
#' @param object A Seurat object
#' @param filename Output filename
#' @param ... Additional arguments passed to \code{\link{writeH5MU}}
#'
#' @return Invisibly returns the filename
#'
#' @seealso \code{\link{writeH5MU}}
#' @export
#'
as.h5mu <- function(object, filename, ...) {
  writeH5MU(object = object, filename = filename, ...)
}
