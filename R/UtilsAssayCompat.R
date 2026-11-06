#' Assay Compatibility Utility Functions
#'
#' This file contains utility functions for handling compatibility between
#' Seurat V4 Assay objects and V5 Assay5 objects. These functions provide
#' unified interfaces for operations that differ between versions.
#'
#' @name UtilsAssayCompat
#' @keywords internal
NULL


#' Get assay data with V4/V5 compatibility
#'
#' Retrieves data from an Assay or Assay5 object, automatically using the
#' appropriate parameter name (layer for V5, slot for V4).
#'
#' @param object Assay or Assay5 object
#' @param layer_or_slot Name of layer (V5) or slot (V4) to retrieve
#' @param verbose Print diagnostic messages
#'
#' @return Matrix or sparse matrix of assay data
#' @keywords internal
#'
#' @seealso V5Compatibility.R for the original implementation
#'
GetAssayDataCompat <- function(object, layer_or_slot, verbose = FALSE) {
  # In SeuratObject 5.0+, the 'slot' argument is defunct - always use 'layer'
  # This works for both Assay and Assay5 objects in modern SeuratObject
  if (verbose) {
    message("Getting layer/slot '", layer_or_slot, "' from ", class(object)[1])
  }
  return(tryCatch(
    expr = GetAssayData(object = object, layer = layer_or_slot),
    error = function(e) {
      if (verbose) {
        message("Error getting '", layer_or_slot, "': ", conditionMessage(e))
      }
      return(NULL)
    }
  ))
}

#' Set assay data with V4/V5 compatibility
#'
#' Sets data in an Assay or Assay5 object, automatically using the appropriate
#' parameter name (layer for V5, slot for V4).
#'
#' @param object Assay or Assay5 object
#' @param layer_or_slot Name of layer (V5) or slot (V4) to set
#' @param new.data New data to set
#' @param verbose Print diagnostic messages
#'
#' @return Modified Assay or Assay5 object
#' @keywords internal
#'
SetAssayDataCompat <- function(object, layer_or_slot, new.data, verbose = FALSE) {
  # In SeuratObject 5.0+, the 'slot' argument is defunct - always use 'layer'
  # This works for both Assay and Assay5 objects in modern SeuratObject
  if (verbose) {
    message("Setting layer/slot '", layer_or_slot, "' in ", class(object)[1])
  }
  return(SetAssayData(
    object = object,
    layer = layer_or_slot,
    new.data = new.data
  ))
}


