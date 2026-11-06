#' @importFrom hdf5r h5attr guess_dtype H5T_STRING H5T_COMPOUND H5S h5const h5types
#' @importFrom utils txtProgressBar
NULL

# Low-level helpers shared by the h5Seurat / loom / h5ad readers and writers.
# All are internal; they work on hdf5r handles (H5File, H5Group, H5D).

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Object and path inspection
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.require_h5 <- function(x, kinds, what) {
  if (!inherits(x, kinds)) {
    stop(sprintf("'%s' must be one of: %s", what, paste(kinds, collapse = ", ")),
         call. = FALSE)
  }
}

# Does `x` carry attribute `name` holding at least one value? (An attribute
# with a zero-length simple dataspace counts as absent.)
AttrExists <- function(x, name) {
  .require_h5(x, c("H5File", "H5Group", "H5D"), "x")
  if (!isTRUE(x$attr_exists(attr_name = name))) return(FALSE)
  attr <- x$attr_open(attr_name = name)
  on.exit(attr$close(), add = TRUE)
  extent <- attr$get_space()$dims
  length(extent) == 0L || all(extent > 0L)
}

# Like `$exists()`, but safe for nested paths whose intermediate groups may
# be missing (hdf5r errors on those instead of returning FALSE).
Exists <- function(x, name) {
  .require_h5(x, c("H5File", "H5Group"), "x")
  parts <- strsplit(name[1], "/", fixed = TRUE)[[1]]
  partial <- NULL
  for (part in parts[nzchar(parts)]) {
    partial <- if (is.null(partial)) part else paste0(partial, "/", part)
    if (!isTRUE(x$exists(name = partial))) return(FALSE)
  }
  TRUE
}

# Build an HDF5 path from components
H5Path <- function(..., collapse = NULL) {
  paste(..., sep = "/", collapse = collapse)
}

# Is the file behind `x` open for writing? Errors (or warns) when it is not.
Writeable <- function(x, error = TRUE) {
  intent <- as.character(x$get_file_id()$get_intent())
  ok <- !identical(intent, "H5F_ACC_RDONLY")
  if (!ok) {
    msg <- sprintf("%s was opened read-only", x$get_filename())
    if (isTRUE(error)) stop(msg, call. = FALSE)
    warning(msg, call. = FALSE, immediate. = TRUE)
  }
  invisible(ok)
}

# Matrix dimensions of a dense dataset or a sparse matrix group
Dims <- function(x) {
  if (inherits(x, "H5D")) return(x$dims)
  if (!inherits(x, "H5Group") || !IsMatrix(x)) {
    stop("Not a matrix", call. = FALSE)
  }
  if (x$attr_exists(attr_name = "dims")) return(h5attr(x, "dims"))
  # anndata records the (rows, cols) of the CSR/CSC matrix as `shape`
  if (x$attr_exists(attr_name = "shape")) return(rev(h5attr(x, "shape")))
  rep(NA_integer_, 2L)
}

# h5Seurat factor: integer `values` (1-based codes) plus string `levels`
IsFactor <- function(x) {
  if (!inherits(x, "H5Group")) return(FALSE)
  one_d <- function(nm, type) {
    x$exists(name = nm) && inherits(x[[nm]], "H5D") &&
      length(x[[nm]]$dims) == 1L && IsDType(x[[nm]], type)
  }
  if (!one_d("levels", "H5T_STRING") || !one_d("values", "H5T_INTEGER")) {
    return(FALSE)
  }
  codes <- x[["values"]]$read()
  length(unique(codes[!is.na(codes)])) <= x[["levels"]]$dims
}

# Dense 2-D dataset, or a group holding compressed-sparse components
IsMatrix <- function(x) {
  if (inherits(x, "H5D")) {
    # Seurat v5 stores per-layer cell membership as a 2-D logical table
    # under assays/<assay>/cells; it is not an expression matrix
    if (basename(x$get_obj_name()) == "cells" &&
        startsWith(dirname(x$get_obj_name()), "/assays/")) {
      return(FALSE)
    }
    return(length(x$dims) == 2L)
  }
  if (inherits(x, "H5Group")) {
    return(all(vapply(c("data", "indices", "indptr"), function(nm) {
      x$exists(name = nm) && inherits(x[[nm]], "H5D")
    }, logical(1))))
  }
  FALSE
}

# Does dataset `x` have (one of) the HDF5 type classes in `dtype`?
IsDType <- function(x, dtype) {
  if (!inherits(x, "H5D")) {
    stop("IsDType() needs an HDF5 dataset", call. = FALSE)
  }
  type_names <- grep("^H5T_", names(h5types), value = TRUE)
  known <- unique(c(vapply(type_names, function(n) class(h5types[[n]])[1], character(1),
                           USE.NAMES = FALSE),
                    "H5T_COMPOUND"))
  bad <- setdiff(dtype, known)
  if (length(bad) == length(dtype)) {
    stop("None of the requested dtypes are HDF5 datatype classes", call. = FALSE)
  }
  if (length(bad)) {
    warning("Ignoring unknown HDF5 datatype classes: ", paste(bad, collapse = ", "),
            call. = FALSE, immediate. = TRUE)
  }
  inherits(x$get_type(), setdiff(dtype, bad))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Datatypes and encodings
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# h5Seurat stores logicals as integers: FALSE = 0, TRUE = 1, NA = 2
BoolToInt <- function(x) {
  codes <- as.integer(x)
  codes[is.na(codes)] <- 2L
  codes
}

# String datatype: variable-length UTF-8, or the fixed 7-byte ASCII used by
# loom 0.1 files
StringType <- function(stype = c("utf8", "ascii7")) {
  stype <- match.arg(stype)
  if (stype == "ascii7") return(H5T_STRING$new(size = 7L))
  st <- H5T_STRING$new(size = Inf)
  st$set_cset(cset = h5const$H5T_CSET_UTF8)
  st$set_strpad(strpad = h5const$H5T_STR_NULLTERM)
  st
}

# hdf5r's guess_dtype(), with every string (also inside compound types)
# switched to `stype` and logicals mapped to integers when
# options(scConvert.dtypes.logical_to_int) is on (the default).
GuessDType <- function(x, stype = "utf8", ...) {
  guessed <- guess_dtype(x = x, ...)
  if (inherits(guessed, "H5T_STRING")) return(StringType(stype = stype))
  if (inherits(guessed, "H5T_LOGICAL")) {
    if (isTRUE(getOption("scConvert.dtypes.logical_to_int", default = TRUE))) {
      return(guess_dtype(x = BoolToInt(x), ...))
    }
    return(guessed)
  }
  if (inherits(guessed, "H5T_COMPOUND")) {
    members <- lapply(guessed$get_cpd_types(), function(t) {
      if (inherits(t, "H5T_STRING")) StringType(stype = stype) else t
    })
    return(H5T_COMPOUND$new(labels = guessed$get_cpd_labels(), dtypes = members,
                            size = guessed$get_size()))
  }
  guessed
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Chunking, versions, names
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Start/end rows (1-based, inclusive) splitting `dsize` items into blocks of
# `csize`; a single block when `csize` is NA.
ChunkPoints <- function(dsize, csize) {
  if (dsize < 1) {
    starts <- ends <- numeric(0)
  } else if (is.na(csize)) {
    starts <- 1
    ends <- dsize
  } else {
    starts <- seq.int(1, dsize, by = csize)
    ends <- pmin(starts + csize - 1, dsize)
  }
  cbind(start = as.numeric(starts), end = as.numeric(ends))
}

# Which axis to iterate over when streaming a dataset. MARGIN is an axis
# number or one of "largest", "smallest", "first", "last".
GetMargin <- function(dims, MARGIN = getOption("scConvert.chunking.MARGIN")) {
  rules <- c("largest", "smallest", "first", "last")
  if (!is.numeric(MARGIN)) {
    rule <- MARGIN[1]
    if (is.null(rule) || !rule %in% rules) {
      warning("Unknown chunking MARGIN '", rule, "', using 'largest'",
              call. = FALSE, immediate. = TRUE)
      rule <- "largest"
    }
    MARGIN <- switch(rule, largest = which.max(dims), smallest = which.min(dims),
                     first = 1L, last = length(dims))
  }
  if (!MARGIN %in% seq_along(dims)) {
    stop("'MARGIN' must be within the dimensions of the dataset", call. = FALSE)
  }
  MARGIN
}

# Pick the supported version matching `query`: with direction "min" the
# highest target not above the query (rounding down), with "max" the lowest
# target not below it (rounding up). `inclusive = FALSE` excludes equality.
ClosestVersion <- function(query, targets, direction = c("min", "max"),
                           inclusive = direction == "min") {
  direction <- match.arg(direction)
  q <- numeric_version(query)
  tv <- sort(numeric_version(targets))
  ok <- if (direction == "min") {
    if (inclusive) tv <= q else tv < q
  } else {
    if (inclusive) tv >= q else tv > q
  }
  if (!any(ok)) {
    stop(sprintf("No supported version is %s %s",
                 if (direction == "min") "at or below" else "at or above",
                 as.character(q)), call. = FALSE)
  }
  pick <- if (direction == "min") max(which(ok)) else min(which(ok))
  as.character(tv[pick])
}

# Seurat feature names may not contain underscores and must be unique
FixFeatures <- function(features) {
  if (any(grepl("_", features, fixed = TRUE))) {
    message("Feature names cannot have underscores ('_'), replacing with dashes ('-')")
    features <- gsub("_", "-", features, fixed = TRUE)
  }
  if (anyDuplicated(features)) {
    message("Non-unique features (rownames) present, making unique")
    features <- make.unique(features)
  }
  features
}

# Class name as stored in h5Seurat `s4class` attributes: "pkg:Class" unless
# the class comes from one of `packages`
GetClass <- function(class, packages = "Seurat") {
  pkg <- methods::getClass(Class = class[1])@package
  if (pkg %in% packages) class[1] else paste0(pkg, ":", class[1])
}

# Empty matrices: no rows/columns at all, or Seurat's 1x1 NA placeholder
IsMatrixEmpty <- function(x) {
  d <- dim(x)
  if (all(d == 0)) return(TRUE)
  all(d == 1) && all(is.na(x))
}

# Name unnamed elements of a list/vector "<prefix><position>"
PadNames <- function(x, prefix = "index") {
  if (!length(x)) return(x)
  nm <- names(x)
  if (is.null(nm)) nm <- character(length(x))
  blank <- !nzchar(nm) | is.na(nm)
  nm[blank] <- paste0(prefix, which(blank))
  names(x) <- nm
  x
}

MakeSpace <- function(n) strrep(" ", n)

PB <- function() txtProgressBar(style = 3, char = "=", file = stderr())

RandomName <- function(length = 5L, ...) {
  paste0(sample(letters, size = length, ...), collapse = "")
}

# Seurat keys are alphanumeric followed by a single underscore
UpdateKey <- function(key) {
  if (grepl("^[[:alnum:]]+_$", key)) return(key)
  core <- gsub("[^[:alnum:]]", "", key)
  if (!nzchar(core)) core <- RandomName(length = 3L)
  paste0(core, "_")
}

# Column index (1-based) of every stored entry of a CSC matrix, from its
# column pointer
PointerToIndex <- function(p) {
  rep.int(seq_len(length(p) - 1L), diff(p))
}
