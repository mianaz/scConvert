#' @include zzz.R
#' @importFrom Matrix sparseMatrix t
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AnnData Encoding/Decoding Helpers (Backend-Agnostic)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Decode AnnData categorical encoding to R factor
#'
#' @param codes Integer vector of 0-based category codes (-1 = NA)
#' @param categories Character vector of category labels
#'
#' @return A factor vector
#'
#' @keywords internal
#'
DecodeCategorical <- function(codes, categories) {
  codes[codes == -1L] <- NA_integer_
  valid <- !is.na(codes) & codes >= 0L & codes < length(categories)
  if (!all(valid[!is.na(codes)])) {
    n_bad <- sum(!valid[!is.na(codes)])
    warning(sprintf(
      "Categorical decoding: %d code(s) out of bounds [0, %d). Setting to NA.",
      n_bad, length(categories)
    ))
    codes[!is.na(codes) & !valid] <- NA_integer_
  }
  factor(categories[codes + 1L], levels = categories)
}

#' Encode R factor as AnnData categorical
#'
#' @param x A factor (or coercible to factor)
#'
#' @return A list with \code{codes} (0-based integer) and \code{categories}
#'   (character vector of levels)
#'
#' @keywords internal
#'
EncodeCategorical <- function(x) {
  x <- as.factor(x)
  codes <- as.integer(x) - 1L
  codes[is.na(x)] <- -1L
  list(codes = codes, categories = levels(x))
}

#' Reconstruct sparse CSR/CSC matrix from array components
#'
#' @param data Numeric vector of non-zero values
#' @param indices Integer vector of column (CSR) or row (CSC) indices (0-based)
#' @param indptr Integer vector of row (CSR) or column (CSC) pointers
#' @param shape Integer vector of length 2: c(n_rows, n_cols)
#' @param transpose If TRUE, transpose result (h5ad stores cells x genes,
#'   Seurat needs genes x cells)
#'
#' @return A dgCMatrix
#'
#' @keywords internal
#'
ReconstructSparseCSR <- function(data, indices, indptr, shape, transpose = TRUE) {
  indices_1based <- indices + 1L
  n_rows <- shape[1]
  n_cols <- shape[2]
  row_indices <- rep(seq_len(n_rows), diff(indptr))
  mat <- sparseMatrix(
    i = row_indices,
    j = indices_1based,
    x = data,
    dims = c(n_rows, n_cols),
    index1 = TRUE
  )
  if (transpose) mat <- t(mat)
  mat
}

#' Deconstruct dgCMatrix to AnnData CSR components
#'
#' AnnData stores expression matrices as cells x genes in CSR format.
#' This function converts a genes x cells dgCMatrix to CSR components
#' for the cells x genes layout.
#'
#' @param mat A dgCMatrix (genes x cells)
#' @param coerce If TRUE (default), coerce mat to dgCMatrix first
#'
#' @return A list with \code{data}, \code{indices} (0-based), \code{indptr},
#'   and \code{shape} (cells, genes)
#'
#' @keywords internal
#'
DeconstructSparseCSR <- function(mat, coerce = TRUE) {
  if (coerce) mat <- as(mat, "dgCMatrix")
  # CSC of (genes x cells) == CSR of (cells x genes) — zero-copy reinterpretation
  # dgCMatrix @p = column pointers, @i = row indices (0-based), @x = values
  # Reinterpreted as CSR: @p = row pointers, @i = column indices, @x = values
  list(
    data = mat@x,
    indices = mat@i,       # 0-based row indices in CSC = column indices in transposed CSR
    indptr = mat@p,        # column pointers in CSC = row pointers in transposed CSR
    shape = as.integer(rev(dim(mat)))  # c(n_cells, n_genes)
  )
}

#' Map AnnData reduction name to Seurat key
#'
#' @param name Reduction name (e.g., "pca", "umap", "tsne")
#'
#' @return A Seurat-style key string (e.g., "PC_", "UMAP_")
#'
#' @keywords internal
#'
AnnDataReductionKey <- function(name) {
  switch(name,
    "pca" = "PC_",
    "tsne" = "tSNE_",
    "umap" = "UMAP_",
    "harmony" = "harmony_",
    "scvi" = "scVI_",
    paste0(toupper(name), "_")
  )
}

#' Map AnnData layer name to Seurat layer name
#'
#' @param name AnnData layer name
#'
#' @return Seurat layer name
#'
#' @keywords internal
#'
AnnDataLayerToSeurat <- function(name) {
  switch(name,
    "counts" = "counts",
    "data" = "data",
    "log_normalized" = "data",
    "scale.data" = "scale.data",
    "scaled" = "scale.data",
    name
  )
}

#' Map Seurat layer name to AnnData layer name
#'
#' @param name Seurat layer name
#'
#' @return AnnData layer name
#'
#' @keywords internal
#'
SeuratLayerToAnnData <- function(name) {
  switch(name,
    "counts" = "counts",
    "data" = "data",
    "scale.data" = "scaled",
    name
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Zarr v2 Low-Level I/O Helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Detect zarr store format version
#'
#' @param store_path Path to zarr store directory
#'
#' @return 2L or 3L
#'
#' @keywords internal
#'
.zarr_store_version <- function(store_path) {
  if (file.exists(file.path(store_path, ".zgroup")) ||
      file.exists(file.path(store_path, ".zarray"))) return(2L)
  if (file.exists(file.path(store_path, "zarr.json"))) return(3L)
  stop("Not a valid zarr store: ", store_path, call. = FALSE)
}

#' Read and parse a JSON file
#'
#' @param path Path to JSON file
#'
#' @return Parsed list
#'
#' @keywords internal
#'
.zarr_read_json <- function(path) {
  if (!file.exists(path)) return(list())
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

#' Read zarr node attributes
#'
#' Reads .zattrs (v2) or zarr.json attributes field (v3).
#'
#' @param store_path Path to zarr store
#' @param rel_path Relative path within store (empty string for root)
#'
#' @return Named list of attributes
#'
#' @keywords internal
#'
.zarr_read_attrs <- function(store_path, rel_path = "") {
  full_path <- if (nchar(rel_path) > 0) {
    file.path(store_path, rel_path)
  } else {
    store_path
  }
  version <- .zarr_store_version(store_path)

  if (version == 2L) {
    .zarr_read_json(file.path(full_path, ".zattrs"))
  } else {
    meta <- .zarr_read_json(file.path(full_path, "zarr.json"))
    meta$attributes %||% list()
  }
}

#' Determine zarr node type
#'
#' @param store_path Path to zarr store
#' @param rel_path Relative path within store
#'
#' @return "group", "array", or "missing"
#'
#' @keywords internal
#'
.zarr_node_type <- function(store_path, rel_path) {
  full_path <- file.path(store_path, rel_path)
  if (!dir.exists(full_path) && !file.exists(full_path)) return("missing")

  version <- .zarr_store_version(store_path)

  if (version == 2L) {
    if (file.exists(file.path(full_path, ".zarray"))) return("array")
    if (file.exists(file.path(full_path, ".zgroup"))) return("group")
    return("missing")
  }

  # v3
  zj <- file.path(full_path, "zarr.json")
  if (!file.exists(zj)) return("missing")
  meta <- .zarr_read_json(zj)
  meta$node_type %||% "missing"
}

#' List children of a zarr group
#'
#' @param store_path Path to zarr store
#' @param rel_path Relative path to group within store
#'
#' @return Character vector of child names
#'
#' @keywords internal
#'
.zarr_list_children <- function(store_path, rel_path = "") {
  full_path <- if (nchar(rel_path) > 0) {
    file.path(store_path, rel_path)
  } else {
    store_path
  }
  if (!dir.exists(full_path)) return(character(0))
  entries <- list.dirs(full_path, full.names = FALSE, recursive = FALSE)
  # Filter out zarr metadata entries
  entries[!grepl("^\\.|^c$|^__", entries)]
}

#' Decompress chunk data
#'
#' @param raw_data Raw vector of compressed data
#' @param compressor Compressor specification list from .zarray metadata
#'
#' @return Decompressed raw vector
#'
#' @keywords internal
#'
.zarr_decompress <- function(raw_data, compressor) {
  if (is.null(compressor)) return(raw_data)
  comp_id <- compressor$id
  if (comp_id %in% c("zlib", "gzip")) {
    memDecompress(raw_data, type = "gzip")
  } else if (comp_id == "blosc") {
    if (!requireNamespace("blosc", quietly = TRUE)) {
      stop("blosc package required for blosc-compressed zarr data. ",
           "Install with: install.packages('blosc')", call. = FALSE)
    }
    blosc::blosc_decompress(raw_data)
  } else {
    stop("Unsupported zarr compressor: ", comp_id,
         ". Supported: zlib, gzip, blosc.", call. = FALSE)
  }
}

#' Compress chunk data
#'
#' @param raw_data Raw vector to compress
#' @param compressor Compressor specification list
#'
#' @return Compressed raw vector
#'
#' @keywords internal
#'
.zarr_compress <- function(raw_data, compressor) {
  if (is.null(compressor)) return(raw_data)
  comp_id <- compressor$id
  if (comp_id %in% c("zlib", "gzip")) {
    memCompress(raw_data, type = "gzip")
  } else if (comp_id == "blosc") {
    if (!requireNamespace("blosc", quietly = TRUE)) {
      stop("blosc package required for blosc compression", call. = FALSE)
    }
    blosc::blosc_compress(raw_data)
  } else {
    stop("Unsupported compressor: ", comp_id, call. = FALSE)
  }
}

#' Parse zarr v2 dtype string to R type info
#'
#' @param dtype Zarr dtype string (e.g., "<f8", "<i4", "|O")
#'
#' @return List with r_type, size, endian, is_object
#'
#' @keywords internal
#'
.zarr_parse_dtype <- function(dtype) {
  endian <- "little"
  if (grepl("^>", dtype)) endian <- "big"
  dtype_clean <- gsub("^[<>|=]", "", dtype)

  type_char <- substr(dtype_clean, 1, 1)
  size_str <- gsub("[^0-9]", "", dtype_clean)
  size <- if (nchar(size_str) > 0) as.integer(size_str) else 1L

  list(
    r_type = switch(type_char,
      "f" = "double",
      "i" = "integer",
      "u" = "integer",
      "b" = "logical",
      "O" = "character",
      "S" = "character",
      "U" = "character",
      "double"
    ),
    size = size,
    endian = endian,
    signed = type_char != "u",
    is_object = dtype_clean == "O"
  )
}

#' Read a zarr v2 numeric array from disk
#'
#' @param store_path Path to zarr store
#' @param rel_path Relative path to array within store
#'
#' @return Numeric vector or matrix
#'
#' @keywords internal
#'
.zarr_read_numeric <- function(store_path, rel_path) {
  full_path <- file.path(store_path, rel_path)
  meta <- .zarr_read_json(file.path(full_path, ".zarray"))

  dtype_info <- .zarr_parse_dtype(meta$dtype)
  shape <- if (is.list(meta$shape)) unlist(meta$shape) else meta$shape
  chunks <- if (is.list(meta$chunks)) unlist(meta$chunks) else meta$chunks
  compressor <- meta$compressor
  order <- meta$order %||% "C"
  n_elements <- prod(shape)

  # Calculate chunk grid
  ndim <- length(shape)
  n_chunks_per_dim <- ceiling(shape / chunks)
  total_chunks <- prod(n_chunks_per_dim)

  if (total_chunks == 1) {
    # Single chunk - most common for AnnData arrays
    chunk_file <- .zarr_chunk_path_v2(full_path, rep(0L, ndim))
    raw_data <- readBin(chunk_file, "raw", file.info(chunk_file)$size)
    raw_data <- .zarr_decompress(raw_data, compressor)
    values <- .raw_to_r_type(raw_data, dtype_info, prod(chunks))
    # Trim to actual shape if chunks > shape
    values <- values[seq_len(n_elements)]
  } else {
    # Multi-chunk: assemble from chunk grid
    values <- .zarr_read_chunked(full_path, meta, dtype_info, shape, chunks,
                                  n_chunks_per_dim, compressor)
  }

  if (ndim == 1) return(values)
  if (ndim == 2) {
    return(matrix(values, nrow = shape[1], ncol = shape[2], byrow = (order == "C")))
  }
  values
}

#' Build chunk file path for zarr v2
#'
#' @param array_dir Array directory path
#' @param chunk_coords Integer vector of chunk coordinates (0-based)
#'
#' @return File path to chunk
#'
#' @keywords internal
#'
.zarr_chunk_path_v2 <- function(array_dir, chunk_coords) {
  if (length(chunk_coords) == 1) {
    file.path(array_dir, as.character(chunk_coords))
  } else {
    file.path(array_dir, paste(chunk_coords, collapse = "."))
  }
}

#' Read multi-chunk zarr v2 numeric array
#'
#' @keywords internal
#'
.zarr_read_chunked <- function(array_dir, meta, dtype_info, shape, chunks,
                                n_chunks_per_dim, compressor) {
  ndim <- length(shape)
  n_elements <- prod(shape)

  if (ndim == 1) {
    values <- numeric(n_elements)
    for (ci in seq_len(n_chunks_per_dim) - 1L) {
      chunk_file <- .zarr_chunk_path_v2(array_dir, ci)
      if (!file.exists(chunk_file)) next
      raw_data <- readBin(chunk_file, "raw", file.info(chunk_file)$size)
      raw_data <- .zarr_decompress(raw_data, compressor)
      chunk_vals <- .raw_to_r_type(raw_data, dtype_info, chunks)
      start_idx <- ci * chunks + 1L
      end_idx <- min(start_idx + chunks - 1L, n_elements)
      n_copy <- end_idx - start_idx + 1L
      values[start_idx:end_idx] <- chunk_vals[seq_len(n_copy)]
    }
    return(values)
  }

  if (ndim == 2) {
    order <- meta$order %||% "C"
    # For 2D, iterate through chunk grid in C order
    values <- numeric(n_elements)
    for (ci in seq_len(n_chunks_per_dim[1]) - 1L) {
      for (cj in seq_len(n_chunks_per_dim[2]) - 1L) {
        chunk_file <- .zarr_chunk_path_v2(array_dir, c(ci, cj))
        if (!file.exists(chunk_file)) next
        raw_data <- readBin(chunk_file, "raw", file.info(chunk_file)$size)
        raw_data <- .zarr_decompress(raw_data, compressor)
        chunk_size <- prod(chunks)
        chunk_vals <- .raw_to_r_type(raw_data, dtype_info, chunk_size)

        # Copy chunk data into the correct position in the flat array
        row_start <- ci * chunks[1]
        col_start <- cj * chunks[2]
        row_end <- min(row_start + chunks[1] - 1L, shape[1] - 1L)
        col_end <- min(col_start + chunks[2] - 1L, shape[2] - 1L)

        if (order == "C") {
          # C order: row-major in chunk, row-major in output
          chunk_mat <- matrix(chunk_vals, nrow = chunks[1], ncol = chunks[2],
                              byrow = TRUE)
          for (r in row_start:row_end) {
            for (cc in col_start:col_end) {
              flat_idx <- r * shape[2] + cc + 1L
              values[flat_idx] <- chunk_mat[r - row_start + 1L, cc - col_start + 1L]
            }
          }
        } else {
          # F order: column-major
          chunk_mat <- matrix(chunk_vals, nrow = chunks[1], ncol = chunks[2])
          for (cc in col_start:col_end) {
            for (r in row_start:row_end) {
              flat_idx <- cc * shape[1] + r + 1L
              values[flat_idx] <- chunk_mat[r - row_start + 1L, cc - col_start + 1L]
            }
          }
        }
      }
    }
    return(values)
  }

  stop("Arrays with > 2 dimensions not supported", call. = FALSE)
}

#' Convert raw bytes to R type
#'
#' @keywords internal
#'
.raw_to_r_type <- function(raw_data, dtype_info, n) {
  if (dtype_info$r_type == "double") {
    readBin(raw_data, "double", n = n, size = dtype_info$size,
            endian = dtype_info$endian)
  } else if (dtype_info$r_type == "integer") {
    readBin(raw_data, "integer", n = n, size = dtype_info$size,
            signed = dtype_info$signed, endian = dtype_info$endian)
  } else if (dtype_info$r_type == "logical") {
    as.logical(readBin(raw_data, "integer", n = n, size = 1L))
  } else {
    stop("Unsupported zarr dtype for numeric reading", call. = FALSE)
  }
}

#' Read a vlen-utf8 encoded string array from zarr v2 store
#'
#' AnnData zarr stores encode string arrays using the vlen-utf8 filter.
#' Each string in the chunk is stored as a 4-byte little-endian int32
#' length prefix followed by that many bytes of UTF-8 text.
#'
#' @param store_path Path to zarr store
#' @param rel_path Relative path to string array within store
#'
#' @return Character vector
#'
#' @keywords internal
#'
.zarr_read_strings <- function(store_path, rel_path) {
  full_path <- file.path(store_path, rel_path)

  # Read array metadata
  meta <- .zarr_read_json(file.path(full_path, ".zarray"))
  shape <- if (is.list(meta$shape)) unlist(meta$shape) else meta$shape
  total_n <- prod(shape)
  chunks <- if (is.list(meta$chunks)) unlist(meta$chunks) else meta$chunks
  compressor <- meta$compressor

  n_chunks <- ceiling(total_n / chunks[1])
  all_strings <- character(0)

  for (chunk_idx in seq_len(n_chunks) - 1L) {
    chunk_file <- file.path(full_path, as.character(chunk_idx))
    if (!file.exists(chunk_file)) {
      # Fill with empty strings for missing chunks
      all_strings <- c(all_strings, rep("", min(chunks[1], total_n - length(all_strings))))
      next
    }
    raw_data <- readBin(chunk_file, "raw", file.info(chunk_file)$size)
    if (length(raw_data) == 0) next
    raw_data <- .zarr_decompress(raw_data, compressor)
    n_in_chunk <- min(chunks[1], total_n - length(all_strings))
    strings <- .decode_vlen_utf8(raw_data, n_in_chunk)
    all_strings <- c(all_strings, strings)
  }

  all_strings[seq_len(total_n)]
}

#' Decode vlen-utf8 encoded raw data to character vector
#'
#' @param raw_data Raw vector containing vlen-utf8 encoded strings
#' @param n Number of strings to decode
#'
#' @return Character vector of length n
#'
#' @keywords internal
#'
.decode_vlen_utf8 <- function(raw_data, n) {
  con <- rawConnection(raw_data, "r")
  on.exit(close(con))
  # numcodecs VLenUTF8 format starts with a 4-byte item count header.
  # Read the first int32: if it equals n, it's the header; otherwise
  # it's the length of the first string (legacy format without header).
  first_int <- readBin(con, "integer", n = 1L, size = 4L, endian = "little")
  if (length(first_int) == 0 || is.na(first_int)) return(character(n))
  strings <- character(n)
  if (first_int == n) {
    # Standard numcodecs format: first_int is the item count header
    for (i in seq_len(n)) {
      len <- readBin(con, "integer", n = 1L, size = 4L, endian = "little")
      if (length(len) == 0 || is.na(len)) break
      if (len < 0) {
        strings[i] <- NA_character_
      } else if (len == 0) {
        strings[i] <- ""
      } else {
        str_bytes <- readBin(con, "raw", n = len)
        strings[i] <- rawToChar(str_bytes)
      }
    }
  } else {
    # Legacy format without header: first_int is the first string length
    if (first_int < 0) {
      strings[1] <- NA_character_
    } else if (first_int == 0) {
      strings[1] <- ""
    } else {
      str_bytes <- readBin(con, "raw", n = first_int)
      strings[1] <- rawToChar(str_bytes)
    }
    for (i in seq_len(n - 1L) + 1L) {
      len <- readBin(con, "integer", n = 1L, size = 4L, endian = "little")
      if (length(len) == 0 || is.na(len)) break
      if (len < 0) {
        strings[i] <- NA_character_
      } else if (len == 0) {
        strings[i] <- ""
      } else {
        str_bytes <- readBin(con, "raw", n = len)
        strings[i] <- rawToChar(str_bytes)
      }
    }
  }
  strings
}

#' Encode character vector to vlen-utf8 raw data
#'
#' @param strings Character vector to encode
#'
#' @return Raw vector in vlen-utf8 format
#'
#' @keywords internal
#'
.encode_vlen_utf8 <- function(strings) {
  n <- length(strings)
  is_na <- is.na(strings)
  strings[is_na] <- ""
  # Vectorized: convert all strings to raw bytes at once
  byte_list <- lapply(strings, charToRaw)
  lens <- vapply(byte_list, length, integer(1L))
  lens[is_na] <- 0L
  # Pre-allocate: 4 (header) + n*4 (length prefixes) + sum(lens) (string bytes)
  total_size <- 4L + n * 4L + sum(lens)
  buf <- raw(total_size)
  # Write item count header
  buf[1:4] <- writeBin(n, raw(), size = 4L, endian = "little")
  # Pre-compute all length prefixes as raw
  len_raw <- writeBin(lens, raw(), size = 4L, endian = "little")
  # Interleave: [len1][bytes1][len2][bytes2]...
  pos <- 5L
  for (i in seq_len(n)) {
    buf[pos:(pos + 3L)] <- len_raw[((i - 1L) * 4L + 1L):(i * 4L)]
    pos <- pos + 4L
    if (lens[i] > 0L) {
      end_pos <- pos + lens[i] - 1L
      buf[pos:end_pos] <- byte_list[[i]]
      pos <- end_pos + 1L
    }
  }
  buf
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Zarr v2 Write Helpers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a zarr v2 group directory with metadata
#'
#' @param dir Directory path for the group
#' @param attrs Named list of attributes (written to .zattrs)
#'
#' @keywords internal
#'
.zarr_create_group <- function(dir, attrs = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines('{"zarr_format":2}', file.path(dir, ".zgroup"))
  if (!is.null(attrs) && length(attrs) > 0) {
    .zarr_write_attrs(dir, attrs)
  }
}

#' Write zarr attributes to .zattrs file
#'
#' @param dir Directory containing the zarr node
#' @param attrs Named list of attributes
#'
#' @keywords internal
#'
.zarr_write_attrs <- function(dir, attrs) {
  json <- jsonlite::toJSON(attrs, auto_unbox = TRUE, null = "null", pretty = TRUE)
  writeLines(json, file.path(dir, ".zattrs"))
}

#' Write a zarr v2 numeric array
#'
#' @param dir Directory for the array
#' @param data Numeric vector or matrix
#' @param dtype Zarr dtype string (default "<f8" for float64)
#' @param compressor Compressor spec list (default gzip level 4)
#' @param attrs Optional attributes for .zattrs
#'
#' @keywords internal
#'
.zarr_write_numeric <- function(dir, data, dtype = NULL, compressor = NULL,
                                 attrs = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  # Determine dtype from data
  if (is.null(dtype)) {
    dtype <- if (is.integer(data)) "<i4" else "<f8"
  }

  # Default compressor
  if (is.null(compressor)) {
    compressor <- list(id = "zlib", level = 4L)
  }

  # Determine shape and flatten data
  if (is.matrix(data)) {
    shape <- dim(data)
    # Convert to C order (row-major) from R's column-major
    flat_data <- as.vector(t(data))
    order <- "C"
  } else {
    shape <- length(data)
    flat_data <- data
    order <- "C"
  }
  # Single chunk per array — our zarr writer only writes one chunk file.
  # Multi-chunk would require writing separate chunk files (0, 1, 2, ...).
  chunks <- shape

  # Write .zarray metadata
  meta <- list(
    zarr_format = 2L,
    shape = as.list(as.integer(shape)),
    chunks = as.list(as.integer(chunks)),
    dtype = dtype,
    compressor = compressor,
    fill_value = if (grepl("b", dtype)) FALSE else if (grepl("f", dtype)) 0.0 else 0L,
    order = order,
    filters = NULL
  )
  json <- jsonlite::toJSON(meta, auto_unbox = TRUE, null = "null", pretty = TRUE)
  writeLines(json, file.path(dir, ".zarray"))

  # Write attributes
  if (!is.null(attrs) && length(attrs) > 0) {
    .zarr_write_attrs(dir, attrs)
  }

  # Convert data to raw bytes
  dtype_info <- .zarr_parse_dtype(dtype)
  if (dtype_info$r_type == "double") {
    raw_data <- writeBin(as.double(flat_data), raw(), size = dtype_info$size,
                         endian = "little")
  } else if (dtype_info$r_type == "integer") {
    raw_data <- writeBin(as.integer(flat_data), raw(), size = dtype_info$size,
                         endian = "little")
  } else if (dtype_info$r_type == "logical") {
    # Boolean: write as single-byte integers (0/1)
    raw_data <- as.raw(as.integer(flat_data))
  } else {
    stop("Unsupported dtype for writing: ", dtype, call. = FALSE)
  }

  # Compress
  raw_data <- .zarr_compress(raw_data, compressor)

  # Write chunk
  chunk_name <- if (length(shape) == 2) "0.0" else "0"
  writeBin(raw_data, file.path(dir, chunk_name))
}

#' Write a zarr v2 vlen-utf8 string array
#'
#' @param dir Directory for the array
#' @param strings Character vector
#' @param compressor Compressor spec list (default gzip level 4)
#' @param attrs Optional attributes for .zattrs
#'
#' @keywords internal
#'
.zarr_write_strings <- function(dir, strings, compressor = NULL, attrs = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(compressor)) {
    compressor <- list(id = "zlib", level = 4L)
  }

  n <- length(strings)

  # Write .zarray metadata
  meta <- list(
    zarr_format = 2L,
    shape = list(as.integer(n)),
    chunks = list(as.integer(n)),
    dtype = "|O",
    compressor = compressor,
    fill_value = 0L,
    order = "C",
    filters = list(list(id = "vlen-utf8"))
  )
  json <- jsonlite::toJSON(meta, auto_unbox = TRUE, null = "null", pretty = TRUE)
  writeLines(json, file.path(dir, ".zarray"))

  if (!is.null(attrs) && length(attrs) > 0) {
    .zarr_write_attrs(dir, attrs)
  }

  # Encode strings as vlen-utf8
  raw_data <- .encode_vlen_utf8(strings)

  # Compress
  raw_data <- .zarr_compress(raw_data, compressor)

  # Write chunk
  writeBin(raw_data, file.path(dir, "0"))
}
