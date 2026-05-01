//! Format-agnostic data model primitives.
//!
//! These types describe *structure*, not *storage*. Concrete IO (B2+)
//! operates on the same vocabulary — `SparseLayout`, `DType`, shape — so
//! converters can reason about metadata without materialising arrays.

use serde::{Deserialize, Serialize};

/// Element types the core understands. Extend by adding variants; never
/// renumber.
#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum DType {
    F32,
    F64,
    I8,
    I16,
    I32,
    I64,
    U8,
    U16,
    U32,
    U64,
    Bool,
    /// Variable-length UTF-8 string.
    Str,
    /// Categorical (codes + levels). See [`CategoricalMeta`].
    Categorical,
}

/// Sparse layout of a 2D matrix. COO is permitted on ingest (some h5ad
/// writers in the wild emit it) but the core canonicalises to CSR/CSC
/// before emitting anything.
#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum SparseLayout {
    Csr,
    Csc,
    Coo,
    Dense,
}

/// Metadata for a categorical column — carried alongside integer codes
/// across the converter. Levels are UTF-8; codes are 0-based with -1
/// reserved for missing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CategoricalMeta {
    pub levels: Vec<String>,
    pub ordered: bool,
}
