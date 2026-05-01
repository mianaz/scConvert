//! Error taxonomy for the scConvert core.
//!
//! Errors that cross the C ABI are serialised via [`ConvertError::code`]
//! plus a JSON payload; callers never see Rust panics.

use thiserror::Error;

pub type Result<T> = std::result::Result<T, ConvertError>;

#[derive(Debug, Error)]
pub enum ConvertError {
    #[error("not yet implemented: {kind}")]
    NotYetImplemented { kind: String },

    #[error("source not found: {path}")]
    SourceMissing { path: String },

    #[error("destination already exists and overwrite is not set: {path}")]
    DestinationExists { path: String },

    #[error("unsupported format pair: {src:?} -> {dst:?}")]
    UnsupportedRoute {
        src: crate::Format,
        dst: crate::Format,
    },

    #[error("malformed input file: {reason}")]
    MalformedInput { reason: String },

    #[error("io error: {0}")]
    Io(#[from] std::io::Error),

    #[error("serialisation error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("memory budget exceeded: needed {needed} bytes, cap is {cap}")]
    MemoryBudgetExceeded { needed: u64, cap: u64 },

    #[error("{0}")]
    Other(String),
}

impl ConvertError {
    /// Stable numeric code used at the C ABI boundary. Do not renumber
    /// after 0.1 ships — bindings rely on these values.
    pub fn code(&self) -> i32 {
        match self {
            ConvertError::NotYetImplemented { .. } => 1,
            ConvertError::SourceMissing { .. } => 2,
            ConvertError::DestinationExists { .. } => 3,
            ConvertError::UnsupportedRoute { .. } => 4,
            ConvertError::MalformedInput { .. } => 5,
            ConvertError::Io(_) => 6,
            ConvertError::Json(_) => 7,
            ConvertError::MemoryBudgetExceeded { .. } => 8,
            ConvertError::Other(_) => 99,
        }
    }
}
