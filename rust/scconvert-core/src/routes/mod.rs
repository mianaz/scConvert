//! Format-to-format conversion routes.
//!
//! Each `pub` function in this module implements one `(src, dst)` pair
//! end-to-end. `crate::convert()` dispatches on the detected formats and
//! calls into here.
//!
//! B2 scope: `h5ad → h5seurat` for the primary sparse `X` matrix. Other
//! fields (obs, var, obsm, raw, uns) land in follow-up commits within
//! this module.

pub mod h5ad_to_h5seurat;
