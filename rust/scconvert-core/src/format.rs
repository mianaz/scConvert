//! Format identifiers.

use serde::{Deserialize, Serialize};
use std::path::Path;

/// Every format scConvert can read or write. Keep variants in lockstep
/// with the R dispatcher so `FileType()` in `R/zzz.R` and this enum agree.
#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Format {
    H5ad,
    H5seurat,
    H5mu,
    Loom,
    Zarr,
    Rds,
    Sce,
    Soma,
    SpatialData,
    /// Format could not be determined from the path.
    Unknown,
}

impl Format {
    pub fn as_str(&self) -> &'static str {
        match self {
            Format::H5ad => "h5ad",
            Format::H5seurat => "h5seurat",
            Format::H5mu => "h5mu",
            Format::Loom => "loom",
            Format::Zarr => "zarr",
            Format::Rds => "rds",
            Format::Sce => "sce",
            Format::Soma => "soma",
            Format::SpatialData => "spatialdata.zarr",
            Format::Unknown => "unknown",
        }
    }
}

/// Mirror of `R/zzz.R::FileType()`. Keep extension logic identical or
/// dispatch will disagree between R and the native CLI.
pub fn detect_format(path: &Path) -> Format {
    let s = path.to_string_lossy().to_lowercase();
    // URI schemes like `soma://...`
    if let Some(scheme_end) = s.find("://") {
        let scheme = &s[..scheme_end];
        return match scheme {
            "soma" => Format::Soma,
            _ => Format::Unknown,
        };
    }
    // Multi-part extension first.
    if s.ends_with(".spatialdata.zarr") || s.ends_with(".spatialdata.zarr/") {
        return Format::SpatialData;
    }
    if s.ends_with(".h5ad") {
        return Format::H5ad;
    }
    if s.ends_with(".h5seurat") || s.ends_with(".h5Seurat") {
        return Format::H5seurat;
    }
    if s.ends_with(".h5mu") {
        return Format::H5mu;
    }
    if s.ends_with(".loom") {
        return Format::Loom;
    }
    if s.ends_with(".zarr") || s.ends_with(".zarr/") {
        return Format::Zarr;
    }
    if s.ends_with(".rds") {
        return Format::Rds;
    }
    if s.ends_with(".sce") {
        return Format::Sce;
    }
    Format::Unknown
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_common_extensions() {
        assert_eq!(detect_format(Path::new("a.h5ad")), Format::H5ad);
        assert_eq!(detect_format(Path::new("a.h5seurat")), Format::H5seurat);
        assert_eq!(detect_format(Path::new("A.h5Seurat")), Format::H5seurat);
        assert_eq!(detect_format(Path::new("a.h5mu")), Format::H5mu);
        assert_eq!(detect_format(Path::new("a.loom")), Format::Loom);
        assert_eq!(detect_format(Path::new("a.zarr")), Format::Zarr);
        assert_eq!(detect_format(Path::new("a.rds")), Format::Rds);
    }

    #[test]
    fn multi_part_extension_beats_plain_zarr() {
        assert_eq!(
            detect_format(Path::new("study.spatialdata.zarr")),
            Format::SpatialData
        );
    }

    #[test]
    fn soma_uri_scheme() {
        assert_eq!(
            detect_format(Path::new("soma://collection/experiment")),
            Format::Soma
        );
    }

    #[test]
    fn unknown_for_bare_names() {
        assert_eq!(detect_format(Path::new("data")), Format::Unknown);
    }
}
