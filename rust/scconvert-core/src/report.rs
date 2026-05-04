//! Fidelity report schema.
//!
//! Every successful conversion emits one. The JSON form is the stable
//! external contract; the Rust types below must stay serde-roundtrippable
//! after field additions. Additive only — consumers (R, Python, CI) pin
//! on field presence, not field count.

use crate::Format;
use serde::{Deserialize, Serialize};

/// Serialization wire format for the report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportFormat {
    Json,
    JsonPretty,
    // Yaml in a future minor version.
}

/// Forward-compat: deserializers ignore unknown fields so an older
/// binding consuming a newer `scconvert-core`'s report does not error.
/// New fields go in additively; renames/removals are breaking changes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FidelityReport {
    /// Schema version of this report. Bump on breaking changes.
    pub schema_version: String,
    /// scconvert-core version that produced the report.
    pub core_version: String,
    pub src_format: Format,
    pub dst_format: Format,
    pub src_format_version: Option<String>,
    pub dst_format_version: Option<String>,
    /// Per-top-level-field record (X, obs, var, obsm, obsp, uns, images, …).
    pub fields: Vec<FieldRecord>,
    /// Machine-readable warnings not tied to a specific field.
    pub warnings: Vec<String>,
    /// Total bytes written to `dst`, best-effort.
    pub bytes_written: Option<u64>,
    /// Monotonic wall-clock conversion time in milliseconds.
    pub duration_ms: Option<u64>,
}

impl FidelityReport {
    pub const SCHEMA_VERSION: &'static str = "0.1";

    pub fn new(src: Format, dst: Format) -> Self {
        FidelityReport {
            schema_version: Self::SCHEMA_VERSION.to_string(),
            core_version: env!("CARGO_PKG_VERSION").to_string(),
            src_format: src,
            dst_format: dst,
            src_format_version: None,
            dst_format_version: None,
            fields: Vec::new(),
            warnings: Vec::new(),
            bytes_written: None,
            duration_ms: None,
        }
    }

    pub fn to_json(&self, pretty: bool) -> crate::Result<String> {
        if pretty {
            Ok(serde_json::to_string_pretty(self)?)
        } else {
            Ok(serde_json::to_string(self)?)
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FieldRecord {
    /// Canonical field path. `/X`, `/obs/cell_type`, `/obsm/X_pca`, …
    pub path: String,
    pub fate: FieldFate,
    /// Free-form explanation when [`FieldFate`] alone is insufficient
    /// (e.g. which categorical encoding was substituted).
    pub note: Option<String>,
    /// Optional digest used to sanity-check numerical identity across
    /// roundtrips; populated when `ConvertOptions::compute_digests = true`.
    pub digest: Option<ArrayDigest>,
}

/// What happened to a field on its way through the converter.
#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FieldFate {
    /// Bit-exact preservation.
    Preserved,
    /// Shape or encoding changed but information content did not (e.g.
    /// CSR→CSC, factor↔categorical level remap).
    Transformed,
    /// Field is representable in `src` but not `dst`; silently dropped.
    Dropped,
    /// Field could be represented in `dst` but was too malformed in `src`
    /// to migrate safely.
    DroppedMalformed,
    /// Field was kept but required lossy conversion.
    Lossy,
}

/// Lightweight numeric digest for per-array verification.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArrayDigest {
    pub element_count: u64,
    /// XxHash64 of the raw little-endian element bytes. Cheap; stable
    /// across endian-normalised ABIs.
    pub xxh64: u64,
    pub dtype: crate::DType,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_report_roundtrips() {
        let r = FidelityReport::new(Format::H5ad, Format::H5seurat);
        let j = r.to_json(false).unwrap();
        let back: FidelityReport = serde_json::from_str(&j).unwrap();
        assert_eq!(r.src_format, back.src_format);
        assert_eq!(r.dst_format, back.dst_format);
        assert_eq!(r.schema_version, back.schema_version);
    }

    #[test]
    fn field_record_roundtrips() {
        let mut r = FidelityReport::new(Format::H5mu, Format::H5seurat);
        r.fields.push(FieldRecord {
            path: "/mod/rna/X".into(),
            fate: FieldFate::Transformed,
            note: Some("CSR -> CSC (zero-copy, shape swap)".into()),
            digest: None,
        });
        let j = r.to_json(true).unwrap();
        let back: FidelityReport = serde_json::from_str(&j).unwrap();
        assert_eq!(back.fields.len(), 1);
        assert_eq!(back.fields[0].fate, FieldFate::Transformed);
    }
}
