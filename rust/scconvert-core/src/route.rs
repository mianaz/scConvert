//! Route planner schema.
//!
//! `scconvert plan src dst` emits a [`Route`]. B1 defines the schema; B5
//! implements the planner that fills it in.

use crate::Format;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Route {
    pub src: Format,
    pub dst: Format,
    pub kind: RouteKind,
    pub estimate: RouteEstimate,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum RouteKind {
    /// HDF5↔HDF5 or Zarr↔Zarr with no materialisation — chunks copied.
    Direct { description: String },
    /// Route goes via an in-memory Seurat/AnnData hub. Memory is bounded
    /// by the input size class.
    Hub { via: String, description: String },
    /// Route requires an external tool (e.g. tiledbsoma Arrow writer).
    ObjectMediated { via: String, description: String },
    /// Planner refuses — reason is human-readable.
    Rejected { reason: String },
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RouteEstimate {
    /// Order-of-magnitude peak RAM in bytes. None = unknown.
    pub peak_memory_bytes: Option<u64>,
    pub runtime_class: Option<RuntimeClass>,
    /// Optional package deps required for this route (R/Python packages,
    /// system libraries). Populated for ObjectMediated.
    pub required_dependencies: Vec<String>,
    /// True when every intermediate can be streamed without holding a
    /// full matrix in memory.
    pub streams: bool,
    /// Fields expected to change fate under this route. Lets users
    /// preview the fidelity report without running the conversion.
    pub expected_losses: Vec<String>,
}

#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum RuntimeClass {
    /// Sub-second on test fixtures.
    Instant,
    /// Seconds on 100K-cell data.
    Fast,
    /// Minutes on 1M-cell data.
    Moderate,
    /// Multi-minute even on moderate data; deserves user warning.
    Slow,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejected_route_roundtrips() {
        let r = Route {
            src: Format::H5ad,
            dst: Format::Soma,
            kind: RouteKind::Rejected {
                reason: "SOMA write requires Python tiledbsoma".into(),
            },
            estimate: RouteEstimate::default(),
        };
        let j = serde_json::to_string_pretty(&r).unwrap();
        let back: Route = serde_json::from_str(&j).unwrap();
        assert_eq!(back.src, Format::H5ad);
        match back.kind {
            RouteKind::Rejected { reason } => assert!(reason.contains("SOMA")),
            _ => panic!("wrong variant"),
        }
    }

    #[test]
    fn direct_route_has_streams_flag() {
        let r = Route {
            src: Format::H5ad,
            dst: Format::H5seurat,
            kind: RouteKind::Direct {
                description: "CSR->CSC shape swap".into(),
            },
            estimate: RouteEstimate {
                streams: true,
                runtime_class: Some(RuntimeClass::Fast),
                ..Default::default()
            },
        };
        let j = serde_json::to_string(&r).unwrap();
        assert!(j.contains("\"streams\":true"));
        let back: Route = serde_json::from_str(&j).unwrap();
        assert!(back.estimate.streams);
    }
}
