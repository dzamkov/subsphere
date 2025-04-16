//! Contains utility functions related to spheres, tessellations and projections.
use crate::{Face, HalfEdge, Sphere, Vertex};

/// Validates the internal consistency of the given [`Sphere`].
/// 
/// This will panic if any of the requirements for a sphere are not met. Since this is a violation
/// of the invariants of the [`Sphere`] trait, this is not a recoverable error.
pub fn validate<S: Sphere>(sphere: S)
where
    S::Face: std::fmt::Debug,
    S::Vertex: std::fmt::Debug,
    S::HalfEdge: std::fmt::Debug,
{
    // Validate faces
    let mut index = 0;
    for f in sphere.faces() {
        assert_eq!(f.index(), index, "face index mismatch");
        index += 1;

        // Validate sides
        let mut side_index = 0;
        let mut prev = None;
        for s in f.sides() {
            assert_eq!(s.inside(), f, "side inside mismatch");
            assert_eq!(s.start(), f.vertex(side_index), "side start mismatch");
            assert_eq!(s.side_index(), side_index, "side index mismatch");
            side_index += 1;
            if let Some(prev) = prev {
                assert_eq!(s.prev(), prev, "side prev mismatch");
                assert_eq!(prev.next(), s, "side next mismatch");
            }
            prev = Some(s);
        }
        if let Some(prev) = prev {
            let s = f.side(0);
            assert_eq!(s.prev(), prev, "side prev mismatch");
            assert_eq!(prev.next(), s, "side next mismatch");
        }
        assert_eq!(side_index, f.num_sides(), "side count mismatch");
    }
    assert_eq!(index, sphere.num_faces(), "face count mismatch");

    // Validate vertices
    let mut index = 0;
    for v in sphere.vertices() {
        assert_eq!(v.index(), index, "vertex index mismatch");
        index += 1;

        // Validate outgoing edges
        let mut edge_index = 0;
        let mut prev = None;
        for e in v.outgoings() {
            assert_eq!(
                e.start(),
                v,
                "outgoing edge with index {} does not start at vertex",
                edge_index
            );
            assert_eq!(e.outgoing_index(), edge_index, "edge index mismatch");
            edge_index += 1;
            if let Some(prev) = prev {
                assert_eq!(e.complement().next(), prev, "outgoing prev mismatch");
                assert_eq!(prev.prev().complement(), e, "outgoing next mismatch");
            }
            prev = Some(e);
        }
        if let Some(prev) = prev {
            let e = v.outgoing(0);
            assert_eq!(e.complement().next(), prev, "outgoing prev mismatch");
            assert_eq!(prev.prev().complement(), e, "outgoing next mismatch");
        }
        assert_eq!(edge_index, v.degree(), "vertex degree mismatch");
    }
    assert_eq!(index, sphere.num_vertices(), "vertex count mismatch");

    // Validate total face area
    let total_area: f64 = sphere.faces().map(|f| f.area()).sum();
    assert!(
        (total_area - 4.0 * std::f64::consts::PI).abs() < 1e-12,
        "total area mismatch: expected {:?}, got {:?}",
        4.0 * std::f64::consts::PI,
        total_area
    );
}