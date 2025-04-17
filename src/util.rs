//! Contains utility functions related to spheres, tessellations and projections.
use crate::{Face, HalfEdge, Sphere, Vertex};
use crate::math::{vec, mat};

/// Computes the signed area of a spherical triangle with the given vertices.
/// 
/// This will be positive if the points are in counter-clockwise order, negative if they
/// are in clockwise order, and zero if they are collinear. If any of the points are not on
/// the unit sphere, the result will be undefined.
pub fn tri_area(points: [[f64; 3]; 3]) -> f64 {
    // https://www.johndcook.com/blog/2021/11/29/area-of-spherical-triangle/
    let d = 1.0
        + vec::dot(points[0], points[1])
        + vec::dot(points[1], points[2])
        + vec::dot(points[2], points[0]);
    (mat::det(points) / d).atan() * 2.0
}

/// Computes the signed area of a spherical polygon with the given vertices.
/// 
/// This will be positive if the points are in counter-clockwise order or negative if they
/// are in clockwise order. If any of the points are not on the unit sphere, the result will be
/// undefined.
pub fn poly_area(mut points: impl Iterator<Item = [f64; 3]>) -> f64 {
    let Some(a) = points.next() else {
        return 0.0;
    };
    let Some(mut b) = points.next() else {
        return 0.0;
    };
    let mut res = 0.0;
    for c in points {
        res += tri_area([a, b, c]);
        b = c;
    }
    res
}

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