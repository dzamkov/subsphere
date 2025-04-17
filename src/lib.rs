#![doc = include_str!("../README.md")]
#![deny(missing_debug_implementations)]
#![deny(missing_docs)]
mod math;

pub mod basetri;
pub mod hex;
pub mod proj;
pub mod tri;
pub mod util;

pub use basetri::BaseTriSphere;
pub use hex::HexSphere;
pub use tri::TriSphere;

use std::num::NonZero;

/// Constructs a tessellation of the unit sphere by projecting an icosahedron onto it.
///
/// The tessellation can be refined by calling methods such as [`TriSphere::subdivide_edge`] or
/// [`TriSphere::truncate`].
pub fn icosphere() -> TriSphere<proj::Fuller> {
    TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(1).unwrap(),
        0,
    )
}

/// Constructs a tessellation of the unit sphere by projecting an octohedron onto it.
///
/// The tessellation can be refined by calling methods such as [`TriSphere::subdivide_edge`] or
/// [`TriSphere::truncate`].
pub fn octosphere() -> TriSphere<proj::Fuller> {
    TriSphere::new(
        BaseTriSphere::Octo,
        proj::Fuller,
        NonZero::new(1).unwrap(),
        0,
    )
}

/// Partitions the surface of the unit sphere into a set of spherical polygons ([`Face`]s).
/// 
/// There are numerous requirements for a valid [`Sphere`] implementation. Custom implementations
/// may use [`util::validate`] to check that these requirements are met.
pub trait Sphere {
    /// The type of [`Face`] on this sphere.
    type Face: Face<Vertex = Self::Vertex, HalfEdge = Self::HalfEdge>;

    /// The type of [`Vertex`] on this sphere.
    type Vertex: Vertex<Face = Self::Face, HalfEdge = Self::HalfEdge>;

    /// The type of [`HalfEdge`] on this sphere.
    type HalfEdge: HalfEdge<Face = Self::Face, Vertex = Self::Vertex>;

    /// The number of faces on this sphere.
    fn num_faces(&self) -> usize;

    /// Gets the [`Face`] with the given [`index`](Face::index).
    fn face(&self, index: usize) -> Self::Face;

    /// Iterates over the faces of this sphere.
    fn faces(&self) -> impl Iterator<Item = Self::Face>;

    /// The number of vertices on this sphere.
    fn num_vertices(&self) -> usize;

    /// Gets the [`Vertex`] with the given [`index`](Vertex::index).
    fn vertex(&self, index: usize) -> Self::Vertex;

    /// Iterates over the vertices of this sphere.
    fn vertices(&self) -> impl Iterator<Item = Self::Vertex>;
}

/// Represents a "face" on a [`Sphere`] of a certain type.
///
/// These are spherical polygons bounded by [`HalfEdge`]s.
pub trait Face: Clone + Eq {
    /// The type of [`Vertex`] on the sphere.
    type Vertex: Vertex<Face = Self, HalfEdge = Self::HalfEdge>;

    /// The type of [`HalfEdge`] on the sphere.
    type HalfEdge: HalfEdge<Face = Self, Vertex = Self::Vertex>;

    /// The index of this vertex within the [`Sphere::faces`] list.
    fn index(&self) -> usize;

    /// The area of this face.
    ///
    /// This is also known as the [solid angle](https://en.wikipedia.org/wiki/Solid_angle)
    /// subtended by the face. The sum of the areas of all faces on a sphere is `4 π`.
    /// 
    /// For faces whose edges are all geodesics, this is equivalent to the [`util::poly_area`]
    /// of the vertices of the face.
    fn area(&self) -> f64;

    /// The number of sides (vertices or edges) that this face has.
    fn num_sides(&self) -> usize;

    /// Gets a [`Vertex`] of this face given its index within [`Face::vertices`].
    ///
    /// This should be equivalent to `side(index).start()`.
    fn vertex(&self, index: usize) -> Self::Vertex {
        self.side(index).start()
    }

    /// Iterates over the [`Vertex`]s of this face in counter-clockwise order.
    ///
    /// Vertices will be returned in an order consistent with [`sides`](Face::sides) and
    /// [`start`](HalfEdge::start). This will iterate over [`num_sides`](Face::num_sides) vertices.
    fn vertices(&self) -> impl Iterator<Item = Self::Vertex> {
        self.sides().map(|h| h.start())
    }

    /// Gets the [`HalfEdge`] which has the given [`index`](HalfEdge::side_index) and this face as
    /// its [`inside`](HalfEdge::inside).
    fn side(&self, index: usize) -> Self::HalfEdge;

    /// Iterates over the [`HalfEdge`]s which have this face as their [`inside`](HalfEdge::inside).
    ///
    /// Edges will be returned in counter-clockwise order around the face, consistent with
    /// [`HalfEdge::side_index`]. The iterator will return [`num_sides`](Face::num_sides) edges.
    fn sides(&self) -> impl Iterator<Item = Self::HalfEdge> {
        (0..self.num_sides()).map(|i| self.side(i))
    }
}

/// Represents a vertex on a [`Sphere`].
pub trait Vertex: Clone + Eq {
    /// The type of [`Face`] on the sphere.
    type Face: Face<Vertex = Self, HalfEdge = Self::HalfEdge>;

    /// The type of [`HalfEdge`] on the sphere.
    type HalfEdge: HalfEdge<Face = Self::Face, Vertex = Self>;

    /// The index of this vertex within the [`Sphere::vertices`] list.
    fn index(&self) -> usize;

    /// The position (or equivalently, the normal) of this vertex.
    fn pos(&self) -> [f64; 3];

    /// The number of edges (or equivalently, [`Face`]s) that are connected to this vertex.
    fn degree(&self) -> usize;

    /// Gets a connected [`Face`] based on its index within the [`Vertex::faces`] list.
    fn face(&self, index: usize) -> Self::Face {
        self.outgoing(index).inside()
    }

    /// Iterates over the [`Face`]s which are connected to this vertex.
    fn faces(&self) -> impl Iterator<Item = Self::Face> {
        self.outgoings().map(|h| h.inside())
    }

    /// Gets an outgoing [`HalfEdge`] based on its index within the [`Vertex::outgoings`] list.
    fn outgoing(&self, index: usize) -> Self::HalfEdge;

    /// Iterates over the outgoing [`HalfEdge`]s which have this vertex as their
    /// [`start`](HalfEdge::start).
    ///
    /// Edges will be returned in counter-clockwise order around the vertex.
    fn outgoings(&self) -> impl Iterator<Item = Self::HalfEdge> {
        (0..self.degree()).map(|i| self.outgoing(i))
    }
}

/// Represents one "side" or direction of an edge on a [`Sphere`].
///
/// Half-edges are oriented such that they go counter-clockwise around their
/// [`inside`](HalfEdge::inside) face.
pub trait HalfEdge: Clone + Eq {
    /// The type of [`Face`] on the sphere.
    type Face: Face<Vertex = Self::Vertex, HalfEdge = Self>;

    /// The type of [`Vertex`] on the sphere.
    type Vertex: Vertex<Face = Self::Face, HalfEdge = Self>;

    /// The index of this half-edge within the [`Face::sides`] list of its
    /// [`inside`](HalfEdge::inside).
    fn side_index(&self) -> usize {
        self.inside().sides().position(|h| h == *self).unwrap()
    }

    /// The index of this half-edge within the [`Vertex::outgoings`] list of its
    /// [`start`](HalfEdge::start).
    fn outgoing_index(&self) -> usize {
        self.start().outgoings().position(|h| h == *self).unwrap()
    }

    /// The length of the edge.
    /// 
    /// For edges that are geodesics, this is equivalent to the [`util::dist`] between its
    /// endpoints.
    fn length(&self) -> f64;

    /// The interior angle between the [previous half-edge](HalfEdge::prev) and this half-edge.
    /// 
    /// The sum of the angles of all [outgoing](Vertex::outgoings) half-edges at a vertex should
    /// be `2 π`.
    /// 
    /// If both this edge, and the previous edge, are geodesics, the angle is equivalent to
    /// [`util::angle`] where the `b` vertex is given by `self.start()`.
    fn angle(&self) -> f64;

    /// Gets the [`Face`] whose interior boundary contains this half-edge.
    fn inside(&self) -> Self::Face;

    /// Gets the [`Vertex`] at the "start" of this half-edge.
    fn start(&self) -> Self::Vertex;

    /// Gets the complementary half-edge on the opposite side of the edge.
    ///
    /// The returned half-edge will go in the opposite direction along the same edge.
    fn complement(&self) -> Self;

    /// Gets the half-edge which shares the [`inside`](HalfEdge::inside) face of this half-edge and
    /// precedes it in counter-clockwise order around the face.
    fn prev(&self) -> Self;

    /// Gets the half-edge which shares the [`inside`](HalfEdge::inside) face of this half-edge and
    /// follows it in counter-clockwise order around the face.
    fn next(&self) -> Self;
}
