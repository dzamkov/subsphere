#![doc = include_str!("../README.md")]
#![deny(missing_debug_implementations)]
#![deny(missing_docs)]
mod math;
mod statictri;

pub mod ico;
pub mod subtri;

pub use subtri::SubTriSphere;
pub use ico::{IcoSphere, BaseIcoSphere};

/// Partitions the surface of the unit sphere into a set of spherical polygons ([`Face`]s).
pub trait Sphere {
    /// The type of [`Face`] on this sphere.
    type Face: Face<Vertex = Self::Vertex, HalfEdge = Self::HalfEdge>;

    /// The type of [`Vertex`] on this sphere.
    type Vertex: Vertex<Face = Self::Face, HalfEdge = Self::HalfEdge>;

    /// The type of [`HalfEdge`] on this sphere.
    type HalfEdge: HalfEdge<Face = Self::Face, Vertex = Self::Vertex>;

    /// The number of faces on this sphere.
    fn num_faces(&self) -> usize;

    /// Iterates over the faces of this sphere.
    fn faces(&self) -> impl Iterator<Item = Self::Face>;

    /// The number of vertices on this sphere.
    fn num_vertices(&self) -> usize;

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
    fn area(&self) -> f64 {
        todo!()
    }

    /// The number of sides (vertices or edges) that this face has.
    fn num_sides(&self) -> usize;

    /// Gets a [`Vertex`] of this face given its index within [`vertices`](Face::vertices).
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

    /// Gets the [`HalfEdge`] which has the given [`index`](HalfEdge::index) and this face as its
    /// [`inside`](HalfEdge::inside).
    fn side(&self, index: usize) -> Self::HalfEdge;

    /// Iterates over the [`HalfEdge`]s which have this face as their [`inside`](HalfEdge::inside).
    ///
    /// Edges will be returned in counter-clockwise order around the face, consistent with
    /// [`HalfEdge::index`]. The iterator will return [`num_sides`](Face::num_sides) edges.
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

    /// Iterates over the [`Face`]s which are connected to this vertex.
    fn faces(&self) -> impl Iterator<Item = Self::Face> {
        self.outgoing().map(|h| h.inside())
    }

    /// The first [`HalfEdge`] in [`Vertex::outgoing`].
    ///
    /// [`Vertex::outgoing`] must not be empty, so this is often a convenient and performant way
    /// of getting an arbitrary outgoing half-edge.
    fn first_outgoing(&self) -> Self::HalfEdge;

    /// Iterates over the [`HalfEdge`]s which have this vertex as their [`HalfEdge::start`].
    ///
    /// Edges will be returned in counter-clockwise order around the vertex, starting with
    /// [`Vertex::first_outgoing`]. The iterator will return [`Vertex::degree`] edges.
    ///
    /// Each returned [`HalfEdge`] corresponding to a [`Face`] which connects to this vertex. These
    /// faces can be obtained using [`HalfEdge::inside`].
    fn outgoing(&self) -> impl Iterator<Item = Self::HalfEdge> {
        let mut h = self.first_outgoing();
        (0..self.degree()).map(move |_| {
            let res = h.clone();
            h = h.prev().complement();
            res
        })
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
    fn index(&self) -> usize;

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
