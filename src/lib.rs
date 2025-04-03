pub mod icosphere;
mod math;

pub use icosphere::IcoSphere;

/// Partitions the surface of a sphere into a set of spherical polygons ([`Face`]s).
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
    fn area(&self) -> f64;

    /// The number of sides (vertices or edges) that this face has.
    fn num_sides(&self) -> usize;

    /// Iterates over the [`Vertex`]s of this face in counter-clockwise order, starting at
    /// an arbitrary vertex.
    ///
    /// This will iterate over [`Face::num_sides`] vertices.
    fn vertices(&self) -> impl Iterator<Item = Self::Vertex> {
        self.boundary().map(|h| h.start())
    }

    /// The first [`HalfEdge`] in [`Face::boundary`].
    ///
    /// [`Face::boundary`] must not be empty, so this is often a convenient and performant way
    /// of getting an arbitrary boundary half-edge.
    fn first_boundary(&self) -> Self::HalfEdge;

    /// Iterates over the [`HalfEdge`]s which have this face as their [`HalfEdge::inside`].
    /// 
    /// Edges will be returned in counter-clockwise order around the face, starting with
    /// [`Face::first_boundary`]. The iterator will return [`Face::num_sides`] edges.
    fn boundary(&self) -> impl Iterator<Item = Self::HalfEdge> {
        let mut h = self.first_boundary();
        (0..self.num_sides()).map(move |_| {
            let res = h.clone();
            h = h.next();
            res
        })
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
/// Half-edges are oriented such that they go counter-clockwise around their [`HalfEdge::inside`]
/// face.
pub trait HalfEdge: Clone + Eq {
    /// The type of [`Face`] on the sphere.
    type Face: Face<Vertex = Self::Vertex, HalfEdge = Self>;

    /// The type of [`Vertex`] on the sphere.
    type Vertex: Vertex<Face = Self::Face, HalfEdge = Self>;

    /// Gets the [`Face`] whose interior boundary contains this half-edge.
    fn inside(&self) -> Self::Face;

    /// Gets the [`Vertex`] at the "start" of this half-edge.
    fn start(&self) -> Self::Vertex;

    /// Gets the complementary half-edge on the opposite side of the edge.
    ///
    /// The returned half-edge will go in the opposite direction along the same edge.
    fn complement(&self) -> Self;

    /// Gets the half-edge which shares the [`HalfEdge::inside`] face of this half-edge and
    /// precedes it in counter-clockwise order around the face.
    fn prev(&self) -> Self;

    /// Gets the half-edge which shares the [`HalfEdge::inside`] face of this half-edge and
    /// follows it in counter-clockwise order around the face.
    fn next(&self) -> Self;
}
