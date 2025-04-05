//! Contains types related to the [`OctoSphere`] tessellation.
use crate::statictri::{self, StaticTriSphere};
use crate::subtri::{self, SubTriSphere};

/// A tessellation of the unit sphere into geodesic triangular [`Face`]s constructed by
/// subdividing an octohedron.
pub type OctoSphere = SubTriSphere<BaseOctoSphere>;

/// Represents a face of an [`OctoSphere`].
pub type Face = subtri::Face<BaseOctoSphere>;

/// Represents a vertex of an [`OctoSphere`].
pub type Vertex = subtri::Vertex<BaseOctoSphere>;

/// Represents a half-edge of an [`OctoSphere`].
pub type HalfEdge = subtri::HalfEdge<BaseOctoSphere>;

/// Partitions the surface of the unit sphere into a set of geodesic triangular [`Face`]s
/// by projecting an octohedron onto it.
#[derive(Default, Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BaseOctoSphere;

/// The number of faces on a [`BaseOctoSphere`].
const NUM_BASE_FACES: usize = 8;

/// The number of vertices on a [`BaseOctoSphere`].
const NUM_BASE_VERTS: usize = 6;

impl AsRef<StaticTriSphere<NUM_BASE_FACES, NUM_BASE_VERTS>> for BaseOctoSphere {
    fn as_ref(&self) -> &StaticTriSphere<NUM_BASE_FACES, NUM_BASE_VERTS> {
        &BASE
    }
}

/// The unsubdivided base sphere of an [`OctoSphere`].
static BASE: StaticTriSphere<NUM_BASE_FACES, NUM_BASE_VERTS> = StaticTriSphere::new(
    [
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, -1.0],
    ],
    [
        // Top cap
        [0, 1, 2],
        [0, 2, 3],
        [0, 3, 4],
        [0, 4, 1],
        // Bottom cap
        [4, 5, 1],
        [1, 5, 2],
        [2, 5, 3],
        [3, 5, 4],
    ],
    EDGE_ANGLE,
    EDGE_COS_ANGLE
);

/// The angle subtended by the edges of the base sphere.
const EDGE_ANGLE: f64 = std::f64::consts::FRAC_PI_2;

/// The cosine of the angle subtended by the edges of the base sphere.
const EDGE_COS_ANGLE: f64 = 0.0;

/// Represents a face of a [`BaseOctoSphere`].
pub type BaseFace = statictri::Face<NUM_BASE_FACES, NUM_BASE_VERTS, BaseOctoSphere>;

/// Represents a vertex of a [`BaseOctoSphere`].
pub type BaseVertex = statictri::Vertex<NUM_BASE_FACES, NUM_BASE_VERTS, BaseOctoSphere>;

/// Represents a half-edge of a [`BaseOctoSphere`].
pub type BaseHalfEdge = statictri::HalfEdge<NUM_BASE_FACES, NUM_BASE_VERTS, BaseOctoSphere>;

impl crate::Sphere for BaseOctoSphere {
    type Face = BaseFace;
    type Vertex = BaseVertex;
    type HalfEdge = BaseHalfEdge;

    fn num_faces(&self) -> usize {
        statictri::Deref(Self).num_faces()
    }

    fn faces(&self) -> impl Iterator<Item = BaseFace> {
        statictri::Deref(Self).faces()
    }

    fn num_vertices(&self) -> usize {
        statictri::Deref(Self).num_vertices()
    }

    fn vertices(&self) -> impl Iterator<Item = BaseVertex> {
        statictri::Deref(Self).vertices()
    }
}

impl crate::subtri::BaseSphere for BaseOctoSphere {}

impl crate::subtri::BaseSphereInternal for BaseOctoSphere {
    type Region = statictri::Region<NUM_BASE_FACES, NUM_BASE_VERTS, BaseOctoSphere>;

    fn edge_angle(&self) -> f64 {
        statictri::Deref(Self).edge_angle()
    }

    fn edge_cos_angle(&self) -> f64 {
        statictri::Deref(Self).edge_cos_angle()
    }

    fn face_owns_vertex_1(&self, face: Self::Face) -> bool {
        statictri::Deref(Self).face_owns_vertex_1(face)
    }

    fn face_owns_edge_2(&self, face: Self::Face) -> bool {
        statictri::Deref(Self).face_owns_edge_2(face)
    }

    fn vertex_owner(&self, vertex: Self::Vertex) -> Self::Face {
        statictri::Deref(Self).vertex_owner(vertex)
    }

    fn num_owned_vertices_before(&self, face: Self::Face) -> usize {
        statictri::Deref(Self).num_owned_vertices_before(face)
    }

    fn num_owned_edges_before(&self, face: Self::Face) -> usize {
        statictri::Deref(Self).num_owned_edges_before(face)
    }

    fn first_region(&self) -> Self::Region {
        statictri::Deref(Self).first_region()
    }

    fn next_region(&self, region: Self::Region) -> Option<Self::Region> {
        statictri::Deref(Self).next_region(region)
    }
}