//! Contains types related to the [`IcoSphere`] tessellation.
use crate::statictri::{self, StaticTriSphere};
use crate::subtri::{self, SubTriSphere};

/// A tessellation of the unit sphere into geodesic triangular [`Face`]s constructed by
/// subdividing an icosahedron.
pub type IcoSphere = SubTriSphere<BaseIcoSphere>;

/// Represents a face of an [`IcoSphere`].
pub type Face = subtri::Face<BaseIcoSphere>;

/// Represents a vertex of an [`IcoSphere`].
pub type Vertex = subtri::Vertex<BaseIcoSphere>;

/// Represents a half-edge of an [`IcoSphere`].
pub type HalfEdge = subtri::HalfEdge<BaseIcoSphere>;

/// Partitions the surface of the unit sphere into a set of geodesic triangular [`Face`]s
/// by projecting an icosahedron onto it.
#[derive(Default, Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BaseIcoSphere;

/// The number of faces on a [`BaseIcoSphere`].
const NUM_BASE_FACES: usize = 20;

/// The number of vertices on a [`BaseIcoSphere`].
const NUM_BASE_VERTS: usize = 12;

/// Represents a face of a [`BaseIcoSphere`].
pub type BaseFace = statictri::Face<NUM_BASE_FACES, NUM_BASE_VERTS, BaseIcoSphere>;

/// Represents a vertex of a [`BaseIcoSphere`].
pub type BaseVertex = statictri::Vertex<NUM_BASE_FACES, NUM_BASE_VERTS, BaseIcoSphere>;

/// Represents a half-edge of a [`BaseIcoSphere`].
pub type BaseHalfEdge = statictri::HalfEdge<NUM_BASE_FACES, NUM_BASE_VERTS, BaseIcoSphere>;

impl AsRef<StaticTriSphere<NUM_BASE_FACES, NUM_BASE_VERTS>> for BaseIcoSphere {
    fn as_ref(&self) -> &StaticTriSphere<NUM_BASE_FACES, NUM_BASE_VERTS> {
        &BASE
    }
}

impl crate::Sphere for BaseIcoSphere {
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

impl crate::subtri::BaseSphere for BaseIcoSphere {}

impl crate::subtri::BaseSphereInternal for BaseIcoSphere {
    type Region = statictri::Region<NUM_BASE_FACES, NUM_BASE_VERTS, BaseIcoSphere>;

    fn face_owns_vertex_0(&self, face: Self::Face) -> bool {
        statictri::Deref(Self).face_owns_vertex_0(face)
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

    /// Interpolation function for a base face.
    ///
    /// Suppose the vertices of a base face are `v₀`, `v₁`, and `v₂`. The interpolation function
    /// is given weights `w₀`, `w₁`, and `w₂` such that `wᵢ ≥ 0` and `w₀ + w₁ + w₂ = 1`. It
    /// will compute an approximation of the "spherical" interpolation of the three vertices
    /// with respect to the given weights, `r = c₀ v₀ + c₁ v₁ + c₂ v₂` and return the coefficients
    /// `[c₀, c₁, c₂]`.
    ///
    /// It is guaranteed that:
    ///  * `r` is on the unit sphere
    ///  * If `wᵢ` is `1`, `cᵢ` will also be `1`.
    ///
    /// See [this paper](https://mathweb.ucsd.edu/~sbuss/ResearchWeb/spheremean/paper.pdf) for
    /// a precise definition of spherical interpolation (which we only approximate here, for
    /// performance reasons).
    fn interpolate(&self, weights: [f64; 3]) -> [f64; 3] {
        let [w_0, w_1, w_2] = weights;

        // Compute initial unnormalized coefficients. `cᵢ = (wᵢ * ANGLE).sin()` assures that this
        // behaves like a perfect SLERP (https://en.wikipedia.org/wiki/Slerp) when one of
        // the weights is `0`.
        // TODO: it might be possible to create a faster implementation of `(x * ANGLE).sin()`.
        // This could also make the result portable, since `sin` isn't.
        let c_0 = (w_0 * ANGLE).sin();
        let c_1 = (w_1 * ANGLE).sin();
        let c_2 = (w_2 * ANGLE).sin();

        // Normalize coefficients such that, when applied to the 3 vertices of the base face,
        // the result will be on the unit sphere
        let norm_sqr =
            c_0 * c_0 + c_1 * c_1 + c_2 * c_2 + 2.0 * DOT * (c_0 * c_1 + c_1 * c_2 + c_2 * c_0);
        let norm = norm_sqr.sqrt();
        let c_0 = c_0 / norm;
        let c_1 = c_1 / norm;
        let c_2 = c_2 / norm;
        return [c_0, c_1, c_2];

        /// Dot product of two adjacent vertices of the base icosahedron.
        const DOT: f64 = C_1;

        /// Angle between two adjacent vertices of the base icosahedron.
        const ANGLE: f64 = 1.1071487177940904;
    }
}

/// The unsubdivided base sphere of an [`IcoSphere`].
static BASE: StaticTriSphere<NUM_BASE_FACES, NUM_BASE_VERTS> = StaticTriSphere::new(
    [
        // Top apex
        [0.0, 0.0, 1.0],
        // Top ring
        [C_0, 0.0, C_1],
        [C_2, C_3, C_1],
        [-C_4, C_5, C_1],
        [-C_4, -C_5, C_1],
        [C_2, -C_3, C_1],
        // Bottom ring
        [C_4, -C_5, -C_1],
        [C_4, C_5, -C_1],
        [-C_2, C_3, -C_1],
        [-C_0, 0.0, -C_1],
        [-C_2, -C_3, -C_1],
        // Bottom apex
        [0.0, 0.0, -1.0],
    ],
    [
        // Top cap
        [0, 1, 2],
        [0, 2, 3],
        [0, 3, 4],
        [0, 4, 5],
        [0, 5, 1],
        // Central ring
        [5, 6, 1],
        [6, 7, 1],
        [1, 7, 2],
        [7, 8, 2],
        [2, 8, 3],
        [8, 9, 3],
        [3, 9, 4],
        [9, 10, 4],
        [4, 10, 5],
        [10, 6, 5],
        // Bottom cap
        [11, 6, 10],
        [11, 10, 9],
        [11, 9, 8],
        [11, 8, 7],
        [11, 7, 6],
    ],
);

/// `sqrt(4 / 5)`
const C_0: f64 = 0.8944271909999159;

/// `sqrt(1 / 5)`
const C_1: f64 = 0.4472135954999579;

/// `(5 - sqrt(5)) / 10`
const C_2: f64 = 0.276393202250021;

/// `sqrt((5 + sqrt(5)) / 10)`
const C_3: f64 = 0.8506508083520399;

/// `(5 + sqrt(5)) / 10`
const C_4: f64 = 0.7236067977499789;

/// `sqrt((5 - sqrt(5)) / 10)`
const C_5: f64 = 0.5257311121191336;