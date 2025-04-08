//! Contains types related to projections.
use crate::basetri;
#[expect(unused_imports)]
use crate::basetri::BaseTriSphere;
use crate::math::{mat, vec};
#[expect(unused_imports)]
use crate::tri::TriSphere;

pub use gnomonic::Gnomonic;
pub use tri::{BaseTriSphereProjection, TriSphereProjection};

/// Maps local coordinates on a planar equilateral triangle to coordinates on the unit sphere.
///
/// The local coordinates are defined like so:
/// ```text
///         [0, a]
///          /  \
///         /    \
///        /      \
///       /        \
///      /          \
///     /            \
/// [0, 0] -------- [a, 0]
/// ```
pub trait TriangleProjection {
    /// Projects a point in the local coordinates of the triangle to a point on the unit sphere.
    ///
    /// This may assume that the given coordinates, `[u, v]` satisfy `0 <= u`, `0 <= v`, and
    /// `u + v <= a`.
    fn to_sphere(&self, coords: [f64; 2]) -> [f64; 3];

    /// Projects a point on the unit sphere to a point in the local coordinates of the triangle.
    #[expect(clippy::wrong_self_convention)]
    fn from_sphere(&self, point: [f64; 3]) -> [f64; 2];
}

/// Maps local coordinates on a planar 60-degree parallelogram to coordinates on the unit sphere.
///
/// The local coordinates are defined like so:
/// ```text
///         [0, c] ------- [b, c]
///          /              /
///         /              /
///        /              /
///       /              /
///      /              /
///     /              /
/// [0, 0] -------- [b, 0]
/// ```
pub trait ParallelogramProjection {
    /// Projects a point in the local coordinates of the parallelogram to a point on the unit
    /// sphere.
    ///
    /// This may assume that the given coordinates, `[u, v]` satisfy `0 <= u <= b` and
    /// `0 <= v <= c`.
    fn to_sphere(&self, coords: [f64; 2]) -> [f64; 3];

    /// Projects a point on the unit sphere to a point in the local coordinates of the
    /// parallelogram.
    #[expect(clippy::wrong_self_convention)]
    fn from_sphere(&self, point: [f64; 3]) -> [f64; 2];
}

/// Contains projection-related types and traits for [`BaseTriSphere`] and [`TriSphere`].
pub mod tri {
    use super::*;

    /// A potential projection function for a [`BaseTriSphere`].
    ///
    /// This defines a mapping from the local planar coordinates of each of its
    /// [`Face`](basetri::Face)s to coordinates on the unit sphere.
    pub trait BaseTriSphereProjection: TriSphereProjection {
        /// The type of [`TriangleProjection`] used for a [`basetri::Face`].
        type Face: TriangleProjection;

        /// Gets the [`TriangleProjection`] for the given [`basetri::Face`], assuming that the
        /// scale of the triangle (`a`) is `1`.
        fn face(&self, face: basetri::Face) -> Self::Face;
    }

    /// A potential projection function for a [`TriSphere`].
    ///
    /// This defines a mapping from the local planar coordinates of each of its regions to
    /// coordinates on the unit sphere. Unlike [`BaseTriSphere`], this actually affects the shape
    /// of its faces, so a projection must be specified to construct a [`TriSphere`].
    pub trait TriSphereProjection {
        /// The type of [`TriangleProjection`] used for the interior of a [`basetri::Face`].
        type Interior: TriangleProjection;

        /// Gets the [`TriangleProjection`] for points on the "interior" region of a particular
        /// [`basetri::Face`].
        ///
        /// The `b` and `c` parameters of the [`TriSphere`] are given. The number of segments along
        /// each edge of the interior triangle is `b - c`.
        fn interior(&self, b: u32, c: u32, face: basetri::Face) -> Self::Interior;

        /// The type of [`ParallelogramProjection`] used for a [`basetri::HalfEdge`].
        type Edge: ParallelogramProjection;

        /// Gets the [`ParallelogramProjection`] for points on an "edge" region identified by a
        /// particular [`basetri::HalfEdge`].
        ///
        /// The `b` and `c` parameters of the [`TriSphere`] are given. These specify the number
        /// of segments along the sides of the parallelogram.
        fn edge(&self, b: u32, c: u32, edge: basetri::HalfEdge) -> Self::Edge;
    }

    /// Implements [`TriangleProjection`] for the interior region of a base face, given the
    /// [`TriangleProjection`] for that face.
    #[derive(Clone, Copy, Debug)]
    pub struct Interior<T> {
        face: T,
        p_0: [f64; 2],
        u: [f64; 2],
        v: [f64; 2],
    }

    impl<T> Interior<T> {
        /// Constructs a [`Interior`] wrapper based on the given face projection and [`TriSphere`]
        /// parameters.
        pub fn new(b: u32, c: u32, face: T) -> Self {
            let (p_0, [u, v]) = interior_trans(b, c);
            Self { face, p_0, u, v }
        }
    }

    impl<T: TriangleProjection> TriangleProjection for Interior<T> {
        fn to_sphere(&self, coords: [f64; 2]) -> [f64; 3] {
            self.face
                .to_sphere(vec::add(self.p_0, mat::apply([self.u, self.v], coords)))
        }

        fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
            todo!()
        }
    }

    /// Implements [`ParallelogramProjection`] for an edge region, given the
    /// [`TriangleProjection`]s for its connected base faces.
    #[derive(Clone, Copy, Debug)]
    pub struct Edge<T> {
        face_left: T,
        face_right: T,
        u: [f64; 2],
        v: [f64; 2],
    }

    impl<T> Edge<T> {
        /// Constructs a [`Edge`] wrapper based on the given connected face projections and
        /// [`TriSphere`] parameters.
        pub fn new(b: u32, c: u32, face_left: T, face_right: T) -> Self {
            let (_, [u, v]) = interior_trans(b, c);
            Self {
                face_left,
                face_right,
                u,
                v,
            }
        }
    }

    /// Computes the transformation from an interior regions coordinates to the coordinates
    /// of the face that contains it.
    fn interior_trans(b: u32, c: u32) -> ([f64; 2], [[f64; 2]; 2]) {
        let b = b as f64;
        let c = c as f64;
        let w_total = b * b + b * c + c * c;
        let w_1 = c * c / w_total;
        let w_2 = b * c / w_total;
        let p_0 = [w_1, w_2];
        let u = vec::div(vec::sub([1.0, 0.0], p_0), b);
        let v = if c > 0.0 {
            vec::div(p_0, c)
        } else {
            [0.0, 1.0 / b]
        };
        (p_0, [u, v])
    }

    impl<T: TriangleProjection> ParallelogramProjection for Edge<T> {
        fn to_sphere(&self, coords: [f64; 2]) -> [f64; 3] {
            let [u, v] = mat::apply([self.u, self.v], coords);
            if v >= 0.0 {
                self.face_left.to_sphere([u, v])
            } else {
                self.face_right.to_sphere([1.0 - u, -v])
            }
        }

        fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
            todo!()
        }
    }
}

/// Contains types related to the [`Gnomonic`] projection.
pub mod gnomonic {
    use super::*;
    use crate::Vertex;

    /// The [gnomonic](https://en.wikipedia.org/wiki/Gnomonic_projection) projection.
    ///
    /// This is the simplest and fastest projection. Given a point on a face of a
    /// [`BaseTriSphere`](super::BaseTriSphere`), it first maps the point to the planar triangle,
    /// and then projects that point to the unit sphere using a perspective projection through the
    /// origin.
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
    pub struct Gnomonic;

    impl Gnomonic {
        /// Constructs a [`Planar`] projection for the inside face of the given edge.
        ///
        /// The projection will be oriented such that the local origin is at the start of the edge.
        pub fn inside(edge: basetri::HalfEdge) -> Planar {
            let p_0 = edge.start().pos();
            let p_1 = edge.next().start().pos();
            let p_2 = edge.prev().start().pos();
            Planar {
                offset: p_0,
                u: vec::sub(p_1, p_0),
                v: vec::sub(p_2, p_0),
            }
        }
    }

    impl BaseTriSphereProjection for Gnomonic {
        type Face = Planar;
        fn face(&self, face: basetri::Face) -> Planar {
            Self::inside(face.side(0))
        }
    }

    impl TriSphereProjection for Gnomonic {
        type Interior = tri::Interior<Planar>;
        fn interior(&self, b: u32, c: u32, face: basetri::Face) -> tri::Interior<Planar> {
            tri::Interior::new(b, c, self.face(face))
        }

        type Edge = tri::Edge<Planar>;
        fn edge(&self, b: u32, c: u32, edge: basetri::HalfEdge) -> Self::Edge {
            let face_left = Self::inside(edge);
            let face_right = Self::inside(edge.complement());
            tri::Edge::new(b, c, face_left, face_right)
        }
    }

    /// A [`TriangleProjection`] or [`ParallelogramProjection`] for a planar region using
    /// the [`Gnomonic`] projection.
    #[derive(Debug)]
    pub struct Planar {
        /// The position of the origin of the local coordinate space of the triangle.
        offset: [f64; 3],

        /// The movement of the unnormalized point for each step in the `u` direction of the
        /// local coordinate space.
        u: [f64; 3],

        /// The movement of the unnormalized point for each step in the `v` direction of the
        /// local coordinate space.
        v: [f64; 3],
    }

    impl Planar {
        /// Projects a point in the local coordinates of the planar region to a point on the unit
        /// sphere.
        fn to_sphere(&self, [u, v]: [f64; 2]) -> [f64; 3] {
            vec::normalize(vec::add(self.offset, mat::apply([self.u, self.v], [u, v])))
        }

        /// Projects a point on the unit sphere to a point in the local coordinates of the
        /// planar region.
        #[expect(clippy::wrong_self_convention)]
        fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
            todo!()
        }
    }

    impl TriangleProjection for Planar {
        fn to_sphere(&self, [u, v]: [f64; 2]) -> [f64; 3] {
            Planar::to_sphere(self, [u, v])
        }

        fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
            Planar::from_sphere(self, point)
        }
    }

    impl ParallelogramProjection for Planar {
        fn to_sphere(&self, [u, v]: [f64; 2]) -> [f64; 3] {
            Planar::to_sphere(self, [u, v])
        }

        fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
            Planar::from_sphere(self, point)
        }
    }
}
