//! Contains types related to projections.
use crate::basetri;
use crate::basetri::BaseTriSphere;
use crate::math::{self, mat, vec};
#[expect(unused_imports)]
use crate::tri::TriSphere;

pub use fuller::Fuller;
pub use gnomonic::Gnomonic;
pub use tri::{BaseTriSphereProjection, TriSphereProjection};

/// Maps local coordinates on a planar triangle to coordinates on the unit sphere.
///
/// The local coordinates are defined like so:
/// ```text
///           p₂
///         [0, a]
///          /  \
///         /    \
///        /      \
///       /        \
///      /          \
///     /            \
/// [0, 0] -------- [a, 0]
///   p₀              p₁
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

/// Maps local coordinates on a planar parallelogram to coordinates on the unit sphere.
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
    ///
    /// This projection has significant area inflation near the centers of faces and has
    /// kinks (discontinuous derivatives) along the edges of faces. These kinks are most
    /// noticeable when the `c` parameter of the [`TriSphere`] is non-zero.
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
    #[derive(Clone, Copy, Debug)]
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
        #[expect(clippy::wrong_self_convention)]
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

/// Contains types related to the [`Fuller`] projection.
pub mod fuller {
    use super::*;
    use crate::{Face, Vertex};

    /// The [Fuller](https://en.wikipedia.org/wiki/Dymaxion_map) projection.
    ///
    /// This projection preserves distances exactly along the boundaries of regions.
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
    pub struct Fuller;

    impl Fuller {
        /// Constructs an instance of the [`Fuller`] projection for a particular [`TriSphere`],
        /// given its [`BaseTriSphere`] and `c / b` ratio.
        pub fn sphere(base: BaseTriSphere, c_over_b: f64) -> Sphere {
            if c_over_b == 0.0 {
                let full_angle_b = base.vertex_angle();
                return Sphere {
                    q: [1.0, 0.0, 0.0],
                    full_angle_b,
                    full_angle_c: 0.0,
                };
            }

            // Compute initial estimate of the first interior point: `p₀ = q₀ v₀ + q₁ v₁ + q₂ v₂`
            let w_total = 1.0 + c_over_b + c_over_b * c_over_b;
            let w_0 = 1.0 / w_total;
            let w_2 = c_over_b * w_0;
            let w_1 = c_over_b * w_2;
            let angle = base.vertex_angle();
            let cos_angle = base.vertex_cos_angle();
            let q_0 = (w_0 * angle).sin();
            let q_1 = (w_1 * angle).sin();
            let q_2 = (w_2 * angle).sin();
            let mut q = vec::normalize_in(cos_angle, [q_0, q_1, q_2]);

            // For continuity between regions, it is very important that `p₂`, `p₀`,
            // and `v₀` all lie on the same plane, and that `angle(p₀, v₀) / angle(p₂, v₀) = c / b`.
            // To ensure this, we will refine our initial estimate above.
            let mut cos_full_angle_b = q[1] + (q[0] + q[2]) * cos_angle;
            let mut sin_full_angle_b = (1.0 - cos_full_angle_b * cos_full_angle_b).sqrt();
            let mut full_angle_b = cos_full_angle_b.acos();

            // TODO: Find a better algorithm that converges faster, or a closed form solution.
            for _ in 0..30 {
                q = vec::add(
                    [
                        ((1.0 - c_over_b) * full_angle_b).sin() / sin_full_angle_b,
                        0.0,
                        0.0,
                    ],
                    vec::mul(
                        [q[1], q[2], q[0]],
                        (c_over_b * full_angle_b).sin() / sin_full_angle_b,
                    ),
                );
                cos_full_angle_b = q[1] + (q[0] + q[2]) * cos_angle;
                sin_full_angle_b = (1.0 - cos_full_angle_b * cos_full_angle_b).sqrt();
                full_angle_b = cos_full_angle_b.acos();
            }

            // Compute angles of region edges
            let full_angle_c = full_angle_b * c_over_b;
            Sphere {
                q,
                full_angle_b,
                full_angle_c,
            }
        }
    }

    impl TriSphereProjection for Fuller {
        type Interior = Triangle;
        fn interior(&self, b: u32, c: u32, face: basetri::Face) -> Triangle {
            Self::sphere(face.sphere(), c as f64 / b as f64).interior(b, c, face)
        }

        type Edge = Parallelogram;
        fn edge(&self, b: u32, c: u32, edge: basetri::HalfEdge) -> Parallelogram {
            Self::sphere(edge.sphere(), c as f64 / b as f64).edge(b, c, edge)
        }
    }

    /// An instance of the [`Fuller`] projection for a particular [`BaseTriSphere`] and
    /// `c / b` ratio.
    ///
    /// For best performance, this should be used instead of [`Fuller`] when performing
    /// operations on the same sphere multiple times.
    #[derive(Clone, Copy, Debug)]
    pub struct Sphere {
        /// The position of the origin (i.e. first point) of an interior region for a base face,
        /// expressed in the basis consisting of the vertices of the face. Due to symmetry, the
        /// other points of the interior region can be computed by rotating this vector:
        ///  * `p₀ = q₀ v₀ + q₁ v₁ + q₂ v₂`
        ///  * `p₁ = q₀ v₁ + q₁ v₂ + q₂ v₀`
        ///  * `p₂ = q₀ v₂ + q₁ v₀ + q₂ v₁`
        q: [f64; 3],

        /// The angle, in radians, subtended by an edge corresponding to the `b` parameter of the
        /// [`TriSphere`].
        full_angle_b: f64,

        /// The angle, in radians, subtended by an edge corresponding to the `c` parameter of the
        /// [`TriSphere`].
        full_angle_c: f64,
    }

    impl TriSphereProjection for Sphere {
        type Interior = Triangle;
        fn interior(&self, b: u32, c: u32, face: basetri::Face) -> Self::Interior {
            debug_assert!(
                ((c as f64 / b as f64) - self.full_angle_c / self.full_angle_b).abs() < 1e-10,
                "sphere-specific projection is not compatible with subdivision parameters"
            );
            let v_0 = face.vertex(0).pos();
            let v_1 = face.vertex(1).pos();
            let v_2 = face.vertex(2).pos();
            Triangle {
                p: [
                    mat::apply([v_0, v_1, v_2], self.q),
                    mat::apply([v_1, v_2, v_0], self.q),
                    mat::apply([v_2, v_0, v_1], self.q),
                ],
                step_angle: self.full_angle_b / b as f64,
                full_angle: self.full_angle_b - self.full_angle_c,
            }
        }

        type Edge = Parallelogram;
        fn edge(&self, b: u32, c: u32, edge: basetri::HalfEdge) -> Self::Edge {
            debug_assert!(
                ((c as f64 / b as f64) - self.full_angle_c / self.full_angle_b).abs() < 1e-10,
                "sphere-specific projection is not compatible with subdivision parameters"
            );
            let comp = edge.complement();
            let v_0 = edge.start().pos();
            let v_1 = comp.start().pos();
            let v_left = edge.prev().start().pos();
            let v_right = comp.prev().start().pos();
            Parallelogram {
                p_0_0: v_0,
                p_1_1: v_1,
                p_0_1: mat::apply([v_0, v_1, v_left], self.q),
                p_1_0: mat::apply([v_1, v_0, v_right], self.q),
                step_angle: self.full_angle_b / b as f64,
                full_angle_u: self.full_angle_b,
                full_angle_v: safe_angle(self.full_angle_c),
            }
        }
    }

    /// A [`TriangleProjection`] for an equilateral triangle using the [`Fuller`] projection.
    #[derive(Clone, Copy, Debug)]
    pub struct Triangle {
        /// The positions of the endpoints of the triangle on the sphere.
        p: [[f64; 3]; 3],

        /// The angle, in radians, subtended by a step of length `1` in the local coordinate space
        /// of this triangle.
        step_angle: f64,

        /// The angle, in radians, between any two endpoints of the triangle.
        full_angle: f64,
    }

    impl TriangleProjection for Triangle {
        fn to_sphere(&self, [u, v]: [f64; 2]) -> [f64; 3] {
            // Special case for zero-size triangle
            if self.full_angle == 0.0 {
                return self.p[0];
            }

            // Adapted from "Exact Equations for Fuller’s Map Projection and Inverse"
            // https://utppublishing.com/doi/pdf/10.3138/carto.43.1.67

            // First, solve the equation:
            // `tan(x + proj_u) + tan(x) + tan(x + proj_v) + t_alpha = 0`
            let t_alpha = (self.full_angle / 2.0).tan();
            let proj_u = (2.0 * u + v) * self.step_angle - self.full_angle;
            let proj_v = (2.0 * v + u) * self.step_angle - self.full_angle;
            let t_half_proj_u = (proj_u / 2.0).tan();
            let t_half_proj_v = (proj_v / 2.0).tan();
            let den_u = 1.0 - t_half_proj_u * t_half_proj_u;
            let den_v = 1.0 - t_half_proj_v * t_half_proj_v;
            let s = t_half_proj_u * den_v + t_half_proj_v * den_u;
            let p = t_half_proj_u * t_half_proj_v;
            let (t_x, _) = math::solve_cubic(
                4.0 * p,
                4.0 * t_alpha * p - 4.0 * s,
                3.0 * den_u * den_v - 8.0 * p - 2.0 * t_alpha * s,
                t_alpha * den_u * den_v + 2.0 * s,
            )
            .map(|t_x| {
                // There may be up to three solutions to the cubic equation, but only one
                // corresponds to the correct projection solution. This is the unique solution
                // which satisfies:
                //   * `-self.full_angle / 2.0 <= x <= self.full_angle / 2.0`
                //   * `-self.full_angle / 2.0 <= x + proj_u <= self.full_angle / 2.0`
                //   * `-self.full_angle / 2.0 <= x + proj_v <= self.full_angle / 2.0`
                //
                // or equivalently:
                //   * `-t_alpha <= tan(x) <= t_alpha`
                //   * `-t_alpha <= tan(x + proj_u) <= t_alpha`
                //   * `-t_alpha <= tan(x + proj_v) <= t_alpha`

                // Assuming there is exactly one such solution, it will be the one with the
                // smallest `max(|tan(x)|, |tan(x + proj_u)|, |tan(x + proj_v)|)`
                let t_x_proj_u =
                    (den_u + 2.0 * t_half_proj_u) / (den_u - 2.0 * t_x * t_half_proj_u);
                let t_x_proj_v =
                    (den_v + 2.0 * t_half_proj_v) / (den_v - 2.0 * t_x * t_half_proj_v);
                (t_x, t_x.abs().max(t_x_proj_u.abs()).max(t_x_proj_v.abs()))
            })
            .min_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
            .unwrap();
            let x = t_x.atan();

            // Use the solution to compute angle offsets of two planes from the second and
            // third edges of the triangle
            let a_0 = self.full_angle / 2.0 - x;
            let a_1 = self.full_angle + a_0 - (2.0 * u + v) * self.step_angle;

            // Construct planes by offsetting the second and third edges of the triangle
            let n_0 = vec::normalize(vec::cross(
                mat::apply(
                    [self.p[0], self.p[1]],
                    [(self.full_angle - a_0).sin(), a_0.sin()],
                ),
                vec::sub(self.p[2], self.p[1]),
            ));
            let n_1 = vec::normalize(vec::cross(
                mat::apply(
                    [self.p[1], self.p[2]],
                    [(self.full_angle - a_1).sin(), a_1.sin()],
                ),
                vec::sub(self.p[0], self.p[2]),
            ));

            // The projected point is the intersection of the two planes
            vec::normalize(vec::cross(n_0, n_1))
        }

        fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
            todo!()
        }
    }

    /// A [`ParallelogramProjection`] for a parallelogram using the [`Fuller`] projection.
    #[derive(Clone, Copy, Debug)]
    pub struct Parallelogram {
        /// The position of the origin of the region on the sphere.
        p_0_0: [f64; 3],

        /// The position of the far (U = 1, V = 1) corner of the region on the sphere.
        p_1_1: [f64; 3],

        /// The position of the left (U = 0, V = 1) corner of the region on the sphere.
        p_0_1: [f64; 3],

        /// The position of the right (U = 1, V = 0) corner of the region on the sphere.
        p_1_0: [f64; 3],

        /// The angle, in radians, subtended by a step of length `1` in the local coordinate space of
        /// this region.
        step_angle: f64,

        /// The angle, in radians, subtended by the U edge of this region.
        full_angle_u: f64,

        /// The angle, in radians, subtended by the V edge of this region.
        full_angle_v: f64,
    }

    impl ParallelogramProjection for Parallelogram {
        fn to_sphere(&self, [u, v]: [f64; 2]) -> [f64; 3] {
            // TODO: This isn't quite right. We need to construct planes and intersect them.
            let c_u_0 = (self.full_angle_u - u * self.step_angle).sin();
            let c_u_1 = (u * self.step_angle).sin();
            let c_v_0 = (self.full_angle_v - v * self.step_angle).sin();
            let c_v_1 = (v * self.step_angle).sin();
            vec::normalize(vec::add(
                vec::add(
                    vec::mul(self.p_0_0, c_u_0 * c_v_0),
                    vec::mul(self.p_1_0, c_u_1 * c_v_0),
                ),
                vec::add(
                    vec::mul(self.p_0_1, c_u_0 * c_v_1),
                    vec::mul(self.p_1_1, c_u_1 * c_v_1),
                ),
            ))
        }

        fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
            todo!()
        }
    }

    /// Checks if `angle` is zero. If so, replaces it with an arbitrary positive value. Otherwise
    /// returns it.
    fn safe_angle(angle: f64) -> f64 {
        if angle == 0.0 { 1.0 } else { angle }
    }
}
