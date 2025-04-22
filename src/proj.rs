//! Contains types related to projections.
use crate::basetri;
#[expect(unused_imports)]
use crate::basetri::BaseTriSphere;
use crate::math::{self, mat, vec};
#[expect(unused_imports)]
use crate::tri::TriSphere;

pub use fuller::Fuller;
pub use gnomonic::Gnomonic;
pub use tri::{BaseTriProjector, TriProjector};

/// Maps two-dimensional "local" coordinates to points on the unit sphere.
/// 
/// The domain and range of the projection function may be restricted, e.g. local coordinates
/// may be restricted to a triangle or rectangle, and only a small section of the unit sphere
/// may be covered.
pub trait Projection {
    /// Projects a point in the local coordinates to a point on the unit sphere.
    /// 
    /// This is the inverse of [`from_sphere`](Projection::from_sphere). If the given local
    /// coordinates are not in the domain of the projection, the result is undefined.
    fn to_sphere(&self, coords: [f64; 2]) -> [f64; 3];

    /// Projects a point on the unit sphere to a point in the local coordinates.
    /// 
    /// This is the inverse of [`to_sphere`](Projection::to_sphere). If the given point is not in
    /// the range of the projection, the result is undefined.
    #[expect(clippy::wrong_self_convention)]
    fn from_sphere(&self, point: [f64; 3]) -> [f64; 2];

    /// The type of [`Projection`] resulting from a call to [`transform`](Projection::transform).
    type Transform: Projection;

    /// Applies an affine transformation to the local coordinates of this projection.
    ///
    /// More concretely, the [`to_sphere`](Projection::transform) method on the returned
    /// projection will compute `self.to_sphere(offset + linear * input)`, where `offset` is
    /// interpreted as a vector and `linear` is interpreted as a 2x2 matrix column-major matrix.
    fn transform(self, offset: [f64; 2], linear: [[f64; 2]; 2]) -> Self::Transform;
}

/// A general implementation of [`Projection::transform`].
#[derive(Debug)]
pub struct Transform<T> {
    source: T,
    offset: [f64; 2],
    linear: [[f64; 2]; 2],
}

impl<T> Transform<T> {
    /// Constructs a transformed view of the given source projection.
    ///
    /// This is always a valid implementation of `source.transform(offset, linear)`.
    pub fn new(source: T, offset: [f64; 2], linear: [[f64; 2]; 2]) -> Self {
        Self {
            source,
            offset,
            linear,
        }
    }
}

impl<T: Projection> Projection for Transform<T> {
    fn to_sphere(&self, coords: [f64; 2]) -> [f64; 3] {
        self.source
            .to_sphere(vec::add(self.offset, mat::apply(self.linear, coords)))
    }

    fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
        todo!()
    }

    type Transform = Self;
    fn transform(self, offset: [f64; 2], linear: [[f64; 2]; 2]) -> Self {
        todo!()
    }
}

/// Contains projection-related types and traits for [`BaseTriSphere`] and [`TriSphere`].
pub mod tri {
    use super::*;

    /// A general projection method which can be used to create a [`Projection`] for any
    /// spherical triangle.
    pub trait TriProjector: BaseTriProjector {
        /// Constructs a [`Projection`] for a specified spherical triangle.
        ///
        /// The points must be on the unit sphere, in counter-clockwise order, and strictly
        /// contained in one hemisphere. The edges of the triangle are geodesics. The length
        /// of each edge in the local coordinate space is assumed to be `1`.
        /// 
        /// The local coordinates are defined like so:
        /// ```text
        ///           p₂
        ///         [0, 1]
        ///          /  \
        ///         /    \
        ///        /      \
        ///       /        \
        ///      /          \
        ///     /            \
        /// [0, 0] -------- [1, 0]
        ///   p₀              p₁
        /// ```
        fn triangle(&self, points: [[f64; 3]; 3]) -> Self::Triangle;
    }

    /// A general projection method which can be used to create a [`Projection`] for any face of a
    /// [`BaseTriSphere`].
    pub trait BaseTriProjector {
        /// The type of [`Projection`] used for a [`basetri::Face`].
        type Triangle: Projection;

        /// Constructs a [`Projection`] for the [`inside`](crate::HalfEdge::inside) face
        /// of the given edge, oriented such that the U axis goes along the edge.
        ///
        /// The length of each edge in the local coordinate space is assumed to be `1`.
        /// 
        /// The local coordinates are defined like so:
        /// ```text
        ///         [0, 1]
        ///          /  \
        ///         /    \
        ///        V      \
        ///       /        \
        ///      /          \
        ///     /            \
        /// [0, 0] --- U -- [1, 0]
        /// ```
        fn inside(&self, edge: basetri::HalfEdge) -> Self::Triangle;
    }
}

/// Contains types related to the [`Gnomonic`] projection.
pub mod gnomonic {
    use super::*;

    /// The [gnomonic](https://en.wikipedia.org/wiki/Gnomonic_projection) projection.
    ///
    /// This is the simplest and fastest projection. Given a point on a face of a
    /// [`BaseTriSphere`](super::BaseTriSphere`), it first maps the point to the planar triangle,
    /// and then projects that point to the unit sphere using a perspective projection through the
    /// origin.
    ///
    /// This projection has significant area inflation near the centers of faces and has
    /// kinks (C₁ discontinuities) along the edges of faces.
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
    pub struct Gnomonic;

    impl TriProjector for Gnomonic {
        fn triangle(&self, [p_0, p_1, p_2]: [[f64; 3]; 3]) -> Planar {
            let u = vec::sub(p_1, p_0);
            let v = vec::sub(p_2, p_0);
            let normal = vec::cross(u, v);
            let normal = vec::div(normal, vec::dot(normal, p_0));
            Planar {
                offset: p_0,
                normal,
                u,
                v,
                inv_linear: inv_linear(u, v),
            }
        }
    }

    impl BaseTriProjector for Gnomonic {
        type Triangle = Planar;
        fn inside(&self, edge: basetri::HalfEdge) -> Planar {
            self.triangle([
                edge.start().pos(),
                edge.next().start().pos(),
                edge.prev().start().pos(),
            ])
        }
    }

    /// A [`Projection`] for a planar region using the [`Gnomonic`] projection.
    #[derive(Clone, Copy, Debug)]
    pub struct Planar {
        /// The position of the origin of the local coordinate space of the triangle.
        offset: [f64; 3],

        /// A vector which is both normal to the plane and on the plane.
        /// 
        /// It follows that every vector `p` which is on the plane satisfies `dot(normal, p) = 1`.
        normal: [f64; 3],

        /// The movement of the unnormalized point for each step in the `u` direction of the
        /// local coordinate space.
        u: [f64; 3],

        /// The movement of the unnormalized point for each step in the `v` direction of the
        /// local coordinate space.
        v: [f64; 3],

        /// The inverse of the matrix:
        /// 
        /// ```text
        /// [dot(u, u), dot(u, v)]
        /// [dot(u, v), dot(v, v)]
        /// ```
        /// 
        /// Used for implementing [`Projection::from_sphere`].
        inv_linear: [[f64; 2]; 2],
    }

    impl Projection for Planar {
        fn to_sphere(&self, [u, v]: [f64; 2]) -> [f64; 3] {
            vec::normalize(vec::add(self.offset, mat::apply([self.u, self.v], [u, v])))
        }

        fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
            let point = vec::div(point, vec::dot(self.normal, point));
            let diff = vec::sub(point, self.offset);
            let dot_u = vec::dot(self.u, diff);
            let dot_v = vec::dot(self.v, diff);
            mat::apply(self.inv_linear, [dot_u, dot_v])
        }

        type Transform = Self;
        fn transform(self, offset: [f64; 2], linear: [[f64; 2]; 2]) -> Self {
            let u = mat::apply([self.u, self.v], linear[0]);
            let v = mat::apply([self.u, self.v], linear[1]);
            Planar {
                offset: vec::add(self.offset, mat::apply([self.u, self.v], offset)),
                normal: self.normal,
                u,
                v,
                inv_linear: inv_linear(u, v),
            }
        }
    }

    /// Computes [`Planar::inv_linear`] based on the vectors `u` and `v`.
    fn inv_linear(u: [f64; 3], v: [f64; 3]) -> [[f64; 2]; 2] {
        let dot_u_u = vec::dot(u, u);
        let dot_v_v = vec::dot(v, v);
        let dot_u_v = vec::dot(u, v);
        let det = dot_u_u * dot_v_v - dot_u_v * dot_u_v;
        [
            [dot_v_v / det, -dot_u_v / det],
            [-dot_u_v / det, dot_u_u / det],
        ]
    }
}

/// Contains types related to the [`Fuller`] projection.
pub mod fuller {
    use super::*;

    /// The [Fuller](https://en.wikipedia.org/wiki/Dymaxion_map) projection method.
    ///
    /// This projection preserves distances exactly along the boundaries of base triangles.
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
    pub struct Fuller;

    impl Fuller {
        /// Constructs a [`Projection`] for a specified spherical triangle. The triangle
        /// must be equilateral, with the angle between any two endpoints equal to `angle`.
        pub fn triangle(self, angle: f64, points: [[f64; 3]; 3]) -> Triangle {
            // TODO: Validate triangle when debug assertions are on
            Triangle { angle, points }
        }
    }

    impl BaseTriProjector for Fuller {
        type Triangle = Triangle;
        fn inside(&self, edge: basetri::HalfEdge) -> Triangle {
            self.triangle(
                edge.sphere().edge_length(),
                [
                    edge.start().pos(),
                    edge.next().start().pos(),
                    edge.prev().start().pos(),
                ],
            )
        }
    }

    /// A [`Projection`] for an equilateral triangle using the [`Fuller`] projection.
    #[derive(Clone, Copy, Debug)]
    pub struct Triangle {
        /// The angle, in radians, between any two endpoints of the triangle.
        angle: f64,

        /// The positions of the endpoints of the triangle on the sphere.
        points: [[f64; 3]; 3],
    }

    impl Projection for Triangle {
        fn to_sphere(&self, [u, v]: [f64; 2]) -> [f64; 3] {
            // Special handling for boundaries to improve accuracy and performance
            if v <= 0.0 {
                debug_assert_eq!(v, 0.0);
                if u <= 0.0 {
                    debug_assert_eq!(u, 0.0);
                    return self.points[0];
                } else {
                    return vec::normalize(vec::add(
                        vec::mul(self.points[0], ((1.0 - u) * self.angle).sin()),
                        vec::mul(self.points[1], (u * self.angle).sin()),
                    ));
                }
            } else if u <= 0.0 {
                debug_assert_eq!(u, 0.0);
                return vec::normalize(vec::add(
                    vec::mul(self.points[0], ((1.0 - v) * self.angle).sin()),
                    vec::mul(self.points[2], (v * self.angle).sin()),
                ));
            }

            // Adapted from "Exact Equations for Fuller’s Map Projection and Inverse"
            // https://utppublishing.com/doi/pdf/10.3138/carto.43.1.67

            // First, solve the equation:
            // `tan(x + proj_u) + tan(x) + tan(x + proj_v) + t_alpha = 0`
            let t_alpha = (self.angle / 2.0).tan();
            let proj_u = (2.0 * u + v - 1.0) * self.angle;
            let proj_v = (2.0 * v + u - 1.0) * self.angle;
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
            let a_0 = self.angle / 2.0 - x;
            let a_1 = a_0 - (2.0 * u + v - 1.0) * self.angle;

            // Construct planes by offsetting the second and third edges of the triangle
            let n_0 = vec::normalize(vec::cross(
                mat::apply(
                    [self.points[0], self.points[1]],
                    [(self.angle - a_0).sin(), a_0.sin()],
                ),
                vec::sub(self.points[2], self.points[1]),
            ));
            let n_1 = vec::normalize(vec::cross(
                mat::apply(
                    [self.points[1], self.points[2]],
                    [(self.angle - a_1).sin(), a_1.sin()],
                ),
                vec::sub(self.points[0], self.points[2]),
            ));

            // The projected point is the intersection of the two planes
            vec::normalize(vec::cross(n_0, n_1))
        }

        fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
            todo!()
        }

        type Transform = Transform<Self>;
        fn transform(self, offset: [f64; 2], linear: [[f64; 2]; 2]) -> Transform<Self> {
            Transform::new(self, offset, linear)
        }
    }
}
