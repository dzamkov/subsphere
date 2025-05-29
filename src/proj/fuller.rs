//! Contains types related to the [`Fuller`] projection.
use super::*;
use std::f64::consts::FRAC_PI_3;

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
        let sin_half_angle = (angle / 2.0).sin();
        let cos_half_angle = (angle / 2.0).cos();
        let tan_half_angle = sin_half_angle / cos_half_angle;
        Triangle {
            angle,
            sin_half_angle,
            cos_half_angle,
            tan_half_angle,
            points,
        }
    }
}

impl BaseTriProjector for Fuller {
    type Triangle = Triangle;
    fn inside(&self, edge: basetri::HalfEdge) -> Triangle {
        let consts = edge.sphere().consts();
        Triangle {
            angle: consts.angle,
            sin_half_angle: consts.sin_half_angle,
            cos_half_angle: consts.cos_half_angle,
            tan_half_angle: consts.tan_half_angle,
            points: [
                edge.start().pos(),
                edge.next().start().pos(),
                edge.prev().start().pos(),
            ],
        }
    }
}

/// A [`Projection`] for an equilateral triangle using the [`Fuller`] projection.
#[derive(Clone, Copy, Debug)]
pub struct Triangle {
    /// The angle, in radians, between any two endpoints of the triangle.
    angle: f64,

    /// `sin(angle / 2)`
    sin_half_angle: f64,

    /// `cos(angle / 2)`
    cos_half_angle: f64,

    /// `tan(angle / 2)`
    tan_half_angle: f64,

    /// The positions of the endpoints of the triangle on the sphere.
    points: [[f64; 3]; 3],
}

impl Projection for Triangle {
    fn to_sphere(&self, [u, v]: [f64; 2]) -> [f64; 3] {
        // Special handling for boundaries to improve accuracy and performance
        if v <= 0.0 {
            debug_assert_eq!(v, 0.0);
            return if u <= 0.0 {
                debug_assert_eq!(u, 0.0);
                self.points[0]
            } else {
                vec::normalize(vec::add(
                    vec::mul(self.points[0], ((1.0 - u) * self.angle).sin()),
                    vec::mul(self.points[1], (u * self.angle).sin()),
                ))
            };
        } else if u <= 0.0 {
            debug_assert_eq!(u, 0.0);
            return vec::normalize(vec::add(
                vec::mul(self.points[0], ((1.0 - v) * self.angle).sin()),
                vec::mul(self.points[2], (v * self.angle).sin()),
            ));
        }

        // Adapted from "Exact Equations for Fuller’s Map Projection and Inverse"
        // https://utppublishing.com/doi/pdf/10.3138/carto.43.1.67
        // "Inverse with Spherical Linear Interpolation"

        // First, solve the equation:
        // `tan(x + proj_u) + tan(x) + tan(x + proj_v) + t_alpha = 0`
        let t_alpha = self.tan_half_angle;
        let proj_u = (2.0 * u + v - 1.0) * self.angle;
        let proj_v = (2.0 * v + u - 1.0) * self.angle;
        let t_half_proj_u = (proj_u / 2.0).tan();
        let t_half_proj_v = (proj_v / 2.0).tan();
        let den_u = 1.0 - t_half_proj_u * t_half_proj_u;
        let den_v = 1.0 - t_half_proj_v * t_half_proj_v;
        let s = t_half_proj_u * den_v + t_half_proj_v * den_u;
        let p = t_half_proj_u * t_half_proj_v;
        let (t_x, _) = solve_cubic_fuller(
            4.0 * p,
            4.0 * t_alpha * p - 4.0 * s,
            3.0 * den_u * den_v - 8.0 * p - 2.0 * t_alpha * s,
            t_alpha * den_u * den_v + 2.0 * s,
        )
        .into_iter()
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
            let t_x_proj_u = (den_u + 2.0 * t_half_proj_u) / (den_u - 2.0 * t_x * t_half_proj_u);
            let t_x_proj_v = (den_v + 2.0 * t_half_proj_v) / (den_v - 2.0 * t_x * t_half_proj_v);
            (t_x, t_x.abs().max(t_x_proj_u.abs()).max(t_x_proj_v.abs()))
        })
        .min_by(|(_, a), (_, b)| a.abs().partial_cmp(&b.abs()).unwrap())
        .unwrap();

        // Use the solution to compute angle offsets of two planes from the second and
        // third edges of the triangle. Construct planes by offsetting the second and third edges
        // of the triangle
        let s_alpha = self.sin_half_angle;
        let c_alpha = self.cos_half_angle;
        let s_k_proj_u = 2.0 * t_half_proj_u;
        let c_k_proj_u = 1.0 - t_half_proj_u * t_half_proj_u;
        let n_0 = vec::cross(
            mat::apply(
                [self.points[0], self.points[1]],
                [s_alpha + c_alpha * t_x, s_alpha - c_alpha * t_x],
            ),
            vec::sub(self.points[2], self.points[1]),
        );
        let n_1 = vec::cross(
            mat::apply(
                [self.points[1], self.points[2]],
                [
                    s_alpha * (c_k_proj_u - t_x * s_k_proj_u)
                        + c_alpha * (t_x * c_k_proj_u + s_k_proj_u),
                    s_alpha * (c_k_proj_u - t_x * s_k_proj_u)
                        - c_alpha * (t_x * c_k_proj_u + s_k_proj_u),
                ],
            ),
            vec::sub(self.points[0], self.points[2]),
        );

        // The projected point is the intersection of the two planes
        vec::normalize(vec::cross(n_0, n_1))
    }

    fn from_sphere(&self, point: [f64; 3]) -> [f64; 2] {
        // Adapted from "Exact Equations for Fuller’s Map Projection and Inverse"
        // https://utppublishing.com/doi/pdf/10.3138/carto.43.1.67
        // "Transformation Equations with Spherical Linear Interpolation"

        let d_0 = mat::det_3([point, self.points[0], self.points[1]]);
        let d_1 = mat::det_3([point, self.points[1], self.points[2]]);
        let d_2 = mat::det_3([point, self.points[2], self.points[0]]);
        let c = vec::dot(self.points[0], self.points[1]);
        let a_0 = (2.0 * (d_0 * c) / (d_0 * c + d_1 + d_2)).atan();
        let a_1 = (2.0 * (d_1 * c) / (d_1 * c + d_2 + d_0)).atan();
        let a_2 = (2.0 * (d_2 * c) / (d_2 * c + d_0 + d_1)).atan();
        let u = (self.angle + 2.0 * a_2 - a_0 - a_1) / (3.0 * self.angle);
        let v = (self.angle + 2.0 * a_0 - a_1 - a_2) / (3.0 * self.angle);
        [u, v]
    }

    type Transform = Transform<Self>;
    fn transform(self, offset: [f64; 2], linear: [[f64; 2]; 2]) -> Transform<Self> {
        Transform::new(self, offset, linear)
    }
}

#[test]
fn test_roundtrip() {
    use crate::basetri::BaseTriSphere;
    use crate::prelude::*;
    const N: usize = 100;
    for v in 0..N {
        for u in 0..(N - v) {
            let u = u as f64 / N as f64;
            let v = v as f64 / N as f64;
            let proj = Fuller.inside(BaseTriSphere::Icosa.face(0).side(0));
            let point = proj.to_sphere([u, v]);
            let [a_u, a_v] = proj.from_sphere(point);
            assert!(
                (a_u - u).abs() < 1e-12 && (a_v - v).abs() < 1e-12,
                "failed roundtrip: got {:?}, expected {:?}",
                [a_u, a_v],
                [u, v]
            );
        }
    }
}

/// Determines the real roots of a cubic polynomial of the form `a x³ + b x² + c x + d` as used in the Fuller
/// projection. This is not a general purpose solver and is optimized for the specific types of cubics encountered in
/// the Fuller projection. As such, multiple branches that are never encountered in that context have been removed for
/// performance purposes. Most notably, this method only handles polynomials with 2 or 3 real roots. Takes 36.9 FLOPs on
/// average in the `test` build.
fn solve_cubic_fuller(a: f64, b: f64, c: f64, d: f64) -> impl IntoIterator<Item = f64> {
    const SMALL: f64 = 1.0e-6;
    if a.abs() < SMALL {
        // Branch needed to avoid loss of precision when `a` is relatively small. Determines the real roots of the
        // quadratic `b x² + c x + d`. In the Fuller projection, this will always have two real roots, so we avoid
        // checking if the discriminant is non-negative in non-debug builds.
        // Branch taken ~10% of the time in testing.
        let j = 2.0 * d;
        let k = 2.0 * b;
        let disc = c * c - j * k;

        // The polynomial discriminant will never be negative so we don't need to branch here.
        // TODO: Can we find a proof?
        debug_assert!(disc >= 0.0, "Quadratic discriminant was negative.");
        let u = -c - disc.sqrt().copysign(c);
        [j / u, u / k, 0.0].into_iter().take(2)
    } else {
        // Convert to a depressed cubic of the form `t³ p t + q` where `t = x - b/3a` and finds the real roots.
        // https://en.wikipedia.org/wiki/Cubic_equation#Depressed_cubic
        // https://mathworld.wolfram.com/CubicFormula.html
        // In a general solver, a case for when `d` is close to 0 could improve performance. In testing this library,
        // `d` was never smaller than 10^-6 so the branch would not be useful.
        let b2 = b * b;
        let a2 = a * a;
        let ac = a * c;

        let p = (3.0 * ac - b2) / (3.0 * a2);
        // If `b` is close to 0, this expression will be very close to `d/a`. In testing this library `|b|` was less
        // than 10^-6 in less than 0.0001% of cases so the branch would cost more performance than it saved.
        let q = (2.0 * b2 * b - 9.0 * ac * b + 27.0 * a2 * d) / (27.0 * a2 * a);

        let s = b / (3.0 * a);

        // Solve depressed cubic
        // In general, a case to just output the negative cube root of `q` when `p` is close to 0 would improve performance.
        // However, in testing, `|p|` was less than 10^-6 in less than 0.0001% of cases, so branching would cost more than
        // it gained.

        // The discriminant will always be negative so we will always have three real roots and don't need to branch here.
        // TODO: Can we find a proof?
        debug_assert!(
            4.0 * p * p * p + 27.0 * q * q < 0.0,
            "Depressed cubic discriminant was not negative."
        );
        let p_3 = -p / 3.0;
        let k = 2.0 * p_3.sqrt();
        // Absolute value not necessary after pulling p_3*p_3 out of sqrt because p is always negative.
        let theta = (-q / (p_3 * k)).acos() / 3.0;
        const PHI: f64 = 2.0 * FRAC_PI_3;

        [
            k * theta.cos() - s,
            k * (theta + PHI).cos() - s,
            k * (theta + 2.0 * PHI).cos() - s,
        ]
        .into_iter()
        .take(3)
    }
}

#[test]
fn test_solve_cubic() {
    // Cubic cases
    assert_similar(
        solve_cubic_fuller(1.0, 4.0, 0.0, -5.0)
            .into_iter()
            .collect(),
        [
            -(5.0f64.sqrt() + 5.0) / 2.0,
            (5.0f64.sqrt() - 5.0) / 2.0,
            1.0,
        ],
    );
    assert_similar(
        solve_cubic_fuller(1.0, -2.0, 0.0, 1.0)
            .into_iter()
            .collect(),
        [
            (1.0 - 5.0f64.sqrt()) / 2.0,
            1.0,
            (1.0 + 5.0f64.sqrt()) / 2.0,
        ],
    );

    // Quadratic cases
    assert_similar(
        solve_cubic_fuller(0.0, 1.0, 0.0, -1.0)
            .into_iter()
            .collect(),
        [-1.0, 1.0],
    );
}

#[cfg(test)]
fn assert_similar<const N: usize>(mut a: Vec<f64>, mut b: [f64; N]) {
    const SMALL: f64 = 1.0e-6;
    a.sort_by(|u, v| u.partial_cmp(v).unwrap());
    b.sort_by(|u, v| u.partial_cmp(v).unwrap());
    if a.len() != b.len() || a.iter().zip(b.iter()).any(|(x, y)| (x - y).abs() > SMALL) {
        panic!("vectors are not similar, a = {:?}, b = {:?}", a, b);
    }
}
