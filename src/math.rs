//! This module contains a minimal set of linear algebra types and functions used by this crate. It
//! has a very specific grab bag of functionality and is not intended to be publicly exposed for
//! general use.
pub use std::f64::consts::PI;

/// Contains functions related to vectors.
pub(crate) mod vec {
    /// Adds `a` to `b`.
    #[inline]
    pub fn add<const N: usize>(a: [f64; N], b: [f64; N]) -> [f64; N] {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = a[i] + b[i];
        }
        res
    }

    /// Subtracts `b` from `a`.
    #[inline]
    pub fn sub<const N: usize>(a: [f64; N], b: [f64; N]) -> [f64; N] {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = a[i] - b[i];
        }
        res
    }

    /// Multiplies a vector by a scalar.
    #[inline]
    pub fn mul<const N: usize>(a: [f64; N], b: f64) -> [f64; N] {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = a[i] * b;
        }
        res
    }

    /// Divides a vector by a scalar.
    #[inline]
    pub fn div<const N: usize>(a: [f64; N], b: f64) -> [f64; N] {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = a[i] / b;
        }
        res
    }

    /// Computes the dot product of two vectors.
    #[inline]
    pub fn dot<const N: usize>(a: [f64; N], b: [f64; N]) -> f64 {
        let mut res = 0.0;
        for i in 0..N {
            res += a[i] * b[i];
        }
        res
    }

    /// Computes the dot product of two vectors in a basis where the dot product between two
    /// basis vectors is `basis_dot`.
    #[inline]
    #[expect(clippy::needless_range_loop)]
    pub fn dot_in<const N: usize>(basis_dot: f64, a: [f64; N], b: [f64; N]) -> f64 {
        let mut res = 0.0;
        for i in 0..N {
            res += a[i] * b[i];
        }
        for i in 0..N {
            for j in 0..N {
                if i != j {
                    res += basis_dot * a[i] * b[j];
                }
            }
        }
        res
    }

    /// Computes the cross product of two vectors.
    #[inline]
    pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    /// Normalizes a vector.
    #[inline]
    pub fn normalize<const N: usize>(a: [f64; N]) -> [f64; N] {
        div(a, dot(a, a).sqrt())
    }

    /// Normalizes a vector in a basis where the dot product between two basis vectors is
    /// `basis_dot`.
    #[inline]
    pub fn normalize_in<const N: usize>(basis_dot: f64, a: [f64; N]) -> [f64; N] {
        div(a, dot_in(basis_dot, a, a).sqrt())
    }
}

/// Contains functions related to matrices.
pub(crate) mod mat {
    use super::*;

    /// Multiplies a matrix with a vector.
    #[inline]
    pub fn apply<const N: usize, const M: usize>(mat: [[f64; N]; M], vec: [f64; M]) -> [f64; N] {
        let mut res = [0.0; N];
        for i in 0..M {
            res = vec::add(res, vec::mul(mat[i], vec[i]));
        }
        res
    }

    /// Computes the determinant of a matrix.
    #[inline]
    pub fn det(m: [[f64; 3]; 3]) -> f64 {
        vec::dot(m[0], vec::cross(m[1], m[2]))
    }
}

/// Computes the area (or equivalently, the solid angle) of a spherical triangle.
pub fn sphere_tri_area(points: [[f64; 3]; 3]) -> f64 {
    // https://www.johndcook.com/blog/2021/11/29/area-of-spherical-triangle/
    let d = 1.0
        + vec::dot(points[0], points[1])
        + vec::dot(points[1], points[2])
        + vec::dot(points[2], points[0]);
    (mat::det(points) / d).atan() * 2.0
}

/// Determines the real roots of a cubic polynomial of the form `a x³ + b x² + c x + d`.
pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> impl Iterator<Item = f64> {
    // TODO: This is an extremely hacky way of handling the quadratic case. Replace with a
    // proper numerically-stable implementation.
    let a = if a == 0.0 {
        1.0e-6
    } else {
        a
    };

    // See https://mathworld.wolfram.com/CubicFormula.html for derivation
    let b_i_3a = (b / 3.0) / a;
    let c_i_a = c / a;
    let d_i_a = d / a;

    // Let `x = t - b / 3 a`. Then the equation becomes `t³ + 3 q t - 2 r = 0` where:
    let q = c_i_a / 3.0 - b_i_3a * b_i_3a;
    let r = (b_i_3a * c_i_a - d_i_a) / 2.0 - b_i_3a * b_i_3a * b_i_3a;
    let disc = q * q * q + r * r;
    if disc >= 0.0 {
        // Equation has one real
        let w = (r + disc.sqrt()).cbrt();
        let t = w - q / w;
        let mut iter = [0.0, 0.0, t - b_i_3a].into_iter();
        iter.next().unwrap();
        iter.next().unwrap();
        iter
    } else {
        // Equation has three real roots
        let h = (r * r - disc).sqrt();
        let s_r = h.cbrt();
        let s_theta_0 = (r / h).acos() / 3.0;
        let t_0 = 2.0 * s_r * (s_theta_0 + 2.0 * PI / 3.0).cos();
        let x_0 = t_0 - b_i_3a;
        let t_1 = 2.0 * s_r * (s_theta_0 - 2.0 * PI / 3.0).cos();
        let x_1 = t_1 - b_i_3a;
        let t_2 = 2.0 * s_r * s_theta_0.cos();
        let x_2 = t_2 - b_i_3a;
        [x_0, x_1, x_2].into_iter()
    }
}

#[test]
fn test_solve_cubic() {
    assert_similar(solve_cubic(1.0, -1.0, 1.0, -1.0).collect(), [1.0]);
    assert_similar(solve_cubic(1.0, 0.0, 0.0, -27.0).collect(), [3.0]);
    assert_similar(solve_cubic(8.0, 8.0, 0.0, -3.0).collect(), [0.5]);
    assert_similar(
        solve_cubic(1.0, 4.0, 0.0, -5.0).collect(),
        [
            -(5.0f64.sqrt() + 5.0) / 2.0,
            (5.0f64.sqrt() - 5.0) / 2.0,
            1.0,
        ],
    );
    assert_similar(
        solve_cubic(1.0, -2.0, 0.0, 1.0).collect(),
        [
            (1.0 - 5.0f64.sqrt()) / 2.0,
            1.0,
            (1.0 + 5.0f64.sqrt()) / 2.0,
        ],
    );
}

#[cfg(test)]
fn assert_similar<const N: usize>(a: Vec<f64>, b: [f64; N]) {
    if a.len() != b.len() || a.iter().zip(b.iter()).any(|(x, y)| (x - y).abs() > 1e-6) {
        panic!("vectors are not similar, a = {:?}, b = {:?}", a, b);
    }
}
