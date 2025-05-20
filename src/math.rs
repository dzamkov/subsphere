//! This module contains a minimal set of linear algebra types and functions used by this crate. It
//! has a very specific grab bag of functionality and is not intended to be publicly exposed for
//! general use.

use std::f64::consts::FRAC_PI_3;

/// Contains functions related to vectors.
pub(crate) mod vec {
    /// Negates a vector.
    #[inline]
    pub const fn neg<const N: usize>(a: [f64; N]) -> [f64; N] {
        let mut res = [0.0; N];
        let mut i = 0;
        while i < N {
            res[i] = -a[i];
            i += 1;
        }
        res
    }

    /// Adds `a` to `b`.
    #[inline]
    pub const fn add<const N: usize>(a: [f64; N], b: [f64; N]) -> [f64; N] {
        let mut res = [0.0; N];
        let mut i = 0;
        while i < N {
            res[i] = a[i] + b[i];
            i += 1;
        }
        res
    }

    /// Subtracts `b` from `a`.
    #[inline]
    pub const fn sub<const N: usize>(a: [f64; N], b: [f64; N]) -> [f64; N] {
        let mut res = [0.0; N];
        let mut i = 0;
        while i < N {
            res[i] = a[i] - b[i];
            i += 1;
        }
        res
    }

    /// Multiplies a vector by a scalar.
    #[inline]
    pub const fn mul<const N: usize>(a: [f64; N], b: f64) -> [f64; N] {
        let mut res = [0.0; N];
        let mut i = 0;
        while i < N {
            res[i] = a[i] * b;
            i += 1;
        }
        res
    }

    /// Divides a vector by a scalar.
    #[inline]
    pub const fn div<const N: usize>(a: [f64; N], b: f64) -> [f64; N] {
        let mut res = [0.0; N];
        let mut i = 0;
        while i < N {
            res[i] = a[i] / b;
            i += 1;
        }
        res
    }

    /// Computes the dot product of two vectors.
    #[inline]
    pub const fn dot<const N: usize>(a: [f64; N], b: [f64; N]) -> f64 {
        let mut res = 0.0;
        let mut i = 0;
        while i < N {
            res += a[i] * b[i];
            i += 1;
        }
        res
    }

    /// Computes the cross product of two vectors.
    #[inline]
    pub const fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
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
}

/// Contains functions related to matrices.
pub(crate) mod mat {
    use super::*;

    /// Multiplies a matrix with a vector.
    #[inline]
    pub const fn apply<const N: usize, const M: usize>(mat: [[f64; N]; M], vec: [f64; M]) -> [f64; N] {
        let mut res = [0.0; N];
        let mut i = 0;
        while i < M {
            res = vec::add(res, vec::mul(mat[i], vec[i]));
            i += 1;
        }
        res
    }

    /// Computes the determinant of a 3x3 matrix.
    #[inline]
    pub const fn det_3(m: [[f64; 3]; 3]) -> f64 {
        vec::dot(m[0], vec::cross(m[1], m[2]))
    }

    /// Computes the determinant of a 2x2 matrix.
    #[inline]
    pub const fn det_2(m: [[f64; 2]; 2]) -> f64 {
        m[0][0] * m[1][1] - m[0][1] * m[1][0]
    }

    /// Computes the adjoint of a 2x2 matrix.
    #[inline]
    pub const fn adjoint_2(m: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
        let [[a, b], [c, d]] = m;
        [[d, -b], [-c, a]]
    }
}
const SMALL: f64 = 1.0e-6;

const fn solve_linear(a: f64, b: f64) -> Option<f64> {
    if a.abs() <= SMALL {
        None
    } else {
        Some(-b / a)
    }
}

fn solve_quadratic(a: f64, b: f64, c: f64) -> [Option<f64>; 2] {
    if a.abs() < SMALL {
        [ solve_linear(b, c), None ]
    } else {
        let disc = b * b - 4.0 * a * c;
        
        if disc.abs() < SMALL {
            [ Some(-b / (2.0 * a)), None ]
        } else if disc < 0.0 {
            [ None, None ]
        } else {
            if b.abs() < SMALL {
                [ Some(c.abs().sqrt()), Some(-c.abs().sqrt()) ]
            } else {
                let u = -b - b.signum() * disc.sqrt();
                [ Some(2.0 * c / u), Some(u / (2.0 * a)) ]
            }
        }
    }
}

// https://en.wikipedia.org/wiki/Cubic_equation#Depressed_cubic
// https://mathworld.wolfram.com/CubicFormula.html
fn solve_depressed_cubic(p: f64, q: f64) -> [Option<f64>; 3] {
    if p.abs() < SMALL {
        [ Some(-q.cbrt()), None, None ]
    } else {
        let p3 = p * p * p;
        let q2 = q * q;
        let disc = 4.0 * p3 + 27.0 * q2;
        
        if disc.abs() < SMALL {
            // Two equal roots
            let r = 3.0 * q / p;
            [ Some(r), Some(-r / 2.0), None ]
        } else {
            if disc < 0.0 {
                // Three distinct real roots
                // 2a, 7m, 3d, 1sqrt, 1acos, 3cos
                let p_3 = -p / 3.0;
                let k = 2.0 * p_3.sqrt();
                // Absolute value not necessary after pulling p_3*p_3 out of sqrt because p is always negative.
                let theta = (-q / (p_3 * k)).acos() / 3.0;
                let phi = 2.0 * FRAC_PI_3;
                
                [
                    Some(k * theta.cos()),
                    Some(k * (theta + phi).cos()),
                    Some(k * (theta + 2.0 * phi).cos())
                ]
            } else {
                // One real root
                let q_2 = -q / 2.0;
                let b = (q2 / 4.0 + p3 / 27.0).sqrt();
                let u1 = q_2 + b;
                let u2 = q_2 - b;
                [ Some(u1.cbrt() + u2.cbrt()), None, None ]
            }
        }
    }
}

/// Determines the real roots of a cubic polynomial of the form `a x³ + b x² + c x + d`. IF a root has multiplicity
/// greater than one, it will only appear once in the output.
pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> [Option<f64>; 3] {
    if a.abs() < SMALL {
        let roots = solve_quadratic(b, c, d);
        
        [ roots[0], roots[1], None ]
    } else if d.abs() < SMALL {
        let roots = solve_quadratic(a, b, c);
        
        [ Some(0.0), roots[0], roots[1] ]
    } else {
        let b2 = b * b;
        let a2 = a * a;
        let ac = a * c;
        
        let p = (3.0 * ac - b2) / (3.0 * a2);
        let q = if b.abs() < SMALL {
            d / a
        } else {
            (2.0 * b2 * b - 9.0 * ac * b + 27.0 * a2 * d) / (27.0 * a2 * a)
        };
        
        let s = b / (3.0 * a);
        
        let roots = solve_depressed_cubic(p, q);
        
        [ roots[0].map(|r| r - s), roots[1].map(|r| r - s), roots[2].map(|r| r - s)]
    }
}

#[test]
fn test_solve_cubic() {
    // Cubic cases
    assert_similar(solve_cubic(1.0, -1.0, 1.0, -1.0).into_iter().filter_map(|a| a).collect(), [1.0]);
    assert_similar(solve_cubic(1.0, 0.0, 0.0, -27.0).into_iter().filter_map(|a| a).collect(), [3.0]);
    assert_similar(solve_cubic(8.0, 8.0, 0.0, -3.0).into_iter().filter_map(|a| a).collect(), [0.5]);
    assert_similar(
        solve_cubic(1.0, 4.0, 0.0, -5.0).into_iter().filter_map(|a| a).collect(),
        [
            -(5.0f64.sqrt() + 5.0) / 2.0,
            (5.0f64.sqrt() - 5.0) / 2.0,
            1.0,
        ],
    );
    assert_similar(
        solve_cubic(1.0, -2.0, 0.0, 1.0).into_iter().filter_map(|a| a).collect(),
        [
            (1.0 - 5.0f64.sqrt()) / 2.0,
            1.0,
            (1.0 + 5.0f64.sqrt()) / 2.0,
        ],
    );

    // Quadratic cases
    assert_similar(solve_cubic(0.0, 1.0, 0.0, -1.0).into_iter().filter_map(|a| a).collect(), [-1.0, 1.0]);
    assert_similar(solve_cubic(0.0, 1.0, 0.0, 1.0).into_iter().filter_map(|a| a).collect(), []);

    // Linear cases
    assert_similar(solve_cubic(0.0, 0.0, 1.0, -1.0).into_iter().filter_map(|a| a).collect(), [1.0]);
}

#[cfg(test)]
fn assert_similar<const N: usize>(mut a: Vec<f64>, mut b: [f64; N]) {
    a.sort_by(|u, v| u.partial_cmp(v).unwrap());
    b.sort_by(|u, v| u.partial_cmp(v).unwrap());
    if a.len() != b.len() || a.iter().zip(b.iter()).any(|(x, y)| (x - y).abs() > SMALL) {
        panic!("vectors are not similar, a = {:?}, b = {:?}", a, b);
    }
}
