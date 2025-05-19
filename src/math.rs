//! This module contains a minimal set of linear algebra types and functions used by this crate. It
//! has a very specific grab bag of functionality and is not intended to be publicly exposed for
//! general use.
pub use std::f64::consts::PI;

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
        let c2 = c * c;
        let ac = a * c;
        let bd = b * d;
        let disc = b2 * c2 - 4.0 * ac * c2 - 4.0 * b2 * bd - 27.0 * a2 * d * d + 18.0 * ac * bd;
        
        // https://en.wikipedia.org/wiki/Cubic_equation#Depressed_cubic
        if disc.abs() < SMALL {
            // Two roots are equal
            if (b2 - 3.0 * ac).abs() < SMALL {
                // Triple root
                [ Some(-b / (3.0 * a)), None, None ]
            } else {
                // Double root
                [
                    Some((4.0 * ac * b - 9.0 * a2 * d - b2 * b) / (a * (b2 - 3.0 * ac))),
                    Some((9.0 * a * d - b * c) / (2.0 * (b2 - 3.0 * ac))),
                    None
                ]
            }
        } else {
            let p = (3.0 * ac - b2) / (3.0 * a2);
            let q = (2.0 * b2 * b - 9.0 * ac * b + 27.0 * a2 * d) / (27.0 * a2 * a);
            let s = b / (3.0 * a);
            
            if p.abs() < SMALL {
                [ Some((-q).cbrt() - s), None, None ]
            } else {
                if disc > 0.0 {
                    // Three distinct real roots
                    let x = 3.0 * q / (2.0 * p) * (-3.0 / p).sqrt();
                    let acos = x.acos() / 3.0;
                    let k = 2.0 * (-p / 3.0).sqrt();
                    let t = [
                        // Unfortunately, it's not possible to reduce this except by solving a second cubic.
                        k * acos.cos(),
                        k * (acos - 2.0 * PI / 3.0).cos(),
                        k * (acos - 4.0 * PI / 3.0).cos(),
                    ];
                    
                    [ Some(t[0] - s), Some(t[1] - s), Some(t[2] - s) ]
                } else {
                    // One real root, two conjugate complex root
                    if p > 0.0 {
                        [
                            Some(-2.0 * (p / 3.0).sqrt() * (((3.0 * q) / (2.0 * p) * (3.0 / p).sqrt()).asinh() / 3.0).sinh() - s),
                            None,
                            None
                        ]
                    } else {
                        [
                            Some(-2.0 * q.signum() * (-p / 3.0).sqrt() * ((-3.0 * q.abs() / (2.0 * p) * (-3.0 / p).sqrt()).acosh() / 3.0).cosh() - s),
                            None,
                            None
                        ]
                    }
                }
            }
        }
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
