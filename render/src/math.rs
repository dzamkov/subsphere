//! This module contains a minimal set of linear algebra types and functions used by this crate. It
//! has a very specific grab bag of functionality and is not intended to be publicly exposed for
//! general use.

/// Contains functions related to vectors.
pub(crate) mod vec {
    /// Negates a vector.
    #[inline]
    pub fn neg<const N: usize>(a: [f64; N]) -> [f64; N] {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = -a[i];
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
    pub fn div<const N: usize>(a: [f64; N], b: f64) -> [f64; N] {
        let mut res = [0.0; N];
        for i in 0..N {
            res[i] = a[i] / b;
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
}
