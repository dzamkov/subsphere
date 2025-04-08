//! This module contains a minimal set of linear algebra types and functions used by this crate.

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
        mul(a, 1.0 / dot(a, a).sqrt())
    }
}

/// Contains functions related to matrices.
pub(crate) mod mat {
    use super::*;

    /// Multiplies a matrix with a vector.
    #[inline]
    pub fn apply<const N: usize, const M: usize>(
        mat: [[f64; N]; M],
        vec: [f64; M],
    ) -> [f64; N] {
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
