//! This module contains a minimal set of linear algebra types and functions used by this crate.

/// A scalar value.
pub type Scalar = f64;

/// A three-dimensional vector of [`Scalar`] values.
pub type Vector3 = [Scalar; 3];

/// A 3x3 matrix of [`Scalar`] values.
pub(crate) type Matrix3 = [Vector3; 3];

/// Contains functions related to [`Vector3`].
pub(crate) mod vec3 {
    use super::*;

    /// Adds `a` to `b`.
    #[inline]
    pub const fn add(a: Vector3, b: Vector3) -> Vector3 {
        [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
    }

    /// Multiplies a vector by a [`Scalar`].
    #[inline]
    pub const fn mul(a: Vector3, b: Scalar) -> Vector3 {
        [a[0] * b, a[1] * b, a[2] * b]
    }

    /// Computes the dot product of two vectors.
    #[inline]
    pub const fn dot(a: Vector3, b: Vector3) -> Scalar {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }

    /// Computes the cross product of two vectors.
    #[inline]
    pub const fn cross(a: Vector3, b: Vector3) -> Vector3 {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ]
    }
}

/// Contains functions related to [`Matrix3`].
pub(crate) mod mat3 {
    use super::*;

    /// Multiplies a matrix with a vector.
    #[inline]
    pub const fn apply(m: Matrix3, v: Vector3) -> Vector3 {
        let x = vec3::mul(m[0], v[0]);
        let y = vec3::mul(m[1], v[1]);
        let z = vec3::mul(m[2], v[2]);
        vec3::add(vec3::add(x, y), z)
    }

    /// Computes the determinant of a matrix.
    #[inline]
    pub fn det(m: Matrix3) -> Scalar {
        vec3::dot(m[0], vec3::cross(m[1], m[2]))
    }
}
