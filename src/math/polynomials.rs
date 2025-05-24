use std::f64::consts::FRAC_PI_3;
use crate::math::SMALL;

/// Represents the distinct real roots of a cubic equation used in the Fuller projection. This library will only ever
/// solve equations with 2 or 3 distinct real roots so we neglect handling other cases in the name of performance.
#[derive(Copy, Clone)]
pub struct Roots(f64, f64, Option<f64>);

impl Roots {
    /// Adds `v` to the value of all roots.
    #[inline]
    const fn add(self, v: f64) -> Self {
        match self.2 {
            None => Self(self.0 + v, self.1 + v, None),
            Some(c) => Self(self.0 + v, self.1 + v, Some(c + v))
        }
    }
}

impl IntoIterator for Roots {
    type Item = f64;
    type IntoIter = std::iter::Take<std::array::IntoIter<f64, 3>>;
    
    fn into_iter(self) -> Self::IntoIter {
        match self.2 {
            None => [self.0, self.1, 0.0].into_iter().take(2),
            Some(c) => [self.0, self.1, c].into_iter().take(3)
        }
    }
}

/// Determines the real roots of a quadratic polynomial of the form `a x² + b x + c`. If a root has multiplicity greater
/// than one, it will only be included once in the output. Will output `NaN` roots if the discriminant is negative, but
/// this did not occur in testing and is thus ignored. Takes 9 FLOPs.
fn solve_quadratic(a: f64, b: f64, c: f64) -> Roots {
    let j = 2.0 * c;
    let k = 2.0 * a;
    let disc = b * b - j * k;
    
    // The polynomial discriminant will never be negative so we don't need to branch here.
    // TODO: Can we find a proof?
    let u = -b - disc.sqrt().copysign(b);
    Roots(j / u, u / k, None)
}

// https://en.wikipedia.org/wiki/Cubic_equation#Depressed_cubic
// https://mathworld.wolfram.com/CubicFormula.html
/// Determines the real roots of a depressed cubic polynomial of the form `x³ + p x + q`. If a root has multiplicity
/// greater than one, it will only appear once in the output. Takes 17 FLOPs.
fn solve_depressed_cubic(p: f64, q: f64) -> Roots {
    // In general, a case to just output the negative cube root of `q` when `p` is close to 0 would improve performance.
    // However, in testing, `|p|` was less than 10^-6 in less than 0.0001% of cases, so branching would cost more than
    // it gained.
    
    // The discriminant will always be negative so we will always have three real roots and don't need to branch here.
    // TODO: Can we find a proof?
    let p_3 = -p / 3.0;
    let k = 2.0 * p_3.sqrt();
    // Absolute value not necessary after pulling p_3*p_3 out of sqrt because p is always negative.
    let theta = (-q / (p_3 * k)).acos() / 3.0;
    const PHI: f64 = 2.0 * FRAC_PI_3;
    
    Roots(
        k * theta.cos(),
        k * (theta + PHI).cos(),
        Some(k * (theta + 2.0 * PHI).cos())
    )
}

/// Determines the real roots of a cubic polynomial of the form `a x³ + b x² + c x + d`. If a root has multiplicity
/// greater than one, it will only appear once in the output. Takes 36.9 FLOPs in this library's typical usage.
/// This is optimized for this library's use case of solving equations with 2 or 3 real roots and will fail when used
/// in other cases.
pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> Roots {
    if a.abs() < SMALL {
        // Branch taken ~10% of the time in testing.
        solve_quadratic(b, c, d)
    } else {
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
        
        solve_depressed_cubic(p, q).add(-s)
    }
}

#[test]
fn test_solve_cubic() {
    // Cubic cases
    assert_similar(
        solve_cubic(1.0, 4.0, 0.0, -5.0).into_iter().collect(),
        [
            -(5.0f64.sqrt() + 5.0) / 2.0,
            (5.0f64.sqrt() - 5.0) / 2.0,
            1.0,
        ],
    );
    assert_similar(
        solve_cubic(1.0, -2.0, 0.0, 1.0).into_iter().collect(),
        [
            (1.0 - 5.0f64.sqrt()) / 2.0,
            1.0,
            (1.0 + 5.0f64.sqrt()) / 2.0,
        ],
    );

    // Quadratic cases
    assert_similar(solve_cubic(0.0, 1.0, 0.0, -1.0).into_iter().collect(), [-1.0, 1.0]);
}

#[cfg(test)]
fn assert_similar<const N: usize>(mut a: Vec<f64>, mut b: [f64; N]) {
    a.sort_by(|u, v| u.partial_cmp(v).unwrap());
    b.sort_by(|u, v| u.partial_cmp(v).unwrap());
    if a.len() != b.len() || a.iter().zip(b.iter()).any(|(x, y)| (x - y).abs() > SMALL) {
        panic!("vectors are not similar, a = {:?}, b = {:?}", a, b);
    }
}
