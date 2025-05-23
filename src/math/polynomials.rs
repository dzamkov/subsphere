use std::f64::consts::FRAC_PI_3;
use crate::math::SMALL;

/// Represents the distinct real roots of a cubic equation.
#[derive(Copy, Clone)]
pub enum Roots {
    Zero,
    One(f64),
    Two(f64, f64),
    Three(f64, f64, f64)
}

impl Roots {
    /// Appends `v` to the end of the list of roots if there are less than three, otherwise does nothing.
    #[inline]
    const fn append(self, v: f64) -> Self {
        match self {
            Roots::Zero => Self::One(v),
            Roots::One(a) => Self::Two(a, v),
            Roots::Two(a, b) => Self::Three(a, b, v),
            Roots::Three(a, b, c) => Self::Three(a, b, c)
        }
    }
    
    /// Adds `v` to the value of all roots.
    #[inline]
    const fn add(self, v: f64) -> Self {
        match self {
            Roots::Zero => Roots::Zero,
            Roots::One(a) => Roots::One(a + v),
            Roots::Two(a, b) => Roots::Two(a + v, b + v),
            Roots::Three(a, b, c) => Roots::Three(a + v, b + v, c + v)
        }
    }
}

impl IntoIterator for Roots {
    type Item = f64;
    type IntoIter = std::iter::Take<std::array::IntoIter<f64, 3>>;
    
    fn into_iter(self) -> Self::IntoIter {
        match self {
            Roots::Zero => [0.0; 3].into_iter().take(0),
            Roots::One(a) => [a, 0.0, 0.0].into_iter().take(1),
            Roots::Two(a, b) => [a, b, 0.0].into_iter().take(2),
            Roots::Three(a, b, c) => [a, b, c].into_iter().take(3)
        }
    }
}

/// Solves the linear system `ax+b=0`.
const fn solve_linear(a: f64, b: f64) -> Roots {
    if a.abs() <= SMALL {
        Roots::Zero
    } else {
        Roots::One(-b / a)
    }
}

/// Determines the real roots of a quadratic polynomial of the form `a x² + b x + c`. If a root has multiplicity greater
/// than one, it will only be including once in the output.
fn solve_quadratic(a: f64, b: f64, c: f64) -> Roots {
    if a.abs() < SMALL {
        solve_linear(b, c)
    } else if b.abs() < SMALL {
        if c > 0.0 {
            Roots::Zero
        } else {
            Roots::Two(c.abs().sqrt(), -c.abs().sqrt())
        }
    } else {
        let disc = b * b - 4.0 * a * c;
        
        if disc.abs() < SMALL {
            Roots::One(-b / (2.0 * a))
        } else if disc < 0.0 {
            Roots::Zero
        } else {
            let u = -b - b.signum() * disc.sqrt();
            Roots::Two(2.0 * c / u, u / (2.0 * a))
        }
    }
}

// https://en.wikipedia.org/wiki/Cubic_equation#Depressed_cubic
// https://mathworld.wolfram.com/CubicFormula.html
/// Determines the real roots of a depressed cubic polynomial of the form `x³ + p x + q`. If a root has multiplicity
/// greater than one, it will only appear once in the output.
fn solve_depressed_cubic(p: f64, q: f64) -> Roots {
    if p.abs() < SMALL {
        Roots::One(-q.cbrt())
    } else {
        let p3 = p * p * p;
        let q2 = q * q;
        let disc = 4.0 * p3 + 27.0 * q2;
        
        if disc.abs() < SMALL {
            // Two equal roots
            let r = 3.0 * q / p;
            Roots::Two(r, -r / 2.0)
        } else if disc < 0.0 {
            // Three distinct real roots
            // 2a, 7m, 3d, 1sqrt, 1acos, 3cos
            let p_3 = -p / 3.0;
            let k = 2.0 * p_3.sqrt();
            // Absolute value not necessary after pulling p_3*p_3 out of sqrt because p is always negative.
            let theta = (-q / (p_3 * k)).acos() / 3.0;
            let phi = 2.0 * FRAC_PI_3;
            
            Roots::Three(
                k * theta.cos(),
                k * (theta + phi).cos(),
                k * (theta + 2.0 * phi).cos()
            )
        } else {
            // One real root
            let q_2 = -q / 2.0;
            let b = (q2 / 4.0 + p3 / 27.0).sqrt();
            let u1 = q_2 + b;
            let u2 = q_2 - b;
            Roots::One(u1.cbrt() + u2.cbrt())
        }
    }
}

/// Determines the real roots of a cubic polynomial of the form `a x³ + b x² + c x + d`. IF a root has multiplicity
/// greater than one, it will only appear once in the output.
pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> Roots {
    if a.abs() < SMALL {
        solve_quadratic(b, c, d)
    } else if d.abs() < SMALL {
        solve_quadratic(a, b, c).append(0.0)
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
        
        solve_depressed_cubic(p, q).add(-s)
    }
}

#[test]
fn test_solve_cubic() {
    // Cubic cases
    assert_similar(solve_cubic(1.0, -1.0, 1.0, -1.0).into_iter().collect(), [1.0]);
    assert_similar(solve_cubic(1.0, 0.0, 0.0, -27.0).into_iter().collect(), [3.0]);
    assert_similar(solve_cubic(8.0, 8.0, 0.0, -3.0).into_iter().collect(), [0.5]);
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
    assert_similar(solve_cubic(0.0, 1.0, 0.0, 1.0).into_iter().collect(), []);

    // Linear cases
    assert_similar(solve_cubic(0.0, 0.0, 1.0, -1.0).into_iter().collect(), [1.0]);
}

#[cfg(test)]
fn assert_similar<const N: usize>(mut a: Vec<f64>, mut b: [f64; N]) {
    a.sort_by(|u, v| u.partial_cmp(v).unwrap());
    b.sort_by(|u, v| u.partial_cmp(v).unwrap());
    if a.len() != b.len() || a.iter().zip(b.iter()).any(|(x, y)| (x - y).abs() > SMALL) {
        panic!("vectors are not similar, a = {:?}, b = {:?}", a, b);
    }
}
