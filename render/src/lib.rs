mod math;

use math::{mat, vec};
use subsphere::prelude::*;

/// Renders a [`Scene`] from the perspective of a [`Camera`] to the given target image.
pub fn render<Sphere: subsphere::Sphere>(
    scene: &Scene<Sphere>,
    camera: &Camera,
    target: &mut image::RgbImage,
) {
    // Raytrace the sphere.
    let size_x = target.width();
    let size_y = target.height();
    for (x, y, pixel) in target.enumerate_pixels_mut() {
        let mut total = [0.0, 0.0, 0.0];
        for [offset_x, offset_y] in SAMPLES {
            let view_x = ((x as f64 + offset_x) / size_x as f64) * 2.0 - 1.0;
            let view_y = 1.0 - ((y as f64 + offset_y) / size_y as f64) * 2.0;
            let (start, dir) = camera.unproject([view_x, view_y]);
            if let Some(t) = trace_sphere(start, dir) {
                let pos = vec::add(start, vec::mul(dir, t));
                let face = scene.sphere.face_at(pos);
                let light =
                    scene.ambient + vec::dot(pos, scene.light_dir).max(0.0) * (1.0 - scene.ambient);
                let color = scene.faces[face.index()];
                total = vec::add(total, vec::mul(color, light));
            } else {
                total = vec::add(total, scene.background);
            }
        }
        *pixel = to_srgb(vec::div(total, SAMPLES.len() as f64));
    }
}

/// The per-pixel sampling positions used for anti-aliasing.
const SAMPLES: &[[f64; 2]] = &[
    [0.375, 0.125],
    [0.875, 0.375],
    [0.125, 0.625],
    [0.625, 0.875],
];

/// Describes a scene to be rendered.
pub struct Scene<'a, Sphere> {
    /// The background color for the scene.
    pub background: Color,

    /// The direction *towards* the light source.
    pub light_dir: [f64; 3],

    /// The amount of ambient light in the scene, between 0.0 and 1.0.
    pub ambient: f64,

    /// The sphere to be rendered.
    ///
    /// This is placed at the origin of the scene.
    pub sphere: &'a Sphere,

    /// The colors of `sphere`'s faces.
    pub faces: Box<[Color]>,
}

/// Describes a perspective a [`Scene`] can be rendered from.
pub struct Camera {
    /// The position of the camera in the scene.
    pos: [f64; 3],

    /// The forward direction of the camera.
    forward: [f64; 3],

    /// The extent of the camera's view plane, along each view axis.
    extent: [[f64; 3]; 2],

    /// The amount that the camera rays diverge from being perpendicular to the view plane.
    ///
    /// This will be zero for othographic cameras, and positive for perspective cameras.
    divergence: f64,
}

impl Camera {
    /// Constructs an orthographic camera which looks at the sphere from the given position.
    ///
    /// The sphere will always be in the center of the view, with the Z axis pointing upwards,
    /// along the positive view Y direction.
    ///
    /// `extent` specifies the scale factor which converts from view coordinates (in the range
    /// [-1.0, 1.0]) to world coordinates along each view axis.
    pub fn ortho(pos: [f64; 3], extent: [f64; 2]) -> Self {
        Self::perspective(pos, extent, 0.0)
    }

    /// Constructs a perspective camera which looks at the sphere from the given position.
    ///
    /// The sphere will always be in the center of the view, with the Z axis pointing upwards,
    /// along the positive view Y direction.
    pub fn perspective(pos: [f64; 3], extent: [f64; 2], divergence: f64) -> Self {
        let up = [0.0, 0.0, 1.0];
        let forward = vec::normalize(vec::neg(pos));
        let extent_x = vec::mul(vec::normalize(vec::cross(forward, up)), extent[0]);
        let extent_y = vec::mul(vec::normalize(vec::cross(extent_x, forward)), extent[1]);
        Self {
            pos,
            forward,
            extent: [extent_x, extent_y],
            divergence,
        }
    }

    /// Gets the ray which corresponds to the given pixel in view coordinates.
    fn unproject(&self, coords: [f64; 2]) -> ([f64; 3], [f64; 3]) {
        let offset = mat::apply(self.extent, coords);
        let start = vec::add(self.pos, offset);
        let dir = vec::normalize(vec::add(self.forward, vec::mul(offset, self.divergence)));
        (start, dir)
    }
}

/// The color type used for rendering.
///
/// Colors are defined in the linear sRGB color space, with components in the range [0.0, 1.0].
pub type Color = [f64; 3];

/// Assigns each face of the given sphere an arbitrary color.
pub fn colorize(sphere: &impl subsphere::Sphere) -> Box<[Color]> {
    use rand::seq::SliceRandom;
    use rand::{Rng, SeedableRng};
    let mut face_colors: Box<[Option<u32>]> = vec![None; sphere.num_faces()].into_boxed_slice();

    // Assign each face a color that is different from its neighbors
    let mut rng = rand::rngs::SmallRng::seed_from_u64(1);
    let mut faces = sphere.faces().collect::<Vec<_>>();
    faces.shuffle(&mut rng);
    let mut num_colors = 6;
    'retry: loop {
        for face in faces.iter() {
            let index = face.index();
            let mut available_colors = (1u64 << num_colors) - 1;
            for side in face.sides() {
                let neighbor = side.twin().inside();
                if let Some(color) = face_colors[neighbor.index()] {
                    available_colors &= !(1 << color);
                }
            }
            if available_colors == 0 {
                // No colors available, try again with more colors
                num_colors += 1;
                continue 'retry;
            }
            let mut select = rng.random_range(0..available_colors.count_ones());
            let mut color = available_colors.trailing_zeros();
            while select > 0 {
                available_colors &= !(1 << color);
                select -= 1;
                color = available_colors.trailing_zeros();
            }
            face_colors[index] = Some(color);
        }
        break;
    }

    // Convert assignments to colors
    face_colors
        .into_iter()
        .map(|i| PALETTE[i.unwrap() as usize])
        .collect::<Vec<_>>()
        .into_boxed_slice()
}

/// The color palette used for coloring faces.
const PALETTE: &[Color] = &[
    [1.0, 0.2, 0.2], // Red
    [0.2, 0.8, 0.2], // Green
    [0.1, 0.4, 1.0], // Blue
    [1.0, 1.0, 0.2], // Yellow
    [1.0, 0.2, 1.0], // Magenta
    [0.2, 1.0, 1.0], // Cyan
    [0.7, 0.7, 0.7], // Gray
];

/// Converts a [`Color`] to an [`image::Rgb<u8>`](image::Rgb).
fn to_srgb(linear: Color) -> image::Rgb<u8> {
    image::Rgb(palette::Srgb::<u8>::from(palette::LinSrgb::from(linear)).into())
}

/// Computes the intersection of a ray with the unit sphere.
///
/// Specifically, this returns the smallest non-negative `t`, such that the point `start + dir * t`
/// is on the sphere, or `None` if the ray does not intersect the sphere.
fn trace_sphere(start: [f64; 3], dir: [f64; 3]) -> Option<f64> {
    let a = 1.0;
    let b = 2.0 * vec::dot(start, dir);
    let c = vec::dot(start, start) - 1.0;
    let disc = b * b - 4.0 * a * c;
    if disc >= 0.0 {
        let u = -b - b.signum() * disc.sqrt();
        let t_0 = u / (2.0 * a);
        let t_1 = 2.0 * c / u;
        [t_0, t_1]
            .into_iter()
            .filter(|&t| t >= 0.0)
            .min_by(|a, b| a.partial_cmp(b).unwrap())
    } else {
        None
    }
}
