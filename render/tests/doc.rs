//! This test module produces renders for the documentation of the `subsphere` crate.

#[test]
fn test_trisphere_icosa_3_1_gnomonic() {
    let sphere = subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Gnomonic,
        std::num::NonZero::new(3).unwrap(),
        1,
    );
    render(&sphere, "trisphere_icosa_3_1_gnomonic");
}

#[test]
fn test_hexsphere_icosa_8_2_gnomonic() {
    let sphere = subsphere::HexSphere::new(subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Gnomonic,
        std::num::NonZero::new(8).unwrap(),
        2,
    ))
    .unwrap();
    render(&sphere, "hexsphere_icosa_8_2_gnomonic");
}

/// Renders the given sphere with the given name.
fn render(sphere: &impl subsphere::Sphere, name: &str) {
    // Render image
    let mut im = image::RgbImage::new(256, 256);
    subsphere_render::render(
        &subsphere_render::Scene {
            background: [1.0, 1.0, 1.0],
            light_dir: [1.0, 0.0, 1.0],
            ambient: 0.4,
            sphere,
            faces: subsphere_render::colorize(sphere),
        },
        &subsphere_render::Camera::perspective([0.0, -1.0, 1.0], [0.5, 0.5], 1.0),
        &mut im,
    );

    // Compare to reference
    let path = format!("out/{}.png", name);
    if std::env::var("OVERWRITE").is_ok() {
        im.save(path).unwrap();
    } else {
        let ref_im = image::open(&path)
            .expect("no reference image found")
            .to_rgb8();
        assert!(
            im == ref_im,
            "rendered image does not match reference: {}",
            name
        );
    }
}
