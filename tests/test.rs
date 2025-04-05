use std::num::NonZero;

use subsphere::{BaseTriSphere, Face, Sphere, SubTriSphere, Vertex};

#[test]
fn test_octosphere_base() {
    let sphere = SubTriSphere::from(BaseTriSphere::Octo);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(sphere), @"1");
}

#[test]
fn test_octosphere_4_0() {
    let sphere = SubTriSphere::new(BaseTriSphere::Octo, NonZero::new(4).unwrap(), 0);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(sphere), @"1.5543100839896127");
}

#[test]
fn test_octosphere_3_1() {
    let sphere = SubTriSphere::new(BaseTriSphere::Octo, NonZero::new(3).unwrap(), 1);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(sphere), @"1.392702192665517");
}

#[test]
fn test_octosphere_2_2() {
    let sphere = SubTriSphere::new(BaseTriSphere::Octo, NonZero::new(2).unwrap(), 2);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(sphere), @"1.2450941184111652");
}

#[test]
fn test_icosphere_base() {
    let sphere = SubTriSphere::from(BaseTriSphere::Icosa);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
}

#[test]
fn test_icosphere_4_0() {
    let sphere = SubTriSphere::new(BaseTriSphere::Icosa, NonZero::new(4).unwrap(), 0);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(sphere), @"1.1643667123353725");
}

#[test]
fn test_icosphere_3_1() {
    let sphere = SubTriSphere::new(BaseTriSphere::Icosa, NonZero::new(3).unwrap(), 1);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(sphere), @"1.16284750155026");
}

#[test]
fn test_icosphere_2_2() {
    let sphere = SubTriSphere::new(BaseTriSphere::Icosa, NonZero::new(2).unwrap(), 2);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(sphere), @"1.073344298185472");
}

/// Validates the internal consistency of the given [`Sphere`].
fn validate(sphere: impl Sphere) {
    // Validate faces
    let mut index = 0;
    for f in sphere.faces() {
        assert_eq!(f.index(), index, "face index mismatch");
        index += 1;
    }
    assert_eq!(index, sphere.num_faces(), "face count mismatch");

    // Validate vertices
    let mut index = 0;
    for v in sphere.vertices() {
        assert_eq!(v.index(), index, "vertex index mismatch");
        index += 1;
    }
    assert_eq!(index, sphere.num_vertices(), "vertex count mismatch");

    // Validate total face area
    let total_area: f64 = sphere.faces().map(|f| f.area()).sum();
    assert!(
        (total_area - 4.0 * std::f64::consts::PI).abs() < 1e-12,
        "total area mismatch: expected {:?}, got {:?}",
        4.0 * std::f64::consts::PI,
        total_area
    );
}

/// The area of the largest face on the given [`Sphere`] divided by the area of the smallest face.
///
/// Smaller values indicate a more uniform distribution of face areas, which is desirable.
fn area_discrepancy(sphere: impl Sphere) -> f64 {
    let mut min_area = f64::MAX;
    let mut max_area = f64::MIN;
    for f in sphere.faces() {
        let area = f.area();
        min_area = min_area.min(area);
        max_area = max_area.max(area);
    }
    max_area / min_area
}

/// Converts the given [`Sphere`] to an OBJ file, returning the UTF-8 bytes of the file.
fn to_obj(sphere: impl Sphere) -> Vec<u8> {
    let mut obj = String::new();
    for v in sphere.vertices() {
        let pos = v.pos();
        obj.push_str(&format!("v {} {} {}\n", pos[0], pos[1], pos[2]));
    }
    for f in sphere.faces() {
        let indices = f
            .vertices()
            .map(|v| format!("{}", v.index() + 1)) // OBJ indices are 1-based
            .collect::<Vec<_>>();
        obj.push_str(&format!("f {}\n", indices.join(" ")));
    }
    obj.into_bytes()
}
