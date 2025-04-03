use std::num::NonZero;

use subsphere::{Face, Sphere, Vertex};

#[test]
fn test_icosphere_base() {
    let sphere = subsphere::IcoSphere::base();
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
}

#[test]
fn test_icosphere_4_0() {
    let sphere = subsphere::IcoSphere::new(NonZero::new(4).unwrap(), 0);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
}

#[test]
fn test_icosphere_3_1() {
    let sphere = subsphere::IcoSphere::new(NonZero::new(3).unwrap(), 1);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
}

#[test]
fn test_icosphere_2_2() {
    let sphere = subsphere::IcoSphere::new(NonZero::new(2).unwrap(), 2);
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
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