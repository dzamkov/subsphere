use std::num::NonZero;

use subsphere::{Face, Sphere, Vertex};

#[test]
fn test_icosphere_1() {
    let sphere = subsphere::IcoSphere::new(NonZero::new(1).unwrap());
    validate_vertex_indices(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
}

#[test]
fn test_icosphere_4() {
    let sphere = subsphere::IcoSphere::new(NonZero::new(4).unwrap());
    validate_vertex_indices(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
}

/// Validates that [`Vertex::index`] is consistent with the [`Sphere::vertices`] iterator.
fn validate_vertex_indices(sphere: impl Sphere) {
    for (i, v) in sphere.vertices().enumerate() {
        assert_eq!(v.index(), i, "vertex index mismatch");
    }
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