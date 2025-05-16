use std::num::NonZero;
use subsphere::util::validate;
use subsphere::prelude::*;
use subsphere::{proj, TriSphere, BaseTriSphere, HexSphere};

#[test]
fn test_octasphere() {
    let sphere = subsphere::octasphere();
    validate_tri(sphere);
    validate(sphere);
}

#[test]
fn test_trisphere_octa_gnomonic_4_0() {
    let sphere = TriSphere::new(
        BaseTriSphere::Octa,
        proj::Gnomonic,
        NonZero::new(4).unwrap(),
        0,
    );
    assert_eq!(
        sphere,
        subsphere::octasphere()
            .subdivide_edge(NonZero::new(4).unwrap())
            .with_projector(proj::Gnomonic)
    );
    validate_tri(sphere);
    validate(sphere);
}

#[test]
fn test_trisphere_octa_gnomonic_3_1() {
    let sphere = TriSphere::new(
        BaseTriSphere::Octa,
        proj::Gnomonic,
        NonZero::new(3).unwrap(),
        1,
    );
    validate_tri(sphere);
    validate(sphere);
}

#[test]
fn test_trisphere_octa_gnomonic_2_2() {
    let sphere = TriSphere::new(
        BaseTriSphere::Octa,
        proj::Gnomonic,
        NonZero::new(2).unwrap(),
        2,
    );
    validate_tri(sphere);
    validate(sphere);
}

#[test]
fn test_icosphere() {
    let sphere = TriSphere::from(BaseTriSphere::Icosa);
    assert_eq!(sphere, subsphere::icosphere());
    validate_tri(sphere);
    validate(sphere);
}

#[test]
fn test_trisphere_icosa_gnomonic_4_0() {
    let sphere = TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Gnomonic,
        NonZero::new(4).unwrap(),
        0,
    );
    validate_tri(sphere);
    validate(sphere);
}

#[test]
fn test_trisphere_icosa_gnomonic_3_1() {
    let sphere = TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Gnomonic,
        NonZero::new(3).unwrap(),
        1,
    );
    validate_tri(sphere);
    validate(sphere);
}

#[test]
fn test_trisphere_icosa_gnomonic_2_2() {
    let sphere = TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Gnomonic,
        NonZero::new(2).unwrap(),
        2,
    );
    validate_tri(sphere);
    validate(sphere);
}

#[test]
fn test_hexsphere_icosa_gnomonic_6_0() {
    let sphere = HexSphere::from_kis(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Gnomonic,
        NonZero::new(6).unwrap(),
        0,
    ))
    .unwrap();
    assert_eq!(
        sphere,
        subsphere::icosphere()
            .subdivide_edge(NonZero::new(2).unwrap())
            .with_projector(proj::Gnomonic)
            .truncate()
    );
    validate(sphere);
}

#[test]
fn test_hexsphere_icosa_gnomonic_4_1() {
    let sphere = HexSphere::from_kis(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Gnomonic,
        NonZero::new(4).unwrap(),
        1,
    ))
    .unwrap();
    validate(sphere);
}

#[test]
fn test_hexsphere_icosa_gnomonic_7_1() {
    let sphere = HexSphere::from_kis(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Gnomonic,
        NonZero::new(7).unwrap(),
        1,
    ))
    .unwrap();
    validate(sphere);
}

#[test]
fn test_hexsphere_icosa_gnomonic_2_2() {
    let sphere = HexSphere::from_kis(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Gnomonic,
        NonZero::new(2).unwrap(),
        2,
    ))
    .unwrap();
    validate(sphere);
}

#[test]
fn test_hexsphere_icosa_gnomonic_8_2() {
    let sphere = HexSphere::from_kis(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Gnomonic,
        NonZero::new(8).unwrap(),
        2,
    ))
    .unwrap();
    validate(sphere);
}

#[test]
fn test_hexsphere_icosa_gnomonic_3_3() {
    let sphere = HexSphere::from_kis(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Gnomonic,
        NonZero::new(3).unwrap(),
        3,
    ))
    .unwrap();
    validate(sphere);
}

/// Performs [`TriSphere`]-specific validation.
fn validate_tri<Proj>(sphere: TriSphere<Proj>)
where
    Proj: Eq + Clone + proj::BaseTriProjector,
    Proj: std::fmt::Debug,
{
    // Verify that faces contain their centers.
    for f in sphere.faces() {
        let center = f.center();
        debug_assert_eq!(sphere.face_at(center), f);
    }
}

/// Converts the given [`subsphere::Sphere`] to an OBJ file, returning the UTF-8 bytes of the file.
fn to_obj(sphere: impl subsphere::Sphere) -> Vec<u8> {
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
