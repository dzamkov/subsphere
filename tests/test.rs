use std::num::NonZero;
use subsphere::{BaseTriSphere, Face, HalfEdge, HexSphere, Sphere, TriSphere, Vertex, proj};

#[test]
fn test_octosphere_base() {
    let sphere = TriSphere::from(BaseTriSphere::Octo);
    assert_eq!(sphere, subsphere::octosphere());
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(3, sphere), @"1");
}

#[test]
fn test_octosphere_4_0() {
    let sphere = TriSphere::new(
        BaseTriSphere::Octo,
        proj::Fuller,
        NonZero::new(4).unwrap(),
        0,
    );
    assert_eq!(
        sphere,
        subsphere::octosphere().subdivide_edge(NonZero::new(4).unwrap())
    );
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(3, sphere), @"1.5184953070560394");
}

#[test]
fn test_octosphere_3_1() {
    let sphere = TriSphere::new(
        BaseTriSphere::Octo,
        proj::Fuller,
        NonZero::new(3).unwrap(),
        1,
    );
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(3, sphere), @"1.8924780020022631");
}

#[test]
fn test_octosphere_2_2() {
    let sphere = TriSphere::new(
        BaseTriSphere::Octo,
        proj::Fuller,
        NonZero::new(2).unwrap(),
        2,
    );
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(3, sphere), @"1.6794205913771127");
}

#[test]
fn test_icosphere_base() {
    let sphere = TriSphere::from(BaseTriSphere::Icosa);
    assert_eq!(sphere, subsphere::icosphere());
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
}

#[test]
fn test_icosphere_4_0() {
    let sphere = TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(4).unwrap(),
        0,
    );
    assert_eq!(
        sphere,
        subsphere::icosphere().subdivide_edge(NonZero::new(4).unwrap())
    );
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(3, sphere), @"1.149592490499728");
}

#[test]
fn test_icosphere_3_1() {
    let sphere = TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(3).unwrap(),
        1,
    );
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(3, sphere), @"1.2142304244815916");
}

#[test]
fn test_icosphere_2_2() {
    let sphere = TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(2).unwrap(),
        2,
    );
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(3, sphere), @"1.17721464511508");
}

#[test]
fn test_hexsphere_6_0() {
    let sphere = HexSphere::new(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(6).unwrap(),
        0,
    ))
    .unwrap();
    assert_eq!(
        sphere,
        subsphere::icosphere()
            .subdivide_edge(NonZero::new(2).unwrap())
            .truncate()
    );
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(6, sphere), @"1.0509668034000925");
}

#[test]
fn test_hexsphere_4_1() {
    let sphere = HexSphere::new(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(4).unwrap(),
        1,
    ))
    .unwrap();
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(6, sphere), @"1.000000000000003");
}

#[test]
fn test_hexsphere_7_1() {
    let sphere = HexSphere::new(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(7).unwrap(),
        1,
    ))
    .unwrap();
    validate(sphere);
}

#[test]
fn test_hexsphere_2_2() {
    let sphere = HexSphere::new(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(2).unwrap(),
        2,
    ))
    .unwrap();
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    insta::assert_snapshot!(area_discrepancy(6, sphere), @"1.0000000000000009");
}

#[test]
fn test_hexsphere_8_2() {
    let sphere = HexSphere::new(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(8).unwrap(),
        2,
    ))
    .unwrap();
    validate(sphere);
}

#[test]
fn test_hexsphere_3_3() {
    let sphere = HexSphere::new(TriSphere::new(
        BaseTriSphere::Icosa,
        proj::Fuller,
        NonZero::new(3).unwrap(),
        3,
    ))
    .unwrap();
    validate(sphere);
}

/// Validates the internal consistency of the given [`Sphere`].
fn validate<S: Sphere>(sphere: S)
where
    S::Face: std::fmt::Debug,
    S::Vertex: std::fmt::Debug,
    S::HalfEdge: std::fmt::Debug,
{
    // Validate faces
    let mut index = 0;
    for f in sphere.faces() {
        assert_eq!(f.index(), index, "face index mismatch");
        index += 1;

        // Validate sides
        let mut side_index = 0;
        let mut prev = None;
        for s in f.sides() {
            assert_eq!(s.inside(), f, "side inside mismatch");
            assert_eq!(s.start(), f.vertex(side_index), "side start mismatch");
            assert_eq!(s.side_index(), side_index, "side index mismatch");
            side_index += 1;
            if let Some(prev) = prev {
                assert_eq!(s.prev(), prev, "side prev mismatch");
                assert_eq!(prev.next(), s, "side next mismatch");
            }
            prev = Some(s);
        }
        if let Some(prev) = prev {
            let s = f.side(0);
            assert_eq!(s.prev(), prev, "side prev mismatch");
            assert_eq!(prev.next(), s, "side next mismatch");
        }
        assert_eq!(side_index, f.num_sides(), "side count mismatch");
    }
    assert_eq!(index, sphere.num_faces(), "face count mismatch");

    // Validate vertices
    let mut index = 0;
    for v in sphere.vertices() {
        assert_eq!(v.index(), index, "vertex index mismatch");
        index += 1;

        // Validate outgoing edges
        let mut edge_index = 0;
        let mut prev = None;
        for e in v.outgoings() {
            assert_eq!(
                e.start(),
                v,
                "outgoing edge with index {} does not start at vertex",
                edge_index
            );
            assert_eq!(e.outgoing_index(), edge_index, "edge index mismatch");
            edge_index += 1;
            if let Some(prev) = prev {
                assert_eq!(e.complement().next(), prev, "outgoing prev mismatch");
                assert_eq!(prev.prev().complement(), e, "outgoing next mismatch");
            }
            prev = Some(e);
        }
        if let Some(prev) = prev {
            let e = v.outgoing(0);
            assert_eq!(e.complement().next(), prev, "outgoing prev mismatch");
            assert_eq!(prev.prev().complement(), e, "outgoing next mismatch");
        }
        assert_eq!(edge_index, v.degree(), "vertex degree mismatch");
    }
    assert_eq!(index, sphere.num_vertices(), "vertex count mismatch");

    // Validate total face area
    let total_area: f64 = sphere.faces().map(|f| f.area()).sum();
    assert!(
        // TODO: Find source of inaccuracy and tighten this tolerance
        (total_area - 4.0 * std::f64::consts::PI).abs() < 1e-3,
        "total area mismatch: expected {:?}, got {:?}",
        4.0 * std::f64::consts::PI,
        total_area
    );
}

/// The area of the largest face on the given [`Sphere`] divided by the area of the smallest face.
///
/// Only counts faces with the given number of sides. Smaller values indicate a more uniform
/// distribution of face areas, which is desirable.
fn area_discrepancy(num_sides: usize, sphere: impl Sphere) -> f64 {
    let mut min_area = f64::MAX;
    let mut max_area = f64::MIN;
    for f in sphere.faces() {
        if f.num_sides() == num_sides {
            let area = f.area();
            min_area = min_area.min(area);
            max_area = max_area.max(area);
        }
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
