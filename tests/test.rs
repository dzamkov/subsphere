use std::num::NonZero;
use subsphere::util::validate;
use subsphere::{BaseTriSphere, Face, HalfEdge, HexSphere, Sphere, TriSphere, Vertex, proj};

#[test]
fn test_octosphere_base() {
    let sphere = TriSphere::from(BaseTriSphere::Octo);
    assert_eq!(sphere, subsphere::octosphere());
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    let metrics = Metrics::new(sphere.faces());
    insta::assert_snapshot!(metrics.area_discrepancy, @"1");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0");
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
    let metrics = Metrics::new(sphere.faces());
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.5184953070560394");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.3955428804408463");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.7458526607730736");
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
    let metrics = Metrics::new(sphere.faces());
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.8924780020022631");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.4647312208585381");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.7426858548903519");
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
    let metrics = Metrics::new(sphere.faces());
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.6794205913771127");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.4938677503333329");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.741758917736975");
}

#[test]
fn test_icosphere_base() {
    let sphere = TriSphere::from(BaseTriSphere::Icosa);
    assert_eq!(sphere, subsphere::icosphere());
    validate(sphere);
    insta::assert_binary_snapshot!(".obj", to_obj(sphere));
    let metrics = Metrics::new(sphere.faces());
    insta::assert_snapshot!(metrics.area_discrepancy, @"1");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.0000000000000002");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.0000000000000006661338147750939");
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
    let metrics = Metrics::new(sphere.faces());
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.149592490499728");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.1656669893755778");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.3077306793682705");
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
    let metrics = Metrics::new(sphere.faces());
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.2142304244815916");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.1668416265723973");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.2924575851995208");
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
    let metrics = Metrics::new(sphere.faces());
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.17721464511508");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.1725986697633959");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.29094312908498365");
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
    let metrics = Metrics::new(sphere.faces().filter(|f| f.num_sides() == 6));
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.0509668034000925");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.171188457818225");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.28834568593235255");
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
    let metrics = Metrics::new(sphere.faces().filter(|f| f.num_sides() == 6));
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.000000000000003");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.1262732071592574");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.16400132885900653");
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
    let metrics = Metrics::new(sphere.faces().filter(|f| f.num_sides() == 6));
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.052006826012178");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.167013900454752");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.32400169018257974");
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
    let metrics = Metrics::new(sphere.faces().filter(|f| f.num_sides() == 6));
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.0000000000000009");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.0606276105914663");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.08150361884566548");
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
    let metrics = Metrics::new(sphere.faces().filter(|f| f.num_sides() == 6));
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.0661630517609373");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.1805389357043357");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.3322464937015097");
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
    let metrics = Metrics::new(sphere.faces().filter(|f| f.num_sides() == 6));
    insta::assert_snapshot!(metrics.area_discrepancy, @"1.0350802014863203");
    insta::assert_snapshot!(metrics.length_discrepancy, @"1.1293976478623056");
    insta::assert_snapshot!(metrics.angle_discrepancy, @"0.16753639503241402");
}

/// Encapsulates measurements for the "quality" of a sphere tessellation.
struct Metrics {
    /// The area of the largest face divided by the area of the smallest face.
    area_discrepancy: f64,

    /// The length of the longest edge divided by the length of the shortest edge.
    length_discrepancy: f64,

    /// The positive difference between the largest interior angle and the smallest interior angle,
    /// in radians.
    angle_discrepancy: f64,
}

impl Metrics {
    /// Computes [`Metrics`] for the given set of faces.
    pub fn new<Face: subsphere::Face>(faces: impl Iterator<Item = Face>) -> Self {
        let mut min_area = f64::MAX;
        let mut max_area = f64::MIN;
        let mut min_length = f64::MAX;
        let mut max_length = f64::MIN;
        let mut min_angle = f64::MAX;
        let mut max_angle = f64::MIN;
        for f in faces {
            let area = f.area();
            min_area = min_area.min(area);
            max_area = max_area.max(area);
            for s in f.sides() {
                let length = s.length();
                min_length = min_length.min(length);
                max_length = max_length.max(length);
                let angle = s.angle();
                min_angle = min_angle.min(angle);
                max_angle = max_angle.max(angle);
            }
        }
        Self {
            area_discrepancy: max_area / min_area,
            length_discrepancy: max_length / min_length,
            angle_discrepancy: max_angle - min_angle,
        }
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
