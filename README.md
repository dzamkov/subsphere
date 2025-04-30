# `subsphere` 
[![](https://img.shields.io/crates/v/subsphere.svg)](https://crates.io/crates/subsphere)
[![](https://docs.rs/subsphere/badge.svg)](https://docs.rs/subsphere/)

**Sphere tessellation toolkit**

This crate provides a general, easy-to-use API for working with tessellated spheres, i.e. spheres
whose surface is partitioned into polygonal cells. It includes implementations for a variety of
such tessellations.

![Example hexagonal tessellation](https://github.com/dzamkov/subsphere/blob/master/render/out/hexsphere_icosa_8_2_fuller.png?raw=true)
![Example triangular tessellation](https://github.com/dzamkov/subsphere/blob/master/render/out/trisphere_icosa_3_1_fuller.png?raw=true)

## Features

* **Implicit Representation:** Instead of storing geometry data directly, the
tessellations, and the elements within them, are represented implicitly. They are compact
zero-allocation `Copy` types that can be used to generate geometry on the fly. This lets you work
with massive tessellations using very little memory.

* **Versatility:** This crate allows you to explore a
[huge variety](https://github.com/dzamkov/subsphere/blob/master/catalog.md) of tessellations, all sharing a
common API. There's a bunch of adjustable parameters and you're free to mix-and-match topologies
and projections to tune the tessellation to your needs.

* **Spherical Geometry:** In the world of `subsphere`, there is only one geometric space: the
sphere. All objects follow the contours of the sphere, and all calculations correctly account for
this.

## Examples

#### Construct a tessellation by refining a base tessellation

```rust
let sphere = subsphere::icosphere()
    .subdivide_edge(NonZero::new(3).unwrap())
    .with_projector(subsphere::proj::Fuller)
    .truncate();
```

#### Construct a tessellation explicitly

```rust
let sphere = subsphere::HexSphere::from_kis(
    subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Fuller,
        NonZero::new(9).unwrap(),
        0,
    )
).unwrap();
```

#### Write a tessellation to an [OBJ file](https://en.wikipedia.org/wiki/Wavefront_.obj_file)

```rust
let mut obj = String::new();
for vert in sphere.vertices() {
    let pos = vert.pos();
    obj.push_str(&format!("v {} {} {}\n", pos[0], pos[1], pos[2]));
}
for face in sphere.faces() {
    let indices = face
        .vertices()
        .map(|vert| format!("{}", vert.index() + 1)) // OBJ indices are 1-based
        .collect::<Vec<_>>();
    obj.push_str(&format!("f {}\n", indices.join(" ")));
}
std::fs::write("sphere.obj", &obj).expect("failed to write to file");
```

#### Toggle neighboring faces when a face is clicked

```rust
// Setup
let mut is_active = vec![false; sphere.num_faces()];

// For each click
let click_face = sphere.face_at(click_point);
for side in click_face.sides() {
    let neighbor_index = side.twin().inside().index();
    is_active[neighbor_index] = !is_active[neighbor_index];
}
```