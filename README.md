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

* **Versatility:** This crate allows you to explore a huge variety of tessellations, all sharing a
common API. There's a bunch of adjustable parameters and you're free to mix-and-match topologies
and projections to tune the tessellation to your needs.

* **Spherical Geometry:** In the world of `subsphere`, there is only one geometric space: the
sphere. All objects follow the contours of the sphere, and all calculations correctly account for
this.

## Usage

### Constructing a tessellation

The first step to using `subsphere` is to specify which tessellation you are working with.
Tessellations implement the [`Sphere`](https://docs.rs/subsphere/latest/subsphere/trait.Sphere.html)
trait. There's two main ways to construct a tessellation:

**Refining a base tessellation**

```rust
let sphere = subsphere::icosphere()
    .subdivide_edge(NonZero::new(3).unwrap())
    .with_projector(subsphere::proj::Fuller)
    .truncate();
```

**Explicitly**

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

These examples both yield the same hexagonal tessellation. Note that there are currently some
tessellations which can only be constructed explicitly.