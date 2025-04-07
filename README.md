# `subsphere` 
[![](https://img.shields.io/crates/v/subsphere.svg)](https://crates.io/crates/subsphere)
[![](https://docs.rs/subsphere/badge.svg)](https://docs.rs/subsphere/)

**Sphere tessellation toolkit**

This crate provides a general, easy-to-use API for working with tessellated spheres, i.e. spheres
whose surface is partitioned into polygonal cells. It includes implementations for a variety of
such tessellations.

## Features

* **Implicit Representation:** Instead of storing geometry data directly, the
tessellations, and the elements within them, are represented implicitly. They are compact
zero-allocation `Copy` types that can be used to generate geometry on the fly. This lets you work
with massive tessellations using very little memory.

* **Spherical Geometry:** We aren't just working with flat meshes and projecting the results
onto a sphere. This crate uses sphere-aware math throughout, giving us less distortion and
more consistent cell shapes than naive implementations.

* **Versatility:** This crate provides a variety of tessellations, including triangular
tessellations based on any [Geodesic Polyhedron](https://en.wikipedia.org/wiki/Geodesic_polyhedron)
and hexagonal tessellations based on any
[Goldberg Polyhedron](https://en.wikipedia.org/wiki/Goldberg_polyhedron), all sharing a common
API. This lets you quickly experiment with them and find the one that works best for your
application.