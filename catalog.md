# Tessellations Catalog

This document provides a listing of all tessellations implemented in this crate, along with
measurements of their "quality" in terms of certain metrics:

* **Area:** The relative difference between the area of the largest face and the smallest face.
This excludes non-hexagonal faces in hexagonal tessellations.

* **Length:** Considering faces whose vertices have the same degree, the maximum relative
difference between the length of the longest edge and the shortest edge of that face.

* **Angle:** Considering faces whose vertices have the same degree, the maximum absolute difference
between the interior angles of the face.

Since most types of tessellations have configurable integer parameters, representative instances
are chosen such that the number of faces is approximately 2000. This provides a "fair" comparison
between tessellations of different types.

## [`TriSphere`](https://docs.rs/subsphere/latest/subsphere/tri/struct.TriSphere.html)

![Example triangular tessellation](https://github.com/dzamkov/subsphere/blob/master/render/out/trisphere_icosa_3_1_fuller.png?raw=true)

| Base    | Projector  | B  | C | # Faces | # Verts | Area    | Length | Angle   |
| ------- | ---------- | -- | - | ------- | ------- | ------- | ------ | ------- |
| `Icosa` | `Gnomonic` | 10 | 0 | 2000    | 1002    | 78.39%  | 14.65% | 14.42°  |
| `Icosa` | `Gnomonic` | 8  | 3 | 1940    | 972     | 93.86%  | 14.51% | 13.83°  |
| `Icosa` | `Gnomonic` | 6  | 6 | 2160    | 1082    | 99.99%  | 12.13% | 11.31°  |
| `Icosa` | `Fuller`   | 10 | 0 | 2000    | 1002    | 14.73%  | 16.96% | 16.72°  |
| `Icosa` | `Fuller`   | 8  | 3 | 1940    | 972     | 24.12%  | 18.30% | 16.67°  |
| `Icosa` | `Fuller`   | 6  | 6 | 2160    | 1082    | 26.27%  | 16.79% | 15.89°  |
| `Octo`  | `Gnomonic` | 16 | 0 | 2048    | 1026    | 356.37% | 40.52% | 40.90°  |
| `Octo`  | `Gnomonic` | 13 | 5 | 2072    | 1038    | 499.91% | 47.67% | 40.03°  |
| `Octo`  | `Gnomonic` | 9  | 9 | 1944    | 974     | 559.14% | 43.09% | 32.14°  |
| `Octo`  | `Fuller`   | 16 | 0 | 2048    | 1026    | 62.60%  | 41.02% | 42.63°  |
| `Octo`  | `Fuller`   | 13 | 5 | 2072    | 1038    | 111.36% | 53.81% | 42.56°  |
| `Octo`  | `Fuller`   | 9  | 9 | 1944    | 974     | 129.43% | 50.51% | 40.73°  |
 |
## [`HexSphere`](https://docs.rs/subsphere/latest/subsphere/hex/struct.HexSphere.html)

![Example hexagonal tessellation](https://github.com/dzamkov/subsphere/blob/master/render/out/hexsphere_icosa_8_2_fuller.png?raw=true)

| Base    | Projector  | B  | C  | # Faces | # Verts | Area    | Length | Angle   |
| ------- | ---------- | -- | -- | ------- | ------- | ------- | ------ | ------- |
| `Icosa` | `Gnomonic` | 24 | 0  | 1922    | 3840    | 73.75%  | 17.58% | 30.21°  |
| `Icosa` | `Gnomonic` | 20 | 8  | 2082    | 4160    | 72.76%  | 20.35% | 18.53°  |
| `Icosa` | `Gnomonic` | 14 | 14 | 1962    | 3920    | 73.64%  | 19.30% | 12.88°  |
| `Icosa` | `Fuller`   | 24 | 0  | 1922    | 3840    | 11.23%  | 17.53% | 31.99°  |
| `Icosa` | `Fuller`   | 20 | 8  | 2082    | 4160    | 10.27%  | 21.02% | 23.10°  |
| `Icosa` | `Fuller`   | 14 | 14 | 1962    | 3920    | 11.07%  | 17.92% | 15.75°  |
| `Octo`  | `Gnomonic` | 39 | 0  | 2030    | 4056    | 343.35% | 44.49% | 83.45°  |
| `Octo`  | `Gnomonic` | 32 | 11 | 1998    | 3992    | 335.61% | 56.14% | 64.48°  |
| `Octo`  | `Gnomonic` | 22 | 22 | 1938    | 3872    | 341.57% | 56.72% | 32.71°  |
| `Octo`  | `Fuller`   | 39 | 0  | 2030    | 4056    | 55.81%  | 41.78% | 84.18°  |
| `Octo`  | `Fuller`   | 32 | 11 | 1998    | 3992    | 57.37%  | 56.55% | 69.86°  |
| `Octo`  | `Fuller`   | 22 | 22 | 1938    | 3872    | 55.44%  | 54.55% | 34.29°  |

