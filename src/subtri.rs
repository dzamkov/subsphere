//! Contains types related to [`SubTriSphere`].
use crate::Face as FaceExt;
use crate::HalfEdge as HalfEdgeExt;
use crate::Sphere as SphereExt;
use crate::Vertex as VertexExt;
use crate::basetri::BaseTriSphere;
use crate::basetri::Face as BaseFace;
use crate::basetri::Vertex as BaseVertex;
use crate::math::{Vector3, mat3, vec3};
use std::num::NonZeroU32;

/// A tessellation of the unit sphere into triangular [`Face`]s constructed by subdividing a
/// triangular platonic solid and projecting it onto the sphere.
///
/// This subdivision of a platonic solid is also known as a
/// [geodesic polyhedron](https://en.wikipedia.org/wiki/Geodesic_polyhedron).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SubTriSphere {
    b: NonZeroU32,
    c_base: u32,
}

impl SubTriSphere {
    /// Constructs a [`SubTriSphere`] from the given base shape and subdivision parameters.
    ///
    /// The `b` and `c` parameters are as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation). For implementation
    /// simplicity and performance, only `b ≥ c` polyhedra are directly supported. `b < c`
    /// polyhedra can be emulated by swapping the `b` and `c` parameters and mirroring the
    /// resulting sphere.
    pub fn new(base: BaseTriSphere, b: NonZeroU32, c: u32) -> Self {
        assert!(
            b.get() < u32::MAX >> 2,
            "b exceeds maxmimum subdivision size. b = {}",
            b
        );
        assert!(
            b.get() >= c,
            "b must be greater than or equal to c. b = {}, c = {}",
            b,
            c
        );
        Self {
            b,
            c_base: (c << 2) | base as u32,
        }
    }

    /// The base shape of this subdivided sphere.
    pub fn base(&self) -> BaseTriSphere {
        unsafe { std::mem::transmute((self.c_base & 0b11) as u8) }
    }

    /// The `b` parameter of the subdivision, as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation).
    pub fn b(&self) -> u32 {
        self.b.get()
    }

    /// The `c` parameter of the subdivision, as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation).
    pub fn c(&self) -> u32 {
        self.c_base >> 2
    }

    /// The number of [`Face`]s each triangle on the base shape is subdivided into.
    ///
    /// This is also known as the "triangulation number" of the shape.
    pub fn num_divisions(&self) -> usize {
        let b = self.b() as usize;
        let c = self.c() as usize;
        b * (b + c) + c * c
    }

    /// The number of faces per edge [`BaseRegion`].
    fn num_faces_per_edge_region(&self) -> usize {
        let b = self.b() as usize;
        let c = self.c() as usize;
        b * c * 2
    }

    /// The number of faces per interior [`BaseRegion`].
    fn num_faces_per_interior_region(&self) -> usize {
        let n = (self.b() - self.c()) as usize;
        n * n
    }

    /// The number of vertices per edge [`BaseRegion`], assuming that the origin vertex is not
    /// counted.
    fn num_vertices_per_edge_region(&self) -> usize {
        let b = self.b() as usize;
        let c = self.c() as usize;
        (b - 1) * (c + 1)
    }

    /// The number of vertices per interior [`BaseRegion`].
    fn num_vertices_per_interior_region(&self) -> usize {
        let n = (self.b() - self.c()) as isize;
        // Yes, this really should return `1` when `n = 0`.
        ((n - 1) * (n - 2) / 2) as usize
    }
}

impl From<BaseTriSphere> for SubTriSphere {
    fn from(base: BaseTriSphere) -> Self {
        Self::new(base, NonZeroU32::new(1).unwrap(), 0)
    }
}

impl crate::Sphere for SubTriSphere {
    type Face = Face;
    type Vertex = Vertex;
    type HalfEdge = HalfEdge;

    fn num_faces(&self) -> usize {
        self.num_divisions() * self.base().num_faces()
    }

    fn face(&self, index: usize) -> Face {
        todo!()
    }

    #[allow(refining_impl_trait)]
    fn faces(&self) -> FaceIter {
        self.faces_from(self.base().first_region())
    }

    fn num_vertices(&self) -> usize {
        self.num_faces() / 2 + 2
    }

    fn vertex(&self, index: usize) -> Vertex {
        todo!()
    }

    #[allow(refining_impl_trait)]
    fn vertices(&self) -> VertexIter {
        VertexIter {
            sphere: *self,
            region: self.base().first_region(),
            u: 0,
            u_end: self.b(),
            v: 0,
        }
    }
}

/// Represents a face in an [`SubTriSphere`].
///
/// This is always a geodesic triangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Face {
    sphere: SubTriSphere,

    /// The [`BaseRegion`] which "owns" this face.
    region: BaseRegion,

    /// The U coordinate of the first vertex of this face, specified in the local coordinate space
    /// of `region`.
    u_0: u32,

    /// The U coordinate of the first vertex of this face, specified in the local coordinate space
    /// of `region`.
    v_0: u32,

    /// If `true`, the first [`HalfEdge`] boundary of this face will go in the +V direction.
    /// Otherwise, it will go in the +U direction.
    boundary_along_v: bool,
}

impl crate::Face for Face {
    type Vertex = Vertex;
    type HalfEdge = HalfEdge;

    fn index(&self) -> usize {
        if self.region.ty().is_edge() {
            self.sphere.base_face_index(self.region)
                + (2 * self.v_0 as usize + self.boundary_along_v as usize)
                    * self.sphere.b() as usize
                + self.u_0 as usize
                - self.boundary_along_v as usize
        } else {
            let n = (self.sphere.b() - self.sphere.c()) as usize;
            self.sphere.base_face_index(self.region)
                + self.v_0 as usize * n * 2
                + self.boundary_along_v as usize * (n - self.v_0 as usize - 1)
                - self.v_0 as usize * self.v_0 as usize
                + self.u_0 as usize
        }
    }

    fn area(&self) -> f64 {
        let [[u_0, v_0], [u_1, v_1], [u_2, v_2]] = self.local_coords();
        let p_0 = self.sphere.project(self.region, u_0 as f64, v_0 as f64);
        let p_1 = self.sphere.project(self.region, u_1 as f64, v_1 as f64);
        let p_2 = self.sphere.project(self.region, u_2 as f64, v_2 as f64);
        let d = 1.0 + vec3::dot(p_0, p_1) + vec3::dot(p_1, p_2) + vec3::dot(p_2, p_0);
        (mat3::det([p_0, p_1, p_2]) / d).atan() * 2.0
    }

    fn num_sides(&self) -> usize {
        3
    }

    fn side(&self, index: usize) -> HalfEdge {
        assert!(index < 3, "side index out of bounds: {}", index);
        let (d_u, d_v, dir) = [
            (0, 0, HalfEdgeDir::Up),
            (0, 0, HalfEdgeDir::Vp),
            (1, 0, HalfEdgeDir::UnVp),
            (0, 1, HalfEdgeDir::Un),
            (0, 1, HalfEdgeDir::Vn),
            (-1, 1, HalfEdgeDir::UpVn),
        ][(index << 1) | self.boundary_along_v as usize];
        HalfEdge {
            sphere: self.sphere,
            region: self.region,
            start_u: self.u_0.wrapping_add_signed(d_u),
            start_v: self.v_0.wrapping_add_signed(d_v),
            dir,
        }
    }
}

impl Face {
    /// Gets the UV coordinates of the vertices of this face, specified in the local coordinate
    /// space of `region`.
    fn local_coords(&self) -> [[u32; 2]; 3] {
        if self.boundary_along_v {
            [
                [self.u_0, self.v_0],
                [self.u_0, self.v_0 + 1],
                [self.u_0 - 1, self.v_0 + 1],
            ]
        } else {
            [
                [self.u_0, self.v_0],
                [self.u_0 + 1, self.v_0],
                [self.u_0, self.v_0 + 1],
            ]
        }
    }
}

/// Represents a vertex in an [`SubTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Vertex {
    sphere: SubTriSphere,

    /// The [`BaseRegion`] which "owns" this vertex.
    region: BaseRegion,

    /// The U coordinate of this vertex, specified in the local coordinate space of `region`.
    u: u32,

    /// The V coordinate of this vertex, specified in the local coordinate space of `region`.
    v: u32,
}

impl Vertex {
    /// Constructs a [`Vertex`] with the given properties.
    ///
    /// This will normalize the vertex `region` to be its proper owner.
    fn new(sphere: SubTriSphere, region: BaseRegion, u: u32, v: u32) -> Self {
        if region.ty().is_edge() {
            if u == 0 {
                if v == 0 {
                    Self::base(sphere, region.owner().vertex(0))
                } else if v == sphere.c() && sphere.b() == sphere.c() {
                    if region.ty() == BaseRegionType::Edge0 {
                        Self::center(sphere, region.owner())
                    } else {
                        Self::center(sphere, region.owner().side(2).complement().inside())
                    }
                } else if region.owner().owns_edge_2() {
                    if region.ty() == BaseRegionType::Edge0 {
                        Self {
                            sphere,
                            region: BaseRegion::new(region.owner(), BaseRegionType::Edge2),
                            u: v,
                            v: 0,
                        }
                    } else {
                        let adj = region.owner().side(2).complement();
                        debug_assert_eq!(adj.index(), 1);
                        Self {
                            sphere,
                            region: BaseRegion::new(adj.inside(), BaseRegionType::Edge0),
                            u: sphere.b() - v,
                            v: sphere.c(),
                        }
                    }
                } else {
                    debug_assert_eq!(region.ty(), BaseRegionType::Edge0);
                    let adj = region.owner().side(2).complement();

                    // With the current base layout, if a base face doesn't own its edge 2,
                    // it will be the edge 0 of the adjacent face that owns it.
                    debug_assert_eq!(adj.index(), 0);
                    Self {
                        sphere,
                        region: BaseRegion::new(adj.inside(), BaseRegionType::Edge0),
                        u: v,
                        v: 0,
                    }
                }
            } else if u == sphere.b() {
                if v == 0 && sphere.b() == sphere.c() {
                    if region.ty() == BaseRegionType::Edge0 {
                        Self::center(sphere, region.owner().side(0).complement().inside())
                    } else {
                        Self::center(sphere, region.owner())
                    }
                } else if region.ty() == BaseRegionType::Edge0 {
                    if v == sphere.c() {
                        Self::base(sphere, region.owner().vertex(1))
                    } else {
                        let adj = region.owner().side(0).complement();
                        if adj.index() == 1 {
                            Self {
                                sphere,
                                region: BaseRegion::new(adj.inside(), BaseRegionType::Edge0),
                                u: sphere.b() - sphere.c() + v,
                                v: sphere.c(),
                            }
                        } else {
                            debug_assert_eq!(adj.index(), 2);
                            Self::adjacent_1(sphere, adj.inside(), (sphere.b() - sphere.c()) + v)
                        }
                    }
                } else if v == sphere.c() {
                    Self::base(sphere, region.owner().vertex(2))
                } else {
                    Self::adjacent_1(sphere, region.owner(), (sphere.b() - sphere.c()) + v)
                }
            } else {
                Self {
                    sphere,
                    region,
                    u,
                    v,
                }
            }
        } else {
            let n = sphere.b() - sphere.c();
            if u + v >= n {
                debug_assert_eq!(u + v, n);
                if n == 0 {
                    todo!()
                } else if v == 0 {
                    if sphere.c() == 0 {
                        Self::base(sphere, region.owner().vertex(1))
                    } else {
                        Self {
                            sphere,
                            region: BaseRegion::new(region.owner(), BaseRegionType::Edge0),
                            u,
                            v: sphere.c(),
                        }
                    }
                } else if u == 0 && sphere.c() == 0 {
                    Self::base(sphere, region.owner().vertex(2))
                } else {
                    Self::adjacent_1(sphere, region.owner(), v)
                }
            } else if u == 0 {
                if v == 0 && sphere.c() == 0 {
                    Self::base(sphere, region.owner().vertex(0))
                } else {
                    let region = if region.owner().owns_edge_2() {
                        BaseRegion::new(region.owner(), BaseRegionType::Edge2)
                    } else {
                        let adj = region.owner().side(2).complement();

                        // With the current base layout, if a base face doesn't own its edge 2,
                        // it will be the edge 0 of the adjacent face that owns it.
                        debug_assert_eq!(adj.index(), 0);
                        BaseRegion::new(adj.inside(), BaseRegionType::Edge0)
                    };
                    Self {
                        sphere,
                        region,
                        u: sphere.c() + v,
                        v: 0,
                    }
                }
            } else if v == 0 {
                Self {
                    sphere,
                    region: BaseRegion::new(region.owner(), BaseRegionType::Edge0),
                    u,
                    v: sphere.c(),
                }
            } else {
                Self {
                    sphere,
                    region,
                    u,
                    v,
                }
            }
        }
    }

    /// Constructs a [`Vertex`] on the edge 1 border of the given base face.
    ///
    /// `r` is the height of the vertex above `face`'s vertex 1. The caller must ensure
    /// `0 < r < b`.
    fn adjacent_1(sphere: SubTriSphere, face: BaseFace, r: u32) -> Self {
        let adj = face.side(1).complement();
        if adj.index() == 0 {
            Self {
                sphere,
                region: BaseRegion::new(adj.inside(), BaseRegionType::Edge0),
                u: sphere.b() - r,
                v: 0,
            }
        } else {
            debug_assert_eq!(adj.index(), 2);
            Self {
                sphere,
                region: BaseRegion::new(adj.inside(), BaseRegionType::Edge2),
                u: r,
                v: sphere.c(),
            }
        }
    }

    /// Gets the [`Vertex`] at the center of the given base face.
    ///
    /// This is only possible when `b == c`.
    fn center(sphere: SubTriSphere, face: BaseFace) -> Self {
        debug_assert_eq!(sphere.b(), sphere.c());
        Self {
            sphere,
            region: BaseRegion::new(face, BaseRegionType::Interior),
            u: 0,
            v: 0,
        }
    }

    /// Constructs a [`Vertex`] for one of the base vertices.
    pub fn base(sphere: SubTriSphere, source: BaseVertex) -> Self {
        if source.index() == 0 {
            Self {
                sphere,
                region: BaseRegion::new(source.owner(), BaseRegionType::Edge0),
                u: 0,
                v: 0,
            }
        } else {
            Self {
                sphere,
                region: BaseRegion::new(source.owner(), BaseRegionType::Edge0),
                u: sphere.b(),
                v: sphere.c(),
            }
        }
    }
}

impl crate::Vertex for Vertex {
    type Face = Face;
    type HalfEdge = HalfEdge;

    fn index(&self) -> usize {
        if self.region.ty().is_edge() {
            self.sphere.base_vertex_index(self.region)
                + self.v as usize * (self.sphere.b() as usize - 1)
                + self.u as usize
        } else if self.sphere.c() < self.sphere.b() {
            let n = (self.sphere.b() - self.sphere.c()) as usize;
            self.sphere.base_vertex_index(self.region)
                + (self.v as usize - 1) * (2 * n - self.v as usize - 2) / 2
                + self.u as usize
        } else {
            self.sphere.base_vertex_index(self.region) + 1
        }
    }

    fn pos(&self) -> [f64; 3] {
        self.sphere
            .project(self.region, self.u as f64, self.v as f64)
    }

    fn degree(&self) -> usize {
        if self.u == 0 && self.v == 0 { 5 } else { 6 }
    }

    fn first_outgoing(&self) -> HalfEdge {
        todo!()
    }
}

/// Represents one "side" of an edge in an [`SubTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfEdge {
    sphere: SubTriSphere,

    /// The [`BaseRegion`] which "owns" this half-edge.
    region: BaseRegion,

    /// The U coordinate of the start vertex, specified in the local coordinate space of `region`.
    start_u: u32,

    /// The V coordinate of the start vertex, specified in the local coordinate space of `region`.
    start_v: u32,

    /// The direction of this half-edge.
    dir: HalfEdgeDir,
}

/// Identifies a possible direction of a half-edge with respect to its [`BaseFace`].
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum HalfEdgeDir {
    /// The +U direction.
    Up = 0,

    /// The +V direction.
    Vp = 1,

    /// The -U/+V direction.
    UnVp = 2,

    /// The -U direction.
    Un = 3,

    /// The -V direction.
    Vn = 4,

    /// The +U/-V direction.
    UpVn = 5,
}

impl crate::HalfEdge for HalfEdge {
    type Vertex = Vertex;
    type Face = Face;

    fn index(&self) -> usize {
        todo!()
    }

    fn inside(&self) -> Face {
        todo!()
    }

    fn start(&self) -> Vertex {
        Vertex::new(self.sphere, self.region, self.start_u, self.start_v)
    }

    fn complement(&self) -> HalfEdge {
        todo!()
    }

    fn prev(&self) -> HalfEdge {
        todo!()
    }

    fn next(&self) -> HalfEdge {
        let mut start_u = self.start_u;
        let mut start_v = self.start_v;
        self.dir.add_to(&mut start_u, &mut start_v);
        let dir = self.dir.rotate_ccw_120();
        HalfEdge {
            sphere: self.sphere,
            region: self.region,
            start_u,
            start_v,
            dir,
        }
    }
}

impl HalfEdgeDir {
    /// Adds this [`HalfEdgeDir`] to the given U and V coordinates.
    pub fn add_to(self, u: &mut u32, v: &mut u32) {
        static TABLE: [(i8, i8); 6] = [
            (1, 0),  // Up
            (0, 1),  // Vp
            (-1, 1), // UnVp
            (-1, 0), // Un
            (0, -1), // Vn
            (1, -1), // UpVn
        ];
        let (d_u, d_v) = TABLE[self as usize];
        *u = u.wrapping_add_signed(d_u as i32);
        *v = v.wrapping_add_signed(d_v as i32);
    }

    /// Rotates this [`HalfEdgeDir`] 120 degrees counter-clockwise.
    pub fn rotate_ccw_120(self) -> Self {
        const TABLE: [HalfEdgeDir; 6] = [
            HalfEdgeDir::UnVp, // Up
            HalfEdgeDir::Un,   // Vp
            HalfEdgeDir::Vn,   // UnVp
            HalfEdgeDir::UpVn, // Un
            HalfEdgeDir::Up,   // Vn
            HalfEdgeDir::Vp,   // UpVn
        ];
        TABLE[self as usize]
    }
}
impl SubTriSphere {
    /// Gets the index of the first face which is owned by the given base region.
    fn base_face_index(&self, region: BaseRegion) -> usize {
        let face = region.owner();
        let num_edge_regions_before = face.num_owned_edges_before()
            + (region.ty() > BaseRegionType::Edge0) as usize;
        let num_interior_regions_before =
            face.index() + (region.ty() > BaseRegionType::Interior) as usize;
        num_edge_regions_before * self.num_faces_per_edge_region()
            + num_interior_regions_before * self.num_faces_per_interior_region()
    }

    /// Gets the index of the first vertex which is owned by the given base region,
    /// not counting the first vertex of the base shape.
    fn base_vertex_index(&self, region: BaseRegion) -> usize {
        let face = region.owner();
        let num_owned_vertices_before = face.num_owned_vertices_before()
            + (face.owns_vertex_1() && region.ty() > BaseRegionType::Edge0) as usize;
        let num_edge_regions_before = face.num_owned_edges_before()
            + (region.ty() > BaseRegionType::Edge0) as usize;
        let num_interior_regions_before =
            face.index() + (region.ty() > BaseRegionType::Interior) as usize;
        num_owned_vertices_before
            + num_edge_regions_before * self.num_vertices_per_edge_region()
            + num_interior_regions_before * self.num_vertices_per_interior_region()
    }

    /// Projects a point in the local coordinate space of this region to a point on the unit
    /// sphere.
    fn project(&self, region: BaseRegion, u: f64, v: f64) -> Vector3 {
        let b = self.b() as f64;
        let c = self.c() as f64;
        let norm_sqr = b * b + b * c + c * c;
        let (face, weights) = if region.ty().is_edge() {
            let w_left = (v * b - u * c) / norm_sqr;
            if region.ty() == BaseRegionType::Edge0 {
                if w_left >= 0.0 {
                    let w_1 = (u * b + (u + v) * c) / norm_sqr;
                    let w_2 = w_left;
                    (region.owner(), [1.0 - w_1 - w_2, w_1, w_2])
                } else {
                    let adj = region.owner().side(0).complement();
                    let w_2 = -w_left;
                    let w_0 = ((u + v) * b + v * c) / norm_sqr;
                    let w_1 = 1.0 - w_0 - w_2;
                    let mut weights = [w_0, w_1, w_2];
                    weights.rotate_right(adj.index());
                    (adj.inside(), weights)
                }
            } else if w_left <= 0.0 {
                let w_1 = -w_left;
                let w_2 = ((u + v) * b + v * c) / norm_sqr;
                (region.owner(), [1.0 - w_1 - w_2, w_1, w_2])
            } else {
                let adj = region.owner().side(2).complement();
                let w_1 = (u * b + (u + v) * c) / norm_sqr;
                let w_2 = w_left;
                let w_0 = 1.0 - w_1 - w_2;
                let mut weights = [w_0, w_1, w_2];
                weights.rotate_right(adj.index());
                (adj.inside(), weights)
            }
        } else {
            let v = v + c;
            let w_1 = (u * b + (u + v) * c) / norm_sqr;
            let w_2 = (v * b - u * c) / norm_sqr;
            (region.owner(), [1.0 - w_1 - w_2, w_1, w_2])
        };
        mat3::apply(
            [
                face.vertex(0).pos(),
                face.vertex(1).pos(),
                face.vertex(2).pos(),
            ],
            interpolate(
                self.base().edge_angle(),
                self.base().edge_cos_angle(),
                weights,
            ),
        )
    }
}

impl SubTriSphere {
    /// Iterates over the faces of this [`SubTriSphere`], starting with the given region.
    fn faces_from(&self, region: BaseRegion) -> FaceIter {
        if region.ty().is_edge() {
            // If `c` is zero, there are no faces in an edge region. This iterator should
            // immediately go to the next region.
            FaceIter {
                sphere: *self,
                region,
                u_0: 0,
                u_0_end: self.b() * (self.c() > 0) as u32,
                v_0: 0,
                boundary_along_v: self.c() == 0,
            }
        } else {
            let n = self.b() - self.c();
            FaceIter {
                sphere: *self,
                region,
                u_0: 0,
                u_0_end: n,
                v_0: 0,
                boundary_along_v: false,
            }
        }
    }
}

/// An iterator over the faces of an [`SubTriSphere`].
#[derive(Clone, Debug)]
pub struct FaceIter {
    sphere: SubTriSphere,
    region: BaseRegion,
    u_0: u32,
    u_0_end: u32,
    v_0: u32,
    boundary_along_v: bool,
}

impl Iterator for FaceIter {
    type Item = Face;
    fn next(&mut self) -> Option<Face> {
        loop {
            if self.u_0 < self.u_0_end {
                let res = Face {
                    sphere: self.sphere,
                    region: self.region,
                    u_0: self.u_0,
                    v_0: self.v_0,
                    boundary_along_v: self.boundary_along_v,
                };
                self.u_0 += 1;
                return Some(res);
            } else if self.region.ty().is_edge() {
                if !self.boundary_along_v {
                    self.u_0 = 1;
                    self.u_0_end = self.sphere.b() + 1;
                    self.boundary_along_v = true;
                    continue;
                } else if self.v_0 + 1 < self.sphere.c() {
                    self.v_0 += 1;
                    self.u_0 = 0;
                    self.u_0_end = self.sphere.b();
                    self.boundary_along_v = false;
                    continue;
                }
            } else if !self.boundary_along_v {
                self.u_0 = 1;
                self.boundary_along_v = true;
                continue;
            } else {
                let n = self.sphere.b() - self.sphere.c();
                if self.v_0 + 1 < n {
                    self.v_0 += 1;
                    self.u_0 = 0;
                    self.u_0_end = n - self.v_0;
                    self.boundary_along_v = false;
                    continue;
                }
            }
            if let Some(region) = self.sphere.base().next_region(self.region) {
                *self = self.sphere.faces_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}

impl SubTriSphere {
    /// Iterates over the vertices of this [`SubTriSphere`], starting with the given region.
    fn vertices_from(&self, region: BaseRegion) -> VertexIter {
        if region.ty().is_edge() {
            VertexIter {
                sphere: *self,
                region,
                u: 1,
                u_end: self.b(),
                v: 0,
            }
        } else if self.c() < self.b() {
            let n = self.b() - self.c();
            VertexIter {
                sphere: *self,
                region,
                u: 1,
                u_end: n - 1,
                v: 1,
            }
        } else {
            // This is a special case. When `b == c` there is exactly one interior vertex,
            // at (0, 0).
            VertexIter {
                sphere: *self,
                region,
                u: 0,
                u_end: 1,
                v: 0,
            }
        }
    }
}

/// An iterator over the vertices of an [`SubTriSphere`].
#[derive(Clone, Debug)]
pub struct VertexIter {
    sphere: SubTriSphere,
    region: BaseRegion,
    u: u32,
    u_end: u32,
    v: u32,
}

impl Iterator for VertexIter {
    type Item = Vertex;
    fn next(&mut self) -> Option<Vertex> {
        loop {
            if self.u < self.u_end {
                let res = Vertex {
                    sphere: self.sphere,
                    region: self.region,
                    u: self.u,
                    v: self.v,
                };
                self.u += 1;
                return Some(res);
            } else if self.region.ty().is_edge() {
                if self.v < self.sphere.c() {
                    self.v += 1;
                    self.u = 1;
                    continue;
                } else if self.u <= self.u_end
                    && self.region.ty() == BaseRegionType::Edge0
                    && self.region.owner().owns_vertex_1()
                {
                    let res = Vertex {
                        sphere: self.sphere,
                        region: self.region,
                        u: self.u,
                        v: self.v,
                    };
                    self.u += 1;
                    return Some(res);
                }
            } else {
                let n = self.sphere.b() - self.sphere.c();
                if self.v + 1 < n {
                    self.v += 1;
                    self.u = 1;
                    self.u_end = n - self.v;
                    continue;
                }
            }
            if let Some(region) = self.sphere.base().next_region(self.region) {
                *self = self.sphere.vertices_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}

/// Identifies a region of a [`BaseTriSphere`] which can contain faces and vertices of a
/// [`SubTriSphere`].
///
/// There is a region corresponding to each face and edge of a [`BaseTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BaseRegion(u8);

impl BaseRegion {
    /// Constructs the [`BaseRegion`] for the given face and type.
    fn new(owner: BaseFace, ty: BaseRegionType) -> Self {
        BaseRegion((owner.0 << 2) | (ty as u8))
    }

    /// Gets the owner of this region.
    fn owner(&self) -> BaseFace {
        BaseFace(self.0 >> 2)
    }

    /// The general type of this [`BaseRegion`].
    fn ty(&self) -> BaseRegionType {
        unsafe { std::mem::transmute(self.0 & 0b11) }
    }
}

/// The general type of a [`BaseRegion`].
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum BaseRegionType {
    /// Corresponds to the first edge of [`BaseRegion::owner`].
    Edge0 = 0,

    /// Corresponds to the interior of [`BaseRegion::owner`].
    Interior = 1,

    /// Corresponds to the third edge of [`BaseRegion::owner`].
    Edge2 = 3,
}

impl BaseRegionType {
    /// Indicates whether this region corresponds to an edge.
    pub fn is_edge(self) -> bool {
        self != Self::Interior
    }
}

impl BaseTriSphere {
    /// Gets the first [`BaseRegion`] of this base shape.
    fn first_region(self) -> BaseRegion {
        BaseRegion::new(self.face(0), BaseRegionType::Edge0)
    }

    /// Gets the [`BaseRegion`] after the given one.
    fn next_region(self, mut region: BaseRegion) -> Option<BaseRegion> {
        if region.ty().is_edge() {
            region.0 += 1;
            if (region.0 >> 2) >= self.last_face_inner() {
                return None;
            }
        } else if region.owner().owns_edge_2() {
            region.0 += 2;
        } else {
            region.0 += 3;
            if (region.0 >> 2) >= self.last_face_inner() {
                return None;
            }
        }
        Some(region)
    }
}

/// Interpolation function for a base face.
///
/// Suppose the vertices of a base face are `v₀`, `v₁`, and `v₂`. The angle between any pair
/// of these vertices is given by `angle` and `cos_angle` is the cosine of that angle.
///
/// The interpolation function is given weights `w₀`, `w₁`, and `w₂` such that `wᵢ ≥ 0` and
/// `w₀ + w₁ + w₂ = 1`. It will compute an approximation of the "spherical" interpolation of the
/// three vertices with respect to the given weights, `r = c₀ v₀ + c₁ v₁ + c₂ v₂` and return the
/// coefficients `[c₀, c₁, c₂]`.
///
/// It is guaranteed that:
///  * `r` is on the unit sphere
///  * If `wᵢ` is `1`, `cᵢ` will also be `1`.
///
/// See [this paper](https://mathweb.ucsd.edu/~sbuss/ResearchWeb/spheremean/paper.pdf) for
/// a precise definition of spherical interpolation (which we only approximate here, for
/// performance reasons).
fn interpolate(angle: f64, cos_angle: f64, weights: [f64; 3]) -> [f64; 3] {
    let [w_0, w_1, w_2] = weights;

    // Compute initial unnormalized coefficients. `cᵢ = (wᵢ * angle).sin()` assures that this
    // behaves like a perfect SLERP (https://en.wikipedia.org/wiki/Slerp) when one of
    // the weights is `0`.
    // TODO: it might be possible to create a faster implementation of `(x * angle).sin()`.
    // This could also make the result portable, since `sin` isn't.
    let c_0 = (w_0 * angle).sin();
    let c_1 = (w_1 * angle).sin();
    let c_2 = (w_2 * angle).sin();

    // Normalize coefficients such that, when applied to the 3 vertices of the base face,
    // the result will be on the unit sphere
    let norm_sqr =
        c_0 * c_0 + c_1 * c_1 + c_2 * c_2 + 2.0 * cos_angle * (c_0 * c_1 + c_1 * c_2 + c_2 * c_0);
    let norm = norm_sqr.sqrt();
    let c_0 = c_0 / norm;
    let c_1 = c_1 / norm;
    let c_2 = c_2 / norm;
    [c_0, c_1, c_2]
}
