//! Contains types related to [`SubTriSphere`].
use crate::Face as FaceExt;
use crate::HalfEdge as HalfEdgeExt;
use crate::Sphere as SphereExt;
use crate::Vertex as VertexExt;
use crate::basetri::BaseTriSphere;
use crate::basetri::Face as BaseFace;
use crate::basetri::HalfEdge as BaseHalfEdge;
use crate::basetri::Vertex as BaseVertex;
use crate::math::{self, Scalar, Vector3, mat3, vec3};
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
        let eval = self.sphere.eval().region(self.region);
        let p_0 = eval.project(u_0 as f64, v_0 as f64);
        let p_1 = eval.project(u_1 as f64, v_1 as f64);
        let p_2 = eval.project(u_2 as f64, v_2 as f64);
        math::sphere_tri_area([p_0, p_1, p_2])
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
            let edge = if region.ty() == BaseRegionType::Edge0 {
                region.owner().side(0)
            } else {
                region.owner().side(2).complement()
            };
            if v == 0 {
                Self::on_edge_u_boundary(sphere, edge, u)
            } else if u == 0 {
                Self::on_edge_v_boundary(sphere, edge, v)
            } else if v == sphere.c() {
                Self::on_edge_u_boundary(sphere, edge.complement(), sphere.b() - u)
            } else if u == sphere.b() {
                Self::on_edge_v_boundary(sphere, edge.complement(), sphere.c() - v)
            } else {
                Self {
                    sphere,
                    region,
                    u,
                    v,
                }
            }
        } else if v == 0 {
            Self::on_interior_boundary(sphere, region.owner().side(0), u)
        } else {
            let n = sphere.b() - sphere.c();
            if u == 0 {
                Self::on_interior_boundary(sphere, region.owner().side(2), n - v)
            } else if u + v >= n {
                debug_assert_eq!(u + v, n);
                Self::on_interior_boundary(sphere, region.owner().side(1), v)
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

    /// Constructs a [`Vertex`] on the bottom (V = 0) boundary of an interior region.
    fn on_interior_boundary(sphere: SubTriSphere, bottom: BaseHalfEdge, u: u32) -> Self {
        Self::on_edge_u_boundary(sphere, bottom.complement(), sphere.b() - u)
    }

    /// Constructs a [`Vertex`] on the V (U = 0) boundary of an edge region.
    fn on_edge_v_boundary(sphere: SubTriSphere, edge: BaseHalfEdge, v: u32) -> Self {
        Self::on_edge_u_boundary(sphere, edge.prev().complement(), v)
    }

    /// Constructs a [`Vertex`] on the U (V = 0) boundary of an edge region.
    fn on_edge_u_boundary(sphere: SubTriSphere, edge: BaseHalfEdge, u: u32) -> Self {
        if u == sphere.b() {
            // In the special case where `b == c`, it's the interior region which owns the
            // central vertex.
            if sphere.b() == sphere.c() {
                Self::center(sphere, edge.complement().inside())
            } else {
                Self::on_edge_u_boundary_exclusive(
                    sphere,
                    edge.complement().prev().complement(),
                    sphere.c(),
                )
            }
        } else {
            Self::on_edge_u_boundary_exclusive(sphere, edge, u)
        }
    }

    /// Constructs a [`Vertex`] on the U (V = 0) boundary of an edge region when it is guaranteed
    /// that `u` is not `b`.
    fn on_edge_u_boundary_exclusive(sphere: SubTriSphere, edge: BaseHalfEdge, u: u32) -> Self {
        if u == 0 {
            Self::base(sphere, edge.start())
        } else {
            let comp = match edge.index() {
                0 => {
                    return Self {
                        sphere,
                        region: BaseRegion::new(edge.inside(), BaseRegionType::Edge0),
                        u,
                        v: 0,
                    };
                }
                1 => edge.complement(),
                _ => {
                    if edge.inside().owns_edge_2() {
                        return Self {
                            sphere,
                            region: BaseRegion::new(edge.inside(), BaseRegionType::Edge2),
                            u: sphere.b() - u,
                            v: sphere.c(),
                        };
                    } else {
                        edge.complement()
                    }
                }
            };
            if comp.index() == 0 {
                Self {
                    sphere,
                    region: BaseRegion::new(comp.inside(), BaseRegionType::Edge0),
                    u: sphere.b() - u,
                    v: sphere.c(),
                }
            } else {
                debug_assert_eq!(comp.index(), 2);
                debug_assert!(comp.inside().owns_edge_2());
                Self {
                    sphere,
                    region: BaseRegion::new(comp.inside(), BaseRegionType::Edge2),
                    u,
                    v: 0,
                }
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
            .eval()
            .region(self.region)
            .project(self.u as f64, self.v as f64)
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
        let num_edge_regions_before =
            face.num_owned_edges_before() + (region.ty() > BaseRegionType::Edge0) as usize;
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
        let num_edge_regions_before =
            face.num_owned_edges_before() + (region.ty() > BaseRegionType::Edge0) as usize;
        let num_interior_regions_before =
            face.index() + (region.ty() > BaseRegionType::Interior) as usize;
        num_owned_vertices_before
            + num_edge_regions_before * self.num_vertices_per_edge_region()
            + num_interior_regions_before * self.num_vertices_per_interior_region()
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

// TODO: For best performance, `SphereEvaluator` should be cached when performing operations
// on the same sphere multiple times. This will require the user's help, so we should figure out a
// way to expose this to them.

/// Encapsulates the information needed to convert between local coordinates and positions
/// on a [`SubTriSphere`].
struct SphereEvaluator {
    /// Suppose the vertices of a base face are `v₀`, `v₁`, and `v₂`. This provides a set of
    /// coefficients `c₀`, `c₁`, and `c₂` such that `p₀ = c₀ v₀ + c₁ v₁ + c₂ v₂` is the
    /// position of the origin of the interior region corresponding to the face. Due to
    /// symmetry, the other points of the interior region can be computed by rotating these
    /// coefficients:
    ///  * `p₀ = c₀ v₀ + c₁ v₁ + c₂ v₂`
    ///  * `p₁ = c₀ v₁ + c₁ v₂ + c₂ v₀`
    ///  * `p₂ = c₀ v₂ + c₁ v₀ + c₂ v₁`
    p_coeffs: [Scalar; 3],

    /// The angle, in radians, subtended by a step of length `1` in the local coordinate space of
    /// a region.
    step_angle: Scalar,

    /// The angle, in radians, subtended by an edge corresponding to the `b` parameter of the
    /// [`SubTriSphere`].
    full_angle_b: Scalar,

    /// The angle, in radians, subtended by an edge corresponding to the `c` parameter of the
    /// [`SubTriSphere`].
    full_angle_c: Scalar,
}

impl SubTriSphere {
    /// Constructs a [`SphereEvaluator`] for this [`SubTriSphere`].
    fn eval(&self) -> SphereEvaluator {
        let b = self.b() as f64;
        let c = self.c() as f64;
        if c == 0.0 {
            let full_angle_b = self.base().edge_angle();
            return SphereEvaluator {
                p_coeffs: [1.0, 0.0, 0.0],
                step_angle: full_angle_b / b,
                full_angle_b,
                full_angle_c: 0.0,
            };
        }

        // Compute initial estimate of the first interior point `p₀ = c₀ v₀ + c₁ v₁ + c₂ v₂`
        let w_total = b * b + b * c + c * c;
        let w_0 = b * b / w_total;
        let w_1 = c * c / w_total;
        let w_2 = b * c / w_total;
        let angle = self.base().edge_angle();
        let cos_angle = self.base().edge_cos_angle();
        let c_0 = (w_0 * angle).sin();
        let c_1 = (w_1 * angle).sin();
        let c_2 = (w_2 * angle).sin();
        let mut p_coeffs = normalize_coeffs(cos_angle, [c_0, c_1, c_2]);

        // For continuity between regions, it is very important that `p₂`, `p₀`,
        // and `v₀` all lie on the same plane, and that `angle(p₀, v₀) / angle(p₂, v₀) = c / b`.
        // To ensure this, we will refine our initial estimate above.
        let mut cos_full_angle_b = p_coeffs[1] + (p_coeffs[0] + p_coeffs[2]) * cos_angle;
        let mut sin_full_angle_b = (1.0 - cos_full_angle_b * cos_full_angle_b).sqrt();
        let mut full_angle_b = cos_full_angle_b.acos();

        // TODO: Find a better algorithm that converges faster, or a closed form solution.
        for _ in 0..30 {
            let t = c / b;
            p_coeffs = vec3::add(
                [
                    ((1.0 - t) * full_angle_b).sin() / sin_full_angle_b,
                    0.0,
                    0.0,
                ],
                vec3::mul(
                    [p_coeffs[1], p_coeffs[2], p_coeffs[0]],
                    (t * full_angle_b).sin() / sin_full_angle_b,
                ),
            );
            cos_full_angle_b = p_coeffs[1] + (p_coeffs[0] + p_coeffs[2]) * cos_angle;
            sin_full_angle_b = (1.0 - cos_full_angle_b * cos_full_angle_b).sqrt();
            full_angle_b = cos_full_angle_b.acos();
        }

        // Compute angles of region edges
        let full_angle_c = full_angle_b * c / b;
        let step_angle = full_angle_b / b;
        SphereEvaluator {
            p_coeffs,
            step_angle,
            full_angle_b,
            full_angle_c,
        }
    }
}

/// Given a set of coefficients which are to be applied to a set of vertices, normalizes them
/// such that the resulting point is on the unit sphere.
///
/// `cos_angle` is the cosine of the angle (or, equivalently, the dot product) between any
/// two of the vertices.
fn normalize_coeffs(cos_angle: Scalar, coeffs: [Scalar; 3]) -> [Scalar; 3] {
    let norm_sqr = vec3::dot(coeffs, coeffs)
        + 2.0 * cos_angle * (coeffs[0] * coeffs[1] + coeffs[1] * coeffs[2] + coeffs[2] * coeffs[0]);
    let norm = norm_sqr.sqrt();
    vec3::mul(coeffs, 1.0 / norm)
}

impl SphereEvaluator {
    /// Constructs an [`RegionEvaluator`] for a particular region of the sphere.
    pub fn region(&self, region: BaseRegion) -> RegionEvaluator {
        match region.ty() {
            BaseRegionType::Edge0 => self.edge(region.owner().side(0)).into(),
            BaseRegionType::Interior => self.interior(region.owner()).into(),
            BaseRegionType::Edge2 => self.edge(region.owner().side(2).complement()).into(),
        }
    }

    /// Constructs an [`InteriorEvaluator`] for a particular interior region of the sphere.
    pub fn interior(&self, face: BaseFace) -> InteriorEvaluator {
        let v_0 = face.vertex(0).pos();
        let v_1 = face.vertex(1).pos();
        let v_2 = face.vertex(2).pos();
        InteriorEvaluator {
            p: [
                mat3::apply([v_0, v_1, v_2], self.p_coeffs),
                mat3::apply([v_1, v_2, v_0], self.p_coeffs),
                mat3::apply([v_2, v_0, v_1], self.p_coeffs),
            ],
            step_angle: self.step_angle,
            full_angle: safe_angle(self.full_angle_b - self.full_angle_c),
        }
    }

    /// Constructs an [`EdgeEvaluator`] for a particular edge region of the sphere.
    pub fn edge(&self, edge: BaseHalfEdge) -> EdgeEvaluator {
        let comp = edge.complement();
        let v_0 = edge.start().pos();
        let v_1 = comp.start().pos();
        let v_left = edge.prev().start().pos();
        let v_right = comp.prev().start().pos();
        EdgeEvaluator {
            p_0_0: v_0,
            p_1_1: v_1,
            p_0_1: mat3::apply([v_0, v_1, v_left], self.p_coeffs),
            p_1_0: mat3::apply([v_1, v_0, v_right], self.p_coeffs),
            step_angle: self.step_angle,
            full_angle_u: self.full_angle_b,
            full_angle_v: safe_angle(self.full_angle_c),
        }
    }
}

/// Encapsulates the information needed to convert between local coordinates on a [`BaseRegion`]
/// and positions on a [`SubTriSphere`].
enum RegionEvaluator {
    Interior(InteriorEvaluator),
    Edge(EdgeEvaluator),
}

impl RegionEvaluator {
    /// Projects a point in the local coordinate space of this region to a point on the sphere.
    pub fn project(&self, u: Scalar, v: Scalar) -> Vector3 {
        match self {
            RegionEvaluator::Interior(eval) => eval.project(u, v),
            RegionEvaluator::Edge(eval) => eval.project(u, v),
        }
    }
}

impl From<InteriorEvaluator> for RegionEvaluator {
    fn from(eval: InteriorEvaluator) -> Self {
        RegionEvaluator::Interior(eval)
    }
}

impl From<EdgeEvaluator> for RegionEvaluator {
    fn from(eval: EdgeEvaluator) -> Self {
        RegionEvaluator::Edge(eval)
    }
}

/// Encapsulates the information needed to convert between local coordinates on an interior
/// [`BaseRegion`] and positions on a [`SubTriSphere`].
struct InteriorEvaluator {
    /// The positions of the endpoints of the region on the sphere.
    p: [Vector3; 3],

    /// The angle, in radians, subtended by a step of length `1` in the local coordinate space of
    /// this region.
    step_angle: Scalar,

    /// The angle, in radians, between any two endpoints of the region.
    full_angle: Scalar,
}

impl InteriorEvaluator {
    /// Projects a point in the local coordinate space of this region to a point on the sphere.
    pub fn project(&self, u: Scalar, v: Scalar) -> Vector3 {
        let c_1 = (u * self.step_angle).sin();
        let c_2 = (v * self.step_angle).sin();
        let c_0 = (self.full_angle - u * self.step_angle - v * self.step_angle).sin();
        vec3::normalize(mat3::apply(self.p, [c_0, c_1, c_2]))
    }
}

/// Encapsulates the information needed to convert between local coordinates on an edge
/// [`BaseRegion`] and positions on a [`SubTriSphere`].
struct EdgeEvaluator {
    /// The position of the origin of the region on the sphere.
    p_0_0: Vector3,

    /// The position of the far (U = 1, V = 1) corner of the region on the sphere.
    p_1_1: Vector3,

    /// The position of the left (U = 0, V = 1) corner of the region on the sphere.
    p_0_1: Vector3,

    /// The position of the right (U = 1, V = 0) corner of the region on the sphere.
    p_1_0: Vector3,

    /// The angle, in radians, subtended by a step of length `1` in the local coordinate space of
    /// this region.
    step_angle: Scalar,

    /// The angle, in radians, subtended by the U edge of this region.
    full_angle_u: Scalar,

    /// The angle, in radians, subtended by the V edge of this region.
    full_angle_v: Scalar,
}

impl EdgeEvaluator {
    /// Projects a point in the local coordinate space of this region to a point on the sphere.
    pub fn project(&self, u: Scalar, v: Scalar) -> Vector3 {
        let c_u_0 = (self.full_angle_u - u * self.step_angle).sin();
        let c_u_1 = (u * self.step_angle).sin();
        let c_v_0 = (self.full_angle_v - v * self.step_angle).sin();
        let c_v_1 = (v * self.step_angle).sin();
        vec3::normalize(vec3::add(
            vec3::add(
                vec3::mul(self.p_0_0, c_u_0 * c_v_0),
                vec3::mul(self.p_1_0, c_u_1 * c_v_0),
            ),
            vec3::add(
                vec3::mul(self.p_0_1, c_u_0 * c_v_1),
                vec3::mul(self.p_1_1, c_u_1 * c_v_1),
            ),
        ))
    }
}

/// Checks if `angle` is zero. If so, replaces it with an arbitrary positive value. Otherwise
/// returns it.
fn safe_angle(angle: Scalar) -> Scalar {
    if angle == 0.0 { 1.0 } else { angle }
}
