//! Contains types related to [`HexSphere`].
use crate::prelude::*;
use crate::proj::BaseTriProjector;
use crate::tri::{self, BaseRegion, BaseRegionType, TriSphere};
use std::num::NonZero;

/// A tessellation of the unit sphere into *mostly* hexagonal [`Face`]s, constructed by grouping
/// triangles of a [`TriSphere`].
///
/// Equivalently, this is a tessellation formed by projecting a
/// [Goldberg polyhedron](https://en.wikipedia.org/wiki/Goldberg_polyhedron) onto
/// the sphere.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HexSphere<Proj> {
    kis: TriSphere<Proj>,
}

impl<Proj: Clone + BaseTriProjector> TriSphere<Proj> {
    /// [Truncates](https://en.wikipedia.org/wiki/Truncation_(geometry)) the vertices of this
    /// [`TriSphere`] to create a [`HexSphere`].
    ///
    /// This will convert all existing faces into hexagons, and create a new face at each vertex.
    pub fn truncate(self) -> HexSphere<Proj> {
        HexSphere {
            kis: self.subdivide_edge(NonZero::new(3).unwrap()),
        }
    }
}

impl<Proj> HexSphere<Proj> {
    /// Attempts to construct a [`HexSphere`] whose [`kis`](HexSphere::kis) is the given
    /// [`TriSphere`].
    ///
    /// For this to succeed, the [`TriSphere`] parameters must satisfy `b % 3 == c % 3`.
    /// This will return [`None`] otherwise.
    pub fn new(kis: TriSphere<Proj>) -> Option<Self> {
        if (kis.b() + 2 * kis.c()) % 3 == 0 {
            Some(Self { kis })
        } else {
            None
        }
    }

    /// Performs the
    /// [kis](https://en.wikipedia.org/wiki/Conway_polyhedron_notation#Original_operations)
    /// operation on this [`HexSphere`].
    ///
    /// This creates a [`Face`](crate::Face) for each [`HalfEdge`] in the [`TriSphere`].
    pub fn kis(self) -> TriSphere<Proj> {
        self.kis
    }

    /// [`TriSphere::num_divisions`] for the dual of this [`HexSphere`].
    const fn dual_num_divisions(&self) -> usize {
        let b = self.kis.b() as usize;
        let c = self.kis.c() as usize;
        let dual_c = (b - c) / 3;
        let dual_b = c + dual_c;
        dual_b * (dual_b + dual_c) + dual_c * dual_c
    }

    /// The number of faces per edge [`BaseRegion`].
    fn num_faces_per_edge_region(&self) -> usize {
        let b = self.kis.b() as usize;
        let c = self.kis.c() as usize;
        num_faces_on_edge(b, c + 1)
    }

    /// The number of faces per interior [`BaseRegion`].
    fn num_faces_per_interior_region(&self) -> usize {
        let b = self.kis.b() as usize;
        let c = self.kis.c() as usize;
        let n = b - c;
        if n == 0 {
            (c % 3 == 0) as usize
        } else {
            num_faces_on_interior(c, n, n)
        }
    }

    /// The number of vertices per edge [`BaseRegion`], assuming that the origin vertex is not
    /// counted.
    fn num_vertices_per_edge_region(&self) -> usize {
        self.kis.num_vertices_per_edge_region() - self.num_faces_per_edge_region()
    }

    /// The number of vertices per interior [`BaseRegion`].
    fn num_vertices_per_interior_region(&self) -> usize {
        self.kis.num_vertices_per_interior_region() - self.num_faces_per_interior_region()
    }
}

/// Counts the number of points `(u, v)` such that `0 < u < b`, `v < d` and `u % 3 == v % 3`.
fn num_faces_on_edge(b: usize, d: usize) -> usize {
    ((b - 1) * d + (b % 3) / 2) / 3
}

#[test]
fn test_num_faces_on_edge() {
    for b in 1..10 {
        for d in 0..10 {
            let mut count = 0;
            for u in 1..b {
                for v in 0..d {
                    if u % 3 == v % 3 {
                        count += 1;
                    }
                }
            }
            assert_eq!(
                num_faces_on_edge(b, d),
                count,
                "failed for b = {}, d = {}",
                b,
                d
            );
        }
    }
}

/// Counts the number of points `(u, v)` such that `0 < u + v < n`, `v < d` and
/// `u % 3 == (v + k) % 3`. Assumes that `n` is a multiple of 3.
fn num_faces_on_interior(k: usize, n: usize, d: usize) -> usize {
    (d - 1) * (2 * n - d - 2) / 6 + ((k % 3 == 0) & (d % 3 != 1)) as usize
}

#[test]
fn test_num_faces_on_interior() {
    for k in 0..2 {
        for n_3 in 1..10 {
            let n = n_3 * 3;
            for d in 1..n {
                let mut count = 0;
                for u in 1..n {
                    for v in 1..d {
                        if u + v < n && u % 3 == (v + k) % 3 {
                            count += 1;
                        }
                    }
                }
                assert_eq!(
                    num_faces_on_interior(k, n, d),
                    count,
                    "failed for k = {}, n = {}, d = {}",
                    k,
                    n,
                    d
                );
            }
        }
    }
}

impl<Proj: Eq + Clone + BaseTriProjector> crate::Sphere for HexSphere<Proj> {
    type Face = Face<Proj>;
    type Vertex = Vertex<Proj>;
    type HalfEdge = HalfEdge<Proj>;

    fn num_faces(&self) -> usize {
        self.num_vertices() / 2 + 2
    }

    fn face(&self, index: usize) -> Face<Proj> {
        // TODO: "Real" implementation with constant-time performance
        self.faces().nth(index).expect("face index out of bounds")
    }

    fn faces(&self) -> impl Iterator<Item = Face<Proj>> {
        FaceIter {
            sphere: self.clone(),
            region: self.kis.base().first_region(),
            u: 0,
            u_end: self.kis.b(),
            v: 0,
        }
    }

    fn face_at(&self, point: [f64; 3]) -> Face<Proj> {
        unsafe { Face::from_kis(self.kis.face_at(point)).unwrap_unchecked() }
    }

    fn num_vertices(&self) -> usize {
        self.dual_num_divisions() * self.kis.base().num_faces()
    }

    fn vertex(&self, index: usize) -> Vertex<Proj> {
        todo!()
    }

    fn vertices(&self) -> impl Iterator<Item = Vertex<Proj>> {
        VertexIter {
            sphere: self.clone(),
            region: self.kis.base().first_region(),
            u: 1,
            u_end: self.kis.b(),
            v: 0,
            skip: false,
        }
    }
}

/// Represents a face on a [`HexSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Face<Proj> {
    center: tri::Vertex<Proj>,
}

impl<Proj: Eq + Clone + BaseTriProjector> Face<Proj> {
    /// Gets the hexsphere face whose central vertex, in the [kised](HexSphere::kis) [`TriSphere`],
    /// is the given [`Vertex`](tri::Vertex).
    ///
    /// This will return [`None`] if the given vertex does not correspond to a face on any
    /// [`HexSphere`].
    ///
    /// # Examples
    /// ```
    /// # use subsphere::prelude::*;
    /// let sphere = subsphere::icosphere().truncate();
    /// let face = sphere.face(0);
    /// let center = face.center();
    /// assert_eq!(subsphere::hex::Face::from_center(center), Some(face));
    /// ```
    pub fn from_center(center: tri::Vertex<Proj>) -> Option<Self> {
        if center.sphere.b() % 3 != center.sphere.c() % 3 {
            return None;
        }
        if center.u % 3 != center.adjusted_v() % 3 {
            return None;
        }
        Some(Self { center })
    }

    /// Gets the hexsphere face which, when [kised](HexSphere::kis), produces the given
    /// triangular [`Face`](tri::Face).
    ///
    /// This will return [`None`] if there is no way to produce the given face by performing a
    /// kis operation.
    ///
    /// # Examples
    /// ```
    /// # use subsphere::prelude::*;
    /// let sphere = subsphere::icosphere().truncate();
    /// let face = sphere.face(0);
    /// let kis_face = face.side(0).kis().inside();
    /// assert_eq!(subsphere::hex::Face::from_kis(kis_face), Some(face));
    /// ```
    pub fn from_kis(kis: tri::Face<Proj>) -> Option<Self> {
        if kis.sphere.b() % 3 != kis.sphere.c() % 3 {
            return None;
        }
        let adj_v = if !kis.region.ty().is_edge() {
            kis.sphere.c()
        } else {
            0
        };
        let [u, v] = match ((kis.u_0 + 2 * (kis.v_0 + adj_v)) % 3, kis.boundary_along_v) {
            (0, _) => [kis.u_0, kis.v_0],
            (1, _) => [kis.u_0, kis.v_0 + 1],
            (_, false) => [kis.u_0 + 1, kis.v_0],
            (_, true) => [kis.u_0 - 1, kis.v_0 + 1],
        };
        Some(unsafe {
            Self::from_center(tri::Vertex::new(kis.sphere, kis.region, u, v)).unwrap_unchecked()
        })
    }

    /// The [`HexSphere`] that this [`Face`] belongs to.
    pub fn sphere(&self) -> HexSphere<Proj> {
        HexSphere {
            kis: self.center.sphere.clone(),
        }
    }

    /// The central vertex of this face on the [kised](HexSphere::kis) [`TriSphere`].
    pub fn center(self) -> tri::Vertex<Proj> {
        self.center
    }

    /// Indicates whether this is a hexagonal face.
    pub fn is_hex(&self) -> bool {
        self.center.as_base().is_none()
    }
}

impl<Proj: Eq + Clone + BaseTriProjector> crate::Face for Face<Proj> {
    type Vertex = Vertex<Proj>;
    type HalfEdge = HalfEdge<Proj>;

    fn index(&self) -> usize {
        let b = self.center.sphere.b() as usize;
        let c = self.center.sphere.c() as usize;
        if self.center.region.ty().is_edge() {
            self.sphere().base_face_index(self.center.region)
                + num_faces_on_edge(b, self.center.v as usize)
                + (self.center.u as usize).div_ceil(3)
        } else if c < b {
            self.sphere().base_face_index(self.center.region)
                + num_faces_on_interior(c, b - c, self.center.v as usize)
                + (self.center.u as usize).div_ceil(3)
        } else {
            self.sphere().base_face_index(self.center.region) + 1
        }
    }

    fn area(&self) -> f64 {
        crate::util::poly_area(self.vertices().map(|v| v.pos()))
    }

    fn num_sides(&self) -> usize {
        if let Some(base) = self.center.as_base() {
            base.degree()
        } else {
            6
        }
    }

    fn side(&self, index: usize) -> HalfEdge<Proj> {
        HalfEdge {
            kis: self.center.outgoing(index).next(),
        }
    }
}

/// Represents a vertex on a [`HexSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Vertex<Proj> {
    kis: tri::Vertex<Proj>,
}

impl<Proj> Vertex<Proj> {
    /// Constructs a [`Vertex`] from the given [`tri::Vertex`].
    fn new(kis: tri::Vertex<Proj>) -> Self {
        debug_assert_ne!(kis.u % 3, kis.adjusted_v() % 3);
        Self { kis }
    }
}

impl<Proj> tri::Vertex<Proj> {
    /// The `v` coordinate of this vertex, measured from the base vertex of its region, rather
    /// than the local origin.
    ///
    /// The combination of `u` and `adjusted_v` can be used to tell whether a [`tri::Vertex`]
    /// corresponds to a face or a vertex on a [`HexSphere`].
    fn adjusted_v(&self) -> u32 {
        let mut v = self.v;
        if !self.region.ty().is_edge() {
            v += self.sphere.c();
        }
        v
    }
}

impl<Proj: Clone> Vertex<Proj> {
    /// The [`HexSphere`] that this [`Vertex`] belongs to.
    pub fn sphere(&self) -> HexSphere<Proj> {
        HexSphere {
            kis: self.kis.sphere.clone(),
        }
    }

    /// The corresponding vertex on the [kised](HexSphere::kis) [`TriSphere`].
    pub fn kis(self) -> tri::Vertex<Proj> {
        self.kis
    }
}

impl<Proj: Eq + Clone + BaseTriProjector> crate::Vertex for Vertex<Proj> {
    type Face = Face<Proj>;
    type HalfEdge = HalfEdge<Proj>;

    fn index(&self) -> usize {
        let b = self.kis.sphere.b() as usize;
        let c = self.kis.sphere.c() as usize;
        if self.kis.region.ty().is_edge() {
            self.sphere().base_vertex_index(self.kis.region) + self.kis.v as usize * (b - 1)
                - num_faces_on_edge(b, self.kis.v as usize)
                + (self.kis.u as usize * 2 - 1 - (self.kis.v % 3 == 1) as usize) / 3
        } else if c < b {
            let n = b - c;
            self.sphere().base_vertex_index(self.kis.region)
                + (self.kis.v as usize - 1) * (2 * n - self.kis.v as usize - 2) / 2
                - num_faces_on_interior(c, n, self.kis.v as usize)
                + (self.kis.u as usize * 2 - 1 - ((self.kis.v as usize + c) % 3 == 1) as usize) / 3
        } else {
            self.sphere().base_vertex_index(self.kis.region)
        }
    }

    fn pos(&self) -> [f64; 3] {
        self.kis.pos()
    }

    fn degree(&self) -> usize {
        3
    }

    fn outgoing(&self, index: usize) -> HalfEdge<Proj> {
        assert!(index < 3, "index out of bounds");
        let phase = (self.kis.u + 2 * self.kis.adjusted_v()) % 3;
        debug_assert_ne!(phase, 0);
        HalfEdge {
            kis: tri::HalfEdge::new(
                self.kis.sphere.clone(),
                self.kis.region,
                self.kis.u,
                self.kis.v,
                tri::HalfEdgeDir::from_index(2 * index + phase as usize - 1),
            ),
        }
    }
}

/// Represents one "side" of an edge on a [`HexSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfEdge<Proj> {
    kis: tri::HalfEdge<Proj>,
}

impl<Proj> HalfEdge<Proj> {
    /// The corresponding half-edge on the [kised](HexSphere::kis) [`TriSphere`].
    pub fn kis(self) -> tri::HalfEdge<Proj> {
        self.kis
    }
}

impl<Proj: Eq + Clone + BaseTriProjector> crate::HalfEdge for HalfEdge<Proj> {
    type Face = Face<Proj>;
    type Vertex = Vertex<Proj>;

    fn side_index(&self) -> usize {
        self.kis.prev().outgoing_index()
    }

    fn length(&self) -> f64 {
        self.kis.length()
    }

    fn angle(&self) -> f64 {
        let v_a = self.prev().start().pos();
        let v_b = self.start().pos();
        let v_c = self.next().start().pos();
        crate::util::angle(v_a, v_b, v_c)
    }

    fn inside(&self) -> Face<Proj> {
        Face {
            center: self.kis.prev().start(),
        }
    }

    fn start(&self) -> Vertex<Proj> {
        Vertex {
            kis: self.kis.start(),
        }
    }

    fn complement(&self) -> Self {
        Self {
            kis: self.kis.complement(),
        }
    }

    fn prev(&self) -> Self {
        Self {
            kis: self.kis.prev().complement().prev(),
        }
    }

    fn next(&self) -> Self {
        Self {
            kis: self.kis.next().complement().next(),
        }
    }
}

impl<Proj> HexSphere<Proj> {
    /// Gets the index of the first face which is owned by the given base region,
    /// not counting the one corresponding to the first vertex of the base shape.
    fn base_face_index(&self, region: BaseRegion) -> usize {
        let face = region.owner();
        let num_owned_vertices_before = face.num_owned_vertices_before()
            + (face.owns_vertex_1() && region.ty() > BaseRegionType::Edge0) as usize;
        let num_edge_regions_before =
            face.num_owned_edges_before() + (region.ty() > BaseRegionType::Edge0) as usize;
        let num_interior_regions_before =
            face.index() + (region.ty() > BaseRegionType::Interior) as usize;
        num_owned_vertices_before
            + num_edge_regions_before * self.num_faces_per_edge_region()
            + num_interior_regions_before * self.num_faces_per_interior_region()
    }

    /// Gets the index of the first vertex which is owned by the given base region.
    fn base_vertex_index(&self, region: BaseRegion) -> usize {
        let face = region.owner();
        let num_edge_regions_before =
            face.num_owned_edges_before() + (region.ty() > BaseRegionType::Edge0) as usize;
        let num_interior_regions_before =
            face.index() + (region.ty() > BaseRegionType::Interior) as usize;
        num_edge_regions_before * self.num_vertices_per_edge_region()
            + num_interior_regions_before * self.num_vertices_per_interior_region()
    }
}

impl<Proj> HexSphere<Proj> {
    /// Iterates over the faces of this [`HexSphere`], starting with the given region.
    fn faces_from(self, region: BaseRegion) -> FaceIter<Proj> {
        let b = self.kis.b();
        let c = self.kis.c();
        if region.ty().is_edge() {
            FaceIter {
                sphere: self,
                region,
                u: 3,
                u_end: b,
                v: 0,
            }
        } else if c < b {
            let n = b - c;
            FaceIter {
                sphere: self,
                region,
                u: 1 + c % 3,
                u_end: n - 1,
                v: 1,
            }
        } else {
            // This is a special case. When `b == c` and `c % 3 == 0`, there is exactly one
            // interior face, at (0, 0).
            FaceIter {
                sphere: self,
                region,
                u: c % 3,
                u_end: 1,
                v: 0,
            }
        }
    }
}

/// An iterator over the faces of an [`TriSphere`].
#[derive(Clone, Debug)]
pub struct FaceIter<Proj> {
    sphere: HexSphere<Proj>,
    region: BaseRegion,
    u: u32,
    u_end: u32,
    v: u32,
}

impl<Proj: Eq + Clone + BaseTriProjector> Iterator for FaceIter<Proj> {
    type Item = Face<Proj>;
    fn next(&mut self) -> Option<Face<Proj>> {
        loop {
            if self.u < self.u_end {
                let res = unsafe {
                    Face::from_center(tri::Vertex {
                        sphere: self.sphere.kis.clone(),
                        region: self.region,
                        u: self.u,
                        v: self.v,
                    })
                    .unwrap_unchecked()
                };
                self.u += 3;
                return Some(res);
            } else if self.region.ty().is_edge() {
                if self.v < self.sphere.kis.c() {
                    self.v += 1;
                    self.u = 1 + (self.v + 2) % 3;
                    continue;
                } else if self.u <= self.u_end
                    && self.region.ty() == BaseRegionType::Edge0
                    && self.region.owner().owns_vertex_1()
                {
                    let res = unsafe {
                        Face::from_center(tri::Vertex {
                            sphere: self.sphere.kis.clone(),
                            region: self.region,
                            u: self.u,
                            v: self.v,
                        })
                        .unwrap_unchecked()
                    };
                    self.u += 3;
                    return Some(res);
                }
            } else {
                let n = self.sphere.kis.b() - self.sphere.kis.c();
                if self.v + 1 < n {
                    self.v += 1;
                    self.u = 1 + (self.v + self.sphere.kis.c() + 2) % 3;
                    self.u_end = n - self.v;
                    continue;
                }
            }
            if let Some(region) = self.sphere.kis.base().next_region(self.region) {
                *self = self.sphere.clone().faces_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}

impl<Proj: Eq + Clone + BaseTriProjector> HexSphere<Proj> {
    /// Iterates over the vertices of this [`HexSphere`], starting with the given region.
    fn vertices_from(self, region: BaseRegion) -> VertexIter<Proj> {
        let b = self.kis.b();
        let c = self.kis.c();
        if region.ty().is_edge() {
            VertexIter {
                sphere: self,
                region,
                u: 1,
                u_end: b,
                v: 0,
                skip: false,
            }
        } else if c < b {
            let n = b - c;
            VertexIter {
                sphere: self,
                region,
                u: 1 + (c % 3 == 0) as u32,
                u_end: n.saturating_sub(1),
                v: 1,
                skip: c % 3 == 1,
            }
        } else {
            // This is a special case. When `b == c` and `c % 3 > 0`, there is exactly one
            // interior vertex, at (0, 0).
            VertexIter {
                sphere: self,
                region,
                u: (c % 3 == 0) as u32,
                u_end: 1,
                v: 0,
                skip: c % 3 == 2,
            }
        }
    }
}

/// An iterator over the vertices of an [`HexSphere`].
#[derive(Clone, Debug)]
pub struct VertexIter<Proj> {
    sphere: HexSphere<Proj>,
    region: BaseRegion,
    u: u32,
    u_end: u32,
    v: u32,
    skip: bool,
}

impl<Proj: Eq + Clone + BaseTriProjector> Iterator for VertexIter<Proj> {
    type Item = Vertex<Proj>;
    fn next(&mut self) -> Option<Vertex<Proj>> {
        loop {
            if self.u < self.u_end {
                let res = Vertex::new(tri::Vertex {
                    sphere: self.sphere.kis.clone(),
                    region: self.region,
                    u: self.u,
                    v: self.v,
                });
                self.u += 1;
                self.u += self.skip as u32;
                self.skip = !self.skip;
                return Some(res);
            } else if self.region.ty().is_edge() {
                if self.v < self.sphere.kis.c() {
                    self.v += 1;
                    (self.u, self.skip) = match self.v % 3 {
                        0 => (1, false),
                        1 => (2, false),
                        _ => (1, true),
                    };
                    continue;
                }
            } else {
                let n = self.sphere.kis.b() - self.sphere.kis.c();
                if self.v + 1 < n {
                    self.v += 1;
                    self.u_end = n - self.v;
                    (self.u, self.skip) = match (self.v + self.sphere.kis.c()) % 3 {
                        0 => (1, false),
                        1 => (2, false),
                        _ => (1, true),
                    };
                    continue;
                }
            }
            if let Some(region) = self.sphere.kis.base().next_region(self.region) {
                *self = self.sphere.clone().vertices_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}
