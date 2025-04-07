//! Contains types related to [`HexSphere`].
use crate::Face as FaceExt;
use crate::Vertex as VertexExt;
use crate::HalfEdge as HalfEdgeExt;
use crate::tri::{self, BaseRegion, BaseRegionType, TriSphere};
use std::num::NonZero;

/// A tessellation of the unit sphere into *mostly* hexagonal [`Face`]s, constructed by grouping
/// triangles of a [`TriSphere`].
/// 
/// Equivalently, this is a tessellation formed by projecting a
/// [Goldberg polyhedron](https://en.wikipedia.org/wiki/Goldberg_polyhedron) onto
/// the sphere.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HexSphere {
    kis: TriSphere,
}

impl TriSphere {
    /// [Truncates](https://en.wikipedia.org/wiki/Truncation_(geometry)) the vertices of this
    /// [`TriSphere`] to create a [`HexSphere`].
    ///
    /// This will convert all existing faces into hexagons, and create a new face at each vertex.
    pub fn truncate(self) -> HexSphere {
        HexSphere {
            kis: self.subdivide_edge(NonZero::new(3).unwrap()),
        }
    }
}

impl HexSphere {
    /// Attempts to construct a [`HexSphere`] whose [`kis`](HexSphere::kis) is the given
    /// [`TriSphere`].
    ///
    /// For this to succeed, the [`TriSphere`] parameters must satisfy `b % 3 == c % 3`.
    /// This will return [`None`] otherwise.
    pub fn new(kis: TriSphere) -> Option<Self> {
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
    pub fn kis(self) -> TriSphere {
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
        // TODO: There's a bug here. The number of faces also depends on the offset of the
        // interior region (i.e. `c`).
        let n_3 = (self.kis.b() - self.kis.c()) as isize / 3;
        (n_3 * n_3 + ((n_3 - 1) * (n_3 - 2) / 2)) as usize
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

/// Counts the number of points `(u, v)` such that `0 < u + v < n`, `v < d` and `u % 3 == v % 3`.
/// Assumes that `n` is a multiple of 3.
fn num_faces_on_interior(n: usize, d: usize) -> usize {
    (d - 1) * (2 * n - d - 2) / 6 + (d % 3 != 1) as usize
}

#[test]
fn test_num_faces_on_interior() {
    for n_3 in 1..10 {
        let n = n_3 * 3;
        for d in 1..n {
            let mut count = 0;
            for u in 1..n {
                for v in 1..d {
                    if u + v < n && u % 3 == v % 3 {
                        count += 1;
                    }
                }
            }
            assert_eq!(
                num_faces_on_interior(n, d),
                count,
                "failed for n = {}, d = {}",
                n,
                d
            );
        }
    }
}

impl crate::Sphere for HexSphere {
    type Face = Face;
    type Vertex = Vertex;
    type HalfEdge = HalfEdge;

    fn num_faces(&self) -> usize {
        self.num_vertices() / 2 + 2
    }

    fn face(&self, index: usize) -> Face {
        todo!()
    }

    fn faces(&self) -> impl Iterator<Item = Face> {
        FaceIter {
            sphere: *self,
            region: self.kis.base().first_region(),
            u: 0,
            u_end: self.kis.b(),
            v: 0,
        }
    }

    fn num_vertices(&self) -> usize {
        self.dual_num_divisions() * self.kis.base().num_faces()
    }

    fn vertex(&self, index: usize) -> Vertex {
        todo!()
    }

    fn vertices(&self) -> impl Iterator<Item = Vertex> {
        VertexIter {
            sphere: *self,
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
pub struct Face {
    center: tri::Vertex,
}

impl Face {
    /// Constructs a [`Face`] from the given central vertex.
    fn new(center: tri::Vertex) -> Self {
        debug_assert!(center.u % 3 == center.v % 3);
        Self { center }
    }

    /// The [`HexSphere`] that this [`Face`] belongs to.
    pub fn sphere(self) -> HexSphere {
        HexSphere {
            kis: self.center.sphere,
        }
    }

    /// Indicates whether this is a hexagonal face.
    pub fn is_hex(&self) -> bool {
        self.center.as_base().is_none()
    }
}

impl crate::Face for Face {
    type Vertex = Vertex;
    type HalfEdge = HalfEdge;

    fn index(&self) -> usize {
        if self.center.region.ty().is_edge() {
            self.sphere().base_face_index(self.center.region)
                + num_faces_on_edge(self.center.sphere.b() as usize, self.center.v as usize)
                + (self.center.u as usize + 2) / 3
        } else if self.center.sphere.c() < self.center.sphere.b() {
            let n = (self.center.sphere.b() - self.center.sphere.c()) as usize;
            self.sphere().base_face_index(self.center.region)
                + num_faces_on_interior(n, self.center.v as usize)
                + (self.center.u as usize + 2) / 3
        } else {
            self.sphere().base_face_index(self.center.region) + 1
        }
    }

    fn num_sides(&self) -> usize {
        if let Some(base) = self.center.as_base() {
            base.degree()
        } else {
            6
        }
    }

    fn side(&self, index: usize) -> HalfEdge {
        HalfEdge {
            kis: self.center.outgoing(index).next()
        }
    }
}

/// Represents a vertex on a [`HexSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Vertex {
    kis: tri::Vertex,
}

impl Vertex {
    /// Constructs a [`Vertex`] from the given [`subtri::Vertex`].
    fn new(kis: tri::Vertex) -> Self {
        debug_assert!((kis.u + kis.v * 2) % 3 != 0);
        Self { kis }
    }

    /// The [`HexSphere`] that this [`Vertex`] belongs to.
    pub fn sphere(self) -> HexSphere {
        HexSphere {
            kis: self.kis.sphere,
        }
    }
}

impl crate::Vertex for Vertex {
    type Face = Face;
    type HalfEdge = HalfEdge;

    fn index(&self) -> usize {
        if self.kis.region.ty().is_edge() {
            self.sphere().base_vertex_index(self.kis.region)
                + self.kis.v as usize * (self.kis.sphere.b() as usize - 1)
                - num_faces_on_edge(self.kis.sphere.b() as usize, self.kis.v as usize)
                + (self.kis.u as usize * 2 - 1 - (self.kis.v % 3 == 1) as usize) / 3
        } else {
            let n = (self.kis.sphere.b() - self.kis.sphere.c()) as usize;
            self.sphere().base_vertex_index(self.kis.region)
                + (self.kis.v as usize - 1) * (2 * n - self.kis.v as usize - 2) / 2
                - num_faces_on_interior(n, self.kis.v as usize)
                + (self.kis.u as usize * 2 - 1 - (self.kis.v % 3 == 1) as usize) / 3
        }
    }

    fn pos(&self) -> [f64; 3] {
        self.kis.pos()
    }

    fn degree(&self) -> usize {
        3
    }
    
    fn outgoing(&self, index: usize) -> HalfEdge {
        todo!()
    }
}

/// Represents one "side" of an edge on a [`HexSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfEdge {
    kis: tri::HalfEdge,
}

impl crate::HalfEdge for HalfEdge {
    type Face = Face;
    type Vertex = Vertex;

    fn side_index(&self) -> usize {
        todo!()
    }

    fn inside(&self) -> Face {
        todo!()
    }

    fn start(&self) -> Vertex {
        Vertex {
            kis: self.kis.start(),
        }
    }

    fn complement(&self) -> Self {
        todo!()
    }

    fn prev(&self) -> Self {
        todo!()
    }

    fn next(&self) -> Self {
        todo!()
    }
}

impl HexSphere {
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

impl HexSphere {
    /// Iterates over the faces of this [`HexSphere`], starting with the given region.
    fn faces_from(&self, region: BaseRegion) -> FaceIter {
        if region.ty().is_edge() {
            FaceIter {
                sphere: *self,
                region,
                u: 3,
                u_end: self.kis.b(),
                v: 0,
            }
        } else if self.kis.c() < self.kis.b() {
            let n = self.kis.b() - self.kis.c();
            FaceIter {
                sphere: *self,
                region,
                u: 1,
                u_end: n - 1,
                v: 1,
            }
        } else {
            // This is a special case. When `b == c` there is exactly one interior vertex,
            // at (0, 0).
            FaceIter {
                sphere: *self,
                region,
                u: 0,
                u_end: 1,
                v: 0,
            }
        }
    }
}

/// An iterator over the faces of an [`TriSphere`].
#[derive(Clone, Debug)]
pub struct FaceIter {
    sphere: HexSphere,
    region: BaseRegion,
    u: u32,
    u_end: u32,
    v: u32,
}

impl Iterator for FaceIter {
    type Item = Face;
    fn next(&mut self) -> Option<Face> {
        loop {
            if self.u < self.u_end {
                let res = Face::new(tri::Vertex {
                    sphere: self.sphere.kis,
                    region: self.region,
                    u: self.u,
                    v: self.v,
                });
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
                    let res = Face::new(tri::Vertex {
                        sphere: self.sphere.kis,
                        region: self.region,
                        u: self.u,
                        v: self.v,
                    });
                    self.u += 3;
                    return Some(res);
                }
            } else {
                let n = self.sphere.kis.b() - self.sphere.kis.c();
                if self.v + 1 < n {
                    self.v += 1;
                    self.u = 1 + (self.v + 2) % 3;
                    self.u_end = n - self.v;
                    continue;
                }
            }
            if let Some(region) = self.sphere.kis.base().next_region(self.region) {
                *self = self.sphere.faces_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}

impl HexSphere {
    /// Iterates over the vertices of this [`HexSphere`], starting with the given region.
    fn vertices_from(&self, region: BaseRegion) -> VertexIter {
        if region.ty().is_edge() {
            VertexIter {
                sphere: *self,
                region,
                u: 1,
                u_end: self.kis.b(),
                v: 0,
                skip: false,
            }
        } else {
            let n = self.kis.b() - self.kis.c();
            VertexIter {
                sphere: *self,
                region,
                u: 2,
                u_end: n.saturating_sub(1),
                v: 1,
                skip: false,
            }
        }
    }
}

/// An iterator over the vertices of an [`HexSphere`].
#[derive(Clone, Debug)]
pub struct VertexIter {
    sphere: HexSphere,
    region: BaseRegion,
    u: u32,
    u_end: u32,
    v: u32,
    skip: bool,
}

impl Iterator for VertexIter {
    type Item = Vertex;
    fn next(&mut self) -> Option<Vertex> {
        loop {
            if self.u < self.u_end {
                let res = Vertex::new(tri::Vertex {
                    sphere: self.sphere.kis,
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
                    (self.u, self.skip) = match (self.v * 2) % 3 {
                        0 => (1, false),
                        1 => (1, true),
                        _ => (2, false),
                    };
                    continue;
                }
            } else {
                let n = self.sphere.kis.b() - self.sphere.kis.c();
                if self.v + 1 < n {
                    self.v += 1;
                    self.u_end = n - self.v;
                    (self.u, self.skip) = match (self.v * 2) % 3 {
                        0 => (1, false),
                        1 => (1, true),
                        _ => (2, false),
                    };
                    continue;
                }
            }
            if let Some(region) = self.sphere.kis.base().next_region(self.region) {
                *self = self.sphere.vertices_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}
