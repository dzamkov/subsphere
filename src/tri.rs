//! Contains types related to [`TriSphere`].
use crate::prelude::*;
use crate::basetri::BaseTriSphere;
use crate::basetri::Face as BaseFace;
use crate::basetri::HalfEdge as BaseHalfEdge;
use crate::basetri::Vertex as BaseVertex;
use crate::math::{mat, vec};
use crate::proj::{BaseTriProjector, Projection};
use std::num::NonZero;

/// A tessellation of the unit sphere into triangular [`Face`]s constructed by subdividing a
/// triangular platonic solid and projecting it onto the sphere.
///
/// This subdivision of a platonic solid is also known as a
/// [geodesic polyhedron](https://en.wikipedia.org/wiki/Geodesic_polyhedron).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct TriSphere<Proj> {
    proj: Proj,
    b: NonZero<u32>,
    c_base: u32,
}

impl<Proj> TriSphere<Proj> {
    /// Constructs a [`TriSphere`] from the given base shape, projection and subdivision
    /// parameters.
    ///
    /// The `b` and `c` parameters are as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation). For implementation
    /// simplicity and performance, only `b ≥ c` polyhedra are directly supported. `b < c`
    /// polyhedra can be emulated by swapping the `b` and `c` parameters and mirroring the
    /// resulting sphere.
    pub const fn new(base: BaseTriSphere, proj: Proj, b: NonZero<u32>, c: u32) -> Self {
        assert!(
            b.get() < u32::MAX >> 2,
            "b exceeds maxmimum subdivision size"
        );
        assert!(b.get() >= c, "b must be greater than or equal to c");
        Self {
            proj,
            b,
            c_base: (c << 2) | base as u32,
        }
    }

    /// The base shape of this subdivided sphere.
    pub const fn base(&self) -> BaseTriSphere {
        unsafe { std::mem::transmute((self.c_base & 0b11) as u8) }
    }

    /// The [`BaseTriProjector`] for this sphere.
    pub const fn proj(&self) -> &Proj {
        &self.proj
    }

    /// The `b` parameter of the subdivision, as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation).
    pub const fn b(&self) -> u32 {
        self.b.get()
    }

    /// The `c` parameter of the subdivision, as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation).
    pub const fn c(&self) -> u32 {
        self.c_base >> 2
    }

    /// The number of [`Face`]s each triangle on the base shape is subdivided into.
    ///
    /// This is also known as the "triangulation number" of the shape.
    pub const fn num_divisions(&self) -> usize {
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
    pub(crate) fn num_vertices_per_edge_region(&self) -> usize {
        let b = self.b() as usize;
        let c = self.c() as usize;
        (b - 1) * (c + 1)
    }

    /// The number of vertices per interior [`BaseRegion`].
    pub(crate) fn num_vertices_per_interior_region(&self) -> usize {
        let n = (self.b() - self.c()) as isize;
        // Yes, this really should return `1` when `n = 0`.
        ((n - 1) * (n - 2) / 2) as usize
    }
}

impl<Proj: Clone + BaseTriProjector> TriSphere<Proj> {
    /// Subdivides each edge of this [`TriSphere`] into the given number of segments.
    pub fn subdivide_edge(&self, divs: NonZero<u32>) -> Self {
        Self::new(
            self.base(),
            self.proj.clone(),
            self.b.checked_mul(divs).expect("arithmetic overflow"),
            self.c() * divs.get(),
        )
    }
}

impl<Proj: BaseTriProjector + Default> From<BaseTriSphere> for TriSphere<Proj> {
    fn from(base: BaseTriSphere) -> Self {
        Self::new(base, Proj::default(), NonZero::new(1).unwrap(), 0)
    }
}

impl<Proj: Eq + Clone + BaseTriProjector> crate::Sphere for TriSphere<Proj> {
    type Face = Face<Proj>;
    type Vertex = Vertex<Proj>;
    type HalfEdge = HalfEdge<Proj>;

    fn num_faces(&self) -> usize {
        self.num_divisions() * self.base().num_faces()
    }

    fn face(&self, index: usize) -> Face<Proj> {
        // TODO: "Real" implementation with constant-time performance
        self.faces().nth(index).expect("face index out of bounds")
    }

    #[allow(refining_impl_trait)]
    fn faces(&self) -> FaceIter<Proj> {
        self.clone().faces_from(self.base().first_region())
    }

    fn face_at(&self, point: [f64; 3]) -> Face<Proj> {
        let proj = self.sphere_proj();
        let (region, [u, v]) = proj.from_sphere_discrete(point);

        // The above only gives us an approximation of which face contains the point. This is
        // because the edges of a `TriSphere` are geodesics, but the boundaries of faces projected
        // according to `proj` usually aren't. We need to find a spherical triangle, whose vertices
        // are vertices of the `TriSphere`, that contains the point.
        let mut v_0 = proj.to_sphere(region, [(u + 1) as f64, v as f64]);
        let mut v_1 = proj.to_sphere(region, [u as f64, (v + 1) as f64]);
        let mut edge_0 = HalfEdge {
            sphere: self.clone(),
            region,
            start_u: (u + 1),
            start_v: v,
            dir: HalfEdgeDir::UnVp,
        };
        if mat::det_3([point, v_0, v_1]) < 0.0 {
            std::mem::swap(&mut v_0, &mut v_1);
            edge_0 = edge_0.complement();
        };

        // Iteratively check if the point is inside the triangle whose first side is `edge_0`.
        // If not, traverse across the edge the point is outside of. Prevent infinite loops by
        // limiting how many times each vertex can be involved in a check.
        let mut limit = 0;
        loop {
            let v_2 = edge_0.prev();
            let v_2 = proj.to_sphere(v_2.region, [v_2.start_u as f64, v_2.start_v as f64]);
            if limit <= 1 && mat::det_3([point, v_1, v_2]) < 0.0 {
                limit = limit.max(0);
                limit += 1;
                v_0 = v_2;
                edge_0 = edge_0.next().complement();
                continue;
            }
            if limit >= -1 && mat::det_3([point, v_2, v_0]) < 0.0 {
                limit = limit.min(0);
                limit -= 1;
                v_1 = v_2;
                edge_0 = edge_0.prev().complement();
                continue;
            }
            return edge_0.inside();
        }
    }

    fn num_vertices(&self) -> usize {
        self.num_faces() / 2 + 2
    }

    fn vertex(&self, index: usize) -> Vertex<Proj> {
        todo!()
    }

    #[allow(refining_impl_trait)]
    fn vertices(&self) -> VertexIter<Proj> {
        VertexIter {
            sphere: self.clone(),
            region: self.base().first_region(),
            u: 0,
            u_end: self.b(),
            v: 0,
        }
    }
}

/// Represents a face on a [`TriSphere`].
///
/// This is always a geodesic triangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Face<Proj> {
    pub(crate) sphere: TriSphere<Proj>,

    /// The [`BaseRegion`] which "owns" this face.
    pub(crate) region: BaseRegion,

    /// The U coordinate of the first vertex of this face, specified in the local coordinate space
    /// of `region`.
    pub(crate) u_0: u32,

    /// The U coordinate of the first vertex of this face, specified in the local coordinate space
    /// of `region`.
    pub(crate) v_0: u32,

    /// If `true`, the first [`HalfEdge`] boundary of this face will go in the +V direction.
    /// Otherwise, it will go in the +U direction.
    pub(crate) boundary_along_v: bool,
}

impl<Proj: Eq + Clone + BaseTriProjector> Face<Proj> {
    /// The approximate position of the center of this face.
    ///
    /// There is no strict definition or guarantees for this function, other than that the
    /// center of a face will be somewhere inside it. This is typically the fastest way to get a
    /// representative position on the face.
    pub fn center(&self) -> [f64; 3] {
        self.sphere.sphere_proj().to_sphere(
            self.region,
            if self.boundary_along_v {
                [self.u_0 as f64 - 1.0 / 3.0, self.v_0 as f64 + 2.0 / 3.0]
            } else {
                [self.u_0 as f64 + 1.0 / 3.0, self.v_0 as f64 + 1.0 / 3.0]
            },
        )
    }
}

#[test]
fn test_center() {
    use crate::util::tri_area;
    let sphere = TriSphere::new(
        BaseTriSphere::Icosa,
        crate::proj::Fuller,
        NonZero::new(4).unwrap(),
        1,
    );
    for face in sphere.faces() {
        let center = face.center();
        assert!(tri_area([face.vertex(0).pos(), face.vertex(1).pos(), center]) > 0.0);
        assert!(tri_area([face.vertex(1).pos(), face.vertex(2).pos(), center]) > 0.0);
        assert!(tri_area([face.vertex(2).pos(), face.vertex(0).pos(), center]) > 0.0);
    }
}

impl<Proj: Eq + Clone + BaseTriProjector> crate::Face for Face<Proj> {
    type Vertex = Vertex<Proj>;
    type HalfEdge = HalfEdge<Proj>;

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
        let proj = self.sphere.sphere_proj();
        let p_0 = proj.to_sphere(self.region, [u_0 as f64, v_0 as f64]);
        let p_1 = proj.to_sphere(self.region, [u_1 as f64, v_1 as f64]);
        let p_2 = proj.to_sphere(self.region, [u_2 as f64, v_2 as f64]);
        crate::util::tri_area([p_0, p_1, p_2])
    }

    fn num_sides(&self) -> usize {
        3
    }

    fn side(&self, index: usize) -> HalfEdge<Proj> {
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
            sphere: self.sphere.clone(),
            region: self.region,
            start_u: self.u_0.wrapping_add_signed(d_u),
            start_v: self.v_0.wrapping_add_signed(d_v),
            dir,
        }
    }
}

impl<Proj> Face<Proj> {
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

/// Represents a vertex in an [`TriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Vertex<Proj> {
    pub(crate) sphere: TriSphere<Proj>,

    /// The [`BaseRegion`] which "owns" this vertex.
    pub(crate) region: BaseRegion,

    /// The U coordinate of this vertex, specified in the local coordinate space of `region`.
    pub(crate) u: u32,

    /// The V coordinate of this vertex, specified in the local coordinate space of `region`.
    pub(crate) v: u32,
}

impl<Proj: Eq + Clone + BaseTriProjector> Vertex<Proj> {
    /// Constructs a [`Vertex`] with the given properties.
    ///
    /// This will normalize the vertex `region` to be its proper owner.
    pub(crate) fn new(sphere: TriSphere<Proj>, region: BaseRegion, u: u32, v: u32) -> Self {
        if region.ty().is_edge() {
            if v == 0 {
                Self::on_edge_u_boundary(sphere, region.as_edge(), u)
            } else if u == 0 {
                Self::on_edge_v_boundary(sphere, region.as_edge(), v)
            } else if v >= sphere.c() {
                debug_assert_eq!(v, sphere.c());
                let u = sphere.b() - u;
                Self::on_edge_u_boundary(sphere, region.as_edge().complement(), u)
            } else if u >= sphere.b() {
                debug_assert_eq!(u, sphere.b());
                let v = sphere.c() - v;
                Self::on_edge_v_boundary(sphere, region.as_edge().complement(), v)
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
    fn on_interior_boundary(sphere: TriSphere<Proj>, bottom: BaseHalfEdge, u: u32) -> Self {
        let u = sphere.b() - u;
        Self::on_edge_u_boundary(sphere, bottom.complement(), u)
    }

    /// Constructs a [`Vertex`] on the V (U = 0) boundary of an edge region.
    fn on_edge_v_boundary(sphere: TriSphere<Proj>, edge: BaseHalfEdge, v: u32) -> Self {
        Self::on_edge_u_boundary(sphere, edge.prev().complement(), v)
    }

    /// Constructs a [`Vertex`] on the U (V = 0) boundary of an edge region.
    fn on_edge_u_boundary(sphere: TriSphere<Proj>, edge: BaseHalfEdge, u: u32) -> Self {
        if u == sphere.b() {
            // In the special case where `b == c`, it's the interior region which owns the
            // central vertex.
            if sphere.b() == sphere.c() {
                Self::center(sphere, edge.complement().inside())
            } else {
                let u = sphere.c();
                Self::on_edge_u_boundary_exclusive(sphere, edge.complement().prev().complement(), u)
            }
        } else {
            Self::on_edge_u_boundary_exclusive(sphere, edge, u)
        }
    }

    /// Constructs a [`Vertex`] on the U (V = 0) boundary of an edge region when it is guaranteed
    /// that `u` is not `b`.
    fn on_edge_u_boundary_exclusive(sphere: TriSphere<Proj>, edge: BaseHalfEdge, u: u32) -> Self {
        if u == 0 {
            Self::base(sphere, edge.start())
        } else {
            let comp = match edge.side_index() {
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
                        let u = sphere.b() - u;
                        let v = sphere.c();
                        return Self {
                            sphere,
                            region: BaseRegion::new(edge.inside(), BaseRegionType::Edge2),
                            u,
                            v,
                        };
                    } else {
                        edge.complement()
                    }
                }
            };
            if comp.side_index() == 0 {
                let u = sphere.b() - u;
                let v = sphere.c();
                Self {
                    sphere,
                    region: BaseRegion::new(comp.inside(), BaseRegionType::Edge0),
                    u,
                    v,
                }
            } else {
                debug_assert_eq!(comp.side_index(), 2);
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
    fn center(sphere: TriSphere<Proj>, face: BaseFace) -> Self {
        debug_assert_eq!(sphere.b(), sphere.c());
        Self {
            sphere,
            region: BaseRegion::new(face, BaseRegionType::Interior),
            u: 0,
            v: 0,
        }
    }

    /// Constructs a [`Vertex`] for one of the base vertices.
    pub fn base(sphere: TriSphere<Proj>, source: BaseVertex) -> Self {
        if source.index() == 0 {
            Self {
                sphere,
                region: BaseRegion::new(source.owner(), BaseRegionType::Edge0),
                u: 0,
                v: 0,
            }
        } else {
            let u = sphere.b();
            let v = sphere.c();
            Self {
                sphere,
                region: BaseRegion::new(source.owner(), BaseRegionType::Edge0),
                u,
                v,
            }
        }
    }

    /// Determines whether this vertex corresponds to a [`BaseVertex`] and if so, returns it.
    pub fn as_base(&self) -> Option<BaseVertex> {
        match self.region.ty() {
            BaseRegionType::Edge0 => {
                if self.u == 0 {
                    debug_assert_eq!(self.v, 0);
                    Some(self.region.owner().vertex(0))
                } else if self.u == self.sphere.b() {
                    debug_assert_eq!(self.v, self.sphere.c());
                    Some(self.region.owner().vertex(1))
                } else {
                    None
                }
            }
            BaseRegionType::Edge2 => {
                if self.u == self.sphere.b() {
                    debug_assert_eq!(self.v, self.sphere.c());
                    Some(self.region.owner().vertex(2))
                } else {
                    None
                }
            }
            _ => None,
        }
    }
}

impl<Proj: Eq + Clone + BaseTriProjector> crate::Vertex for Vertex<Proj> {
    type Face = Face<Proj>;
    type HalfEdge = HalfEdge<Proj>;

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
            .sphere_proj()
            .to_sphere(self.region, [self.u as f64, self.v as f64])
    }

    fn degree(&self) -> usize {
        if let Some(base) = self.as_base() {
            base.degree()
        } else {
            6
        }
    }

    fn outgoing(&self, index: usize) -> HalfEdge<Proj> {
        if let Some(base) = self.as_base() {
            HalfEdge::base(self.sphere.clone(), base.outgoing(index))
        } else {
            HalfEdge::new(
                self.sphere.clone(),
                self.region,
                self.u,
                self.v,
                HalfEdgeDir::from_index(index),
            )
        }
    }
}

/// Represents one "side" of an edge in an [`TriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfEdge<Proj> {
    sphere: TriSphere<Proj>,

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
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum HalfEdgeDir {
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

impl<Proj> HalfEdge<Proj> {
    /// Attempts to construct a [`HalfEdge`] with the given properties, normalizing
    /// `region` to be its proper owner.
    ///
    /// This requires `(start_u, start_v)` to be a valid coordinate in the given `region`.
    pub(crate) fn new(
        sphere: TriSphere<Proj>,
        region: BaseRegion,
        start_u: u32,
        start_v: u32,
        dir: HalfEdgeDir,
    ) -> Self {
        if region.ty().is_edge() {
            if start_v == 0 {
                Self::on_edge_u_boundary(sphere, region.as_edge(), start_u, dir)
            } else if start_u == 0 {
                Self::on_edge_v_boundary(sphere, region.as_edge(), start_v, dir)
            } else if start_v >= sphere.c() {
                debug_assert_eq!(start_v, sphere.c());
                let start_u = sphere.b() - start_u;
                Self::on_edge_u_boundary(
                    sphere,
                    region.as_edge().complement(),
                    start_u,
                    dir.rotate_ccw(3),
                )
            } else if start_u >= sphere.b() {
                debug_assert_eq!(start_u, sphere.b());
                let start_v = sphere.c() - start_v;
                Self::on_edge_v_boundary(
                    sphere,
                    region.as_edge().complement(),
                    start_v,
                    dir.rotate_ccw(3),
                )
            } else {
                Self {
                    sphere,
                    region,
                    start_u,
                    start_v,
                    dir,
                }
            }
        } else if start_v == 0 {
            Self::on_interior_boundary(sphere, region.owner().side(0), start_u, dir)
        } else {
            let n = sphere.b() - sphere.c();
            if start_u == 0 {
                Self::on_interior_boundary(
                    sphere,
                    region.owner().side(2),
                    n - start_v,
                    dir.rotate_ccw(2),
                )
            } else if start_u + start_v >= n {
                debug_assert_eq!(start_u + start_v, n);
                Self::on_interior_boundary(
                    sphere,
                    region.owner().side(1),
                    start_v,
                    dir.rotate_ccw(4),
                )
            } else {
                Self {
                    sphere,
                    region,
                    start_u,
                    start_v,
                    dir,
                }
            }
        }
    }

    /// Constructs a [`HalfEdge`] on the bottom (V = 0) boundary of an interior region.
    fn on_interior_boundary(
        sphere: TriSphere<Proj>,
        bottom: BaseHalfEdge,
        start_u: u32,
        dir: HalfEdgeDir,
    ) -> Self {
        if start_u == 0 {
            let c = sphere.c();
            if dir > HalfEdgeDir::Up {
                return Self::on_edge_u_boundary(
                    sphere,
                    bottom.prev().complement(),
                    c + start_u,
                    dir.rotate_ccw(5),
                );
            }
        }
        let n = sphere.b() - sphere.c();
        let in_region = if start_u >= n {
            debug_assert_eq!(start_u, n);
            dir == HalfEdgeDir::UnVp
        } else {
            dir < HalfEdgeDir::Un
        };
        if !in_region {
            let start_u = sphere.b() - start_u;
            return Self::on_edge_u_boundary(
                sphere,
                bottom.complement(),
                start_u,
                dir.rotate_ccw(3),
            );
        }
        let (start_u, start_v, dir) = match bottom.side_index() {
            0 => (start_u, 0, dir),
            1 => (n - start_u, start_u, dir.rotate_ccw(2)),
            _ => (0, n - start_u, dir.rotate_ccw(4)),
        };
        Self {
            sphere,
            region: BaseRegion::new(bottom.inside(), BaseRegionType::Interior),
            start_u,
            start_v,
            dir,
        }
    }

    /// Constructs a [`HalfEdge`] on the V (U = 0) boundary of an edge region.
    ///
    /// Assumes that `start_v > 0`.
    fn on_edge_v_boundary(
        sphere: TriSphere<Proj>,
        edge: BaseHalfEdge,
        start_v: u32,
        dir: HalfEdgeDir,
    ) -> Self {
        debug_assert!(start_v > 0);
        if dir >= HalfEdgeDir::Vn {
            Self::on_edge_exclusive(sphere, edge, 0, start_v, dir)
        } else if start_v >= sphere.c() {
            debug_assert_eq!(start_v, sphere.c());
            Self::on_interior_boundary(sphere, edge, 0, dir)
        } else if dir == HalfEdgeDir::Up {
            Self::on_edge_exclusive(sphere, edge, 0, start_v, dir)
        } else {
            Self::on_edge_u_boundary(sphere, edge.prev().complement(), start_v, dir.rotate_ccw(5))
        }
    }

    /// Constructs a [`HalfEdge`] on the U (V = 0) boundary of an edge region.
    fn on_edge_u_boundary(
        sphere: TriSphere<Proj>,
        edge: BaseHalfEdge,
        start_u: u32,
        dir: HalfEdgeDir,
    ) -> Self {
        if sphere.c() == 0 {
            Self::on_interior_boundary(sphere, edge, start_u, dir)
        } else if start_u == 0 {
            let mut edge = edge;
            let mut dir = dir;
            while dir > HalfEdgeDir::Up {
                edge = edge.prev().complement();
                dir = dir.rotate_ccw(5);
            }
            Self::on_edge_exclusive(sphere, edge, 0, 0, HalfEdgeDir::Up)
        } else if start_u == sphere.b() && dir == HalfEdgeDir::Up {
            let c = sphere.c();
            Self::on_edge_u_boundary(
                sphere,
                edge.complement().prev().complement(),
                c,
                HalfEdgeDir::UnVp,
            )
        } else if dir < HalfEdgeDir::Un {
            Self::on_edge_exclusive(sphere, edge, start_u, 0, dir)
        } else if start_u <= sphere.c() {
            Self::on_edge_v_boundary(sphere, edge.complement().next(), start_u, dir.rotate_ccw(1))
        } else {
            let start_u = sphere.b() - start_u;
            Self::on_interior_boundary(sphere, edge.complement(), start_u, dir.rotate_ccw(3))
        }
    }

    /// Constructs a [`HalfEdge`] on an edge region when it is guaranteed that the resulting
    /// half-edge will belong to the region.
    fn on_edge_exclusive(
        sphere: TriSphere<Proj>,
        edge: BaseHalfEdge,
        start_u: u32,
        start_v: u32,
        dir: HalfEdgeDir,
    ) -> Self {
        debug_assert!(sphere.c() > 0);
        let comp = match edge.side_index() {
            0 => {
                return Self {
                    sphere,
                    region: BaseRegion::new(edge.inside(), BaseRegionType::Edge0),
                    start_u,
                    start_v,
                    dir,
                };
            }
            1 => edge.complement(),
            _ => {
                if edge.inside().owns_edge_2() {
                    let start_u = sphere.b() - start_u;
                    let start_v = sphere.c() - start_v;
                    return Self {
                        sphere,
                        region: BaseRegion::new(edge.inside(), BaseRegionType::Edge2),
                        start_u,
                        start_v,
                        dir: dir.rotate_ccw(3),
                    };
                } else {
                    edge.complement()
                }
            }
        };
        if comp.side_index() == 0 {
            let start_u = sphere.b() - start_u;
            let start_v = sphere.c() - start_v;
            Self {
                sphere,
                region: BaseRegion::new(comp.inside(), BaseRegionType::Edge0),
                start_u,
                start_v,
                dir: dir.rotate_ccw(3),
            }
        } else {
            debug_assert_eq!(comp.side_index(), 2);
            debug_assert!(comp.inside().owns_edge_2());
            Self {
                sphere,
                region: BaseRegion::new(comp.inside(), BaseRegionType::Edge2),
                start_u,
                start_v,
                dir,
            }
        }
    }

    /// Constructs a [`HalfEdge`] corresponding to a [`BaseHalfEdge`].
    fn base(sphere: TriSphere<Proj>, source: BaseHalfEdge) -> Self {
        if sphere.c() == 0 {
            // In this case, edge regions do not own any faces, and therefore no half-edges.
            let (start_u, start_v, dir) = match source.side_index() {
                0 => (0, 0, HalfEdgeDir::Up),
                1 => (sphere.b(), 0, HalfEdgeDir::UnVp),
                _ => (0, sphere.b(), HalfEdgeDir::Vn),
            };
            Self {
                sphere,
                region: BaseRegion::new(source.inside(), BaseRegionType::Interior),
                start_u,
                start_v,
                dir,
            }
        } else {
            Self::on_edge_exclusive(sphere, source, 0, 0, HalfEdgeDir::Up)
        }
    }
}

impl<Proj: Eq + Clone + BaseTriProjector> crate::HalfEdge for HalfEdge<Proj> {
    type Vertex = Vertex<Proj>;
    type Face = Face<Proj>;

    fn side_index(&self) -> usize {
        self.dir as usize / 2
    }

    fn length(&self) -> f64 {
        crate::util::dist(self.start().pos(), self.next().start().pos())
    }

    fn angle(&self) -> f64 {
        let v_a = self.prev().start().pos();
        let v_b = self.start().pos();
        let v_c = self.next().start().pos();
        crate::util::angle(v_a, v_b, v_c)
    }

    fn inside(&self) -> Face<Proj> {
        let (d_u, d_v) = match self.dir {
            HalfEdgeDir::Up => (0, 0),
            HalfEdgeDir::Vp => (0, 0),
            HalfEdgeDir::UnVp => (-1, 0),
            HalfEdgeDir::Un => (0, -1),
            HalfEdgeDir::Vn => (0, -1),
            HalfEdgeDir::UpVn => (1, -1),
        };
        Face {
            sphere: self.sphere.clone(),
            region: self.region,
            u_0: self.start_u.wrapping_add_signed(d_u),
            v_0: self.start_v.wrapping_add_signed(d_v),
            boundary_along_v: (self.dir as usize) % 2 == 1,
        }
    }

    fn start(&self) -> Vertex<Proj> {
        Vertex::new(self.sphere.clone(), self.region, self.start_u, self.start_v)
    }

    fn complement(&self) -> HalfEdge<Proj> {
        let [d_u, d_v] = self.dir.into_vec();
        HalfEdge::new(
            self.sphere.clone(),
            self.region,
            self.start_u.wrapping_add_signed(d_u),
            self.start_v.wrapping_add_signed(d_v),
            self.dir.rotate_ccw(3),
        )
    }

    fn prev(&self) -> HalfEdge<Proj> {
        let dir = self.dir.rotate_ccw(4);
        let [d_u, d_v] = dir.into_vec();
        HalfEdge {
            sphere: self.sphere.clone(),
            region: self.region,
            start_u: self.start_u.wrapping_add_signed(-d_u),
            start_v: self.start_v.wrapping_add_signed(-d_v),
            dir,
        }
    }

    fn next(&self) -> HalfEdge<Proj> {
        let [d_u, d_v] = self.dir.into_vec();
        HalfEdge {
            sphere: self.sphere.clone(),
            region: self.region,
            start_u: self.start_u.wrapping_add_signed(d_u),
            start_v: self.start_v.wrapping_add_signed(d_v),
            dir: self.dir.rotate_ccw(2),
        }
    }
}

impl HalfEdgeDir {
    /// Constructs a [`HalfEdgeDir`] from the given index.
    ///
    /// This corresponds to an angle of `index * 60` degrees.
    pub fn from_index(index: usize) -> Self {
        assert!(index < 6, "index out of bounds: {}", index);
        unsafe { std::mem::transmute(index as u8) }
    }

    /// Converts this direction into a vector representation.
    pub fn into_vec(self) -> [i32; 2] {
        let u = (0x020100000102u64 >> (self as usize * 8)) as u8 as i32 - 1;
        let v = (0x000001020201u64 >> (self as usize * 8)) as u8 as i32 - 1;
        [u, v]
    }

    /// Rotates this [`HalfEdgeDir`] counter-clockwise by `amount * 60` degrees.
    pub fn rotate_ccw(self, amount: usize) -> Self {
        let index = (self as usize + amount) % 6;
        unsafe { std::mem::transmute(index as u8) }
    }
}

impl<Proj> TriSphere<Proj> {
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

impl<Proj> TriSphere<Proj> {
    /// Iterates over the faces of this [`TriSphere`], starting with the given region.
    fn faces_from(self, region: BaseRegion) -> FaceIter<Proj> {
        if region.ty().is_edge() {
            // If `c` is zero, there are no faces in an edge region. This iterator should
            // immediately go to the next region.
            let u_0_end = self.b() * (self.c() > 0) as u32;
            let boundary_along_v = self.c() == 0;
            FaceIter {
                sphere: self,
                region,
                u_0: 0,
                u_0_end,
                v_0: 0,
                boundary_along_v,
            }
        } else {
            let n = self.b() - self.c();
            FaceIter {
                sphere: self,
                region,
                u_0: 0,
                u_0_end: n,
                v_0: 0,
                boundary_along_v: false,
            }
        }
    }
}

/// An iterator over the faces of an [`TriSphere`].
#[derive(Clone, Debug)]
pub struct FaceIter<Proj> {
    sphere: TriSphere<Proj>,
    region: BaseRegion,
    u_0: u32,
    u_0_end: u32,
    v_0: u32,
    boundary_along_v: bool,
}

impl<Proj: Eq + Clone + BaseTriProjector> Iterator for FaceIter<Proj> {
    type Item = Face<Proj>;
    fn next(&mut self) -> Option<Face<Proj>> {
        loop {
            if self.u_0 < self.u_0_end {
                let res = Face {
                    sphere: self.sphere.clone(),
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
                *self = self.sphere.clone().faces_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}

impl<Proj> TriSphere<Proj> {
    /// Iterates over the vertices of this [`TriSphere`], starting with the given region.
    fn vertices_from(self, region: BaseRegion) -> VertexIter<Proj> {
        if region.ty().is_edge() {
            let u_end = self.b();
            VertexIter {
                sphere: self,
                region,
                u: 1,
                u_end,
                v: 0,
            }
        } else if self.c() < self.b() {
            let n = self.b() - self.c();
            VertexIter {
                sphere: self,
                region,
                u: 1,
                u_end: n - 1,
                v: 1,
            }
        } else {
            // This is a special case. When `b == c` there is exactly one interior vertex,
            // at (0, 0).
            VertexIter {
                sphere: self,
                region,
                u: 0,
                u_end: 1,
                v: 0,
            }
        }
    }
}

/// An iterator over the vertices of an [`TriSphere`].
#[derive(Clone, Debug)]
pub struct VertexIter<Proj> {
    sphere: TriSphere<Proj>,
    region: BaseRegion,
    u: u32,
    u_end: u32,
    v: u32,
}

impl<Proj: Eq + Clone + BaseTriProjector> Iterator for VertexIter<Proj> {
    type Item = Vertex<Proj>;
    fn next(&mut self) -> Option<Vertex<Proj>> {
        loop {
            if self.u < self.u_end {
                let res = Vertex {
                    sphere: self.sphere.clone(),
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
                        sphere: self.sphere.clone(),
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
                *self = self.sphere.clone().vertices_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}

/// Identifies a region of a [`BaseTriSphere`] which can contain faces and vertices of a
/// [`TriSphere`].
///
/// There is a region corresponding to each face and edge of a [`BaseTriSphere`].
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct BaseRegion(u8);

impl BaseRegion {
    /// Constructs the [`BaseRegion`] for the given face and type.
    pub fn new(owner: BaseFace, ty: BaseRegionType) -> Self {
        BaseRegion((owner.0 << 2) | (ty as u8))
    }

    /// Gets the owner of this region.
    pub fn owner(&self) -> BaseFace {
        BaseFace(self.0 >> 2)
    }

    /// The general type of this [`BaseRegion`].
    pub fn ty(&self) -> BaseRegionType {
        unsafe { std::mem::transmute(self.0 & 0b11) }
    }

    /// Assuming that this region is for an edge, gets the edge corresponding to this region. The
    /// start of the returned edge is the origin of this region.
    pub fn as_edge(self) -> BaseHalfEdge {
        if self.ty() == BaseRegionType::Edge0 {
            self.owner().side(0)
        } else {
            self.owner().side(2).complement()
        }
    }
}

impl std::fmt::Debug for BaseRegion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("BaseRegion::new")
            .field(&self.owner())
            .field(&self.ty())
            .finish()
    }
}

/// The general type of a [`BaseRegion`].
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum BaseRegionType {
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
    pub(crate) fn first_region(self) -> BaseRegion {
        BaseRegion::new(self.face(0), BaseRegionType::Edge0)
    }

    /// Gets the [`BaseRegion`] after the given one.
    pub(crate) fn next_region(self, mut region: BaseRegion) -> Option<BaseRegion> {
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

/// A piecewise projection which maps local coordinates on a [`BaseRegion`] to points on the
/// unit sphere, according to the parameters of a [`TriSphere`].
struct SphereProjection<'a, Proj> {
    /// The [`TriSphere`] this projection is based on.
    sphere: &'a TriSphere<Proj>,

    /// The local coordinates of the origin of an interior region on its corresponding base
    /// face.
    p_0: [f64; 2],

    /// The linear transformation which maps local coordinates on a region to coordinates on
    /// its corresponding base face.
    linear: [[f64; 2]; 2],

    /// The reciprocal of the determinant of `linear`.
    inv_linear_det: f64,
}

impl<Proj> TriSphere<Proj> {
    /// Constructs a [`SphereProjection`] based on the parameters of this [`TriSphere`].
    fn sphere_proj(&self) -> SphereProjection<Proj> {
        let b = self.b() as f64;
        let c = self.c() as f64;
        let w_total = b * b + b * c + c * c;
        let w_1 = c * c / w_total;
        let w_2 = b * c / w_total;
        let p_0 = [w_1, w_2];
        let u = vec::div(vec::sub([1.0, 0.0], p_0), b);
        let v = if c > 0.0 {
            vec::div(p_0, c)
        } else {
            [0.0, 1.0 / b]
        };
        let linear = [u, v];
        SphereProjection {
            sphere: self,
            p_0,
            linear,
            inv_linear_det: 1.0 / mat::det_2(linear),
        }
    }
}

impl<Proj: BaseTriProjector> SphereProjection<'_, Proj> {
    /// Projects a point in the local coordinates of a given region to a point on the unit sphere.
    pub fn to_sphere(&self, region: BaseRegion, coords: [f64; 2]) -> [f64; 3] {
        if region.ty().is_edge() {
            let coords = mat::apply(self.linear, coords);
            if coords[1] >= 0.0 {
                self.sphere
                    .proj()
                    .inside(region.as_edge())
                    .to_sphere(coords)
            } else {
                self.sphere
                    .proj()
                    .inside(region.as_edge().complement())
                    .to_sphere(vec::add([1.0, 0.0], vec::neg(coords)))
            }
        } else {
            self.sphere
                .proj()
                .inside(region.owner().side(0))
                .to_sphere(vec::add(self.p_0, mat::apply(self.linear, coords)))
        }
    }

    /// Projects a point on the unit sphere to "discrete" local coordinates of some region.
    #[expect(clippy::wrong_self_convention)]
    pub fn from_sphere_discrete(&self, point: [f64; 3]) -> (BaseRegion, [u32; 2]) {
        let face = self.sphere.base().face_at(point);
        let coords = self.sphere.proj().inside(face.side(0)).from_sphere(point);
        let b = self.sphere.b();
        let c = self.sphere.c();
        if c == 0 {
            let [u, v] = vec::mul(coords, b as f64);
            let u = (u.max(0.0) as u32).min(b - 1);
            let v = (v.max(0.0) as u32).min(b - u - 1);
            debug_assert!(u + v < b);
            return (BaseRegion::new(face, BaseRegionType::Interior), [u, v]);
        }
        let coords = vec::sub(coords, self.p_0);
        let coords = vec::mul(coords, self.inv_linear_det);
        let [u, v] = mat::apply(mat::adjoint_2(self.linear), coords);
        let n = b - c;
        let (comp, [r_u, r_v]) = if u + v < n as f64 {
            if u >= 0.0 {
                let u = u as u32;
                if v >= 0.0 {
                    let v = v as u32;
                    debug_assert!(u + v < n);
                    return (BaseRegion::new(face, BaseRegionType::Interior), [u, v]);
                } else {
                    debug_assert!(u < b);
                    let r_v = ((v + c as f64).max(0.0) as u32).min(c - 1);
                    return (BaseRegion::new(face, BaseRegionType::Edge0), [u, r_v]);
                }
            } else {
                let r_u = ((u + (v + c as f64)).max(0.0) as u32).min(b - 1);
                let r_v = ((-u) as u32).min(c - 1);
                if face.owns_edge_2() {
                    return (BaseRegion::new(face, BaseRegionType::Edge2), [r_u, r_v]);
                } else {
                    (face.side(2).complement(), [r_u, r_v])
                }
            }
        } else if v >= 0.0 {
            let r_u = ((b as f64 - v).max(0.0) as u32).min(b - 1);
            let r_v = ((u - (n as f64 - v)) as u32).min(c - 1);
            (face.side(1).complement(), [r_u, r_v])
        } else {
            let r_u = (u as u32).min(b - 1);
            let r_v = ((v + c as f64).max(0.0) as u32).min(c - 1);
            return (BaseRegion::new(face, BaseRegionType::Edge0), [r_u, r_v]);
        };
        if comp.side_index() == 0 {
            (BaseRegion::new(comp.inside(), BaseRegionType::Edge0), [r_u, r_v])
        } else {
            debug_assert_eq!(comp.side_index(), 2);
            debug_assert!(comp.inside().owns_edge_2());
            (BaseRegion::new(comp.inside(), BaseRegionType::Edge2), [b - r_u - 1, c - r_v - 1])
        }
    }
}
