//! Contains types related to [`SubTriSphere`].
use crate::Face as BaseFace;
use crate::HalfEdge as BaseHalfEdge;
use crate::Vertex as BaseVertex;
use crate::math::{Vector3, mat3, vec3};
use std::num::NonZeroU32;

/// A tessellation of the unit sphere into triangular [`Face`]s constructed by subdividing a
/// base polyhedron and projecting it onto the sphere.
///
/// This can be used to construct a sphere tessellation based on any
/// [geodesic polyhedron](https://en.wikipedia.org/wiki/Geodesic_polyhedron).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SubTriSphere<Base: BaseSphere> {
    base: Base,
    b: NonZeroU32,
    c: u32,
}

/// Describes a potential base shape for a [`SubTriSphere`].
/// 
/// This is a [`crate::Sphere`] with equilateral triangular faces. There are effectively
/// three topologically distinct implementations of this trait:
///  * `BaseTetraSphere` (**TODO**)
///  * `BaseOctoSphere` (**TODO**)
///  * [`BaseIcoSphere`](crate::ico::BaseIcoSphere)
#[allow(private_bounds)]
pub trait BaseSphere: crate::Sphere + Eq + Clone + BaseSphereInternal {

}

/// Internal trait for [`BaseSphere`].
/// 
/// This is a [`crate::Sphere`] whose faces are triangular, which also assigns an "owner" (face) to
/// each vertex and edge, subject to the following rules:
///  * Every vertex and edge is owned by exactly one face.
///  * Faces must own their first edge.
///  * Faces must not own their second edge.
///  * Faces may own their third edge.
///  * Faces may only own their first vertex.
pub(crate) trait BaseSphereInternal: crate::Sphere + Eq + Clone {
    /// The type of [`BaseRegion`] for this base shape.
    type Region: BaseRegion<Face = Self::Face>;

    /// Indicates whether the given face owns its first vertex.
    fn face_owns_vertex_0(&self, face: Self::Face) -> bool;

    /// Indicates whether the given face owns its third edge.
    fn face_owns_edge_2(&self, face: Self::Face) -> bool;

    /// Gets the face which owns the given vertex.
    fn vertex_owner(&self, vertex: Self::Vertex) -> Self::Face;

    /// Gets the number of vertices that are owned by the faces preceding the given face in
    /// iteration order.
    fn num_owned_vertices_before(&self, face: Self::Face) -> usize;

    /// Gets the number of edges that are owned by the faces preceding the given face in iteration
    /// order.
    fn num_owned_edges_before(&self, face: Self::Face) -> usize;

    /// Gets the first [`BaseRegion`] of this base shape.
    fn first_region(&self) -> Self::Region;

    /// Gets the [`BaseRegion`] after the given one.
    fn next_region(&self, region: Self::Region) -> Option<Self::Region>;

    // TODO: Remove
    #[allow(missing_docs)]
    fn interpolate(&self, weights: [f64; 3]) -> [f64; 3];
}

/// Identifies a region of a [`BaseSphere`] which can contain faces and vertices of a
/// [`SubTriSphere`].
///
/// There is a region corresponding to each face and edge of a [`BaseSphere`].
pub trait BaseRegion: Clone + Eq {
    /// The type of [`crate::Face`] for the [`BaseSphere`].
    type Face: crate::Face;

    /// Constructs the [`BaseRegion`] for the given face and type.
    fn new(face: Self::Face, ty: BaseRegionType) -> Self;

    /// The owner of this region.
    fn owner(&self) -> Self::Face;

    /// The general type of this [`BaseRegion`].
    fn ty(&self) -> BaseRegionType;
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

impl<Base: BaseSphere + Default> SubTriSphere<Base> {
    /// Constructs a [`SubTriSphere`] wrapper over the base shape, without any subdivision.
    pub fn base() -> Self {
        Self::new_internal(Base::default(), NonZeroU32::new(1).unwrap(), 0)
    }

    /// Constructs a [`SubTriSphere`] with the given subdivision parameters.
    ///
    /// The `b` and `c` parameters are as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation). For implementation
    /// simplicity and performance, only `b ≥ c` polyhedra are directly supported. `b < c`
    /// polyhedra can be emulated by simply swapping the `b` and `c` parameters and mirroring the
    /// resulting sphere.
    pub fn new(b: NonZeroU32, c: u32) -> Self {
        Self::new_internal(Base::default(), b, c)
    }
}

impl<Base: BaseSphere> SubTriSphere<Base> {
    /// Constructs a [`SubTriSphere`] by subdividing the given base shape.
    pub(crate) fn new_internal(base: Base, b: NonZeroU32, c: u32) -> Self {
        assert!(
            b.get() >= c,
            "b must be greater than or equal to c. b = {}, c = {}",
            b,
            c
        );
        Self { base, b, c }
    }

    /// The `b` parameter of the subdivision, as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation).
    pub fn b(&self) -> u32 {
        self.b.get()
    }

    /// The `c` parameter of the subdivision, as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation).
    pub fn c(&self) -> u32 {
        self.c
    }

    /// The number of [`Face`]s each triangle on the base shape is subdivided into.
    ///
    /// This is also known as the "triangulation number" of the shape.
    pub fn num_divisions(&self) -> usize {
        let b = self.b() as usize;
        let c = self.c as usize;
        b * (b + c) + c * c
    }

    /// The number of faces per edge [`BaseRegion`].
    fn num_faces_per_edge_region(&self) -> usize {
        let b = self.b() as usize;
        let c = self.c as usize;
        b * c * 2
    }

    /// The number of faces per interior [`BaseRegion`].
    fn num_faces_per_interior_region(&self) -> usize {
        let n = (self.b() - self.c) as usize;
        n * n
    }

    /// The number of vertices per edge [`BaseRegion`], assuming that the origin vertex is not
    /// counted.
    fn num_vertices_per_edge_region(&self) -> usize {
        let b = self.b() as usize;
        let c = self.c as usize;
        (b - 1) * (c + 1)
    }

    /// The number of vertices per interior [`BaseRegion`].
    fn num_vertices_per_interior_region(&self) -> usize {
        let n = (self.b() - self.c) as isize;
        // Yes, this really should return `1` when `n = 0`.
        ((n - 1) * (n - 2) / 2) as usize
    }
}

impl<Base: BaseSphere> crate::Sphere for SubTriSphere<Base> {
    type Face = Face<Base>;
    type Vertex = Vertex<Base>;
    type HalfEdge = HalfEdge<Base>;

    fn num_faces(&self) -> usize {
        self.num_divisions() * self.base.num_faces()
    }

    #[allow(refining_impl_trait)]
    fn faces(&self) -> FaceIter<Base> {
        self.faces_from(self.base.first_region())
    }

    fn num_vertices(&self) -> usize {
        self.num_faces() / 2 + 2
    }

    #[allow(refining_impl_trait)]
    fn vertices(&self) -> VertexIter<Base> {
        self.vertices_from(self.base.first_region())
    }
}

/// Represents a face in an [`SubTriSphere`].
///
/// This is always a geodesic triangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Face<Base: BaseSphere> {
    sphere: SubTriSphere<Base>,

    /// The [`BaseRegion`] which "owns" this face.
    region: Base::Region,

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

impl<Base: BaseSphere> crate::Face for Face<Base> {
    type Vertex = Vertex<Base>;
    type HalfEdge = HalfEdge<Base>;

    fn index(&self) -> usize {
        if self.region.ty().is_edge() {
            self.sphere.base_face_index(self.region.clone())
                + (2 * self.v_0 as usize + self.boundary_along_v as usize)
                    * self.sphere.b() as usize
                + self.u_0 as usize
                - self.boundary_along_v as usize
        } else {
            let n = (self.sphere.b() - self.sphere.c) as usize;
            self.sphere.base_face_index(self.region.clone())
                + self.v_0 as usize * n * 2
                + self.boundary_along_v as usize * (n - self.v_0 as usize - 1)
                - self.v_0 as usize * self.v_0 as usize
                + self.u_0 as usize
        }
    }

    fn area(&self) -> f64 {
        let [[u_0, v_0], [u_1, v_1], [u_2, v_2]] = self.local_coords();
        let p_0 = self
            .sphere
            .project(self.region.clone(), u_0 as f64, v_0 as f64);
        let p_1 = self
            .sphere
            .project(self.region.clone(), u_1 as f64, v_1 as f64);
        let p_2 = self
            .sphere
            .project(self.region.clone(), u_2 as f64, v_2 as f64);
        let d = 1.0 + vec3::dot(p_0, p_1) + vec3::dot(p_1, p_2) + vec3::dot(p_2, p_0);
        (mat3::det([p_0, p_1, p_2]) / d).atan() * 2.0
    }

    fn num_sides(&self) -> usize {
        3
    }

    fn side(&self, index: usize) -> HalfEdge<Base> {
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
            region: self.region.clone(),
            start_u: self.u_0.wrapping_add_signed(d_u),
            start_v: self.v_0.wrapping_add_signed(d_v),
            dir,
        }
    }
}

impl<Base: BaseSphere> Face<Base> {
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
pub struct Vertex<Base: BaseSphere> {
    sphere: SubTriSphere<Base>,

    /// The [`BaseRegion`] which "owns" this vertex.
    region: Base::Region,

    /// The U coordinate of this vertex, specified in the local coordinate space of `region`.
    u: u32,

    /// The V coordinate of this vertex, specified in the local coordinate space of `region`.
    v: u32,
}

impl<Base: BaseSphere> Vertex<Base> {
    /// Constructs a [`Vertex`] with the given properties.
    ///
    /// This will normalize the vertex `region` to be its proper owner.
    fn new(sphere: SubTriSphere<Base>, region: Base::Region, u: u32, v: u32) -> Self {
        if region.ty().is_edge() {
            if u == 0 {
                if v == 0 {
                    Self::base(sphere, region.owner().vertex(0))
                } else if v == sphere.c && sphere.b() == sphere.c {
                    if region.ty() == BaseRegionType::Edge0 {
                        Self::center(sphere, region.owner())
                    } else {
                        Self::center(sphere, region.owner().side(2).complement().inside())
                    }
                } else if sphere.base.face_owns_edge_2(region.owner()) {
                    if region.ty() == BaseRegionType::Edge0 {
                        Self {
                            sphere,
                            region: Base::Region::new(region.owner(), BaseRegionType::Edge2),
                            u: v,
                            v: 0,
                        }
                    } else {
                        let adj = region.owner().side(2).complement();
                        debug_assert_eq!(adj.index(), 1);
                        Self {
                            sphere: sphere.clone(),
                            region: Base::Region::new(adj.inside(), BaseRegionType::Edge0),
                            u: sphere.b() - v,
                            v: sphere.c,
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
                        region: Base::Region::new(adj.inside(), BaseRegionType::Edge0),
                        u: v,
                        v: 0,
                    }
                }
            } else if u == sphere.b() {
                if v == 0 && sphere.b() == sphere.c {
                    if region.ty() == BaseRegionType::Edge0 {
                        Self::center(sphere, region.owner().side(0).complement().inside())
                    } else {
                        Self::center(sphere, region.owner())
                    }
                } else if region.ty() == BaseRegionType::Edge0 {
                    if v == sphere.c {
                        Self::base(sphere, region.owner().vertex(1))
                    } else {
                        let adj = region.owner().side(0).complement();
                        if adj.index() == 1 {
                            Self {
                                sphere: sphere.clone(),
                                region: Base::Region::new(adj.inside(), BaseRegionType::Edge0),
                                u: sphere.b() - sphere.c + v,
                                v: sphere.c,
                            }
                        } else {
                            debug_assert_eq!(adj.index(), 2);
                            Self::adjacent_1(
                                sphere.clone(),
                                adj.inside(),
                                (sphere.b() - sphere.c) + v,
                            )
                        }
                    }
                } else if v == sphere.c {
                    Self::base(sphere, region.owner().vertex(2))
                } else {
                    Self::adjacent_1(sphere.clone(), region.owner(), (sphere.b() - sphere.c) + v)
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
            let n = sphere.b() - sphere.c;
            if u + v >= n {
                debug_assert_eq!(u + v, n);
                if n == 0 {
                    todo!()
                } else if v == 0 {
                    if sphere.c == 0 {
                        Self::base(sphere, region.owner().vertex(1))
                    } else {
                        Self {
                            sphere: sphere.clone(),
                            region: Base::Region::new(region.owner(), BaseRegionType::Edge0),
                            u,
                            v: sphere.c,
                        }
                    }
                } else if u == 0 && sphere.c == 0 {
                    Self::base(sphere, region.owner().vertex(2))
                } else {
                    Self::adjacent_1(sphere, region.owner(), v)
                }
            } else if u == 0 {
                if v == 0 && sphere.c == 0 {
                    Self::base(sphere, region.owner().vertex(0))
                } else {
                    let region = if sphere.base.face_owns_edge_2(region.owner()) {
                        Base::Region::new(region.owner(), BaseRegionType::Edge2)
                    } else {
                        let adj = region.owner().side(2).complement();

                        // With the current base layout, if a base face doesn't own its edge 2,
                        // it will be the edge 0 of the adjacent face that owns it.
                        debug_assert_eq!(adj.index(), 0);
                        Base::Region::new(adj.inside(), BaseRegionType::Edge0)
                    };
                    Self {
                        sphere: sphere.clone(),
                        region,
                        u: sphere.c + v,
                        v: 0,
                    }
                }
            } else if v == 0 {
                Self {
                    sphere: sphere.clone(),
                    region: Base::Region::new(region.owner(), BaseRegionType::Edge0),
                    u,
                    v: sphere.c,
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
    fn adjacent_1(sphere: SubTriSphere<Base>, face: Base::Face, r: u32) -> Self {
        let adj = face.side(1).complement();
        if adj.index() == 0 {
            Self {
                sphere: sphere.clone(),
                region: Base::Region::new(adj.inside(), BaseRegionType::Edge0),
                u: sphere.b() - r,
                v: 0,
            }
        } else {
            debug_assert_eq!(adj.index(), 2);
            Self {
                sphere: sphere.clone(),
                region: Base::Region::new(adj.inside(), BaseRegionType::Edge2),
                u: r,
                v: sphere.c,
            }
        }
    }

    /// Gets the [`Vertex`] at the center of the given base face.
    ///
    /// This is only possible when `b == c`.
    fn center(sphere: SubTriSphere<Base>, face: Base::Face) -> Self {
        debug_assert_eq!(sphere.b(), sphere.c);
        Self {
            sphere,
            region: Base::Region::new(face, BaseRegionType::Interior),
            u: 0,
            v: 0,
        }
    }

    /// Constructs a [`Vertex`] for one of the base vertices.
    pub fn base(sphere: SubTriSphere<Base>, source: Base::Vertex) -> Self {
        let owner = sphere.base.vertex_owner(source);
        Self {
            sphere,
            region: Base::Region::new(owner, BaseRegionType::Edge0),
            u: 0,
            v: 0,
        }
    }
}

impl<Base: BaseSphere> crate::Vertex for Vertex<Base> {
    type Face = Face<Base>;
    type HalfEdge = HalfEdge<Base>;

    fn index(&self) -> usize {
        if self.region.ty().is_edge() {
            let owns_origin = self.sphere.base.face_owns_vertex_0(self.region.owner())
                && self.region.ty() == BaseRegionType::Edge0;
            self.sphere.base_vertex_index(self.region.clone())
                + self.v as usize * (self.sphere.b() as usize - 1)
                + self.u as usize
                - !owns_origin as usize
        } else if self.sphere.c < self.sphere.b() {
            let n = (self.sphere.b() - self.sphere.c) as usize;
            self.sphere.base_vertex_index(self.region.clone())
                + (self.v as usize - 1) * (2 * n - self.v as usize - 2) / 2
                + (self.u as usize - 1)
        } else {
            self.sphere.base_vertex_index(self.region.clone())
        }
    }

    fn pos(&self) -> [f64; 3] {
        self.sphere
            .project(self.region.clone(), self.u as f64, self.v as f64)
    }

    fn degree(&self) -> usize {
        if self.u == 0 && self.v == 0 { 5 } else { 6 }
    }

    fn first_outgoing(&self) -> HalfEdge<Base> {
        todo!()
    }
}

/// Represents one "side" of an edge in an [`SubTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfEdge<Base: BaseSphere> {
    sphere: SubTriSphere<Base>,

    /// The [`BaseRegion`] which "owns" this half-edge.
    region: Base::Region,

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

impl<Base: BaseSphere> crate::HalfEdge for HalfEdge<Base> {
    type Vertex = Vertex<Base>;
    type Face = Face<Base>;

    fn index(&self) -> usize {
        todo!()
    }

    fn inside(&self) -> Face<Base> {
        todo!()
    }

    fn start(&self) -> Vertex<Base> {
        Vertex::new(
            self.sphere.clone(),
            self.region.clone(),
            self.start_u,
            self.start_v,
        )
    }

    fn complement(&self) -> HalfEdge<Base> {
        todo!()
    }

    fn prev(&self) -> HalfEdge<Base> {
        todo!()
    }

    fn next(&self) -> HalfEdge<Base> {
        let mut start_u = self.start_u;
        let mut start_v = self.start_v;
        self.dir.add_to(&mut start_u, &mut start_v);
        let dir = self.dir.rotate_ccw_120();
        HalfEdge {
            sphere: self.sphere.clone(),
            region: self.region.clone(),
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
impl<Base: BaseSphere> SubTriSphere<Base> {
    /// Gets the index of the first face which is owned by the given base region.
    fn base_face_index(&self, region: Base::Region) -> usize {
        let face = region.owner();
        let num_edge_regions_before = self.base.num_owned_edges_before(face.clone())
            + (region.ty() > BaseRegionType::Edge0) as usize;
        let num_interior_regions_before =
            face.index() + (region.ty() > BaseRegionType::Interior) as usize;
        num_edge_regions_before * self.num_faces_per_edge_region()
            + num_interior_regions_before * self.num_faces_per_interior_region()
    }

    /// Gets the index of the first vertex which is owned by the given base region.
    fn base_vertex_index(&self, region: Base::Region) -> usize {
        let face = region.owner();
        let num_origin_vertices_before = self.base.num_owned_vertices_before(face.clone())
            + (self.base.face_owns_vertex_0(face.clone()) && region.ty() > BaseRegionType::Edge0)
                as usize;
        let num_edge_regions_before = self.base.num_owned_edges_before(face.clone())
            + (region.ty() > BaseRegionType::Edge0) as usize;
        let num_interior_regions_before =
            face.index() + (region.ty() > BaseRegionType::Interior) as usize;
        num_origin_vertices_before
            + num_edge_regions_before * self.num_vertices_per_edge_region()
            + num_interior_regions_before * self.num_vertices_per_interior_region()
    }

    /// Projects a point in the local coordinate space of this region to a point on the unit
    /// sphere.
    fn project(&self, region: Base::Region, u: f64, v: f64) -> Vector3 {
        let b = self.b() as f64;
        let c = self.c as f64;
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
            self.base.interpolate(weights),
        )
    }
}

impl<Base: BaseSphere> SubTriSphere<Base> {
    /// Iterates over the faces of this [`SubTriSphere`], starting with the given region.
    fn faces_from(&self, region: Base::Region) -> FaceIter<Base> {
        if region.ty().is_edge() {
            // If `c` is zero, there are no faces in an edge region. This iterator should
            // immediately go to the next region.
            FaceIter {
                sphere: self.clone(),
                region,
                u_0: 0,
                u_0_end: self.b() * (self.c > 0) as u32,
                v_0: 0,
                boundary_along_v: self.c == 0,
            }
        } else {
            let n = self.b() - self.c;
            FaceIter {
                sphere: self.clone(),
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
pub struct FaceIter<Base: BaseSphere> {
    sphere: SubTriSphere<Base>,
    region: Base::Region,
    u_0: u32,
    u_0_end: u32,
    v_0: u32,
    boundary_along_v: bool,
}

impl<Base: BaseSphere> Iterator for FaceIter<Base> {
    type Item = Face<Base>;
    fn next(&mut self) -> Option<Face<Base>> {
        loop {
            if self.u_0 < self.u_0_end {
                let res = Face {
                    sphere: self.sphere.clone(),
                    region: self.region.clone(),
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
                } else if self.v_0 + 1 < self.sphere.c {
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
                let n = self.sphere.b() - self.sphere.c;
                if self.v_0 + 1 < n {
                    self.v_0 += 1;
                    self.u_0 = 0;
                    self.u_0_end = n - self.v_0;
                    self.boundary_along_v = false;
                    continue;
                }
            }
            if let Some(region) = self.sphere.base.next_region(self.region.clone()) {
                *self = self.sphere.faces_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}

impl<Base: BaseSphere> SubTriSphere<Base> {
    /// Iterates over the vertices of this [`SubTriSphere`], starting with the given region.
    fn vertices_from(&self, region: Base::Region) -> VertexIter<Base> {
        if region.ty().is_edge() {
            let has_origin = self.base.face_owns_vertex_0(region.owner())
                && region.ty() == BaseRegionType::Edge0;
            VertexIter {
                sphere: self.clone(),
                region,
                u: !has_origin as u32,
                u_end: self.b(),
                v: 0,
            }
        } else if self.c < self.b() {
            let n = self.b() - self.c;
            VertexIter {
                sphere: self.clone(),
                region,
                u: 1,
                u_end: n - 1,
                v: 1,
            }
        } else {
            // This is a special case. When `b == c` there is exactly one interior vertex,
            // at (0, 0).
            VertexIter {
                sphere: self.clone(),
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
pub struct VertexIter<Base: BaseSphere> {
    sphere: SubTriSphere<Base>,
    region: Base::Region,
    u: u32,
    u_end: u32,
    v: u32,
}

impl<Base: BaseSphere> Iterator for VertexIter<Base> {
    type Item = Vertex<Base>;
    fn next(&mut self) -> Option<Vertex<Base>> {
        loop {
            if self.u < self.u_end {
                let res = Vertex {
                    sphere: self.sphere.clone(),
                    region: self.region.clone(),
                    u: self.u,
                    v: self.v,
                };
                self.u += 1;
                return Some(res);
            } else if self.region.ty().is_edge() {
                if self.v < self.sphere.c {
                    self.v += 1;
                    self.u = 1;
                    continue;
                }
            } else {
                let n = self.sphere.b() - self.sphere.c;
                if self.v + 1 < n {
                    self.v += 1;
                    self.u = 1;
                    self.u_end = n - self.v;
                    continue;
                }
            }
            if let Some(region) = self.sphere.base.next_region(self.region.clone()) {
                *self = self.sphere.vertices_from(region);
                continue;
            } else {
                return None;
            }
        }
    }
}
