use crate::math::{mat3, vec3, Scalar, Vector3};
use std::num::NonZeroU32;

/// A partition of the unit sphere into geodesic triangular [`Face`]s using an icosahedron
/// subdivision scheme.
///
/// More technically, this partitions the sphere by projecting the
/// geodesic polyhedron, `{3, 5+}_(b, c)`, onto it. The `b` and `c` parameters are as described
/// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct IcoSphere {
    b: NonZeroU32,
    c: u32,
}

impl IcoSphere {
    /// The base icosahedron, with no subdivisions.
    pub const fn base() -> Self {
        Self {
            b: const { NonZeroU32::new(1).unwrap() },
            c: 0,
        }
    }

    /// Constructs an [`IcoSphere`] with the given subdivision parameters.
    ///
    /// The `b` and `c` parameters are as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation). For implementation
    /// simplicity and performance, only `b ≥ c` polyhedra are directly supported. `b < c`
    /// polyhedra can be emulated by simply swapping the `b` and `c` parameters and mirroring the
    /// resulting sphere.
    pub fn new(b: NonZeroU32, c: u32) -> Self {
        assert!(
            b.get() >= c,
            "b must be greater than or equal to c. b = {}, c = {}",
            b,
            c
        );
        Self { b, c }
    }

    /// The `b` parameter of this [`IcoSphere`], as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation).
    pub fn b(&self) -> u32 {
        self.b.get()
    }

    /// The `c` parameter of this [`IcoSphere`], as described
    /// [here](https://en.wikipedia.org/wiki/Geodesic_polyhedron#Notation).
    pub fn c(&self) -> u32 {
        self.c
    }

    /// The number of [`Face`]s each triangle on the base icosahedron is subdivided into.
    ///
    /// This is also known as the "triangulation number" of the icosphere.
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

impl crate::Sphere for IcoSphere {
    type Face = Face;
    type Vertex = Vertex;
    type HalfEdge = HalfEdge;

    fn num_faces(&self) -> usize {
        self.num_divisions() * 20
    }

    #[allow(refining_impl_trait)]
    fn faces(&self) -> FaceIter {
        self.faces_from(BaseRegion::FIRST)
    }

    fn num_vertices(&self) -> usize {
        self.num_divisions() * 10 + 2
    }

    #[allow(refining_impl_trait)]
    fn vertices(&self) -> VertexIter {
        self.vertices_from(BaseRegion::FIRST)
    }
}

/// Represents a face in an [`IcoSphere`].
///
/// This is always a geodesic triangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Face {
    sphere: IcoSphere,

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
            let n = (self.sphere.b() - self.sphere.c) as usize;
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

    fn first_boundary(&self) -> HalfEdge {
        HalfEdge {
            sphere: self.sphere,
            region: self.region,
            start_u: self.u_0,
            start_v: self.v_0,
            dir: if self.boundary_along_v {
                HalfEdgeDir::Vp
            } else {
                HalfEdgeDir::Up
            },
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

/// Represents a vertex in an [`IcoSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Vertex {
    sphere: IcoSphere,

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
    fn new(sphere: IcoSphere, region: BaseRegion, u: u32, v: u32) -> Self {
        if region.ty().is_edge() {
            if u == 0 {
                if v == 0 {
                    Self::base(sphere, region.owner().vertex_owner(0))
                } else if v == sphere.c && sphere.b() == sphere.c {
                    if region.ty() == BaseRegionType::Edge0 {
                        Self::center(sphere, region.owner())
                    } else {
                        Self::center(sphere, region.owner().adjacent(2).0)
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
                        let (adj, adj_index) = region.owner().adjacent(2);
                        debug_assert_eq!(adj_index, 1);
                        Self {
                            sphere,
                            region: BaseRegion::new(adj, BaseRegionType::Edge0),
                            u: sphere.b() - v,
                            v: sphere.c,
                        }
                    }
                } else {
                    debug_assert_eq!(region.ty(), BaseRegionType::Edge0);
                    let (adj, adj_index) = region.owner().adjacent(2);

                    // With the current base layout, if a base face doesn't own its edge 2,
                    // it will be the edge 0 of the adjacent face that owns it.
                    debug_assert_eq!(adj_index, 0);
                    Self {
                        sphere,
                        region: BaseRegion::new(adj, BaseRegionType::Edge0),
                        u: v,
                        v: 0,
                    }
                }
            } else if u == sphere.b() {
                if v == 0 && sphere.b() == sphere.c {
                    if region.ty() == BaseRegionType::Edge0 {
                        Self::center(sphere, region.owner().adjacent(0).0)
                    } else {
                        Self::center(sphere, region.owner())
                    }
                } else if region.ty() == BaseRegionType::Edge0 {
                    if v == sphere.c {
                        Self::base(sphere, region.owner().vertex_owner(1))
                    } else {
                        let (adj, adj_index) = region.owner().adjacent(0);
                        if adj_index == 1 {
                            Self {
                                sphere,
                                region: BaseRegion::new(adj, BaseRegionType::Edge0),
                                u: sphere.b() - sphere.c + v,
                                v: sphere.c,
                            }
                        } else {
                            debug_assert_eq!(adj_index, 2);
                            Self::adjacent_1(sphere, adj, (sphere.b() - sphere.c) + v)
                        }
                    }
                } else if v == sphere.c {
                    Self::base(sphere, region.owner().vertex_owner(2))
                } else {
                    Self::adjacent_1(sphere, region.owner(), (sphere.b() - sphere.c) + v)
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
                        Self::base(sphere, region.owner().vertex_owner(1))
                    } else {
                        Self {
                            sphere,
                            region: BaseRegion::new(region.owner(), BaseRegionType::Edge0),
                            u,
                            v: sphere.c,
                        }
                    }
                } else if u == 0 && sphere.c == 0 {
                    Self::base(sphere, region.owner().vertex_owner(2))
                } else {
                    Self::adjacent_1(sphere, region.owner(), v)
                }
            } else if u == 0 {
                if v == 0 && sphere.c == 0 {
                    Self::base(sphere, region.owner().vertex_owner(0))
                } else {
                    let region = if region.owner().owns_edge_2() {
                        BaseRegion::new(region.owner(), BaseRegionType::Edge2)
                    } else {
                        let (adj, adj_index) = region.owner().adjacent(2);

                        // With the current base layout, if a base face doesn't own its edge 2,
                        // it will be the edge 0 of the adjacent face that owns it.
                        debug_assert_eq!(adj_index, 0);
                        BaseRegion::new(adj, BaseRegionType::Edge0)
                    };
                    Self {
                        sphere,
                        region,
                        u: sphere.c + v,
                        v: 0,
                    }
                }
            } else if v == 0 {
                Self {
                    sphere,
                    region: BaseRegion::new(region.owner(), BaseRegionType::Edge0),
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

    /// Constructs a [`Vertex`] on the edge 1 border of the given [`BaseFace`].
    ///
    /// `r` is the height of the vertex above `face`'s vertex 1. The caller must ensure
    /// `0 < r < b`.
    fn adjacent_1(sphere: IcoSphere, face: BaseFace, r: u32) -> Self {
        let (adj, adj_index) = face.adjacent(1);
        if adj_index == 0 {
            Self {
                sphere,
                region: BaseRegion::new(adj, BaseRegionType::Edge0),
                u: sphere.b() - r,
                v: 0,
            }
        } else {
            debug_assert_eq!(adj_index, 2);
            Self {
                sphere,
                region: BaseRegion::new(adj, BaseRegionType::Edge2),
                u: r,
                v: sphere.c,
            }
        }
    }

    /// Gets the [`Vertex`] at the center of the given [`BaseFace`].
    ///
    /// This is only possible when `b == c`.
    fn center(sphere: IcoSphere, face: BaseFace) -> Self {
        debug_assert_eq!(sphere.b(), sphere.c);
        Self {
            sphere,
            region: BaseRegion::new(face, BaseRegionType::Interior),
            u: 0,
            v: 0,
        }
    }

    /// Constructs a [`Vertex`] for one of the 12 base vertices, identified by its owning
    /// [`BaseFace`].
    const fn base(sphere: IcoSphere, owner: BaseFace) -> Self {
        Self {
            sphere,
            region: BaseRegion::new(owner, BaseRegionType::Edge0),
            u: 0,
            v: 0,
        }
    }
}

impl crate::Vertex for Vertex {
    type Face = Face;
    type HalfEdge = HalfEdge;

    fn index(&self) -> usize {
        if self.region.ty().is_edge() {
            let owns_origin =
                self.region.owner().owns_vertex_0() && self.region.ty() == BaseRegionType::Edge0;
            self.sphere.base_vertex_index(self.region)
                + self.v as usize * (self.sphere.b() as usize - 1)
                + self.u as usize
                - !owns_origin as usize
        } else if self.sphere.c < self.sphere.b() {
            let n = (self.sphere.b() - self.sphere.c) as usize;
            self.sphere.base_vertex_index(self.region)
                + (self.v as usize - 1) * (2 * n - self.v as usize - 2) / 2
                + (self.u as usize - 1)
        } else {
            self.sphere.base_vertex_index(self.region)
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

/// Represents one "side" of an edge in an [`IcoSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfEdge {
    sphere: IcoSphere,

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

/// Represents one of the 20 faces of the base icosphere.
///
/// For any particular subdivision of the icosahedron, every vertex and face is "owned" by some
/// [`BaseFace`]. This determines the order in which they are indexed/iterated.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct BaseFace(u8);

impl BaseFace {
    /// Gets the positions of the vertices of this base face.
    pub fn vertices_pos(&self) -> [Vector3; 3] {
        unsafe {
            let indices = BASE_INDICES.get_unchecked(self.0 as usize);
            [
                *BASE_VERTS.get_unchecked(indices[0] as usize),
                *BASE_VERTS.get_unchecked(indices[1] as usize),
                *BASE_VERTS.get_unchecked(indices[2] as usize),
            ]
        }
    }

    /// Identifies the [owner](Vertex::base) of a specified vertex of this base face.
    pub fn vertex_owner(self, index: u8) -> BaseFace {
        unsafe { Self::VERTEX_OWNER.get_unchecked(self.0 as usize)[index as usize] }
    }

    /// Table used to implement [`Self::vertex_owner`].
    const VERTEX_OWNER: &'static [[BaseFace; 3]; NUM_BASE_FACES] = &const {
        // Assign ownership of each base vertex to the first base face that contains it as it's 0
        // vertex.
        let mut vert_owner: [Option<BaseFace>; NUM_BASE_VERTS] = [None; NUM_BASE_VERTS];
        let mut i = 0;
        while i < NUM_BASE_FACES {
            let v_0 = BASE_INDICES[i][0] as usize;
            if vert_owner[v_0].is_none() {
                vert_owner[v_0] = Some(BaseFace(i as u8));
            }
            i += 1;
        }

        // Verify that every base vertex has an owner
        let mut j = 0;
        while j < NUM_BASE_VERTS {
            assert!(
                vert_owner[j].is_some(),
                "base vertex does not have an owner"
            );
            j += 1;
        }

        // Convert into ownership table
        let mut table: [[BaseFace; 3]; NUM_BASE_FACES] = [[BaseFace(0); 3]; NUM_BASE_FACES];
        let mut i = 0;
        while i < NUM_BASE_FACES {
            table[i] = [
                vert_owner[BASE_INDICES[i][0] as usize].unwrap(),
                vert_owner[BASE_INDICES[i][1] as usize].unwrap(),
                vert_owner[BASE_INDICES[i][2] as usize].unwrap(),
            ];
            i += 1;
        }
        table
    };

    /// Identifies the [`BaseFace`] adjacent to this one along the specified edge. Also returns the
    /// index of that edge relative to the *returned* [`BaseFace`].
    pub fn adjacent(self, index: u8) -> (BaseFace, u8) {
        unsafe { Self::ADJACENT.get_unchecked(self.0 as usize)[index as usize] }
    }

    /// Table used to implement [`Self::adjacent`].
    const ADJACENT: &'static [[(BaseFace, u8); 3]; NUM_BASE_FACES] = &const {
        // Build adjacent mapping for potential edges
        let mut adjacent: [Option<(BaseFace, u8)>; NUM_BASE_VERTS * NUM_BASE_VERTS] =
            [None; NUM_BASE_VERTS * NUM_BASE_VERTS];
        let mut i = 0;
        while i < NUM_BASE_FACES {
            let v_0 = BASE_INDICES[i][0] as usize;
            let v_1 = BASE_INDICES[i][1] as usize;
            let v_2 = BASE_INDICES[i][2] as usize;
            let edge_0 = v_0 * NUM_BASE_VERTS + v_1;
            let edge_1 = v_1 * NUM_BASE_VERTS + v_2;
            let edge_2 = v_2 * NUM_BASE_VERTS + v_0;
            assert!(adjacent[edge_0].is_none(), "duplicate edge detected");
            assert!(adjacent[edge_1].is_none(), "duplicate edge detected");
            assert!(adjacent[edge_2].is_none(), "duplicate edge detected");
            adjacent[edge_0] = Some((BaseFace(i as u8), 0));
            adjacent[edge_1] = Some((BaseFace(i as u8), 1));
            adjacent[edge_2] = Some((BaseFace(i as u8), 2));
            i += 1;
        }

        // Convert to adjacency table
        let mut table: [[(BaseFace, u8); 3]; NUM_BASE_FACES] =
            [[(BaseFace(0), 0); 3]; NUM_BASE_FACES];
        let mut i = 0;
        while i < NUM_BASE_FACES {
            let v_0 = BASE_INDICES[i][0] as usize;
            let v_1 = BASE_INDICES[i][1] as usize;
            let v_2 = BASE_INDICES[i][2] as usize;
            table[i] = [
                adjacent[v_1 * NUM_BASE_VERTS + v_0].unwrap(),
                adjacent[v_2 * NUM_BASE_VERTS + v_1].unwrap(),
                adjacent[v_0 * NUM_BASE_VERTS + v_2].unwrap(),
            ];
            i += 1;
        }
        table
    };

    /// Indicates whether this face "owns" its first vertex.
    ///
    /// This is the only vertex which may be owned by the base face.
    pub fn owns_vertex_0(self) -> bool {
        (Self::OWNS_VERTEX_0 >> self.0) & 1 == 1
    }

    /// Bitset used to implement [`Self::owns_origin`].
    const OWNS_VERTEX_0: u32 = const {
        let mut owns_vertex_0 = 0;
        let mut i = 0;
        while i < NUM_BASE_FACES {
            if Self::VERTEX_OWNER[i][0].0 as usize == i {
                owns_vertex_0 |= 1 << i;
            }
            i += 1;
        }
        owns_vertex_0
    };

    /// Indicates whether this face "owns" the region corresponding to its second edge.
    pub fn owns_edge_2(self) -> bool {
        (Self::OWNS_EDGE_2 >> self.0) & 1 == 1
    }

    /// Bitset used to implement [`Self::owns_edge_2`].
    const OWNS_EDGE_2: u32 = const {
        // Assign ownership of each edge. First by edge 0, then optionally by edge 2.
        let mut edge_has_owner: [bool; NUM_BASE_VERTS * NUM_BASE_VERTS] =
            [false; NUM_BASE_VERTS * NUM_BASE_VERTS];
        let mut num_owned_edges = 0;
        let mut i = 0;
        while i < NUM_BASE_FACES {
            let v_0 = BASE_INDICES[i][0] as usize;
            let v_1 = BASE_INDICES[i][1] as usize;
            let edge_0 = if v_0 < v_1 {
                v_0 * NUM_BASE_VERTS + v_1
            } else {
                v_1 * NUM_BASE_VERTS + v_0
            };
            assert!(!edge_has_owner[edge_0], "edge already has an owner");
            edge_has_owner[edge_0] = true;
            num_owned_edges += 1;
            i += 1;
        }
        let mut owns_edge_2 = 0;
        let mut i = 0;
        while i < NUM_BASE_FACES {
            let v_2 = BASE_INDICES[i][2] as usize;
            let v_0 = BASE_INDICES[i][0] as usize;
            let edge_2 = if v_0 < v_2 {
                v_0 * NUM_BASE_VERTS + v_2
            } else {
                v_2 * NUM_BASE_VERTS + v_0
            };
            if !edge_has_owner[edge_2] {
                edge_has_owner[edge_2] = true;
                num_owned_edges += 1;
                owns_edge_2 |= 1 << i;
            }
            i += 1;
        }

        // Verify that every edge has an owner
        assert!(
            num_owned_edges == NUM_BASE_EDGES,
            "not all edges have an owner"
        );
        owns_edge_2
    };

    /// Interpolation function for a base face.
    ///
    /// Suppose the vertices of a base face are `v₀`, `v₁`, and `v₂`. The interpolation function
    /// is given weights `w₀`, `w₁`, and `w₂` such that `wᵢ ≥ 0` and `w₀ + w₁ + w₂ = 1`. It
    /// will compute an approximation of the "spherical" interpolation of the three vertices
    /// with respect to the given weights, `r = c₀ v₀ + c₁ v₁ + c₂ v₂` and return the coefficients
    /// `[c₀, c₁, c₂]`.
    ///
    /// It is guaranteed that:
    ///  * `r` is on the unit sphere
    ///  * If `wᵢ` is `1`, `cᵢ` will also be `1`.
    ///
    /// See [this paper](https://mathweb.ucsd.edu/~sbuss/ResearchWeb/spheremean/paper.pdf) for
    /// a precise definition of spherical interpolation (which we only approximate here, for
    /// performance reasons).
    pub fn interpolate(weights: [Scalar; 3]) -> [Scalar; 3] {
        let [w_0, w_1, w_2] = weights;

        // Compute initial unnormalized coefficients. `cᵢ = (wᵢ * ANGLE).sin()` assures that this
        // behaves like a perfect SLERP (https://en.wikipedia.org/wiki/Slerp) when one of
        // the weights is `0`.
        // TODO: it might be possible to create a faster implementation of `(x * ANGLE).sin()`.
        // This could also make the result portable, since `sin` isn't.
        let c_0 = (w_0 * ANGLE).sin();
        let c_1 = (w_1 * ANGLE).sin();
        let c_2 = (w_2 * ANGLE).sin();

        // Normalize coefficients such that, when applied to the 3 vertices of the base face,
        // the result will be on the unit sphere
        let norm_sqr =
            c_0 * c_0 + c_1 * c_1 + c_2 * c_2 + 2.0 * DOT * (c_0 * c_1 + c_1 * c_2 + c_2 * c_0);
        let norm = norm_sqr.sqrt();
        let c_0 = c_0 / norm;
        let c_1 = c_1 / norm;
        let c_2 = c_2 / norm;
        return [c_0, c_1, c_2];

        /// Dot product of two adjacent vertices of the base icosahedron.
        const DOT: Scalar = C_1;

        /// Angle between two adjacent vertices of the base icosahedron.
        const ANGLE: Scalar = 1.1071487177940904;
    }
}

/// Represents one of the 50 "regions" on an icosphere.
///
/// There is a region corresponding to every base face and base edge. Each region defines a local
/// UV coordinate space, shaped like a parallelogram for edges, and a triangle for faces.
///
/// Every base region is "owned" by exactly one base face. Base faces always own the region
/// corresponding to their first edge, their interior, and possibly the region corresponding to
/// their third edge.
///
/// Every vertex, face and half-edge is "owned" by exactly one base region.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct BaseRegion(u8);

impl BaseRegion {
    /// The first [`BaseRegion`].
    pub const FIRST: BaseRegion = BaseRegion(0);

    /// Gets the [`BaseRegion`] with the given owner and type.
    pub const fn new(owner: BaseFace, ty: BaseRegionType) -> Self {
        BaseRegion((owner.0 << 2) | (ty as u8))
    }

    /// The [`BaseFace`] which owns this region.
    pub const fn owner(self) -> BaseFace {
        BaseFace(self.0 >> 2)
    }

    /// The type of this region.
    pub const fn ty(self) -> BaseRegionType {
        unsafe { std::mem::transmute(self.0 & 3) }
    }
}

/// The general type of a [`BaseRegion`].
#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum BaseRegionType {
    /// Corresponds to the first edge of [`BaseRegion::owner`].
    Edge0 = 0,

    /// Corresponds to the interior of [`BaseRegion::owner`].
    Interior = 1,

    /// Corresponds to the third edge of [`BaseRegion::owner`].
    Edge2 = 3,
}

impl BaseRegionType {
    /// Indicates whether this region corresponds to an edge.
    fn is_edge(self) -> bool {
        self != Self::Interior
    }
}

impl IcoSphere {
    /// Gets the index of the first face which is owned by the given base region.
    fn base_face_index(&self, region: BaseRegion) -> usize {
        let face = region.owner();
        let before_mask = (1 << face.0) - 1;
        let num_edge_regions_before = face.0 as usize
            + (BaseFace::OWNS_EDGE_2 & before_mask).count_ones() as usize
            + (region.ty() > BaseRegionType::Edge0) as usize;
        let num_interior_regions_before =
            face.0 as usize + (region.ty() > BaseRegionType::Interior) as usize;
        num_edge_regions_before * self.num_faces_per_edge_region()
            + num_interior_regions_before * self.num_faces_per_interior_region()
    }

    /// Gets the index of the first vertex which is owned by the given base region.
    fn base_vertex_index(&self, region: BaseRegion) -> usize {
        let face = region.owner();
        let before_mask = (1 << face.0) - 1;
        let num_origin_vertices_before = (BaseFace::OWNS_VERTEX_0 & before_mask).count_ones()
            as usize
            + (face.owns_vertex_0() && region.ty() > BaseRegionType::Edge0) as usize;
        let num_edge_regions_before = face.0 as usize
            + (BaseFace::OWNS_EDGE_2 & before_mask).count_ones() as usize
            + (region.ty() > BaseRegionType::Edge0) as usize;
        let num_interior_regions_before =
            face.0 as usize + (region.ty() > BaseRegionType::Interior) as usize;
        num_origin_vertices_before
            + num_edge_regions_before * self.num_vertices_per_edge_region()
            + num_interior_regions_before * self.num_vertices_per_interior_region()
    }

    /// Projects a point in the local coordinate space of this region to a point on the unit
    /// sphere.
    fn project(&self, region: BaseRegion, u: f64, v: f64) -> Vector3 {
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
                    let (adj, adj_index) = region.owner().adjacent(0);
                    let w_2 = -w_left;
                    let w_0 = ((u + v) * b + v * c) / norm_sqr;
                    let w_1 = 1.0 - w_0 - w_2;
                    let mut weights = [w_0, w_1, w_2];
                    weights.rotate_right(adj_index as usize);
                    (adj, weights)
                }
            } else if w_left <= 0.0 {
                let w_1 = -w_left;
                let w_2 = ((u + v) * b + v * c) / norm_sqr;
                (region.owner(), [1.0 - w_1 - w_2, w_1, w_2])
            } else {
                let (adj, adj_index) = region.owner().adjacent(2);
                let w_1 = (u * b + (u + v) * c) / norm_sqr;
                let w_2 = w_left;
                let w_0 = 1.0 - w_1 - w_2;
                let mut weights = [w_0, w_1, w_2];
                weights.rotate_right(adj_index as usize);
                (adj, weights)
            }
        } else {
            let v = v + c;
            let w_1 = (u * b + (u + v) * c) / norm_sqr;
            let w_2 = (v * b - u * c) / norm_sqr;
            (region.owner(), [1.0 - w_1 - w_2, w_1, w_2])
        };
        mat3::apply(face.vertices_pos(), BaseFace::interpolate(weights))
    }
}

impl IcoSphere {
    /// Iterates over the faces of this [`IcoSphere`], starting with the given region.
    fn faces_from(&self, region: BaseRegion) -> FaceIter {
        if region.ty().is_edge() {
            // If `c` is zero, there are no faces in an edge region. This iterator should
            // immediately go to the next region.
            FaceIter {
                sphere: *self,
                region,
                u_0: 0,
                u_0_end: self.b() * (self.c > 0) as u32,
                v_0: 0,
                boundary_along_v: self.c == 0,
            }
        } else {
            let n = self.b() - self.c;
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

/// An iterator over the faces of an [`IcoSphere`].
pub struct FaceIter {
    sphere: IcoSphere,
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
            let next_region;
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
                } else if self.v_0 + 1 < self.sphere.c {
                    self.v_0 += 1;
                    self.u_0 = 0;
                    self.u_0_end = self.sphere.b();
                    self.boundary_along_v = false;
                    continue;
                } else {
                    next_region = self.region.0 + 1;
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
                } else if self.region.owner().owns_edge_2() {
                    next_region = self.region.0 + 2;
                } else {
                    next_region = self.region.0 + 3;
                }
            }
            if (next_region as usize) < NUM_BASE_FACES << 2 {
                *self = self.sphere.faces_from(BaseRegion(next_region));
                continue;
            } else {
                return None;
            }
        }
    }
}

impl IcoSphere {
    /// Iterates over the vertices of this [`IcoSphere`], starting with the given region.
    fn vertices_from(&self, region: BaseRegion) -> VertexIter {
        if region.ty().is_edge() {
            let has_origin = region.owner().owns_vertex_0() && region.ty() == BaseRegionType::Edge0;
            VertexIter {
                sphere: *self,
                region,
                u: !has_origin as u32,
                u_end: self.b(),
                v: 0,
            }
        } else if self.c < self.b() {
            let n = self.b() - self.c;
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

/// An iterator over the vertices of an [`IcoSphere`].
pub struct VertexIter {
    sphere: IcoSphere,
    region: BaseRegion,
    u: u32,
    u_end: u32,
    v: u32,
}

impl Iterator for VertexIter {
    type Item = Vertex;
    fn next(&mut self) -> Option<Vertex> {
        loop {
            let next_region;
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
                if self.v < self.sphere.c {
                    self.v += 1;
                    self.u = 1;
                    continue;
                } else {
                    next_region = self.region.0 + 1;
                }
            } else {
                let n = self.sphere.b() - self.sphere.c;
                if self.v + 1 < n {
                    self.v += 1;
                    self.u = 1;
                    self.u_end = n - self.v;
                    continue;
                } else if self.region.owner().owns_edge_2() {
                    next_region = self.region.0 + 2;
                } else {
                    next_region = self.region.0 + 3;
                }
            }
            if (next_region as usize) < NUM_BASE_FACES << 2 {
                *self = self.sphere.vertices_from(BaseRegion(next_region));
                continue;
            } else {
                return None;
            }
        }
    }
}

/// The number of faces on the base icosahedron.
const NUM_BASE_FACES: usize = 20;

/// The number of vertices on the base icosahedron.
const NUM_BASE_VERTS: usize = 12;

/// The number of edges on the base icosahedron.
const NUM_BASE_EDGES: usize = 30;

/// The positions of the vertices of the base icosahedron.
static BASE_VERTS: [[f64; 3]; NUM_BASE_VERTS] = [
    // Top apex
    [0.0, 0.0, 1.0],
    // Top ring
    [C_0, 0.0, C_1],
    [C_2, C_3, C_1],
    [-C_4, C_5, C_1],
    [-C_4, -C_5, C_1],
    [C_2, -C_3, C_1],
    // Bottom ring
    [C_4, -C_5, -C_1],
    [C_4, C_5, -C_1],
    [-C_2, C_3, -C_1],
    [-C_0, 0.0, -C_1],
    [-C_2, -C_3, -C_1],
    // Bottom apex
    [0.0, 0.0, -1.0],
];

/// The indices of the face vertices of the base icosahedron.
static BASE_INDICES: [[u32; 3]; NUM_BASE_FACES] = [
    // Top cap
    [0, 1, 2],
    [0, 2, 3],
    [0, 3, 4],
    [0, 4, 5],
    [0, 5, 1],
    // Central ring
    [5, 6, 1],
    [6, 7, 1],
    [1, 7, 2],
    [7, 8, 2],
    [2, 8, 3],
    [8, 9, 3],
    [3, 9, 4],
    [9, 10, 4],
    [4, 10, 5],
    [10, 6, 5],
    // Bottom cap
    [11, 6, 10],
    [11, 10, 9],
    [11, 9, 8],
    [11, 8, 7],
    [11, 7, 6],
];

/// `sqrt(4 / 5)`
const C_0: f64 = 0.8944271909999159;

/// `sqrt(1 / 5)`
const C_1: f64 = 0.4472135954999579;

/// `(5 - sqrt(5)) / 10`
const C_2: f64 = 0.276393202250021;

/// `sqrt((5 + sqrt(5)) / 10)`
const C_3: f64 = 0.8506508083520399;

/// `(5 + sqrt(5)) / 10`
const C_4: f64 = 0.7236067977499789;

/// `sqrt((5 - sqrt(5)) / 10)`
const C_5: f64 = 0.5257311121191336;
