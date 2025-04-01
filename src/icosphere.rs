use crate::math::{Matrix3, Scalar, Vector3, mat3, vec3};
use std::num::NonZeroU32;

/// A partition of the unit sphere into geodesic triangular [`Face`]s using an icosahedron-like
/// subdivision scheme.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct IcoSphere {
    divisions: NonZeroU32,
}

impl IcoSphere {
    /// Constructs an [`IcoSphere`] with the given number of subdivisions.
    ///
    /// The `divisions` parameters specifies how many segments each edge on the base icosahedron
    /// should be divided into.
    pub fn new(divisions: NonZeroU32) -> Self {
        Self { divisions }
    }

    /// The number of segments each edge on the base icosahedron is divided into.
    pub fn divisions(&self) -> NonZeroU32 {
        self.divisions
    }
}

impl crate::Sphere for IcoSphere {
    type Face = Face;
    type Vertex = Vertex;
    type HalfEdge = HalfEdge;

    #[allow(refining_impl_trait)]
    fn faces(&self) -> FaceIter {
        FaceIter {
            sphere: *self,
            base: 0,
            u_0: 0,
            v_0: 0,
            boundary_along_v: false,
        }
    }

    #[allow(refining_impl_trait)]
    fn vertices(&self) -> VertexIter {
        VertexIter {
            sphere: *self,
            base: 0,
            u: 0,
            v: 0,
        }
    }
}

/// An iterator over the faces of an [`IcoSphere`].
pub struct FaceIter {
    sphere: IcoSphere,
    base: u8,
    u_0: u32,
    v_0: u32,
    boundary_along_v: bool,
}

impl Iterator for FaceIter {
    type Item = Face;
    fn next(&mut self) -> Option<Face> {
        if self.base as usize >= NUM_BASE_FACES {
            return None;
        }
        let base = BaseFace(self.base);
        let res = Face {
            sphere: self.sphere,
            base,
            u_0: self.u_0,
            v_0: self.v_0,
            boundary_along_v: self.boundary_along_v,
        };
        self.u_0 += 1;
        if self.u_0 + self.v_0 >= self.sphere.divisions().get() {
            self.v_0 += u32::from(self.boundary_along_v);
            self.boundary_along_v = !self.boundary_along_v;
            self.u_0 = u32::from(self.boundary_along_v);
            if self.u_0 + self.v_0 >= self.sphere.divisions().get() {
                self.base += 1;
                self.u_0 = 0;
                self.v_0 = 0;
                self.boundary_along_v = false;
            }
        }
        Some(res)
    }
}

/// An iterator over the vertices of an [`IcoSphere`].
pub struct VertexIter {
    sphere: IcoSphere,
    base: u8,
    u: u32,
    v: u32,
}

impl Iterator for VertexIter {
    type Item = Vertex;
    fn next(&mut self) -> Option<Vertex> {
        if self.base as usize >= NUM_BASE_FACES {
            return None;
        }
        let base = BaseFace(self.base);
        let res = Vertex {
            sphere: self.sphere,
            base,
            u: self.u,
            v: self.v,
        };
        self.u += 1;
        if self.u + self.v >= self.sphere.divisions().get() {
            self.u = u32::from(!base.owns_v_edge());
            self.v += 1;
            while self.u + self.v >= self.sphere.divisions().get() {
                self.base += 1;
                if self.base as usize >= NUM_BASE_FACES {
                    break;
                }
                self.u = u32::from(!BaseFace(self.base).owns_origin());
                self.v = 0;
            }
        }
        Some(res)
    }
}

/// Represents a face in an [`IcoSphere`].
///
/// This is always a geodesic triangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Face {
    sphere: IcoSphere,

    /// The [`BaseFace`] this face belongs to.
    base: BaseFace,

    /// The U coordinate of the first vertex of this face on `base`.
    u_0: u32,

    /// The V coordinate of the first vertex of this face on `base`.
    v_0: u32,

    /// If `true`, the first [`HalfEdge`] boundary of this face will go in the +V direction.
    /// Otherwise, it will go in the +U direction.
    boundary_along_v: bool,
}

impl crate::Face for Face {
    type Vertex = Vertex;
    type HalfEdge = HalfEdge;

    fn index(&self) -> usize {
        todo!()
    }

    fn area(&self) -> f64 {
        todo!()
    }

    fn center(&self) -> Vector3 {
        todo!()
    }

    fn num_sides(&self) -> usize {
        3
    }

    fn first_boundary(&self) -> HalfEdge {
        HalfEdge {
            sphere: self.sphere,
            base: self.base,
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

/// Represents a vertex in an [`IcoSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Vertex {
    sphere: IcoSphere,

    /// The [`BaseFace`] which "owns" this vertex.
    ///
    /// The order that vertices are indexed/iterated is determined by their "owning" base face.
    /// Each base face:
    ///   * Sometimes owns its local origin vertex (the vertex at `U = 0` and `V = 0`).
    ///   * Always owns the vertices along its U edge (i.e. vertices with  `U = 1..` and `V = 0`).
    ///   * Sometimes owns the vertices along its V edge (i.e. vertices with `U = 0` and `V = 1..`).
    ///   * Never owns the vertices along its far edge (i.e. vertices with `U + V = N`).
    ///   * Always owns the vertices in its interior (i.e. vertices with `U > 0`, `V > 0` and
    ///     `U + V < N`).
    base: BaseFace,

    /// The U coordinate of this vertex on `base`.
    u: u32,

    /// The V coordinate of this vertex on `base`.
    v: u32,
}

impl Vertex {
    /// Constructs a [`Vertex`] with the given properties.
    ///
    /// This will normalize the vertex `base` to be its proper owner.
    fn new(sphere: IcoSphere, base: BaseFace, u: u32, v: u32) -> Self {
        if u == 0 {
            if v == 0 {
                Self {
                    sphere,
                    base: base.vertex_owner(0),
                    u: 0,
                    v: 0,
                }
            } else if v == sphere.divisions().get() {
                Self {
                    sphere,
                    base: base.vertex_owner(2),
                    u: 0,
                    v: 0,
                }
            } else if base.owns_v_edge() {
                Self { sphere, base, u, v }
            } else {
                let (adj, adj_index) = base.adjacent(2);

                // With the current configuration of the base icosahedron. Every base face which
                // doesn't own its V edge is adjacent to the U edge of the base face that does.
                debug_assert_eq!(adj_index, 0);
                Self {
                    sphere,
                    base: adj,
                    u: v,
                    v: 0,
                }
            }
        } else if u + v == sphere.divisions().get() {
            if v == 0 {
                Self {
                    sphere,
                    base: base.vertex_owner(1),
                    u: 0,
                    v: 0,
                }
            } else {
                let (adj, adj_index) = base.adjacent(1);
                if adj_index == 0 {
                    Self {
                        sphere,
                        base: adj,
                        u,
                        v: 0,
                    }
                } else {
                    debug_assert_eq!(adj_index, 2);
                    Self {
                        sphere,
                        base: adj,
                        u: 0,
                        v,
                    }
                }
            }
        } else {
            debug_assert!(u + v < sphere.divisions().get());
            Self { sphere, base, u, v }
        }
    }
}

impl crate::Vertex for Vertex {
    type Face = Face;
    type HalfEdge = HalfEdge;

    fn index(&self) -> usize {
        let divs = self.sphere.divisions().get() as usize;
        let u = self.u as usize;
        let v = self.v as usize;
        let not_owns_origin = !self.base.owns_origin() as usize;
        let not_owns_v_edge = !self.base.owns_v_edge() as usize;
        self.base.base_vertex_index(divs) + u + v * (2 * divs - v + 1) / 2
            - not_owns_origin
            - not_owns_v_edge * v
    }

    fn pos(&self) -> [f64; 3] {
        let rel_u = self.u as f64 / self.sphere.divisions().get() as f64;
        let rel_v = self.v as f64 / self.sphere.divisions().get() as f64;
        let weights = [1.0 - rel_u - rel_v, rel_u, rel_v];
        let coeffs = BaseFace::interpolate(weights);
        mat3::apply(self.base.vertices_pos(), coeffs)
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

    /// The [`BaseFace`] this face belongs to.
    base: BaseFace,

    /// The U coordinate of the start vertex on `base`.
    start_u: u32,

    /// The V coordinate of the start vertex on `base`.
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
        Vertex::new(self.sphere, self.base, self.start_u, self.start_v)
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
            base: self.base,
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

/// Represents one of the 20 faces of the base (0-subdivision) icosphere.
///
/// All [`Face`]s of any subdivision level belong to some [`BaseFace`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct BaseFace(u8);

impl BaseFace {
    /// Gets the index of the first vertex owned by this base face.
    pub fn base_vertex_index(&self, divs: usize) -> usize {
        let index = self.0 as usize;
        let before_mask = (1 << index) - 1;
        (Self::OWNS_ORIGIN & before_mask).count_ones() as usize
            + (Self::OWNS_V_EDGE & before_mask).count_ones() as usize * (divs - 1)
            + divs * (divs - 1) * index / 2
    }

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

    /// Indicates whether this face [owns](Vertex::base) the vertex at its local origin (i.e. the
    /// vertex with `U = 0` and `V = 0`).
    pub fn owns_origin(self) -> bool {
        (Self::OWNS_ORIGIN >> self.0) & 1 == 1
    }

    /// Bitset used to implement [`Self::owns_origin`].
    const OWNS_ORIGIN: u32 = const {
        let mut owns_origin = 0;
        let mut i = 0;
        while i < NUM_BASE_FACES {
            if Self::VERTEX_OWNER[i][0].0 as usize == i {
                owns_origin |= 1 << i;
            }
            i += 1;
        }
        owns_origin
    };

    /// Indicates whether this face [owns](Vertex::base) the vertices along its V edge (i.e.
    /// vertices with `U = 0` and `V = 1..`).
    pub fn owns_v_edge(self) -> bool {
        (Self::OWNS_V_EDGE >> self.0) & 1 == 1
    }

    /// Bitset used to implement [`Self::owns_v_edge`].
    const OWNS_V_EDGE: u32 = const {
        // Assign ownership of each edge. First by U edge, then optionally by V edge.
        let mut edge_has_owner: [bool; NUM_BASE_VERTS * NUM_BASE_VERTS] =
            [false; NUM_BASE_VERTS * NUM_BASE_VERTS];
        let mut num_owned_edges = 0;
        let mut i = 0;
        while i < NUM_BASE_FACES {
            let v_0 = BASE_INDICES[i][0] as usize;
            let v_1 = BASE_INDICES[i][1] as usize;
            let edge = if v_0 < v_1 {
                v_0 * NUM_BASE_VERTS + v_1
            } else {
                v_1 * NUM_BASE_VERTS + v_0
            };
            assert!(!edge_has_owner[edge], "edge already has an owner");
            edge_has_owner[edge] = true;
            num_owned_edges += 1;
            i += 1;
        }
        let mut owns_v_edge = 0;
        let mut i = 0;
        while i < NUM_BASE_FACES {
            let v_2 = BASE_INDICES[i][2] as usize;
            let v_0 = BASE_INDICES[i][0] as usize;
            let edge = if v_0 < v_2 {
                v_0 * NUM_BASE_VERTS + v_2
            } else {
                v_2 * NUM_BASE_VERTS + v_0
            };
            if !edge_has_owner[edge] {
                edge_has_owner[edge] = true;
                num_owned_edges += 1;
                owns_v_edge |= 1 << i;
            }
            i += 1;
        }

        // Verify that every edge has an owner
        assert!(
            num_owned_edges == NUM_BASE_EDGES,
            "not all edges have an owner"
        );
        owns_v_edge
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
