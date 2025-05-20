//! Contains types related to [`BaseTriSphere`].
use crate::math::vec;

/// A tessellation of the unit sphere constructed by projecting a triangular platonic solid
/// onto it.
#[repr(u8)]
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq, Hash)]
pub enum BaseTriSphere {
    /// A tessellation of the unit sphere constructed by projecting an icosahedron onto it.
    #[default]
    Icosa = 0,

    /// A tessellation of the unit sphere constructed by projecting an octahedron onto it.
    Octa = 1,

    /// A tessellation of the unit sphere constructed by projecting a tetrahedron onto it.
    Tetra = 2,
}

impl BaseTriSphere {
    #[inline]
    pub(crate) const fn from_u8(value: u8) -> Self {
        [ Self::Icosa, Self::Octa, Self::Tetra ][value as usize]
    }
    
    /// The degree of the vertices in this base shape.
    pub const fn vertex_degree(self) -> usize {
        self.lookup::<5, 4, 3>() as usize
    }

    /// The length of any edge on this base shape, or equivalently, the angle between any
    /// two adjacent vertices.
    pub const fn edge_length(self) -> f64 {
        [1.1071487177940904, std::f64::consts::FRAC_PI_2][self as usize]
    }

    /// The cosine of [`edge_length`](BaseTriSphere::edge_length`), or equivalently, the dot
    /// product between any two adjacent vertices.
    pub const fn edge_cos_length(self) -> f64 {
        [C_1, 0.0][self as usize]
    }

    /// The internal representation of the first face of this base shape.
    pub(crate) const fn first_face_inner(self) -> u8 {
        self.lookup::<0, 20, 28>()
    }

    /// One more than the internal representation of the last face of this base shape.
    pub(crate) const fn last_face_inner(self) -> u8 {
        self.lookup::<20, 28, 32>()
    }

    /// The internal representation of the first vertex of this base shape.
    pub(crate) const fn first_vertex_inner(self) -> u8 {
        self.lookup::<0, 12, 18>()
    }

    /// One more than the internal representation of the last vertex of this base shape.
    pub(crate) const fn last_vertex_inner(self) -> u8 {
        self.lookup::<12, 18, 22>()
    }

    /// Returns the constant value corresponding to this base shape.
    const fn lookup<const ICOSA: u8, const OCTA: u8, const TETRA: u8>(self) -> u8 {
        ((((TETRA as u32) << 16) | ((OCTA as u32) << 8) | ICOSA as u32) >> (self as u32 * 8)) as u8
    }

    /// Determines which face contains the given point on the unit sphere.
    pub const fn face_at(self, point: [f64; 3]) -> Face {
        match self {
            BaseTriSphere::Icosa => {
                let index = icosa_point_index(point);
                ICOSA_FACE_AT[index as usize]
            }
            BaseTriSphere::Octa => {
                let index = octa_point_index(point);
                Face((OCTA_FACE_AT >> (index * 8)) as u8)
            }
            BaseTriSphere::Tetra => todo!(),
        }
    }

    /// The number of vertices on the sphere.
    pub const fn num_vertices(self) -> usize {
        self.lookup::<12, 6, 4>() as usize
    }

    /// Gets the [`Vertex`] with the specified index.
    pub const fn vertex(self, index: usize) -> Vertex {
        assert!(index < self.num_vertices(), "index out of bounds");
        Vertex(self.first_vertex_inner() + index as u8)
    }
}

impl crate::Sphere for BaseTriSphere {
    type Face = Face;
    type Vertex = Vertex;
    type HalfEdge = HalfEdge;

    fn num_faces(&self) -> usize {
        self.lookup::<20, 8, 4>() as usize
    }

    fn face(&self, index: usize) -> Face {
        assert!(index < self.num_faces(), "index out of bounds");
        Face(self.first_face_inner() + index as u8)
    }

    fn faces(&self) -> impl Iterator<Item = Face> {
        (self.first_face_inner()..self.last_face_inner()).map(Face)
    }

    fn face_at(&self, point: [f64; 3]) -> Face {
        (*self).face_at(point)
    }

    fn num_vertices(&self) -> usize {
        (*self).num_vertices()
    }

    fn vertex(&self, index: usize) -> Vertex {
        (*self).vertex(index)
    }

    fn vertices(&self) -> impl Iterator<Item = Vertex> {
        (self.first_vertex_inner()..self.last_vertex_inner()).map(Vertex)
    }
}

/// A face of a [`BaseTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Face(pub(crate) u8);

impl Face {
    /// Gets the [`BaseTriSphere`] this face belongs to.
    pub const fn sphere(self) -> BaseTriSphere {
        BaseTriSphere::from_u8(self.0.saturating_sub(12) / 8)
    }

    /// Indicates whether this face [owns](OwnershipInfo) its second vertex.
    pub(crate) fn owns_vertex_1(self) -> bool {
        OWNERSHIP.owns_vert_1 & (1 << self.0) != 0
    }

    /// Indicates whether this face [owns](OwnershipInfo) its third edge.
    pub(crate) fn owns_edge_2(self) -> bool {
        OWNERSHIP.owns_edge_2 & (1 << self.0) != 0
    }

    /// Gets the number of vertices that are owned by the faces preceding this face in
    /// iteration order, not counting the first vertex of the shape.
    pub(crate) fn num_owned_vertices_before(self) -> usize {
        let start = self.sphere().first_face_inner();
        let before_mask = (1u32 << self.0) - 1;
        ((OWNERSHIP.owns_vert_1 & before_mask) >> start).count_ones() as usize
    }

    /// Gets the number of edges that are owned by the faces preceding this face in iteration
    /// order.
    pub(crate) fn num_owned_edges_before(self) -> usize {
        let start = self.sphere().first_face_inner();
        let index = self.0 - start;
        let before_mask = (1u32 << self.0) - 1;
        index as usize + ((OWNERSHIP.owns_edge_2 & before_mask) >> start).count_ones() as usize
    }

    /// Gets the [`HalfEdge`] which has the given [`index`](HalfEdge::side_index) and this face as
    /// its [`inside`](HalfEdge::inside).
    pub const fn side(self, index: usize) -> HalfEdge {
        assert!(index < 3, "index out of bounds");
        HalfEdge((self.0 << 2) | index as u8)
    }

    /// Gets the point at the center of this face.
    pub const fn center(self) -> [f64; 3] {
        let v_0 = self.side(0).start().pos();
        let v_1 = self.side(1).start().pos();
        let v_2 = self.side(2).start().pos();
        let mul = [0.4194695241216063, 0.5773502691896257][self.sphere() as usize]; // TODO
        vec::mul(vec::add(vec::add(v_0, v_1), v_2), mul)
    }
}

#[test]
fn test_center_face_at() {
    use crate::Sphere;
    // TODO: Extend to other spheres
    for sphere in [BaseTriSphere::Icosa, BaseTriSphere::Octa] {
        for face in sphere.faces() {
            let center = face.center();
            assert!((vec::dot(center, center) - 1.0).abs() < 1.0e-12);
            assert_eq!(sphere.face_at(center), face);
        }
    }
}

impl crate::Face for Face {
    type Vertex = Vertex;
    type HalfEdge = HalfEdge;

    fn index(&self) -> usize {
        (self.0 - self.sphere().first_face_inner()) as usize
    }

    fn area(&self) -> f64 {
        // TODO: Replace with lookup table
        let v_0 = self.side(0).start().pos();
        let v_1 = self.side(1).start().pos();
        let v_2 = self.side(2).start().pos();
        crate::util::tri_area([v_0, v_1, v_2])
    }

    fn num_sides(&self) -> usize {
        3
    }

    fn side(&self, index: usize) -> HalfEdge {
        (*self).side(index)
    }
}

/// A vertex of a [`BaseTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Vertex(u8);

impl Vertex {
    /// Gets the [`BaseTriSphere`] this vertex belongs to.
    pub const fn sphere(self) -> BaseTriSphere {
        if self.0 < 12 {
            BaseTriSphere::Icosa
        } else if self.0 < 18 {
            BaseTriSphere::Octa
        } else {
            BaseTriSphere::Tetra
        }
    }

    /// Gets the face which [owns](OwnershipInfo) this vertex.
    pub(crate) const fn owner(self) -> Face {
        OWNERSHIP.vert_owner[self.0 as usize]
    }

    /// The position of this vertex.
    pub const fn pos(self) -> [f64; 3] {
        VERTS[self.0 as usize]
    }
}

impl crate::Vertex for Vertex {
    type Face = Face;
    type HalfEdge = HalfEdge;

    fn index(&self) -> usize {
        (self.0 - self.sphere().first_vertex_inner()) as usize
    }

    fn pos(&self) -> [f64; 3] {
        (*self).pos()
    }

    fn degree(&self) -> usize {
        self.sphere().vertex_degree()
    }

    fn outgoing(&self, index: usize) -> HalfEdge {
        assert!(index < self.degree(), "index out of bounds");
        OUTGOING[self.0 as usize][index]
    }
}

/// A half-edge of a [`BaseTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfEdge(u8);

impl HalfEdge {
    /// Gets the [`BaseTriSphere`] this face belongs to.
    pub const fn sphere(self) -> BaseTriSphere {
        self.inside().sphere()
    }

    /// The index of this half-edge within the [`sides`](Face::sides) list of its
    /// [`inside`](HalfEdge::inside).
    pub const fn side_index(self) -> usize {
        (self.0 & 0b11) as usize
    }

    /// Gets the [`Face`] whose interior boundary contains this half-edge.
    pub const fn inside(self) -> Face {
        Face(self.0 >> 2)
    }

    /// Gets the [`Vertex`] at the "start" of this half-edge.
    pub const fn start(self) -> Vertex {
        Vertex(INDICES[self.inside().0 as usize][self.side_index()])
    }

    /// Gets the complementary half-edge on the opposite side of the edge.
    ///
    /// The returned half-edge will go in the opposite direction along the same edge.
    pub const fn twin(self) -> Self {
        TWIN[self.inside().0 as usize][self.side_index()]
    }

    /// Gets the half-edge which shares the [`inside`](HalfEdge::inside) face of this half-edge and
    /// precedes it in counter-clockwise order around the face.
    pub const fn prev(self) -> Self {
        HalfEdge(if self.0 & 0b11 == 0 {
            self.0 + 2
        } else {
            self.0 - 1
        })
    }

    /// Gets the half-edge which shares the [`inside`](HalfEdge::inside) face of this half-edge and
    /// follows it in counter-clockwise order around the face.
    pub const fn next(self) -> Self {
        HalfEdge(if self.0 & 0b11 == 2 {
            self.0 - 2
        } else {
            self.0 + 1
        })
    }
}

impl crate::HalfEdge for HalfEdge {
    type Face = Face;
    type Vertex = Vertex;

    fn side_index(&self) -> usize {
        (*self).side_index()
    }

    fn length(&self) -> f64 {
        self.sphere().edge_length()
    }

    fn angle(&self) -> f64 {
        // TODO: Replace with lookup table
        let v_a = self.prev().start().pos();
        let v_b = self.start().pos();
        let v_c = self.next().start().pos();
        crate::util::angle(v_a, v_b, v_c)
    }

    fn inside(&self) -> Face {
        (*self).inside()
    }

    fn start(&self) -> Vertex {
        (*self).start()
    }

    fn twin(&self) -> Self {
        (*self).twin()
    }

    fn prev(&self) -> Self {
        (*self).prev()
    }

    fn next(&self) -> Self {
        (*self).next()
    }
}

/// Table used to implement [`crate::Vertex::outgoing`].
static OUTGOING: [[HalfEdge; 5]; NUM_VERTS] = const {
    let mut res: [[HalfEdge; 5]; NUM_VERTS] = [[HalfEdge(u8::MAX); 5]; NUM_VERTS];
    let mut j = 0;
    while j < NUM_VERTS {
        let vert = Vertex(j as u8);
        let owner = vert.owner();

        // Determine a unique first outgoing edge for each vertex
        let first_outgoing = {
            let mut k = 0;
            loop {
                let edge = owner.side(k);
                if edge.start().0 == vert.0 {
                    break edge;
                }
                k += 1;
                if k == 3 {
                    panic!("owner does not contain vertex");
                }
            }
        };

        // Build outgoing edge list
        let mut k = 0;
        let mut outgoing = first_outgoing;
        loop {
            res[j][k] = outgoing;
            k += 1;
            outgoing = outgoing.prev().twin();
            assert!(
                outgoing.start().0 == vert.0,
                "outgoing edge does not start at vertex"
            );
            if outgoing.0 == first_outgoing.0 {
                break;
            }
        }
        assert!(
            k == Vertex(j as u8).sphere().vertex_degree(),
            "degree mismatch"
        );
        j += 1;
    }
    res
};

/// Provides information about the "ownership" relationships between faces, edges, and vertices.
///
/// The rules for ownership are as follows:
///  * Every vertex and edge is owned by exactly one face.
///  * Faces must own their first edge.
///  * Faces may own their third edge.
///  * The first face in a [`BaseTriSphere`] must own its first vertex, which must be the first
///    vertex in the [`BaseTriSphere`].
///  * Faces may own their second vertex.
struct OwnershipInfo {
    vert_owner: [Face; NUM_VERTS],
    owns_vert_1: u32,
    owns_edge_2: u32,
}

/// Provides information about the "ownership" relationships between faces, edges, and vertices.
const OWNERSHIP: OwnershipInfo = const {
    // Assign ownership of first vertex in each shape
    let mut vert_owner: [Face; NUM_VERTS] = [Face(u8::MAX); NUM_VERTS];
    vert_owner[0] = Face(0);
    vert_owner[12] = Face(20);

    // Assign ownership of each vertex to the first face that contains it as it's 1 vertex
    let mut owns_vert_1 = 0;
    let mut i = 0;
    while i < NUM_FACES {
        let v_1 = INDICES[i][1] as usize;
        if vert_owner[v_1].0 == u8::MAX {
            vert_owner[v_1] = Face(i as u8);
            owns_vert_1 |= 1 << i;
        }
        i += 1;
    }

    // Verify that every vertex has an owner
    let mut j = 0;
    while j < NUM_VERTS {
        assert!(vert_owner[j].0 != u8::MAX, "vertex does not have an owner");
        j += 1;
    }

    // Assign ownership of each edge. First by edge 0, then optionally by edge 2.
    let mut edge_has_owner: [[bool; NUM_VERTS]; NUM_VERTS] = [[false; NUM_VERTS]; NUM_VERTS];
    let mut num_owned_edges = 0;
    let mut i = 0;
    while i < NUM_FACES {
        let v_0 = INDICES[i][0] as usize;
        let v_1 = INDICES[i][1] as usize;
        let edge_has_owner = if v_0 < v_1 {
            &mut edge_has_owner[v_0][v_1]
        } else {
            &mut edge_has_owner[v_1][v_0]
        };
        assert!(!*edge_has_owner, "edge already has an owner");
        *edge_has_owner = true;
        num_owned_edges += 1;
        i += 1;
    }
    let mut owns_edge_2 = 0;
    let mut i = 0;
    while i < NUM_FACES {
        let v_2 = INDICES[i][2] as usize;
        let v_0 = INDICES[i][0] as usize;
        let edge_has_owner = if v_0 < v_2 {
            &mut edge_has_owner[v_0][v_2]
        } else {
            &mut edge_has_owner[v_2][v_0]
        };
        if !*edge_has_owner {
            *edge_has_owner = true;
            num_owned_edges += 1;
            owns_edge_2 |= 1 << i;
        }
        i += 1;
    }

    // Verify that every edge has an owner
    assert!(
        num_owned_edges == NUM_FACES * 3 / 2,
        "not all edges have an owner"
    );

    // Finalize ownership info
    OwnershipInfo {
        vert_owner,
        owns_vert_1,
        owns_edge_2,
    }
};

/// Table used to implement [`crate::HalfEdge::twin`].
const TWIN: [[HalfEdge; 3]; NUM_FACES] = const {
    // Build adjacent mapping for potential edges
    let mut adjacent: [[u8; NUM_VERTS]; NUM_VERTS] = [[u8::MAX; NUM_VERTS]; NUM_VERTS];
    let mut i = 0;
    while i < NUM_FACES {
        let v_0 = INDICES[i][0] as usize;
        let v_1 = INDICES[i][1] as usize;
        let v_2 = INDICES[i][2] as usize;
        assert!(adjacent[v_0][v_1] == u8::MAX, "duplicate edge detected");
        assert!(adjacent[v_1][v_2] == u8::MAX, "duplicate edge detected");
        assert!(adjacent[v_2][v_0] == u8::MAX, "duplicate edge detected");
        adjacent[v_0][v_1] = (i as u8) << 2;
        adjacent[v_1][v_2] = ((i as u8) << 2) | 1;
        adjacent[v_2][v_0] = ((i as u8) << 2) | 2;
        i += 1;
    }

    // Convert to adjacency table for faces
    let mut res: [[HalfEdge; 3]; NUM_FACES] = [[HalfEdge(0); 3]; NUM_FACES];
    i = 0;
    while i < NUM_FACES {
        let v_0 = INDICES[i][0] as usize;
        let v_1 = INDICES[i][1] as usize;
        let v_2 = INDICES[i][2] as usize;
        assert!(adjacent[v_1][v_0] != u8::MAX, "hanging edge detected");
        assert!(adjacent[v_2][v_1] != u8::MAX, "hanging edge detected");
        assert!(adjacent[v_0][v_2] != u8::MAX, "hanging edge detected");
        res[i][0] = HalfEdge(adjacent[v_1][v_0]);
        res[i][1] = HalfEdge(adjacent[v_2][v_1]);
        res[i][2] = HalfEdge(adjacent[v_0][v_2]);
        i += 1;
    }
    res
};

/// Given a point on the unit sphere, gets an index which can be used to identify which
/// icosahedron [`Face`] contains it.
///
/// The indexing scheme is arbitrary, but points on different [`Face`]s must have different
/// indices.
const fn icosa_point_index(point: [f64; 3]) -> u8 {
    // The index is constructed by testing the dot product of `point` with 5 distinct (and
    // non-antipodal) vertices. For each test, we get one of 3 results. This gives us 243
    // "possible" indices, with each face corresponding to exactly one of these.
    let mut res = 0;
    let mut i = 0;
    while i < 5 {
        let dot = vec::dot(point, BaseTriSphere::Icosa.vertex(i).pos());
        let comp = (dot.is_sign_positive() as u8 + 1) * ((dot.abs() > C_1) as u8);
        res = res * 3 + comp;
        i += 1;
    }
    res
}

/// A lookup table which identifies which [`Face`] corresponds to a particular result from
/// [`icosa_point_index`].
///
/// For indices that don't correspond exactly to a face, this will provide an arbitrary
/// nearby [`Face`].
const ICOSA_FACE_AT: [Face; 243] = const {
    let mut res = [Face(u8::MAX); 243];

    // Assign each face to its proper index
    let sphere = BaseTriSphere::Icosa;
    let mut i = sphere.first_face_inner();
    while i < sphere.last_face_inner() {
        let face = Face(i);
        let j = icosa_point_index(face.center());
        assert!(res[j as usize].0 == u8::MAX, "index already assigned");
        res[j as usize] = face;
        i += 1;
    }

    // Fill remaining indices by iteratively copying the contents of a "nearby" index which is
    // already filled.
    let mut next = res;
    let mut all_filled = false;
    while !all_filled {
        let mut index = 0;
        all_filled = true;
        while index < 243 {
            if res[index].0 == u8::MAX {
                all_filled = false;

                // Try to find a component we can change slightly to get a filled index.
                let mut mul = 1;
                while mul <= 83 {
                    let comp = (index / mul) % 3;
                    if comp == 0 {
                        if res[index + mul].0 != u8::MAX {
                            next[index] = res[index + mul];
                            break;
                        }
                        if res[index + 2 * mul].0 != u8::MAX {
                            next[index] = res[index + 2 * mul];
                            break;
                        }
                    } else if res[index - comp * mul].0 != u8::MAX {
                        next[index] = res[index - comp * mul];
                        break;
                    }
                    mul *= 3;
                }
            }
            index += 1;
        }
        res = next;
    }
    res
};

/// Given a point on the unit sphere, gets an index which can be used to identify which
/// octahedron [`Face`] contains it.
///
/// The indexing scheme is arbitrary, but points on different [`Face`]s must have different
/// indices.
const fn octa_point_index([x, y, z]: [f64; 3]) -> u8 {
    (((x >= 0.0) as u8) << 2) | (((y >= 0.0) as u8) << 1) | ((z >= 0.0) as u8)
}

/// A compact lookup table which identifies which [`Face`] corresponds to a particular result
/// from [`octa_point_index`].
const OCTA_FACE_AT: u64 = const {
    let sphere = BaseTriSphere::Octa;
    let mut i = sphere.first_face_inner();
    let mut res = 0;
    while i < sphere.last_face_inner() {
        let face = Face(i);
        let j = octa_point_index(face.center());
        res |= (i as u64) << (j * 8);
        i += 1;
    }
    res
};

/// The total number of vertices across all [`BaseTriSphere`]s.
const NUM_VERTS: usize = 12 + 6;

/// The total number of faces across all [`BaseTriSphere`]s.
const NUM_FACES: usize = 20 + 8;

/// The vertex position data for all potential vertices on a [`BaseTriSphere`].
const VERTS: [[f64; 3]; NUM_VERTS] = [
    // Icosahedron top apex
    [0.0, 0.0, 1.0],
    // Icosahedron top ring
    [C_0, 0.0, C_1],
    [C_2, C_3, C_1],
    [-C_4, C_5, C_1],
    [-C_4, -C_5, C_1],
    [C_2, -C_3, C_1],
    // Icosahedron bottom ring
    [C_4, -C_5, -C_1],
    [C_4, C_5, -C_1],
    [-C_2, C_3, -C_1],
    [-C_0, 0.0, -C_1],
    [-C_2, -C_3, -C_1],
    // Icosahedron bottom apex
    [0.0, 0.0, -1.0],
    // Octahedron
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [-1.0, 0.0, 0.0],
    [0.0, -1.0, 0.0],
    [0.0, 0.0, -1.0],
];

/// The face index data for all potential faces on a [`BaseTriSphere`].
const INDICES: [[u8; 3]; NUM_FACES] = [
    // Icosahedron top cap
    [0, 1, 2],
    [0, 2, 3],
    [0, 3, 4],
    [0, 4, 5],
    [0, 5, 1],
    // Icosahedron central ring
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
    // Icosahedron bottom cap
    [10, 11, 6],
    [6, 11, 7],
    [7, 11, 8],
    [8, 11, 9],
    [9, 11, 10],
    // Octahedron top cap
    [12, 13, 14],
    [12, 14, 15],
    [12, 15, 16],
    [12, 16, 13],
    // Octahedron bottom cap
    [16, 17, 13],
    [13, 17, 14],
    [14, 17, 15],
    [15, 17, 16],
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
