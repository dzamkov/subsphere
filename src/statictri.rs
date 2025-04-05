//! Contains types related to [`StaticTriSphere`].
use crate::subtri::{BaseRegion, BaseRegionType, BaseSphere, BaseSphereInternal};

/// A [`BaseSphere`] defined explicitly using inline arrays.
///
/// `F` is the number of faces, and `V` is the number of vertices. Since the size of the
/// structure is given, this requires no dynamic memory allocation. This also allows it to
/// constructed in a `const` context.
#[derive(Clone, Copy, Debug)]
pub struct StaticTriSphere<const F: usize, const V: usize> {
    verts: [[f64; 3]; V],
    indices: [[u8; 3]; F],
    edge_angle: f64,
    edge_cos_angle: f64,
    adjacent: [[u8; 3]; F],
    vert_owner: [u8; V],
    owns_vert_1: u64,
    owns_edge_2: u64,
}

impl<const F: usize, const V: usize> StaticTriSphere<F, V> {
    /// Creates a new [`StaticTriSphere`] with the given vertex and face data.
    ///
    /// This will validate that the data describes a valid [`BaseSphere`]. If it doesn't, this will
    /// panic.
    pub const fn new(
        verts: [[f64; 3]; V],
        indices: [[u8; 3]; F],
        edge_angle: f64,
        edge_cos_angle: f64,
    ) -> Self {
        assert!(F < 64, "too many faces");
        assert!(V < 64, "too many vertices");

        // Build adjacent mapping for potential edges
        let mut s_adjacent: [[u8; V]; V] = [[u8::MAX; V]; V];
        let mut i = 0;
        while i < F {
            let v_0 = indices[i][0] as usize;
            let v_1 = indices[i][1] as usize;
            let v_2 = indices[i][2] as usize;
            assert!(s_adjacent[v_0][v_1] == u8::MAX, "duplicate edge detected");
            assert!(s_adjacent[v_1][v_2] == u8::MAX, "duplicate edge detected");
            assert!(s_adjacent[v_2][v_0] == u8::MAX, "duplicate edge detected");
            s_adjacent[v_0][v_1] = (i as u8) << 2;
            s_adjacent[v_1][v_2] = ((i as u8) << 2) | 1;
            s_adjacent[v_2][v_0] = ((i as u8) << 2) | 2;
            i += 1;
        }

        // Convert to adjacency table for faces
        let mut adjacent: [[u8; 3]; F] = [[u8::MAX; 3]; F];
        i = 0;
        while i < F {
            let v_0 = indices[i][0] as usize;
            let v_1 = indices[i][1] as usize;
            let v_2 = indices[i][2] as usize;
            assert!(s_adjacent[v_1][v_0] != u8::MAX, "hanging edge detected");
            assert!(s_adjacent[v_2][v_1] != u8::MAX, "hanging edge detected");
            assert!(s_adjacent[v_0][v_2] != u8::MAX, "hanging edge detected");
            adjacent[i][0] = s_adjacent[v_1][v_0];
            adjacent[i][1] = s_adjacent[v_2][v_1];
            adjacent[i][2] = s_adjacent[v_0][v_2];
            i += 1;
        }

        // Assign ownership of each vertex to the first face that contains it as it's 0 vertex
        let mut vert_owner: [u8; V] = [u8::MAX; V];
        vert_owner[0] = 0;
        let mut owns_vert_1 = 0;
        let mut i = 0;
        while i < F {
            let v_1 = indices[i][1] as usize;
            if vert_owner[v_1] == u8::MAX {
                vert_owner[v_1] = i as u8;
                owns_vert_1 |= 1 << i;
            }
            i += 1;
        }

        // Verify that every vertex has an owner
        let mut j = 0;
        while j < V {
            assert!(vert_owner[j] != u8::MAX, "vertex does not have an owner");
            j += 1;
        }

        // Assign ownership of each edge. First by edge 0, then optionally by edge 2.
        let mut edge_has_owner: [[bool; V]; V] = [[false; V]; V];
        let mut num_owned_edges = 0;
        let mut i = 0;
        while i < F {
            let v_0 = indices[i][0] as usize;
            let v_1 = indices[i][1] as usize;
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
        while i < F {
            let v_2 = indices[i][2] as usize;
            let v_0 = indices[i][0] as usize;
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
        assert!(num_owned_edges == F * 3 / 2, "not all edges have an owner");

        // Finalize sphere
        Self {
            verts,
            indices,
            edge_angle,
            edge_cos_angle,
            adjacent,
            vert_owner,
            owns_vert_1,
            owns_edge_2,
        }
    }
}

/// A wrapper over a [`StaticTriSphere`] reference which implements [`Sphere`](crate::Sphere)
/// and [`BaseSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Deref<const F: usize, const V: usize, S: AsRef<StaticTriSphere<F, V>>>(pub S);

impl<const F: usize, const V: usize, S: Eq + Clone + AsRef<StaticTriSphere<F, V>>> crate::Sphere
    for Deref<F, V, S>
{
    type Face = Face<F, V, S>;
    type Vertex = Vertex<F, V, S>;
    type HalfEdge = HalfEdge<F, V, S>;

    fn num_faces(&self) -> usize {
        F
    }

    fn faces(&self) -> impl Iterator<Item = Self::Face> {
        (0..F).map(move |i| Face {
            sphere: self.0.clone(),
            index: i as u8,
        })
    }

    fn num_vertices(&self) -> usize {
        V
    }

    fn vertices(&self) -> impl Iterator<Item = Self::Vertex> {
        (0..V).map(move |i| Vertex {
            sphere: self.0.clone(),
            index: i as u8,
        })
    }
}

impl<const F: usize, const V: usize, S: Eq + Clone + AsRef<StaticTriSphere<F, V>>> BaseSphere
    for Deref<F, V, S>
{
}

impl<const F: usize, const V: usize, S: Eq + Clone + AsRef<StaticTriSphere<F, V>>>
    BaseSphereInternal for Deref<F, V, S>
{
    type Region = Region<F, V, S>;

    fn edge_angle(&self) -> f64 {
        self.0.as_ref().edge_angle
    }

    fn edge_cos_angle(&self) -> f64 {
        self.0.as_ref().edge_cos_angle
    }

    fn face_owns_vertex_1(&self, face: Self::Face) -> bool {
        let sphere = self.0.as_ref();
        (sphere.owns_vert_1 >> face.index) & 1 == 1
    }

    fn face_owns_edge_2(&self, face: Self::Face) -> bool {
        let sphere = self.0.as_ref();
        (sphere.owns_edge_2 >> face.index) & 1 == 1
    }

    fn vertex_owner(&self, vertex: Self::Vertex) -> Self::Face {
        let sphere = self.0.as_ref();
        Face {
            sphere: self.0.clone(),
            index: sphere.vert_owner[vertex.index as usize],
        }
    }

    fn num_owned_vertices_before(&self, face: Self::Face) -> usize {
        let sphere = self.0.as_ref();
        let before_mask = (1 << face.index) - 1;
        (sphere.owns_vert_1 & before_mask).count_ones() as usize
    }

    fn num_owned_edges_before(&self, face: Self::Face) -> usize {
        let sphere = self.0.as_ref();
        let before_mask = (1 << face.index) - 1;
        face.index as usize + (sphere.owns_edge_2 & before_mask).count_ones() as usize
    }

    fn first_region(&self) -> Self::Region {
        Region {
            sphere: self.0.clone(),
            data: 0,
        }
    }

    fn next_region(&self, mut region: Self::Region) -> Option<Self::Region> {
        if region.ty().is_edge() {
            region.data += 1;
            if (region.data >> 2) >= F as u8 {
                return None;
            }
        } else if self.face_owns_edge_2(region.owner()) {
            region.data += 2;
        } else {
            region.data += 3;
            if (region.data >> 2) >= F as u8 {
                return None;
            }
        }
        Some(region)
    }
}

/// A face of a [`StaticTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Face<const F: usize, const V: usize, S: AsRef<StaticTriSphere<F, V>>> {
    sphere: S,
    index: u8,
}

impl<const F: usize, const V: usize, S: Eq + Clone + AsRef<StaticTriSphere<F, V>>> crate::Face
    for Face<F, V, S>
{
    type Vertex = Vertex<F, V, S>;
    type HalfEdge = HalfEdge<F, V, S>;

    fn index(&self) -> usize {
        self.index as usize
    }

    fn num_sides(&self) -> usize {
        3
    }

    fn side(&self, index: usize) -> HalfEdge<F, V, S> {
        HalfEdge {
            sphere: self.sphere.clone(),
            data: (self.index << 2) | (index as u8),
        }
    }
}

/// A vertex of a [`StaticTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Vertex<const F: usize, const V: usize, S: AsRef<StaticTriSphere<F, V>>> {
    sphere: S,
    index: u8,
}

impl<const F: usize, const V: usize, S: Eq + Clone + AsRef<StaticTriSphere<F, V>>> crate::Vertex
    for Vertex<F, V, S>
{
    type Face = Face<F, V, S>;
    type HalfEdge = HalfEdge<F, V, S>;

    fn index(&self) -> usize {
        self.index as usize
    }

    fn pos(&self) -> [f64; 3] {
        let sphere = self.sphere.as_ref();
        unsafe { *sphere.verts.get_unchecked(self.index as usize) }
    }

    fn degree(&self) -> usize {
        todo!()
    }

    fn first_outgoing(&self) -> Self::HalfEdge {
        todo!()
    }
}

/// A half-edge of a [`StaticTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct HalfEdge<const F: usize, const V: usize, S: AsRef<StaticTriSphere<F, V>>> {
    sphere: S,
    data: u8,
}

impl<const F: usize, const V: usize, S: Eq + Clone + AsRef<StaticTriSphere<F, V>>> crate::HalfEdge
    for HalfEdge<F, V, S>
{
    type Face = Face<F, V, S>;
    type Vertex = Vertex<F, V, S>;

    fn index(&self) -> usize {
        (self.data & 3) as usize
    }

    fn inside(&self) -> Self::Face {
        Face {
            sphere: self.sphere.clone(),
            index: self.data >> 2,
        }
    }

    fn start(&self) -> Self::Vertex {
        let sphere = self.sphere.as_ref();
        let face_index = self.data >> 2;
        let side_index = self.data & 3;
        Vertex {
            sphere: self.sphere.clone(),
            index: unsafe {
                *sphere
                    .indices
                    .get_unchecked(face_index as usize)
                    .get_unchecked(side_index as usize)
            },
        }
    }

    fn complement(&self) -> Self {
        let sphere = self.sphere.as_ref();
        let face_index = self.data >> 2;
        let side_index = self.data & 3;
        HalfEdge {
            sphere: self.sphere.clone(),
            data: unsafe {
                *sphere
                    .adjacent
                    .get_unchecked(face_index as usize)
                    .get_unchecked(side_index as usize)
            },
        }
    }

    fn prev(&self) -> Self {
        todo!()
    }

    fn next(&self) -> Self {
        todo!()
    }
}

/// A [`BaseRegion`] for a [`StaticTriSphere`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Region<const F: usize, const V: usize, S: AsRef<StaticTriSphere<F, V>>> {
    sphere: S,
    data: u8,
}

impl<const F: usize, const V: usize, S: Eq + Clone + AsRef<StaticTriSphere<F, V>>> BaseRegion
    for Region<F, V, S>
{
    type Face = Face<F, V, S>;

    fn new(face: Face<F, V, S>, ty: BaseRegionType) -> Self {
        Region {
            sphere: face.sphere,
            data: (face.index << 2) | (ty as u8),
        }
    }

    fn owner(&self) -> Face<F, V, S> {
        Face {
            sphere: self.sphere.clone(),
            index: self.data >> 2,
        }
    }

    fn ty(&self) -> BaseRegionType {
        unsafe { std::mem::transmute(self.data & 3) }
    }
}
