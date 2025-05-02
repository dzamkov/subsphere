use std::hint::black_box;
use std::num::NonZero;
use subsphere::prelude::*;

fn main() {
    divan::main();
}

fn tessellation_vertex_benchmark<Tessellation: subsphere::Sphere>(sphere: Tessellation) {
    black_box(
        sphere.vertices()
            .map(|v| v.pos())
            .collect::<Vec<_>>()
    );
}

fn tessellation_face_benchmark<Tessellation: subsphere::Sphere>(sphere: Tessellation) {
    let mut indices = Vec::new();
    for face in sphere.faces() {
        indices.extend(face.vertices().map(|v| v.index()));
    }
    black_box(indices);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn trisphere_icosa_fuller_b_0_vertex(b: u32) {
    let sphere = subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Fuller,
        NonZero::new(b).unwrap(),
        0,
    );
    
    tessellation_vertex_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn trisphere_icosa_gnomonic_b_0_vertex(b: u32) {
    let sphere = subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Gnomonic,
        NonZero::new(b).unwrap(),
        0,
    );
    
    tessellation_vertex_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn trisphere_octo_fuller_b_0_vertex(b: u32) {
    let sphere = subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Octo,
        subsphere::proj::Fuller,
        NonZero::new(b).unwrap(),
        0,
    );
    
    tessellation_vertex_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn trisphere_octo_gnomonic_b_0_vertex(b: u32) {
    let sphere = subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Octo,
        subsphere::proj::Gnomonic,
        NonZero::new(b).unwrap(),
        0,
    );
    
    tessellation_vertex_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn trisphere_icosa_fuller_b_0_face(b: u32) {
    let sphere = subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Fuller,
        NonZero::new(b).unwrap(),
        0,
    );
    
    tessellation_face_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn trisphere_icosa_gnomonic_b_0_face(b: u32) {
    let sphere = subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Gnomonic,
        NonZero::new(b).unwrap(),
        0,
    );
    
    tessellation_face_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn trisphere_octo_fuller_b_0_face(b: u32) {
    let sphere = subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Octo,
        subsphere::proj::Fuller,
        NonZero::new(b).unwrap(),
        0,
    );
    
    tessellation_face_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn trisphere_octo_gnomonic_b_0_face(b: u32) {
    let sphere = subsphere::TriSphere::new(
        subsphere::BaseTriSphere::Octo,
        subsphere::proj::Gnomonic,
        NonZero::new(b).unwrap(),
        0,
    );
    
    tessellation_face_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn hexsphere_icosa_fuller_b_bmod3_vertex(b: u32) {
    let c = b % 3;
    
    let sphere = subsphere::HexSphere::from_kis(
        subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Icosa,
            subsphere::proj::Fuller,
            NonZero::new(b).unwrap(),
            c,
        )
    ).unwrap();
    
    tessellation_vertex_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn hexsphere_icosa_gnomonic_b_bmod3_vertex(b: u32) {
    let c = b % 3;
    
    let sphere = subsphere::HexSphere::from_kis(
        subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Icosa,
            subsphere::proj::Gnomonic,
            NonZero::new(b).unwrap(),
            c,
        )
    ).unwrap();
    
    tessellation_vertex_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn hexsphere_octo_fuller_b_bmod3_vertex(b: u32) {
    let c = b % 3;
    
    let sphere = subsphere::HexSphere::from_kis(
        subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Octo,
            subsphere::proj::Fuller,
            NonZero::new(b).unwrap(),
            c,
        )
    ).unwrap();
    
    tessellation_vertex_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn hexsphere_octo_gnomonic_b_bmod3_vertex(b: u32) {
    let c = b % 3;
    
    let sphere = subsphere::HexSphere::from_kis(
        subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Octo,
            subsphere::proj::Gnomonic,
            NonZero::new(b).unwrap(),
            c,
        )
    ).unwrap();
    
    tessellation_vertex_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn hexsphere_icosa_fuller_b_bmod3_face(b: u32) {
    let c = b % 3;
    
    let sphere = subsphere::HexSphere::from_kis(
        subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Icosa,
            subsphere::proj::Fuller,
            NonZero::new(b).unwrap(),
            c,
        )
    ).unwrap();
    
    tessellation_face_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn hexsphere_icosa_gnomonic_b_bmod3_face(b: u32) {
    let c = b % 3;
    
    let sphere = subsphere::HexSphere::from_kis(
        subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Icosa,
            subsphere::proj::Gnomonic,
            NonZero::new(b).unwrap(),
            c,
        )
    ).unwrap();
    
    tessellation_face_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn hexsphere_octo_fuller_b_bmod3_face(b: u32) {
    let c = b % 3;
    
    let sphere = subsphere::HexSphere::from_kis(
        subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Octo,
            subsphere::proj::Fuller,
            NonZero::new(b).unwrap(),
            c,
        )
    ).unwrap();
    
    tessellation_face_benchmark(sphere);
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn hexsphere_octo_gnomonic_b_bmod3_face(b: u32) {
    let c = b % 3;
    
    let sphere = subsphere::HexSphere::from_kis(
        subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Octo,
            subsphere::proj::Gnomonic,
            NonZero::new(b).unwrap(),
            c,
        )
    ).unwrap();
    
    tessellation_face_benchmark(sphere);
}