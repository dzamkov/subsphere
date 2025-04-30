use std::hint::black_box;
use std::num::NonZero;
use subsphere::prelude::*;

fn main() {
    divan::main();
}

fn general_benchmark<Proj: Eq + Clone + subsphere::proj::BaseTriProjector>(base: subsphere::BaseTriSphere, proj: Proj, b: NonZero<u32>, c: u32) {
    let sphere = subsphere::TriSphere::new(
        base,
        proj,
        b,
        c,
    );
    
    black_box(
        sphere.vertices()
            .map(|v| black_box(v.pos()))
            .collect::<Vec<_>>()
    );
    
    black_box(
        sphere.faces()
            .collect::<Vec<_>>()
    );
}

fn general_benchmark_kis<Proj: Eq + Clone + subsphere::proj::BaseTriProjector>(base: subsphere::BaseTriSphere, proj: Proj, b: NonZero<u32>, c: u32) {
    let sphere = subsphere::HexSphere::from_kis(
        subsphere::TriSphere::new(
            base,
            proj,
            b,
            c,
        )
    ).unwrap();
    
    black_box(
        sphere.vertices()
            .map(|v| black_box(v.pos()))
            .collect::<Vec<_>>()
    );
    
    black_box(
        sphere.faces()
            .collect::<Vec<_>>()
    );
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn icosahedron_fuller(b: u32) {
    general_benchmark(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Fuller,
        NonZero::new(b).unwrap(),
        0,
    );
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn icosahedron_gnomonic(b: u32) {
    general_benchmark(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Gnomonic,
        NonZero::new(b).unwrap(),
        0,
    );
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn octahedron_fuller(b: u32) {
    general_benchmark(
        subsphere::BaseTriSphere::Octo,
        subsphere::proj::Fuller,
        NonZero::new(b).unwrap(),
        0,
    );
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn octahedron_gnomonic(b: u32) {
    general_benchmark(
        subsphere::BaseTriSphere::Octo,
        subsphere::proj::Gnomonic,
        NonZero::new(b).unwrap(),
        0,
    );
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn icosahedron_fuller_kis(b: u32) {
    let c = b % 3;
    general_benchmark_kis(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Fuller,
        NonZero::new(b).unwrap(),
        c,
    );
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn icosahedron_gnomonic_kis(b: u32) {
    let c = b % 3;
    general_benchmark_kis(
        subsphere::BaseTriSphere::Icosa,
        subsphere::proj::Gnomonic,
        NonZero::new(b).unwrap(),
        c,
    );
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn octahedron_fuller_kis(b: u32) {
    let c = b % 3;
    general_benchmark_kis(
        subsphere::BaseTriSphere::Octo,
        subsphere::proj::Fuller,
        NonZero::new(b).unwrap(),
        c,
    );
}

#[divan::bench(args = [1, 2, 4, 8, 16, 32, 64, 128])]
fn octahedron_gnomonic_kis(b: u32) {
    let c = b % 3;
    general_benchmark_kis(
        subsphere::BaseTriSphere::Octo,
        subsphere::proj::Gnomonic,
        NonZero::new(b).unwrap(),
        c,
    );
}