use std::hint::black_box;
use std::num::NonZero;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use subsphere::prelude::*;

criterion_group!(
    trisphere_octa,
    trisphere_octa_fuller_b_0_face,
    trisphere_octa_fuller_b_0_vertex,
    trisphere_octa_gnomonic_b_0_face,
    trisphere_octa_gnomonic_b_0_vertex
);

criterion_group!(
    trisphere_icosa,
    trisphere_icosa_fuller_b_0_face,
    trisphere_icosa_fuller_b_0_vertex,
    trisphere_icosa_gnomonic_b_0_face,
    trisphere_icosa_gnomonic_b_0_vertex,
);

criterion_group!(
    hexsphere_octa,
    hexsphere_octa_fuller_b_bmod3_face,
    hexsphere_octa_fuller_b_bmod3_vertex,
    hexsphere_octa_gnomonic_b_bmod3_face,
    hexsphere_octa_gnomonic_b_bmod3_vertex
);

criterion_group!(
    hexsphere_icosa,
    hexsphere_icosa_fuller_b_bmod3_face,
    hexsphere_icosa_fuller_b_bmod3_vertex,
    hexsphere_icosa_gnomonic_b_bmod3_face,
    hexsphere_icosa_gnomonic_b_bmod3_vertex,
);

criterion_main!(
    trisphere_octa,
    trisphere_icosa,
    hexsphere_icosa,
    hexsphere_octa
);

fn tessellation_vertex_benchmark<Tessellation: subsphere::Sphere>(sphere: Tessellation, n: usize) {
    let mut vertices = Vec::with_capacity(n);
    
    vertices.extend(sphere.vertices()
        .map(|v| v.pos()));
    
    black_box(vertices);
}

fn tessellation_face_benchmark<Tessellation: subsphere::Sphere>(sphere: Tessellation, n: usize) {
    let mut indices = Vec::with_capacity(n);
    
    for face in sphere.faces() {
        indices.extend(face.vertices().map(|v| v.index()));
    }
    black_box(indices);
}

fn trisphere_icosa_fuller_b_0_vertex(c: &mut Criterion) {
    let mut group = c.benchmark_group("trisphere_icosa_fuller_b_0_vertex");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let sphere = subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Icosa,
            subsphere::proj::Fuller,
            NonZero::new(b).unwrap(),
            0,
        );
        let n = sphere.vertices().count();
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_vertex_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn trisphere_icosa_gnomonic_b_0_vertex(c: &mut Criterion) {
    let mut group = c.benchmark_group("trisphere_icosa_gnomonic_b_0_vertex");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let sphere = subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Icosa,
            subsphere::proj::Gnomonic,
            NonZero::new(b).unwrap(),
            0,
        );
        let n = sphere.vertices().count();
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_vertex_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn trisphere_octa_fuller_b_0_vertex(c: &mut Criterion) {
    let mut group = c.benchmark_group("trisphere_octa_fuller_b_0_vertex");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let sphere = subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Octa,
            subsphere::proj::Fuller,
            NonZero::new(b).unwrap(),
            0,
        );
        let n = sphere.vertices().count();
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_vertex_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn trisphere_octa_gnomonic_b_0_vertex(c: &mut Criterion) {
    let mut group = c.benchmark_group("trisphere_octa_gnomonic_b_0_vertex");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let sphere = subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Octa,
            subsphere::proj::Gnomonic,
            NonZero::new(b).unwrap(),
            0,
        );
        let n = sphere.vertices().count();
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_vertex_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn trisphere_icosa_fuller_b_0_face(c: &mut Criterion) {
    let mut group = c.benchmark_group("trisphere_icosa_fuller_b_0_face");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let sphere = subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Icosa,
            subsphere::proj::Fuller,
            NonZero::new(b).unwrap(),
            0,
        );
        let n = sphere.num_faces() * (sphere.face(0).vertices().count() + 1);
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_face_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn trisphere_icosa_gnomonic_b_0_face(c: &mut Criterion) {
    let mut group = c.benchmark_group("trisphere_icosa_gnomonic_b_0_face");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let sphere = subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Icosa,
            subsphere::proj::Gnomonic,
            NonZero::new(b).unwrap(),
            0,
        );
        let n = sphere.num_faces() * (sphere.face(0).vertices().count() + 1);
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_face_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn trisphere_octa_fuller_b_0_face(c: &mut Criterion) {
    let mut group = c.benchmark_group("trisphere_octa_fuller_b_0_face");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let sphere = subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Octa,
            subsphere::proj::Fuller,
            NonZero::new(b).unwrap(),
            0,
        );
        let n = sphere.num_faces() * (sphere.face(0).vertices().count() + 1);
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_face_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn trisphere_octa_gnomonic_b_0_face(c: &mut Criterion) {
    let mut group = c.benchmark_group("trisphere_octa_gnomonic_b_0_face");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let sphere = subsphere::TriSphere::new(
            subsphere::BaseTriSphere::Octa,
            subsphere::proj::Gnomonic,
            NonZero::new(b).unwrap(),
            0,
        );
        let n = sphere.num_faces() * (sphere.face(0).vertices().count() + 1);
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_face_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn hexsphere_icosa_fuller_b_bmod3_vertex(c: &mut Criterion) {
    let mut group = c.benchmark_group("hexsphere_icosa_fuller_b_bmod3_vertex");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let c = b % 3;
        
        let sphere = subsphere::HexSphere::from_kis(
            subsphere::TriSphere::new(
                subsphere::BaseTriSphere::Icosa,
                subsphere::proj::Fuller,
                NonZero::new(b).unwrap(),
                c,
            )
        ).unwrap();
        let n = sphere.vertices().count();
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_vertex_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn hexsphere_icosa_gnomonic_b_bmod3_vertex(c: &mut Criterion) {
    let mut group = c.benchmark_group("hexsphere_icosa_gnomonic_b_bmod3_vertex");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let c = b % 3;
        
        let sphere = subsphere::HexSphere::from_kis(
            subsphere::TriSphere::new(
                subsphere::BaseTriSphere::Icosa,
                subsphere::proj::Gnomonic,
                NonZero::new(b).unwrap(),
                c,
            )
        ).unwrap();
        let n = sphere.vertices().count();
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_vertex_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn hexsphere_octa_fuller_b_bmod3_vertex(c: &mut Criterion) {
    let mut group = c.benchmark_group("hexsphere_octa_fuller_b_bmod3_vertex");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let c = b % 3;
        
        let sphere = subsphere::HexSphere::from_kis(
            subsphere::TriSphere::new(
                subsphere::BaseTriSphere::Octa,
                subsphere::proj::Fuller,
                NonZero::new(b).unwrap(),
                c,
            )
        ).unwrap();
        let n = sphere.vertices().count();
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_vertex_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn hexsphere_octa_gnomonic_b_bmod3_vertex(c: &mut Criterion) {
    let mut group = c.benchmark_group("hexsphere_octa_gnomonic_b_bmod3_vertex");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let c = b % 3;
        
        let sphere = subsphere::HexSphere::from_kis(
            subsphere::TriSphere::new(
                subsphere::BaseTriSphere::Octa,
                subsphere::proj::Gnomonic,
                NonZero::new(b).unwrap(),
                c,
            )
        ).unwrap();
        let n = sphere.vertices().count();
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_vertex_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn hexsphere_icosa_fuller_b_bmod3_face(c: &mut Criterion) {
    let mut group = c.benchmark_group("hexsphere_icosa_fuller_b_bmod3_face");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let c = b % 3;
        
        let sphere = subsphere::HexSphere::from_kis(
            subsphere::TriSphere::new(
                subsphere::BaseTriSphere::Icosa,
                subsphere::proj::Fuller,
                NonZero::new(b).unwrap(),
                c,
            )
        ).unwrap();
        let n = sphere.num_faces() * (sphere.face(0).vertices().count() + 1);
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_face_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn hexsphere_icosa_gnomonic_b_bmod3_face(c: &mut Criterion) {
    let mut group = c.benchmark_group("hexsphere_icosa_gnomonic_b_bmod3_face");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let c = b % 3;
        
        let sphere = subsphere::HexSphere::from_kis(
            subsphere::TriSphere::new(
                subsphere::BaseTriSphere::Icosa,
                subsphere::proj::Gnomonic,
                NonZero::new(b).unwrap(),
                c,
            )
        ).unwrap();
        let n = sphere.num_faces() * (sphere.face(0).vertices().count() + 1);
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_face_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn hexsphere_octa_fuller_b_bmod3_face(c: &mut Criterion) {
    let mut group = c.benchmark_group("hexsphere_octa_fuller_b_bmod3_face");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let c = b % 3;
        
        let sphere = subsphere::HexSphere::from_kis(
            subsphere::TriSphere::new(
                subsphere::BaseTriSphere::Octa,
                subsphere::proj::Fuller,
                NonZero::new(b).unwrap(),
                c,
            )
        ).unwrap();
        let n = sphere.num_faces() * (sphere.face(0).vertices().count() + 1);
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_face_benchmark(sphere, n));
        });
    }
    group.finish();
}

fn hexsphere_octa_gnomonic_b_bmod3_face(c: &mut Criterion) {
    let mut group = c.benchmark_group("hexsphere_octa_gnomonic_b_bmod3_face");
    for b in [1, 2, 4, 8, 16, 32, 64, 128] {
        let c = b % 3;
        
        let sphere = subsphere::HexSphere::from_kis(
            subsphere::TriSphere::new(
                subsphere::BaseTriSphere::Octa,
                subsphere::proj::Gnomonic,
                NonZero::new(b).unwrap(),
                c,
            )
        ).unwrap();
        let n = sphere.num_faces() * (sphere.face(0).vertices().count() + 1);
        
        group.bench_with_input(BenchmarkId::from_parameter(b), &b, |bencher, _| {
            bencher.iter(|| tessellation_face_benchmark(sphere, n));
        });
    }
    group.finish();
}
