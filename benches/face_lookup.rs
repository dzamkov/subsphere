use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use std::{hint::black_box, num::NonZero};
use subsphere::prelude::*;

fn face_lookup(c: &mut Criterion) {
    let mut group = c.benchmark_group("hexsphere_face_lookup");
    for frequency in [4, 44, 444] {
        let sphere = subsphere::icosphere()
            .subdivide_edge(NonZero::new(frequency).unwrap())
            .with_projector(subsphere::proj::Fuller)
            .truncate();
        let count = sphere.num_faces();
        let mut seed = 1234567_u64;
        let random: Vec<_> = (0..4096)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                (seed % count as u64) as usize
            })
            .collect();
        for (name, indices) in [
            ("first", vec![0]),
            ("middle", vec![count / 2]),
            ("last", vec![count - 1]),
            ("random", random),
        ] {
            for direct in [true, false] {
                let method = if direct { "direct" } else { "iterator" };
                let id = BenchmarkId::new(format!("{method}/{name}"), count);
                group.bench_with_input(id, &indices, |bencher, indices| {
                    let sphere = black_box(sphere);
                    let mut cursor = 0;
                    if direct {
                        bencher.iter(|| {
                            let index = black_box(indices[cursor]);
                            cursor = (cursor + 1) % indices.len();
                            black_box(sphere.face(index))
                        });
                    } else {
                        // The implementation replaced by direct lookup.
                        bencher.iter(|| {
                            let index = black_box(indices[cursor]);
                            cursor = (cursor + 1) % indices.len();
                            black_box(sphere.faces().nth(index).unwrap())
                        });
                    }
                });
            }
        }
    }
    group.finish();
}

criterion_group!(benches, face_lookup);
criterion_main!(benches);
