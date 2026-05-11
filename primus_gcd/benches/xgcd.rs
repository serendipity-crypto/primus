use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use primus_gcd::Xgcd;
use rand::distr::{Distribution, Uniform};

type ValueT = u64;

fn bench_gcd(c: &mut Criterion) {
    let distr = Uniform::new(0, ValueT::MAX >> 1).unwrap();
    let mut rng = rand::rng();
    // (larger, smaller) so x >= y for xgcd
    let data: Vec<(ValueT, ValueT)> = (0..1000)
        .map(|_| {
            let x = distr.sample(&mut rng);
            let y = distr.sample(&mut rng);
            (x.max(y), x.min(y))
        })
        .collect();

    c.bench_function("gcd", |b| {
        b.iter(|| {
            for &(x, y) in &data {
                black_box(black_box(x).gcd(black_box(y)));
            }
        })
    });

    c.bench_function("is_coprime", |b| {
        b.iter(|| {
            for &(x, y) in &data {
                black_box(black_box(x).is_coprime(black_box(y)));
            }
        })
    });

    c.bench_function("xgcd", |b| {
        b.iter(|| {
            for &(x, y) in &data {
                black_box(ValueT::xgcd(black_box(x), black_box(y)));
            }
        })
    });

    c.bench_function("gcdinv", |b| {
        b.iter(|| {
            // gcdinv requires x < y, so swap the pair
            for &(larger, smaller) in &data {
                black_box(ValueT::gcdinv(black_box(smaller), black_box(larger)));
            }
        })
    });
}

criterion_group!(benches, bench_gcd);
criterion_main!(benches);
