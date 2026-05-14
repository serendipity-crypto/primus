use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use primus_modulus::{BarrettModulus, MontgomeryModulus, PowOf2Modulus};
use primus_reduce::ops::*;
use rand::distr::{Distribution, Uniform};

type ValueT = u64;

fn bench_reduce_mul(c: &mut Criterion) {
    // Barrett
    let barrett_mod = BarrettModulus::<ValueT>::new((1u64 << 62) - 1);
    let m = barrett_mod.value();

    let distr = Uniform::new(0, m).unwrap();
    let mut rng = rand::rng();
    let data: Vec<(ValueT, ValueT)> = (0..1000)
        .map(|_| (distr.sample(&mut rng), distr.sample(&mut rng)))
        .collect();

    c.bench_function("Barrett_reduce_mul", |b| {
        b.iter(|| {
            for &(a, b_val) in &data {
                black_box(barrett_mod.reduce_mul(black_box(a), black_box(b_val)));
            }
        })
    });

    // Montgomery
    let mont_mod = MontgomeryModulus::<ValueT>::new(m | 1); // ensure odd
    let m = mont_mod.value();
    let distr = Uniform::new(0, m).unwrap();
    let mut rng = rand::rng();
    let data: Vec<(ValueT, ValueT)> = (0..1000)
        .map(|_| (distr.sample(&mut rng), distr.sample(&mut rng)))
        .collect();

    c.bench_function("Montgomery_reduce_mul", |b| {
        b.iter(|| {
            for &(a, b_val) in &data {
                let a_m = mont_mod.to_montgomery(a);
                let b_m = mont_mod.to_montgomery(b_val);
                black_box(mont_mod.reduce_mul(black_box(a_m), black_box(b_m)));
            }
        })
    });

    // PowOf2
    let pow2_mod = PowOf2Modulus::<ValueT>::new(1u64 << 50);
    let m = pow2_mod.value();
    let distr = Uniform::new(0, m).unwrap();
    let mut rng = rand::rng();
    let data: Vec<(ValueT, ValueT)> = (0..1000)
        .map(|_| (distr.sample(&mut rng), distr.sample(&mut rng)))
        .collect();

    c.bench_function("PowOf2_reduce_mul", |b| {
        b.iter(|| {
            for &(a, b_val) in &data {
                black_box(pow2_mod.reduce_mul(black_box(a), black_box(b_val)));
            }
        })
    });
}

criterion_group!(benches, bench_reduce_mul);
criterion_main!(benches);
