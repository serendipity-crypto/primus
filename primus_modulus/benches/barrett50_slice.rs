use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use primus_modulus::{Barrett50Modulus, BarrettModulus};
use primus_reduce::{LazyReduceMulSlice, ReduceDotProduct, ReduceMulAddSlice, ReduceMulSlice};
use rand::distr::{Distribution, Uniform};

type V = u64;

/// FHE-style 50-bit prime: 2^50 − 27. Comfortably inside `[2^48, 2^50)`.
const MODULUS: V = (1u64 << 50) - 27;

fn rand_vec(n: usize, m: V) -> Vec<V> {
    let mut rng = rand::rng();
    let distr = Uniform::new(0, m).unwrap();
    distr.sample_iter(&mut rng).take(n).collect()
}

fn bench_mul_to(c: &mut Criterion) {
    let m50 = Barrett50Modulus::new(MODULUS);
    let mref = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("barrett50_mul_slice_to");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let mut out = vec![0; n];

        group.bench_function(format!("barrett64/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                mref.reduce_mul_slice_to(a, b, out);
            });
        });

        group.bench_function(format!("barrett50/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                m50.reduce_mul_slice_to(a, b, out);
            });
        });
    }
    group.finish();
}

fn bench_lazy_mul_to(c: &mut Criterion) {
    let m50 = Barrett50Modulus::new(MODULUS);
    let mref = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("barrett50_lazy_mul_slice_to");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let mut out = vec![0; n];

        group.bench_function(format!("barrett64/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                mref.lazy_reduce_mul_slice_to(a, b, out);
            });
        });

        group.bench_function(format!("barrett50/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                m50.lazy_reduce_mul_slice_to(a, b, out);
            });
        });
    }
    group.finish();
}

fn bench_mul_add_to(c: &mut Criterion) {
    let m50 = Barrett50Modulus::new(MODULUS);
    let mref = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("barrett50_mul_add_slice_to");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let c_in = rand_vec(n, MODULUS);
        let mut out = vec![0; n];

        group.bench_function(format!("barrett64/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                let c = black_box(&c_in);
                mref.reduce_mul_add_slice_to(a, b, c, out);
            });
        });

        group.bench_function(format!("barrett50/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                let c = black_box(&c_in);
                m50.reduce_mul_add_slice_to(a, b, c, out);
            });
        });
    }
    group.finish();
}

fn bench_dot_product(c: &mut Criterion) {
    let m50 = Barrett50Modulus::new(MODULUS);
    let mref = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("barrett50_dot_product");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);

        group.bench_function(format!("barrett64/n={n}"), |bencher| {
            bencher.iter(|| {
                let a = black_box(&a);
                let b = black_box(&b);
                mref.reduce_dot_product(a, b)
            });
        });

        group.bench_function(format!("barrett50/n={n}"), |bencher| {
            bencher.iter(|| {
                let a = black_box(&a);
                let b = black_box(&b);
                m50.reduce_dot_product(a, b)
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_mul_to,
    bench_lazy_mul_to,
    bench_mul_add_to,
    bench_dot_product,
);
criterion_main!(benches);
