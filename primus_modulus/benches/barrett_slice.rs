use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use primus_modulus::BarrettModulus;
use primus_reduce::{
    LazyReduceMulAdd, LazyReduceMulAddSlice, ReduceAddAssign, ReduceAddSlice, ReduceMul,
    ReduceMulAdd, ReduceMulAddSlice, ReduceMulSlice,
};
use rand::distr::{Distribution, Uniform};

type V = u64;

// 62-bit NTT-friendly prime, well within Barrett's `bit_count < T::BITS - 1` bound.
const MODULUS: V = 4_611_686_018_427_322_369;

fn rand_vec(n: usize, m: V) -> Vec<V> {
    let mut rng = rand::rng();
    let distr = Uniform::new(0, m).unwrap();
    distr.sample_iter(&mut rng).take(n).collect()
}

fn bench_add(c: &mut Criterion) {
    let m = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("add_slice");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let mut buf = a.clone();

        group.bench_function(format!("per_elem/n={n}"), |bencher| {
            bencher.iter(|| {
                let dst = black_box(&mut buf);
                let src = black_box(&b);
                dst.iter_mut()
                    .zip(src.iter())
                    .for_each(|(d, &s)| m.reduce_add_assign(d, s));
            });
        });

        group.bench_function(format!("slice/n={n}"), |bencher| {
            bencher.iter(|| {
                let dst = black_box(&mut buf);
                let src = black_box(&b);
                m.reduce_add_slice_assign(dst, src);
            });
        });
    }
    group.finish();
}

fn bench_mul(c: &mut Criterion) {
    let m = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("mul_slice_to");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let mut out = vec![0; n];

        group.bench_function(format!("per_elem/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                out.iter_mut()
                    .zip(a.iter().zip(b.iter()))
                    .for_each(|(o, (&a, &b))| *o = m.reduce_mul(a, b));
            });
        });

        group.bench_function(format!("slice/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                m.reduce_mul_slice_to(a, b, out);
            });
        });
    }
    group.finish();
}

fn bench_add_mul_assign(c: &mut Criterion) {
    let m = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("add_mul_slice_assign");
    for &n in &[1024usize, 4096, 16384] {
        let acc_init = rand_vec(n, MODULUS);
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let mut acc = acc_init.clone();

        group.bench_function(format!("per_elem/n={n}"), |bencher| {
            bencher.iter(|| {
                let acc = black_box(&mut acc);
                let a = black_box(&a);
                let b = black_box(&b);
                acc.iter_mut()
                    .zip(a.iter().zip(b.iter()))
                    .for_each(|(z, (&x, &y))| *z = m.reduce_mul_add(x, y, *z));
            });
        });

        group.bench_function(format!("slice/n={n}"), |bencher| {
            bencher.iter(|| {
                let acc = black_box(&mut acc);
                let a = black_box(&a);
                let b = black_box(&b);
                m.reduce_add_mul_slice_assign(acc, a, b);
            });
        });
    }
    group.finish();
}

fn bench_lazy_add_mul_assign(c: &mut Criterion) {
    let m = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("lazy_add_mul_slice_assign");
    for &n in &[1024usize, 4096, 16384] {
        let acc_init = rand_vec(n, MODULUS);
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let mut acc = acc_init.clone();

        group.bench_function(format!("per_elem/n={n}"), |bencher| {
            bencher.iter(|| {
                let acc = black_box(&mut acc);
                let a = black_box(&a);
                let b = black_box(&b);
                acc.iter_mut()
                    .zip(a.iter().zip(b.iter()))
                    .for_each(|(z, (&x, &y))| *z = m.lazy_reduce_mul_add(x, y, *z));
            });
        });

        group.bench_function(format!("slice/n={n}"), |bencher| {
            bencher.iter(|| {
                let acc = black_box(&mut acc);
                let a = black_box(&a);
                let b = black_box(&b);
                m.lazy_reduce_add_mul_slice_assign(acc, a, b);
            });
        });
    }
    group.finish();
}

fn bench_mul_add_to(c: &mut Criterion) {
    let m = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("mul_add_slice_to");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let c_in = rand_vec(n, MODULUS);
        let mut out = vec![0; n];

        group.bench_function(format!("per_elem/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                let c = black_box(&c_in);
                out.iter_mut()
                    .zip(a.iter().zip(b.iter()).zip(c.iter()))
                    .for_each(|(d, ((&a, &b), &c))| *d = m.reduce_mul_add(a, b, c));
            });
        });

        group.bench_function(format!("slice/n={n}"), |bencher| {
            bencher.iter(|| {
                let out = black_box(&mut out);
                let a = black_box(&a);
                let b = black_box(&b);
                let c = black_box(&c_in);
                m.reduce_mul_add_slice_to(a, b, c, out);
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_add,
    bench_mul,
    bench_add_mul_assign,
    bench_lazy_add_mul_assign,
    bench_mul_add_to,
);
criterion_main!(benches);
