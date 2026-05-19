//! Step 9 验证 bench：`BarrettModulus<u64>` portable_simd 路径在 `N = 4 / 8 / 16`
//! 三个 lane width 下的吞吐对比。
//!
//! plan.md Step 9 想问的问题：在 AVX-512 主机上，调用 simd_kernel 时把 lane
//! width 从默认的 `default_lanes` (= 8 on AVX-512) 显式覆盖成 `16`，是否真能
//! 给 NTT 大尺寸场景多 +10~15%。这里 bench 4 个 mul-class 算子来回答。
//!
//! 用 ~2^62 NTT-friendly prime（与 `barrett_slice.rs` 一致），所以走的是
//! `BarrettModulus<u64>` 的 portable_simd 通路，**不是 Barrett50 IFMA**。

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use primus_modulus::{BarrettModulus, barrett_simd_kernel as kern};
use rand::distr::{Distribution, Uniform};

type V = u64;

const MODULUS: V = 4_611_686_018_427_322_369;

fn rand_vec(n: usize, m: V) -> Vec<V> {
    let mut rng = rand::rng();
    let distr = Uniform::new(0, m).unwrap();
    distr.sample_iter(&mut rng).take(n).collect()
}

macro_rules! bench_lane {
    ($group:expr, $kern:ident, $lanes:literal, $n:expr, $modulus:expr, $($arg:expr),* $(,)?) => {
        $group.bench_function(format!("N={}/n={}", $lanes, $n), |bencher| {
            bencher.iter(|| {
                kern::$kern::<u64, $lanes>(
                    $modulus,
                    $(black_box($arg)),*
                );
            });
        });
    };
}

fn bench_mul_slice_to(c: &mut Criterion) {
    let m = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("lane_width_mul_slice_to");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let mut out = vec![0u64; n];
        bench_lane!(group, reduce_mul_slice_to, 4, n, m, &a, &b, &mut out);
        bench_lane!(group, reduce_mul_slice_to, 8, n, m, &a, &b, &mut out);
        bench_lane!(group, reduce_mul_slice_to, 16, n, m, &a, &b, &mut out);
    }
    group.finish();
}

fn bench_lazy_mul_slice_to(c: &mut Criterion) {
    let m = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("lane_width_lazy_mul_slice_to");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let mut out = vec![0u64; n];
        bench_lane!(group, lazy_reduce_mul_slice_to, 4, n, m, &a, &b, &mut out);
        bench_lane!(group, lazy_reduce_mul_slice_to, 8, n, m, &a, &b, &mut out);
        bench_lane!(group, lazy_reduce_mul_slice_to, 16, n, m, &a, &b, &mut out);
    }
    group.finish();
}

fn bench_mul_add_slice_to(c: &mut Criterion) {
    let m = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("lane_width_mul_add_slice_to");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        let c_in = rand_vec(n, MODULUS);
        let mut out = vec![0u64; n];
        bench_lane!(group, reduce_mul_add_slice_to, 4, n, m, &a, &b, &c_in, &mut out);
        bench_lane!(group, reduce_mul_add_slice_to, 8, n, m, &a, &b, &c_in, &mut out);
        bench_lane!(group, reduce_mul_add_slice_to, 16, n, m, &a, &b, &c_in, &mut out);
    }
    group.finish();
}

fn bench_dot_product(c: &mut Criterion) {
    let m = BarrettModulus::<V>::new(MODULUS);
    let mut group = c.benchmark_group("lane_width_dot_product");
    for &n in &[1024usize, 4096, 16384] {
        let a = rand_vec(n, MODULUS);
        let b = rand_vec(n, MODULUS);
        group.bench_function(format!("N=4/n={n}"), |bencher| {
            bencher.iter(|| {
                let a = black_box(&a);
                let b = black_box(&b);
                kern::reduce_dot_product::<u64, 4>(m, a, b)
            });
        });
        group.bench_function(format!("N=8/n={n}"), |bencher| {
            bencher.iter(|| {
                let a = black_box(&a);
                let b = black_box(&b);
                kern::reduce_dot_product::<u64, 8>(m, a, b)
            });
        });
        group.bench_function(format!("N=16/n={n}"), |bencher| {
            bencher.iter(|| {
                let a = black_box(&a);
                let b = black_box(&b);
                kern::reduce_dot_product::<u64, 16>(m, a, b)
            });
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_mul_slice_to,
    bench_lazy_mul_slice_to,
    bench_mul_add_slice_to,
    bench_dot_product,
);
criterion_main!(benches);
