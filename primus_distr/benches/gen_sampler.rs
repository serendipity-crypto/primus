// cargo bench -p primus_distr --bench gen_sampler

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use primus_distr::DiscreteZiggurat;

const MODULUS_MINUS_ONE: u64 = 1125899906826241 - 1;

const TAIL_CUT: f64 = 12.0;

fn bench_different_sampler(c: &mut Criterion) {
    let mut group = c.benchmark_group("GenSampler");

    for sigma in [1.0, 3.0, 10.0, 15.0, 20.0, 25.0] {
        if sigma >= 10.0 {
            group.bench_function(format!("DiscreteZiggurat({sigma})"), |b| {
                b.iter(|| {
                    black_box(DiscreteZiggurat::new(sigma, TAIL_CUT, MODULUS_MINUS_ONE));
                })
            });
        }

        if sigma <= 20.0 {
            use primus_distr::CDTSampler;

            group.bench_function(format!("CDTSampler({sigma})"), |b| {
                b.iter(|| black_box(CDTSampler::new(sigma, TAIL_CUT, MODULUS_MINUS_ONE)))
            });
        }

        #[cfg(all(target_os = "linux", feature = "high_precision"))]
        if sigma <= 20.0 {
            use primus_distr::UnixCDTSampler;
            group.bench_function(format!("UnixCDTSampler({sigma})"), |b| {
                b.iter(|| black_box(UnixCDTSampler::new(sigma, TAIL_CUT, MODULUS_MINUS_ONE)))
            });
        }
    }

    group.finish();
}

criterion_group!(benches, bench_different_sampler);

criterion_main!(benches);
