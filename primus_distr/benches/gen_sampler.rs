use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
#[cfg(target_os = "linux")]
use primus_distr::UnixCDTSampler;
use primus_distr::{CDTSamplerLogSpaceDD, DiscreteZiggurat};

const MODULUS_MINUS_ONE: u64 = 1125899906826241 - 1;

const TAIL_CUT: f64 = 12.0;

fn bench_different_sampler(c: &mut Criterion) {
    let mut group = c.benchmark_group("GenSampler");

    for sigma in [2.0, 2.5, 3.0] {
        if sigma >= 1.5 {
            group.bench_function(format!("DiscreteZiggurat({sigma})"), |b| {
                b.iter(|| {
                    black_box(DiscreteZiggurat::new(sigma, TAIL_CUT, MODULUS_MINUS_ONE));
                })
            });
        }

        {
            use primus_distr::CDTSamplerLogSpace;

            group.bench_function(format!("CDTSamplerLogSpace({sigma})"), |b| {
                b.iter(|| black_box(CDTSamplerLogSpace::new(sigma, TAIL_CUT, MODULUS_MINUS_ONE)))
            });
        }

        {
            group.bench_function(format!("CDTSamplerLogSpaceDD({sigma})"), |b| {
                b.iter(|| {
                    black_box(CDTSamplerLogSpaceDD::new(
                        sigma,
                        TAIL_CUT,
                        MODULUS_MINUS_ONE,
                    ))
                })
            });
        }

        // #[cfg(not(target_os = "linux"))]
        // {
        //     use primus_distr::CDTSampler;
        //     group.bench_function(format!("CDTSampler({sigma})"), |b| {
        //         b.iter(|| black_box(CDTSampler::new(sigma, TAIL_CUT, MODULUS_MINUS_ONE)))
        //     });
        // }

        #[cfg(target_os = "linux")]
        {
            group.bench_function(format!("UnixCDTSampler_std={std_dev}"), |b| {
                b.iter(|| black_box(UnixCDTSampler::new(sigma, TAIL_CUT, MODULUS_MINUS_ONE)))
            });
        }
    }

    group.finish();
}

criterion_group!(benches, bench_different_sampler,);

criterion_main!(benches);
