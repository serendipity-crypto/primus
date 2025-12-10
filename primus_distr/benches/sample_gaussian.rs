// cargo bench -p primus_distr --bench sample_gaussian

use criterion::{Criterion, criterion_group, criterion_main};

use primus_distr::DiscreteZiggurat;
use rand::distr::Distribution;

const N: usize = 100_000;
const MODULUS_MINUS_ONE: u64 = 1125899906826241 - 1;

fn bench_sample(c: &mut Criterion) {
    let mut group = c.benchmark_group("Sample");

    for sigma in [3.19, 10.0, 16.0, 17.0, 20.0, 25.0] {
        if sigma >= 10.0 {
            let mut rng = rand::rng();
            let ziggurat = DiscreteZiggurat::new(sigma, 12.0, MODULUS_MINUS_ONE);
            group.bench_function(format!("DiscreteZiggurat({sigma})"), |b| {
                b.iter(|| {
                    for _ in 0..N {
                        ziggurat.sample(&mut rng);
                    }
                })
            });
        }

        if sigma <= 20.0 {
            use primus_distr::CDTSampler;

            let mut rng = rand::rng();
            let cdt = CDTSampler::new(sigma, 12.0, MODULUS_MINUS_ONE);
            group.bench_function(format!("CDTSampler({sigma})"), |b| {
                b.iter(|| {
                    for _ in 0..N {
                        cdt.sample(&mut rng);
                    }
                })
            });
        }

        #[cfg(all(target_os = "linux", feature = "high_precision"))]
        if sigma <= 20.0 {
            use primus_distr::UnixCDTSampler;

            let mut rng = rand::rng();
            let unix_cdt = UnixCDTSampler::new(sigma, 12.0, MODULUS_MINUS_ONE);
            group.bench_function(format!("UnixCDTSampler({sigma})"), |b| {
                b.iter(|| {
                    for _ in 0..N {
                        unix_cdt.sample(&mut rng);
                    }
                })
            });
        }
    }

    group.finish();
}

criterion_group!(benches, bench_sample);

criterion_main!(benches);
