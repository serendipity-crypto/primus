use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
#[cfg(not(target_os = "linux"))]
use primus_distr::CDTSampler;
use primus_distr::DiscreteZiggurat;
#[cfg(target_os = "linux")]
use primus_distr::UnixCDTSampler;
use rand::distr::Distribution;

const N: usize = 100_000;
const MODULUS_MINUS_ONE: u64 = 1125899906826241 - 1;

fn bench_discrete_ziggurat(c: &mut Criterion) {
    let mut group = c.benchmark_group("Discrete Ziggurat");

    for std_dev in [1.0, 2.0, 2.5, 3.0, 5.0, 10.0, 20.0] {
        let mut rng = rand::rng();
        let sampler = DiscreteZiggurat::new(std_dev, 12.0, MODULUS_MINUS_ONE);

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("std_dev={}", std_dev)),
            &std_dev,
            |b, _| {
                b.iter(|| {
                    for _ in 0..N {
                        sampler.sample(&mut rng);
                    }
                })
            },
        );
    }

    group.finish();
}

#[cfg(not(target_os = "linux"))]
fn bench_cdt_sampler(c: &mut Criterion) {
    let mut group = c.benchmark_group("CDT Sampler");

    for std_dev in [1.0, 2.0, 2.5, 3.0, 5.0, 10.0, 20.0] {
        let tail_cut = if std_dev < 3.0 { 12.0 } else { 6.0 };
        let mut rng = rand::rng();
        let sampler = CDTSampler::new(std_dev, tail_cut, MODULUS_MINUS_ONE);

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("std_dev={}", std_dev)),
            &std_dev,
            |b, _| {
                b.iter(|| {
                    for _ in 0..N {
                        sampler.sample(&mut rng);
                    }
                })
            },
        );
    }

    group.finish();
}

#[cfg(target_os = "linux")]
fn bench_unix_cdt_sampler(c: &mut Criterion) {
    let mut group = c.benchmark_group("Unix CDT Sampler");

    for std_dev in [1.0, 2.0, 2.5, 3.0, 5.0, 10.0, 20.0] {
        let tail_cut = if std_dev < 3.0 { 12.0 } else { 6.0 };
        let mut rng = rand::rng();
        let sampler = UnixCDTSampler::new(std_dev, tail_cut, MODULUS_MINUS_ONE);

        group.bench_with_input(
            BenchmarkId::from_parameter(format!("std_dev={}", std_dev)),
            &std_dev,
            |b, _| {
                b.iter(|| {
                    for _ in 0..N {
                        sampler.sample(&mut rng);
                    }
                })
            },
        );
    }

    group.finish();
}

fn bench_different_sampler(c: &mut Criterion) {
    let mut group = c.benchmark_group("Compare_Same_Sigma");

    for std_dev in [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1] {
        let mut rng = rand::rng();
        let ziggurat = DiscreteZiggurat::new(std_dev, 12.0, MODULUS_MINUS_ONE);
        group.bench_function(format!("DiscreteZiggurat_std={std_dev}"), |b| {
            b.iter(|| {
                for _ in 0..N {
                    ziggurat.sample(&mut rng);
                }
            })
        });

        #[cfg(not(target_os = "linux"))]
        {
            let mut rng = rand::rng();
            let cdt = CDTSampler::new(std_dev, 12.0, MODULUS_MINUS_ONE);
            group.bench_function(format!("CDTSampler_std={std_dev}"), |b| {
                b.iter(|| {
                    for _ in 0..N {
                        cdt.sample(&mut rng);
                    }
                })
            });
        }

        #[cfg(target_os = "linux")]
        {
            let mut rng = rand::rng();
            let unix_cdt = UnixCDTSampler::new(std_dev, 12.0, MODULUS_MINUS_ONE);
            group.bench_function(format!("UnixCDTSampler_std={std_dev}"), |b| {
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

#[cfg(not(target_os = "linux"))]
criterion_group!(
    benches,
    bench_discrete_ziggurat,
    bench_cdt_sampler,
    bench_different_sampler,
);

#[cfg(target_os = "linux")]
criterion_group!(
    benches,
    bench_discrete_ziggurat,
    bench_unix_cdt_sampler,
    bench_different_sampler,
);

criterion_main!(benches);
