// cargo r -r -p primus_distr --example compare_samplers
//
// Comparative analysis of multiple discrete Gaussian samplers.
// Compares accuracy and characteristics across different implementations:
// - CDTSamplerLogSpace (f64 precision, 128-bit sampling)
// - CDTSamplerLogSpaceDD (double-double ~106-bit precision, 128-bit sampling)

use comfy_table::{Attribute, Cell, Color, ContentArrangement, Table, presets::UTF8_FULL};
use rand::distr::Distribution;
use std::time::Instant;

type ValueT = u64;

// Modulus for discrete Gaussian sampling
const Q: ValueT = 1125899906826241;
const HALF_Q: ValueT = Q >> 1;

// Number of samples for statistical analysis (2^20 = 1,048,576)
const N: usize = 1 << 20;

// Tail cut for discrete Gaussian
const TAIL_CUT: f64 = 12.0;

// Sigma ranges to test for cumulative probability
const SIGMA_RANGES: [f64; 6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

fn main() {
    let sigmas: Vec<f64> = vec![0.8, 1.5, 3.19, 9.0, 13.0, 15.0, 16.0, 17.0];

    println!("\n{}", "═".repeat(100));
    println!("Discrete Gaussian Sampler Comparison");
    println!("Samples per test: {}", N);
    println!("Testing {} sigma values: {:?}", sigmas.len(), sigmas);
    println!("{}\n", "═".repeat(100));

    for sigma in sigmas {
        compare_samplers_at_sigma(sigma);
    }

    println!("\n{}", "═".repeat(100));
    println!("All comparisons completed!");
    println!("{}", "═".repeat(100));
}

/// Helper function to convert modular value to signed value
#[inline]
fn to_signed_i128(x: ValueT) -> i128 {
    if x <= HALF_Q {
        x as i128
    } else {
        (x as i128) - (Q as i128)
    }
}

/// Helper function to convert modular value to signed f64
#[inline]
fn to_signed_f64(x: ValueT) -> f64 {
    if x <= HALF_Q {
        x as f64
    } else {
        (x as f64) - (Q as f64)
    }
}

struct SamplerStats {
    name: String,
    actual_std: f64,
    std_error: f64,
    std_error_pct: f64,
    cumulative_probs: Vec<f64>,
    sample_time_ms: f64,
}

fn compare_samplers_at_sigma(sigma: f64) {
    println!("\n{}", "═".repeat(100));
    println!("Comparing samplers at σ = {:.2}", sigma);
    println!("{}", "═".repeat(100));

    let mut rng = rand::rng();
    let mut all_stats = Vec::new();

    let mut i = 1;

    // Test CDTSampler (f64 precision)
    {
        println!("\n[{i}] Testing CDTSampler ...");
        let sampler = primus_distr::CDTSampler::<ValueT>::new(sigma, TAIL_CUT, Q - 1);

        let start = Instant::now();
        let data: Vec<ValueT> = sampler.sample_iter(&mut rng).take(N).collect();
        let sample_time_ms = start.elapsed().as_secs_f64() * 1000.0;

        let stats = compute_stats("CDTSampler (f64)", &data, sigma, sample_time_ms);
        all_stats.push(stats);
        i += 1;
    }

    #[allow(unused_assignments)]
    {
        println!("[{i}] Testing Discrete Ziggurat ...");
        let sampler = primus_distr::DiscreteZiggurat::<ValueT>::new(sigma, TAIL_CUT, Q - 1);

        let start = Instant::now();
        let data: Vec<ValueT> = sampler.sample_iter(&mut rng).take(N).collect();
        let sample_time_ms = start.elapsed().as_secs_f64() * 1000.0;

        let stats = compute_stats("Discrete Ziggurat", &data, sigma, sample_time_ms);
        all_stats.push(stats);
        i += 1;
    }

    #[cfg(all(target_os = "linux", feature = "high_precision"))]
    {
        println!("[{i}] Testing UnixCDTSampler ...");
        let sampler = primus_distr::UnixCDTSampler::<ValueT>::new(sigma, TAIL_CUT, Q - 1);

        let start = Instant::now();
        let data: Vec<ValueT> = sampler.sample_iter(&mut rng).take(N).collect();
        let sample_time_ms = start.elapsed().as_secs_f64() * 1000.0;

        let stats = compute_stats("UnixCDTSampler", &data, sigma, sample_time_ms);
        all_stats.push(stats);
    }

    // Display comparison table
    display_comparison_table(sigma, &all_stats);
}

fn compute_stats(name: &str, data: &[ValueT], sigma: f64, sample_time_ms: f64) -> SamplerStats {
    // Pre-compute sigma limits for cumulative probability counting
    let sigma_limits: [u64; 6] = SIGMA_RANGES
        .iter()
        .map(|n| (n * sigma).floor() as u64)
        .collect::<Vec<_>>()
        .try_into()
        .expect("SIGMA_RANGES should have exactly 6 elements");

    // Single-pass computation: mean, cumulative counts
    let mut sum = 0i128;
    let mut counts = [0usize; 6];

    for &x in data.iter() {
        let signed_val = to_signed_i128(x);
        sum += signed_val;

        // Count for cumulative probabilities
        let abs_val = signed_val.unsigned_abs() as u64;
        for (i, &limit) in sigma_limits.iter().enumerate() {
            if abs_val <= limit {
                counts[i] += 1;
            }
        }
    }

    // Calculate mean using f64
    let mean_f64 = sum as f64 / N as f64;

    // Calculate variance
    let variance_sum: f64 = data
        .iter()
        .map(|&x| {
            let diff = to_signed_f64(x) - mean_f64;
            diff * diff
        })
        .sum();

    let actual_std = (variance_sum / N as f64).sqrt();
    let std_error = actual_std - sigma;
    let std_error_pct = (std_error / sigma) * 100.0;

    // Compute cumulative probabilities
    let cumulative_probs: Vec<f64> = counts.iter().map(|&c| c as f64 / N as f64).collect();

    SamplerStats {
        name: name.to_string(),
        actual_std,
        std_error,
        std_error_pct,
        cumulative_probs,
        sample_time_ms,
    }
}

fn display_comparison_table(sigma: f64, all_stats: &[SamplerStats]) {
    // Compute theoretical probabilities
    let gaussian_pdf = |k: i64| -> f64 {
        let k_f = k as f64;
        (-k_f * k_f / (2.0 * sigma * sigma)).exp()
    };

    let norm_limit = (sigma * TAIL_CUT).ceil() as i64;
    let mut z = gaussian_pdf(0);
    for k in 1..=norm_limit {
        z += 2.0 * gaussian_pdf(k);
    }

    let mut theoretical_probs = Vec::with_capacity(6);
    for &n_sigma in SIGMA_RANGES.iter() {
        let limit = (n_sigma * sigma).floor() as i64;
        let mut prob = gaussian_pdf(0);
        for k in 1..=limit {
            prob += 2.0 * gaussian_pdf(k);
        }
        theoretical_probs.push(prob / z);
    }

    // Standard Deviation Comparison Table
    println!("\n{}", "━".repeat(100));
    println!("Standard Deviation Comparison (σ = {:.2})", sigma);
    println!("{}", "━".repeat(100));

    let mut std_table = Table::new();
    std_table
        .load_preset(UTF8_FULL)
        .set_content_arrangement(ContentArrangement::Dynamic);

    std_table.set_header(vec![
        Cell::new("Sampler").add_attribute(Attribute::Bold),
        Cell::new("Expected σ").add_attribute(Attribute::Bold),
        Cell::new("Actual σ").add_attribute(Attribute::Bold),
        Cell::new("Error (abs)").add_attribute(Attribute::Bold),
        Cell::new("Error (%)").add_attribute(Attribute::Bold),
        Cell::new("Sample Time (ms)").add_attribute(Attribute::Bold),
    ]);

    for stats in all_stats {
        let error_color = if stats.std_error_pct.abs() < 0.1 {
            Color::Green
        } else if stats.std_error_pct.abs() < 0.5 {
            Color::Yellow
        } else {
            Color::Red
        };

        std_table.add_row(vec![
            Cell::new(&stats.name),
            Cell::new(format!("{:.10}", sigma)),
            Cell::new(format!("{:.10}", stats.actual_std)),
            Cell::new(format!("{:+.10}", stats.std_error)).fg(error_color),
            Cell::new(format!("{:+.4}%", stats.std_error_pct)).fg(error_color),
            Cell::new(format!("{:.2}", stats.sample_time_ms)),
        ]);
    }

    println!("{}", std_table);

    // Cumulative Probability Comparison Table
    println!("\n{}", "━".repeat(100));
    println!("Cumulative Probability Comparison (σ = {:.2})", sigma);
    println!("{}", "━".repeat(100));

    for (idx, &n_sigma) in SIGMA_RANGES.iter().enumerate() {
        let mut prob_table = Table::new();
        prob_table
            .load_preset(UTF8_FULL)
            .set_content_arrangement(ContentArrangement::Dynamic);

        let expected = theoretical_probs[idx];

        prob_table.set_header(vec![
            Cell::new(format!("±{}σ (±{:.2})", n_sigma, n_sigma * sigma))
                .add_attribute(Attribute::Bold),
            Cell::new("Actual P").add_attribute(Attribute::Bold),
            Cell::new("Expected P").add_attribute(Attribute::Bold),
            Cell::new("Diff (abs)").add_attribute(Attribute::Bold),
            Cell::new("Diff (%)").add_attribute(Attribute::Bold),
        ]);

        for stats in all_stats {
            let actual = stats.cumulative_probs[idx];
            let diff = actual - expected;
            let diff_pct = if expected > 0.0 {
                diff / expected * 100.0
            } else {
                0.0
            };

            let diff_color = if diff.abs() > 0.01 {
                Color::Red
            } else if diff.abs() > 0.005 {
                Color::Yellow
            } else {
                Color::Green
            };

            let diff_pct_color = if diff_pct.abs() > 1.0 {
                Color::Red
            } else if diff_pct.abs() > 0.5 {
                Color::Yellow
            } else {
                Color::Green
            };

            prob_table.add_row(vec![
                Cell::new(&stats.name),
                Cell::new(format!("{:.6} ({:.2}%)", actual, actual * 100.0)),
                Cell::new(format!("{:.6} ({:.2}%)", expected, expected * 100.0)),
                Cell::new(format!("{:+.6}", diff)).fg(diff_color),
                Cell::new(format!("{:+.2}%", diff_pct)).fg(diff_pct_color),
            ]);
        }

        println!("{}", prob_table);
    }

    // Summary comparison
    println!("\n{}", "━".repeat(100));
    println!("Summary Comparison (σ = {:.2})", sigma);
    println!("{}", "━".repeat(100));

    let mut summary_table = Table::new();
    summary_table
        .load_preset(UTF8_FULL)
        .set_content_arrangement(ContentArrangement::Dynamic);

    summary_table.set_header(vec![
        Cell::new("Sampler").add_attribute(Attribute::Bold),
        Cell::new("Avg Prob Error (%)").add_attribute(Attribute::Bold),
        Cell::new("Max Prob Error (%)").add_attribute(Attribute::Bold),
        Cell::new("σ Error (%)").add_attribute(Attribute::Bold),
        Cell::new("Overall Quality").add_attribute(Attribute::Bold),
    ]);

    for stats in all_stats {
        // Compute average and max probability errors
        let prob_errors: Vec<f64> = stats
            .cumulative_probs
            .iter()
            .zip(theoretical_probs.iter())
            .map(|(&actual, &expected)| {
                if expected > 0.0 {
                    ((actual - expected) / expected * 100.0).abs()
                } else {
                    0.0
                }
            })
            .collect();

        let avg_prob_error = prob_errors.iter().sum::<f64>() / prob_errors.len() as f64;
        let max_prob_error = prob_errors
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);

        let quality = if stats.std_error_pct.abs() < 0.1 && avg_prob_error < 0.1 {
            Cell::new("★★★ Excellent").fg(Color::Green)
        } else if stats.std_error_pct.abs() < 0.5 && avg_prob_error < 0.5 {
            Cell::new("★★☆ Very Good").fg(Color::Cyan)
        } else if stats.std_error_pct.abs() < 1.0 && avg_prob_error < 1.0 {
            Cell::new("★☆☆ Good").fg(Color::Yellow)
        } else {
            Cell::new("☆☆☆ Acceptable").fg(Color::Red)
        };

        summary_table.add_row(vec![
            Cell::new(&stats.name),
            Cell::new(format!("{:.4}%", avg_prob_error)),
            Cell::new(format!("{:.4}%", max_prob_error)),
            Cell::new(format!("{:+.4}%", stats.std_error_pct.abs())),
            quality,
        ]);
    }

    println!("{}", summary_table);
    println!("{}", "━".repeat(100));
}
