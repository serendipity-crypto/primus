// cargo r -r -p primus_distr --example check_gaussian
//
// Test harness for discrete Gaussian samplers.
// Validates:
// - Standard deviation accuracy
// - Cumulative probability distribution across sigma ranges

use comfy_table::{Attribute, Cell, Color, ContentArrangement, Table, presets::UTF8_FULL};
use num_traits::ConstZero;
use rand::distr::Distribution;

type ValueT = u64;

// Modulus for discrete Gaussian sampling
const Q: ValueT = 1125899906826241;
const HALF_Q: ValueT = Q >> 1;

// Number of samples for statistical analysis (2^20 = 1,048,576)
const N: usize = 1 << 20;

// Number of distributions to sum for convolution test
const CHUNK_SIZE: usize = 10;

// Tail cut for discrete Gaussian (must match DiscreteZiggurat::new parameter)
const TAIL_CUT: f64 = 12.0;

// Sigma ranges to test for cumulative probability
const SIGMA_RANGES: [f64; 6] = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

// Quality thresholds for standard deviation error (percentage)
const QUALITY_EXCELLENT: f64 = 0.1;
const QUALITY_VERY_GOOD: f64 = 0.5;
const QUALITY_GOOD: f64 = 1.0;
const QUALITY_ACCEPTABLE: f64 = 2.0;

// Probability difference thresholds for color coding
const PROB_DIFF_RED_THRESHOLD: f64 = 0.01;      // 1% absolute difference
const PROB_DIFF_YELLOW_THRESHOLD: f64 = 0.005;  // 0.5% absolute difference
const PROB_DIFF_PCT_RED_THRESHOLD: f64 = 1.0;   // 1% relative difference
const PROB_DIFF_PCT_YELLOW_THRESHOLD: f64 = 0.5; // 0.5% relative difference

fn main() {
    check_standard_deviation();
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

/// Determine quality level based on error percentage
fn quality_level(error_pct: f64) -> &'static str {
    let abs_error = error_pct.abs();
    if abs_error < QUALITY_EXCELLENT {
        "✓ Excellent"
    } else if abs_error < QUALITY_VERY_GOOD {
        "✓ Very Good"
    } else if abs_error < QUALITY_GOOD {
        "○ Good"
    } else if abs_error < QUALITY_ACCEPTABLE {
        "△ Acceptable"
    } else {
        "✗ Poor"
    }
}

/// Create colored cell based on difference value
fn colored_diff_cell(diff: f64) -> Cell {
    let color = if diff.abs() > PROB_DIFF_RED_THRESHOLD {
        Color::Red
    } else if diff.abs() > PROB_DIFF_YELLOW_THRESHOLD {
        Color::Yellow
    } else {
        Color::Green
    };
    Cell::new(format!("{:+.6}", diff)).fg(color)
}

/// Create colored cell based on percentage difference
fn colored_pct_cell(diff_pct: f64) -> Cell {
    let color = if diff_pct.abs() > PROB_DIFF_PCT_RED_THRESHOLD {
        Color::Red
    } else if diff_pct.abs() > PROB_DIFF_PCT_YELLOW_THRESHOLD {
        Color::Yellow
    } else {
        Color::Green
    };
    Cell::new(format!("{:+.2}%", diff_pct)).fg(color)
}

fn check_standard_deviation() {
    let mut rng = rand::rng();

    // Test sigma values (including the previously problematic 3.19)
    let sigmas: Vec<f64> = vec![1.9, 3.19, 5.0, 10.0];

    println!("\n{}", "═".repeat(80));
    println!("Discrete Gaussian Sampler Validation");
    println!("Samples per test: {}", N);
    println!("Testing {} sigma values: {:?}", sigmas.len(), sigmas);
    println!("{}\n", "═".repeat(80));

    let mut data: Vec<ValueT> = vec![ValueT::ZERO; N];
    for (idx, sigma) in sigmas.iter().enumerate() {
        println!("[{}/{}] Testing σ = {:.2}...", idx + 1, sigmas.len(), sigma);
        let distr = <primus_distr::DiscreteZiggurat<ValueT>>::new(*sigma, TAIL_CUT, Q - 1);

        // Sample data
        data.iter_mut()
            .zip(distr.clone().sample_iter(&mut rng))
            .for_each(|(d, v)| *d = v);

        check(&data, *sigma, "");

        // Test convolution property: sum of N(0, σ²) should be N(0, n*σ²)
        println!("  Testing convolution property (sum of {} distributions)...", CHUNK_SIZE);
        let datas: Vec<Vec<ValueT>> = (0..CHUNK_SIZE)
            .map(|_| distr.clone().sample_iter(&mut rng).take(N).collect())
            .collect();

        let new_data = datas
            .into_iter()
            .reduce(|mut acc, x| {
                for (a, b) in acc.iter_mut().zip(x) {
                    *a = (*a + b) % Q;
                }
                acc
            })
            .expect("reduce should never fail with non-empty iterator");

        check(&new_data, (CHUNK_SIZE as f64).sqrt() * sigma, "Convolution ");
    }

    println!("\n{}", "═".repeat(80));
    println!("All tests completed!");
    println!("{}", "═".repeat(80));
}

fn check(data: &[ValueT], sigma: f64, info: &str) {
    // Pre-compute sigma limits for cumulative probability counting
    let sigma_limits: [u64; 6] = SIGMA_RANGES
        .iter()
        .map(|n| (n * sigma).floor() as u64)
        .collect::<Vec<_>>()
        .try_into()
        .expect("SIGMA_RANGES should have exactly 6 elements");

    // Single-pass computation: mean, cumulative counts
    // Using f64 for mean/variance is 100-1000x faster than BigDecimal
    let mut sum = 0i128; // Use i128 to avoid overflow
    let mut counts = [0usize; 6]; // Counts for each sigma range

    for &x in data.iter() {
        let signed_val = to_signed_i128(x);
        sum += signed_val;

        // Count for cumulative probabilities (use absolute value)
        let abs_val = signed_val.unsigned_abs() as u64;
        for (i, &limit) in sigma_limits.iter().enumerate() {
            if abs_val <= limit {
                counts[i] += 1;
            }
        }
    }

    // Calculate mean using f64 (much faster than BigDecimal)
    let mean_f64 = sum as f64 / N as f64;

    // Calculate variance in second pass
    let variance_sum: f64 = data
        .iter()
        .map(|&x| {
            let diff = to_signed_f64(x) - mean_f64;
            diff * diff
        })
        .sum();

    let real_std = (variance_sum / N as f64).sqrt();

    // Standard deviation comparison
    let std_error = real_std - sigma;
    let std_error_pct = (std_error / sigma) * 100.0;

    println!("\n{}", "━".repeat(80));
    println!("  {}Standard Deviation Analysis (σ)", info);
    println!("{}", "━".repeat(80));
    println!("  Expected:  {:.10}", sigma);
    println!("  Actual:    {:.10}", real_std);
    println!("  Error:     {:+.10} (absolute)", std_error);
    println!("             {:+.4}% (relative)", std_error_pct);
    println!("  Quality:   {}", quality_level(std_error_pct));
    println!("{}", "━".repeat(80));

    // Pre-compute theoretical probabilities (avoid recomputation)
    let gaussian_pdf = |k: i64| -> f64 {
        let k_f = k as f64;
        (-k_f * k_f / (2.0 * sigma * sigma)).exp()
    };

    // Compute normalization constant once (use TAIL_CUT to match sampler truncation)
    let norm_limit = (sigma * TAIL_CUT).ceil() as i64;
    let mut z = gaussian_pdf(0);
    for k in 1..=norm_limit {
        z += 2.0 * gaussian_pdf(k);
    }

    // Pre-compute cumulative probabilities for all sigma ranges (avoid repeated calculation)
    let mut theoretical_probs = Vec::with_capacity(6);
    for &n_sigma in SIGMA_RANGES.iter() {
        let limit = (n_sigma * sigma).floor() as i64;
        let mut prob = gaussian_pdf(0);
        for k in 1..=limit {
            prob += 2.0 * gaussian_pdf(k);
        }
        theoretical_probs.push(prob / z);
    }

    // Actual cumulative probabilities already computed in counts array
    let actual_probs: Vec<f64> = counts.iter().map(|&c| c as f64 / N as f64).collect();

    // Create sigma range comparison table
    let mut table = Table::new();
    table
        .load_preset(UTF8_FULL)
        .set_content_arrangement(ContentArrangement::Dynamic);

    table.set_header(vec![
        Cell::new("Range").add_attribute(Attribute::Bold),
        Cell::new("σ Value").add_attribute(Attribute::Bold),
        Cell::new("Actual P(|X| ≤ nσ)").add_attribute(Attribute::Bold),
        Cell::new("Expected P(|X| ≤ nσ)").add_attribute(Attribute::Bold),
        Cell::new("Diff").add_attribute(Attribute::Bold),
        Cell::new("Diff %").add_attribute(Attribute::Bold),
    ]);

    // Use pre-computed probabilities
    for i in 0..6 {
        let n_sigma = SIGMA_RANGES[i];
        let actual = actual_probs[i];
        let expected = theoretical_probs[i];
        let diff = actual - expected;
        let diff_pct = if expected > 0.0 {
            diff / expected * 100.0
        } else {
            0.0
        };

        table.add_row(vec![
            Cell::new(format!("±{}σ", n_sigma)),
            Cell::new(format!("±{:.2}", n_sigma * sigma)),
            Cell::new(format!("{:.6} ({:.2}%)", actual, actual * 100.0)),
            Cell::new(format!("{:.6} ({:.2}%)", expected, expected * 100.0)),
            colored_diff_cell(diff),
            colored_pct_cell(diff_pct),
        ]);
    }

    println!("\n{}", "=".repeat(80));
    println!("Cumulative Probability Distribution (σ = {:.2}):", sigma);
    println!("{}", "=".repeat(80));
    println!("{}", table);
    println!("{}\n", "=".repeat(80));
}
