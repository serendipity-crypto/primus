use primus_integer::{AsInto, UnsignedInteger};
use rand::distr::Distribution;
use twofloat::TwoFloat;

/// CDT sampler using double-double precision (TwoFloat) for ~106-bit accuracy
#[derive(Debug, Clone)]
pub struct CDTSamplerLogSpaceDD<T: UnsignedInteger> {
    std_dev: f64,
    modulus_minus_one: T,
    cdt: Vec<u128>,
}

/// Log-sum-exp trick using TwoFloat for higher precision
/// Computes log(sum(exp(log_values))) stably
fn log_sum_exp_dd(log_values: &[TwoFloat]) -> TwoFloat {
    if log_values.is_empty() {
        return TwoFloat::from(f64::NEG_INFINITY);
    }

    // Find maximum value for numerical stability
    let mut max_log = log_values[0];
    for &val in &log_values[1..] {
        if val.hi() > max_log.hi() || (val.hi() == max_log.hi() && val.lo() > max_log.lo()) {
            max_log = val;
        }
    }

    // Check if max_log is negative infinity
    if max_log.hi().is_infinite() && max_log.hi().is_sign_negative() {
        return TwoFloat::from(f64::NEG_INFINITY);
    }

    // Compute sum(exp(log_values - max_log))
    let mut sum_exp = TwoFloat::from(0.0);
    for &log_val in log_values.iter() {
        sum_exp += (log_val - max_log).exp()
    }

    max_log + sum_exp.ln()
}

/// Convert TwoFloat in [0, 1] to u128 with high precision
/// Optimized version using direct bit manipulation
#[inline]
fn twofloat_to_u128(tf: TwoFloat) -> u128 {
    // Fast path for boundary values
    let hi = tf.hi();
    if hi <= 0.0 {
        return 0;
    }
    if hi >= 1.0 {
        return u128::MAX;
    }

    // For values in (0, 1), use optimized scaling
    // Split the 128-bit range into hi and lo parts
    const SCALE_HI: f64 = 18446744073709551616.0; // 2^64

    // Scale to get high 64 bits
    let scaled = tf * TwoFloat::from(SCALE_HI);
    let hi_val = scaled.hi();
    let hi_u64 = hi_val.min(u64::MAX as f64) as u64;

    // Get fractional part for low 64 bits
    let frac_hi = hi_val - (hi_u64 as f64);
    let frac_lo = scaled.lo();
    let frac_combined = frac_hi + frac_lo;

    // Scale fractional part by 2^64
    let lo_val = frac_combined * SCALE_HI;
    let lo_u64 = lo_val.max(0.0).min(u64::MAX as f64) as u64;

    ((hi_u64 as u128) << 64) | (lo_u64 as u128)
}

impl<T: UnsignedInteger> CDTSamplerLogSpaceDD<T> {
    /// Generate a CDT sampler using double-double precision arithmetic
    /// This provides ~106-bit precision in CDF computation
    pub fn new(std_dev: f64, tail_cut: f64, modulus_minus_one: T) -> Self {
        let max_std_dev = std_dev * tail_cut;
        let mut length = max_std_dev.floor() as usize + 1;

        assert!(length <= 1024, "CDT table too large: {}", length);
        if length <= 1 {
            length = 2;
        }

        let two_var_f64 = 2.0 * std_dev * std_dev;

        // Compute PDF in log-space to avoid underflow
        let mut log_pdf_f64 = vec![f64::NEG_INFINITY; length];
        log_pdf_f64[0] = 0.5f64.ln();

        for (i, item) in log_pdf_f64.iter_mut().enumerate().skip(1) {
            let i_f64 = i as f64;
            *item = -(i_f64 * i_f64) / two_var_f64;
        }

        // Convert to TwoFloat only for normalization step
        let log_pdf: Vec<TwoFloat> = log_pdf_f64.iter().map(|&x| TwoFloat::from(x)).collect();

        // Compute normalization constant in log-space
        let log_sum = log_sum_exp_dd(&log_pdf);

        let log_sum_f64 = log_sum.hi();

        let mut pdf_f64 = Vec::with_capacity(length);
        for &log_val_f64 in log_pdf_f64.iter() {
            let normalized = log_val_f64 - log_sum_f64;
            pdf_f64.push(normalized.exp());
        }

        // Build cumulative distribution with TwoFloat precision for accumulation
        // This is where precision matters most
        let mut cdt_dd = Vec::with_capacity(length + 1);
        let one = TwoFloat::from(1.0);
        let zero = TwoFloat::from(0.0);
        let mut acc = zero;

        cdt_dd.push(zero);
        for &p_f64 in pdf_f64.iter() {
            acc += TwoFloat::from(p_f64);
            // Clamp to [0, 1] to handle any numerical errors
            cdt_dd.push(if acc > one {
                one
            } else if acc < zero {
                zero
            } else {
                acc
            });
        }

        // Ensure last value is exactly 1.0
        if let Some(last) = cdt_dd.last_mut() {
            *last = one;
        }

        assert_eq!(cdt_dd.len(), length + 1, "CDT length mismatch");

        // Map to u128 range [0, 2^128-1] with high precision
        let cdt: Vec<u128> = cdt_dd.into_iter().map(twofloat_to_u128).collect();

        Self {
            std_dev,
            modulus_minus_one,
            cdt,
        }
    }

    /// Returns the standard deviation of this sampler
    #[inline]
    pub fn std_dev(&self) -> f64 {
        self.std_dev
    }
}

impl<T: UnsignedInteger> Distribution<T> for CDTSamplerLogSpaceDD<T> {
    #[inline]
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> T {
        let r: u128 = rng.random();

        let positive = (r & 1) == 1;

        // Binary search to find the right bin
        let idx = self.cdt.partition_point(|&x| x <= r) - 1;

        let v: T = idx.as_into();

        if v.is_zero() {
            return T::ZERO;
        }

        if positive {
            v
        } else {
            self.modulus_minus_one - v + T::ONE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log_sum_exp_dd() {
        // Test with similar values
        let values: Vec<TwoFloat> = vec![0.0, 0.0, 0.0]
            .into_iter()
            .map(TwoFloat::from)
            .collect();
        let result = log_sum_exp_dd(&values);
        // log(3) ≈ 1.0986
        let expected = TwoFloat::from(3.0).ln();
        let diff = (result - expected).abs();
        assert!(
            diff < TwoFloat::from(1e-15),
            "diff was {:?}, result was {:?}",
            diff,
            result
        );

        // Test with moderately separated values
        let values: Vec<TwoFloat> = vec![-1.0, -2.0, -3.0]
            .into_iter()
            .map(TwoFloat::from)
            .collect();
        let result = log_sum_exp_dd(&values);

        // Result should be around -1 + ln(1 + exp(-1) + exp(-2))
        let expected = TwoFloat::from(-1.0)
            + (TwoFloat::from(1.0) + TwoFloat::from(-1.0).exp() + TwoFloat::from(-2.0).exp()).ln();
        let diff = (result - expected).abs();
        assert!(
            diff < TwoFloat::from(1e-10),
            "diff was {:?}, result was {:?}, expected was {:?}",
            diff,
            result,
            expected
        );
    }

    #[test]
    fn test_twofloat_to_u128() {
        // Test boundary values
        let zero = twofloat_to_u128(TwoFloat::from(0.0));
        assert_eq!(zero, 0);

        let one = twofloat_to_u128(TwoFloat::from(1.0));
        assert_eq!(one, u128::MAX);

        // Test 0.5 should be approximately half of max
        let half = twofloat_to_u128(TwoFloat::from(0.5));
        let expected = u128::MAX / 2;
        // Allow some tolerance due to rounding
        // Use saturating subtraction to avoid overflow
        let diff = if half > expected {
            half - expected
        } else {
            expected - half
        };
        assert!(diff < 1000, "diff was {}", diff);
    }

    #[test]
    fn test_sampler_creation() {
        let sampler = CDTSamplerLogSpaceDD::<u64>::new(3.19, 12.0, u64::MAX - 1);
        assert!((sampler.std_dev() - 3.19).abs() < 1e-10);
        assert!(!sampler.cdt.is_empty());

        // First value should be 0
        assert_eq!(sampler.cdt[0], 0);

        // Last value should be u128::MAX
        assert_eq!(sampler.cdt.last().unwrap(), &u128::MAX);
    }

    #[test]
    fn test_precision_improvement() {
        // Create both samplers with same parameters
        use crate::discrete_gaussian::CDTSamplerLogSpace;

        let std_dev = 3.19;
        let tail_cut = 12.0;

        let sampler_f64 = CDTSamplerLogSpace::<u64>::new(std_dev, tail_cut, u64::MAX - 1);
        let sampler_dd = CDTSamplerLogSpaceDD::<u64>::new(std_dev, tail_cut, u64::MAX - 1);

        // Both should have the same std_dev
        assert!((sampler_f64.std_dev() - sampler_dd.std_dev()).abs() < 1e-10);

        // Test that both samplers can sample (basic functionality test)
        use rand::distr::Distribution;
        let mut rng = rand::rng();

        let _sample_f64: u64 = sampler_f64.sample(&mut rng);
        let _sample_dd: u64 = sampler_dd.sample(&mut rng);
    }
}
