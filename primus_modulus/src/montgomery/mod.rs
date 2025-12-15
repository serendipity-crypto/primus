use core::fmt::Display;

use primus_integer::DivRemScalar;

use crate::integer::UnsignedInteger;

mod ops;

/// A modulus using Montgomery reduction algorithm.
///
/// Montgomery reduction is efficient for performing many modular multiplications
/// with the same modulus. It transforms numbers into Montgomery form where
/// multiplication followed by reduction can be done more efficiently.
///
/// For a modulus N and R = 2^T::BITS, a number x is represented in Montgomery
/// form as x' = x * R mod N. Montgomery reduction efficiently computes
/// x * y * R^(-1) mod N from x' and y'.
#[derive(Debug, Clone, Copy)]
pub struct MontgomeryModulus<T: UnsignedInteger> {
    /// The modulus value (must be odd)
    value: T,
    /// R mod N
    r: T,
    /// R^2 mod N, used to convert numbers into Montgomery form
    r2: T,
    /// N' = -N^(-1) mod R, used in Montgomery reduction
    n_prime: T,
}

impl<T: UnsignedInteger> MontgomeryModulus<T> {
    /// Creates a new [`MontgomeryModulus<T>`] with the given value.
    ///
    /// # Panics
    ///
    /// Panics if:
    /// - `value` is 0 or 1
    /// - `value` is even (Montgomery reduction requires odd modulus)
    pub fn new(value: T) -> Self {
        assert!(value > T::ONE, "modulus can't be 0 or 1.");
        assert!(
            value & T::ONE == T::ONE,
            "Montgomery reduction requires odd modulus."
        );

        // Compute N' = -N^(-1) mod R
        // We need n_prime such that N * n_prime ≡ -1 (mod R)
        // Since R = 2^BITS, we compute this iteratively
        let n_prime = compute_n_prime(value);

        // Compute R^2 mod N
        // R = 2^BITS, so R mod N = 2^BITS mod N
        // Then R^2 mod N = (R mod N)^2 mod N
        let (r, r2) = compute_r2(value);

        Self {
            value,
            r,
            r2,
            n_prime,
        }
    }

    /// Returns the value of this [`MontgomeryModulus<T>`].
    #[inline]
    pub const fn value(&self) -> T {
        self.value
    }

    /// Returns R mod N.
    #[inline]
    pub const fn r(&self) -> T {
        self.r
    }

    /// Returns R^2 mod N.
    #[inline]
    pub const fn r2(&self) -> T {
        self.r2
    }

    /// Returns N' = -N^(-1) mod R.
    #[inline]
    pub const fn n_prime(&self) -> T {
        self.n_prime
    }

    /// Converts a value to Montgomery form.
    ///
    /// Computes `value * R mod N` where R = 2^T::BITS.
    #[inline]
    pub fn to_montgomery(&self, value: T) -> T {
        // To convert to Montgomery form: x' = x * R mod N
        // We compute this as REDC(x * R^2 mod N)
        // Since we have R^2 mod N precomputed, we can use montgomery_mul
        let (lo, hi) = value.widening_mul(self.r2);
        self.montgomery_reduce([lo, hi])
    }

    /// Converts a value from Montgomery form.
    ///
    /// Computes `value * R^(-1) mod N` where R = 2^T::BITS.
    #[inline]
    pub fn from_montgomery(&self, value: T) -> T {
        let m = value.wrapping_mul(self.n_prime);

        let (m_lo, m_hi) = m.widening_mul(self.value);

        let (_, carry1) = value.overflowing_add(m_lo);

        let result = m_hi.wrapping_add(T::as_from(carry1));

        if result >= self.value {
            result.wrapping_sub(self.value)
        } else {
            result
        }
    }

    /// Montgomery reduction (REDC algorithm).
    ///
    /// Given T = [t_hi, t_lo] representing t_hi * 2^BITS + t_lo,
    /// computes T * R^(-1) mod N efficiently.
    #[inline]
    pub fn montgomery_reduce(&self, value: [T; 2]) -> T {
        // REDC algorithm:
        // 1. m = (T mod R) * N' mod R = t_lo * N' mod R
        // 2. t = (T + m * N) / R
        // 3. if t >= N then return t - N else return t

        let m = value[0].wrapping_mul(self.n_prime);

        // Compute (T + m * N) / R
        // T = value[1] * R + value[0]
        // m * N gives us a double-word result
        let (m_lo, m_hi) = m.widening_mul(self.value);

        // Add m * N to T
        // (value[0] + m_lo) may overflow, producing a carry
        let (_sum_lo, carry1) = value[0].overflowing_add(m_lo);

        // Add the high word and the carry
        let (mut result, carry2) = value[1].overflowing_add(m_hi);
        if carry1 {
            let (r, c) = result.overflowing_add(T::ONE);
            result = r;
            if c {
                // This means we had an overflow
                // After division by R, we need to handle this
                // But with proper m, sum_lo should be 0, so we just take result
                return result.wrapping_sub(self.value);
            }
        }

        // Division by R means taking the high word (result)
        // Reduce if necessary
        if carry2 || result >= self.value {
            result.wrapping_sub(self.value)
        } else {
            result
        }
    }
}

/// Computes N' = -N^(-1) mod R where R = 2^T::BITS
fn compute_n_prime<T: UnsignedInteger>(n: T) -> T {
    // We use Newton's method: x_{i+1} = x_i * (2 - N * x_i)
    // Starting with x_0 = 1 (works for odd N)
    let mut x = T::ONE;

    // Iterate log2(BITS) times
    let mut i = 1;
    while i < T::BITS {
        x = x.wrapping_mul(T::TWO.wrapping_sub(n.wrapping_mul(x)));
        i *= 2;
    }

    // We want -N^(-1), so negate the result
    x.wrapping_neg()
}

fn compute_r2<T: UnsignedInteger>(n: T) -> (T, T) {
    let mut quotient = [T::ZERO; 2];

    let r = DivRemScalar::div_rem_scalar(&[T::ZERO, T::ONE], n, &mut quotient);

    let (lo, hi) = r.widening_mul(r);

    let r2 = DivRemScalar::div_rem_scalar(&[lo, hi], n, &mut quotient);

    (r, r2)
}

impl<T: UnsignedInteger> Display for MontgomeryModulus<T> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl<T: UnsignedInteger> primus_reduce::Modulus for MontgomeryModulus<T> {
    type ValueT = T;

    #[inline]
    fn value(self) -> Option<Self::ValueT> {
        Some(self.value)
    }

    #[inline(always)]
    fn value_unchecked(self) -> Self::ValueT {
        self.value
    }

    #[inline(always)]
    fn minus_one(self) -> Self::ValueT {
        self.value - T::ONE
    }
}

#[cfg(test)]
mod tests {
    use primus_reduce::ops::*;
    use rand::distr::{Distribution, Uniform};

    use super::*;

    #[test]
    fn test_to_from_montgomery() {
        let mut m = rand::random_range(3..u32::MAX);
        if m & 1 == 0 {
            m |= 1;
        }
        let modulus = MontgomeryModulus::<u32>::new(m);
        let value = rand::random_range(0..m);

        let mont_form = modulus.to_montgomery(value);
        let back = modulus.from_montgomery(mont_form);

        assert_eq!(back, value);
    }

    #[test]
    fn test_n_prime_property() {
        let mut m = rand::random_range(3..u32::MAX);
        if m & 1 == 0 {
            m |= 1;
        }

        let modulus = MontgomeryModulus::<u32>::new(m);

        // N * N' ≡ -1 (mod R)
        // Since R = 2^32 for u32, this wraps around
        let product = modulus.value().wrapping_mul(modulus.n_prime());
        assert_eq!(product, u32::MAX); // -1 in two's complement
    }

    #[test]
    fn test_ops() {
        let mut rng = rand::rng();

        let mut m = rand::random_range(3..u32::MAX);
        if m & 1 == 0 {
            m |= 1;
        }

        let modulus = MontgomeryModulus::<u32>::new(m);

        let distr = Uniform::new(0, m).unwrap();

        let a = distr.sample(&mut rng);
        let b = distr.sample(&mut rng);
        let c = distr.sample(&mut rng);

        let a_m = modulus.to_montgomery(a);
        let b_m = modulus.to_montgomery(b);
        let c_m = modulus.to_montgomery(c);

        if m < u32::MAX >> 1 {
            assert_eq!(modulus.from_montgomery(a_m.wrapping_add(m)), a);
        }

        assert_eq!(
            modulus.from_montgomery(modulus.reduce_add(a_m, b_m)),
            ((a as u64 + b as u64) % m as u64) as u32
        );

        assert_eq!(
            modulus.from_montgomery(modulus.reduce_sub(a_m, b_m)),
            ((m as u64 + a as u64 - b as u64) % m as u64) as u32
        );

        let p_m = modulus.reduce_mul(a_m, b_m);
        assert_eq!(
            modulus.from_montgomery(p_m),
            ((a as u64 * b as u64) % m as u64) as u32
        );

        let p_m = modulus.reduce_mul_add(a_m, b_m, c_m);
        assert_eq!(
            modulus.from_montgomery(p_m),
            ((a as u64 * b as u64 + c as u64) % m as u64) as u32
        );
    }
}
