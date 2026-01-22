/// Pre-computes a Barrett factor with which modular multiplication can
/// be performed more efficiently
#[derive(Clone, Copy)]
pub struct MultiplyFactor {
    operand: u64,
    barrett_factor: u64,
}

impl MultiplyFactor {
    /// Computes and stores the Barrett factor `floor((operand << bit_shift) / modulus)`.
    ///
    /// This is useful when modular multiplications of the form `(x * operand) mod modulus`
    /// are performed many times with the same `modulus` and `operand`.
    ///
    /// Note: passing `operand = 1` can be used to precompute a Barrett factor for
    /// multiplications of the form `(x * y) mod modulus`, where only the `modulus`
    /// is reused across calls.
    pub fn new(operand: u64, bit_shift: u32, modulus: u64) -> Self {
        assert!(
            operand <= modulus,
            "operand {operand} must be less than modulus {modulus}",
        );
        assert!(
            bit_shift == 32 || bit_shift == 52 || bit_shift == 64,
            "Unsupported BitShift {bit_shift}"
        );

        let op_hi = operand >> (64u32 - bit_shift);
        let op_lo = if bit_shift == 64 {
            0
        } else {
            operand << bit_shift
        };

        let barrett_factor = divide_u128_u64_lo(op_hi, op_lo, modulus);

        Self {
            operand,
            barrett_factor,
        }
    }

    pub fn operand(&self) -> u64 {
        self.operand
    }

    pub fn barrett_factor(&self) -> u64 {
        self.barrett_factor
    }

    #[inline]
    pub fn lazy_mul_modulo<const BIT_SHIFT: u32>(self, b: u64, modulus: u64) -> u64 {
        let hw = if BIT_SHIFT == 32 {
            (self.barrett_factor * b) >> BIT_SHIFT
        } else {
            ((self.barrett_factor as u128 * b as u128) >> BIT_SHIFT) as u64
        };
        self.operand
            .wrapping_mul(b)
            .wrapping_sub(modulus.wrapping_mul(hw))
    }

    #[inline]
    pub fn mul_modulo<const BIT_SHIFT: u32>(self, b: u64, modulus: u64) -> u64 {
        let r = self.lazy_mul_modulo::<BIT_SHIFT>(b, modulus);
        if r >= modulus { r - modulus } else { r }
    }
}

// Returns low 64bit of 128b/64b where x1=high 64b, x0=low 64b
fn divide_u128_u64_lo(x1: u64, x0: u64, y: u64) -> u64 {
    let n = ((x1 as u128) << 64) | (x0 as u128);
    let q = n / y as u128;
    q as u64
}

#[cfg(test)]
mod tests {
    use rand::{Rng, distr::Uniform};

    use super::*;

    #[test]
    fn test_32() {
        let mut rng = rand::rng();

        let q = 536813569;
        let distr = Uniform::new(0, q).unwrap();

        for _ in 0..1000 {
            let a = rng.sample(distr);
            let b = rng.sample(distr);

            let mf_a = MultiplyFactor::new(a, 32, q);
            let c = mf_a.mul_modulo::<32>(b, q);

            assert_eq!(c, (a * b) % q);
        }
    }

    #[test]
    fn test_52() {
        let mut rng = rand::rng();

        let q = 562949953392641;
        let distr = Uniform::new(0, q).unwrap();

        for _ in 0..1000 {
            let a = rng.sample(distr);
            let b = rng.sample(distr);

            let mf_a = MultiplyFactor::new(a, 52, q);
            let c = mf_a.mul_modulo::<52>(b, q);

            assert_eq!(c, ((a as u128 * b as u128) % q as u128) as u64);
        }
    }

    #[test]
    fn test_64() {
        let mut rng = rand::rng();

        let q = 1152921504606830593;
        let distr = Uniform::new(0, q).unwrap();

        for _ in 0..1000 {
            let a = rng.sample(distr);
            let b = rng.sample(distr);

            let mf_a = MultiplyFactor::new(a, 64, q);
            let c = mf_a.mul_modulo::<64>(b, q);

            assert_eq!(c, ((a as u128 * b as u128) % q as u128) as u64);
        }
    }
}
