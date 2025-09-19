use integer::UnsignedInteger;

mod ops;

/// Power of 2 modulus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(transparent)]
pub struct PowOf2Modulus<T: UnsignedInteger> {
    /// The special value for performing `reduce`.
    ///
    /// It's equal to modulus value sub one.
    mask: T,
}

impl<T: UnsignedInteger> PowOf2Modulus<T> {
    /// Creates a [`PowOf2Modulus<T>`].
    ///
    /// - `value`: The value of the modulus.
    #[inline]
    pub fn new(value: T) -> Self {
        assert!(
            value > T::ONE && value.is_power_of_two(),
            "The value is not a power of 2."
        );
        Self {
            mask: value - T::ONE,
        }
    }

    /// Creates a [`PowOf2Modulus<T>`].
    ///
    /// - `mask`: modulus value minus one.
    #[inline]
    pub fn with_mask(mask: T) -> Self {
        assert!(mask.count_zeros() == mask.leading_zeros() && !mask.is_zero());
        assert!(
            mask.leading_zeros() > 0,
            "NativeModulus<T> supports modulus value such as 2⁸, 2¹⁶, 2³², 2⁶⁴, 2¹²⁸"
        );
        Self { mask }
    }

    /// Returns the value of this [`PowOf2Modulus<T>`].
    #[inline]
    pub fn value(self) -> T {
        self.mask + T::ONE
    }

    /// Returns the mask of this [`PowOf2Modulus<T>`],
    /// which is equal to modulus value sub one.
    #[inline]
    pub const fn mask(self) -> T {
        self.mask
    }
}

impl<T: UnsignedInteger> reduce::Modulus for PowOf2Modulus<T> {
    type ValueT = T;

    #[inline]
    fn value(self) -> Option<Self::ValueT> {
        Some(self.mask + T::ONE)
    }

    #[inline(always)]
    fn minus_one(self) -> Self::ValueT {
        self.mask
    }
}

#[cfg(test)]
mod tests {
    use rand::{distr::Uniform, prelude::*, rng};

    use reduce::ops::*;

    use super::*;

    #[test]
    fn test_modulus_create() {
        let mut rng = rng();

        let _m = <PowOf2Modulus<u8>>::new(rng.random_range(2..=(u8::MAX >> 2)).next_power_of_two());
        let _m =
            <PowOf2Modulus<u16>>::new(rng.random_range(2..=(u16::MAX >> 2)).next_power_of_two());
        let _m =
            <PowOf2Modulus<u32>>::new(rng.random_range(2..=(u32::MAX >> 2)).next_power_of_two());
        let _m =
            <PowOf2Modulus<u64>>::new(rng.random_range(2..=(u64::MAX >> 2)).next_power_of_two());
        let _m =
            <PowOf2Modulus<u128>>::new(rng.random_range(2..=(u128::MAX >> 2)).next_power_of_two());

        let _m = <PowOf2Modulus<u8>>::with_mask(
            rng.random_range(2..=(u8::MAX >> 2)).next_power_of_two() - 1,
        );
        let _m = <PowOf2Modulus<u16>>::with_mask(
            rng.random_range(2..=(u16::MAX >> 2)).next_power_of_two() - 1,
        );
        let _m = <PowOf2Modulus<u32>>::with_mask(
            rng.random_range(2..=(u32::MAX >> 2)).next_power_of_two() - 1,
        );
        let _m = <PowOf2Modulus<u64>>::with_mask(
            rng.random_range(2..=(u64::MAX >> 2)).next_power_of_two() - 1,
        );
        let _m = <PowOf2Modulus<u128>>::with_mask(
            rng.random_range(2..=(u128::MAX >> 2)).next_power_of_two() - 1,
        );
    }

    #[test]
    #[should_panic]
    fn test_modulus_create_panic() {
        let mut rng = rng();

        let m = loop {
            let r = rng.random_range(0..=(u64::MAX >> 2));
            if !r.is_power_of_two() {
                break r;
            }
        };

        let _m = PowOf2Modulus::<u64>::new(m);
    }

    #[test]
    fn test_reduce() {
        let mut rng = rng();

        let m: u64 = rng.random_range(2..=(u64::MAX >> 2)).next_power_of_two();
        let modulus = PowOf2Modulus::<u64>::new(m);
        let dis = Uniform::new_inclusive(0, modulus.mask()).unwrap();

        let v: u64 = rng.sample(dis);
        assert_eq!(modulus.reduce(v), v % m);

        let a: u64 = rng.sample(dis);
        let b: u64 = rng.sample(dis);
        assert_eq!(modulus.reduce_add(a, b), (a + b) % m);

        let a: u64 = rng.sample(dis);
        let b: u64 = rng.sample(dis);
        assert_eq!(modulus.reduce_sub(a, b), (m + a - b) % m);

        let a: u64 = rng.sample(dis);
        let b: u64 = rng.sample(dis);
        assert_eq!(
            modulus.reduce_mul(a, b),
            ((a as u128 * b as u128) % m as u128) as u64
        );

        let a: u64 = rng.sample(dis);
        let a_neg = modulus.reduce_neg(a);
        assert_eq!(modulus.reduce_add(a, a_neg), 0);

        assert_eq!(modulus.reduce_neg(0), 0);
    }
}
