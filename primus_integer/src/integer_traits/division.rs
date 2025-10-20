use num_traits::Zero;

const FAST_DIV_WIDE: bool = cfg!(any(target_arch = "x86", target_arch = "x86_64"));

/// A trait for types that support combined division and remainder operation.
pub trait DivRem: Sized {
    fn div_rem(self, divisor: Self) -> (Self, Self);
}

macro_rules! impl_div_rem {
    ($($T:ty)*) => {$(
        impl DivRem for $T {
            #[inline]
            fn div_rem(self, rhs: Self) -> (Self, Self) {
                (self / rhs, self % rhs)
            }
        })*
    };
}

impl_div_rem! { u8 u16 u32 u64 u128 usize }

/// A trait for array types that support division and remainder operation by a scalar.
pub trait DivRemScalar: Sized {
    fn div_rem_scalar(dividend: &[Self], divisor: Self, quotient: &mut [Self]) -> Self;
}

macro_rules! impl_div_rem_scalar {
    ($T:ty, $W:ty, $HALF_BITS:ident, $LO_MASK:ident, $div_half:ident, $div_wide:ident) => {
        const $HALF_BITS: u32 = <$T>::BITS >> 1;
        const $LO_MASK: $T = <$T>::MAX >> $HALF_BITS;

        /// Perform (rem * 2^{<$T>::BITS} + digit) / divisor, returns the quotient and the remainder.
        ///
        /// # Correctness
        ///
        /// * rem < divisor
        /// * divisor < $LO_MASK
        #[inline]
        fn $div_half(rem: $T, digit: $T, divisor: $T) -> ($T, $T) {
            debug_assert!(rem < divisor && divisor <= $LO_MASK);
            let (hi, rem) = ((rem << $HALF_BITS) | (digit >> $HALF_BITS)).div_rem(divisor);
            let (lo, rem) = ((rem << $HALF_BITS) | (digit & $LO_MASK)).div_rem(divisor);
            ((hi << $HALF_BITS) | lo, rem)
        }

        /// Perform (hi * 2^{<$T>::BITS} + lo) / divisor, returns the quotient and the remainder.
        ///
        /// # Correctness
        ///
        /// * hi < divisor
        #[inline]
        fn $div_wide(hi: $T, lo: $T, divisor: $T) -> ($T, $T) {
            debug_assert!(hi < divisor);

            let lhs = lo as $W | ((hi as $W) << <$T>::BITS);
            let rhs = divisor as $W;
            ((lhs / rhs) as $T, (lhs % rhs) as $T)
        }

        impl DivRemScalar for $T {
            fn div_rem_scalar(dividend: &[Self], divisor: Self, quotient: &mut [Self]) -> $T {
                debug_assert_eq!(dividend.len(), quotient.len());
                if divisor.is_zero() {
                    panic!("attempt to divide by zero")
                }

                if divisor == 1 {
                    quotient.copy_from_slice(dividend);
                    return 0;
                }

                quotient.fill(0);

                // if dividend.iter().all(|&v| v == 0) {
                //     return 0;
                // }

                // if dividend.len() == 1 {
                //     let (q, r) = dividend[0].div_rem(divisor);
                //     quotient[0] = q;
                //     return r;
                // }

                if dividend.iter().skip(1).all(|&v| v.is_zero()) {
                    if dividend[0] < divisor {
                        return dividend[0];
                    } else if dividend[0] == divisor {
                        quotient[0] = 1;
                        return 0;
                    } else {
                        let (q, r) = dividend[0].div_rem(divisor);
                        quotient[0] = q;
                        return r;
                    }
                }

                let mut rem = 0;

                if !FAST_DIV_WIDE && divisor <= $LO_MASK {
                    for (&d_elem, q_elem) in dividend.iter().rev().zip(quotient.iter_mut().rev()) {
                        let (q, r) = $div_half(rem, d_elem, divisor);
                        *q_elem = q;
                        rem = r;
                    }
                } else {
                    for (&d_elem, q_elem) in dividend.iter().rev().zip(quotient.iter_mut().rev()) {
                        let (q, r) = $div_wide(rem, d_elem, divisor);
                        *q_elem = q;
                        rem = r;
                    }
                }

                rem
            }
        }
    };
}

impl_div_rem_scalar!(u8, u16, HALF_BITS_U8, LO_MASK_U8, div_half_u8, div_wide_u8);
impl_div_rem_scalar!(
    u16,
    u32,
    HALF_BITS_U16,
    LO_MASK_U16,
    div_half_u16,
    div_wide_u16
);
impl_div_rem_scalar!(
    u32,
    u64,
    HALF_BITS_U32,
    LO_MASK_U32,
    div_half_u32,
    div_wide_u32
);
impl_div_rem_scalar!(
    u64,
    u128,
    HALF_BITS_U64,
    LO_MASK_U64,
    div_half_u64,
    div_wide_u64
);

#[cfg(target_pointer_width = "64")]
impl_div_rem_scalar!(
    usize,
    u128,
    HALF_BITS_USIZE,
    LO_MASK_USIZE,
    div_half_usize,
    div_wide_usize
);

#[cfg(target_pointer_width = "32")]
impl_div_rem_scalar!(
    usize,
    u64,
    HALF_BITS_USIZE,
    LO_MASK_USIZE,
    div_half_usize,
    div_wide_usize
);

mod big_digit {
    use std::num::NonZeroU128;

    use super::{DivRem, DivRemScalar};

    const HALF_BITS_U128: u32 = 64;
    const LO_MASK_U128: u128 = u64::MAX as u128;

    impl DivRemScalar for u128 {
        fn div_rem_scalar(dividend: &[u128], divisor: u128, quotient: &mut [u128]) -> u128 {
            debug_assert!(divisor != 0);
            debug_assert_eq!(dividend.len(), quotient.len());

            let len = dividend.len();

            if len == 1 {
                quotient[0] = dividend[0] / divisor;
                dividend[0] % divisor
            } else if divisor <= LO_MASK_U128 {
                let mut rem = 0;
                for (&d_elem, q_elem) in dividend.iter().rev().zip(quotient.iter_mut().rev()) {
                    let (q, r) = div_half_u128(rem, d_elem, divisor);
                    *q_elem = q;
                    rem = r;
                }
                rem
            } else {
                let mut remainder = 0;
                let non_zero_divisor = unsafe { NonZeroU128::new_unchecked(divisor) };

                if dividend[len - 1] < divisor {
                    // The result fits in 128 bits.
                    quotient[len - 1] = 0;
                    quotient[len - 2] = udiv256_by_128_to_128(
                        dividend[len - 1],
                        dividend[len - 2],
                        non_zero_divisor,
                        &mut remainder,
                    );
                } else {
                    // First, divide with the high part to get the remainder in dividend.s.high.
                    // After that dividend.s.high < divisor.s.low.
                    let (q, r) = dividend[len - 1].div_rem(divisor);
                    quotient[len - 1] = q;
                    quotient[len - 2] = udiv256_by_128_to_128(
                        r,
                        dividend[len - 2],
                        non_zero_divisor,
                        &mut remainder,
                    );
                }

                for i in (0..len - 2).rev() {
                    quotient[i] = udiv256_by_128_to_128(
                        remainder,
                        dividend[i],
                        non_zero_divisor,
                        &mut remainder,
                    );
                }

                remainder
            }
        }
    }

    /// Perform (rem * 2¹²⁸ + digit) / divisor, returns the quotient and the remainder.
    ///
    /// # Correctness
    ///
    /// * rem < divisor
    /// * divisor < 2⁶⁴
    #[inline]
    fn div_half_u128(rem: u128, digit: u128, divisor: u128) -> (u128, u128) {
        debug_assert!(rem < divisor && divisor <= LO_MASK_U128);
        let (hi, rem) = ((rem << HALF_BITS_U128) | (digit >> HALF_BITS_U128)).div_rem(divisor);
        let (lo, rem) = ((rem << HALF_BITS_U128) | (digit & LO_MASK_U128)).div_rem(divisor);
        ((hi << HALF_BITS_U128) | lo, rem)
    }

    #[inline(always)]
    fn udiv256_by_128_to_128(u1: u128, u0: u128, mut v: NonZeroU128, r: &mut u128) -> u128 {
        const N_UDWORD_BITS: u32 = 128;

        #[inline]
        unsafe fn shl_nz(x: NonZeroU128, n: u32) -> NonZeroU128 {
            debug_assert!(n < N_UDWORD_BITS);
            let res: u128 = x.get() << n;
            debug_assert_ne!(res, 0);
            unsafe { NonZeroU128::new_unchecked(res) }
        }

        #[inline]
        unsafe fn shr_nz(x: NonZeroU128, n: u32) -> NonZeroU128 {
            debug_assert!(n < N_UDWORD_BITS);
            let res: u128 = x.get() >> n;
            debug_assert_ne!(res, 0);
            unsafe { NonZeroU128::new_unchecked(res) }
        }

        const B: u128 = 1 << (N_UDWORD_BITS / 2); // Number base (128 bits)
        let (un1, un0): (u128, u128); // Norm. dividend LSD's
        let (vn1, vn0): (NonZeroU128, u128); // Norm. divisor digits
        let (mut q1, mut q0): (u128, u128); // Quotient digits
        let (un128, un21, un10): (u128, u128, u128); // Dividend digit pairs

        debug_assert!(v.get() > u1);

        let s = v.leading_zeros();
        debug_assert_ne!(s, N_UDWORD_BITS);
        if s > 0 {
            // Normalize the divisor.
            v = unsafe { shl_nz(v, s) };
            un128 = (u1 << s) | (u0 >> (N_UDWORD_BITS - s));
            un10 = u0 << s; // Shift dividend left
        } else {
            // Avoid undefined behavior of (u0 >> 128).
            un128 = u1;
            un10 = u0;
        }

        // Break divisor up into two 64-bit digits.
        vn1 = unsafe { shr_nz(v, N_UDWORD_BITS / 2) };
        vn0 = v.get() & 0xFFFF_FFFF_FFFF_FFFF;

        // Break right half of dividend into two digits.
        un1 = un10 >> (N_UDWORD_BITS / 2);
        un0 = un10 & 0xFFFF_FFFF_FFFF_FFFF;

        // Compute the first quotient digit, q1.
        q1 = un128 / vn1;
        let mut rhat = un128 - q1 * vn1.get();

        // q1 has at most error 2. No more than 2 iterations.
        while q1 >= B || q1 * vn0 > B * rhat + un1 {
            q1 -= 1;
            rhat += vn1.get();
            if rhat >= B {
                break;
            }
        }

        un21 = un128
            .wrapping_mul(B)
            .wrapping_add(un1)
            .wrapping_sub(q1.wrapping_mul(v.get()));

        // Compute the second quotient digit.
        q0 = un21 / vn1;
        rhat = un21 - q0 * vn1.get();

        // q0 has at most error 2. No more than 2 iterations.
        while q0 >= B || q0 * vn0 > B * rhat + un0 {
            q0 -= 1;
            rhat += vn1.get();
            if rhat >= B {
                break;
            }
        }

        *r = (un21
            .wrapping_mul(B)
            .wrapping_add(un0)
            .wrapping_sub(q0.wrapping_mul(v.get())))
            >> s;
        q1 * B + q0
    }
}
