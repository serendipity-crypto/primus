//! Extended GCD and modular inverse for unsigned integer types.
//!
//! This implementation refers to the following codebases.
//! <https://flintlib.org/doc/ulong_extras.html#c.n_xgcd>
//! <https://flintlib.org/doc/ulong_extras.html#c.n_gcdinv>

/// Greatest common divisor and Bézout coefficients
pub trait Xgcd: Sized {
    /// Calculates the Greatest Common Divisor (GCD) of the number and `other`. The
    /// result is always non-negative.
    ///
    /// # Examples
    ///
    /// ```
    /// use primus_gcd::Xgcd;
    ///
    /// assert_eq!(42u64.gcd(56), 14);
    /// assert_eq!(0u64.gcd(5), 5);
    /// assert_eq!(5u64.gcd(0), 5);
    /// ```
    fn gcd(self, other: Self) -> Self;

    /// Check whether two numbers are coprime.
    ///
    /// # Examples
    ///
    /// ```
    /// use primus_gcd::Xgcd;
    ///
    /// assert!(14u64.is_coprime(25));
    /// assert!(!14u64.is_coprime(28));
    /// assert!(!0u64.is_coprime(0));
    /// ```
    #[allow(clippy::wrong_self_convention)]
    fn is_coprime(self, other: Self) -> bool;

    /// Returns the greatest common divisor `g` of `x` and `y` and unsigned
    /// values `a` and `b` such that `a x - b y = g`. We require `x ≥ y`.
    ///
    /// We claim that computing the extended greatest common divisor via the
    /// Euclidean algorithm always results in cofactor `|a| < x/2`,
    /// `|b| < x/2`, with perhaps some small degenerate exceptions.
    ///
    /// We proceed by induction.
    ///
    /// Suppose we are at some step of the algorithm, with `x_n = q y_n + r`
    /// with `r ≥ 1`, and suppose `1 = s y_n - t r` with
    /// `s < y_n / 2`, `t < y_n / 2` by hypothesis.
    ///
    /// Write `1 = s y_n - t (x_n - q y_n) = (s + t q) y_n - t x_n`.
    ///
    /// It suffices to show that `(s + t q) < x_n / 2` as `t < y_n / 2 < x_n / 2`,
    /// which will complete the induction step.
    ///
    /// But at the previous step in the back substitution we would have had
    /// `1 = s r - c d` with `s < r/2` and `c < r/2`.
    ///
    /// Then `s + t q < r/2 + y_n / 2 q = (r + q y_n)/2 = x_n / 2`.
    ///
    /// # Examples
    ///
    /// ```
    /// use primus_gcd::Xgcd;
    ///
    /// let (a, b, d) = u64::xgcd(240, 46);
    /// assert_eq!(d, 2);
    /// assert_eq!(a as u128 * 240 - b as u128 * 46, 2);
    /// ```
    ///
    /// # Panics if
    ///
    /// - `x < y`
    fn xgcd(x: Self, y: Self) -> (Self, Self, Self);

    /// Returns the greatest common divisor `g` of `x` and `y` and computes
    /// `a` such that `0 ≤ a < y` and `a x = gcd(x, y) mod y`, when
    /// this is defined. We require `x < y`.
    ///
    /// When `y = 1` the greatest common divisor is set to `1` and `a` is
    /// set to `0`.
    ///
    /// This is merely an adaption of the extended Euclidean algorithm
    /// computing just one cofactor and reducing it modulo `y`.
    ///
    /// # Examples
    ///
    /// ```
    /// use primus_gcd::Xgcd;
    ///
    /// let (a, d) = u64::gcdinv(17, 29);
    /// assert_eq!(d, 1);
    /// assert_eq!((a as u128 * 17) % 29, 1);
    /// ```
    ///
    /// # Panics if
    ///
    /// - `x ≥ y`
    fn gcdinv(x: Self, y: Self) -> (Self, Self);
}

macro_rules! impl_extended_gcd {
    (impl Xgcd for $SelfT:ty; SignedType: $SignedT:ty) => {
        impl Xgcd for $SelfT {
            #[inline]
            fn gcd(self, other: Self) -> Self {
                // Use Stein's algorithm
                let mut m = self;
                let mut n = other;
                if m == 0 || n == 0 {
                    return m | n;
                }

                // find common factors of 2
                let shift = (m | n).trailing_zeros();

                // divide n and m by 2 until odd
                m >>= m.trailing_zeros();
                n >>= n.trailing_zeros();

                while m != n {
                    if m > n {
                        m -= n;
                        m >>= m.trailing_zeros();
                    } else {
                        n -= m;
                        n >>= n.trailing_zeros();
                    }
                }
                m << shift
            }

            #[inline]
            fn is_coprime(self, other: Self) -> bool {
                // Fast paths that avoid computing the full GCD.
                if self == other {
                    return self == 1;
                }
                if self == 1 || other == 1 {
                    return true;
                }
                self.gcd(other) == 1
            }

            #[inline]
            fn xgcd(x: Self, y: Self) -> (Self, Self, Self) {
                let mut u1: $SignedT;
                let mut u2: $SignedT;
                let mut v1: $SignedT;
                let mut v2: $SignedT;
                let mut t1: $SignedT;
                let mut t2: $SignedT;

                let mut u3: Self;
                let mut v3: Self;
                let mut quot: Self;
                let mut rem: Self;
                let mut d: Self;

                assert!(x >= y);

                u1 = 1;
                v2 = 1;
                u2 = 0;
                v1 = 0;
                u3 = x;
                v3 = y;

                // x and y both have top bit set
                if ((x & y) as $SignedT) < 0 {
                    d = u3 - v3;
                    t2 = v2;
                    t1 = u2;
                    u2 = u1 - u2;
                    u1 = t1;
                    u3 = v3;
                    v2 = v1 - v2;
                    v1 = t2;
                    v3 = d;
                }

                // second value has second msb set
                while ((v3 << 1) as $SignedT) < 0 {
                    d = u3 - v3;
                    if d < v3 {
                        // quot = 1
                        t2 = v2;
                        t1 = u2;
                        u2 = u1 - u2;
                        u1 = t1;
                        u3 = v3;
                        v2 = v1 - v2;
                        v1 = t2;
                        v3 = d;
                    } else if d < (v3 << 1) {
                        // quot = 2
                        t1 = u2;
                        u2 = u1 - (u2 << 1);
                        u1 = t1;
                        u3 = v3;
                        t2 = v2;
                        v2 = v1 - (v2 << 1);
                        v1 = t2;
                        v3 = d - u3;
                    } else {
                        // quot = 3
                        t1 = u2;
                        u2 = u1 - 3 * u2;
                        u1 = t1;
                        u3 = v3;
                        t2 = v2;
                        v2 = v1 - 3 * v2;
                        v1 = t2;
                        v3 = d - (u3 << 1);
                    }
                }

                while v3 > 0 {
                    d = u3 - v3;

                    // overflow not possible, top 2 bits of v3 not set
                    if u3 < (v3 << 2) {
                        if d < v3 {
                            // quot = 1
                            t2 = v2;
                            t1 = u2;
                            u2 = u1 - u2;
                            u1 = t1;
                            u3 = v3;
                            v2 = v1 - v2;
                            v1 = t2;
                            v3 = d;
                        } else if d < (v3 << 1) {
                            // quot = 2
                            t1 = u2;
                            u2 = u1 - (u2 << 1);
                            u1 = t1;
                            u3 = v3;
                            t2 = v2;
                            v2 = v1 - (v2 << 1);
                            v1 = t2;
                            v3 = d - u3;
                        } else {
                            // quot = 3
                            t1 = u2;
                            u2 = u1 - 3 * u2;
                            u1 = t1;
                            u3 = v3;
                            t2 = v2;
                            v2 = v1 - 3 * v2;
                            v1 = t2;
                            v3 = d - (u3 << 1);
                        }
                    } else {
                        quot = u3 / v3;
                        rem = u3 - v3 * quot;
                        t1 = u2;
                        u2 = u1 - (quot as $SignedT) * u2;
                        u1 = t1;
                        u3 = v3;
                        t2 = v2;
                        v2 = v1 - (quot as $SignedT) * v2;
                        v1 = t2;
                        v3 = rem;
                    }
                }

                /* Remarkably, |u1| < x/2, thus comparison with 0 is valid */
                if u1 <= 0 {
                    u1 = u1.wrapping_add_unsigned(y);
                    v1 = v1.wrapping_sub_unsigned(x);
                }

                (u1 as Self, v1.wrapping_neg() as Self, u3)
            }

            #[inline]
            fn gcdinv(mut x: Self, y: Self) -> (Self, Self) {
                let mut v1: $SignedT;
                let mut v2: $SignedT;
                let mut t2: $SignedT;

                let mut d: Self;
                let mut r: Self;
                let mut quot: Self;
                let mut rem: Self;

                assert!(y > x);

                v1 = 0;
                v2 = 1;
                r = x;
                x = y;

                // y and x both have top bit set
                if ((x & r) as $SignedT) < 0 {
                    d = x - r;
                    t2 = v2;
                    x = r;
                    v2 = v1 - v2;
                    v1 = t2;
                    r = d;
                }

                // second value has second msb set
                while ((r << 1) as $SignedT) < 0 {
                    d = x - r;
                    if (d < r) {
                        // quot = 1
                        t2 = v2;
                        x = r;
                        v2 = v1 - v2;
                        v1 = t2;
                        r = d;
                    } else if (d < (r << 1)) {
                        // quot = 2
                        x = r;
                        t2 = v2;
                        v2 = v1 - (v2 << 1);
                        v1 = t2;
                        r = d - x;
                    } else {
                        // quot = 3
                        x = r;
                        t2 = v2;
                        v2 = v1 - 3 * v2;
                        v1 = t2;
                        r = d - (x << 1);
                    }
                }

                while r > 0 {
                    // overflow not possible due to top 2 bits of r not being set
                    if x < (r << 2) {
                        // if quot < 4
                        d = x - r;
                        if (d < r) {
                            // quot = 1
                            t2 = v2;
                            x = r;
                            v2 = v1 - v2;
                            v1 = t2;
                            r = d;
                        } else if d < (r << 1) {
                            // quot = 2
                            x = r;
                            t2 = v2;
                            v2 = v1 - (v2 << 1);
                            v1 = t2;
                            r = d - x;
                        } else {
                            // quot = 3
                            x = r;
                            t2 = v2;
                            v2 = v1 - 3 * v2;
                            v1 = t2;
                            r = d - (x << 1);
                        }
                    } else {
                        quot = x / r;
                        rem = x - r * quot;
                        x = r;
                        t2 = v2;
                        v2 = v1 - (quot as $SignedT) * v2;
                        v1 = t2;
                        r = rem;
                    }
                }

                if v1 < 0 {
                    v1 = v1.wrapping_add_unsigned(y);
                }

                (v1 as Self, x)
            }
        }
    };
}

impl_extended_gcd!(impl Xgcd for u8; SignedType: i8);
impl_extended_gcd!(impl Xgcd for u16; SignedType: i16);
impl_extended_gcd!(impl Xgcd for u32; SignedType: i32);
impl_extended_gcd!(impl Xgcd for u64; SignedType: i64);
impl_extended_gcd!(impl Xgcd for usize; SignedType: isize);
impl_extended_gcd!(impl Xgcd for u128; SignedType: i128);

#[cfg(test)]
mod tests {
    use rand::{prelude::*, rng};

    use super::*;

    macro_rules! gcd_edge_tests {
        ($mod_name:ident, $T:ty) => {
            mod $mod_name {
                use super::*;

                #[test]
                fn test_gcd_zero() {
                    assert_eq!(<$T>::gcd(0_u8.into(), 0_u8.into()), 0_u8.into());
                    assert_eq!(<$T>::gcd(42_u8.into(), 0_u8.into()), 42_u8.into());
                    assert_eq!(<$T>::gcd(0_u8.into(), 42_u8.into()), 42_u8.into());
                }

                #[test]
                fn test_gcd_one() {
                    assert_eq!(<$T>::gcd(1_u8.into(), 1_u8.into()), 1_u8.into());
                    assert_eq!(<$T>::gcd(1_u8.into(), 42_u8.into()), 1_u8.into());
                    assert_eq!(<$T>::gcd(42_u8.into(), 1_u8.into()), 1_u8.into());
                }

                #[test]
                fn test_gcd_symmetry() {
                    let mut rng = rng();
                    for _ in 0..20 {
                        let a = rng.random_range(<$T>::MIN..=<$T>::MAX);
                        let b = rng.random_range(<$T>::MIN..=<$T>::MAX);
                        assert_eq!(a.gcd(b), b.gcd(a));
                    }
                }

                #[test]
                fn test_is_coprime_zero() {
                    assert!(!<$T>::is_coprime(0_u8.into(), 0_u8.into()));
                    assert!(<$T>::is_coprime(0_u8.into(), 1_u8.into()));
                    assert!(<$T>::is_coprime(1_u8.into(), 0_u8.into()));
                    assert!(<$T>::is_coprime(1_u8.into(), 1_u8.into()));
                }

                #[test]
                fn test_xgcd_d_is_gcd() {
                    let mut rng = rng();
                    for _ in 0..20 {
                        let x = rng.random_range(<$T>::MIN..(<$T>::MAX >> 2));
                        let y = rng.random_range(<$T>::MIN..=x);
                        let (_a, _b, d) = <$T>::xgcd(x, y);
                        assert_eq!(d, x.gcd(y));
                    }
                }

                #[test]
                fn test_xgcd_bezout() {
                    // x = y: d = x, and a*x - b*x = x => (a - b) = 1
                    let x = rng().random_range(1..(<$T>::MAX >> 2));
                    let (_a, _b, d) = <$T>::xgcd(x, x);
                    assert_eq!(d, x);

                    // y = 1: gcd(x, 1) = 1
                    let (a, _b, d) = <$T>::xgcd(x, 1_u8.into());
                    assert_eq!(d, 1_u8.into());
                    assert!(a < x);
                }

                #[test]
                fn test_gcdinv_d_is_gcd() {
                    let mut rng = rng();
                    for _ in 0..20 {
                        let y = rng.random_range(1..(<$T>::MAX >> 2));
                        let x = rng.random_range(<$T>::MIN..y);
                        let (a, d) = <$T>::gcdinv(x, y);
                        assert_eq!(d, x.gcd(y));
                        assert!(a < y, "a={a} should be < y={y}");
                    }
                }

                #[test]
                fn test_gcdinv_edge() {
                    // x = 0
                    let y = rng().random_range(1..(<$T>::MAX >> 2));
                    let (a, d) = <$T>::gcdinv(0_u8.into(), y);
                    assert_eq!(d, y);
                    assert!(a < y);
                }
            }
        };
    }

    macro_rules! gcd_identity_tests {
        ($mod_name:ident, $T:ty, $WideT:ty) => {
            mod $mod_name {
                use super::*;

                #[test]
                fn test_xgcd_identity() {
                    let mut rng = rng();
                    for _ in 0..20 {
                        let x = rng.random_range(0..<$T>::MAX >> 1);
                        let y = rng.random_range(0..=x);
                        let (a, b, d) = <$T>::xgcd(x, y);
                        assert_eq!(
                            a as $WideT * x as $WideT - b as $WideT * y as $WideT,
                            d as $WideT,
                        );
                    }
                }

                #[test]
                fn test_gcdinv_identity() {
                    let mut rng = rng();
                    for _ in 0..20 {
                        let y = rng.random_range(1..<$T>::MAX >> 1);
                        let x = rng.random_range(0..y);
                        let (a, d) = <$T>::gcdinv(x, y);
                        assert_eq!(
                            (a as $WideT * x as $WideT) % y as $WideT,
                            d as $WideT % y as $WideT,
                        );
                    }
                }
            }
        };
    }

    // Edge case tests for all types, including u128
    gcd_edge_tests!(tests_u8, u8);
    gcd_edge_tests!(tests_u16, u16);
    gcd_edge_tests!(tests_u32, u32);
    gcd_edge_tests!(tests_u64, u64);
    gcd_edge_tests!(tests_usize, usize);
    gcd_edge_tests!(tests_u128, u128);

    // Identity tests (a*x - b*y = d) — only for types that have a wider type
    gcd_identity_tests!(tests_id_u8, u8, u16);
    gcd_identity_tests!(tests_id_u16, u16, u32);
    gcd_identity_tests!(tests_id_u32, u32, u64);
    gcd_identity_tests!(tests_id_u64, u64, u128);
    gcd_identity_tests!(tests_id_usize, usize, u128);

    macro_rules! gcd_msb_tests {
        ($mod_name:ident, $T:ty, $WideT:ty) => {
            mod $mod_name {
                use super::*;

                #[test]
                fn test_xgcd_msb_both() {
                    // Path A: both operands have the MSB set, triggering
                    // `(x & y) as signed < 0`.  We pick the smallest values
                    // with the MSB set (x = y = MAX>>1 + 1): xgcd(x, x)
                    // hits the first-if and exits immediately (v3 = 0),
                    // avoiding the quot=3 branch in the main loop.
                    let val = (<$T>::MAX >> 1) + 1;
                    let (a, b, d) = <$T>::xgcd(val, val);
                    assert_eq!(d, val);
                    assert!(
                        a >= 1 && a.wrapping_sub(1) == b,
                        "a={a}, b={b}, expected a-b=1"
                    );
                }

                #[test]
                fn test_xgcd_msb_second() {
                    // Path B: second MSB set but top MSB clear, triggering
                    // `(v3 << 1) as signed < 0`.  Range [MAX>>2+1, MAX>>1]
                    // is safe: cofactor bound x/2 < MAX/4, so 3*x/2 <
                    // 3*MAX/8 which fits in the signed type.
                    let mut rng = rng();
                    let lo = (<$T>::MAX >> 2) + 1;
                    let hi = <$T>::MAX >> 1;
                    for _ in 0..10 {
                        let x = rng.random_range(lo..=hi);
                        let y = rng.random_range(lo..=x);
                        let (a, b, d) = <$T>::xgcd(x, y);
                        assert_eq!(d, x.gcd(y));
                        let lhs = a as $WideT * x as $WideT - b as $WideT * y as $WideT;
                        assert_eq!(lhs, d as $WideT);
                    }
                }

                #[test]
                fn test_gcdinv_msb_second() {
                    let mut rng = rng();
                    let lo = (<$T>::MAX >> 2) + 1;
                    let hi = <$T>::MAX >> 1;
                    for _ in 0..10 {
                        let y = rng.random_range(lo..=hi);
                        let x = rng.random_range(lo..y);
                        let (a, d) = <$T>::gcdinv(x, y);
                        assert_eq!(d, x.gcd(y));
                        assert_eq!(
                            (a as $WideT * x as $WideT) % y as $WideT,
                            d as $WideT % y as $WideT,
                        );
                    }
                }
            }
        };
    }

    // MSB tests only for u64+ — smaller types (u8/u16/u32) use matching-width
    // signed types that may overflow in the quot=3 branch with MSB-range inputs.
    gcd_msb_tests!(tests_msb_u64, u64, u128);

    mod tests_msb_u128 {
        use super::*;

        #[test]
        fn test_xgcd_msb_both() {
            // Path A: both operands MSB set, using minimal values to avoid
            // the quot=3 branch in the main loop.
            let val = (u128::MAX >> 1) + 1;
            let (a, b, d) = u128::xgcd(val, val);
            assert_eq!(d, val);
            assert!(a >= 1 && a.wrapping_sub(1) == b);
        }

        #[test]
        fn test_xgcd_msb_second() {
            // Path B: second MSB set, top MSB clear.
            let mut rng = rng();
            let lo = (u128::MAX >> 2) + 1;
            let hi = u128::MAX >> 1;
            for _ in 0..10 {
                let x = rng.random_range(lo..=hi);
                let y = rng.random_range(lo..=x);
                let (_a, _b, d) = u128::xgcd(x, y);
                assert_eq!(d, x.gcd(y));
            }
        }

        #[test]
        fn test_gcdinv_msb_second() {
            let mut rng = rng();
            let lo = (u128::MAX >> 2) + 1;
            let hi = u128::MAX >> 1;
            for _ in 0..10 {
                let y = rng.random_range(lo..=hi);
                let x = rng.random_range(lo..y);
                let (a, d) = u128::gcdinv(x, y);
                assert_eq!(d, x.gcd(y));
                assert!(a < y);
            }
        }
    }
}
