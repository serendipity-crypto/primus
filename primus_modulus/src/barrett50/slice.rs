//! Slice trait wiring for [`Barrett50Modulus`].
//!
//! The mul-family slice traits dispatch to the IFMA kernels in
//! [`super::simd_ifma`] when the target feature combo
//! `avx512f + avx512dq + avx512ifma` is enabled at compile time;
//! otherwise (and for the non-mul traits in all configurations) they
//! delegate to the wrapped [`crate::BarrettModulus<u64>`], whose own
//! slice impls already SIMD-vectorize via portable_simd or fall back to
//! scalar.

use primus_reduce::{
    LazyReduceMulAddSlice, LazyReduceMulSlice, ReduceAddSlice, ReduceDotProduct, ReduceMulAddSlice,
    ReduceMulSlice, ReduceNegSlice, ReduceOnceSlice, ReduceSubSlice,
};

use super::Barrett50Modulus;

// ---------------------------------------------------------------------------
// Non-mul slice traits: always delegate to inner. IFMA gives no speedup on
// add/sub/neg/once, so don't duplicate kernels.
// ---------------------------------------------------------------------------

impl ReduceOnceSlice<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_once_slice_assign(self, values: &mut [u64]) {
        self.inner.reduce_once_slice_assign(values)
    }
    #[inline]
    fn reduce_once_slice_to(self, input: &[u64], output: &mut [u64]) {
        self.inner.reduce_once_slice_to(input, output)
    }
}

impl ReduceNegSlice<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_neg_slice_assign(self, values: &mut [u64]) {
        self.inner.reduce_neg_slice_assign(values)
    }
    #[inline]
    fn reduce_neg_slice_to(self, input: &[u64], output: &mut [u64]) {
        self.inner.reduce_neg_slice_to(input, output)
    }
}

impl ReduceAddSlice<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_add_slice_assign(self, a: &mut [u64], b: &[u64]) {
        self.inner.reduce_add_slice_assign(a, b)
    }
    #[inline]
    fn reduce_add_slice_to(self, a: &[u64], b: &[u64], output: &mut [u64]) {
        self.inner.reduce_add_slice_to(a, b, output)
    }
}

impl ReduceSubSlice<u64> for Barrett50Modulus {
    #[inline]
    fn reduce_sub_slice_assign(self, a: &mut [u64], b: &[u64]) {
        self.inner.reduce_sub_slice_assign(a, b)
    }
    #[inline]
    fn reduce_sub_slice_to(self, a: &[u64], b: &[u64], output: &mut [u64]) {
        self.inner.reduce_sub_slice_to(a, b, output)
    }
}

// ---------------------------------------------------------------------------
// Mul-family slice traits: IFMA fast path vs fallback.
//
// We use two mutually exclusive `#[cfg(...)]`-gated impl modules. Rust does
// not allow macros inside `#[cfg(...)]`, so the gate condition is repeated
// verbatim on both branches.
// ---------------------------------------------------------------------------

#[cfg(all(
    feature = "nightly",
    feature = "simd",
    target_feature = "avx512f",
    target_feature = "avx512dq",
    target_feature = "avx512ifma",
))]
mod ifma_impls {
    use super::*;

    impl ReduceMulSlice<u64> for Barrett50Modulus {
        #[inline]
        fn reduce_mul_slice_assign(self, a: &mut [u64], b: &[u64]) {
            unsafe { super::super::simd_ifma::reduce_mul_slice_assign(self, a, b) }
        }
        #[inline]
        fn reduce_mul_slice_to(self, a: &[u64], b: &[u64], output: &mut [u64]) {
            unsafe { super::super::simd_ifma::reduce_mul_slice_to(self, a, b, output) }
        }
    }

    impl LazyReduceMulSlice<u64> for Barrett50Modulus {
        #[inline]
        fn lazy_reduce_mul_slice_assign(self, a: &mut [u64], b: &[u64]) {
            unsafe { super::super::simd_ifma::lazy_reduce_mul_slice_assign(self, a, b) }
        }
        #[inline]
        fn lazy_reduce_mul_slice_to(self, a: &[u64], b: &[u64], output: &mut [u64]) {
            unsafe { super::super::simd_ifma::lazy_reduce_mul_slice_to(self, a, b, output) }
        }
    }

    impl ReduceMulAddSlice<u64> for Barrett50Modulus {
        #[inline]
        fn reduce_add_mul_slice_assign(self, acc: &mut [u64], a: &[u64], b: &[u64]) {
            unsafe { super::super::simd_ifma::reduce_add_mul_slice_assign(self, acc, a, b) }
        }
        #[inline]
        fn reduce_sub_mul_slice_assign(self, acc: &mut [u64], a: &[u64], b: &[u64]) {
            unsafe { super::super::simd_ifma::reduce_sub_mul_slice_assign(self, acc, a, b) }
        }
        #[inline]
        fn reduce_mul_add_slice_to(self, a: &[u64], b: &[u64], c: &[u64], output: &mut [u64]) {
            unsafe { super::super::simd_ifma::reduce_mul_add_slice_to(self, a, b, c, output) }
        }
        #[inline]
        fn reduce_scalar_mul_add_slice_to(
            self,
            scalar: u64,
            b: &[u64],
            c: &[u64],
            output: &mut [u64],
        ) {
            unsafe {
                super::super::simd_ifma::reduce_scalar_mul_add_slice_to(self, scalar, b, c, output)
            }
        }
    }

    impl LazyReduceMulAddSlice<u64> for Barrett50Modulus {
        #[inline]
        fn lazy_reduce_add_mul_slice_assign(self, acc: &mut [u64], a: &[u64], b: &[u64]) {
            unsafe { super::super::simd_ifma::lazy_reduce_add_mul_slice_assign(self, acc, a, b) }
        }
        #[inline]
        fn lazy_reduce_sub_mul_slice_assign(self, acc: &mut [u64], a: &[u64], b: &[u64]) {
            unsafe { super::super::simd_ifma::lazy_reduce_sub_mul_slice_assign(self, acc, a, b) }
        }
        #[inline]
        fn lazy_reduce_mul_add_slice_to(self, a: &[u64], b: &[u64], c: &[u64], output: &mut [u64]) {
            unsafe { super::super::simd_ifma::lazy_reduce_mul_add_slice_to(self, a, b, c, output) }
        }
        #[inline]
        fn lazy_reduce_scalar_mul_add_slice_to(
            self,
            scalar: u64,
            b: &[u64],
            c: &[u64],
            output: &mut [u64],
        ) {
            unsafe {
                super::super::simd_ifma::lazy_reduce_scalar_mul_add_slice_to(
                    self, scalar, b, c, output,
                )
            }
        }
    }

    impl ReduceDotProduct<u64> for Barrett50Modulus {
        type Output = u64;

        #[inline]
        fn reduce_dot_product(self, a: impl AsRef<[u64]>, b: impl AsRef<[u64]>) -> u64 {
            unsafe { super::super::simd_ifma::reduce_dot_product(self, a.as_ref(), b.as_ref()) }
        }

        #[inline]
        fn reduce_dot_product_iter(
            self,
            a: impl IntoIterator<Item = u64>,
            b: impl IntoIterator<Item = u64>,
        ) -> u64 {
            self.inner.reduce_dot_product_iter(a, b)
        }
    }
}

#[cfg(not(all(
    feature = "nightly",
    feature = "simd",
    target_feature = "avx512f",
    target_feature = "avx512dq",
    target_feature = "avx512ifma",
)))]
mod fallback_impls {
    use super::*;

    impl ReduceMulSlice<u64> for Barrett50Modulus {
        #[inline]
        fn reduce_mul_slice_assign(self, a: &mut [u64], b: &[u64]) {
            self.inner.reduce_mul_slice_assign(a, b)
        }
        #[inline]
        fn reduce_mul_slice_to(self, a: &[u64], b: &[u64], output: &mut [u64]) {
            self.inner.reduce_mul_slice_to(a, b, output)
        }
    }

    impl LazyReduceMulSlice<u64> for Barrett50Modulus {
        #[inline]
        fn lazy_reduce_mul_slice_assign(self, a: &mut [u64], b: &[u64]) {
            self.inner.lazy_reduce_mul_slice_assign(a, b)
        }
        #[inline]
        fn lazy_reduce_mul_slice_to(self, a: &[u64], b: &[u64], output: &mut [u64]) {
            self.inner.lazy_reduce_mul_slice_to(a, b, output)
        }
    }

    impl ReduceMulAddSlice<u64> for Barrett50Modulus {
        #[inline]
        fn reduce_add_mul_slice_assign(self, acc: &mut [u64], a: &[u64], b: &[u64]) {
            self.inner.reduce_add_mul_slice_assign(acc, a, b)
        }
        #[inline]
        fn reduce_sub_mul_slice_assign(self, acc: &mut [u64], a: &[u64], b: &[u64]) {
            self.inner.reduce_sub_mul_slice_assign(acc, a, b)
        }
        #[inline]
        fn reduce_mul_add_slice_to(self, a: &[u64], b: &[u64], c: &[u64], output: &mut [u64]) {
            self.inner.reduce_mul_add_slice_to(a, b, c, output)
        }
        #[inline]
        fn reduce_scalar_mul_add_slice_to(
            self,
            scalar: u64,
            b: &[u64],
            c: &[u64],
            output: &mut [u64],
        ) {
            self.inner
                .reduce_scalar_mul_add_slice_to(scalar, b, c, output)
        }
    }

    impl LazyReduceMulAddSlice<u64> for Barrett50Modulus {
        #[inline]
        fn lazy_reduce_add_mul_slice_assign(self, acc: &mut [u64], a: &[u64], b: &[u64]) {
            self.inner.lazy_reduce_add_mul_slice_assign(acc, a, b)
        }
        #[inline]
        fn lazy_reduce_sub_mul_slice_assign(self, acc: &mut [u64], a: &[u64], b: &[u64]) {
            self.inner.lazy_reduce_sub_mul_slice_assign(acc, a, b)
        }
        #[inline]
        fn lazy_reduce_mul_add_slice_to(self, a: &[u64], b: &[u64], c: &[u64], output: &mut [u64]) {
            self.inner.lazy_reduce_mul_add_slice_to(a, b, c, output)
        }
        #[inline]
        fn lazy_reduce_scalar_mul_add_slice_to(
            self,
            scalar: u64,
            b: &[u64],
            c: &[u64],
            output: &mut [u64],
        ) {
            self.inner
                .lazy_reduce_scalar_mul_add_slice_to(scalar, b, c, output)
        }
    }

    impl ReduceDotProduct<u64> for Barrett50Modulus {
        type Output = u64;

        #[inline]
        fn reduce_dot_product(self, a: impl AsRef<[u64]>, b: impl AsRef<[u64]>) -> u64 {
            self.inner.reduce_dot_product(a, b)
        }

        #[inline]
        fn reduce_dot_product_iter(
            self,
            a: impl IntoIterator<Item = u64>,
            b: impl IntoIterator<Item = u64>,
        ) -> u64 {
            self.inner.reduce_dot_product_iter(a, b)
        }
    }
}

#[cfg(test)]
mod tests {
    use rand::distr::{Distribution, Uniform};

    use super::*;
    use crate::BarrettModulus;

    /// FHE-style 50-bit prime: 2^50 − 27.
    const MODULUS: u64 = (1u64 << 50) - 27;

    fn rand_slice(len: usize) -> Vec<u64> {
        let mut rng = rand::rng();
        let distr = Uniform::new(0, MODULUS).unwrap();
        distr.sample_iter(&mut rng).take(len).collect()
    }

    fn mul_mod(a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % MODULUS as u128) as u64
    }

    fn mmod(x: u64) -> u64 {
        x % MODULUS
    }

    fn reference() -> BarrettModulus<u64> {
        BarrettModulus::<u64>::new(MODULUS)
    }

    fn test_lengths() -> &'static [usize] {
        &[
            0usize, 1, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 1000,
        ]
    }

    #[test]
    fn reduce_mul_slice_to_matches_reference() {
        let m50 = Barrett50Modulus::new(MODULUS);
        for &len in test_lengths() {
            let a = rand_slice(len);
            let b = rand_slice(len);
            let mut got = vec![0u64; len];
            let mut want = vec![0u64; len];
            m50.reduce_mul_slice_to(&a, &b, &mut got);
            reference().reduce_mul_slice_to(&a, &b, &mut want);
            assert_eq!(got, want, "len={len}");
            for (&g, (&x, &y)) in got.iter().zip(a.iter().zip(b.iter())) {
                assert_eq!(g, mul_mod(x, y), "len={len}, x={x}, y={y}");
            }
        }
    }

    #[test]
    fn reduce_mul_slice_assign_matches_reference() {
        let m50 = Barrett50Modulus::new(MODULUS);
        for &len in test_lengths() {
            let a = rand_slice(len);
            let b = rand_slice(len);
            let mut got = a.clone();
            let mut want = a.clone();
            m50.reduce_mul_slice_assign(&mut got, &b);
            reference().reduce_mul_slice_assign(&mut want, &b);
            assert_eq!(got, want, "len={len}");
        }
    }

    #[test]
    fn lazy_reduce_mul_slice_to_canonical_matches() {
        let m50 = Barrett50Modulus::new(MODULUS);
        for &len in test_lengths() {
            let a = rand_slice(len);
            let b = rand_slice(len);
            let mut got = vec![0u64; len];
            m50.lazy_reduce_mul_slice_to(&a, &b, &mut got);
            for (i, (&g, (&x, &y))) in got.iter().zip(a.iter().zip(b.iter())).enumerate() {
                assert!(g < 2 * MODULUS, "lazy output {g} >= 2m at len={len}, i={i}");
                let canonical = if g >= MODULUS { g - MODULUS } else { g };
                assert_eq!(canonical, mul_mod(x, y), "len={len}, i={i}");
            }
        }
    }

    #[test]
    fn reduce_mul_add_slice_to_matches_reference() {
        let m50 = Barrett50Modulus::new(MODULUS);
        for &len in test_lengths() {
            let a = rand_slice(len);
            let b = rand_slice(len);
            let c = rand_slice(len);
            let mut got = vec![0u64; len];
            m50.reduce_mul_add_slice_to(&a, &b, &c, &mut got);
            for (i, (&g, ((&x, &y), &z))) in got
                .iter()
                .zip(a.iter().zip(b.iter()).zip(c.iter()))
                .enumerate()
            {
                assert_eq!(g, mmod(mul_mod(x, y) + z), "len={len}, i={i}");
            }
        }
    }

    #[test]
    fn lazy_reduce_mul_add_slice_to_canonical_matches() {
        let m50 = Barrett50Modulus::new(MODULUS);
        for &len in test_lengths() {
            let a = rand_slice(len);
            let b = rand_slice(len);
            let c = rand_slice(len);
            let mut got = vec![0u64; len];
            m50.lazy_reduce_mul_add_slice_to(&a, &b, &c, &mut got);
            for (i, (&g, ((&x, &y), &z))) in got
                .iter()
                .zip(a.iter().zip(b.iter()).zip(c.iter()))
                .enumerate()
            {
                assert!(
                    g < 2 * MODULUS,
                    "lazy mul_add {g} >= 2m at len={len}, i={i}"
                );
                let canonical = if g >= MODULUS { g - MODULUS } else { g };
                assert_eq!(canonical, mmod(mul_mod(x, y) + z), "len={len}, i={i}");
            }
        }
    }

    #[test]
    fn reduce_add_mul_slice_assign_matches_reference() {
        let m50 = Barrett50Modulus::new(MODULUS);
        for &len in test_lengths() {
            let acc0 = rand_slice(len);
            let a = rand_slice(len);
            let b = rand_slice(len);
            let mut got = acc0.clone();
            let mut want = acc0.clone();
            m50.reduce_add_mul_slice_assign(&mut got, &a, &b);
            reference().reduce_add_mul_slice_assign(&mut want, &a, &b);
            assert_eq!(got, want, "len={len}");
        }
    }

    #[test]
    fn reduce_sub_mul_slice_assign_matches_reference() {
        let m50 = Barrett50Modulus::new(MODULUS);
        for &len in test_lengths() {
            let acc0 = rand_slice(len);
            let a = rand_slice(len);
            let b = rand_slice(len);
            let mut got = acc0.clone();
            let mut want = acc0.clone();
            m50.reduce_sub_mul_slice_assign(&mut got, &a, &b);
            reference().reduce_sub_mul_slice_assign(&mut want, &a, &b);
            assert_eq!(got, want, "len={len}");
        }
    }

    #[test]
    fn reduce_scalar_mul_add_slice_to_matches_reference() {
        let m50 = Barrett50Modulus::new(MODULUS);
        let mut rng = rand::rng();
        let scalar = Uniform::new(0, MODULUS).unwrap().sample(&mut rng);
        for &len in test_lengths() {
            let b = rand_slice(len);
            let c = rand_slice(len);
            let mut got = vec![0u64; len];
            let mut want = vec![0u64; len];
            m50.reduce_scalar_mul_add_slice_to(scalar, &b, &c, &mut got);
            reference().reduce_scalar_mul_add_slice_to(scalar, &b, &c, &mut want);
            assert_eq!(got, want, "len={len}");
        }
    }

    #[test]
    fn reduce_dot_product_matches_reference() {
        let m50 = Barrett50Modulus::new(MODULUS);
        for &len in test_lengths() {
            let a = rand_slice(len);
            let b = rand_slice(len);
            let got = m50.reduce_dot_product(&a, &b);
            let want = reference().reduce_dot_product(&a, &b);
            assert_eq!(got, want, "len={len}");
        }
    }

    /// Exercise boundary moduli on the most demanding path (mul_slice).
    #[test]
    fn boundary_moduli() {
        for &m in &[
            1u64 << 48,
            (1u64 << 48) + 1,
            (1u64 << 49) + 17,
            (1u64 << 50) - 27,
            (1u64 << 50) - 9,
            (1u64 << 50) - 1,
        ] {
            let m50 = Barrett50Modulus::new(m);
            let bref = BarrettModulus::<u64>::new(m);
            let mut rng = rand::rng();
            let distr = Uniform::new(0, m).unwrap();
            let a: Vec<u64> = (0..256).map(|_| distr.sample(&mut rng)).collect();
            let b: Vec<u64> = (0..256).map(|_| distr.sample(&mut rng)).collect();
            // Worst-case inputs: m−1.
            let mut a2 = a.clone();
            let mut b2 = b.clone();
            a2[0] = m - 1;
            b2[0] = m - 1;
            let mut got = vec![0u64; 256];
            let mut want = vec![0u64; 256];
            m50.reduce_mul_slice_to(&a2, &b2, &mut got);
            bref.reduce_mul_slice_to(&a2, &b2, &mut want);
            assert_eq!(got, want, "modulus={m}");
        }
    }
}
