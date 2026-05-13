//! Slice-level (bulk) modular operations.
//!
//! These traits mirror the scalar `Reduce*` traits in [`crate::ops`] but
//! operate on whole slices, so that implementations can dispatch to a
//! SIMD kernel internally and amortize the per-call overhead.
//!
//! Each trait bundles an in-place (`*_assign`) form and an out-of-place
//! (`*_to`) form. There are no default impls: every modulus type provides
//! its own body, which is typically a thin wrapper around a hand-written
//! scalar / SIMD kernel.
//!
//! # Length and value-range invariants
//!
//! All slice traits use `debug_assert*!` to check length agreement and
//! value-range pre-conditions. In release builds the checks are stripped;
//! callers (typically the polynomial / NTT layer) are expected to uphold
//! them at higher-level boundaries.

/// Slice form of [`crate::ReduceOnce`].
pub trait ReduceOnceSlice<T> {
    /// For each `v` in `values`: `v -= modulus` if `v >= modulus`, where
    /// `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - Each `values[i] < 2 * modulus`
    /// - Each result is `< modulus`
    fn reduce_once_slice_assign(self, values: &mut [T]);

    /// For each `v` in `input`: writes `v - modulus` if `v >= modulus`,
    /// otherwise `v`, into `output`, where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `input.len() == output.len()`
    /// - Each `input[i] < 2 * modulus`
    /// - Each result is `< modulus`
    fn reduce_once_slice_to(self, input: &[T], output: &mut [T]);
}

/// Slice form of [`crate::ReduceNeg`].
pub trait ReduceNegSlice<T> {
    /// Calculates `v = -v (mod modulus)` for each element in-place, where
    /// `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - Each `values[i] < modulus`
    fn reduce_neg_slice_assign(self, values: &mut [T]);

    /// Writes `-input[i] (mod modulus)` into `output[i]` for each element,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `input.len() == output.len()`
    /// - Each `input[i] < modulus`
    fn reduce_neg_slice_to(self, input: &[T], output: &mut [T]);
}

/// Slice form of [`crate::ReduceAdd`].
pub trait ReduceAddSlice<T, B = T> {
    /// Calculates `a[i] = (a[i] + b[i]) (mod modulus)` element-wise,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len()`
    /// - Each `a[i] < modulus` and `b[i] < modulus`
    fn reduce_add_slice_assign(self, a: &mut [T], b: &[B]);

    /// Writes `a[i] + b[i] (mod modulus)` into `output[i]` element-wise,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len() == output.len()`
    /// - Each `a[i] < modulus` and `b[i] < modulus`
    fn reduce_add_slice_to(self, a: &[T], b: &[B], output: &mut [T]);
}

/// Slice form of [`crate::ReduceSub`].
pub trait ReduceSubSlice<T, B = T> {
    /// Calculates `a[i] = (a[i] - b[i]) (mod modulus)` element-wise,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len()`
    /// - Each `a[i] < modulus` and `b[i] < modulus`
    fn reduce_sub_slice_assign(self, a: &mut [T], b: &[B]);

    /// Writes `a[i] - b[i] (mod modulus)` into `output[i]` element-wise,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len() == output.len()`
    /// - Each `a[i] < modulus` and `b[i] < modulus`
    fn reduce_sub_slice_to(self, a: &[T], b: &[B], output: &mut [T]);
}

/// Slice form of [`crate::ReduceMul`].
pub trait ReduceMulSlice<T, B = T> {
    /// Calculates `a[i] = (a[i] * b[i]) (mod modulus)` element-wise,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len()`
    /// - Each `a[i] * b[i] < modulus²`
    fn reduce_mul_slice_assign(self, a: &mut [T], b: &[B]);

    /// Writes `a[i] * b[i] (mod modulus)` into `output[i]` element-wise,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len() == output.len()`
    /// - Each `a[i] * b[i] < modulus²`
    fn reduce_mul_slice_to(self, a: &[T], b: &[B], output: &mut [T]);
}

/// Slice form of [`crate::ReduceMulAdd`].
///
/// Provides the four fused multiply-add shapes that the polynomial /
/// NTT layer needs:
///
/// 1. `acc[i] += a[i] * b[i]`              — FMAC accumulate
/// 2. `acc[i] -= a[i] * b[i]`              — fused multiply-subtract
/// 3. `out[i]  = a[i] * b[i] + c[i]`       — three-input one-output
/// 4. `out[i]  = scalar * b[i] + c[i]`     — scalar × slice plus addend
pub trait ReduceMulAddSlice<T> {
    /// Calculates `acc[i] = (acc[i] + a[i] * b[i]) (mod modulus)`
    /// element-wise, where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `acc.len() == a.len() == b.len()`
    /// - Each `acc[i] < modulus`, `a[i] < modulus`, `b[i] < modulus`
    fn reduce_add_mul_slice_assign(self, acc: &mut [T], a: &[T], b: &[T]);

    /// Calculates `acc[i] = (acc[i] - a[i] * b[i]) (mod modulus)`
    /// element-wise, where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `acc.len() == a.len() == b.len()`
    /// - Each `acc[i] < modulus`, `a[i] < modulus`, `b[i] < modulus`
    fn reduce_sub_mul_slice_assign(self, acc: &mut [T], a: &[T], b: &[T]);

    /// Writes `a[i] * b[i] + c[i] (mod modulus)` into `output[i]`,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len() == c.len() == output.len()`
    /// - Each `a[i] < modulus`, `b[i] < modulus`, `c[i] < modulus`
    fn reduce_mul_add_slice_to(self, a: &[T], b: &[T], c: &[T], output: &mut [T]);

    /// Writes `scalar * b[i] + c[i] (mod modulus)` into `output[i]`,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `b.len() == c.len() == output.len()`
    /// - `scalar < modulus`, each `b[i] < modulus`, `c[i] < modulus`
    fn reduce_scalar_mul_add_slice_to(self, scalar: T, b: &[T], c: &[T], output: &mut [T]);
}
