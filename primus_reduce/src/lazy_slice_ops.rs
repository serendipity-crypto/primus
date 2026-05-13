//! Lazy slice-level modular operations.
//!
//! These traits mirror [`crate::lazy_ops`] but operate on whole slices.
//! Results are only guaranteed to be in `[0, 2 * modulus)`; callers must
//! perform a final reduction (e.g. via [`crate::ReduceOnceSlice`]) when a
//! canonical representative is required.
//!
//! See [`crate::slice_ops`] for the conventions on length checks
//! (`debug_assert*!`) and the lack of default impls.

/// Lazy slice form of [`crate::ReduceMul`] / [`crate::LazyReduceMul`].
pub trait LazyReduceMulSlice<T, B = T> {
    /// Calculates `a[i] = (a[i] * b[i]) (mod 2 * modulus)` element-wise,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len()`
    /// - Each `a[i] * b[i] < modulus²`
    /// - Each result is in `[0, 2 * modulus)`
    fn lazy_reduce_mul_slice_assign(self, a: &mut [T], b: &[B]);

    /// Writes `a[i] * b[i] (mod 2 * modulus)` into `output[i]`
    /// element-wise, where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len() == output.len()`
    /// - Each `a[i] * b[i] < modulus²`
    /// - Each result is in `[0, 2 * modulus)`
    fn lazy_reduce_mul_slice_to(self, a: &[T], b: &[B], output: &mut [T]);
}

/// Lazy slice form of [`crate::ReduceMulAdd`] / [`crate::LazyReduceMulAdd`].
///
/// Same four shapes as [`crate::ReduceMulAddSlice`]; results are in
/// `[0, 2 * modulus)`.
pub trait LazyReduceMulAddSlice<T> {
    /// Calculates `acc[i] = (acc[i] + a[i] * b[i]) (mod 2 * modulus)`
    /// element-wise, where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `acc.len() == a.len() == b.len()`
    /// - Each `acc[i] < modulus`, `a[i] < modulus`, `b[i] < modulus`
    /// - Each result is in `[0, 2 * modulus)`
    fn lazy_reduce_add_mul_slice_assign(self, acc: &mut [T], a: &[T], b: &[T]);

    /// Calculates `acc[i] = (acc[i] - a[i] * b[i]) (mod 2 * modulus)`
    /// element-wise, where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `acc.len() == a.len() == b.len()`
    /// - Each `acc[i] < modulus`, `a[i] < modulus`, `b[i] < modulus`
    /// - Each result is in `[0, 2 * modulus)`
    fn lazy_reduce_sub_mul_slice_assign(self, acc: &mut [T], a: &[T], b: &[T]);

    /// Writes `a[i] * b[i] + c[i] (mod 2 * modulus)` into `output[i]`,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `a.len() == b.len() == c.len() == output.len()`
    /// - Each `a[i] < modulus`, `b[i] < modulus`, `c[i] < modulus`
    /// - Each result is in `[0, 2 * modulus)`
    fn lazy_reduce_mul_add_slice_to(self, a: &[T], b: &[T], c: &[T], output: &mut [T]);

    /// Writes `scalar * b[i] + c[i] (mod 2 * modulus)` into `output[i]`,
    /// where `self` is the modulus.
    ///
    /// # Correctness
    ///
    /// - `b.len() == c.len() == output.len()`
    /// - `scalar < modulus`, each `b[i] < modulus`, `c[i] < modulus`
    /// - Each result is in `[0, 2 * modulus)`
    fn lazy_reduce_scalar_mul_add_slice_to(self, scalar: T, b: &[T], c: &[T], output: &mut [T]);
}
