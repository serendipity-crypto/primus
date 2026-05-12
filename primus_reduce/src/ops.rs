use primus_integer::UnsignedInteger;

use crate::ReduceError;

/// The reduce operation.
pub trait Reduce<T> {
    /// Output type.
    type Output;

    /// Calculates `value (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// The legal input range for `value` is **implementation-defined** and
    /// depends on the modulus form: Barrett admits up to `value < modulus²`,
    /// Montgomery requires `value < N·R`, native power-of-two admits any
    /// `T`. Consult the concrete modulus type's docs.
    ///
    /// When `T` is a slice (`&[U]`) or iterator-like input, implementations
    /// interpret the call as a sum-reduce (`Σ value_i mod modulus`); see
    /// [`FieldContext`](crate::FieldContext) for the bound shape.
    #[must_use]
    fn reduce(self, value: T) -> Self::Output;
}

/// The reduce assignment operation.
pub trait ReduceAssign<T> {
    /// Calculates `value (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// Input range is implementation-defined; see [`Reduce::reduce`].
    fn reduce_assign(self, value: &mut T);
}

/// At most one minus operation.
pub trait ReduceOnce<T> {
    /// Output type.
    type Output;

    /// Calculates `value - modulus` if `value >= modulus`.
    ///
    /// # Correctness
    ///
    /// - `value < 2 * modulus`
    /// - Result is `< modulus`
    #[must_use]
    fn reduce_once(self, value: T) -> Self::Output;
}

/// At most one minus operation assignment.
pub trait ReduceOnceAssign<T> {
    /// Calculates `value - modulus` if `value >= modulus`.
    ///
    /// # Correctness
    ///
    /// - `value < 2 * modulus`
    /// - Result is `< modulus`
    fn reduce_once_assign(self, value: &mut T);
}

/// The modular addition.
pub trait ReduceAdd<T, B = T> {
    /// Output type.
    type Output;

    /// Calculates `a + b (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a < modulus`
    /// - `b < modulus`
    #[must_use]
    fn reduce_add(self, a: T, b: B) -> Self::Output;
}

/// The modular addition assignment.
pub trait ReduceAddAssign<T, B = T> {
    /// Calculates `a += b (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a < modulus`
    /// - `b < modulus`
    fn reduce_add_assign(self, a: &mut T, b: B);
}

/// The modular double.
pub trait ReduceDouble<T> {
    /// Output type.
    type Output;

    /// Calculates `2*value (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `value < modulus`
    #[must_use]
    fn reduce_double(self, value: T) -> Self::Output;
}

/// The modular double assignment.
pub trait ReduceDoubleAssign<T> {
    /// Calculates `value = 2*value (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `value < modulus`
    fn reduce_double_assign(self, value: &mut T);
}

/// The modular subtraction.
pub trait ReduceSub<T, B = T> {
    /// Output type.
    type Output;

    /// Calculates `a - b (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a < modulus`
    /// - `b < modulus`
    #[must_use]
    fn reduce_sub(self, a: T, b: B) -> Self::Output;
}

/// The modular subtraction assignment.
pub trait ReduceSubAssign<T, B = T> {
    /// Calculates `a -= b (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a < modulus`
    /// - `b < modulus`
    fn reduce_sub_assign(self, a: &mut T, b: B);
}

/// The modular negation.
pub trait ReduceNeg<T> {
    /// Output type.
    type Output;

    /// Calculates `-value (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `value < modulus`
    #[must_use]
    fn reduce_neg(self, value: T) -> Self::Output;
}

/// The modular negation assignment.
pub trait ReduceNegAssign<T> {
    /// Calculates `-value (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `value < modulus`
    fn reduce_neg_assign(self, value: &mut T);
}

/// The modular multiplication.
pub trait ReduceMul<T, B = T> {
    /// Output type.
    type Output;

    /// Calculates `a * b (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a*b < modulus²`
    #[must_use]
    fn reduce_mul(self, a: T, b: B) -> Self::Output;
}

/// The modular multiplication assignment.
pub trait ReduceMulAssign<T, B = T> {
    /// Calculates `a *= b (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a*b < modulus²`
    fn reduce_mul_assign(self, a: &mut T, b: B);
}

/// The modular square.
pub trait ReduceSquare<T> {
    /// Output type.
    type Output;

    /// Calculates `value² (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `value < modulus`
    #[must_use]
    fn reduce_square(self, value: T) -> Self::Output;
}

/// The modular square assignment.
pub trait ReduceSquareAssign<T> {
    /// Calculates `value = value² (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `value < modulus`
    fn reduce_square_assign(self, value: &mut T);
}

/// The modular multiply-add.
pub trait ReduceMulAdd<T, B = T, C = T> {
    /// Output type.
    type Output;

    /// Calculates `(a * b) + c (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a < modulus`
    /// - `b < modulus`
    /// - `c < modulus`
    #[must_use]
    fn reduce_mul_add(self, a: T, b: B, c: C) -> Self::Output;
}

/// The modular multiply-add assignment.
pub trait ReduceMulAddAssign<T, B = T, C = T> {
    /// Calculates `(a * b) + c (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a < modulus`
    /// - `b < modulus`
    /// - `c < modulus`
    fn reduce_mul_add_assign(self, a: &mut T, b: B, c: C);
}

/// Calculate the inverse element for a field.
pub trait ReduceInv<T> {
    /// Output type.
    type Output;

    /// Calculate the multiplicative inverse of `value (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `value < modulus`
    /// - `value` and `modulus` must be coprime
    ///
    /// # Panics
    ///
    /// Panics if `value` has no inverse modulo `modulus`. Use
    /// [`TryReduceInv`] for a non-panicking variant.
    #[must_use]
    fn reduce_inv(self, value: T) -> Self::Output;
}

/// The modular inversion assignment for a field.
pub trait ReduceInvAssign<T> {
    /// Calculates `value^(-1) (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `value < modulus`
    /// - `value` and `modulus` must be coprime
    ///
    /// # Panics
    ///
    /// Panics if `value` has no inverse modulo `modulus`. Use
    /// [`TryReduceInv`] for a non-panicking variant.
    fn reduce_inv_assign(self, value: &mut T);
}

/// Try to calculate the inverse element when there may be not a field.
pub trait TryReduceInv<T> {
    /// Output type.
    type Output;

    /// Try to calculate the multiplicative inverse of `value modulo modulus` where `self` is modulus.
    ///
    /// # Errors
    ///
    /// If there dose not exist the such inverse, a [`ReduceError`] will be returned.
    fn try_reduce_inv(self, value: T) -> Result<Self::Output, ReduceError<T>>;
}

/// The modular division.
pub trait ReduceDiv<T, B = T> {
    /// Output type.
    type Output;

    /// Calculates `a / b (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a < modulus`
    /// - `b < modulus`
    /// - `b` and `modulus` must be coprime
    ///
    /// # Panics
    ///
    /// Panics if `b` has no inverse modulo `modulus`.
    #[must_use]
    fn reduce_div(self, a: T, b: B) -> Self::Output;
}

/// The modular division assignment.
pub trait ReduceDivAssign<T, B = T> {
    /// Calculates `a /= b (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `a < modulus`
    /// - `b < modulus`
    /// - `b` and `modulus` must be coprime
    ///
    /// # Panics
    ///
    /// Panics if `b` has no inverse modulo `modulus`.
    fn reduce_div_assign(self, a: &mut T, b: B);
}

/// The modular exponentiation.
pub trait ReduceExp<T> {
    /// Calculates `base^exp (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `base < modulus`
    #[must_use]
    fn reduce_exp<Exponent: UnsignedInteger>(self, base: T, exp: Exponent) -> T;
}

/// The modular exponentiation.
pub trait ReduceExpPowOf2<T> {
    /// Calculates `base^(2^exp_log) (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - `base < modulus`
    #[must_use]
    fn reduce_exp_power_of_2(self, base: T, exp_log: u32) -> T;
}

/// The modular dot product.
///
/// This is always used for slice. For example, `u64` slice `[u64]`.
///
/// For two same length slice `a = (a₀, a₁, ..., an)` and `b = (b₀, b₁, ..., bn)`.
///
/// This trait will calculate `a₀×b₀ + a₁×b₁ + ... + an×bn mod modulus`.
pub trait ReduceDotProduct<T> {
    /// Output type.
    type Output;

    /// Calculate `∑a_i×b_i (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - Each `a_i < modulus` and `b_i < modulus`
    ///
    /// # Panics
    ///
    /// Panics if `a.as_ref().len() != b.as_ref().len()`.
    #[must_use]
    fn reduce_dot_product(self, a: impl AsRef<[T]>, b: impl AsRef<[T]>) -> Self::Output;

    /// Calculate `∑a_i×b_i (mod modulus)` where `self` is modulus.
    ///
    /// # Correctness
    ///
    /// - Each `a_i < modulus` and `b_i < modulus`
    /// - If the iterators yield different numbers of elements, iteration
    ///   stops at the shorter (standard `zip` semantics); callers that
    ///   require equal length should use [`reduce_dot_product`](Self::reduce_dot_product).
    #[must_use]
    fn reduce_dot_product_iter(
        self,
        a: impl IntoIterator<Item = T>,
        b: impl IntoIterator<Item = T>,
    ) -> Self::Output;
}
