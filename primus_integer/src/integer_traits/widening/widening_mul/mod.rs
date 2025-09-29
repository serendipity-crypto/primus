mod primitive;
#[cfg(all(feature = "nightly", feature = "simd"))]
mod simd;

/// Widening mul operation trait.
pub trait WideningMul: Sized {
    /// Calculates the complete product `self` * `rhs` without the possibility to overflow.
    ///
    /// This returns the low-order (wrapping) bits and the high-order (overflow) bits
    /// of the result as two separate values, in that order.
    fn widening_mul(self, rhs: Self) -> (Self, Self);

    /// Calculates the complete product `self` * `rhs` without the possibility to overflow.
    ///
    /// This returns only the high-order (overflow) bits of the result.
    fn widening_mul_hw(self, rhs: Self) -> Self;
}
