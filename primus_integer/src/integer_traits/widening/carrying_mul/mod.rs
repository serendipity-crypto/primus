mod primitive;
#[cfg(all(feature = "nightly", feature = "simd"))]
mod simd;

/// Carrying mul operation trait.
pub trait CarryingMul: Sized {
    /// Calculates the "full multiplication" `self` * `rhs` + `carry` without
    /// the possibility to overflow.
    ///
    /// This returns the low-order (wrapping) bits and the high-order (overflow) bits
    /// of the result as two separate values, in that order.
    ///
    /// Performs "long multiplication" which takes in an extra amount to add, and may return
    /// an additional amount of overflow. This allows for chaining together multiple multiplications
    /// to create "big integers" which represent larger values.
    fn carrying_mul(self, rhs: Self, carry: Self) -> (Self, Self);

    /// Calculates the "full multiplication" `self * rhs + carry1 + carry2`.
    ///
    /// This returns the low-order (wrapping) bits and the high-order (overflow) bits
    /// of the result as two separate values, in that order.
    ///
    /// This cannot overflow, as the double-width result has exactly enough
    /// space for the largest possible result. This is equivalent to how, in
    /// decimal, 9 × 9 + 9 + 9 = 81 + 18 = 99 = 9×10⁰ + 9×10¹ = 10² - 1.
    ///
    /// Performs "long multiplication" which takes in an extra amount to add, and may return an
    /// additional amount of overflow. This allows for chaining together multiple
    /// multiplications to create "big integers" which represent larger values.
    fn carrying_mul_add(self, rhs: Self, carry: Self, add: Self) -> (Self, Self);

    /// Calculates the "full multiplication" `self * rhs + carry` without
    /// the possibility to overflow.
    ///
    /// This returns only the high-order (overflow) bits of the result.
    ///
    /// Performs "long multiplication" which takes in an extra amount to add, and may return
    /// an additional amount of overflow. This allows for chaining together multiple multiplications
    /// to create "big integers" which represent larger values.
    fn carrying_mul_hw(self, rhs: Self, carry: Self) -> Self;

    /// Calculates the "full multiplication" `self * rhs + carry1 + carry2`.
    ///
    /// This returns only the high-order (overflow) bits of the result.
    ///
    /// This cannot overflow, as the double-width result has exactly enough
    /// space for the largest possible result. This is equivalent to how, in
    /// decimal, 9 × 9 + 9 + 9 = 81 + 18 = 99 = 9×10⁰ + 9×10¹ = 10² - 1.
    ///
    /// Performs "long multiplication" which takes in an extra amount to add, and may return an
    /// additional amount of overflow. This allows for chaining together multiple
    /// multiplications to create "big integers" which represent larger values.
    fn carrying_mul_add_hw(self, rhs: Self, carry: Self, add: Self) -> Self;
}
