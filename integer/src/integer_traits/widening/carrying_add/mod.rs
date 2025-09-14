mod primitive;
#[cfg(all(feature = "nightly", feature = "simd"))]
mod simd;

/// Carrying add operation trait
pub trait CarryingAdd: Sized {
    /// The type of `carry`.
    type CarryT;

    /// Calculates `self` + `rhs` + `carry` and checks for overflow.
    ///
    /// Performs “ternary addition” of two integer operands and a carry-in bit,
    /// and returns a tuple of the sum along with a boolean indicating
    /// whether an arithmetic overflow would occur. On overflow, the wrapped value is returned.
    ///
    /// This allows chaining together multiple additions to create a wider addition,
    /// and can be useful for bignum addition.
    /// This method should only be used for the most significant word.
    ///
    /// The output boolean returned by this method is not a carry flag,
    /// and should not be added to a more significant word.
    ///
    /// If the input carry is false, this method is equivalent to `overflowing_add`.
    fn carrying_add(self, rhs: Self, carry: Self::CarryT) -> (Self, Self::CarryT);
}
