use num_traits::Zero;
use primus_integer::UnsignedInteger;
use primus_reduce::ops::ReduceInvAssign;

use super::DcrtPolynomial;

impl<T: UnsignedInteger> DcrtPolynomial<T> {
    /// Try to calculate the inverse of the polynomial.
    #[inline]
    pub fn try_inv<M>(mut self, moduli: &[M]) -> Result<Self, Self>
    where
        M: Copy + ReduceInvAssign<T>,
    {
        if self.iter().any(|poly| poly.iter().any(Zero::is_zero)) {
            Err(self)
        } else {
            self.iter_mut().zip(moduli).for_each(|(poly, modulus)| {
                poly.iter_mut().for_each(|v| modulus.reduce_inv_assign(v))
            });
            Ok(self)
        }
    }
}
