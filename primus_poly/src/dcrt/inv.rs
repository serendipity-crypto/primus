use num_traits::Zero;
use primus_integer::UnsignedInteger;
use primus_reduce::ops::ReduceInvAssign;

use crate::{DataOwned, PolyLength, RawData};

use super::DcrtPolynomial;

impl<S, T> DcrtPolynomial<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Try to calculate the inverse of the polynomial.
    #[inline]
    pub fn try_inv<M>(
        mut self,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) -> Result<Self, Self>
    where
        M: Copy + ReduceInvAssign<T>,
    {
        if self.0.iter().any(Zero::is_zero) {
            Err(self)
        } else {
            self.iter_mut(poly_length)
                .zip(moduli)
                .for_each(|(poly, modulus)| {
                    poly.iter_mut().for_each(|v| modulus.reduce_inv_assign(v))
                });
            Ok(self)
        }
    }
}
