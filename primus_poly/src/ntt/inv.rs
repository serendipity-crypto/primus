use num_traits::Zero;
use primus_integer::{DataMut, RawData, UnsignedInteger};
use primus_reduce::ops::ReduceInvAssign;

use super::NttPolynomial;

impl<S, T> NttPolynomial<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Try to calculate the inverse of the polynomial.
    #[inline]
    pub fn try_inv<M>(mut self, modulus: M) -> Result<Self, Self>
    where
        M: Copy + ReduceInvAssign<T>,
    {
        if self.iter().any(Zero::is_zero) {
            Err(self)
        } else {
            self.iter_mut().for_each(|v| modulus.reduce_inv_assign(v));
            Ok(self)
        }
    }
}
