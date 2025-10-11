use primus_integer::UnsignedInteger;
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};

use crate::{CrtRlwe, CrtRlweInfo};

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone)]
pub struct DcrtRlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T: UnsignedInteger> DcrtRlwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`DcrtRlwe<S>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }

    /// Creates a [`DcrtRlwe<S>`] with all entries equal to zero.
    #[inline]
    pub fn zero(info: CrtRlweInfo) -> Self {
        let len = info.moduli_count.0 * info.poly_length.0 * 2;
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }
}

impl<S, T: UnsignedInteger> DcrtRlwe<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.set_zero();
    }
}

impl<S, T: UnsignedInteger> DcrtRlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
}
