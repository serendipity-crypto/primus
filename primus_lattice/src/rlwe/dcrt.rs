use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{
    ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;

use crate::CrtRlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
#[derive(Clone)]
pub struct DcrtRlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(DcrtRlwe<S, T>);
impl_bytes_conversion!(DcrtRlwe<S, T>);
impl_zero!(DcrtRlwe<S, T>);
impl_basic_operation_multiple_modulus!(DcrtRlwe<S, T>);
impl_crt_intt!(DcrtRlwe<S, T>, CrtRlwe);

impl<S, T> DcrtRlwe<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
}

impl<S, T> DcrtRlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
}

impl<S, T> DcrtRlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs a multiplication on the `self` [`DcrtRlwe<S>`] with another `dcrt_polynomial` [`DcrtPolynomial<A>`],
    /// store the result into `result` [`DcrtRlwe<B>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, A, B>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<A>,
        result: &mut DcrtRlwe<B>,
        moduli: &[M],
        poly_length: usize,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.data.chunks_exact(poly_length * 2),
            result.data.chunks_exact_mut(poly_length * 2),
            dcrt_polynomial.iter_each_modulus(poly_length),
            moduli
        )
        .for_each(|(rlwe0, rlwe1, poly, modulus)| {
            rlwe0
                .chunks_exact(poly_length)
                .zip(rlwe1.chunks_exact_mut(poly_length))
                .for_each(|(a0, a1)| {
                    NttPolynomial(ArrayBase(a0)).mul_inplace(
                        &NttPolynomial(ArrayBase(poly)),
                        &mut NttPolynomial(ArrayBase(a1)),
                        *modulus,
                    );
                });
        });
    }
}
