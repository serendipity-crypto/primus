use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{
    ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;

use crate::CrtGlwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
#[derive(Clone)]
pub struct DcrtGlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(DcrtGlwe<S, T>);
impl_bytes_conversion!(DcrtGlwe<S, T>);
impl_zero!(DcrtGlwe<S, T>);
impl_basic_operation_multiple_modulus!(DcrtGlwe<S, T>);
impl_crt_intt!(DcrtGlwe<S, T>, CrtGlwe);

impl<S, T> DcrtGlwe<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
}

impl<S, T> DcrtGlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
}

impl<S, T> DcrtGlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs a multiplication on the `self` [`DcrtGlwe<S>`] with another `dcrt_polynomial` [`DcrtPolynomial<A>`],
    /// store the result into `result` [`DcrtGlwe<B>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, A, B>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<A>,
        result: &mut DcrtGlwe<B>,
        moduli: &[M],
        cipher_single_modulus_len: usize,
        poly_length: usize,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.data.chunks_exact(cipher_single_modulus_len),
            result.data.chunks_exact_mut(cipher_single_modulus_len),
            dcrt_polynomial.iter(poly_length),
            moduli
        )
        .for_each(|(glwe0, glwe1, poly, modulus)| {
            glwe0
                .chunks_exact(poly_length)
                .zip(glwe1.chunks_exact_mut(poly_length))
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
