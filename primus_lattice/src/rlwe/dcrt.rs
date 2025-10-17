use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData, dcrt::DcrtPolynomial};
use primus_reduce::FieldContext;

use super::CrtRlwe;

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
        poly_length: usize,
        moduli: &[M],
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let crt_poly_length = dcrt_polynomial.dcrt_poly_length();

        izip!(
            self.data.chunks_exact(crt_poly_length),
            result.data.chunks_exact_mut(crt_poly_length),
        )
        .for_each(|(a, b)| {
            DcrtPolynomial(ArrayBase(a)).mul_inplace(
                dcrt_polynomial,
                &mut DcrtPolynomial(ArrayBase(b)),
                poly_length,
                moduli,
            );
        });
    }
}
