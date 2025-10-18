use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData, dcrt::DcrtPolynomial};
use primus_reduce::FieldContext;

use super::CrtGlwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, GLWE).
///
/// ## Structure of the `data`
///
/// |--a1--|....|--ak--|--b--|
///
/// where `a1`...`ak` and `b` are [`DcrtPolynomial`] with same poly length and moduli count, `k` is the dimension.
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
impl_iter_sub_structure!(DcrtGlwe<S, T>, dcrt_poly);
impl_basic_operation_multiple_modulus!(DcrtGlwe<S, T>);
impl_crt_intt!(DcrtGlwe<S, T>, CrtGlwe);

impl<S, T> DcrtGlwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    pub fn add_dcrt_glwe_mul_dcrt_polynomial_assign<M, A, B>(
        &mut self,
        dcrt_glwe: &DcrtGlwe<A>,
        dcrt_polynomial: &DcrtPolynomial<B>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + Data,
    {
        let dcrt_poly_length = dcrt_polynomial.dcrt_poly_length();
        izip!(
            self.data.chunks_exact_mut(dcrt_poly_length),
            dcrt_glwe.data.chunks_exact(dcrt_poly_length)
        )
        .for_each(|(x, y)| {
            DcrtPolynomial(ArrayBase(x)).add_mul_assign(
                &DcrtPolynomial(ArrayBase(y)),
                dcrt_polynomial,
                poly_length,
                moduli,
            );
        });
    }
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
