use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{
    ArrayBase, Data, DataMut, DataOwned, DcrtPolynomial, DcrtPolynomialIter, DcrtPolynomialIterMut,
    RawData,
};
use primus_reduce::FieldContext;

use super::CrtRlwe;

pub type DcrtRlweOwned<T> = DcrtRlwe<Vec<T>>;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
///
/// ## Structure of the `data`
///
/// |------a------|------b------|
///
/// where `a` and `b` are [`DcrtPolynomial`] with same poly length and moduli count.
#[derive(Clone)]
pub struct DcrtRlwe<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(DcrtRlwe<S, T>);
impl_bytes_conversion!(DcrtRlwe<S, T>);
impl_zero!(DcrtRlwe<S, T>);
impl_iters!(DcrtRlwe);
impl_iter_sub_structure!(DcrtRlwe<S, T>, DcrtPolynomial, dcrt_poly);
impl_basic_operation_multiple_modulus!(DcrtRlwe<S, T>);
impl_crt_intt!(DcrtRlwe<S, T>, CrtRlwe);

impl<S, T> DcrtRlwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs a multiplication on the `self` [`DcrtRlwe<S>`] with another `dcrt_poly` [`DcrtPolynomial<A>`],
    /// store the result into `result` [`DcrtRlwe<B>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, A, B>(
        &self,
        dcrt_poly: &DcrtPolynomial<A>,
        result: &mut DcrtRlwe<B>,
        poly_length: usize,
        moduli: &[M],
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let dcrt_poly_len = dcrt_poly.dcrt_poly_len();

        self.iter_dcrt_poly(dcrt_poly_len)
            .zip(result.iter_dcrt_poly_mut(dcrt_poly_len))
            .for_each(|(a, mut b)| {
                a.mul_inplace(dcrt_poly, &mut b, poly_length, moduli);
            });
    }
}
