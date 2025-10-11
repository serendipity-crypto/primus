use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{
    ArrayBase, Data, DataMut, DataOwned, NttPolynomial, PolyLength, RawData, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;

use crate::{CrtRlweInfo, DcrtRlwe};

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone)]
pub struct CrtRlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T: UnsignedInteger> CrtRlwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`CrtRlwe<T>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }

    /// Creates a [`CrtRlwe<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(info: CrtRlweInfo) -> Self {
        let len = info.moduli_count.0 * info.poly_length.0 * 2;
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }

    /// Perform element-wise modular addition of two [`CrtRlwe<S>`].
    #[inline]
    pub fn add_element_wise<M, A>(
        mut self,
        rhs: &CrtRlwe<A>,
        moduli: &[M],
        poly_length: PolyLength,
    ) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_assign_element_wise(rhs, moduli, poly_length);
        self
    }

    /// Perform element-wise modular subtraction of two [`CrtRlwe<S>`].
    #[inline]
    pub fn sub_element_wise<M, A>(
        mut self,
        rhs: &CrtRlwe<A>,
        moduli: &[M],
        poly_length: PolyLength,
    ) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.sub_assign_element_wise(rhs, moduli, poly_length);
        self
    }
}

impl<S, T: UnsignedInteger> CrtRlwe<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.set_zero();
    }

    /// Performs an in-place element-wise modular addition
    /// on the `self` [`CrtRlwe<T>`] with another `rhs` [`CrtRlwe<T>`].
    #[inline]
    pub fn add_assign_element_wise<M, A>(
        &mut self,
        rhs: &CrtRlwe<A>,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.data.chunks_exact_mut(poly_length * 2),
            rhs.data.chunks_exact(poly_length * 2),
            moduli
        )
        .for_each(|(x, y, m)| {
            ArrayBase(x).add_assign(&ArrayBase(y), *m);
        });
    }

    /// Performs an in-place element-wise modular subtraction
    /// on the `self` [`CrtRlwe<T>`] with another `rhs` [`CrtRlwe<T>`].
    #[inline]
    pub fn sub_assign_element_wise<M, A>(
        &mut self,
        rhs: &CrtRlwe<A>,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        izip!(
            self.data.chunks_exact_mut(poly_length * 2),
            rhs.data.chunks_exact(poly_length * 2),
            moduli
        )
        .for_each(|(x, y, m)| {
            ArrayBase(x).add_assign(&ArrayBase(y), *m);
        });
    }
}

impl<S, T: UnsignedInteger> CrtRlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `result`.
    #[inline]
    pub fn add_inplace<M, A, B>(
        &self,
        rhs: &CrtRlwe<A>,
        result: &mut CrtRlwe<B>,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.data.chunks_exact(poly_length * 2),
            rhs.data.chunks_exact(poly_length * 2),
            result.data.chunks_exact_mut(poly_length * 2),
            moduli
        )
        .for_each(|(x, y, z, m)| {
            ArrayBase(x).add_inplace(&ArrayBase(y), &mut ArrayBase(z), *m);
        });
    }

    /// Performs subtraction operation:`self - rhs`,
    /// and put the result to the `result`.
    #[inline]
    pub fn sub_inplace<M, A, B>(
        &self,
        rhs: &CrtRlwe<A>,
        result: &mut CrtRlwe<B>,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.data.chunks_exact(poly_length * 2),
            rhs.data.chunks_exact(poly_length * 2),
            result.data.chunks_exact_mut(poly_length * 2),
            moduli
        )
        .for_each(|(x, y, z, m)| {
            ArrayBase(x).sub_inplace(&ArrayBase(y), &mut ArrayBase(z), *m);
        });
    }

    /// Performs a multiplication on the `self` [`CrtRlwe<T>`] with another `dcrt_polynomial` [`DcrtPolynomial<T>`],
    /// store the result into `result` [`DcrtRlwe<T>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, Table, A, B>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<A>,
        result: &mut DcrtRlwe<B>,
        moduli: &[M],
        table: &Table,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        result.data.copy_from_slice(self.data.as_ref());

        let poly_length = table.poly_length();

        izip!(
            result.data.chunks_exact_mut(poly_length * 2),
            dcrt_polynomial.iter(poly_length),
            table.iter(),
            moduli
        )
        .for_each(|(x, p, t, m)| {
            let (a, b) = unsafe { x.split_at_mut_unchecked(poly_length) };
            t.transform_slice(a);
            t.transform_slice(b);
            NttPolynomial(ArrayBase(a)).mul_assign(&NttPolynomial(ArrayBase(p)), *m);
            NttPolynomial(ArrayBase(b)).mul_assign(&NttPolynomial(ArrayBase(p)), *m);
        });
    }
}
