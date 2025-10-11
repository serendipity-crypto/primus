use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{
    ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;

use crate::{CrtGlweInfo, DcrtGlwe};

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone)]
pub struct CrtGlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T> CrtGlwe<S>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`CrtGlwe<S>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T> CrtGlwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a [`CrtGlwe<S>`] with all entries equal to zero.
    #[inline]
    pub fn zero(info: CrtGlweInfo) -> Self {
        let len = info.moduli_count * (info.dimension + 1) * info.poly_length;
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }

    /// Perform element-wise modular addition of two [`CrtGlwe<S>`].
    #[inline]
    pub fn add_element_wise<M, A>(
        mut self,
        rhs: &CrtGlwe<A>,
        moduli: &[M],
        info: CrtGlweInfo,
    ) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_assign_element_wise(rhs, moduli, info);
        self
    }

    /// Perform element-wise modular subtraction of two [`CrtGlwe<S>`].
    #[inline]
    pub fn sub_element_wise<M, A>(
        mut self,
        rhs: &CrtGlwe<A>,
        moduli: &[M],
        info: CrtGlweInfo,
    ) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.sub_assign_element_wise(rhs, moduli, info);
        self
    }

    /// ntt transform
    #[inline]
    pub fn into_ntt_form<Table>(self, table: &Table, dimension: usize) -> DcrtGlwe<S>
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let poly_length = table.poly_length();

        let Self { mut data } = self;

        data.chunks_exact_mut(poly_length * (dimension + 1))
            .zip(table.iter())
            .for_each(|(glwe, ntt_table)| {
                glwe.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.transform_slice(poly);
                });
            });

        DcrtGlwe::new(data)
    }
}

impl<S, T> CrtGlwe<S>
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
    /// on the `self` [`CrtGlwe<S>`] with another `rhs` [`CrtGlwe<A>`].
    #[inline]
    pub fn add_assign_element_wise<M, A>(
        &mut self,
        rhs: &CrtGlwe<A>,
        moduli: &[M],
        info: CrtGlweInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        let glwe_len = info.glwe_len;
        izip!(
            self.data.chunks_exact_mut(glwe_len),
            rhs.data.chunks_exact(glwe_len),
            moduli
        )
        .for_each(|(x, y, m)| {
            ArrayBase(x).add_assign(&ArrayBase(y), *m);
        });
    }

    /// Performs an in-place element-wise modular subtraction
    /// on the `self` [`CrtGlwe<S>`] with another `rhs` [`CrtGlwe<A>`].
    #[inline]
    pub fn sub_assign_element_wise<M, A>(
        &mut self,
        rhs: &CrtGlwe<A>,
        moduli: &[M],
        info: CrtGlweInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        let glwe_len = info.glwe_len;
        izip!(
            self.data.chunks_exact_mut(glwe_len),
            rhs.data.chunks_exact(glwe_len),
            moduli
        )
        .for_each(|(x, y, m)| {
            ArrayBase(x).add_assign(&ArrayBase(y), *m);
        });
    }
}

impl<S, T> CrtGlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `result`.
    #[inline]
    pub fn add_inplace<M, A, B>(
        &self,
        rhs: &CrtGlwe<A>,
        result: &mut CrtGlwe<B>,
        moduli: &[M],
        info: CrtGlweInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let glwe_len = info.glwe_len;
        izip!(
            self.data.chunks_exact(glwe_len),
            rhs.data.chunks_exact(glwe_len),
            result.data.chunks_exact_mut(glwe_len),
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
        rhs: &CrtGlwe<A>,
        result: &mut CrtGlwe<B>,
        moduli: &[M],
        info: CrtGlweInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let glwe_len = info.glwe_len;
        izip!(
            self.data.chunks_exact(glwe_len),
            rhs.data.chunks_exact(glwe_len),
            result.data.chunks_exact_mut(glwe_len),
            moduli
        )
        .for_each(|(x, y, z, m)| {
            ArrayBase(x).sub_inplace(&ArrayBase(y), &mut ArrayBase(z), *m);
        });
    }

    /// Performs a multiplication on the `self` [`CrtGlwe<S>`] with another `dcrt_polynomial` [`DcrtPolynomial<A>`],
    /// store the result into `result` [`DcrtGlwe<T>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, Table, A, B>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<A>,
        result: &mut DcrtGlwe<B>,
        moduli: &[M],
        table: &Table,
        info: CrtGlweInfo,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        result.data.copy_from_slice(self.data.as_ref());

        let poly_length = table.poly_length();
        let glwe_len = info.glwe_len;

        izip!(
            result.data.chunks_exact_mut(glwe_len),
            dcrt_polynomial.iter(poly_length),
            table.iter(),
            moduli
        )
        .for_each(|(glwe, poly, ntt_table, modulus)| {
            glwe.chunks_exact_mut(poly_length).for_each(|a| {
                ntt_table.transform_slice(a);
                NttPolynomial(ArrayBase(a)).mul_assign(&NttPolynomial(ArrayBase(poly)), *modulus);
            });
        });
    }

    /// ntt transform
    #[inline]
    pub fn to_ntt_form_inplace<Table, A>(
        &self,
        result: &mut DcrtGlwe<A>,
        table: &Table,
        info: CrtGlweInfo,
    ) where
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + DataMut,
    {
        result.data.copy_from_slice(self.data.as_ref());

        let poly_length = table.poly_length();
        let glwe_len = info.glwe_len;

        result
            .data
            .chunks_exact_mut(glwe_len)
            .zip(table.iter())
            .for_each(|(glwe, ntt_table)| {
                glwe.chunks_exact_mut(poly_length).for_each(|a| {
                    ntt_table.transform_slice(a);
                });
            });
    }
}
