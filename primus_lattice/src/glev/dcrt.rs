use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::{CrtGlev, CrtGlevInfo};

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct DcrtGlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T> DcrtGlev<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`DcrtGlev<S, T>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T> DcrtGlev<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a [`DcrtGlev<S>`] with all entries equal to zero.
    #[inline]
    pub fn zero(info: CrtGlevInfo) -> Self {
        let len =
            info.moduli_count * info.decompose_length * info.poly_length * (info.dimension + 1);
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }

    /// Perform element-wise modular addition of two [`DcrtGlev<S>`].
    #[inline]
    pub fn add_element_wise<M, A>(
        mut self,
        rhs: &DcrtGlev<A>,
        moduli: &[M],
        info: CrtGlevInfo,
    ) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_assign_element_wise(rhs, moduli, info);
        self
    }

    /// Perform element-wise modular subtraction of two [`DcrtGlev<S>`].
    #[inline]
    pub fn sub_element_wise<M, A>(
        mut self,
        rhs: &DcrtGlev<A>,
        moduli: &[M],
        info: CrtGlevInfo,
    ) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.sub_assign_element_wise(rhs, moduli, info);
        self
    }

    /// ntt inverse transform
    #[inline]
    pub fn into_coeff_form<Table>(self, table: &Table, glev_len: usize) -> CrtGlev<S>
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let poly_length = table.poly_length();

        let Self { mut data } = self;

        data.chunks_exact_mut(glev_len)
            .zip(table.iter())
            .for_each(|(glev, ntt_table)| {
                glev.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.inverse_transform_slice(poly);
                });
            });

        CrtGlev::new(data)
    }
}

impl<S, T> DcrtGlev<S, T>
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
    /// on the `self` [`DcrtGlev<S>`] with another `rhs` [`DcrtGlev<A>`].
    #[inline]
    pub fn add_assign_element_wise<M, A>(
        &mut self,
        rhs: &DcrtGlev<A>,
        moduli: &[M],
        info: CrtGlevInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        let glev_len = info.glev_len;
        izip!(
            self.data.chunks_exact_mut(glev_len),
            rhs.data.chunks_exact(glev_len),
            moduli
        )
        .for_each(|(x, y, m)| {
            ArrayBase(x).add_assign(&ArrayBase(y), *m);
        });
    }

    /// Performs an in-place element-wise modular subtraction
    /// on the `self` [`DcrtGlev<S>`] with another `rhs` [`DcrtGlev<A>`].
    #[inline]
    pub fn sub_assign_element_wise<M, A>(
        &mut self,
        rhs: &DcrtGlev<A>,
        moduli: &[M],
        info: CrtGlevInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        let glev_len = info.glev_len;
        izip!(
            self.data.chunks_exact_mut(glev_len),
            rhs.data.chunks_exact(glev_len),
            moduli
        )
        .for_each(|(x, y, m)| {
            ArrayBase(x).add_assign(&ArrayBase(y), *m);
        });
    }
}

impl<S, T> DcrtGlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `result`.
    #[inline]
    pub fn add_inplace<M, A, B>(
        &self,
        rhs: &DcrtGlev<A>,
        result: &mut DcrtGlev<B>,
        moduli: &[M],
        info: CrtGlevInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let glev_len = info.glev_len;
        izip!(
            self.data.chunks_exact(glev_len),
            rhs.data.chunks_exact(glev_len),
            result.data.chunks_exact_mut(glev_len),
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
        rhs: &DcrtGlev<A>,
        result: &mut DcrtGlev<B>,
        moduli: &[M],
        info: CrtGlevInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let glev_len = info.glev_len;
        izip!(
            self.data.chunks_exact(glev_len),
            rhs.data.chunks_exact(glev_len),
            result.data.chunks_exact_mut(glev_len),
            moduli
        )
        .for_each(|(x, y, z, m)| {
            ArrayBase(x).sub_inplace(&ArrayBase(y), &mut ArrayBase(z), *m);
        });
    }

    /// ntt inverse transform
    #[inline]
    pub fn to_coeff_form_inplace<Table, A>(
        &self,
        result: &mut CrtGlev<A>,
        table: &Table,
        glev_len: usize,
    ) where
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + DataMut,
    {
        result.data.copy_from_slice(self.data.as_ref());

        let poly_length = table.poly_length();

        result
            .data
            .chunks_exact_mut(glev_len)
            .zip(table.iter())
            .for_each(|(glev, ntt_table)| {
                glev.chunks_exact_mut(poly_length).for_each(|a| {
                    ntt_table.inverse_transform_slice(a);
                });
            });
    }
}
