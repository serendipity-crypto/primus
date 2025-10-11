use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::{CrtGlevInfo, DcrtGlev};

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct CrtGlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T> CrtGlev<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`CrtGlev<S, T>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T> CrtGlev<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a [`CrtGlev<S>`] with all entries equal to zero.
    #[inline]
    pub fn zero(info: CrtGlevInfo) -> Self {
        let len =
            info.moduli_count * info.decompose_length * info.poly_length * (info.dimension + 1);
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }

    /// Perform element-wise modular addition of two [`CrtGlev<S>`].
    #[inline]
    pub fn add_element_wise<M, A>(
        mut self,
        rhs: &CrtGlev<A>,
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

    /// Perform element-wise modular subtraction of two [`CrtGlev<S>`].
    #[inline]
    pub fn sub_element_wise<M, A>(
        mut self,
        rhs: &CrtGlev<A>,
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

    /// ntt transform
    #[inline]
    pub fn into_ntt_form<Table>(self, table: &Table, glev_len: usize) -> DcrtGlev<S>
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let poly_length = table.poly_length();

        let Self { mut data } = self;

        data.chunks_exact_mut(glev_len)
            .zip(table.iter())
            .for_each(|(rlev, ntt_table)| {
                rlev.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.transform_slice(poly);
                });
            });

        DcrtGlev::new(data)
    }
}

impl<S, T> CrtGlev<S, T>
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
    /// on the `self` [`CrtGlev<S>`] with another `rhs` [`CrtGlev<A>`].
    #[inline]
    pub fn add_assign_element_wise<M, A>(
        &mut self,
        rhs: &CrtGlev<A>,
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
    /// on the `self` [`CrtGlev<S>`] with another `rhs` [`CrtGlev<A>`].
    #[inline]
    pub fn sub_assign_element_wise<M, A>(
        &mut self,
        rhs: &CrtGlev<A>,
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

impl<S, T> CrtGlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `result`.
    #[inline]
    pub fn add_inplace<M, A, B>(
        &self,
        rhs: &CrtGlev<A>,
        result: &mut CrtGlev<B>,
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
        rhs: &CrtGlev<A>,
        result: &mut CrtGlev<B>,
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

    /// ntt transform
    #[inline]
    pub fn to_ntt_form_inplace<Table, A>(
        &self,
        result: &mut DcrtGlev<A>,
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
            .for_each(|(rlev, ntt_table)| {
                rlev.chunks_exact_mut(poly_length).for_each(|a| {
                    ntt_table.transform_slice(a);
                });
            });
    }
}
