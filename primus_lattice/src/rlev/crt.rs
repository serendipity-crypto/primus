use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, RawData};
use primus_reduce::FieldContext;

use crate::{CrtRlevInfo, DcrtRlev};

/// A representation of Ring Learning with Errors (RLWE) ciphertexts with respect to different base,
/// used to control noise growth in polynomial multiplications.
#[derive(Clone)]
pub struct CrtRlev<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T> CrtRlev<S, T>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`CrtRlev<S, T>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T> CrtRlev<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a [`CrtRlev<S>`] with all entries equal to zero.
    #[inline]
    pub fn zero(info: CrtRlevInfo) -> Self {
        let len = info.moduli_count * info.decompose_length * info.poly_length * 2;
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }

    /// Perform element-wise modular addition of two [`CrtRlev<S>`].
    #[inline]
    pub fn add_element_wise<M, A>(
        mut self,
        rhs: &CrtRlev<A>,
        moduli: &[M],
        info: CrtRlevInfo,
    ) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_assign_element_wise(rhs, moduli, info);
        self
    }

    /// Perform element-wise modular subtraction of two [`CrtRlev<S>`].
    #[inline]
    pub fn sub_element_wise<M, A>(
        mut self,
        rhs: &CrtRlev<A>,
        moduli: &[M],
        info: CrtRlevInfo,
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
    pub fn into_ntt_form<Table>(self, table: &Table, rlev_len: usize) -> DcrtRlev<S>
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let poly_length = table.poly_length();

        let Self { mut data } = self;

        data.chunks_exact_mut(rlev_len)
            .zip(table.iter())
            .for_each(|(rlev, ntt_table)| {
                rlev.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.transform_slice(poly);
                });
            });

        DcrtRlev::new(data)
    }
}

impl<S, T> CrtRlev<S, T>
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
    /// on the `self` [`CrtRlev<S>`] with another `rhs` [`CrtRlev<S>`].
    #[inline]
    pub fn add_assign_element_wise<M, A>(
        &mut self,
        rhs: &CrtRlev<A>,
        moduli: &[M],
        info: CrtRlevInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        let rlev_len = info.rlev_len;
        izip!(
            self.data.chunks_exact_mut(rlev_len),
            rhs.data.chunks_exact(rlev_len),
            moduli
        )
        .for_each(|(x, y, m)| {
            ArrayBase(x).add_assign(&ArrayBase(y), *m);
        });
    }

    /// Performs an in-place element-wise modular subtraction
    /// on the `self` [`CrtRlev<S>`] with another `rhs` [`CrtRlev<S>`].
    #[inline]
    pub fn sub_assign_element_wise<M, A>(
        &mut self,
        rhs: &CrtRlev<A>,
        moduli: &[M],
        info: CrtRlevInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        let rlev_len = info.rlev_len;
        izip!(
            self.data.chunks_exact_mut(rlev_len),
            rhs.data.chunks_exact(rlev_len),
            moduli
        )
        .for_each(|(x, y, m)| {
            ArrayBase(x).add_assign(&ArrayBase(y), *m);
        });
    }
}

impl<S, T> CrtRlev<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `result`.
    #[inline]
    pub fn add_inplace<M, A, B>(
        &self,
        rhs: &CrtRlev<A>,
        result: &mut CrtRlev<B>,
        moduli: &[M],
        info: CrtRlevInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let rlev_len = info.rlev_len;
        izip!(
            self.data.chunks_exact(rlev_len),
            rhs.data.chunks_exact(rlev_len),
            result.data.chunks_exact_mut(rlev_len),
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
        rhs: &CrtRlev<A>,
        result: &mut CrtRlev<B>,
        moduli: &[M],
        info: CrtRlevInfo,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let rlev_len = info.rlev_len;
        izip!(
            self.data.chunks_exact(rlev_len),
            rhs.data.chunks_exact(rlev_len),
            result.data.chunks_exact_mut(rlev_len),
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
        result: &mut DcrtRlev<A>,
        table: &Table,
        rlev_len: usize,
    ) where
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + DataMut,
    {
        result.data.copy_from_slice(self.data.as_ref());

        let poly_length = table.poly_length();

        result
            .data
            .chunks_exact_mut(rlev_len)
            .zip(table.iter())
            .for_each(|(rlev, ntt_table)| {
                rlev.chunks_exact_mut(poly_length).for_each(|a| {
                    ntt_table.transform_slice(a);
                });
            });
    }
}
