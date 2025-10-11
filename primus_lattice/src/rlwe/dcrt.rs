use primus_integer::{UnsignedInteger, izip};
use primus_ntt::{Dcrt, DcrtTable, Ntt};
use primus_poly::{
    ArrayBase, Data, DataMut, DataOwned, NttPolynomial, PolyLength, RawData, dcrt::DcrtPolynomial,
};
use primus_reduce::FieldContext;

use crate::{CrtRlwe, CrtRlweInfo};

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone)]
pub struct DcrtRlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T> DcrtRlwe<S>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`DcrtRlwe<S>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T> DcrtRlwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a [`DcrtRlwe<S>`] with all entries equal to zero.
    #[inline]
    pub fn zero(info: CrtRlweInfo) -> Self {
        let len = info.moduli_count.0 * info.poly_length.0 * 2;
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }

    /// Perform element-wise modular addition of two [`DcrtRlwe<S>`].
    #[inline]
    pub fn add_element_wise<M, A>(
        mut self,
        rhs: &DcrtRlwe<A>,
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

    /// Perform element-wise modular subtraction of two [`DcrtRlwe<S>`].
    #[inline]
    pub fn sub_element_wise<M, A>(
        mut self,
        rhs: &DcrtRlwe<A>,
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

    /// inverse ntt transform
    #[inline]
    pub fn into_coeff_form<Table>(self, table: &Table) -> CrtRlwe<S>
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let poly_length = table.poly_length();

        let Self { mut data } = self;

        data.chunks_exact_mut(poly_length * 2)
            .zip(table.iter())
            .for_each(|(rlwe, ntt_table)| {
                rlwe.chunks_exact_mut(poly_length).for_each(|poly| {
                    ntt_table.inverse_transform_slice(poly);
                });
            });

        CrtRlwe::new(data)
    }
}

impl<S, T> DcrtRlwe<S>
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
    /// on the `self` [`DcrtRlwe<S>`] with another `rhs` [`DcrtRlwe<A>`].
    #[inline]
    pub fn add_assign_element_wise<M, A>(
        &mut self,
        rhs: &DcrtRlwe<A>,
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
    /// on the `self` [`DcrtRlwe<S>`] with another `rhs` [`DcrtRlwe<A>`].
    #[inline]
    pub fn sub_assign_element_wise<M, A>(
        &mut self,
        rhs: &DcrtRlwe<A>,
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

impl<S, T> DcrtRlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `result`.
    #[inline]
    pub fn add_inplace<M, A, B>(
        &self,
        rhs: &DcrtRlwe<A>,
        result: &mut DcrtRlwe<B>,
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
        rhs: &DcrtRlwe<A>,
        result: &mut DcrtRlwe<B>,
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

    /// Performs a multiplication on the `self` [`DcrtRlwe<S>`] with another `dcrt_polynomial` [`DcrtPolynomial<A>`],
    /// store the result into `result` [`DcrtRlwe<B>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, A, B>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<A>,
        result: &mut DcrtRlwe<B>,
        moduli: &[M],
        PolyLength(poly_length): PolyLength,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        izip!(
            self.data.chunks_exact(poly_length * 2),
            result.data.chunks_exact_mut(poly_length * 2),
            dcrt_polynomial.iter(poly_length),
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

    /// inverse ntt transform
    #[inline]
    pub fn to_coeff_form_inplace<Table, A>(&self, result: &mut CrtRlwe<A>, table: &Table)
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
        A: RawData<Elem = T> + DataMut,
    {
        result.data.copy_from_slice(self.data.as_ref());

        let poly_length = table.poly_length();

        result
            .data
            .chunks_exact_mut(poly_length * 2)
            .zip(table.iter())
            .for_each(|(rlwe, t)| {
                let (a, b) = unsafe { rlwe.split_at_mut_unchecked(poly_length) };
                t.inverse_transform_slice(a);
                t.inverse_transform_slice(b);
            });
    }
}
