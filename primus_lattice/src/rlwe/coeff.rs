use primus_distr::DiscreteGaussian;
use primus_factor::ShoupFactor;
use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, Polynomial, RawData};
use primus_reduce::{
    FieldContext,
    ops::{ReduceNeg, ReduceNegAssign},
};

use crate::lwe::{Lwe, MultiMsgLwe};

use super::NttRlwe;

pub type RlweOwned<T> = Rlwe<Vec<T>>;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
///
/// ## Structure of the `data`
///
/// |------a------|------b------|
///
/// where `a` and `b` are [`Polynomial`] with same poly length.
#[derive(Clone)]
pub struct Rlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl_common!(Rlwe<S, T>);
impl_bytes_conversion!(Rlwe<S, T>);
impl_zero!(Rlwe<S, T>);
impl_iter_sub_structure!(Rlwe<S, T>, poly);
impl_basic_operation_single_modulus!(Rlwe<S, T>);
impl_ntt!(Rlwe<S, T>, NttRlwe);

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`Rlwe<S>`] with reference of [`Polynomial<A>`].
    #[inline]
    pub fn from_ref<A>(a: &Polynomial<A>, b: &Polynomial<A>) -> Self
    where
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(a.poly_length(), b.poly_length());
        Self {
            data: ArrayBase::from_vec([a.0.as_ref(), b.0.as_ref()].concat()),
        }
    }
}

impl<T> Rlwe<Vec<T>, T>
where
    T: UnsignedInteger,
{
    /// Extract an LWE sample from RLWE.
    #[inline]
    pub fn extract_lwe_locally<M>(self, modulus: M) -> Lwe<T>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let Self { data } = self;
        let ArrayBase(mut data) = data;
        let len = data.len() / 2;
        let b = data[len];
        data.truncate(len);

        data[1..].reverse();
        data[1..]
            .iter_mut()
            .for_each(|v| modulus.reduce_neg_assign(v));

        Lwe::new(data, b)
    }

    /// Sample extract a [`MultiMsgLwe<T>`] with several encrypted messages.
    pub fn extract_first_few_lwe_locally<M>(self, count: usize, modulus: M) -> MultiMsgLwe<T>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let Self { data } = self;
        let ArrayBase(mut data) = data;
        let len = data.len() / 2;

        let b = data[len..len + count].to_vec();

        data.truncate(len);

        data[1..].reverse();
        data[1..]
            .iter_mut()
            .for_each(|v| modulus.reduce_neg_assign(v));

        MultiMsgLwe::new(data, b)
    }

    /// Generate a [`Rlwe<Vec<T>>`] sample which encrypts `0`.
    pub fn generate_random_zero_sample<R, Table, M, A>(
        secret_key: &NttPolynomial<A>,
        gaussian: &DiscreteGaussian<T>,
        ntt_table: &Table,
        modulus: M,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        Table: NttTable<ValueT = T> + Ntt,
        A: RawData<Elem = T> + Data,
        M: FieldContext<T>,
    {
        let poly_length = secret_key.poly_length();

        let mut data = Rlwe::zero(poly_length * 2);

        let (a, b) = data.a_b_mut_slices();

        Polynomial(ArrayBase(&mut *a)).random_assign(modulus, rng);

        b.copy_from_slice(a);
        ntt_table.transform_slice(b);
        NttPolynomial(ArrayBase(&mut *b)).mul_assign(secret_key, modulus);
        ntt_table.inverse_transform_slice(b);

        Polynomial(ArrayBase(b)).add_random_gaussian_assign(gaussian, modulus, rng);

        data
    }
}

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Extracts mutable slice of `a` and `b` of this [`Rlwe<S>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self) -> (&mut [T], &mut [T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.0.split_at_mut_unchecked(mid) }
    }

    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, modulus: M)
    where
        M: FieldContext<T>,
    {
        ArrayBase(self.as_mut()).mul_scalar_assign(scalar, modulus);
    }

    #[inline]
    pub fn mul_factor_assign(&mut self, scalar: ShoupFactor<T>, modulus: T) {
        ArrayBase(self.as_mut()).mul_factor_assign(scalar, modulus);
    }
}

impl<S, T> Rlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts slice of `a` and `b` of this [`Rlwe<S>`].
    #[inline]
    pub fn a_b_slices(&self) -> (&[T], &[T]) {
        let mid = self.data.len() >> 1;
        unsafe { self.data.split_at_unchecked(mid) }
    }

    /// Performs a multiplication on the `self` [`Rlwe<S>`] with another `ntt_polynomial` [`NttPolynomial<A>`],
    /// store the result into `result` [`NttRlwe<B>`].
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, Table, A, B>(
        &self,
        ntt_polynomial: &NttPolynomial<A>,
        result: &mut NttRlwe<B>,
        modulus: M,
        ntt_table: &Table,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = ntt_table.poly_length();

        result.data.copy_from_slice(self.data.as_ref());

        result.data.chunks_exact_mut(poly_length).for_each(|poly| {
            ntt_table.transform_slice(poly);
            NttPolynomial(ArrayBase(poly)).mul_assign(ntt_polynomial, modulus);
        });
    }

    /// Extract an LWE sample from RLWE.
    #[inline]
    pub fn extract_lwe_with_index<M>(&self, index: usize, modulus: M) -> Lwe<T>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let split = index + 1;

        let (a, b) = self.a_b_slices();

        let mut a: Vec<T> = a.to_vec();

        a[..split].reverse();
        a[split..].reverse();
        a[split..]
            .iter_mut()
            .for_each(|x| modulus.reduce_neg_assign(x));

        Lwe::new(a, b[index])
    }

    /// Extract an LWE sample from RLWE.
    #[inline]
    pub fn extract_first_few_lwe<M>(&self, count: usize, modulus: M) -> MultiMsgLwe<T>
    where
        M: Copy + ReduceNeg<T, Output = T> + ReduceNegAssign<T>,
    {
        let (a, b) = self.a_b_slices();

        let mut a: Vec<_> = a.iter().map(|&x| modulus.reduce_neg(x)).collect();
        a[1..].reverse();
        modulus.reduce_neg_assign(&mut a[0]);

        MultiMsgLwe::new(a, b[..count].to_vec())
    }

    /// Extract an LWE sample from RLWE.
    #[inline]
    pub fn extract_lwe<M>(&self, modulus: M) -> Lwe<T>
    where
        M: Copy + ReduceNeg<T, Output = T> + ReduceNegAssign<T>,
    {
        let (a, b) = self.a_b_slices();

        let mut a: Vec<_> = a.iter().map(|&x| modulus.reduce_neg(x)).collect();
        a[1..].reverse();
        modulus.reduce_neg_assign(&mut a[0]);

        Lwe::new(a, b[0])
    }
}
