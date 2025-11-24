use primus_distr::DiscreteGaussian;
use primus_integer::UnsignedInteger;
use primus_ntt::NttTable;
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, Polynomial, RawData};
use primus_reduce::{FieldContext, ops::ReduceNegAssign};
use rand::distr::Uniform;

use crate::lwe::{Lwe, MultiMsgLwe};

use super::Rlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
#[derive(Clone)]
pub struct TruncatedRlwe<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl_common!(TruncatedRlwe<S, T>);
impl_bytes_conversion!(TruncatedRlwe<S, T>);
impl_zero!(TruncatedRlwe<S, T>);
impl_basic_operation_single_modulus!(TruncatedRlwe<S, T>);

impl<T> TruncatedRlwe<Vec<T>, T>
where
    T: UnsignedInteger,
{
    /// Extract an LWE sample from RLWE.
    #[inline]
    pub fn extract_lwe_locally<M>(self, poly_length: usize, modulus: M) -> Lwe<T>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let Self(mut data) = self;
        let b = data[poly_length];
        data.truncate(poly_length);

        data[1..].reverse();
        data[1..]
            .iter_mut()
            .for_each(|v| modulus.reduce_neg_assign(v));

        Lwe::new(data, b)
    }

    /// Sample extract a [`MultiMsgLwe<T>`] with several encrypted messages.
    pub fn extract_first_few_lwe_locally<M>(
        self,
        count: usize,
        poly_length: usize,
        modulus: M,
    ) -> MultiMsgLwe<T>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let Self(mut data) = self;

        let b = data[poly_length..poly_length + count].to_vec();

        data.truncate(poly_length);

        data[1..].reverse();
        data[1..]
            .iter_mut()
            .for_each(|v| modulus.reduce_neg_assign(v));

        MultiMsgLwe::new(data, b)
    }

    /// Generate a [`TruncatedRlwe<Vec<T>>`] sample which encrypts `0`.
    pub fn generate_random_zero_sample<R, Table, M, A>(
        msg_count: usize,
        secret_key: &NttPolynomial<A>,
        uniform: Uniform<T>,
        gaussian: &DiscreteGaussian<T>,
        ntt_table: &Table,
        modulus: M,
        rng: &mut R,
    ) -> Self
    where
        R: rand::Rng + rand::CryptoRng,
        Table: NttTable<ValueT = T>,
        A: RawData<Elem = T> + Data,
        M: FieldContext<T>,
    {
        let poly_length = secret_key.poly_length();

        let mut cipher = Rlwe::zero(poly_length * 2);

        let (a, b) = cipher.a_b_mut_slices();

        Polynomial(&mut *a).random_with_distribution_assign(&uniform, rng);

        b.copy_from_slice(a);
        ntt_table.transform_slice(b);
        NttPolynomial(&mut *b).mul_assign(secret_key, modulus);
        ntt_table.inverse_transform_slice(b);

        Polynomial(b).add_random_gaussian_assign(gaussian, modulus, rng);

        let mut data: Vec<T> = cipher.0;

        data.truncate(poly_length + msg_count);

        TruncatedRlwe(data)
    }
}

impl<S, T> TruncatedRlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Extracts slice of `a` and `b` of this [`TruncatedRlwe<S>`].
    #[inline]
    pub fn a_b_slices(&self, poly_length: usize) -> (&[T], &[T]) {
        unsafe { self.0.split_at_unchecked(poly_length) }
    }
}
