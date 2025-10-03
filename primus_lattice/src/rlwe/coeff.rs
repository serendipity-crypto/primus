use primus_distr::DiscreteGaussian;
use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{NttPolynomial, Polynomial};
use primus_reduce::{
    FieldContext,
    ops::{ReduceNeg, ReduceNegAssign},
};
use primus_utils::{ByteCount, Size};
use rand::{CryptoRng, Rng};
use serde::{Deserialize, Serialize};

use crate::{Lwe, MultiMsgLwe, NttRlwe};

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone, Default, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct Rlwe<T: UnsignedInteger> {
    /// Represents the first component in the RLWE structure.
    /// It is a polynomial where the coefficients are elements of the field `F`.
    pub(crate) a: Polynomial<T>,
    /// Represents the second component in the RLWE structure.
    /// It's also a polynomial with coefficients in the field `F`.
    pub(crate) b: Polynomial<T>,
}

impl<T: UnsignedInteger> Rlwe<T> {
    /// Creates a new [`Rlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (a, b) = converted_data.split_at(converted_data.len() >> 1);

        Self {
            a: Polynomial::from_slice(a),
            b: Polynomial::from_slice(b),
        }
    }

    /// Creates a new [`Rlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (a, b) = converted_data.split_at(converted_data.len() >> 1);

        self.a.copy_from(a);
        self.b.copy_from(b);
    }

    /// Converts [`Rlwe<T>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let data_a: &[u8] = bytemuck::cast_slice(self.a.as_slice());
        let data_b: &[u8] = bytemuck::cast_slice(self.b.as_slice());

        [data_a, data_b].concat()
    }

    /// Converts [`Rlwe<T>`] into bytes, stored in `data``.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let data_a: &[u8] = bytemuck::cast_slice(self.a.as_slice());
        let data_b: &[u8] = bytemuck::cast_slice(self.b.as_slice());

        assert_eq!(data.len(), data_a.len() + data_b.len());

        let (a, b) = unsafe { data.split_at_mut_unchecked(data_a.len()) };

        a.copy_from_slice(data_a);
        b.copy_from_slice(data_b);
    }

    /// Returns the bytes count of [`Rlwe<T>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        (self.a.poly_length() << 1) * <T as ByteCount>::BYTES_COUNT
    }
}

impl<T: UnsignedInteger> Rlwe<T> {
    /// Creates a new [`Rlwe<T>`].
    #[inline]
    pub fn new(a: Polynomial<T>, b: Polynomial<T>) -> Self {
        Self { a, b }
    }

    /// Creates a new [`Rlwe<T>`] with reference of [`Polynomial<T>`].
    #[inline]
    pub fn from_ref(a: &Polynomial<T>, b: &Polynomial<T>) -> Self {
        assert_eq!(a.poly_length(), b.poly_length());
        Self {
            a: a.clone(),
            b: b.clone(),
        }
    }

    /// Creates a new [`Rlwe<T>`] that is initialized to zero,
    /// both `a` and `b` polynomials are initialized to zero.
    #[inline]
    pub fn zero(poly_length: usize) -> Self {
        Self {
            a: Polynomial::zero(poly_length),
            b: Polynomial::zero(poly_length),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.set_zero();
        self.b.set_zero();
    }

    /// Returns a reference to the a of this [`Rlwe<T>`].
    #[inline]
    pub fn a(&self) -> &Polynomial<T> {
        &self.a
    }

    /// Returns a reference to the b of this [`Rlwe<T>`].
    #[inline]
    pub fn b(&self) -> &Polynomial<T> {
        &self.b
    }

    /// Returns a mutable reference to the a of this [`Rlwe<T>`].
    #[inline]
    pub fn a_mut(&mut self) -> &mut Polynomial<T> {
        &mut self.a
    }

    /// Returns a mutable reference to the b of this [`Rlwe<T>`].
    #[inline]
    pub fn b_mut(&mut self) -> &mut Polynomial<T> {
        &mut self.b
    }

    /// Returns mutable references to the `a` and `b` of this [`Rlwe<T>`].
    #[inline]
    pub fn a_b_mut(&mut self) -> (&mut Polynomial<T>, &mut Polynomial<T>) {
        (&mut self.a, &mut self.b)
    }

    /// Extracts a slice of `a` of this [`Rlwe<T>`].
    #[inline]
    pub fn a_slice(&self) -> &[T] {
        self.a.as_slice()
    }

    /// Extracts a slice of `b` of this [`Rlwe<T>`].
    #[inline]
    pub fn b_slice(&self) -> &[T] {
        self.b.as_slice()
    }

    /// Extracts a mutable slice of `a` of this [`Rlwe<T>`].
    #[inline]
    pub fn a_mut_slice(&mut self) -> &mut [T] {
        self.a.as_mut_slice()
    }

    /// Extracts a mutable slice of `b` of this [`Rlwe<T>`].
    #[inline]
    pub fn b_mut_slice(&mut self) -> &mut [T] {
        self.b.as_mut_slice()
    }

    /// Extracts mutable slice of `a` and `b` of this [`Rlwe<T>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self) -> (&mut [T], &mut [T]) {
        (self.a.as_mut_slice(), self.b.as_mut_slice())
    }
}

impl<T: UnsignedInteger> Rlwe<T> {
    /// ntt transform
    #[inline]
    pub fn into_ntt_form<Table>(self, ntt_table: &Table) -> NttRlwe<T>
    where
        Table: NttTable<ValueT = T> + Ntt<CoeffPoly = Polynomial<T>, NttPoly = NttPolynomial<T>>,
    {
        let Self { a, b } = self;

        let a = ntt_table.transform_inplace(a);
        let b = ntt_table.transform_inplace(b);

        NttRlwe::new(a, b)
    }

    /// ntt transform
    #[inline]
    pub fn transform_inplace<Table>(&self, ntt_table: &Table, result: &mut NttRlwe<T>)
    where
        Table: NttTable<ValueT = T> + Ntt<CoeffPoly = Polynomial<T>, NttPoly = NttPolynomial<T>>,
    {
        let (a, b) = result.a_b_mut_slices();

        a.copy_from_slice(self.a_slice());
        b.copy_from_slice(self.b_slice());

        ntt_table.transform_slice(a);
        ntt_table.transform_slice(b);
    }
}

impl<T: UnsignedInteger> Rlwe<T> {
    /// Perform element-wise modular addition of two [`Rlwe<T>`].
    #[inline]
    pub fn add_element_wise<M>(self, rhs: &Self, modulus: M) -> Self
    where
        M: FieldContext<T>,
    {
        Self {
            a: self.a.add(rhs.a(), modulus),
            b: self.b.add(rhs.b(), modulus),
        }
    }

    /// Perform element-wise modular subtraction of two [`Rlwe<T>`].
    #[inline]
    pub fn sub_element_wise<M>(self, rhs: &Self, modulus: M) -> Self
    where
        M: FieldContext<T>,
    {
        Self {
            a: self.a.sub(rhs.a(), modulus),
            b: self.b.sub(rhs.b(), modulus),
        }
    }

    /// Performs an in-place element-wise modular addition
    /// on the `self` [`Rlwe<T>`] with another `rhs` [`Rlwe<T>`].
    #[inline]
    pub fn add_assign_element_wise<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        self.a.add_assign(rhs.a(), modulus);
        self.b.add_assign(rhs.b(), modulus);
    }

    /// Performs an in-place element-wise modular subtraction
    /// on the `self` [`Rlwe<T>`] with another `rhs` [`Rlwe<T>`].
    #[inline]
    pub fn sub_assign_element_wise<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        self.a.sub_assign(rhs.a(), modulus);
        self.b.sub_assign(rhs.b(), modulus);
    }

    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `result`.
    #[inline]
    pub fn add_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        self.a.add_inplace(rhs.a(), result.a_mut(), modulus);
        self.b.add_inplace(rhs.b(), result.b_mut(), modulus);
    }

    /// Performs subtraction operation:`self - rhs`,
    /// and put the result to the `result`.
    #[inline]
    pub fn sub_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        self.a.sub_inplace(rhs.a(), result.a_mut(), modulus);
        self.b.sub_inplace(rhs.b(), result.b_mut(), modulus);
    }

    /// Performs a multiplication on the `self` [`Rlwe<T>`] with another `ntt_polynomial` [`NttPolynomial<T>`],
    /// store the result into `result` [`NttRlwe<T>`].
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, Table>(
        &self,
        ntt_polynomial: &NttPolynomial<T>,
        result: &mut NttRlwe<T>,
        modulus: M,
        ntt_table: &Table,
    ) where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt<CoeffPoly = Polynomial<T>, NttPoly = NttPolynomial<T>>,
    {
        let (a, b) = result.a_b_mut();

        a.copy_from(self.a());
        b.copy_from(self.b());

        ntt_table.transform_slice(a.as_mut_slice());
        ntt_table.transform_slice(b.as_mut_slice());

        a.mul_assign(ntt_polynomial, modulus);
        b.mul_assign(ntt_polynomial, modulus);
    }
}

impl<T: UnsignedInteger> Rlwe<T> {
    /// Extract an LWE sample from RLWE.
    #[inline]
    pub fn extract_lwe_with_index<M>(&self, index: usize, modulus: M) -> Lwe<T>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let split = index + 1;

        let mut a: Vec<_> = self.a_slice().to_vec();
        a[..split].reverse();
        a[split..].reverse();
        a[split..]
            .iter_mut()
            .for_each(|x| modulus.reduce_neg_assign(x));

        Lwe::new(a, self.b[index])
    }

    /// Extract an LWE sample from RLWE.
    #[inline]
    pub fn extract_first_few_lwe<M>(&self, count: usize, modulus: M) -> MultiMsgLwe<T>
    where
        M: Copy + ReduceNeg<T, Output = T> + ReduceNegAssign<T>,
    {
        let mut a: Vec<_> = self.a.iter().map(|&x| modulus.reduce_neg(x)).collect();
        a[1..].reverse();
        modulus.reduce_neg_assign(&mut a[0]);

        MultiMsgLwe::new(a, self.b[..count].to_vec())
    }

    /// Extract an LWE sample from RLWE.
    #[inline]
    pub fn extract_lwe<M>(&self, modulus: M) -> Lwe<T>
    where
        M: Copy + ReduceNeg<T, Output = T> + ReduceNegAssign<T>,
    {
        let mut a: Vec<_> = self.a.iter().map(|&x| modulus.reduce_neg(x)).collect();
        a[1..].reverse();
        modulus.reduce_neg_assign(&mut a[0]);

        Lwe::new(a, self.b[0])
    }

    /// Extract an LWE sample from RLWE.
    #[inline]
    pub fn extract_lwe_locally<M>(self, modulus: M) -> Lwe<T>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let Self { a, b } = self;
        let mut a = a.into_vec();

        a[1..].reverse();
        a[1..].iter_mut().for_each(|v| modulus.reduce_neg_assign(v));

        Lwe::new(a, b[0])
    }

    /// Sample extract a [`MultiMsgLwe<T>`] with several encrypted messages.
    pub fn extract_first_few_lwe_locally<M>(self, count: usize, modulus: M) -> MultiMsgLwe<T>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let Self { a, b } = self;
        let mut a = a.into_vec();
        let mut b = b.into_vec();

        b.truncate(count);

        a[1..].reverse();
        a[1..].iter_mut().for_each(|v| modulus.reduce_neg_assign(v));

        MultiMsgLwe::new(a, b)
    }
}

impl<T: UnsignedInteger> Rlwe<T> {
    /// Generate a [`Rlwe<T>`] sample which encrypts `0`.
    pub fn generate_random_zero_sample<M, Table, R>(
        secret_key: &NttPolynomial<T>,
        gaussian: &DiscreteGaussian<T>,
        ntt_table: &Table,
        modulus: M,
        rng: &mut R,
    ) -> Self
    where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt<CoeffPoly = Polynomial<T>, NttPoly = NttPolynomial<T>>,
        R: Rng + CryptoRng,
    {
        let rlwe_dimension = secret_key.poly_length();
        let a = <Polynomial<T>>::random(rlwe_dimension, modulus, rng);

        let mut a_ntt = ntt_table.transform(&a);
        a_ntt.mul_assign(secret_key, modulus);

        let mut e = <Polynomial<T>>::random_gaussian(rlwe_dimension, gaussian, rng);
        e.add_assign(&ntt_table.inverse_transform_inplace(a_ntt), modulus);

        Self { a, b: e }
    }
}

impl<T: UnsignedInteger> Size for Rlwe<T> {
    #[inline]
    fn size(&self) -> usize {
        self.a.size() * 2
    }
}
