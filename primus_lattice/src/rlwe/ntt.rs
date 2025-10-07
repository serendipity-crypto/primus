use primus_distr::DiscreteGaussian;
use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{NttPolynomial, Polynomial};
use primus_reduce::FieldContext;
use primus_utils::{ByteCount, Size};
use rand::{CryptoRng, Rng};
use serde::{Deserialize, Serialize};

use crate::Rlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone, Default, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct NttRlwe<T: UnsignedInteger> {
    /// Represents the first component in the RLWE structure.
    pub(crate) a: NttPolynomial<T>,
    /// Represents the second component in the RLWE structure.
    pub(crate) b: NttPolynomial<T>,
}

impl<T: UnsignedInteger> NttRlwe<T> {
    /// Creates a new [`NttRlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (a, b) = converted_data.split_at(converted_data.len() >> 1);

        Self {
            a: NttPolynomial::from_slice(a),
            b: NttPolynomial::from_slice(b),
        }
    }

    /// Creates a new [`NttRlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (a, b) = converted_data.split_at(converted_data.len() >> 1);

        self.a.copy_from(a);
        self.b.copy_from(b);
    }

    /// Converts [`NttRlwe<T>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let data_a: &[u8] = bytemuck::cast_slice(self.a.as_slice());
        let data_b: &[u8] = bytemuck::cast_slice(self.b.as_slice());

        [data_a, data_b].concat()
    }

    /// Converts [`NttRlwe<T>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let data_a: &[u8] = bytemuck::cast_slice(self.a.as_slice());
        let data_b: &[u8] = bytemuck::cast_slice(self.b.as_slice());

        assert_eq!(data.len(), data_a.len() + data_b.len());

        let (a, b) = unsafe { data.split_at_mut_unchecked(data_a.len()) };

        a.copy_from_slice(data_a);
        b.copy_from_slice(data_b);
    }

    /// Returns the bytes count of [`NttRlwe<T>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        (self.a.poly_length() << 1) * <T as ByteCount>::BYTES_COUNT
    }
}

impl<T: UnsignedInteger> NttRlwe<T> {
    /// Creates a new [`NttRlwe<T>`].
    #[inline]
    pub fn new(a: NttPolynomial<T>, b: NttPolynomial<T>) -> Self {
        assert_eq!(a.poly_length(), b.poly_length());
        Self { a, b }
    }

    /// Creates a new [`NttRlwe<T>`] with reference of [`NttPolynomial<T>`].
    #[inline]
    pub fn from_ref(a: &NttPolynomial<T>, b: &NttPolynomial<T>) -> Self {
        assert_eq!(a.poly_length(), b.poly_length());
        Self {
            a: a.clone(),
            b: b.clone(),
        }
    }

    /// Creates a [`NttRlwe<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(coeff_count: usize) -> NttRlwe<T> {
        Self {
            a: <NttPolynomial<T>>::zero(coeff_count),
            b: <NttPolynomial<T>>::zero(coeff_count),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.set_zero();
        self.b.set_zero();
    }

    /// Returns a reference to the a of this [`NttRlwe<T>`].
    #[inline]
    pub fn a(&self) -> &NttPolynomial<T> {
        &self.a
    }

    /// Returns a mutable reference to the a of this [`NttRlwe<T>`].
    #[inline]
    pub fn a_mut(&mut self) -> &mut NttPolynomial<T> {
        &mut self.a
    }

    /// Returns a reference to the b of this [`NttRlwe<T>`].
    #[inline]
    pub fn b(&self) -> &NttPolynomial<T> {
        &self.b
    }

    /// Returns a mutable reference to the b of this [`NttRlwe<T>`].
    #[inline]
    pub fn b_mut(&mut self) -> &mut NttPolynomial<T> {
        &mut self.b
    }

    /// Returns a mutable reference to the `a` and `b` of this [`NttRlwe<T>`].
    #[inline]
    pub fn a_b_mut(&mut self) -> (&mut NttPolynomial<T>, &mut NttPolynomial<T>) {
        (&mut self.a, &mut self.b)
    }

    /// Extracts a slice of `a` of this [`NttRlwe<T>`].
    #[inline]
    pub fn a_slice(&self) -> &[T] {
        self.a.as_slice()
    }

    /// Extracts a mutable slice of `a` of this [`NttRlwe<T>`].
    #[inline]
    pub fn a_mut_slice(&mut self) -> &mut [T] {
        self.a.as_mut_slice()
    }

    /// Extracts a slice of `b` of this [`NttRlwe<T>`].
    #[inline]
    pub fn b_slice(&self) -> &[T] {
        self.b.as_slice()
    }

    /// Extracts a mutable slice of `b` of this [`NttRlwe<T>`].
    #[inline]
    pub fn b_mut_slice(&mut self) -> &mut [T] {
        self.b.as_mut_slice()
    }

    /// Extracts mutable slice of `a` and `b` of this [`NttRlwe<T>`].
    #[inline]
    pub fn a_b_mut_slices(&mut self) -> (&mut [T], &mut [T]) {
        (self.a.as_mut_slice(), self.b.as_mut_slice())
    }
}

impl<T: UnsignedInteger> NttRlwe<T> {
    /// ntt inverse transform
    #[inline]
    pub fn into_coeff_form<Table>(self, ntt_table: &Table) -> Rlwe<T>
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let Self { a, b } = self;

        let a = ntt_table.inverse_transform_inplace(a);
        let b = ntt_table.inverse_transform_inplace(b);

        Rlwe::new(a, b)
    }

    /// ntt inverse transform
    #[inline]
    pub fn to_coeff_form_inplace<Table>(&self, ntt_table: &Table, result: &mut Rlwe<T>)
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let (a, b) = result.a_b_mut_slices();

        a.copy_from_slice(self.a_slice());
        b.copy_from_slice(self.b_slice());

        ntt_table.inverse_transform_slice(a);
        ntt_table.inverse_transform_slice(b);
    }
}

impl<T: UnsignedInteger> NttRlwe<T> {
    /// Perform element-wise modular addition of two [`NttRlwe<T>`].
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

    /// Perform element-wise modular subtraction of two [`NttRlwe<T>`].
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
    /// on the `self` [`NttRlwe<T>`] with another `rhs` [`NttRlwe<T>`].
    #[inline]
    pub fn add_assign_element_wise<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        self.a.add_assign(rhs.a(), modulus);
        self.b.add_assign(rhs.b(), modulus);
    }

    /// Performs an in-place element-wise modular subtraction
    /// on the `self` [`NttRlwe<T>`] with another `rhs` [`NttRlwe<T>`].
    #[inline]
    pub fn sub_assign_element_wise<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        self.a.sub_assign(rhs.a(), modulus);
        self.b.sub_assign(rhs.b(), modulus);
    }

    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `destination`.
    #[inline]
    pub fn add_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        self.a.add_inplace(rhs.a(), result.a_mut(), modulus);
        self.b.add_inplace(rhs.b(), result.b_mut(), modulus);
    }

    /// Performs subtraction operation:`self - rhs`,
    /// and put the result to the `destination`.
    #[inline]
    pub fn sub_inplace<M>(&self, rhs: &Self, result: &mut Self, modulus: M)
    where
        M: FieldContext<T>,
    {
        self.a.sub_inplace(rhs.a(), result.a_mut(), modulus);
        self.b.sub_inplace(rhs.b(), result.b_mut(), modulus);
    }

    /// Performs a modular multiplication on the `self` [`NttRlwe<T>`] with another `polynomial` [`NttPolynomial<T>`].
    #[inline]
    pub fn mul_ntt_polynomial_assign<M>(&mut self, polynomial: &NttPolynomial<T>, modulus: M)
    where
        M: FieldContext<T>,
    {
        self.a.mul_assign(polynomial, modulus);
        self.b.mul_assign(polynomial, modulus);
    }

    /// Performs a modular multiplication on the `self` [`NttRlwe<T>`] with another `polynomial` [`NttPolynomial<T>`],
    /// stores the result into `result`.
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M>(
        &self,
        polynomial: &NttPolynomial<T>,
        result: &mut Self,
        modulus: M,
    ) where
        M: FieldContext<T>,
    {
        self.a.mul_inplace(polynomial, result.a_mut(), modulus);
        self.b.mul_inplace(polynomial, result.b_mut(), modulus);
    }

    /// Performs `self = self + ntt_rlwe * ntt_polynomial`.
    #[inline]
    pub fn add_ntt_rlwe_mul_ntt_polynomial_assign<M>(
        &mut self,
        ntt_rlwe: &Self,
        ntt_polynomial: &NttPolynomial<T>,
        modulus: M,
    ) where
        M: FieldContext<T>,
    {
        self.a_mut()
            .add_mul_assign(ntt_rlwe.a(), ntt_polynomial, modulus);
        self.b_mut()
            .add_mul_assign(ntt_rlwe.b(), ntt_polynomial, modulus);
    }

    /// Performs `self = self + ntt_rlwe * ntt_polynomial`.
    ///
    /// The result coefficients may be in [0, 2*modulus) for some case,
    /// and fall back to [0, modulus) for normal case.
    #[inline]
    pub fn add_ntt_rlwe_mul_ntt_polynomial_assign_fast<M>(
        &mut self,
        ntt_rlwe: &Self,
        ntt_polynomial: &NttPolynomial<T>,
        modulus: M,
    ) where
        M: FieldContext<T>,
    {
        self.a_mut()
            .add_mul_assign_fast(ntt_rlwe.a(), ntt_polynomial, modulus);
        self.b_mut()
            .add_mul_assign_fast(ntt_rlwe.b(), ntt_polynomial, modulus);
    }

    /// Performs `destination = self + ntt_rlwe * ntt_polynomial`.
    #[inline]
    pub fn add_ntt_rlwe_mul_ntt_polynomial_inplace<M>(
        &self,
        ntt_rlwe: &Self,
        ntt_polynomial: &NttPolynomial<T>,
        result: &mut Self,
        modulus: M,
    ) where
        M: FieldContext<T>,
    {
        ntt_rlwe
            .a()
            .mul_add_inplace(ntt_polynomial, self.a(), result.a_mut(), modulus);
        ntt_rlwe
            .b()
            .mul_add_inplace(ntt_polynomial, self.b(), result.b_mut(), modulus);
    }
}

impl<T: UnsignedInteger> NttRlwe<T> {
    /// Generate a [`NttRlwe<T>`] sample which encrypts `0`.
    pub fn generate_random_zero_sample<M, Table, R>(
        secret_key: &NttPolynomial<T>,
        gaussian: &DiscreteGaussian<T>,
        ntt_table: &Table,
        modulus: M,
        rng: &mut R,
    ) -> Self
    where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        R: Rng + CryptoRng,
    {
        let rlwe_dimension = secret_key.poly_length();
        let a = <NttPolynomial<T>>::random(rlwe_dimension, modulus, rng);

        let e = <Polynomial<T>>::random_gaussian(rlwe_dimension, gaussian, rng);
        let mut e = ntt_table.transform_inplace(e);
        e.add_mul_assign(&a, secret_key, modulus);

        Self { a, b: e }
    }

    /// Generate a [`NttRlwe<T>`] sample which encrypts `value`.
    pub fn generate_random_value_sample<M, Table, R>(
        secret_key: &NttPolynomial<T>,
        value: T,
        gaussian: &DiscreteGaussian<T>,
        ntt_table: &Table,
        modulus: M,
        rng: &mut R,
    ) -> Self
    where
        M: FieldContext<T>,
        Table: NttTable<ValueT = T> + Ntt,
        R: Rng + CryptoRng,
    {
        let rlwe_dimension = secret_key.poly_length();
        let a = <NttPolynomial<T>>::random(rlwe_dimension, modulus, rng);

        let mut e = <Polynomial<T>>::random_gaussian(rlwe_dimension, gaussian, rng);
        modulus.reduce_add_assign(&mut e[0], value);

        let mut b = ntt_table.transform_inplace(e);
        b.add_mul_assign(&a, secret_key, modulus);

        Self { a, b }
    }
}

impl<T: UnsignedInteger> Size for NttRlwe<T> {
    #[inline]
    fn size(&self) -> usize {
        self.a.size() * 2
    }
}
