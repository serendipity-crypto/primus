use bytemuck::Pod;
use primus_distr::DiscreteGaussian;
use primus_integer::{ByteCount, UnsignedInteger, size::Size};
use primus_modulo::ops::*;
use primus_reduce::{Modulus, ops::*};
use rand::distr::{Distribution, Uniform};
use serde::{Deserialize, Serialize};

/// Represents a cryptographic structure based on the Learning with Errors (LWE) problem.
/// The LWE problem is a fundamental component in modern cryptography, often used to build
/// secure cryptographic systems that are considered hard to crack by quantum computers.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Lwe<T: Copy> {
    /// A vector of elements of `T`, representing the public vector part of the LWE instance.
    a: Vec<T>,
    /// An element of `T`, representing the value which is computed as
    /// the dot product of `a` with a secret vector, plus message and some noise.
    b: T,
}

impl<T: Copy + Pod + ByteCount> Lwe<T> {
    /// Creates a new [`Lwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (&b, a) = converted_data.split_last().unwrap();

        Self { a: a.to_vec(), b }
    }

    /// Creates a new [`Lwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (&b, a) = converted_data.split_last().unwrap();

        self.a.copy_from_slice(a);
        self.b = b;
    }

    /// Converts [`Lwe<T>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let b = &[self.b];

        let data_a: &[u8] = bytemuck::cast_slice(&self.a);
        let data_b: &[u8] = bytemuck::cast_slice(b);

        [data_a, data_b].concat()
    }

    /// Converts [`Lwe<T>`] into bytes, stored in `data``.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let b = &[self.b];

        let data_a: &[u8] = bytemuck::cast_slice(&self.a);
        let data_b: &[u8] = bytemuck::cast_slice(b);

        assert_eq!(data.len(), data_a.len() + data_b.len());

        let (a, b) = unsafe { data.split_at_mut_unchecked(data_a.len()) };

        a.copy_from_slice(data_a);
        b.copy_from_slice(data_b);
    }

    /// Returns the bytes count of [`Lwe<T>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        (self.a.len() + 1) * T::BYTES_COUNT
    }
}

impl<T: Copy> Lwe<T> {
    /// Creates a new [`Lwe<T>`].
    #[inline]
    pub fn new(a: Vec<T>, b: T) -> Self {
        Self { a, b }
    }

    /// Creates a new [`Lwe<T>`] with reference.
    #[inline]
    pub fn from_ref(a: &[T], b: T) -> Self {
        Self { a: a.to_vec(), b }
    }

    /// Returns a reference to the `a` of this [`Lwe<T>`].
    #[inline]
    pub fn a(&self) -> &[T] {
        self.a.as_ref()
    }

    /// Returns a mutable reference to the `a` of this [`Lwe<T>`].
    #[inline]
    pub fn a_mut(&mut self) -> &mut Vec<T> {
        &mut self.a
    }

    /// Returns a slice reference to the a of this [`Lwe<T>`].
    #[inline]
    pub fn a_slice(&self) -> &[T] {
        self.a.as_slice()
    }

    /// Returns a mutable slice reference to the a of this [`Lwe<T>`].
    #[inline]
    pub fn a_mut_slice(&mut self) -> &mut [T] {
        self.a.as_mut_slice()
    }

    /// Returns the `b` of this [`Lwe<T>`].
    #[inline]
    pub fn b(&self) -> T {
        self.b
    }

    /// Returns a mutable reference to the `b` of this [`Lwe<T>`].
    #[inline]
    pub fn b_mut(&mut self) -> &mut T {
        &mut self.b
    }

    /// Returns the dimension of this [`Lwe<T>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.a.len()
    }
}

impl<T: UnsignedInteger> Lwe<T> {
    /// Generates a [`Lwe<T>`] with all values are `0`.
    #[inline]
    pub fn zero(dimension: usize) -> Self {
        Self {
            a: vec![T::ZERO; dimension],
            b: T::ZERO,
        }
    }

    /// Sets all values to `0`.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.fill(T::ZERO);
        self.b = T::ZERO;
    }

    /// Performs component-wise modular addition of two [`Lwe<T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is a reference.
    /// If your `self` is not a reference, you can use function `add_component_wise`.
    #[inline]
    pub fn add_component_wise_ref<M>(&self, rhs: &Self, modulus: M) -> Self
    where
        M: Copy + ReduceAdd<T, Output = T>,
    {
        debug_assert_eq!(self.a.len(), rhs.a.len());
        Self::new(
            self.a
                .iter()
                .zip(rhs.a())
                .map(|(&a, &b)| a.add_modulo(b, modulus))
                .collect(),
            self.b.add_modulo(rhs.b, modulus),
        )
    }

    /// Performs component-wise modular addition of two [`Lwe<T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is not a reference.
    /// If your `self` is a reference, you can use function `add_component_wise_ref`.
    #[inline]
    pub fn add_component_wise<M>(mut self, rhs: &Self, modulus: M) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
    {
        self.add_component_wise_assign(rhs, modulus);
        self
    }

    /// Performs an in-place component-wise modular addition
    /// on the `self` [`Lwe<T>`] with another `rhs` [`Lwe<T>`].
    #[inline]
    pub fn add_component_wise_assign<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: Copy + ReduceAddAssign<T>,
    {
        debug_assert_eq!(self.a.len(), rhs.a.len());
        self.a
            .iter_mut()
            .zip(rhs.a())
            .for_each(|(a, &b)| a.add_modulo_assign(b, modulus));
        self.b.add_modulo_assign(rhs.b, modulus);
    }

    /// Performs component-wise modular subtraction of two [`Lwe<T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is a reference.
    /// If your `self` is not a reference, you can use function `sub_component_wise`.
    #[inline]
    pub fn sub_component_wise_ref<M>(&self, rhs: &Self, modulus: M) -> Self
    where
        M: Copy + ReduceSub<T, Output = T>,
    {
        debug_assert_eq!(self.a.len(), rhs.a.len());
        Self::new(
            self.a
                .iter()
                .zip(rhs.a())
                .map(|(&a, &b)| a.sub_modulo(b, modulus))
                .collect(),
            self.b.sub_modulo(rhs.b, modulus),
        )
    }

    /// Performs component-wise modular subtraction of two [`Lwe<T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is not a reference.
    /// If your `self` is a reference, you can use function `sub_component_wise_ref`.
    #[inline]
    pub fn sub_component_wise<M>(mut self, rhs: &Self, modulus: M) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
    {
        self.sub_component_wise_assign(rhs, modulus);
        self
    }

    /// Performs an in-place component-wise modular subtraction
    /// on the `self` [`Lwe<T>`] with another `rhs` [`Lwe<T>`].
    #[inline]
    pub fn sub_component_wise_assign<M>(&mut self, rhs: &Self, modulus: M)
    where
        M: Copy + ReduceSubAssign<T>,
    {
        debug_assert_eq!(self.a.len(), rhs.a.len());
        self.a
            .iter_mut()
            .zip(rhs.a())
            .for_each(|(a, &b)| a.sub_modulo_assign(b, modulus));
        self.b.sub_modulo_assign(rhs.b, modulus)
    }

    /// Performs an in-place modular scalar multiplication
    /// on the `self` [`Lwe<T>`] with scalar `T`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.a
            .iter_mut()
            .for_each(|v| v.mul_modulo_assign(scalar, modulus));
        self.b.mul_modulo_assign(scalar, modulus);
    }

    /// Performs an in-place modular scalar multiplication
    /// on the `rhs` [`Lwe<T>`] with `scalar` `T`,
    /// then add to `self`.
    #[inline]
    pub fn add_rhs_mul_scalar_assign<M>(&mut self, rhs: &Self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
    {
        self.a
            .iter_mut()
            .zip(rhs.a())
            .for_each(|(v, &r)| *v = modulus.reduce_mul_add(r, scalar, *v));
        self.b = modulus.reduce_mul_add(rhs.b, scalar, self.b);
    }

    /// Performs a modular negation on the `self` [`Lwe<T>`].
    #[inline]
    pub fn neg<M>(&self, modulus: M) -> Self
    where
        M: Copy + ReduceNeg<T, Output = T>,
    {
        let a = self.a.iter().map(|&v| modulus.reduce_neg(v)).collect();
        Self::new(a, modulus.reduce_neg(self.b))
    }

    /// Performs an negation on the `self` [`Lwe<T>`].
    #[inline]
    pub fn neg_assign<M>(&mut self, modulus: M)
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.a.iter_mut().for_each(|v| modulus.reduce_neg_assign(v));
        modulus.reduce_neg_assign(&mut self.b)
    }

    /// Generate a [`Lwe<T>`] sample which encrypts `0`.
    #[inline]
    pub fn generate_random_zero_sample<M, R>(
        secret_key: &[T],
        modulus: M,
        uniform: Uniform<T>,
        gaussian: &DiscreteGaussian<T>,
        rng: &mut R,
    ) -> Self
    where
        M: Copy + Modulus<ValueT = T> + ReduceDotProduct<T> + ReduceAdd<T, Output = T>,
        R: rand::Rng + rand::CryptoRng,
    {
        let len = secret_key.len();

        let a: Vec<T> = uniform.sample_iter(&mut *rng).take(len).collect();
        let e = gaussian.sample(rng);

        let b = modulus.reduce_dot_product(a.as_slice(), secret_key);
        let b = modulus.reduce_add(b, e);

        Lwe { a, b }
    }
}

impl<T: Copy + ByteCount> Size for Lwe<T> {
    #[inline]
    fn byte_count(&self) -> usize {
        (self.a.len() + 1) * T::BYTES_COUNT
    }
}
