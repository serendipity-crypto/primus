use std::mem::MaybeUninit;

use primus_distr::DiscreteGaussian;
use primus_integer::{UnsignedInteger, size::Size};
use primus_modulo::ops::*;
use primus_poly::{Data, DataMut, DataOwned, RawData};
use primus_reduce::{Modulus, ops::*};
use rand::distr::{Distribution, Uniform};
use serde::{Deserialize, Serialize};

/// Represents a cryptographic structure based on the Learning with Errors (LWE) problem.
/// The LWE problem is a fundamental component in modern cryptography, often used to build
/// secure cryptographic systems that are considered hard to crack by quantum computers.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Lwe<S, T = <S as RawData>::Elem>(pub S)
where
    S: RawData<Elem = T>,
    T: UnsignedInteger;

impl<S, T> Lwe<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`Lwe<S, T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        Self(DataOwned::from_slice(converted_data))
    }
}

impl<S, T> Lwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Creates a new [`Lwe<S, T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        self.0.copy_from_slice(converted_data);
    }
}

impl<S, T> Lwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Converts [`Lwe<S, T>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let data: &[u8] = bytemuck::cast_slice(self.0.as_ref());
        data.to_vec()
    }

    /// Converts [`Lwe<S, T>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let src: &[u8] = bytemuck::cast_slice(self.0.as_ref());

        assert_eq!(data.len(), src.len());

        data.copy_from_slice(src);
    }
}

impl<S, T> Lwe<S, T>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`Lwe<S, T>`].
    #[inline]
    pub fn new(data: S) -> Self {
        Self(data)
    }

    // /// Creates a new [`Lwe<S, T>`] with reference.
    // #[inline]
    // pub fn from_ref(a: &[T], b: T) -> Self {
    //     Self { a: a.to_vec(), b }
    // }

    /// Generates a [`Lwe<S, T>`] with all values are `0`.
    #[inline]
    pub fn zero(dimension: usize) -> Self {
        Self(S::zero(dimension + 1))
    }
}

impl<T> Lwe<Vec<T>, T>
where
    T: UnsignedInteger,
{
    /// Generate a [`Lwe<S, T>`] sample which encrypts `0`.
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

        let mut data: Vec<MaybeUninit<T>> = Vec::with_capacity(len + 1);
        unsafe {
            data.set_len(len + 1);
        }
        data[0..len]
            .iter_mut()
            .zip(uniform.sample_iter(&mut *rng))
            .for_each(|(x, y)| {
                x.write(y);
            });
        data[len].write(gaussian.sample(rng));

        let mut data = unsafe { std::mem::transmute::<_, Vec<T>>(data) };

        let b = modulus.reduce_dot_product(&data[0..len], secret_key);
        data[len] = modulus.reduce_add(b, data[len]);

        Lwe(data)
    }
}

impl<S, T> Lwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Returns a mutable reference to the `a` of this [`Lwe<S, T>`].
    #[inline]
    pub fn a_mut(&mut self) -> &mut [T] {
        self.0.as_mut().split_last_mut().unwrap().1
    }

    /// Returns a mutable reference to the `b` of this [`Lwe<S, T>`].
    #[inline]
    pub fn b_mut(&mut self) -> &mut T {
        self.0.as_mut().last_mut().unwrap()
    }

    /// Returns mutable references to `a` and `b` of this [`Lwe<S, T>`].
    #[inline]
    pub fn a_b_mut(&mut self) -> (&mut [T], &mut T) {
        let (b, a) = self.0.as_mut().split_last_mut().unwrap();
        (a, b)
    }

    /// Sets all values to `0`.
    #[inline]
    pub fn set_zero(&mut self) {
        self.0.fill(T::ZERO);
    }
}

impl<S, T> Lwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Returns a reference to the `a` of this [`Lwe<S, T>`].
    #[inline]
    pub fn a(&self) -> &[T] {
        self.0.as_ref().split_last().unwrap().1
    }

    /// Returns the `b` of this [`Lwe<S, T>`].
    #[inline]
    pub fn b(&self) -> T {
        *self.0.as_ref().last().unwrap()
    }

    pub fn a_b(&self) -> (&[T], T) {
        let (b, a) = self.0.as_ref().split_last().unwrap();
        (a, *b)
    }

    /// Returns the dimension of this [`Lwe<S, T>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.0.len() - 1
    }
}

impl<S, T> Lwe<S, T>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Performs component-wise modular addition of two [`Lwe<S, T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is not a reference.
    /// If your `self` is a reference, you can use function `add_component_wise_ref`.
    #[inline]
    pub fn add_component_wise<M, A>(mut self, rhs: &Lwe<A>, modulus: M) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_component_wise_assign(rhs, modulus);
        self
    }

    /// Performs an in-place component-wise modular addition
    /// on the `self` [`Lwe<S, T>`] with another `rhs` [`Lwe<S, T>`].
    #[inline]
    pub fn add_component_wise_assign<M, A>(&mut self, rhs: &Lwe<A>, modulus: M)
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        self.0
            .iter_mut()
            .zip(rhs.0.iter())
            .for_each(|(x, &y)| x.add_modulo_assign(y, modulus));
    }

    /// Performs component-wise modular subtraction of two [`Lwe<S, T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is not a reference.
    /// If your `self` is a reference, you can use function `sub_component_wise_ref`.
    #[inline]
    pub fn sub_component_wise<M, A>(mut self, rhs: &Lwe<A>, modulus: M) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.sub_component_wise_assign(rhs, modulus);
        self
    }

    /// Performs an in-place component-wise modular subtraction
    /// on the `self` [`Lwe<S, T>`] with another `rhs` [`Lwe<S, T>`].
    #[inline]
    pub fn sub_component_wise_assign<M, A>(&mut self, rhs: &Lwe<A>, modulus: M)
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        self.0
            .iter_mut()
            .zip(rhs.0.iter())
            .for_each(|(x, &y)| x.sub_modulo_assign(y, modulus));
    }

    /// Performs an in-place modular scalar multiplication
    /// on the `self` [`Lwe<S, T>`] with scalar `T`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.0
            .iter_mut()
            .for_each(|v| v.mul_modulo_assign(scalar, modulus));
    }

    /// Performs an in-place modular scalar multiplication
    /// on the `rhs` [`Lwe<S, T>`] with `scalar` `T`,
    /// then add to `self`.
    #[inline]
    pub fn add_rhs_mul_scalar_assign<M, A>(&mut self, rhs: &Lwe<A>, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        self.0
            .iter_mut()
            .zip(rhs.0.iter())
            .for_each(|(v, &r)| *v = modulus.reduce_mul_add(r, scalar, *v));
    }

    /// Performs an negation on the `self` [`Lwe<S, T>`].
    #[inline]
    pub fn neg_assign<M>(&mut self, modulus: M)
    where
        M: Copy + ReduceNegAssign<T>,
    {
        self.0.iter_mut().for_each(|v| modulus.reduce_neg_assign(v));
    }
}

impl<S, T> Lwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Performs component-wise modular addition of two [`Lwe<S, T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is a reference.
    /// If your `self` is not a reference, you can use function `add_component_wise`.
    #[inline]
    pub fn add_component_wise_ref<M, A, B>(&self, rhs: &Lwe<A>, modulus: M) -> Lwe<B>
    where
        M: Copy + ReduceAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataOwned,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        Lwe::new(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(&a, &b)| a.add_modulo(b, modulus))
                .collect(),
        )
    }

    /// Performs component-wise modular subtraction of two [`Lwe<S, T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is a reference.
    /// If your `self` is not a reference, you can use function `sub_component_wise`.
    #[inline]
    pub fn sub_component_wise_ref<M, A, B>(&self, rhs: &Lwe<A>, modulus: M) -> Lwe<B>
    where
        M: Copy + ReduceSub<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataOwned,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        Lwe::new(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(&a, &b)| a.sub_modulo(b, modulus))
                .collect(),
        )
    }

    /// Performs a modular negation on the `self` [`Lwe<S, T>`].
    #[inline]
    pub fn neg<M, A>(&self, modulus: M) -> Lwe<A>
    where
        M: Copy + ReduceNeg<T, Output = T>,
        A: RawData<Elem = T> + DataOwned,
    {
        Lwe::new(self.0.iter().map(|&v| modulus.reduce_neg(v)).collect())
    }
}

impl<S, T> Size for Lwe<S, T>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn byte_count(&self) -> usize {
        self.0.len() * T::BYTES
    }
}
