use primus_integer::{Data, DataMut, DataOwned, RawData, Size, UnsignedInteger};
use primus_reduce::ops::*;
use serde::{Deserialize, Serialize};

use super::Lwe;

/// Represents a cryptographic structure based on the Learning with Errors (LWE) problem.
///
/// This structure encrypts several messages like a rlwe but truncated `b`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MultiMsgLwe<S>(pub S)
where
    S: RawData,
    <S as RawData>::Elem: UnsignedInteger;

impl<S, T> MultiMsgLwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`MultiMsgLwe<S, T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        Self(S::from_slice(converted_data))
    }

    /// Creates a new [`MultiMsgLwe<S, T>`].
    #[inline]
    pub fn new(data: S) -> Self {
        Self(data)
    }

    /// Generates a [`MultiMsgLwe<S, T>`] with all values are `0`.
    #[inline]
    pub fn zero(dimension: usize, msg_count: usize) -> Self {
        Self(S::from_vec(vec![T::ZERO; dimension + msg_count]))
    }
}

impl<S, T> MultiMsgLwe<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Creates a new [`MultiMsgLwe<S, T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        self.0.copy_from_slice(converted_data);
    }

    /// Returns mutable references to `a` and `b` of this [`MultiMsgLwe<S, T>`].
    #[inline]
    pub fn a_b_mut(&mut self, mid: usize) -> (&mut [T], &mut [T]) {
        self.0.split_at_mut(mid)
    }

    /// Sets all values to `0`.
    #[inline]
    pub fn set_zero(&mut self) {
        self.0.fill(T::ZERO);
    }

    /// Perform component-wise modular addition of two [`MultiMsgLwe<S, T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is not a reference.
    /// If your `self` is a reference, you can use function `add_component_wise_ref`.
    #[inline]
    pub fn add_component_wise<M, A>(mut self, rhs: &MultiMsgLwe<A>, modulus: M) -> Self
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.add_component_wise_assign(rhs, modulus);
        self
    }

    /// Performs an in-place component-wise modular addition
    /// on the `self` [`MultiMsgLwe<S, T>`] with another `rhs` [`MultiMsgLwe<S, T>`].
    #[inline]
    pub fn add_component_wise_assign<M, A>(&mut self, rhs: &MultiMsgLwe<A>, modulus: M)
    where
        M: Copy + ReduceAddAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        self.0
            .iter_mut()
            .zip(rhs.0.iter())
            .for_each(|(a, &b)| modulus.reduce_add_assign(a, b));
    }

    /// Perform component-wise modular subtraction of two [`MultiMsgLwe<S, T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is not a reference.
    /// If your `self` is a reference, you can use function `sub_component_wise_ref`.
    #[inline]
    pub fn sub_component_wise<M, A>(mut self, rhs: &MultiMsgLwe<A>, modulus: M) -> Self
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        self.sub_component_wise_assign(rhs, modulus);
        self
    }

    /// Performs an in-place component-wise modular subtraction
    /// on the `self` [`MultiMsgLwe<S, T>`] with another `rhs` [`MultiMsgLwe<S, T>`].
    #[inline]
    pub fn sub_component_wise_assign<M, A>(&mut self, rhs: &MultiMsgLwe<A>, modulus: M)
    where
        M: Copy + ReduceSubAssign<T>,
        A: RawData<Elem = T> + Data,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        self.0
            .iter_mut()
            .zip(rhs.0.iter())
            .for_each(|(a, &b)| modulus.reduce_sub_assign(a, b));
    }

    /// Performs an in-place modular scalar multiplication
    /// on the `self` [`MultiMsgLwe<S, T>`] with scalar `T`.
    #[inline]
    pub fn mul_scalar_assign<M>(&mut self, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAssign<T>,
    {
        self.0
            .iter_mut()
            .for_each(|v| modulus.reduce_mul_assign(v, scalar));
    }

    /// Performs an in-place modular scalar multiplication
    /// on the `rhs` [`MultiMsgLwe<S, T>`] with `scalar` `T`,
    /// then add to `self`.
    #[inline]
    pub fn add_rhs_mul_scalar_assign<M, A>(&mut self, rhs: &MultiMsgLwe<A>, scalar: T, modulus: M)
    where
        M: Copy + ReduceMulAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
    {
        self.0
            .iter_mut()
            .zip(rhs.0.iter())
            .for_each(|(v, &r)| *v = modulus.reduce_mul_add(r, scalar, *v));
    }
}

impl<S, T> MultiMsgLwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Converts [`MultiMsgLwe<S, T>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let data: &[u8] = bytemuck::cast_slice(self.0.as_slice());

        data.to_vec()
    }

    /// Converts [`MultiMsgLwe<S, T>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let src: &[u8] = bytemuck::cast_slice(self.0.as_slice());

        assert_eq!(data.len(), src.len());

        data.copy_from_slice(src);
    }

    /// Returns references to `a` and `b` of this [`MultiMsgLwe<S, T>`].
    #[inline]
    pub fn a_b(&self, mid: usize) -> (&[T], &[T]) {
        self.0.split_at(mid)
    }

    /// Perform component-wise modular addition of two [`MultiMsgLwe<S, T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is a reference.
    /// If your `self` is not a reference, you can use function `add_component_wise`.
    #[inline]
    pub fn add_component_wise_ref<M, A, B>(
        &self,
        rhs: &MultiMsgLwe<A>,
        modulus: M,
    ) -> MultiMsgLwe<B>
    where
        M: Copy + ReduceAdd<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataOwned,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        MultiMsgLwe::new(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(&a, &b)| modulus.reduce_add(a, b))
                .collect(),
        )
    }

    /// Perform component-wise modular subtraction of two [`MultiMsgLwe<S, T>`].
    ///
    /// # Attention
    ///
    /// In this function, `self` is a reference.
    /// If your `self` is not a reference, you can use function `sub_component_wise`.
    #[inline]
    pub fn sub_component_wise_ref<M, A, B>(
        &self,
        rhs: &MultiMsgLwe<A>,
        modulus: M,
    ) -> MultiMsgLwe<B>
    where
        M: Copy + ReduceSub<T, Output = T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataOwned,
    {
        debug_assert_eq!(self.0.len(), rhs.0.len());
        MultiMsgLwe::new(
            self.0
                .iter()
                .zip(rhs.0.iter())
                .map(|(&a, &b)| modulus.reduce_sub(a, b))
                .collect(),
        )
    }
}

impl<T: UnsignedInteger> MultiMsgLwe<Vec<T>> {
    /// Sample extract [`Lwe<Vec<T>>`].
    #[inline]
    pub fn extract_rlwe_mode<M>(&self, dimension: usize, index: usize, modulus: M) -> Lwe<Vec<T>>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let mut data = self.0[..dimension + 1].to_vec();
        if index == 0 {
            Lwe::new(data)
        } else {
            data[..dimension].rotate_right(index);
            data[..index]
                .iter_mut()
                .for_each(|v| modulus.reduce_neg_assign(v));
            data[dimension] = self.0[dimension + index];
            Lwe::new(data)
        }
    }

    /// Sample extract all [`Lwe<T>`].
    #[inline]
    pub fn extract_all<M>(&self, msg_count: usize, modulus: M) -> Vec<Lwe<Vec<T>>>
    where
        M: Copy + ReduceNegAssign<T>,
    {
        let dimension = self.0.len() - msg_count;
        let mut result = Vec::with_capacity(msg_count);

        let mut data = self.0[..dimension].to_vec();
        self.0[dimension + 1..].iter().for_each(|&b| {
            let lwe = Lwe::new(data.clone());
            result.push(lwe);

            data[..dimension].rotate_right(1);
            modulus.reduce_neg_assign(&mut data[0]);
            data[dimension] = b;
        });
        result.push(Lwe::new(data));

        result
    }
}

impl<S, T> Size for MultiMsgLwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    #[inline]
    fn byte_count(&self) -> usize {
        self.0.len() * T::BYTES
    }
}
