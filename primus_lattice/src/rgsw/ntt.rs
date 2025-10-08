use primus_integer::UnsignedInteger;
use primus_poly::NttPolynomial;
use primus_utils::ByteCount;
use serde::{Deserialize, Serialize};

use crate::{NttRlev, NttRlwe};

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct NttRgsw<T: UnsignedInteger> {
    a: NttRlev<T>,
    b: NttRlev<T>,
}

impl<T: UnsignedInteger> NttRgsw<T> {
    /// Creates a new [`NttRgsw<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8], poly_length: usize) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (data_minus_s_m, data_m) = converted_data.split_at(converted_data.len() >> 1);

        let m: Vec<NttRlwe<T>> = data_m
            .chunks_exact(poly_length << 1)
            .map(|chunk| {
                let (a, b) = unsafe { chunk.split_at_unchecked(poly_length) };
                NttRlwe {
                    a: NttPolynomial::from_slice(a),
                    b: NttPolynomial::from_slice(b),
                }
            })
            .collect();

        let minus_s_m: Vec<NttRlwe<T>> = data_minus_s_m
            .chunks_exact(poly_length << 1)
            .map(|chunk| {
                let (a, b) = unsafe { chunk.split_at_unchecked(poly_length) };
                NttRlwe {
                    a: NttPolynomial::from_slice(a),
                    b: NttPolynomial::from_slice(b),
                }
            })
            .collect();

        Self {
            a: NttRlev::new(minus_s_m),
            b: NttRlev::new(m),
        }
    }

    /// Creates a new [`NttRgsw<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8], poly_length: usize) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        self.a
            .iter_mut()
            .chain(self.b.iter_mut())
            .zip(converted_data.chunks_exact(poly_length << 1))
            .for_each(|(rlwe, chunk): (&mut NttRlwe<T>, &[T])| {
                let (a, b) = unsafe { chunk.split_at_unchecked(poly_length) };
                rlwe.a.copy_from(a);
                rlwe.b.copy_from(b);
            });
    }

    /// Converts [`NttRgsw<T>`] into bytes.
    #[inline]
    pub fn into_bytes(&self, poly_length: usize) -> Vec<u8> {
        let size = (self.b.data().len() << 2) * poly_length * <T as ByteCount>::BYTES_COUNT;
        let mut result: Vec<u8> = Vec::with_capacity(size);

        self.a
            .iter()
            .chain(self.b.iter())
            .for_each(|rlwe: &NttRlwe<T>| {
                result.extend_from_slice(bytemuck::cast_slice(rlwe.a_slice()));
                result.extend_from_slice(bytemuck::cast_slice(rlwe.b_slice()));
            });

        result
    }

    /// Converts [`NttRgsw<T>`] into bytes, stored in `data``.
    #[inline]
    pub fn into_bytes_inplace(&self, data: &mut [u8], poly_length: usize) {
        let poly_bytes_count = poly_length * <T as ByteCount>::BYTES_COUNT;

        data.chunks_exact_mut(poly_bytes_count << 1)
            .zip(self.a.iter().chain(self.b.iter()))
            .for_each(|(chunk, rlwe): (&mut [u8], &NttRlwe<T>)| {
                let (a, b) = unsafe { chunk.split_at_mut_unchecked(poly_bytes_count) };
                a.copy_from_slice(bytemuck::cast_slice(rlwe.a_slice()));
                b.copy_from_slice(bytemuck::cast_slice(rlwe.b_slice()));
            });
    }

    // /// Returns the bytes count of [`NttRgsw<T>`].
    // #[inline]
    // pub fn bytes_count(&self) -> usize {
    //     self.m.bytes_count() << 1
    // }
}

impl<T: UnsignedInteger> NttRgsw<T> {
    /// Creates a new [`NttRgsw<T>`].
    #[inline]
    pub fn new(a: NttRlev<T>, b: NttRlev<T>) -> Self {
        Self { a, b }
    }

    /// Creates a new [`NttRgsw<T>`] with reference.
    #[inline]
    pub fn from_ref(a: &NttRlev<T>, b: &NttRlev<T>) -> Self {
        Self {
            a: a.clone(),
            b: b.clone(),
        }
    }

    /// Creates a [`NttRgsw<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(decompose_length: usize, poly_length: usize) -> Self {
        Self {
            a: NttRlev::zero(decompose_length, poly_length),
            b: NttRlev::zero(decompose_length, poly_length),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.set_zero();
        self.b.set_zero();
    }

    /// Returns a reference to the `-s*m` of this [`NttRgsw<T>`].
    #[inline]
    pub fn a(&self) -> &NttRlev<T> {
        &self.a
    }

    /// Returns a mutable reference to the `-s*m` of this [`NttRgsw<T>`].
    #[inline]
    pub fn a_mut(&mut self) -> &mut NttRlev<T> {
        &mut self.a
    }

    /// Returns a reference to the `m` of this [`NttRgsw<T>`].
    #[inline]
    pub fn b(&self) -> &NttRlev<T> {
        &self.b
    }

    /// Returns a mutable reference to the `m` of this [`NttRgsw<T>`].
    #[inline]
    pub fn b_mut(&mut self) -> &mut NttRlev<T> {
        &mut self.b
    }
}
