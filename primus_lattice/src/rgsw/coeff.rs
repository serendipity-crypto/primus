use primus_integer::UnsignedInteger;
use primus_poly::Polynomial;
use primus_utils::ByteCount;
use serde::{Deserialize, Serialize};

use crate::{Rlev, Rlwe};

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
///
/// The [`Rgsw`] struct holds two components, `m` and `minus_s_m`, each of type [`Rlev`]. These components are
/// integral to the RGSW encryption scheme, which is a variant of GSW encryption that operates over polynomial
/// rings for efficiency. This scheme is often used in lattice-based cryptography for constructing fully
/// homomorphic encryption systems.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct Rgsw<T: UnsignedInteger> {
    /// The first part of the RGSW ciphertext, which is often used for homomorphic operations
    /// and can represent the encrypted data multiplied by some secret value.
    minus_s_m: Rlev<T>,
    /// The second part of the RGSW ciphertext, typically representing the encrypted data.
    m: Rlev<T>,
}

impl<T: UnsignedInteger> Rgsw<T> {
    /// Creates a new [`Rgsw<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8], poly_length: usize) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (data_minus_s_m, data_m) = converted_data.split_at(converted_data.len() >> 1);

        let m: Vec<Rlwe<T>> = data_m
            .chunks_exact(poly_length << 1)
            .map(|chunk| {
                let (a, b) = unsafe { chunk.split_at_unchecked(poly_length) };
                Rlwe {
                    a: Polynomial::from_slice(a),
                    b: Polynomial::from_slice(b),
                }
            })
            .collect();

        let minus_s_m: Vec<Rlwe<T>> = data_minus_s_m
            .chunks_exact(poly_length << 1)
            .map(|chunk| {
                let (a, b) = unsafe { chunk.split_at_unchecked(poly_length) };
                Rlwe {
                    a: Polynomial::from_slice(a),
                    b: Polynomial::from_slice(b),
                }
            })
            .collect();

        Self {
            minus_s_m: Rlev::new(minus_s_m),
            m: Rlev::new(m),
        }
    }

    /// Creates a new [`Rgsw<F>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8], poly_length: usize) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        self.minus_s_m
            .iter_mut()
            .chain(self.m.iter_mut())
            .zip(converted_data.chunks_exact(poly_length << 1))
            .for_each(|(rlwe, chunk): (&mut Rlwe<T>, &[T])| {
                let (a, b) = unsafe { chunk.split_at_unchecked(poly_length) };
                rlwe.a.copy_from(a);
                rlwe.b.copy_from(b);
            });
    }

    /// Converts [`Rgsw<F>`] into bytes.
    #[inline]
    pub fn into_bytes(&self, poly_length: usize) -> Vec<u8> {
        let size = (self.m.data().len() << 2) * poly_length * <T as ByteCount>::BYTES_COUNT;
        let mut result: Vec<u8> = Vec::with_capacity(size);

        self.minus_s_m
            .iter()
            .chain(self.m.iter())
            .for_each(|rlwe: &Rlwe<T>| {
                result.extend_from_slice(bytemuck::cast_slice(rlwe.a_slice()));
                result.extend_from_slice(bytemuck::cast_slice(rlwe.b_slice()));
            });

        result
    }

    /// Converts [`Rgsw<F>`] into bytes, stored in `data``.
    #[inline]
    pub fn into_bytes_inplace(&self, data: &mut [u8], poly_length: usize) {
        let poly_bytes_count = poly_length * <T as ByteCount>::BYTES_COUNT;

        data.chunks_exact_mut(poly_bytes_count << 1)
            .zip(self.minus_s_m.iter().chain(self.m.iter()))
            .for_each(|(chunk, rlwe): (&mut [u8], &Rlwe<T>)| {
                let (a, b) = unsafe { chunk.split_at_mut_unchecked(poly_bytes_count) };
                a.copy_from_slice(bytemuck::cast_slice(rlwe.a_slice()));
                b.copy_from_slice(bytemuck::cast_slice(rlwe.b_slice()));
            });
    }

    // /// Returns the bytes count of [`Rgsw<T>`].
    // #[inline]
    // pub fn bytes_count(&self) -> usize {
    //     self.m.bytes_count() << 1
    // }
}

impl<T: UnsignedInteger> Rgsw<T> {
    /// Creates a new [`Rgsw<T>`].
    #[inline]
    pub fn new(minus_s_m: Rlev<T>, m: Rlev<T>) -> Self {
        Self { minus_s_m, m }
    }

    /// Creates a new [`Rgsw<T>`] with reference.
    #[inline]
    pub fn from_ref(minus_s_m: &Rlev<T>, m: &Rlev<T>) -> Self {
        Self {
            minus_s_m: minus_s_m.clone(),
            m: m.clone(),
        }
    }

    /// Creates a [`Rgsw<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(decompose_length: usize, poly_length: usize) -> Self {
        Self {
            minus_s_m: Rlev::zero(decompose_length, poly_length),
            m: Rlev::zero(decompose_length, poly_length),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.minus_s_m.set_zero();
        self.m.set_zero();
    }

    /// Returns a reference to the `-s*m` of this [`Rgsw<F>`].
    #[inline]
    pub fn minus_s_m(&self) -> &Rlev<T> {
        &self.minus_s_m
    }

    /// Returns a mutable reference to the `-s*m` of this [`Rgsw<F>`].
    #[inline]
    pub fn minus_s_m_mut(&mut self) -> &mut Rlev<T> {
        &mut self.minus_s_m
    }

    /// Returns a reference to the `m` of this [`Rgsw<F>`].
    #[inline]
    pub fn m(&self) -> &Rlev<T> {
        &self.m
    }

    /// Returns a mutable reference to the `m` of this [`Rgsw<F>`].
    #[inline]
    pub fn m_mut(&mut self) -> &mut Rlev<T> {
        &mut self.m
    }
}
