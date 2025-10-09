use std::iter::once;

use primus_integer::{ByteCount, UnsignedInteger};
use primus_poly::NttPolynomial;
use serde::{Deserialize, Serialize};

use crate::{NttGlev, NttGlwe};

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct NttGgsw<T: UnsignedInteger> {
    a: Vec<NttGlev<T>>,
    b: NttGlev<T>,
}

impl<T: UnsignedInteger> NttGgsw<T> {
    /// Creates a new [`NttGgsw<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(
        data: &[u8],
        decompose_length: usize,
        dimension: usize,
        poly_length: usize,
    ) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let glwe_a_len = dimension * poly_length;
        let glwe_len = glwe_a_len + poly_length;
        let glev_len = decompose_length * glwe_len;

        let (data_a, data_b) = converted_data.split_at(dimension * glev_len);

        let f = |chunk: &[T]| -> Vec<NttGlwe<T>> {
            chunk
                .chunks_exact(glwe_len)
                .map(|chunk| {
                    let (a, b) = unsafe { chunk.split_at_unchecked(glwe_a_len) };
                    NttGlwe {
                        a: a.chunks_exact(poly_length)
                            .map(NttPolynomial::from_slice)
                            .collect(),
                        b: NttPolynomial::from_slice(b),
                    }
                })
                .collect()
        };

        let a: Vec<NttGlev<T>> = data_a
            .chunks_exact(glev_len)
            .map(|chunk| NttGlev::new(f(chunk)))
            .collect();

        let b: Vec<NttGlwe<T>> = f(data_b);

        Self {
            a,
            b: NttGlev::new(b),
        }
    }

    /// Creates a new [`NttGgsw<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(
        &mut self,
        data: &[u8],
        decompose_length: usize,
        dimension: usize,
        poly_length: usize,
    ) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let glwe_a_len = dimension * poly_length;
        let glwe_len = glwe_a_len + poly_length;
        let glev_len = decompose_length * glwe_len;

        self.a
            .iter_mut()
            .chain(once(&mut self.b))
            .zip(converted_data.chunks_exact(glev_len))
            .for_each(|(glev, chunk): (&mut NttGlev<T>, &[T])| {
                glev.iter_mut()
                    .zip(chunk.chunks_exact(glwe_len))
                    .for_each(|(glwe, inchunk)| {
                        let (a, b) = unsafe { inchunk.split_at_unchecked(glwe_a_len) };
                        glwe.a
                            .iter_mut()
                            .zip(a.chunks_exact(poly_length))
                            .for_each(|(p, c)| p.copy_from(c));
                        glwe.b.copy_from(b);
                    })
            });
    }

    /// Converts [`Ggsw<T>`] into bytes.
    #[inline]
    pub fn into_bytes(
        &self,
        decompose_length: usize,
        dimension: usize,
        poly_length: usize,
    ) -> Vec<u8> {
        let size = (dimension + 1)
            * decompose_length
            * (dimension + 1)
            * poly_length
            * <T as ByteCount>::BYTES_COUNT;
        let mut result: Vec<u8> = Vec::with_capacity(size);

        self.a.iter().chain(once(&self.b)).for_each(|glev| {
            glev.iter().for_each(|glwe| {
                glwe.a.iter().for_each(|p| {
                    result.extend_from_slice(bytemuck::cast_slice(p.as_slice()));
                });
                result.extend_from_slice(bytemuck::cast_slice(glwe.b_slice()));
            })
        });

        result
    }

    /// Converts [`Ggsw<T>`] into bytes, stored in `data`.
    #[inline]
    pub fn into_bytes_inplace(
        &self,
        data: &mut [u8],
        decompose_length: usize,
        dimension: usize,
        poly_length: usize,
    ) {
        let poly_bytes_count = poly_length * <T as ByteCount>::BYTES_COUNT;
        let glwe_a_bytes_count = dimension * poly_bytes_count;
        let glwe_bytes_count = glwe_a_bytes_count + poly_bytes_count;
        let glev_bytes_count = decompose_length * glwe_bytes_count;

        self.a
            .iter()
            .chain(once(&self.b))
            .zip(data.chunks_exact_mut(glev_bytes_count))
            .for_each(|(glev, chunk): (&NttGlev<T>, &mut [u8])| {
                glev.iter()
                    .zip(chunk.chunks_exact_mut(glwe_bytes_count))
                    .for_each(|(glwe, glwe_chunk)| {
                        let (a, b) =
                            unsafe { glwe_chunk.split_at_mut_unchecked(glwe_a_bytes_count) };
                        glwe.a
                            .iter()
                            .zip(a.chunks_exact_mut(poly_bytes_count))
                            .for_each(|(ai, ai_chunk)| {
                                ai_chunk.copy_from_slice(bytemuck::cast_slice(ai.as_slice()));
                            });
                        b.copy_from_slice(bytemuck::cast_slice(glwe.b_slice()));
                    })
            });
    }
}
