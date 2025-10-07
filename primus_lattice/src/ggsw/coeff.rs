use std::iter::once;

use primus_integer::UnsignedInteger;
use primus_poly::Polynomial;
use serde::{Deserialize, Serialize};

use crate::{Glev, Glwe};

/// Represents a ciphertext in the Ring-GSW (Ring Learning With Errors) homomorphic encryption scheme.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound = "T: UnsignedInteger")]
pub struct Ggsw<T: UnsignedInteger> {
    a: Vec<Glev<T>>,
    b: Glev<T>,
}

impl<T: UnsignedInteger> Ggsw<T> {
    /// Creates a new [`Ggsw<T>`] from bytes `data`.
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

        let f = |chunk: &[T]| -> Vec<Glwe<T>> {
            chunk
                .chunks_exact(glwe_len)
                .map(|chunk| {
                    let (a, b) = unsafe { chunk.split_at_unchecked(glwe_a_len) };
                    Glwe {
                        a: a.chunks_exact(poly_length)
                            .map(Polynomial::from_slice)
                            .collect(),
                        b: Polynomial::from_slice(b),
                    }
                })
                .collect()
        };

        let a: Vec<Glev<T>> = data_a
            .chunks_exact(glev_len)
            .map(|chunk| Glev::new(f(chunk)))
            .collect();

        let b: Vec<Glwe<T>> = f(data_b);

        Self { a, b: Glev::new(b) }
    }

    /// Creates a new [`NttRgsw<F>`] from bytes `data`.
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
            .for_each(|(glev, chunk): (&mut Glev<T>, &[T])| {
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
}
