use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::NttPolynomial;
use primus_utils::ByteCount;
use serde::{Deserialize, Serialize};

use crate::Glwe;

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, Glwe).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone, Default, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct NttGlwe<T: UnsignedInteger> {
    pub(crate) a: Vec<NttPolynomial<T>>,
    pub(crate) b: NttPolynomial<T>,
}

impl<T: UnsignedInteger> NttGlwe<T> {
    /// Creates a new [`NttGlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8], dimension: usize, poly_length: usize) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (a, b) = converted_data.split_at(dimension * poly_length);

        Self {
            a: a.chunks_exact(poly_length)
                .map(|s| NttPolynomial::from_slice(s))
                .collect(),
            b: NttPolynomial::from_slice(b),
        }
    }

    /// Creates a new [`NttGlwe<T>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8], dimension: usize, poly_length: usize) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        let (a, b) = converted_data.split_at(dimension * poly_length);

        self.a
            .iter_mut()
            .zip(a.chunks_exact(poly_length))
            .for_each(|(x, y)| x.copy_from(y));
        self.b.copy_from(b);
    }

    /// Converts [`NttGlwe<T>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let data_b: &[u8] = bytemuck::cast_slice(self.b.as_slice());

        let dimension = self.a.len();
        let poly_bytes_len = data_b.len();
        let mid = poly_bytes_len * dimension;

        let mut result: Vec<u8> = vec![0; mid + poly_bytes_len];

        let (a, b) = result.split_at_mut(mid);

        a.chunks_exact_mut(poly_bytes_len)
            .zip(self.a.iter())
            .for_each(|(x, y)| {
                x.copy_from_slice(bytemuck::cast_slice(y.as_slice()));
            });

        b.copy_from_slice(data_b);

        result
    }

    /// Converts [`NttGlwe<T>`] into bytes, stored in `data``.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let data_b: &[u8] = bytemuck::cast_slice(self.b.as_slice());

        let dimension = self.a.len();
        let poly_bytes_len = data_b.len();
        let mid = poly_bytes_len * dimension;

        assert_eq!(data.len(), mid + poly_bytes_len);

        let (a, b) = unsafe { data.split_at_mut_unchecked(mid) };

        a.chunks_exact_mut(poly_bytes_len)
            .zip(self.a.iter())
            .for_each(|(x, y)| {
                x.copy_from_slice(bytemuck::cast_slice(y.as_slice()));
            });
        b.copy_from_slice(data_b);
    }

    /// Returns the bytes count of [`NttGlwe<T>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        (self.a.len() + 1) * self.b.poly_length() * <T as ByteCount>::BYTES_COUNT
    }
}

impl<T: UnsignedInteger> NttGlwe<T> {
    /// Creates a new [`NttGlwe<T>`].
    #[inline]
    pub fn new(a: Vec<NttPolynomial<T>>, b: NttPolynomial<T>) -> Self {
        Self { a, b }
    }

    /// Creates a new [`NttGlwe<T>`] with reference of [`NttPolynomial<T>`].
    #[inline]
    pub fn from_ref(a: &[NttPolynomial<T>], b: &NttPolynomial<T>) -> Self {
        Self {
            a: a.to_vec(),
            b: b.clone(),
        }
    }

    /// Creates a new [`NttGlwe<T>`] that is initialized to zero,
    /// both `a` and `b` polynomials are initialized to zero.
    #[inline]
    pub fn zero(dimension: usize, poly_length: usize) -> Self {
        Self {
            a: (0..dimension)
                .map(|_| NttPolynomial::zero(poly_length))
                .collect(),
            b: NttPolynomial::zero(poly_length),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.iter_mut().for_each(|s| s.set_zero());
        self.b.set_zero();
    }

    /// Returns a reference to the a of this [`NttGlwe<T>`].
    pub fn a(&self) -> &[NttPolynomial<T>] {
        &self.a
    }

    /// Returns a mutable reference to the a of this [`NttGlwe<T>`].
    pub fn a_mut(&mut self) -> &mut [NttPolynomial<T>] {
        &mut self.a
    }

    /// Returns a reference to the b of this [`NttGlwe<T>`].
    pub fn b(&self) -> &NttPolynomial<T> {
        &self.b
    }

    /// Returns a mutable reference to the b of this [`NttGlwe<T>`].
    pub fn b_mut(&mut self) -> &mut NttPolynomial<T> {
        &mut self.b
    }

    /// Returns mutable references to the `a` and `b` of this [`NttGlwe<T>`].
    #[inline]
    pub fn a_b_mut(&mut self) -> (&mut [NttPolynomial<T>], &mut NttPolynomial<T>) {
        (&mut self.a, &mut self.b)
    }

    /// Extracts a slice of `b` of this [`NttGlwe<T>`].
    #[inline]
    pub fn b_slice(&self) -> &[T] {
        self.b.as_slice()
    }

    /// Extracts a mutable slice of `b` of this [`NttGlwe<T>`].
    #[inline]
    pub fn b_mut_slice(&mut self) -> &mut [T] {
        self.b.as_mut_slice()
    }
}

impl<T: UnsignedInteger> NttGlwe<T> {
    /// ntt inverse transform
    #[inline]
    pub fn into_coeff_form<Table>(self, ntt_table: &Table) -> Glwe<T>
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let Self { a, b } = self;

        let a = a
            .into_iter()
            .map(|p| ntt_table.inverse_transform_inplace(p))
            .collect();
        let b = ntt_table.inverse_transform_inplace(b);

        Glwe::new(a, b)
    }

    /// ntt inverse transform
    #[inline]
    pub fn inverse_transform_inplace<Table>(&self, ntt_table: &Table, result: &mut Glwe<T>)
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let (a, b) = result.a_b_mut();

        a.iter_mut().zip(self.a()).for_each(|(x, y)| x.copy_from(y));
        b.copy_from(self.b());

        a.iter_mut()
            .for_each(|p| ntt_table.inverse_transform_slice(p.as_mut_slice()));
        ntt_table.inverse_transform_slice(b.as_mut_slice());
    }
}
