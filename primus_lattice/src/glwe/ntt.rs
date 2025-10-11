use primus_integer::{UnsignedInteger, size::Size};
use primus_ntt::{Ntt, NttTable};
use primus_poly::{ArrayBase, Data, DataMut, DataOwned, NttPolynomial, RawData};
use primus_reduce::FieldContext;

use crate::{Glwe, GlweInfo};

/// A cryptographic structure for Module(General) Learning with Errors (MLWE, Glwe).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone)]
pub struct NttGlwe<S, T = <S as RawData>::Elem>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    pub data: ArrayBase<S>,
}

impl<S, T> NttGlwe<S>
where
    S: RawData<Elem = T>,
    T: UnsignedInteger,
{
    /// Creates a new [`NttGlwe<S>`].
    #[inline]
    pub fn new(data: ArrayBase<S>) -> Self {
        Self { data }
    }
}

impl<S, T> NttGlwe<S>
where
    S: RawData<Elem = T> + DataOwned,
    T: UnsignedInteger,
{
    /// Creates a new [`NttGlwe<S>`] that is initialized to zero,
    /// both `a` and `b` polynomials are initialized to zero.
    #[inline]
    pub fn zero(info: GlweInfo) -> Self {
        let len = (info.dimension + 1) * info.poly_length;
        Self {
            data: ArrayBase::from_vec(vec![T::ZERO; len]),
        }
    }

    /// Creates a new [`NttGlwe<S>`] from bytes `data`.
    #[inline]
    pub fn from_bytes(data: &[u8]) -> Self {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        Self {
            data: ArrayBase::from_slice(converted_data),
        }
    }

    /// Perform element-wise modular addition `self + rhs`.
    #[inline]
    pub fn add_element_wise<M, A>(mut self, rhs: &NttGlwe<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
        self
    }

    /// Perform element-wise modular subtraction `self - rhs`.
    #[inline]
    pub fn sub_element_wise<M, A>(mut self, rhs: &NttGlwe<A>, modulus: M) -> Self
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
        self
    }

    /// ntt inverse transform
    #[inline]
    pub fn into_coeff_form<Table>(mut self, ntt_table: &Table) -> Glwe<S>
    where
        Table: NttTable<ValueT = T> + Ntt,
    {
        let poly_length = ntt_table.poly_length();
        self.data.chunks_exact_mut(poly_length).for_each(|values| {
            ntt_table.inverse_transform_slice(values);
        });

        Glwe::new(self.data)
    }
}

impl<S, T> NttGlwe<S>
where
    S: RawData<Elem = T> + DataMut,
    T: UnsignedInteger,
{
    /// Creates a new [`NttGlwe<S>`] from bytes `data`.
    #[inline]
    pub fn from_bytes_assign(&mut self, data: &[u8]) {
        let converted_data: &[T] = bytemuck::cast_slice(data);

        self.data.copy_from_slice(converted_data);
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.data.set_zero();
    }

    /// Performs an in-place element-wise modular addition `self += rhs`.
    #[inline]
    pub fn add_assign_element_wise<M, A>(&mut self, rhs: &NttGlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.add_assign(&rhs.data, modulus);
    }

    /// Performs an in-place element-wise modular subtraction `self -= rhs`
    #[inline]
    pub fn sub_assign_element_wise<M, A>(&mut self, rhs: &NttGlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        self.data.sub_assign(&rhs.data, modulus);
    }

    /// Performs a modular multiplication on the `self` [`NttGlwe<S>`] with another `polynomial` [`NttPolynomial<A>`].
    #[inline]
    pub fn mul_ntt_polynomial_assign<M, A>(&mut self, polynomial: &NttPolynomial<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
    {
        let poly_length = polynomial.poly_length();

        self.data.chunks_exact_mut(poly_length).for_each(|p| {
            NttPolynomial(ArrayBase(p)).mul_assign(polynomial, modulus);
        });
    }
}

impl<S, T> NttGlwe<S>
where
    S: RawData<Elem = T> + Data,
    T: UnsignedInteger,
{
    /// Converts [`NttGlwe<S>`] into bytes.
    #[inline]
    pub fn to_bytes(&self) -> Vec<u8> {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        converted_data.to_vec()
    }

    /// Converts [`NttGlwe<S>`] into bytes, stored in `data`.
    #[inline]
    pub fn to_bytes_inplace(&self, data: &mut [u8]) {
        let converted_data: &[u8] = bytemuck::cast_slice(self.data.as_ref());

        data.copy_from_slice(converted_data);
    }

    /// Returns the bytes count of [`NttGlwe<S>`].
    #[inline]
    pub fn bytes_count(&self) -> usize {
        self.data.byte_count()
    }

    /// Performs element-wise modular addition:`result = self + rhs`,
    #[inline]
    pub fn add_inplace<M, A>(&self, rhs: &Self, result: &mut NttGlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.add_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// Performs element-wise modular addition:`result = self - rhs`,
    #[inline]
    pub fn sub_inplace<M, A>(&self, rhs: &Self, result: &mut NttGlwe<A>, modulus: M)
    where
        M: FieldContext<T>,
        A: RawData<Elem = T> + DataMut,
    {
        self.data.sub_inplace(&rhs.data, &mut result.data, modulus)
    }

    /// Performs a modular multiplication on the `self` [`NttGlwe<S>`] with another `polynomial` [`NttPolynomial`],
    /// stores the result into `result`.
    #[inline]
    pub fn mul_ntt_polynomial_inplace<M, A, B>(
        &self,
        polynomial: &NttPolynomial<A>,
        result: &mut NttGlwe<B>,
        modulus: M,
    ) where
        M: FieldContext<T>,
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        let poly_length = polynomial.poly_length();

        self.data
            .chunks_exact(poly_length)
            .zip(result.data.chunks_exact_mut(poly_length))
            .for_each(|(x, y)| {
                NttPolynomial(ArrayBase(x)).mul_inplace(
                    polynomial,
                    &mut NttPolynomial(ArrayBase(y)),
                    modulus,
                );
            });
    }

    /// ntt inverse transform
    #[inline]
    pub fn to_coeff_form_inplace<Table, A>(&self, result: &mut Glwe<A>, ntt_table: &Table)
    where
        A: RawData<Elem = T> + DataMut,
        Table: NttTable<ValueT = T> + Ntt,
    {
        result.data.copy_from_slice(self.data.as_ref());

        let poly_length = ntt_table.poly_length();
        result
            .data
            .chunks_exact_mut(poly_length)
            .for_each(|values| {
                ntt_table.inverse_transform_slice(values);
            });
    }
}
