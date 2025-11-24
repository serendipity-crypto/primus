use primus_poly::{CrtPolynomial, DataMut, DcrtPolynomial, RawData};
use primus_reduce::FieldContext;

use crate::{NttError, NttTable, PrimitiveRoot};

#[cfg(feature = "concrete-ntt")]
mod concrete;
mod primitive;

#[cfg(feature = "concrete-ntt")]
pub use concrete::prime32::CrtConcrete32Table;
#[cfg(feature = "concrete-ntt")]
pub use concrete::prime64::CrtConcrete64Table;
pub use primitive::UintCrtNttTable;

pub trait DcrtTable: Sized {
    /// The value type.
    type ValueT: PrimitiveRoot;

    type NttTables: NttTable<ValueT = Self::ValueT>;

    /// Creates a new [`DcrtTable`].
    fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, NttError<Self::ValueT>>
    where
        M: FieldContext<Self::ValueT>;

    /// Returns a reference to the ntt tables.
    fn ntt_tables(&self) -> &[Self::NttTables];

    /// Returns an iterator over the ntt tables.
    fn iter(&self) -> std::slice::Iter<'_, Self::NttTables>;

    fn poly_length(&self) -> usize;

    fn moduli_count(&self) -> usize;

    fn crt_poly_length(&self) -> usize;

    /// Perform a fast number theory transform in place.
    ///
    /// This function transforms a crt polynomial to a dcrt polynomial.
    ///
    /// # Arguments
    ///
    /// * `crt_poly` - inputs in normal order, outputs in bit-reversed order
    fn transform_inplace<S: RawData<Elem = Self::ValueT> + DataMut>(
        &self,
        crt_poly: CrtPolynomial<S, Self::ValueT>,
    ) -> DcrtPolynomial<S, Self::ValueT>;

    /// Perform a fast inverse number theory transform in place.
    ///
    /// This function transforms a dcrt polynomial to a crt polynomial.
    ///
    /// # Arguments
    ///
    /// * `dcrt_poly` - inputs in bit-reversed order, outputs in normal order
    fn inverse_transform_inplace<S: RawData<Elem = Self::ValueT> + DataMut>(
        &self,
        dcrt_poly: DcrtPolynomial<S, Self::ValueT>,
    ) -> CrtPolynomial<S, Self::ValueT>;

    /// Perform a fast number theory transform in place.
    ///
    /// This function transforms a crt polynomial slice with coefficient in `[0, 4*modulus)`
    /// to a dcrt polynomial slice with coefficient in `[0, 4*modulus)`.
    ///
    /// # Arguments
    ///
    /// * `poly` - inputs in normal order, outputs in bit-reversed order
    fn lazy_transform_slice(&self, poly: &mut [Self::ValueT]);

    /// Perform a fast number theory transform in place.
    ///
    /// This function transforms a polynomial slice with coefficient in `[0, 4*modulus)`
    /// to a ntt polynomial slice with coefficient in `[0, modulus)`.
    ///
    /// # Arguments
    ///
    /// * `poly` - inputs in normal order, outputs in bit-reversed order
    fn transform_slice(&self, poly: &mut [Self::ValueT]);

    /// Perform a fast inverse number theory transform in place.
    ///
    /// This function transforms a ntt polynomial slice with coefficient in `[0, 2*modulus)`
    /// to a polynomial slice with coefficient in `[0, 2*modulus)`.
    ///
    /// # Arguments
    ///
    /// * `values` - inputs in bit-reversed order, outputs in normal order
    fn lazy_inverse_transform_slice(&self, poly: &mut [Self::ValueT]);

    /// Perform a fast inverse number theory transform in place.
    ///
    /// This function transforms a ntt polynomial slice with coefficient in `[0, 2*modulus)`
    /// to a polynomial slice with coefficient in `[0, modulus)`.
    ///
    /// # Arguments
    ///
    /// * `values` - inputs in bit-reversed order, outputs in normal order
    fn inverse_transform_slice(&self, poly: &mut [Self::ValueT]);
}
