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

pub trait CrtNttTable: Sized {
    /// The value type.
    type ValueT: PrimitiveRoot;

    type Table: NttTable<ValueT = Self::ValueT>;

    /// Creates a new [`CrtNttTable`].
    fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, NttError<Self::ValueT>>
    where
        M: FieldContext<Self::ValueT>;

    /// Returns a reference to the ntt tables.
    fn ntt_tables(&self) -> &[Self::Table];

    /// Returns an iterator over the ntt tables.
    fn iter(&self) -> std::slice::Iter<'_, Self::Table>;
}
