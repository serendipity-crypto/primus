use reduce::FieldAdapter;

use crate::{NttError, PrimitiveRoot};

#[cfg(feature = "concrete-ntt")]
mod concrete;
mod primitive;

pub trait CrtNttTable: Sized {
    /// The value type.
    type ValueT: PrimitiveRoot;

    /// Creates a new [`CrtNttTable`].
    fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, NttError<Self::ValueT>>
    where
        M: FieldAdapter<Self::ValueT>;
}
