use reduce::FieldAdapter;

use crate::{NttError, root::PrimitiveRoot};

mod primitive;

/// An abstract for ntt table generation.
pub trait NttTable: Sized {
    /// The value type.
    type ValueT: PrimitiveRoot;

    /// Creates a new [`NttTable`].
    fn new<M>(log_n: u32, modulus: M) -> Result<Self, NttError<Self::ValueT>>
    where
        M: FieldAdapter<Self::ValueT>;

    /// Creates a new [`NttTable`] with a given primitive root.
    fn with_root(
        log_n: u32,
        modulus: Self::ValueT,
        root: Self::ValueT,
    ) -> Result<Self, NttError<Self::ValueT>>;

    /// Get the polynomial modulus degree.
    fn dimension(&self) -> usize;
}

pub trait CrtNttTable: Sized {
    /// The value type.
    type ValueT: PrimitiveRoot;
}
