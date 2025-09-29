use primus_integer::UnsignedInteger;
use reduce::FieldContext;

use crate::{NttError, NttTable, UintNttTable};

use super::CrtNttTable;

pub struct UintCrtNttTable<T: UnsignedInteger> {
    ntt_tables: Vec<UintNttTable<T>>,
}

impl<T: UnsignedInteger> UintCrtNttTable<T> {}

impl<T: UnsignedInteger> CrtNttTable for UintCrtNttTable<T> {
    type ValueT = T;

    type Table = UintNttTable<T>;

    #[inline]
    fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, NttError<Self::ValueT>>
    where
        M: FieldContext<Self::ValueT>,
    {
        let mut ntt_tables = Vec::with_capacity(moduli.len());
        for modulus in moduli {
            ntt_tables.push(UintNttTable::new(log_n, *modulus)?);
        }
        Ok(Self { ntt_tables })
    }

    /// Returns a reference to the ntt tables of this [`UintCrtNttTable<T>`].
    #[inline]
    fn ntt_tables(&self) -> &[Self::Table] {
        &self.ntt_tables
    }

    /// Returns an iterator over this [`UintCrtNttTable<T>`].
    #[inline]
    fn iter(&self) -> std::slice::Iter<'_, Self::Table> {
        self.ntt_tables.iter()
    }
}
