use integer::UnsignedInteger;
use reduce::FieldAdapter;

use crate::{NttTable, UintNttTable};

use super::CrtNttTable;

pub struct UintCrtNttTable<T: UnsignedInteger> {
    ntt_tables: Vec<UintNttTable<T>>,
}

impl<T: UnsignedInteger> UintCrtNttTable<T> {}

impl<T: UnsignedInteger> CrtNttTable for UintCrtNttTable<T> {
    type ValueT = T;

    fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, crate::NttError<Self::ValueT>>
    where
        M: FieldAdapter<Self::ValueT>,
    {
        let mut ntt_tables = Vec::with_capacity(moduli.len());
        for modulus in moduli {
            ntt_tables.push(UintNttTable::new(log_n, *modulus)?);
        }
        Ok(Self { ntt_tables })
    }
}
