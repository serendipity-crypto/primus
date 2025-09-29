/// ntt for 32bits
pub mod prime32 {

    use crate::{Concrete32Table, CrtNttTable, NttTable};

    /// Wrapping crt concrete NTT for 32bit primes.
    pub struct CrtConcrete32Table {
        tables: Vec<Concrete32Table>,
    }

    impl CrtNttTable for CrtConcrete32Table {
        type ValueT = u32;

        type Table = Concrete32Table;

        #[inline]
        fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, crate::NttError<Self::ValueT>>
        where
            M: primus_reduce::FieldContext<Self::ValueT>,
        {
            let mut tables = Vec::with_capacity(moduli.len());
            for modulus in moduli {
                tables.push(Self::Table::new(log_n, *modulus)?);
            }
            Ok(Self { tables })
        }

        #[inline]
        fn ntt_tables(&self) -> &[Self::Table] {
            &self.tables
        }

        #[inline]
        fn iter(&self) -> std::slice::Iter<'_, Self::Table> {
            self.tables.iter()
        }
    }
}

/// ntt for 64bits
pub mod prime64 {

    use crate::{Concrete64Table, CrtNttTable, NttTable};

    /// Wrapping crt concrete NTT for 64bit primes.
    pub struct CrtConcrete64Table {
        tables: Vec<Concrete64Table>,
    }

    impl CrtNttTable for CrtConcrete64Table {
        type ValueT = u64;

        type Table = Concrete64Table;

        #[inline]
        fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, crate::NttError<Self::ValueT>>
        where
            M: primus_reduce::FieldContext<Self::ValueT>,
        {
            let mut tables = Vec::with_capacity(moduli.len());
            for modulus in moduli {
                tables.push(Self::Table::new(log_n, *modulus)?);
            }
            Ok(Self { tables })
        }

        #[inline]
        fn ntt_tables(&self) -> &[Self::Table] {
            &self.tables
        }

        #[inline]
        fn iter(&self) -> std::slice::Iter<'_, Self::Table> {
            self.tables.iter()
        }
    }
}
