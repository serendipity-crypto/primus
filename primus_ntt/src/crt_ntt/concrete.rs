/// ntt for 32bits
pub mod prime32 {

    use crate::{Concrete32Table, CrtNttTable};

    /// Wrapping concrete NTT for 32bit primes.
    pub struct CrtConcrete32Table {
        tables: Vec<Concrete32Table>,
    }

    impl CrtNttTable for CrtConcrete32Table {
        type ValueT = u32;

        fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, crate::NttError<Self::ValueT>>
        where
            M: reduce::FieldAdapter<Self::ValueT>,
        {
            todo!()
        }
    }
}
