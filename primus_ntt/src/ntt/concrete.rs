/// ntt for 32bits
pub mod prime32 {

    use concrete_ntt::prime32::Plan;
    use reduce::FieldAdapter;

    use crate::{NttError, NttTable};

    /// Wrapping concrete NTT for 32bit primes.
    pub struct Concrete32Table {
        plan: Plan,
        root: u32,
        modulus: u32,
    }

    impl Concrete32Table {
        /// Create a new NTT table for 32bit prime.
        #[inline]
        pub fn new(modulus: u32, log_n: u32) -> Result<Self, NttError<u32>> {
            let plan = Plan::try_new(1 << log_n, modulus).ok_or(NttError::NttTableErr)?;
            let root = plan.root();

            Ok(Self {
                plan,
                modulus,
                root,
            })
        }

        /// Get the root of unity.
        #[inline]
        pub fn root(&self) -> u32 {
            self.root
        }

        /// Returns the modulus.
        #[inline]
        pub fn modulus(&self) -> u32 {
            self.modulus
        }
    }

    impl NttTable for Concrete32Table {
        type ValueT = u32;

        #[inline]
        fn new<M>(log_n: u32, modulus: M) -> Result<Self, NttError<Self::ValueT>>
        where
            M: FieldAdapter<Self::ValueT>,
        {
            let modulus = modulus.value_unchecked();
            let plan = Plan::try_new(1 << log_n, modulus).ok_or(NttError::NttTableErr)?;
            let root = plan.root();

            Ok(Self {
                plan,
                modulus,
                root,
            })
        }

        #[inline]
        fn poly_length(&self) -> usize {
            self.plan.ntt_size()
        }
    }
}

/// ntt for 64bits
pub mod prime64 {
    use concrete_ntt::prime64::Plan;

    use crate::{NttError, NttTable};

    /// Wrapping concrete NTT for 64bit primes.
    pub struct Concrete64Table {
        plan: Plan,
        root: u64,
        modulus: u64,
    }

    impl Concrete64Table {
        /// Get the root of unity.
        #[inline]
        pub fn root(&self) -> u64 {
            self.root
        }

        /// Returns the modulus.
        #[inline]
        pub fn modulus(&self) -> u64 {
            self.modulus
        }
    }

    impl NttTable for Concrete64Table {
        type ValueT = u64;

        fn new<M>(log_n: u32, modulus: M) -> Result<Self, NttError<Self::ValueT>>
        where
            M: reduce::FieldAdapter<Self::ValueT>,
        {
            let modulus = modulus.value_unchecked();
            let plan = Plan::try_new(1 << log_n, modulus).ok_or(NttError::NttTableErr)?;
            let root = plan.root();

            Ok(Self {
                plan,
                root,
                modulus,
            })
        }

        fn poly_length(&self) -> usize {
            self.plan.ntt_size()
        }
    }
}
