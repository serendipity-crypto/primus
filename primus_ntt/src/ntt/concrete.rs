/// ntt for 32bits
pub mod prime32 {

    use concrete_ntt::prime32::Plan;
    use primus_poly::{NttPolynomial, Polynomial};
    use reduce::FieldContext;

    use crate::{Ntt, NttError, NttTable};

    /// Wrapping concrete NTT for 32bit primes.
    pub struct Concrete32Table {
        plan: Plan,
        root: u32,
        modulus: u32,
    }

    impl Concrete32Table {
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
            M: FieldContext<Self::ValueT>,
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

    impl Ntt for Concrete32Table {
        type CoeffPoly = Polynomial<u32>;

        type NttPoly = NttPolynomial<u32>;

        #[inline]
        fn transform_inplace(&self, mut poly: Self::CoeffPoly) -> Self::NttPoly {
            self.transform_slice(poly.as_mut_slice());
            Self::NttPoly::new(poly.into_vec())
        }

        #[inline]
        fn inverse_transform_inplace(&self, mut values: Self::NttPoly) -> Self::CoeffPoly {
            self.inverse_transform_slice(values.as_mut_slice());
            Self::CoeffPoly::new(values.into_vec())
        }

        #[inline]
        fn lazy_transform_slice(&self, poly: &mut [<Self as NttTable>::ValueT]) {
            self.plan.fwd(poly);
        }

        #[inline]
        fn transform_slice(&self, poly: &mut [<Self as NttTable>::ValueT]) {
            self.plan.fwd(poly);
        }

        #[inline]
        fn lazy_inverse_transform_slice(&self, values: &mut [<Self as NttTable>::ValueT]) {
            self.plan.inv(values);
            self.plan.normalize(values);
        }

        #[inline]
        fn inverse_transform_slice(&self, values: &mut [<Self as NttTable>::ValueT]) {
            self.plan.inv(values);
            self.plan.normalize(values);
        }

        #[inline]
        fn transform_monomial(
            &self,
            coeff: Self::ValueT,
            degree: usize,
            values: &mut [<Self as NttTable>::ValueT],
        ) {
            self.plan.fwd_monomial(coeff, degree, values);
        }

        #[inline]
        fn transform_coeff_one_monomial(
            &self,
            degree: usize,
            values: &mut [<Self as NttTable>::ValueT],
        ) {
            self.plan.fwd_coeff_one_monomial(degree, values);
        }

        #[inline]
        fn transform_coeff_minus_one_monomial(
            &self,
            degree: usize,
            values: &mut [<Self as NttTable>::ValueT],
        ) {
            self.plan.fwd_coeff_minus_one_monomial(degree, values);
        }
    }
}

/// ntt for 64bits
pub mod prime64 {
    use concrete_ntt::prime64::Plan;
    use primus_poly::{NttPolynomial, Polynomial};

    use crate::{Ntt, NttError, NttTable};

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
            M: reduce::FieldContext<Self::ValueT>,
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

    impl Ntt for Concrete64Table {
        type CoeffPoly = Polynomial<u64>;

        type NttPoly = NttPolynomial<u64>;

        #[inline]
        fn transform_inplace(&self, mut poly: Self::CoeffPoly) -> Self::NttPoly {
            self.transform_slice(poly.as_mut_slice());
            Self::NttPoly::new(poly.into_vec())
        }

        #[inline]
        fn inverse_transform_inplace(&self, mut values: Self::NttPoly) -> Self::CoeffPoly {
            self.inverse_transform_slice(values.as_mut_slice());
            Self::CoeffPoly::new(values.into_vec())
        }

        #[inline]
        fn lazy_transform_slice(&self, poly: &mut [<Self as NttTable>::ValueT]) {
            self.plan.fwd(poly);
        }

        #[inline]
        fn transform_slice(&self, poly: &mut [<Self as NttTable>::ValueT]) {
            self.plan.fwd(poly);
        }

        #[inline]
        fn lazy_inverse_transform_slice(&self, values: &mut [<Self as NttTable>::ValueT]) {
            self.plan.inv(values);
            self.plan.normalize(values);
        }

        #[inline]
        fn inverse_transform_slice(&self, values: &mut [<Self as NttTable>::ValueT]) {
            self.plan.inv(values);
            self.plan.normalize(values);
        }

        #[inline]
        fn transform_monomial(
            &self,
            coeff: Self::ValueT,
            degree: usize,
            values: &mut [<Self as NttTable>::ValueT],
        ) {
            self.plan.fwd_monomial(coeff, degree, values);
        }

        #[inline]
        fn transform_coeff_one_monomial(
            &self,
            degree: usize,
            values: &mut [<Self as NttTable>::ValueT],
        ) {
            self.plan.fwd_coeff_one_monomial(degree, values);
        }

        #[inline]
        fn transform_coeff_minus_one_monomial(
            &self,
            degree: usize,
            values: &mut [<Self as NttTable>::ValueT],
        ) {
            self.plan.fwd_coeff_minus_one_monomial(degree, values);
        }
    }
}
