/// ntt for 32bits
pub mod prime32 {

    use primus_poly::{NttPolynomial, Polynomial, crt::CrtPolynomial, dcrt::DcrtPolynomial};

    use crate::{Concrete32Table, Dcrt, DcrtTable, Ntt, NttTable};

    /// Wrapping crt concrete NTT for 32bit primes.
    pub struct CrtConcrete32Table {
        ntt_tables: Vec<Concrete32Table>,
    }

    impl DcrtTable for CrtConcrete32Table {
        type ValueT = u32;

        type NttTables = Concrete32Table;

        #[inline]
        fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, crate::NttError<Self::ValueT>>
        where
            M: primus_reduce::FieldContext<Self::ValueT>,
        {
            let mut ntt_tables = Vec::with_capacity(moduli.len());
            for modulus in moduli {
                ntt_tables.push(Self::NttTables::new(log_n, *modulus)?);
            }
            Ok(Self { ntt_tables })
        }

        #[inline]
        fn ntt_tables(&self) -> &[Self::NttTables] {
            &self.ntt_tables
        }

        #[inline]
        fn iter(&self) -> std::slice::Iter<'_, Self::NttTables> {
            self.ntt_tables.iter()
        }
    }

    impl Dcrt for CrtConcrete32Table {
        #[inline]
        fn transform_inplace(&self, crt_poly: CrtPolynomial<u32>) -> DcrtPolynomial<u32> {
            let r: Vec<NttPolynomial<u32>> = self
                .iter()
                .zip(crt_poly)
                .map(|(t, p)| t.transform_inplace(p))
                .collect();

            DcrtPolynomial::new(r)
        }

        #[inline]
        fn inverse_transform_inplace(&self, dcrt_poly: DcrtPolynomial<u32>) -> CrtPolynomial<u32> {
            let r: Vec<Polynomial<u32>> = self
                .iter()
                .zip(dcrt_poly)
                .map(|(t, p)| t.inverse_transform_inplace(p))
                .collect();

            CrtPolynomial::new(r)
        }

        #[inline]
        fn lazy_transform_slice<P: AsMut<[Self::ValueT]>>(&self, poly: &mut [P]) {
            self.iter()
                .zip(poly)
                .for_each(|(ntt_table, s)| ntt_table.lazy_transform_slice(s.as_mut()))
        }

        #[inline]
        fn transform_slice<P: AsMut<[Self::ValueT]>>(&self, poly: &mut [P]) {
            self.iter()
                .zip(poly)
                .for_each(|(ntt_table, s)| ntt_table.transform_slice(s.as_mut()))
        }

        #[inline]
        fn lazy_inverse_transform_slice<P: AsMut<[Self::ValueT]>>(&self, poly: &mut [P]) {
            self.iter()
                .zip(poly)
                .for_each(|(ntt_table, s)| ntt_table.lazy_inverse_transform_slice(s.as_mut()))
        }

        #[inline]
        fn inverse_transform_slice<P: AsMut<[Self::ValueT]>>(&self, poly: &mut [P]) {
            self.iter()
                .zip(poly)
                .for_each(|(ntt_table, s)| ntt_table.inverse_transform_slice(s.as_mut()))
        }
    }
}

/// ntt for 64bits
pub mod prime64 {

    use primus_poly::{NttPolynomial, Polynomial, crt::CrtPolynomial, dcrt::DcrtPolynomial};

    use crate::{Concrete64Table, Dcrt, DcrtTable, Ntt, NttTable};

    /// Wrapping crt concrete NTT for 64bit primes.
    pub struct CrtConcrete64Table {
        ntt_tables: Vec<Concrete64Table>,
    }

    impl DcrtTable for CrtConcrete64Table {
        type ValueT = u64;

        type NttTables = Concrete64Table;

        #[inline]
        fn new<M>(log_n: u32, moduli: &[M]) -> Result<Self, crate::NttError<Self::ValueT>>
        where
            M: primus_reduce::FieldContext<Self::ValueT>,
        {
            let mut ntt_tables = Vec::with_capacity(moduli.len());
            for modulus in moduli {
                ntt_tables.push(Self::NttTables::new(log_n, *modulus)?);
            }
            Ok(Self { ntt_tables })
        }

        #[inline]
        fn ntt_tables(&self) -> &[Self::NttTables] {
            &self.ntt_tables
        }

        #[inline]
        fn iter(&self) -> std::slice::Iter<'_, Self::NttTables> {
            self.ntt_tables.iter()
        }
    }

    impl Dcrt for CrtConcrete64Table {
        #[inline]
        fn transform_inplace(&self, crt_poly: CrtPolynomial<u64>) -> DcrtPolynomial<u64> {
            let r: Vec<NttPolynomial<u64>> = self
                .iter()
                .zip(crt_poly)
                .map(|(t, p)| t.transform_inplace(p))
                .collect();

            DcrtPolynomial::new(r)
        }

        #[inline]
        fn inverse_transform_inplace(&self, dcrt_poly: DcrtPolynomial<u64>) -> CrtPolynomial<u64> {
            let r: Vec<Polynomial<u64>> = self
                .iter()
                .zip(dcrt_poly)
                .map(|(t, p)| t.inverse_transform_inplace(p))
                .collect();

            CrtPolynomial::new(r)
        }

        #[inline]
        fn lazy_transform_slice<P: AsMut<[Self::ValueT]>>(&self, poly: &mut [P]) {
            self.iter()
                .zip(poly)
                .for_each(|(ntt_table, s)| ntt_table.lazy_transform_slice(s.as_mut()))
        }

        #[inline]
        fn transform_slice<P: AsMut<[Self::ValueT]>>(&self, poly: &mut [P]) {
            self.iter()
                .zip(poly)
                .for_each(|(ntt_table, s)| ntt_table.transform_slice(s.as_mut()))
        }

        #[inline]
        fn lazy_inverse_transform_slice<P: AsMut<[Self::ValueT]>>(&self, poly: &mut [P]) {
            self.iter()
                .zip(poly)
                .for_each(|(ntt_table, s)| ntt_table.lazy_inverse_transform_slice(s.as_mut()))
        }

        #[inline]
        fn inverse_transform_slice<P: AsMut<[Self::ValueT]>>(&self, poly: &mut [P]) {
            self.iter()
                .zip(poly)
                .for_each(|(ntt_table, s)| ntt_table.inverse_transform_slice(s.as_mut()))
        }
    }
}
