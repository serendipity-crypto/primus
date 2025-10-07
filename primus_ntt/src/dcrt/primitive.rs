use primus_integer::UnsignedInteger;
use primus_poly::{NttPolynomial, Polynomial, crt::CrtPolynomial, dcrt::DcrtPolynomial};
use primus_reduce::FieldContext;

use crate::{Dcrt, Ntt, NttError, NttTable, UintNttTable};

use super::DcrtTable;

pub struct UintCrtNttTable<T: UnsignedInteger> {
    ntt_tables: Vec<UintNttTable<T>>,
}

impl<T: UnsignedInteger> UintCrtNttTable<T> {}

impl<T: UnsignedInteger> DcrtTable for UintCrtNttTable<T> {
    type ValueT = T;

    type NttTables = UintNttTable<T>;

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
    fn ntt_tables(&self) -> &[Self::NttTables] {
        &self.ntt_tables
    }

    /// Returns an iterator over this [`UintCrtNttTable<T>`].
    #[inline]
    fn iter(&self) -> std::slice::Iter<'_, Self::NttTables> {
        self.ntt_tables.iter()
    }
}

impl<T: UnsignedInteger> Dcrt for UintCrtNttTable<T> {
    #[inline]
    fn transform_inplace(&self, crt_poly: CrtPolynomial<T>) -> DcrtPolynomial<T> {
        let r: Vec<NttPolynomial<T>> = self
            .iter()
            .zip(crt_poly)
            .map(|(t, p)| t.transform_inplace(p))
            .collect();

        DcrtPolynomial::new(r)
    }

    #[inline]
    fn inverse_transform_inplace(&self, dcrt_poly: DcrtPolynomial<T>) -> CrtPolynomial<T> {
        let r: Vec<Polynomial<T>> = self
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
