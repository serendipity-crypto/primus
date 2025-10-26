use primus_integer::UnsignedInteger;
use primus_poly::{DataMut, RawData, crt::CrtPolynomial, dcrt::DcrtPolynomial};
use primus_reduce::FieldContext;

use crate::{Dcrt, Ntt, NttError, NttTable, UintNttTable};

use super::DcrtTable;

pub struct UintCrtNttTable<T: UnsignedInteger> {
    ntt_tables: Vec<UintNttTable<T>>,
    poly_length: usize,
    moduli_count: usize,
    crt_poly_length: usize,
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
        let moduli_count = moduli.len();
        let poly_length = 1 << log_n;
        let crt_poly_length = moduli_count * poly_length;

        let mut ntt_tables = Vec::with_capacity(moduli_count);
        for modulus in moduli {
            ntt_tables.push(UintNttTable::new(log_n, *modulus)?);
        }

        Ok(Self {
            ntt_tables,
            poly_length,
            moduli_count,
            crt_poly_length,
        })
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

    #[inline]
    fn poly_length(&self) -> usize {
        self.poly_length
    }

    #[inline]
    fn moduli_count(&self) -> usize {
        self.moduli_count
    }

    #[inline]
    fn crt_poly_length(&self) -> usize {
        self.crt_poly_length
    }
}

impl<T: UnsignedInteger> Dcrt for UintCrtNttTable<T> {
    #[inline]
    fn transform_inplace<S: RawData<Elem = Self::ValueT> + DataMut>(
        &self,
        mut crt_poly: CrtPolynomial<S, T>,
    ) -> DcrtPolynomial<S, T> {
        let poly_length = self.poly_length();

        debug_assert_eq!(poly_length * self.moduli_count, crt_poly.crt_poly_length());

        self.iter()
            .zip(crt_poly.iter_each_modulus_mut(poly_length))
            .for_each(|(t, p)| t.transform_slice(p));

        DcrtPolynomial::new(crt_poly.0)
    }

    #[inline]
    fn inverse_transform_inplace<S: RawData<Elem = Self::ValueT> + DataMut>(
        &self,
        mut dcrt_poly: DcrtPolynomial<S, T>,
    ) -> CrtPolynomial<S, T> {
        let poly_length = self.poly_length();

        debug_assert_eq!(poly_length * self.moduli_count, dcrt_poly.dcrt_poly_length());

        self.iter()
            .zip(dcrt_poly.iter_each_modulus_mut(poly_length))
            .for_each(|(t, p)| t.inverse_transform_slice(p));

        CrtPolynomial::new(dcrt_poly.0)
    }

    #[inline]
    fn lazy_transform_slice(&self, poly: &mut [Self::ValueT]) {
        let poly_length = self.poly_length();
        debug_assert_eq!(poly_length * self.moduli_count, poly.len());
        self.iter()
            .zip(poly.chunks_exact_mut(poly_length))
            .for_each(|(ntt_table, poly)| ntt_table.lazy_transform_slice(poly))
    }

    #[inline]
    fn transform_slice(&self, poly: &mut [Self::ValueT]) {
        let poly_length = self.poly_length();
        debug_assert_eq!(poly_length * self.moduli_count, poly.len());
        self.iter()
            .zip(poly.chunks_exact_mut(poly_length))
            .for_each(|(ntt_table, poly)| ntt_table.transform_slice(poly))
    }

    #[inline]
    fn lazy_inverse_transform_slice(&self, poly: &mut [Self::ValueT]) {
        let poly_length = self.poly_length();
        debug_assert_eq!(poly_length * self.moduli_count, poly.len());
        self.iter()
            .zip(poly.chunks_exact_mut(poly_length))
            .for_each(|(ntt_table, values)| ntt_table.lazy_inverse_transform_slice(values))
    }

    #[inline]
    fn inverse_transform_slice(&self, poly: &mut [Self::ValueT]) {
        let poly_length = self.poly_length();
        debug_assert_eq!(poly_length * self.moduli_count, poly.len());
        self.iter()
            .zip(poly.chunks_exact_mut(poly_length))
            .for_each(|(ntt_table, values)| ntt_table.inverse_transform_slice(values))
    }
}
