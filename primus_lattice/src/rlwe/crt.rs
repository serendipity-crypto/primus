use primus_integer::{UnsignedInteger, size::Size};
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::{crt::CrtPolynomial, dcrt::DcrtPolynomial};
use primus_reduce::FieldContext;
use serde::{Deserialize, Serialize};

use crate::DcrtRlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct CrtRlwe<T: UnsignedInteger> {
    pub(crate) a: CrtPolynomial<T>,
    pub(crate) b: CrtPolynomial<T>,
}

impl<T: UnsignedInteger> CrtRlwe<T> {
    /// Creates a new [`CrtRlwe<T>`].
    #[inline]
    pub fn new(a: CrtPolynomial<T>, b: CrtPolynomial<T>) -> Self {
        Self { a, b }
    }

    /// Creates a [`CrtRlwe<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(moduli_count: usize, poly_length: usize) -> Self {
        Self {
            a: CrtPolynomial::zero(moduli_count, poly_length),
            b: CrtPolynomial::zero(moduli_count, poly_length),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.set_zero();
        self.b.set_zero();
    }

    pub fn a(&self) -> &CrtPolynomial<T> {
        &self.a
    }

    pub fn a_mut(&mut self) -> &mut CrtPolynomial<T> {
        &mut self.a
    }

    pub fn b(&self) -> &CrtPolynomial<T> {
        &self.b
    }

    pub fn b_mut(&mut self) -> &mut CrtPolynomial<T> {
        &mut self.b
    }

    pub fn a_b_mut(&mut self) -> (&mut CrtPolynomial<T>, &mut CrtPolynomial<T>) {
        (&mut self.a, &mut self.b)
    }
}

impl<T: UnsignedInteger> CrtRlwe<T> {
    /// ntt transform
    #[inline]
    pub fn into_ntt_form<Table>(self, table: &Table) -> DcrtRlwe<T>
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let Self { a, b } = self;

        let a = table.transform_inplace(a);
        let b = table.transform_inplace(b);

        DcrtRlwe::new(a, b)
    }

    /// ntt transform
    #[inline]
    pub fn to_ntt_form_inplace<Table>(&self, table: &Table, result: &mut DcrtRlwe<T>)
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let (a, b) = result.a_b_mut();

        a.copy_from(&self.a);
        b.copy_from(&self.b);

        table.transform_slice(a.as_mut());
        table.transform_slice(b.as_mut());
    }
}

impl<T: UnsignedInteger> CrtRlwe<T> {
    /// Perform element-wise modular addition of two [`CrtRlwe<T>`].
    #[inline]
    pub fn add_element_wise<M>(self, rhs: &Self, moduli: &[M]) -> Self
    where
        M: FieldContext<T>,
    {
        Self {
            a: self.a.add(rhs.a(), moduli),
            b: self.b.add(rhs.b(), moduli),
        }
    }

    /// Perform element-wise modular subtraction of two [`CrtRlwe<T>`].
    #[inline]
    pub fn sub_element_wise<M>(self, rhs: &Self, moduli: &[M]) -> Self
    where
        M: FieldContext<T>,
    {
        Self {
            a: self.a.sub(rhs.a(), moduli),
            b: self.b.sub(rhs.b(), moduli),
        }
    }

    /// Performs an in-place element-wise modular addition
    /// on the `self` [`CrtRlwe<T>`] with another `rhs` [`CrtRlwe<T>`].
    #[inline]
    pub fn add_assign_element_wise<M>(&mut self, rhs: &Self, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        self.a.add_assign(rhs.a(), moduli);
        self.b.add_assign(rhs.b(), moduli);
    }

    /// Performs an in-place element-wise modular subtraction
    /// on the `self` [`CrtRlwe<T>`] with another `rhs` [`CrtRlwe<T>`].
    #[inline]
    pub fn sub_assign_element_wise<M>(&mut self, rhs: &Self, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        self.a.sub_assign(rhs.a(), moduli);
        self.b.sub_assign(rhs.b(), moduli);
    }

    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `result`.
    #[inline]
    pub fn add_inplace<M>(&self, rhs: &Self, result: &mut Self, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        self.a.add_inplace(rhs.a(), result.a_mut(), moduli);
        self.b.add_inplace(rhs.b(), result.b_mut(), moduli);
    }

    /// Performs subtraction operation:`self - rhs`,
    /// and put the result to the `result`.
    #[inline]
    pub fn sub_inplace<M>(&self, rhs: &Self, result: &mut Self, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        self.a.sub_inplace(rhs.a(), result.a_mut(), moduli);
        self.b.sub_inplace(rhs.b(), result.b_mut(), moduli);
    }

    /// Performs a multiplication on the `self` [`CrtRlwe<T>`] with another `dcrt_polynomial` [`DcrtPolynomial<T>`],
    /// store the result into `result` [`DcrtRlwe<T>`].
    #[inline]
    pub fn mul_dcrt_polynomial_inplace<M, Table>(
        &self,
        dcrt_polynomial: &DcrtPolynomial<T>,
        result: &mut DcrtRlwe<T>,
        moduli: &[M],
        table: &Table,
    ) where
        M: FieldContext<T>,
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let (a, b) = result.a_b_mut();

        a.copy_from(self.a());
        b.copy_from(self.b());

        table.transform_slice(a.as_mut());
        table.transform_slice(b.as_mut());

        a.mul_assign(dcrt_polynomial, moduli);
        b.mul_assign(dcrt_polynomial, moduli);
    }
}

impl<T: UnsignedInteger> Size for CrtRlwe<T> {
    #[inline]
    fn byte_count(&self) -> usize {
        self.a.byte_count() * 2
    }
}
