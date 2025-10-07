use primus_integer::UnsignedInteger;
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::dcrt::DcrtPolynomial;
use primus_reduce::FieldContext;
use serde::{Deserialize, Serialize};

use crate::CrtRlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct DcrtRlwe<T: UnsignedInteger> {
    pub(crate) a: DcrtPolynomial<T>,
    pub(crate) b: DcrtPolynomial<T>,
}

impl<T: UnsignedInteger> DcrtRlwe<T> {
    /// Creates a new [`DcrtRlwe<T>`].
    #[inline]
    pub fn new(a: DcrtPolynomial<T>, b: DcrtPolynomial<T>) -> Self {
        Self { a, b }
    }

    pub fn a(&self) -> &DcrtPolynomial<T> {
        &self.a
    }

    pub fn a_mut(&mut self) -> &mut DcrtPolynomial<T> {
        &mut self.a
    }

    pub fn b(&self) -> &DcrtPolynomial<T> {
        &self.b
    }

    pub fn b_mut(&mut self) -> &mut DcrtPolynomial<T> {
        &mut self.b
    }

    pub fn a_b_mut(&mut self) -> (&mut DcrtPolynomial<T>, &mut DcrtPolynomial<T>) {
        (&mut self.a, &mut self.b)
    }
}

impl<T: UnsignedInteger> DcrtRlwe<T> {
    /// ntt transform
    #[inline]
    pub fn into_coeff_form<Table>(self, table: &Table) -> CrtRlwe<T>
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let Self { a, b } = self;

        let a = table.inverse_transform_inplace(a);
        let b = table.inverse_transform_inplace(b);

        CrtRlwe::new(a, b)
    }

    /// ntt transform
    #[inline]
    pub fn to_coeff_form_inplace<Table>(&self, table: &Table, result: &mut CrtRlwe<T>)
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let (a, b) = result.a_b_mut();

        a.copy_from(&self.a);
        b.copy_from(&self.b);

        table.inverse_transform_slice(a.as_mut());
        table.inverse_transform_slice(b.as_mut());
    }
}

impl<T: UnsignedInteger> DcrtRlwe<T> {
    /// Perform element-wise modular addition of two [`DcrtRlwe<T>`].
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

    /// Perform element-wise modular subtraction of two [`DcrtRlwe<T>`].
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
    /// on the `self` [`DcrtRlwe<T>`] with another `rhs` [`DcrtRlwe<T>`].
    #[inline]
    pub fn add_assign_element_wise<M>(&mut self, rhs: &Self, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        self.a.add_assign(rhs.a(), moduli);
        self.b.add_assign(rhs.b(), moduli);
    }

    /// Performs an in-place element-wise modular subtraction
    /// on the `self` [`DcrtRlwe<T>`] with another `rhs` [`DcrtRlwe<T>`].
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
}
