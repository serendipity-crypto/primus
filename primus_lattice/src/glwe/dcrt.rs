use primus_integer::UnsignedInteger;
use primus_ntt::{Dcrt, DcrtTable};
use primus_poly::dcrt::DcrtPolynomial;
use primus_reduce::FieldContext;
use primus_utils::{Size, izip};
use serde::{Deserialize, Serialize};

use super::CrtGlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Clone, Serialize, Deserialize)]
#[serde(bound(deserialize = "T: UnsignedInteger"))]
pub struct DcrtGlwe<T: UnsignedInteger> {
    pub(crate) a: Vec<DcrtPolynomial<T>>,
    pub(crate) b: DcrtPolynomial<T>,
}

impl<T: UnsignedInteger> DcrtGlwe<T> {
    /// Creates a new [`DcrtGlwe<T>`].
    #[inline]
    pub fn new(a: Vec<DcrtPolynomial<T>>, b: DcrtPolynomial<T>) -> Self {
        Self { a, b }
    }

    /// Creates a [`DcrtGlwe<T>`] with all entries equal to zero.
    #[inline]
    pub fn zero(dimension: usize, moduli_count: usize, poly_length: usize) -> Self {
        Self {
            a: (0..dimension)
                .map(|_| DcrtPolynomial::zero(moduli_count, poly_length))
                .collect(),
            b: DcrtPolynomial::zero(moduli_count, poly_length),
        }
    }

    /// Set all entries equal to zero.
    #[inline]
    pub fn set_zero(&mut self) {
        self.a.iter_mut().for_each(|v| v.set_zero());
        self.b.set_zero();
    }

    pub fn a(&self) -> &[DcrtPolynomial<T>] {
        &self.a
    }

    pub fn a_mut(&mut self) -> &mut [DcrtPolynomial<T>] {
        &mut self.a
    }

    pub fn b(&self) -> &DcrtPolynomial<T> {
        &self.b
    }

    pub fn b_mut(&mut self) -> &mut DcrtPolynomial<T> {
        &mut self.b
    }

    pub fn a_b_mut(&mut self) -> (&mut [DcrtPolynomial<T>], &mut DcrtPolynomial<T>) {
        (&mut self.a, &mut self.b)
    }
}

impl<T: UnsignedInteger> DcrtGlwe<T> {
    /// ntt transform
    #[inline]
    pub fn into_coeff_form<Table>(self, table: &Table) -> CrtGlwe<T>
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let Self { a, b } = self;

        let a = a
            .into_iter()
            .map(|v| table.inverse_transform_inplace(v))
            .collect();
        let b = table.inverse_transform_inplace(b);

        CrtGlwe::new(a, b)
    }

    /// ntt transform
    #[inline]
    pub fn to_coeff_form_inplace<Table>(&self, table: &Table, result: &mut CrtGlwe<T>)
    where
        Table: DcrtTable<ValueT = T> + Dcrt,
    {
        let (a, b) = result.a_b_mut();

        a.iter_mut().zip(&self.a).for_each(|(x, y)| {
            x.copy_from(y);
            table.inverse_transform_slice(x.as_mut());
        });

        b.copy_from(&self.b);
        table.inverse_transform_slice(b.as_mut());
    }
}

impl<T: UnsignedInteger> DcrtGlwe<T> {
    /// Perform element-wise modular addition of two [`DcrtGlwe<T>`].
    #[inline]
    pub fn add_element_wise<M>(self, rhs: &Self, moduli: &[M]) -> Self
    where
        M: FieldContext<T>,
    {
        Self {
            a: self
                .a
                .into_iter()
                .zip(rhs.a())
                .map(|(x, y)| x.add(y, moduli))
                .collect(),
            b: self.b.add(rhs.b(), moduli),
        }
    }

    /// Perform element-wise modular subtraction of two [`DcrtGlwe<T>`].
    #[inline]
    pub fn sub_element_wise<M>(self, rhs: &Self, moduli: &[M]) -> Self
    where
        M: FieldContext<T>,
    {
        Self {
            a: self
                .a
                .into_iter()
                .zip(rhs.a())
                .map(|(x, y)| x.sub(y, moduli))
                .collect(),
            b: self.b.sub(rhs.b(), moduli),
        }
    }

    /// Performs an in-place element-wise modular addition
    /// on the `self` [`DcrtGlwe<T>`] with another `rhs` [`DcrtGlwe<T>`].
    #[inline]
    pub fn add_assign_element_wise<M>(&mut self, rhs: &Self, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        self.a
            .iter_mut()
            .zip(rhs.a())
            .for_each(|(x, y)| x.add_assign(y, moduli));
        self.b.add_assign(rhs.b(), moduli);
    }

    /// Performs an in-place element-wise modular subtraction
    /// on the `self` [`DcrtGlwe<T>`] with another `rhs` [`DcrtGlwe<T>`].
    #[inline]
    pub fn sub_assign_element_wise<M>(&mut self, rhs: &Self, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        self.a
            .iter_mut()
            .zip(rhs.a())
            .for_each(|(x, y)| x.sub_assign(y, moduli));
        self.b.sub_assign(rhs.b(), moduli);
    }

    /// Performs addition operation:`self + rhs`,
    /// and puts the result to the `result`.
    #[inline]
    pub fn add_inplace<M>(&self, rhs: &Self, result: &mut Self, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        izip!(self.a(), rhs.a(), result.a_mut()).for_each(|(x, y, z)| {
            x.add_inplace(y, z, moduli);
        });
        self.b.add_inplace(rhs.b(), result.b_mut(), moduli);
    }

    /// Performs subtraction operation:`self - rhs`,
    /// and put the result to the `result`.
    #[inline]
    pub fn sub_inplace<M>(&self, rhs: &Self, result: &mut Self, moduli: &[M])
    where
        M: FieldContext<T>,
    {
        izip!(self.a(), rhs.a(), result.a_mut()).for_each(|(x, y, z)| {
            x.sub_inplace(y, z, moduli);
        });
        self.b.sub_inplace(rhs.b(), result.b_mut(), moduli);
    }
}

impl<T: UnsignedInteger> Size for DcrtGlwe<T> {
    #[inline]
    fn size(&self) -> usize {
        self.b.size() * (self.a.len() + 1)
    }
}
