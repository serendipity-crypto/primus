use primus_integer::UnsignedInteger;
use primus_ntt::{Ntt, NttTable};
use primus_poly::{NttPolynomial, Polynomial, crt::CrtPolynomial, dcrt::DcrtPolynomial};
use serde::{Deserialize, Serialize};

use crate::DcrtRlwe;

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Serialize, Deserialize)]
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

    // /// ntt transform
    // #[inline]
    // pub fn into_ntt_form<Table>(self, ntt_table: &Table) -> DcrtRlwe<T>
    // where
    //     Table: NttTable<ValueT = T> + Ntt<CoeffPoly = Polynomial<T>, NttPoly = NttPolynomial<T>>,
    // {
    //     let Self { a, b } = self;

    //     let a = DcrtPolynomial::new(
    //         a.into_iter()
    //             .map(|p| ntt_table.transform_inplace(p))
    //             .collect(),
    //     );
    //     let b = DcrtPolynomial::new(
    //         b.into_iter()
    //             .map(|p| ntt_table.transform_inplace(p))
    //             .collect(),
    //     );

    //     DcrtRlwe::new(a, b)
    // }

    // /// ntt transform
    // #[inline]
    // pub fn transform_inplace<Table>(&self, ntt_table: &Table, result: &mut DcrtRlwe<T>)
    // where
    //     Table: NttTable<ValueT = T> + Ntt<CoeffPoly = Polynomial<T>, NttPoly = NttPolynomial<T>>,
    // {
    //     let (a, b) = result.a_b_mut();

    //     a.iter_mut().zip(self.a.iter()).for_each(|(x, y)| {
    //         x.copy_from(y);
    //         ntt_table.transform_slice(x.as_mut_slice());
    //     });
    //     b.iter_mut().zip(self.b.iter()).for_each(|(x, y)| {
    //         x.copy_from(y);
    //         ntt_table.transform_slice(x.as_mut_slice());
    //     });
    // }
}
