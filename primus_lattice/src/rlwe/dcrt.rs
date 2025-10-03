use primus_integer::UnsignedInteger;
use primus_poly::dcrt::DcrtPolynomial;
use serde::{Deserialize, Serialize};

/// A cryptographic structure for Ring Learning with Errors (RLWE).
/// This structure is used in advanced cryptographic systems and protocols, particularly
/// those that require efficient homomorphic encryption properties.
#[derive(Serialize, Deserialize)]
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
