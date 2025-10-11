mod array;

mod big_uint_poly;
pub mod crt;
pub mod dcrt;
mod ntt;
mod poly;

pub use array::{Array, ArrayBase, ArrayMut, ArrayRef, Data, DataMut, DataOwned, RawData};

pub use big_uint_poly::BigUintPolynomial;

pub use ntt::{NttPolynomial, NttPolynomialMut, NttPolynomialOwned, NttPolynomialRef};
pub use poly::{Polynomial, PolynomialMut, PolynomialOwned, PolynomialRef};

#[derive(Debug, Clone, Copy)]
pub struct PolyLength(pub usize);
