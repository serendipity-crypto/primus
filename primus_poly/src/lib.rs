#[macro_use]
mod macros;

mod array;

mod big_uint_poly;
pub mod crt;
pub mod dcrt;
mod ntt;
mod poly;

pub use array::{Array, ArrayBase, ArrayMut, ArrayRef, Data, DataMut, DataOwned, RawData};

pub use big_uint_poly::{BigUintPolynomial, BigUintPolynomialIter, BigUintPolynomialIterMut};

pub use ntt::{
    NttPolynomial, NttPolynomialIter, NttPolynomialIterMut, NttPolynomialMut, NttPolynomialOwned,
    NttPolynomialRef,
};
pub use poly::{
    Polynomial, PolynomialIter, PolynomialIterMut, PolynomialMut, PolynomialOwned, PolynomialRef,
};
