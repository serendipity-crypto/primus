//! Defines some lattice cryptographic structure.

mod lwe;
mod rlwe;

pub use lwe::{ExLwe, Lwe};
pub use rlwe::{NttRlwe, Rlwe};
