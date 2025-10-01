//! Defines some lattice cryptographic structure.

mod lwe;
mod rlwe;

pub use lwe::{Lwe, MultiMsgLwe};
pub use rlwe::{NttRlwe, Rlwe};
