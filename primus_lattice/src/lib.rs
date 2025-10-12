//! Defines some lattice cryptographic structure.

#[macro_use]
mod macros;

mod glwe;
mod lwe;
mod rlwe;

mod glev;
mod rlev;

mod ggsw;
mod rgsw;

pub use glwe::{CrtGlwe, DcrtGlwe, Glwe, NttGlwe};
pub use lwe::{Lwe, MultiMsgLwe};
pub use rlwe::{CrtRlwe, DcrtRlwe, NttRlwe, Rlwe};

pub use glev::{CrtGlev, DcrtGlev, Glev, NttGlev};
pub use rlev::{CrtRlev, DcrtRlev, NttRlev, Rlev};

pub use ggsw::{CrtGgsw, DcrtGgsw, Ggsw, NttGgsw};
pub use rgsw::{CrtRgsw, DcrtRgsw, NttRgsw, Rgsw};
