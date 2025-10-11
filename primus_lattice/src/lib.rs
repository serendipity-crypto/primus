//! Defines some lattice cryptographic structure.

mod glwe;
mod lwe;
mod rlwe;

mod glev;
mod rlev;

// mod ggsw;
// mod rgsw;

pub use glwe::{CrtGlwe, CrtGlweInfo, DcrtGlwe, Glwe, GlweInfo, NttGlwe};
pub use lwe::{Lwe, MultiMsgLwe};
pub use rlwe::{CrtRlwe, CrtRlweInfo, DcrtRlwe, NttRlwe, Rlwe};

pub use glev::{CrtGlev, CrtGlevInfo, DcrtGlev, Glev, GlevInfo, NttGlev};
pub use rlev::{CrtRlev, CrtRlevInfo, DcrtRlev, NttRlev, Rlev, RlevInfo};

// pub use ggsw::{CrtGgsw, DcrtGgsw, Ggsw, NttGgsw};
// pub use rgsw::{CrtRgsw, DcrtRgsw, NttRgsw, Rgsw};

#[derive(Debug, Clone, Copy)]
pub struct ModuliCount(pub usize);
