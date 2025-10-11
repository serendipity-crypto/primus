mod coeff;
mod ntt;

mod crt;
mod dcrt;

pub use coeff::Rlev;
pub use ntt::NttRlev;

pub use crt::CrtRlev;
pub use dcrt::DcrtRlev;

#[derive(Debug, Clone, Copy)]
pub struct RlevInfo {
    pub(crate) decompose_length: usize,
    pub(crate) poly_length: usize,
}

impl RlevInfo {
    /// Creates a new [`RlevInfo`].
    #[inline]
    pub fn new(decompose_length: usize, poly_length: usize) -> Self {
        Self {
            decompose_length,
            poly_length,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CrtRlevInfo {
    pub(crate) moduli_count: usize,
    pub(crate) decompose_length: usize,
    pub(crate) poly_length: usize,
    pub(crate) rlev_len: usize,
}

impl CrtRlevInfo {
    /// Creates a new [`CrtRlevInfo`].
    #[inline]
    pub fn new(moduli_count: usize, decompose_length: usize, poly_length: usize) -> Self {
        let rlev_len = decompose_length * poly_length * 2;
        Self {
            moduli_count,
            decompose_length,
            poly_length,
            rlev_len,
        }
    }
}
