mod coeff;
mod ntt;

mod crt;
mod dcrt;

pub use coeff::Glev;
pub use ntt::NttGlev;

pub use crt::CrtGlev;
pub use dcrt::DcrtGlev;

#[derive(Debug, Clone, Copy)]
pub struct GlevInfo {
    pub(crate) decompose_length: usize,
    pub(crate) dimension: usize,
    pub(crate) poly_length: usize,
}

impl GlevInfo {
    /// Creates a new [`GlevInfo`].
    #[inline]
    pub fn new(decompose_length: usize, dimension: usize, poly_length: usize) -> Self {
        Self {
            decompose_length,
            dimension,
            poly_length,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CrtGlevInfo {
    pub(crate) moduli_count: usize,
    pub(crate) decompose_length: usize,
    pub(crate) dimension: usize,
    pub(crate) poly_length: usize,
    pub(crate) glev_len: usize,
}

impl CrtGlevInfo {
    /// Creates a new [`CrtGlevInfo`].
    #[inline]
    pub fn new(
        moduli_count: usize,
        decompose_length: usize,
        dimension: usize,
        poly_length: usize,
    ) -> Self {
        let glev_len = decompose_length * (dimension + 1) * poly_length;
        Self {
            moduli_count,
            decompose_length,
            dimension,
            poly_length,
            glev_len,
        }
    }
}
