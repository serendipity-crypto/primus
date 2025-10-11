mod coeff;
mod ntt;

mod crt;
mod dcrt;

pub use coeff::Glwe;
pub use ntt::NttGlwe;

pub use crt::CrtGlwe;
pub use dcrt::DcrtGlwe;

#[derive(Debug, Clone, Copy)]
pub struct GlweInfo {
    pub(crate) dimension: usize,
    pub(crate) poly_length: usize,
}

impl GlweInfo {
    /// Creates a new [`GlweInfo`].
    #[inline]
    pub fn new(dimension: usize, poly_length: usize) -> Self {
        Self {
            dimension,
            poly_length,
        }
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn poly_length(&self) -> usize {
        self.poly_length
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CrtGlweInfo {
    moduli_count: usize,
    dimension: usize,
    poly_length: usize,
    glwe_len: usize,
}

impl CrtGlweInfo {
    /// Creates a new [`CrtGlweInfo`].
    pub fn new(moduli_count: usize, dimension: usize, poly_length: usize) -> Self {
        let glwe_len = (dimension + 1) * poly_length;
        Self {
            moduli_count,
            dimension,
            poly_length,
            glwe_len,
        }
    }

    pub fn moduli_count(&self) -> usize {
        self.moduli_count
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn poly_length(&self) -> usize {
        self.poly_length
    }
}
