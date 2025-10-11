mod coeff;
mod ntt;

mod crt;
mod dcrt;

pub use coeff::Rlwe;
pub use ntt::NttRlwe;

pub use crt::CrtRlwe;
pub use dcrt::DcrtRlwe;
use primus_poly::PolyLength;

use crate::ModuliCount;

#[derive(Debug, Clone, Copy)]
pub struct CrtRlweInfo {
    pub(crate) moduli_count: ModuliCount,
    pub(crate) poly_length: PolyLength,
}

impl CrtRlweInfo {
    /// Creates a new [`CrtRlweInfo`].
    #[inline]
    pub fn new(moduli_count: ModuliCount, poly_length: PolyLength) -> Self {
        Self {
            moduli_count,
            poly_length,
        }
    }

    /// Returns the moduli count of this [`CrtRlweInfo`].
    #[inline]
    pub fn moduli_count(&self) -> ModuliCount {
        self.moduli_count
    }

    /// Returns the poly length of this [`CrtRlweInfo`].
    #[inline]
    pub fn poly_length(&self) -> PolyLength {
        self.poly_length
    }
}
