use primus_decompose::{big_integer::BigUintApproxSignedBasis, primitive::ApproxSignedBasis};
use primus_distr::DiscreteGaussian;
use primus_integer::UnsignedInteger;
use primus_reduce::FieldContext;

use crate::RingSecretKeyType;

/// Glwe Parameters.
#[derive(Debug, Clone, Copy)]
pub struct GlweParameters<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> {
    /// The dimension, refers to **k** in the paper.
    pub dimension: usize,
    /// The polynomial length, refers to **N** in the paper.
    pub poly_length: usize,
    /// **RLWE** message modulus, refers to **t** in the paper.
    pub plain_modulus_value: ValueT,
    /// **RLWE** cipher modulus minus one, refers to **Q-1**.
    pub modulus_minus_one: ValueT,
    /// The modulus, refers to **Q** in the paper.
    pub modulus: ModulusT,
    /// The distribution type of the secret key.
    pub secret_key_type: RingSecretKeyType,
    /// The noise error's standard deviation.
    pub noise_standard_deviation: f64,
}

impl<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> GlweParameters<ValueT, ModulusT> {
    /// Returns the dimension of this [`GlweParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the poly length of this [`GlweParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the cipher modulus of this [`GlweParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn cipher_modulus(&self) -> ValueT {
        self.modulus.value_unchecked()
    }

    /// Returns the noise distribution.
    #[inline]
    pub fn noise_distribution(&self) -> DiscreteGaussian<ValueT> {
        DiscreteGaussian::new(0.0, self.noise_standard_deviation, self.modulus_minus_one).unwrap()
    }

    /// Returns the noise distribution.
    #[inline]
    pub fn noise_distribution_div_count(&self, count: u32) -> DiscreteGaussian<ValueT> {
        let var = self.noise_standard_deviation * self.noise_standard_deviation;
        let sigma = (var / count as f64).sqrt();
        DiscreteGaussian::new(0.0, sigma, self.modulus_minus_one).unwrap()
    }
}

/// Glwe Parameters.
#[derive(Debug, Clone)]
pub struct CrtGlweParameters<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> {
    /// The dimension, refers to **k** in the paper.
    pub dimension: usize,
    /// The polynomial length, refers to **N** in the paper.
    pub poly_length: usize,
    /// **RLWE** message modulus, refers to **t** in the paper.
    pub plain_modulus_value: ValueT,
    /// **RLWE** cipher modulus minus one, refers to **Q-1**.
    pub modulus_minus_one: Vec<ValueT>,
    pub modulus: Vec<ValueT>,
    /// The modulus, refers to **Q** in the paper.
    pub moduli: Vec<ModulusT>,
    /// The distribution type of the secret key.
    pub secret_key_type: RingSecretKeyType,
    /// The noise error's standard deviation.
    pub noise_standard_deviation: f64,
}

impl<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> CrtGlweParameters<ValueT, ModulusT> {
    /// Returns the dimension of this [`CrtGlweParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the poly length of this [`CrtGlweParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns a reference to the moduli of this [`CrtGlweParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn moduli(&self) -> &[ModulusT] {
        &self.moduli
    }
}

/// Ggsw Parameters.
#[derive(Debug, Clone, Copy)]
pub struct GlevParameters<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> {
    /// The dimension, refers to **k** in the paper.
    pub dimension: usize,
    /// The dimension, refers to **N** in the paper.
    pub poly_length: usize,
    /// cipher modulus minus one, refers to **Q-1**.
    pub modulus_minus_one: ValueT,
    /// The modulus, refers to **Q** in the paper.
    pub modulus: ModulusT,
    /// The distribution type of the secret key.
    pub secret_key_type: RingSecretKeyType,
    /// The noise error's standard deviation.
    pub noise_standard_deviation: f64,
    /// Decompose basis for `Q`.
    pub basis: ApproxSignedBasis<ValueT>,
}

/// Ggsw Parameters.
#[derive(Debug, Clone)]
pub struct CrtGlevParameters<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> {
    /// The dimension, refers to **k** in the paper.
    pub dimension: usize,
    /// The dimension, refers to **N** in the paper.
    pub poly_length: usize,
    /// cipher modulus minus one, refers to **Q-1**.
    pub modulus_minus_one: Vec<ValueT>,
    /// The modulus, refers to **Q** in the paper.
    pub modulus: Vec<ValueT>,
    pub moduli: Vec<ModulusT>,
    /// The distribution type of the secret key.
    pub secret_key_type: RingSecretKeyType,
    /// The noise error's standard deviation.
    pub noise_standard_deviation: f64,
    /// Decompose basis for `Q`.
    pub basis: BigUintApproxSignedBasis<ValueT>,
}
