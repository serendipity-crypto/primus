use primus_decompose::primitive::ApproxSignedBasis;
use primus_distr::DiscreteGaussian;
use primus_integer::UnsignedInteger;
use primus_reduce::FieldContext;

use crate::RingSecretKeyType;

/// Rlwe Parameters.
#[derive(Clone)]
pub struct RlweParameters<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> {
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
    /// The noise's distribution.
    pub noise_distribution: DiscreteGaussian<ValueT>,
}

impl<ValueT, ModulusT> RlweParameters<ValueT, ModulusT>
where
    ValueT: UnsignedInteger,
    ModulusT: FieldContext<ValueT>,
{
    /// Returns the poly length of this [`RlweParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the cipher modulus of this [`RlweParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn cipher_modulus(&self) -> ValueT {
        self.modulus.value_unchecked()
    }

    /// Returns the noise distribution.
    #[inline]
    pub fn noise_distribution(&self) -> &DiscreteGaussian<ValueT> {
        // DiscreteGaussian::new(0.0, self.noise_standard_deviation, self.modulus_minus_one).unwrap()
        &self.noise_distribution
    }

    /// Returns the noise distribution.
    #[inline]
    pub fn noise_distribution_div_count(&self, count: u32) -> DiscreteGaussian<ValueT> {
        let var = self.noise_standard_deviation * self.noise_standard_deviation;
        let sigma = (var / count as f64).sqrt();
        DiscreteGaussian::new(0.0, sigma, self.modulus_minus_one).unwrap()
    }
}

/// Rlev Parameters.
#[derive(Clone)]
pub struct RlevParameters<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> {
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
    /// The noise's distribution.
    pub noise_distribution: DiscreteGaussian<ValueT>,
}

impl<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> RlevParameters<ValueT, ModulusT> {
    /// Returns the decompose basis.
    #[inline]
    pub fn basis(&self) -> &ApproxSignedBasis<ValueT> {
        &self.basis
    }

    /// Returns the poly length of this [`GadgetRlweParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the noise distribution.
    #[inline]
    pub fn noise_distribution(&self) -> &DiscreteGaussian<ValueT> {
        // DiscreteGaussian::new(0.0, self.noise_standard_deviation, self.modulus_minus_one).unwrap()
        &self.noise_distribution
    }

    /// Returns the noise distribution.
    #[inline]
    pub fn noise_distribution_div_count(&self, count: u32) -> DiscreteGaussian<ValueT> {
        let var = self.noise_standard_deviation * self.noise_standard_deviation;
        let sigma = (var / count as f64).sqrt();
        DiscreteGaussian::new(0.0, sigma, self.modulus_minus_one).unwrap()
    }
}

/// Rgsw Parameters.
pub type RgswParameters<ValueT, ModulusT> = RlevParameters<ValueT, ModulusT>;
