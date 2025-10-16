use primus_decompose::{big_integer::BigUintApproxSignedBasis, primitive::ApproxSignedBasis};
use primus_distr::{DiscreteGaussian, SignedDiscreteGaussian};
use primus_integer::{UnsignedInteger, multiply_many_values};
use primus_reduce::FieldContext;

use crate::RingSecretKeyType;

/// Glwe Parameters.
#[derive(Clone)]
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

/// Big Unsigned Integer Glwe Parameters.
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

impl<ValueT, ModulusT> CrtGlweParameters<ValueT, ModulusT>
where
    ValueT: UnsignedInteger,
    ModulusT: FieldContext<ValueT>,
{
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

    #[inline]
    pub fn value_len(&self) -> usize {
        self.modulus.len()
    }
}

/// Glev Parameters.
#[derive(Clone)]
pub struct GlevParameters<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> {
    /// The dimension, refers to **k** in the paper.
    dimension: usize,
    /// The dimension, refers to **N** in the paper.
    poly_length: usize,
    /// cipher modulus minus one, refers to **Q-1**.
    modulus_minus_one: ValueT,
    /// The modulus, refers to **Q** in the paper.
    modulus: ModulusT,
    /// The distribution type of the secret key.
    secret_key_type: RingSecretKeyType,
    /// The noise error's standard deviation.
    noise_standard_deviation: f64,
    /// The noise's distribution.
    noise_distribution: DiscreteGaussian<ValueT>,
    /// Decompose basis for `Q`.
    basis: ApproxSignedBasis<ValueT>,
}

impl<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> GlevParameters<ValueT, ModulusT> {
    /// Creates a new [`GlevParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn new(
        dimension: usize,
        poly_length: usize,
        modulus: ModulusT,
        secret_key_type: RingSecretKeyType,
        noise_standard_deviation: f64,
        basis: ApproxSignedBasis<ValueT>,
    ) -> Self {
        let modulus_minus_one = modulus.value_unchecked() - ValueT::ONE;

        let noise_distribution =
            DiscreteGaussian::new(0.0, noise_standard_deviation, modulus_minus_one).unwrap();

        Self {
            dimension,
            poly_length,
            modulus_minus_one,
            modulus,
            secret_key_type,
            noise_standard_deviation,
            basis,
            noise_distribution,
        }
    }

    /// Returns the dimension of this [`GlevParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the poly length of this [`GlevParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the modulus of this [`GlevParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn modulus(&self) -> ModulusT {
        self.modulus
    }

    /// Returns a reference to the noise distribution of this [`GlevParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn noise_distribution(&self) -> &DiscreteGaussian<ValueT> {
        &self.noise_distribution
    }

    /// Returns the basis of this [`GlevParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn basis(&self) -> ApproxSignedBasis<ValueT> {
        self.basis
    }
}

/// Ggsw Parameters.
pub type GgswParameters<ValueT, ModulusT> = GlevParameters<ValueT, ModulusT>;

/// Big Unsigned Integer Ggsw Parameters.
#[derive(Clone)]
pub struct CrtGlevParameters<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> {
    /// The dimension, refers to **k** in the paper.
    dimension: usize,
    /// The dimension, refers to **N** in the paper.
    poly_length: usize,
    /// cipher modulus minus one, refers to **Q-1**.
    modulus_minus_one: Vec<ValueT>,
    /// The modulus, refers to **Q** in the paper.
    modulus: Vec<ValueT>,
    moduli: Vec<ModulusT>,
    /// The distribution type of the secret key.
    secret_key_type: RingSecretKeyType,
    /// The noise error's standard deviation.
    noise_standard_deviation: f64,
    /// The noise's distribution.
    noise_distribution: SignedDiscreteGaussian<ValueT::SignedInteger>,
    /// Decompose basis for `Q`.
    basis: BigUintApproxSignedBasis<ValueT>,
}

impl<ValueT: UnsignedInteger, ModulusT: FieldContext<ValueT>> CrtGlevParameters<ValueT, ModulusT> {
    /// Creates a new [`CrtGlevParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn new(
        dimension: usize,
        poly_length: usize,
        moduli: &[ModulusT],
        secret_key_type: RingSecretKeyType,
        noise_standard_deviation: f64,
        basis: BigUintApproxSignedBasis<ValueT>,
    ) -> Self {
        let modulus_values: Vec<ValueT> = moduli.iter().map(|m| m.value_unchecked()).collect();
        let modulus = multiply_many_values(&modulus_values);
        let modulus_minus_one = {
            let mut temp = modulus.clone();
            temp[0] -= ValueT::ONE;
            temp
        };

        let noise_distribution =
            SignedDiscreteGaussian::new(0.0, noise_standard_deviation).unwrap();

        Self {
            dimension,
            poly_length,
            modulus_minus_one,
            modulus,
            moduli: moduli.to_vec(),
            secret_key_type,
            noise_standard_deviation,
            basis,
            noise_distribution,
        }
    }

    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    #[inline]
    pub fn modulus(&self) -> &[ValueT] {
        &self.modulus
    }

    #[inline]
    pub fn value_len(&self) -> usize {
        self.modulus.len()
    }

    #[inline]
    pub fn moduli(&self) -> &[ModulusT] {
        &self.moduli
    }

    #[inline]
    pub fn moduli_count(&self) -> usize {
        self.moduli.len()
    }

    #[inline]
    pub fn secret_key_type(&self) -> RingSecretKeyType {
        self.secret_key_type
    }

    /// Returns a reference to the noise distribution of this [`CrtGlevParameters<ValueT, ModulusT>`].
    #[inline]
    pub fn noise_distribution(&self) -> &SignedDiscreteGaussian<ValueT::SignedInteger> {
        &self.noise_distribution
    }

    #[inline]
    pub fn basis(&self) -> &BigUintApproxSignedBasis<ValueT> {
        &self.basis
    }
}

/// Big Unsigned Integer Ggsw Parameters.
pub type CrtGgswParameters<ValueT, ModulusT> = CrtGlevParameters<ValueT, ModulusT>;
