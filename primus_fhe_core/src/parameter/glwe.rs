use primus_decompose::{big_integer::BigUintApproxSignedBasis, primitive::ApproxSignedBasis};
use primus_distr::{DiscreteGaussian, SignedDiscreteGaussian};
use primus_integer::{BigIntegerOps, DivRemScalar, UnsignedInteger, multiply_many_values};
use primus_reduce::FieldContext;
use rand::distr::Uniform;

use crate::RingSecretKeyType;

/// Glwe Parameters.
#[derive(Clone)]
pub struct GlweParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// The dimension, refers to **k** in the paper.
    dimension: usize,
    /// The polynomial length, refers to **N** in the paper.
    poly_length: usize,
    /// **RLWE** message modulus, refers to **t** in the paper.
    plain_modulus_value: T,
    /// **RLWE** cipher modulus minus one, refers to **Q-1**.
    cipher_modulus_minus_one: T,
    /// The modulus, refers to **Q** in the paper.
    cipher_modulus: M,
    /// The distribution type of the secret key.
    secret_key_type: RingSecretKeyType,
    /// The noise's distribution.
    noise_distribution: DiscreteGaussian<T>,
}

impl<T, M> GlweParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// Creates a new [`GlweParameters<T, M>`].
    pub fn new(
        dimension: usize,
        poly_length: usize,
        plain_modulus_value: T,
        cipher_modulus: M,
        secret_key_type: RingSecretKeyType,
        noise_standard_deviation: f64,
    ) -> Self {
        let cipher_modulus_minus_one = cipher_modulus.minus_one();

        let noise_distribution =
            DiscreteGaussian::new(0.0, noise_standard_deviation, cipher_modulus_minus_one).unwrap();

        Self {
            dimension,
            poly_length,
            plain_modulus_value,
            cipher_modulus_minus_one,
            cipher_modulus,
            secret_key_type,
            noise_distribution,
        }
    }

    /// Returns the dimension of this [`GlweParameters<T, M>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the poly length of this [`GlweParameters<T, M>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the plain modulus value of this [`GlweParameters<T, M>`].
    pub fn plain_modulus_value(&self) -> T {
        self.plain_modulus_value
    }

    /// Returns the cipher modulus of this [`GlweParameters<T, M>`].
    pub fn cipher_modulus(&self) -> M {
        self.cipher_modulus
    }

    /// Returns the cipher modulus of this [`GlweParameters<T, M>`].
    #[inline]
    pub fn cipher_modulus_value(&self) -> T {
        self.cipher_modulus.value_unchecked()
    }

    /// Returns the cipher modulus minus one of this [`GlweParameters<T, M>`].
    pub fn cipher_modulus_minus_one(&self) -> T {
        self.cipher_modulus_minus_one
    }

    /// Returns the secret key type of this [`GlweParameters<T, M>`].
    pub fn secret_key_type(&self) -> RingSecretKeyType {
        self.secret_key_type
    }

    /// Returns a reference to the noise distribution of this [`GlweParameters<T, M>`].
    #[inline]
    pub fn noise_distribution(&self) -> &DiscreteGaussian<T> {
        &self.noise_distribution
    }

    /// Returns the noise distribution.
    #[inline]
    pub fn noise_distribution_div_count(&self, count: u32) -> DiscreteGaussian<T> {
        let noise_standard_deviation = self.noise_distribution.standard_deviation();
        let var = noise_standard_deviation * noise_standard_deviation;
        let sigma = (var / count as f64).sqrt();
        DiscreteGaussian::new(0.0, sigma, self.cipher_modulus_minus_one).unwrap()
    }
}

/// Big Unsigned Integer Glwe Parameters.
#[derive(Clone)]
pub struct CrtGlweParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// The dimension, refers to **k** in the paper.
    dimension: usize,
    /// The polynomial length, refers to **N** in the paper.
    poly_length: usize,
    /// **RLWE** message modulus, refers to **t** in the paper.
    plain_modulus_value: T,
    /// **RLWE** cipher modulus minus one, refers to **Q-1**.
    cipher_modulus_minus_one: Vec<T>,
    /// **RLWE** cipher modulus, refers to **Q**.
    cipher_modulus: Vec<T>,
    /// The moduli, refers to **Q=Q1*Q2*...** in the paper.
    cipher_moduli: Vec<M>,
    /// The moduli, refers to **Q=Q1*Q2*...** in the paper.
    cipher_moduli_value: Vec<T>,
    /// Refers to `Q1-1`, `Q2-1` ...
    cipher_moduli_minus_one: Vec<T>,

    cipher_moduli_uniform_distr: Vec<Uniform<T>>,

    /// Refers to `Q/t`.
    delta: Vec<T>,
    delta_residues: Vec<T>,
    inverse_delta_residues: Vec<T>,
    /// The distribution type of the secret key.
    secret_key_type: RingSecretKeyType,
    /// The noise distribution
    noise_distribution: SignedDiscreteGaussian<T::SignedInteger>,
}

impl<T, M> CrtGlweParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// Creates a new [`CrtGlweParameters<T, M>`].
    pub fn new(
        dimension: usize,
        poly_length: usize,
        plain_modulus_value: T,
        cipher_moduli: &[M],
        secret_key_type: RingSecretKeyType,
        noise_standard_deviation: f64,
    ) -> Self {
        let cipher_moduli_value: Vec<T> =
            cipher_moduli.iter().map(|m| m.value_unchecked()).collect();
        let cipher_moduli_minus_one = cipher_moduli_value.iter().map(|&m| m - T::ONE).collect();
        let cipher_modulus = multiply_many_values(&cipher_moduli_value);
        let cipher_modulus_minus_one = {
            let mut temp = cipher_modulus.clone();
            temp[0] -= T::ONE;
            temp
        };

        let cipher_moduli_uniform_distr = cipher_moduli
            .iter()
            .map(|m| m.uniform_distribution())
            .collect();

        let noise_distribution =
            SignedDiscreteGaussian::new(0.0, noise_standard_deviation).unwrap();

        let mut delta = vec![T::ZERO; cipher_modulus.len()];

        let rem = DivRemScalar::div_rem_scalar(&cipher_modulus, plain_modulus_value, &mut delta);

        if rem * T::TWO >= plain_modulus_value {
            let _ = delta.slice_add_value_assign(T::ONE);
        }

        let delta_residues: Vec<T> = cipher_moduli
            .iter()
            .map(|modulus| modulus.reduce(delta.as_ref()))
            .collect();

        let inverse_delta_residues: Vec<T> = delta_residues
            .iter()
            .zip(cipher_moduli)
            .map(|(&v, modulus)| modulus.reduce_inv(v))
            .collect();

        Self {
            dimension,
            poly_length,
            plain_modulus_value,
            cipher_modulus,
            cipher_modulus_minus_one,
            cipher_moduli: cipher_moduli.to_vec(),
            cipher_moduli_value,
            cipher_moduli_minus_one,
            cipher_moduli_uniform_distr,
            delta,
            delta_residues,
            inverse_delta_residues,
            secret_key_type,
            noise_distribution,
        }
    }

    /// Returns the dimension of this [`CrtGlweParameters<T, M>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the poly length of this [`CrtGlweParameters<T, M>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the plain modulus value of this [`CrtGlweParameters<T, M>`].
    pub fn plain_modulus_value(&self) -> T {
        self.plain_modulus_value
    }

    /// Returns a reference to the cipher modulus of this [`CrtGlweParameters<T, M>`].
    pub fn cipher_modulus(&self) -> &[T] {
        &self.cipher_modulus
    }

    /// Returns a reference to the modulus minus one of this [`CrtGlweParameters<T, M>`].
    pub fn cipher_modulus_minus_one(&self) -> &[T] {
        &self.cipher_modulus_minus_one
    }

    /// Returns a reference to the moduli of this [`CrtGlweParameters<T, M>`].
    #[inline]
    pub fn cipher_moduli(&self) -> &[M] {
        &self.cipher_moduli
    }

    /// Returns a reference to the cipher moduli value of this [`CrtGlweParameters<T, M>`].
    pub fn cipher_moduli_value(&self) -> &[T] {
        &self.cipher_moduli_value
    }

    /// Returns a reference to the cipher moduli minus one of this [`CrtGlweParameters<T, M>`].
    pub fn cipher_moduli_minus_one(&self) -> &[T] {
        &self.cipher_moduli_minus_one
    }

    /// Returns the moduli count of this [`CrtGlweParameters<T, M>`].
    pub fn cipher_moduli_count(&self) -> usize {
        self.cipher_moduli.len()
    }

    /// Returns a reference to the cipher moduli uniform distr of this [`CrtGlweParameters<T, M>`].
    pub fn cipher_moduli_uniform_distr(&self) -> &[Uniform<T>] {
        &self.cipher_moduli_uniform_distr
    }

    /// Returns the big uint value len of this [`CrtGlweParameters<T, M>`].
    #[inline]
    pub fn big_uint_value_len(&self) -> usize {
        self.cipher_modulus.len()
    }

    /// Returns the secret key type of this [`CrtGlweParameters<T, M>`].
    pub fn secret_key_type(&self) -> RingSecretKeyType {
        self.secret_key_type
    }

    /// Returns a reference to the noise distribution of this [`CrtGlweParameters<T, M>`].
    pub fn noise_distribution(&self) -> &SignedDiscreteGaussian<T::SignedInteger> {
        &self.noise_distribution
    }

    /// Returns the noise standard deviation of this [`CrtGlweParameters<T, M>`].
    pub fn noise_standard_deviation(&self) -> f64 {
        self.noise_distribution.standard_deviation()
    }

    /// Returns a reference to the delta of this [`CrtGlweParameters<T, M>`].
    pub fn delta(&self) -> &[T] {
        &self.delta
    }

    /// Returns a reference to the delta residues of this [`CrtGlweParameters<T, M>`].
    pub fn delta_residues(&self) -> &[T] {
        &self.delta_residues
    }

    /// Returns a reference to the inverse delta residues of this [`CrtGlweParameters<T, M>`].
    pub fn inverse_delta_residues(&self) -> &[T] {
        &self.inverse_delta_residues
    }
}

/// Glev Parameters.
#[derive(Clone)]
pub struct GlevParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// The dimension, refers to **k** in the paper.
    dimension: usize,
    /// The dimension, refers to **N** in the paper.
    poly_length: usize,
    /// cipher modulus minus one, refers to **Q-1**.
    cipher_modulus_minus_one: T,
    /// The modulus, refers to **Q** in the paper.
    cipher_modulus: M,
    /// The distribution type of the secret key.
    secret_key_type: RingSecretKeyType,
    /// The noise's distribution.
    noise_distribution: DiscreteGaussian<T>,
    /// Decompose basis for `Q`.
    basis: ApproxSignedBasis<T>,
}

impl<T, M> GlevParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// Creates a new [`GlevParameters<T, M>`].
    #[inline]
    pub fn new(
        dimension: usize,
        poly_length: usize,
        cipher_modulus: M,
        secret_key_type: RingSecretKeyType,
        noise_standard_deviation: f64,
        basis: ApproxSignedBasis<T>,
    ) -> Self {
        let cipher_modulus_minus_one = cipher_modulus.minus_one();

        let noise_distribution =
            DiscreteGaussian::new(0.0, noise_standard_deviation, cipher_modulus_minus_one).unwrap();

        Self {
            dimension,
            poly_length,
            cipher_modulus_minus_one,
            cipher_modulus,
            secret_key_type,
            basis,
            noise_distribution,
        }
    }

    /// Returns the dimension of this [`GlevParameters<T, M>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the poly length of this [`GlevParameters<T, M>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns the cipher modulus minus one of this [`GlevParameters<T, M>`].
    pub fn cipher_modulus_minus_one(&self) -> T {
        self.cipher_modulus_minus_one
    }

    /// Returns the modulus of this [`GlevParameters<T, M>`].
    #[inline]
    pub fn cipher_modulus(&self) -> M {
        self.cipher_modulus
    }

    /// Returns the secret key type of this [`GlevParameters<T, M>`].
    pub fn secret_key_type(&self) -> RingSecretKeyType {
        self.secret_key_type
    }

    /// Returns a reference to the noise distribution of this [`GlevParameters<T, M>`].
    #[inline]
    pub fn noise_distribution(&self) -> &DiscreteGaussian<T> {
        &self.noise_distribution
    }

    /// Returns the noise standard deviation of this [`GlevParameters<T, M>`].
    pub fn noise_standard_deviation(&self) -> f64 {
        self.noise_distribution.standard_deviation()
    }

    /// Returns a reference to the basis of this [`GlevParameters<T, M>`].
    #[inline]
    pub fn basis(&self) -> &ApproxSignedBasis<T> {
        &self.basis
    }
}

/// Ggsw Parameters.
pub type GgswParameters<T, M> = GlevParameters<T, M>;

/// Big Unsigned Integer Ggsw Parameters.
#[derive(Clone)]
pub struct CrtGlevParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// The dimension, refers to **k** in the paper.
    dimension: usize,
    /// The dimension, refers to **N** in the paper.
    poly_length: usize,
    /// cipher modulus minus one, refers to **Q-1**.
    cipher_modulus_minus_one: Vec<T>,
    /// The modulus, refers to **Q** in the paper.
    cipher_modulus: Vec<T>,
    /// The moduli, refers to **Q=Q1*Q2*...** in the paper.
    cipher_moduli: Vec<M>,
    /// The moduli, refers to **Q=Q1*Q2*...** in the paper.
    cipher_moduli_value: Vec<T>,
    /// Refers to `Q1-1`, `Q2-1` ...
    cipher_moduli_minus_one: Vec<T>,

    cipher_moduli_uniform_distr: Vec<Uniform<T>>,
    /// The distribution type of the secret key.
    secret_key_type: RingSecretKeyType,
    /// The noise's distribution.
    noise_distribution: SignedDiscreteGaussian<<T as UnsignedInteger>::SignedInteger>,
    /// Decompose basis for `Q`.
    basis: BigUintApproxSignedBasis<T>,
}

impl<T, M> CrtGlevParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// Creates a new [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn new(
        dimension: usize,
        poly_length: usize,
        cipher_moduli: &[M],
        secret_key_type: RingSecretKeyType,
        noise_standard_deviation: f64,
        basis: BigUintApproxSignedBasis<T>,
    ) -> Self {
        let cipher_moduli_value: Vec<T> =
            cipher_moduli.iter().map(|m| m.value_unchecked()).collect();
        let cipher_moduli_minus_one = cipher_moduli_value.iter().map(|&m| m - T::ONE).collect();
        let cipher_modulus = multiply_many_values(&cipher_moduli_value);
        let cipher_modulus_minus_one = {
            let mut temp = cipher_modulus.clone();
            temp[0] -= T::ONE;
            temp
        };

        let cipher_moduli_uniform_distr = cipher_moduli
            .iter()
            .map(|m| m.uniform_distribution())
            .collect();

        let noise_distribution =
            SignedDiscreteGaussian::new(0.0, noise_standard_deviation).unwrap();

        Self {
            dimension,
            poly_length,
            cipher_modulus,
            cipher_modulus_minus_one,
            cipher_moduli: cipher_moduli.to_vec(),
            cipher_moduli_value,
            cipher_moduli_minus_one,
            cipher_moduli_uniform_distr,
            secret_key_type,
            noise_distribution,
            basis,
        }
    }

    /// Returns the dimension of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Returns the poly length of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    /// Returns a reference to the cipher modulus of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn cipher_modulus(&self) -> &[T] {
        &self.cipher_modulus
    }

    /// Returns a reference to the cipher modulus minus one of this [`CrtGlevParameters<T, M>`].
    pub fn cipher_modulus_minus_one(&self) -> &[T] {
        &self.cipher_modulus_minus_one
    }

    /// Returns the big uint value len of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn big_uint_value_len(&self) -> usize {
        self.cipher_modulus.len()
    }

    /// Returns a reference to the moduli of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn cipher_moduli(&self) -> &[M] {
        &self.cipher_moduli
    }

    /// Returns the moduli count of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn cipher_moduli_count(&self) -> usize {
        self.cipher_moduli.len()
    }

    /// Returns a reference to the cipher moduli value of this [`CrtGlevParameters<T, M>`].
    pub fn cipher_moduli_value(&self) -> &[T] {
        &self.cipher_moduli_value
    }

    /// Returns a reference to the cipher moduli minus one of this [`CrtGlevParameters<T, M>`].
    pub fn cipher_moduli_minus_one(&self) -> &[T] {
        &self.cipher_moduli_minus_one
    }

    /// Returns a reference to the cipher moduli uniform distr of this [`CrtGlevParameters<T, M>`].
    pub fn cipher_moduli_uniform_distr(&self) -> &[Uniform<T>] {
        &self.cipher_moduli_uniform_distr
    }

    /// Returns the secret key type of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn secret_key_type(&self) -> RingSecretKeyType {
        self.secret_key_type
    }

    /// Returns a reference to the noise distribution of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn noise_distribution(&self) -> &SignedDiscreteGaussian<T::SignedInteger> {
        &self.noise_distribution
    }

    /// Returns the noise standard deviation of this  [`CrtGlevParameters<T, M>`].
    pub fn noise_standard_deviation(&self) -> f64 {
        self.noise_distribution.standard_deviation()
    }

    /// Returns a reference to the basis of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn basis(&self) -> &BigUintApproxSignedBasis<T> {
        &self.basis
    }
}

/// Big Unsigned Integer Ggsw Parameters.
pub type CrtGgswParameters<T, M> = CrtGlevParameters<T, M>;
