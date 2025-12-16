use primus_decompose::{big_integer::BigUintApproxSignedBasis, primitive::ApproxSignedBasis};
use primus_distr::{DiscreteGaussian, SignedDiscreteGaussian};
use primus_factor::ShoupFactor;
use primus_integer::{BigUint, DivRemScalar, UnsignedInteger, multiply_many_values};
use primus_modulo::ops::*;
use primus_reduce::FieldContext;
use primus_rns::{BaseConverter, RNSBase};
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
            DiscreteGaussian::new(noise_standard_deviation, cipher_modulus_minus_one).unwrap();

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
        DiscreteGaussian::new(sigma, self.cipher_modulus_minus_one).unwrap()
    }
}

#[derive(Clone, Copy)]
pub struct RNSGlweCommonSize {
    /// The dimension, refers to **k** in the paper.
    dimension: usize,
    /// The polynomial length, refers to **N** in the paper.
    poly_length: usize,
    /// The moduli count, refers to *Q1, Q2, ...* in the paper.
    moduli_count: usize,
    /// The modulus size, refers to *Q* in the paper.
    big_uint_value_len: usize,
    rns_poly_len: usize,
    big_uint_poly_len: usize,
    rns_glwe_mid: usize,
    rns_glwe_len: usize,
}

impl RNSGlweCommonSize {
    pub fn new(
        dimension: usize,
        poly_length: usize,
        moduli_count: usize,
        big_uint_value_len: usize,
    ) -> Self {
        assert!(poly_length.is_power_of_two());

        let big_uint_poly_len = poly_length * big_uint_value_len;
        let rns_poly_len = poly_length * moduli_count;
        let rns_glwe_mid = dimension * rns_poly_len;
        let rns_glwe_len = rns_glwe_mid + rns_poly_len;

        Self {
            dimension,
            poly_length,
            moduli_count,
            big_uint_value_len,
            rns_poly_len,
            big_uint_poly_len,
            rns_glwe_mid,
            rns_glwe_len,
        }
    }

    pub fn dimension(&self) -> usize {
        self.dimension
    }

    pub fn poly_length(&self) -> usize {
        self.poly_length
    }

    pub fn moduli_count(&self) -> usize {
        self.moduli_count
    }

    pub fn big_uint_value_len(&self) -> usize {
        self.big_uint_value_len
    }

    pub fn rns_poly_len(&self) -> usize {
        self.rns_poly_len
    }

    pub fn big_uint_poly_len(&self) -> usize {
        self.big_uint_poly_len
    }

    pub fn rns_glwe_mid(&self) -> usize {
        self.rns_glwe_mid
    }

    pub fn secret_key_len(&self) -> usize {
        self.rns_glwe_mid
    }

    pub fn rns_glwe_len(&self) -> usize {
        self.rns_glwe_len
    }

    pub fn public_key_len(&self) -> usize {
        self.rns_glwe_len
    }
}

/// Big Unsigned Integer Glwe Parameters.
#[derive(Clone)]
pub struct CrtGlweParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    common_size: RNSGlweCommonSize,
    /// The message modulus, refers to **t** in the paper.
    plain_modulus_value: T,
    /// The message modulus, refers to **t** in the paper.
    plain_modulus: M,
    /// The cipher modulus minus one, refers to **Q-1**.
    cipher_modulus_minus_one: Vec<T>,
    /// The moduli, refers to *Q1, Q2, ...* in the paper.
    cipher_moduli: Vec<M>,
    /// The moduli, refers to *Q1, Q2, ...* in the paper.
    cipher_moduli_value: Vec<T>,
    /// Refers to `Q1-1`, `Q2-1` ...
    cipher_moduli_minus_one: Vec<T>,
    /// The uniform distribuition to sample values over `Q1`, `Q2` ...
    cipher_moduli_uniform_distr: Vec<Uniform<T>>,
    /// Residue Number System for *Q*.
    base_q: RNSBase<T, M>,
    /// Refers to `Q/t`.
    delta: Vec<T>,
    delta_mod_q: Vec<T>,
    inv_delta_mod_q: Vec<T>,
    gamma: T,
    base_t_gamma: RNSBase<T, M>,
    t_gamma_mod_q: Vec<T>,
    minus_inv_q_mod_t_gamma: Vec<T>,
    inv_gamma_mod_t: ShoupFactor<T>,
    converter: BaseConverter<T, M>,
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
        plain_modulus: M,
        gamma_modulus: M,
        cipher_moduli: &[M],
        secret_key_type: RingSecretKeyType,
        noise_standard_deviation: f64,
    ) -> Self {
        let t = plain_modulus.value_unchecked();
        let gamma = gamma_modulus.value_unchecked();

        let cipher_moduli_value: Vec<T> = cipher_moduli
            .iter()
            .map(|qi| qi.value_unchecked())
            .collect();
        let cipher_moduli_minus_one = cipher_moduli_value.iter().map(|&qi| qi - T::ONE).collect();
        let base_q = RNSBase::new(cipher_moduli).unwrap();
        let cipher_modulus = base_q.moduli_product();
        let cipher_modulus_minus_one = {
            let mut temp = BigUint(cipher_modulus.0.to_vec());
            let _ = temp.sub_value_assign(T::ONE);
            temp
        };

        let cipher_moduli_uniform_distr = cipher_moduli
            .iter()
            .map(|qi| qi.uniform_distribution())
            .collect();

        let noise_distribution = SignedDiscreteGaussian::new(noise_standard_deviation).unwrap();

        let mut delta = BigUint(vec![T::ZERO; cipher_modulus.len()]);

        let rem = DivRemScalar::div_rem_scalar(cipher_modulus.0, t, delta.digits_mut());
        if rem * T::TWO >= t {
            let _ = delta.add_value_assign(T::ONE);
        }

        let delta_mod_q: Vec<T> = base_q.decompose(delta.digits());

        let inv_delta_mod_q: Vec<T> = delta_mod_q
            .iter()
            .zip(cipher_moduli)
            .map(|(&v, modulus)| modulus.reduce_inv(v))
            .collect();

        let t_gamma = [plain_modulus, gamma_modulus];
        let base_t_gamma = RNSBase::new(&t_gamma).unwrap();
        let q_mod_t_gamma = base_t_gamma.decompose(cipher_modulus.digits());
        let minus_inv_q_mod_t_gamma: Vec<T> = q_mod_t_gamma
            .iter()
            .zip(&t_gamma)
            .map(|(&x, modulus)| modulus.reduce_neg(modulus.reduce_inv(x)))
            .collect();
        let inv_gamma_mod_t =
            ShoupFactor::new(gamma.modulo(plain_modulus).inv_modulo(plain_modulus), t);
        let t_gamma_value = multiply_many_values(&[t, gamma]);
        let t_gamma_mod_q: Vec<T> = base_q.decompose(t_gamma_value.as_ref());

        let converter = BaseConverter::new(&base_q, &base_t_gamma);

        let common_size = RNSGlweCommonSize::new(
            dimension,
            poly_length,
            base_q.moduli_count(),
            base_q.big_uint_value_len(),
        );

        Self {
            common_size,
            plain_modulus_value: t,
            plain_modulus,
            cipher_modulus_minus_one: cipher_modulus_minus_one.0,
            cipher_moduli: cipher_moduli.to_vec(),
            cipher_moduli_value,
            cipher_moduli_minus_one,
            cipher_moduli_uniform_distr,
            delta: delta.0,
            delta_mod_q,
            inv_delta_mod_q,
            gamma,
            t_gamma_mod_q,
            minus_inv_q_mod_t_gamma,
            inv_gamma_mod_t,
            base_q,
            base_t_gamma,
            converter,
            secret_key_type,
            noise_distribution,
        }
    }

    /// Returns the dimension of this [`CrtGlweParameters<T, M>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.common_size.dimension()
    }

    /// Returns the poly length of this [`CrtGlweParameters<T, M>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.common_size.poly_length()
    }

    /// Returns the plain modulus value of this [`CrtGlweParameters<T, M>`].
    pub fn plain_modulus_value(&self) -> T {
        self.plain_modulus_value
    }

    pub fn plain_modulus(&self) -> M {
        self.plain_modulus
    }

    /// Returns a reference to the cipher modulus of this [`CrtGlweParameters<T, M>`].
    pub fn cipher_modulus(&self) -> BigUint<&[T]> {
        self.base_q.moduli_product()
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
        self.base_q.big_uint_value_len()
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
    pub fn delta_mod_q(&self) -> &[T] {
        &self.delta_mod_q
    }

    /// Returns a reference to the inverse delta residues of this [`CrtGlweParameters<T, M>`].
    pub fn inv_delta_mod_q(&self) -> &[T] {
        &self.inv_delta_mod_q
    }

    pub fn t_gamma_mod_q(&self) -> &[T] {
        &self.t_gamma_mod_q
    }

    pub fn converter(&self) -> &BaseConverter<T, M> {
        &self.converter
    }

    pub fn minus_inv_q_mod_t_gamma(&self) -> &[T] {
        &self.minus_inv_q_mod_t_gamma
    }

    pub fn t_gamma(&self) -> &[M] {
        self.base_t_gamma.moduli()
    }

    pub fn gamma(&self) -> T {
        self.gamma
    }

    pub fn inv_gamma_mod_t(&self) -> ShoupFactor<T> {
        self.inv_gamma_mod_t
    }

    pub fn base_q(&self) -> &RNSBase<T, M> {
        &self.base_q
    }

    pub fn common_size(&self) -> RNSGlweCommonSize {
        self.common_size
    }

    pub fn big_uint_poly_len(&self) -> usize {
        self.common_size.big_uint_poly_len()
    }

    pub fn rns_poly_len(&self) -> usize {
        self.common_size.rns_poly_len()
    }

    pub fn rns_glwe_mid(&self) -> usize {
        self.common_size.rns_glwe_mid()
    }

    pub fn rns_glwe_len(&self) -> usize {
        self.common_size.rns_glwe_len()
    }

    pub fn secret_key_len(&self) -> usize {
        self.common_size.secret_key_len()
    }

    pub fn public_key_len(&self) -> usize {
        self.common_size.public_key_len()
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
            DiscreteGaussian::new(noise_standard_deviation, cipher_modulus_minus_one).unwrap();

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

#[derive(Clone, Copy)]
pub struct RNSGlevCommonSize {
    rns_glwe_common_size: RNSGlweCommonSize,
    decompose_length: usize,
    rns_glev_len: usize,
    rns_ggsw_len: usize,
}

impl RNSGlevCommonSize {
    pub fn new(rns_glwe_common_size: RNSGlweCommonSize, decompose_length: usize) -> Self {
        let rns_glev_len = decompose_length * rns_glwe_common_size.rns_glwe_len();
        let rns_ggsw_len = rns_glev_len * (rns_glwe_common_size.dimension() + 1);
        Self {
            rns_glwe_common_size,
            decompose_length,
            rns_glev_len,
            rns_ggsw_len,
        }
    }

    pub fn decompose_length(&self) -> usize {
        self.decompose_length
    }

    pub fn rns_glev_len(&self) -> usize {
        self.rns_glev_len
    }

    pub fn rns_ggsw_len(&self) -> usize {
        self.rns_ggsw_len
    }

    pub fn dimension(&self) -> usize {
        self.rns_glwe_common_size.dimension()
    }

    pub fn moduli_count(&self) -> usize {
        self.rns_glwe_common_size.moduli_count()
    }

    pub fn poly_length(&self) -> usize {
        self.rns_glwe_common_size.poly_length()
    }

    pub fn rns_poly_len(&self) -> usize {
        self.rns_glwe_common_size.rns_poly_len()
    }

    pub fn rns_glwe_mid(&self) -> usize {
        self.rns_glwe_common_size.rns_glwe_mid()
    }

    pub fn rns_glwe_len(&self) -> usize {
        self.rns_glwe_common_size.rns_glwe_len()
    }

    pub fn secret_key_len(&self) -> usize {
        self.rns_glwe_common_size.secret_key_len()
    }

    pub fn public_key_len(&self) -> usize {
        self.rns_glwe_common_size.public_key_len()
    }

    pub fn big_uint_poly_len(&self) -> usize {
        self.rns_glwe_common_size.big_uint_poly_len()
    }
}

/// Big Unsigned Integer Ggsw Parameters.
#[derive(Clone)]
pub struct CrtGlevParameters<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    common_size: RNSGlevCommonSize,
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
        let cipher_moduli_value: Vec<T> = cipher_moduli
            .iter()
            .map(|qi| qi.value_unchecked())
            .collect();
        let cipher_moduli_minus_one = cipher_moduli_value.iter().map(|&qi| qi - T::ONE).collect();
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

        let noise_distribution = SignedDiscreteGaussian::new(noise_standard_deviation).unwrap();

        let rns_glwe_common_size = RNSGlweCommonSize::new(
            dimension,
            poly_length,
            cipher_moduli.len(),
            cipher_modulus.len(),
        );

        let common_size = RNSGlevCommonSize::new(rns_glwe_common_size, basis.decompose_length());

        Self {
            common_size,
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

    pub fn with_glwe_params(
        glwe_params: &CrtGlweParameters<T, M>,
        basis: BigUintApproxSignedBasis<T>,
    ) -> Self {
        let decompose_length = basis.decompose_length();
        Self {
            cipher_modulus_minus_one: glwe_params.cipher_modulus_minus_one().to_vec(),
            cipher_modulus: glwe_params.cipher_modulus().0.to_vec(),
            cipher_moduli: glwe_params.cipher_moduli().to_vec(),
            cipher_moduli_value: glwe_params.cipher_moduli_value().to_vec(),
            cipher_moduli_minus_one: glwe_params.cipher_moduli_minus_one().to_vec(),
            cipher_moduli_uniform_distr: glwe_params.cipher_moduli_uniform_distr().to_vec(),
            secret_key_type: glwe_params.secret_key_type,
            noise_distribution: glwe_params.noise_distribution().clone(),
            basis,
            common_size: RNSGlevCommonSize::new(glwe_params.common_size(), decompose_length),
        }
    }

    /// Returns the dimension of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn dimension(&self) -> usize {
        self.common_size.dimension()
    }

    /// Returns the poly length of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn poly_length(&self) -> usize {
        self.common_size.poly_length()
    }

    /// Returns a reference to the cipher modulus of this [`CrtGlevParameters<T, M>`].
    #[inline]
    pub fn cipher_modulus(&self) -> BigUint<&[T]> {
        BigUint(&self.cipher_modulus)
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

    pub fn common_size(&self) -> RNSGlevCommonSize {
        self.common_size
    }

    pub fn rns_glev_len(&self) -> usize {
        self.common_size.rns_glev_len()
    }

    pub fn rns_ggsw_len(&self) -> usize {
        self.common_size.rns_ggsw_len()
    }

    pub fn rns_poly_len(&self) -> usize {
        self.common_size.rns_poly_len()
    }

    pub fn rns_glwe_mid(&self) -> usize {
        self.common_size.rns_glwe_mid()
    }

    pub fn rns_glwe_len(&self) -> usize {
        self.common_size.rns_glwe_len()
    }

    pub fn decompose_length(&self) -> usize {
        self.basis.decompose_length()
    }

    pub fn big_uint_poly_len(&self) -> usize {
        self.common_size.big_uint_poly_len()
    }
}

/// Big Unsigned Integer Ggsw Parameters.
pub type CrtGgswParameters<T, M> = CrtGlevParameters<T, M>;
