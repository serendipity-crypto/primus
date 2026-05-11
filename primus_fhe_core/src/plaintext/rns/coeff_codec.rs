//! Pipeline:
//!
//!   plaintext m: Polynomial<T> in Z_t
//!     │
//!     │  encode_coeffs / add_encode_coeffs_assign
//!     ▼
//!   m·Δ: CrtPolynomial<T>  (coefficient domain, ready to NTT)
//!     │
//!     │  caller: NTT each modulus chunk
//!     ▼
//!   (encryption operates in DcrtPolynomial)
//!     ⋮
//!   phase output: DcrtPolynomial<T>
//!     │
//!     │  caller: inverse-NTT each modulus chunk
//!     ▼
//!   coefficient-domain DcrtPolynomial<T>
//!     │
//!     │  decode_coeffs (consumes msg_mod_q as workspace + scratch buffers)
//!     ▼
//!   recovered m: Polynomial<T>

use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{
    BigUint, Data, DataMut, DivRemScalar, RawData, UnsignedInteger, multiply_many_values,
};
use primus_modulo::ops::*;
use primus_poly::{CrtPolynomial, DcrtPolynomial, Polynomial};
use primus_reduce::FieldContext;
use primus_rns::{BaseConverter, RNSBase};

pub struct RnsCoeffDecodeContext<T: UnsignedInteger> {
    msg_mod_q: DcrtPolynomial<Vec<T>>,
    fast_convert_buffer: Vec<T>,
}

pub struct RnsCoeffDecodeContextRefMut<'a, T: UnsignedInteger> {
    pub msg_mod_q: &'a mut DcrtPolynomial<Vec<T>>,
    pub fast_convert_buffer: &'a mut [T],
}

impl<T: UnsignedInteger> RnsCoeffDecodeContext<T> {
    /// Creates a new [`RnsCoeffDecodeContext<T>`].
    #[inline]
    pub fn new(moduli_count: usize, poly_length: usize) -> Self {
        let msg_mod_q: DcrtPolynomial<Vec<T>> = DcrtPolynomial::zero(moduli_count * poly_length);
        let fast_convert_buffer = vec![T::ZERO; moduli_count * poly_length];

        Self {
            msg_mod_q,
            fast_convert_buffer,
        }
    }

    #[inline]
    pub fn as_mut(&mut self) -> RnsCoeffDecodeContextRefMut<'_, T> {
        RnsCoeffDecodeContextRefMut {
            msg_mod_q: &mut self.msg_mod_q,
            fast_convert_buffer: &mut self.fast_convert_buffer,
        }
    }
}

/// BFV-style RNS coefficient codec.
///
/// Encodes `m ∈ Z_t` as RNS residues of `m · Δ mod Q` where `Δ = round(Q/t)`,
/// and decodes via the HPS / Bajard et al. fast base extension to `{t, γ}`.
#[derive(Clone)]
pub struct RnsCoeffCodec<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    // Plaintext / ciphertext bases
    t: T,
    base_q: RNSBase<T, M>,
    moduli_values: Vec<T>,

    // Encode 用：Δ = round(Q / t)（BigUint），分解为 Q 各素数下的 Shoup 因子
    delta: BigUint<Vec<T>>,
    delta_factor_mod_q: Vec<ShoupFactor<T>>,

    // Decode 用：HPS γ-trick，输出 m mod t
    gamma: T,
    base_t_gamma: RNSBase<T, M>,               // {t, γ}
    t_gamma_factor_mod_q: Vec<ShoupFactor<T>>, // [(t·γ) mod q_i]
    minus_inv_q_mod_t_gamma: Vec<T>,           // [(−Q^{-1}) mod m_j], m_j ∈ {t, γ}
    inv_gamma_mod_t: ShoupFactor<T>,           // (γ^{-1}) mod t
    converter_q_to_t_gamma: BaseConverter<T, M>,
}

impl<T, M> RnsCoeffCodec<T, M>
where
    T: UnsignedInteger,
    M: FieldContext<T>,
{
    /// Builds a BFV-style RNS codec.
    ///
    /// # Panics
    /// - `t` not coprime with `gamma` or any `q_i`
    /// - `gamma <= t` (HPS requires that γ > t)
    pub fn new(t_modulus: M, base_q: RNSBase<T, M>, gamma_modulus: M) -> Self {
        let t = t_modulus.value_unchecked();
        let gamma = gamma_modulus.value_unchecked();

        let moduli_values: Vec<T> = base_q
            .moduli()
            .iter()
            .map(|m| m.value_unchecked())
            .collect();

        Self::validate_moduli(t, &moduli_values, gamma);

        let cipher_modulus = base_q.moduli_product();

        let mut delta = BigUint(vec![T::ZERO; cipher_modulus.len()]);

        let rem = DivRemScalar::div_rem_scalar(cipher_modulus.digits(), t, delta.digits_mut());
        if rem >= ((t >> 1u32) + (t & T::ONE)) {
            let _ = delta.add_value_assign(T::ONE);
        }

        let delta_factor_mod_q = base_q.decompose_to_rns_factor(delta.view());

        let t_gamma = [t_modulus, gamma_modulus];
        let base_t_gamma = RNSBase::new(&t_gamma).unwrap();
        let q_mod_t_gamma = base_t_gamma.decompose(cipher_modulus.view());
        let minus_inv_q_mod_t_gamma: Vec<T> = q_mod_t_gamma
            .iter()
            .zip(&t_gamma)
            .map(|(&x, modulus)| modulus.reduce_neg(modulus.reduce_inv(x)))
            .collect();
        let inv_gamma_mod_t = ShoupFactor::new(gamma.modulo(t_modulus).inv_modulo(t_modulus), t);
        let t_gamma_value = multiply_many_values(&[t, gamma]);
        let t_gamma_factor_mod_q = base_q.decompose_to_rns_factor(t_gamma_value.view());

        let converter_q_to_t_gamma = BaseConverter::new(&base_q, &base_t_gamma);

        Self {
            t,
            base_q,
            moduli_values,
            delta,
            delta_factor_mod_q,
            gamma,
            base_t_gamma,
            t_gamma_factor_mod_q,
            minus_inv_q_mod_t_gamma,
            inv_gamma_mod_t,
            converter_q_to_t_gamma,
        }
    }

    fn validate_moduli(plain_modulus_value: T, cipher_moduli_value: &[T], gamma: T) {
        assert!(
            plain_modulus_value >= T::TWO,
            "plain modulus must be at least 2"
        );
        assert!(
            gamma > plain_modulus_value,
            "gamma modulus must be greater than the plain modulus for HPS decoding"
        );
        assert!(
            plain_modulus_value.is_coprime(gamma),
            "plain modulus and gamma modulus must be coprime"
        );
        assert!(
            gamma <= T::MAX / T::TWO,
            "gamma too large: HPS final step may overflow"
        );

        for &qi in cipher_moduli_value {
            assert!(
                qi.is_coprime(plain_modulus_value),
                "cipher moduli must be coprime with the plain modulus"
            );
            assert!(
                qi > plain_modulus_value,
                "each RNS ciphertext modulus must be greater than the plain modulus for centered coefficient lifting"
            );
            assert!(
                qi.is_coprime(gamma),
                "cipher moduli must be coprime with the gamma modulus"
            );
        }
    }

    pub fn t(&self) -> T {
        self.t
    }

    pub fn base_q(&self) -> &RNSBase<T, M> {
        &self.base_q
    }

    pub fn moduli_count(&self) -> usize {
        self.base_q.moduli_count()
    }

    pub fn delta(&self) -> BigUint<&[T]> {
        self.delta.view()
    }

    pub fn gamma(&self) -> T {
        self.gamma
    }

    pub fn centered_encode_coeffs<A, B>(
        &self,
        message: &Polynomial<A>,
        crt_message: &mut CrtPolynomial<B>,
        poly_length: usize,
    ) where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        self.base_q.wrapping_decompose_small_polynomial_inplace(
            message,
            crt_message,
            poly_length,
            self.t,
        );

        crt_message.mul_factor_assign(&self.delta_factor_mod_q, poly_length, &self.moduli_values);
    }

    pub fn add_centered_encode_coeffs_assign<A, B>(
        &self,
        message: &Polynomial<A>,
        destination: &mut CrtPolynomial<B>,
        poly_length: usize,
    ) where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
    {
        self.base_q.add_wrapping_decompose_small_polynomial_scaled(
            message,
            destination,
            poly_length,
            self.t,
            &self.delta_factor_mod_q,
        );
    }

    pub fn decode_coeffs<A, B>(
        &self,
        msg_mod_q: &mut DcrtPolynomial<A>,
        msg: &mut Polynomial<B>,
        poly_length: usize,
        fast_convert_buffer: &mut [T],
    ) where
        A: RawData<Elem = T> + DataMut,
        B: RawData<Elem = T> + DataMut,
    {
        let rns_poly_len = self
            .moduli_count()
            .checked_mul(poly_length)
            .expect("RNS polynomial length overflow");
        assert_eq!(msg_mod_q.as_ref().len(), rns_poly_len);
        assert_eq!(msg.as_ref().len(), poly_length);
        assert_eq!(fast_convert_buffer.len(), rns_poly_len);

        let t = self.t;
        let gamma = self.gamma;
        let t_modulus = self.base_t_gamma.moduli()[0];
        let gamma_modulus = self.base_t_gamma.moduli()[1];
        let minus_inv_q_mod_t = self.minus_inv_q_mod_t_gamma[0];
        let minus_inv_q_mod_gamma = self.minus_inv_q_mod_t_gamma[1];
        let inv_gamma_mod_t = self.inv_gamma_mod_t;
        let msg = msg.as_mut();

        // HPS γ-trick decode (Bajard et al. 2017):
        //  1. Multiply by t·γ mod q_i: msg_mod_q := t·γ·c mod q_i
        //  2. Fast base-extend Q -> {t, γ}: msg_mod_t_gamma := round(t·γ·c/Q) mod {t,γ}
        //  3. Multiply by -Q^{-1} mod {t,γ}: yields (y_t, y_γ) = (m·γ + ε, ε) approx
        //  4. Centered-lift y_γ and recover m·γ mod t = y_t - centered(y_γ)
        //  5. Multiply by γ^{-1} mod t to get m

        msg_mod_q.mul_factor_assign(&self.t_gamma_factor_mod_q, poly_length, &self.moduli_values);

        self.converter_q_to_t_gamma
            .fast_convert_array_to_pair_iter(msg_mod_q.as_ref(), poly_length, fast_convert_buffer)
            .zip(msg.iter_mut())
            .for_each(|((y_t, y_gamma), m)| {
                let y_t = t_modulus.reduce_mul(y_t, minus_inv_q_mod_t);
                let y_gamma = gamma_modulus.reduce_mul(y_gamma, minus_inv_q_mod_gamma);

                let mut temp = gamma - y_gamma + y_t;
                if temp >= gamma {
                    temp -= gamma;
                }
                *m = inv_gamma_mod_t.factor_mul_modulo(temp, t);
            });
    }
}
