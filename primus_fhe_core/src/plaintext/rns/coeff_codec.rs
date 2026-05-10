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
//!

use primus_factor::{FactorMul, ShoupFactor};
use primus_integer::{
    BigUint, Data, DataMut, DivRemScalar, RawData, UnsignedInteger, izip, multiply_many_values,
};
use primus_modulo::ops::*;
use primus_poly::{CrtPolynomial, DcrtPolynomial, Polynomial};
use primus_reduce::FieldContext;
use primus_rns::{BaseConverter, RNSBase};

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
    moduli_value: Vec<T>,

    // Encode 用：Δ = round(Q / t)（BigUint），分解为 Q 各素数下的 Shoup 因子
    delta: BigUint<Vec<T>>,
    delta_mod_q: Vec<T>,
    delta_factor_mod_q: Vec<ShoupFactor<T>>,

    // Decode 用：HPS γ-trick，输出 m mod t
    gamma: T,
    base_t_gamma: RNSBase<T, M>,     // {t, γ}
    t_gamma_mod_q: Vec<T>,           // [(t·γ) mod q_i]
    minus_inv_q_mod_t_gamma: Vec<T>, // [(−Q^{-1}) mod m_j], m_j ∈ {t, γ}
    inv_gamma_mod_t: ShoupFactor<T>, // (γ^{-1}) mod t
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
    /// - `gamma <= t` (HPS 要求 γ > t)
    pub fn new(t_modulus: M, base_q: RNSBase<T, M>, gamma_modulus: M) -> Self {
        let t = t_modulus.value_unchecked();
        let gamma = gamma_modulus.value_unchecked();

        assert!(gamma > t, "gamma must exceed t for HPS decoding");
        debug_assert!(
            gamma <= T::MAX / T::TWO,
            "gamma too large: HPS final step may overflow"
        );

        let cipher_modulus = base_q.moduli_product();
        let moduli_value: Vec<T> = base_q
            .moduli()
            .iter()
            .map(|m| m.value_unchecked())
            .collect();

        let mut delta = BigUint(vec![T::ZERO; cipher_modulus.len()]);

        let rem = DivRemScalar::div_rem_scalar(cipher_modulus.digits(), t, delta.digits_mut());
        if rem * T::TWO >= t {
            let _ = delta.add_value_assign(T::ONE);
        }

        let delta_mod_q: Vec<T> = base_q.decompose(delta.view());

        let delta_factor_mod_q: Vec<ShoupFactor<T>> = delta_mod_q
            .iter()
            .zip(moduli_value.iter())
            .map(|(&value, &modulus)| ShoupFactor::new(value, modulus))
            .collect();

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
        let t_gamma_mod_q: Vec<T> = base_q.decompose(t_gamma_value.view());

        let converter_q_to_t_gamma = BaseConverter::new(&base_q, &base_t_gamma);

        Self {
            t,
            base_q,
            moduli_value,
            delta,
            delta_mod_q,
            delta_factor_mod_q,
            gamma,
            base_t_gamma,
            t_gamma_mod_q,
            minus_inv_q_mod_t_gamma,
            inv_gamma_mod_t,
            converter_q_to_t_gamma,
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

    pub fn delta_mod_q(&self) -> &[T] {
        &self.delta_mod_q
    }

    pub fn gamma(&self) -> T {
        self.gamma
    }

    pub fn encode_coeffs<A, B>(
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

        crt_message.mul_factor_assign(&self.delta_factor_mod_q, poly_length, &self.moduli_value);
    }

    pub fn add_encode_coeffs_assign<A, B, C>(
        &self,
        message: &Polynomial<A>,
        crt_message: &mut CrtPolynomial<B>,
        destination: &mut CrtPolynomial<C>,
        poly_length: usize,
    ) where
        A: RawData<Elem = T> + Data,
        B: RawData<Elem = T> + DataMut,
        C: RawData<Elem = T> + DataMut,
    {
        self.base_q.wrapping_decompose_small_polynomial_inplace(
            message,
            crt_message,
            poly_length,
            self.t,
        );

        destination.add_mul_factor_assign(
            crt_message,
            &self.delta_factor_mod_q,
            poly_length,
            &self.moduli_value,
        );
    }

    pub fn decode_coeffs<A, B, C>(
        &self,
        msg_mod_q: &mut DcrtPolynomial<A>,
        msg: &mut Polynomial<B>,
        poly_length: usize,
        msg_mod_t_gamma: &mut CrtPolynomial<C>,
        fast_convert_buffer: &mut [T],
    ) where
        A: RawData<Elem = T> + DataMut,
        B: RawData<Elem = T> + DataMut,
        C: RawData<Elem = T> + DataMut,
    {
        let gamma = self.gamma;

        // HPS γ-trick decode (Bajard et al. 2017):
        //  1. Multiply by t·γ mod q_i: msg_mod_q := t·γ·c mod q_i
        //  2. Fast base-extend Q -> {t, γ}: msg_mod_t_gamma := round(t·γ·c/Q) mod {t,γ}
        //  3. Multiply by -Q^{-1} mod {t,γ}: yields (y_t, y_γ) = (m·γ + ε, ε) approx
        //  4. Centered-lift y_γ and recover m·γ mod t = y_t - centered(y_γ)
        //  5. Multiply by γ^{-1} mod t to get m

        msg_mod_q.mul_scalar_assign(&self.t_gamma_mod_q, poly_length, self.base_q.moduli());

        self.converter_q_to_t_gamma.fast_convert_array(
            msg_mod_q.as_ref(),
            msg_mod_t_gamma.as_mut(),
            poly_length,
            fast_convert_buffer,
        );

        msg_mod_t_gamma.mul_scalar_assign(
            &self.minus_inv_q_mod_t_gamma,
            poly_length,
            self.base_t_gamma.moduli(),
        );

        let (y_t_slices, y_gamma_slices) = msg_mod_t_gamma.as_ref().split_at(poly_length);

        izip!(msg.iter_mut(), y_t_slices, y_gamma_slices).for_each(|(res, &y_t, &y_gamma)| {
            let mut temp = gamma - y_gamma + y_t;
            if temp >= gamma {
                temp -= gamma;
            }
            *res = self.inv_gamma_mod_t.factor_mul_modulo(temp, self.t);
        });
    }
}
