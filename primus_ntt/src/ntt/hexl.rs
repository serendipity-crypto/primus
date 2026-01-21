#![cfg(target_arch = "x86_64")]

use aligned_vec::{AVec, avec};
use primus_factor::{FactorMul, MultiplyFactor, ShoupFactor};
use primus_reduce::FieldContext;

use crate::{
    NttError, NttTable, PrimitiveRoot,
    ntt::hexl::{
        fwd::forward_transform_to_bit_reverse_avx512,
        internal::{
            DEFAULT_SHIFT_BITS, IFMA_SHIFT_BITS, MAX_FWD_32_MODULUS, MAX_FWD_IFMA_MODULUS,
            MAX_INV_32_MODULUS, MAX_INV_IFMA_MODULUS, check_arguments,
        },
        inv::inverse_transform_from_bit_reverse_avx512,
        radix2::forward_transform_to_bit_reverse_radix2_inplace,
    },
    reverse::ReverseLsbs,
};

mod fwd;
mod internal;
mod inv;
mod number_theory;
mod radix2;
mod utils;

type Factor = ShoupFactor<u64>;

/// Performs negacyclic forward and inverse number-theoretic transform
/// (NTT), commonly used in RLWE cryptography.
///
/// The number-theoretic transform (NTT) specializes the discrete
/// Fourier transform (DFT) to the finite field \f$ \mathbb{Z}_q[X] / (X^N + 1)
/// \f$.
pub struct HexlNttTable {
    /// size of NTT transform, should be power of 2
    n: usize,
    /// prime modulus. Must satisfy q == 1 mod 2n
    q: u64,
    /// log_2(n)
    log_n: u32,
    inv_n: u64,
    /// A 2N'th root of unity
    root: u64,
    /// Inverse of minimal root of unity
    inv_root: u64,
    /// powers of the minimal root of unity
    root_of_unity_powers: AVec<u64>,
    /// vector of floor(W * 2**32 / q), with W the root of unity powers
    precon32_root_of_unity_powers: AVec<u64>,
    /// vector of floor(W * 2**64 / q), with W the root of unity powers
    precon64_root_of_unity_powers: AVec<u64>,

    /// powers of the minimal root of unity adjusted for use in AVX512
    avx512_root_of_unity_powers: AVec<u64>,
    /// vector of floor(W * 2**32 / q), with W the AVX512 root of unity powers
    avx512_precon32_root_of_unity_powers: AVec<u64>,
    /// vector of floor(W * 2**52 / q), with W the AVX512 root of unity powers
    avx512_precon52_root_of_unity_powers: AVec<u64>,
    /// vector of floor(W * 2**64 / q), with W the AVX512 root of unity powers
    avx512_precon64_root_of_unity_powers: AVec<u64>,

    /// vector of floor(W * 2**32 / q), with W the inverse root of unity powers
    precon32_inv_root_of_unity_powers: AVec<u64>,
    /// vector of floor(W * 2**52 / q), with W the inverse root of unity powers
    precon52_inv_root_of_unity_powers: AVec<u64>,
    /// vector of floor(W * 2**64 / q), with W the inverse root of unity powers
    precon64_inv_root_of_unity_powers: AVec<u64>,

    inv_root_of_unity_powers: AVec<u64>,
}

impl HexlNttTable {
    pub fn n(&self) -> usize {
        self.n
    }

    pub fn q(&self) -> u64 {
        self.q
    }

    pub fn log_n(&self) -> u32 {
        self.log_n
    }

    pub fn root(&self) -> u64 {
        self.root
    }

    pub fn inv_root(&self) -> u64 {
        self.inv_root
    }

    pub fn root_of_unity_powers(&self) -> &[u64] {
        &self.root_of_unity_powers
    }

    pub fn precon32_root_of_unity_powers(&self) -> &[u64] {
        &self.precon32_root_of_unity_powers
    }

    pub fn precon64_root_of_unity_powers(&self) -> &[u64] {
        &self.precon64_root_of_unity_powers
    }

    pub fn avx512_root_of_unity_powers(&self) -> &[u64] {
        &self.avx512_root_of_unity_powers
    }

    pub fn avx512_precon32_root_of_unity_powers(&self) -> &[u64] {
        &self.avx512_precon32_root_of_unity_powers
    }

    pub fn avx512_precon52_root_of_unity_powers(&self) -> &[u64] {
        &self.avx512_precon52_root_of_unity_powers
    }

    pub fn avx512_precon64_root_of_unity_powers(&self) -> &[u64] {
        &self.avx512_precon64_root_of_unity_powers
    }

    pub fn precon32_inv_root_of_unity_powers(&self) -> &[u64] {
        &self.precon32_inv_root_of_unity_powers
    }

    pub fn precon52_inv_root_of_unity_powers(&self) -> &[u64] {
        &self.precon52_inv_root_of_unity_powers
    }

    pub fn precon64_inv_root_of_unity_powers(&self) -> &[u64] {
        &self.precon64_inv_root_of_unity_powers
    }

    pub fn inv_root_of_unity_powers(&self) -> &[u64] {
        &self.inv_root_of_unity_powers
    }
}

impl NttTable for HexlNttTable {
    type ValueT = u64;

    /// Initializes an NTT object with degree `2^log_n` and modulus `q`.
    ///
    /// ## Parameters
    /// - `log_n`: Also known as log(n) where n is size of the NTT transform.
    /// - `q`: Prime modulus. Must satisfy `q ≡ 1 (mod 2n)`.
    ///
    /// Performs pre-computation necessary for forward and inverse transforms.
    fn new<M>(log_n: u32, modulus: M) -> Result<Self, NttError<Self::ValueT>>
    where
        M: FieldContext<Self::ValueT>,
    {
        let q = modulus.value_unchecked();
        let n = 1usize << log_n;
        check_arguments(log_n, q);

        let root = <u64 as PrimitiveRoot>::try_minimal_primitive_root(log_n + 1, modulus)?;

        let to_root_type = |x: u64| -> Factor { Factor::new(x, q) };

        let root_factor = to_root_type(root);

        let mut power = root;

        let mut ordinal_root_powers = vec![0; n * 2];
        let mut iter = ordinal_root_powers.iter_mut();
        *iter.next().unwrap() = 1;
        *iter.next().unwrap() = root;
        for root_power in iter {
            power = root_factor.factor_mul_modulo(power, q);
            *root_power = power;
        }

        let inv_root = *ordinal_root_powers.last().unwrap();

        debug_assert_eq!(root_factor.factor_mul_modulo(inv_root, q), 1);

        let reverse_lsbs: Vec<usize> = (0..n).map(|i| i.reverse_lsbs(log_n)).collect();

        let mut root_of_unity_powers = avec![0; n];
        root_of_unity_powers[0] = 1;
        for (&root_power, &i) in ordinal_root_powers[0..n].iter().zip(reverse_lsbs.iter()) {
            root_of_unity_powers[i] = root_power;
        }

        let mut inv_root_of_unity_powers = avec![0; n];
        inv_root_of_unity_powers[0] = 1;
        for (&inv_root_power, &i) in ordinal_root_powers[n + 1..]
            .iter()
            .rev()
            .zip(reverse_lsbs.iter())
        {
            inv_root_of_unity_powers[i + 1] = inv_root_power;
        }

        let mut avx512_root_of_unity_powers = AVec::with_capacity(64, n / 8 + 3 * n / 2);

        // Duplicate each root of unity at indices [N/4, N/2].
        // These are the roots of unity used in the FwdNTT FwdT2 function
        // By creating these duplicates, we avoid extra permutations while loading the
        // roots of unity;
        let w2_roots: Vec<u64> = root_of_unity_powers[n / 4..n / 2]
            .iter()
            .flat_map(|&x| std::iter::repeat(x).take(2))
            .collect();

        // Duplicate each root of unity at indices [N/8, N/4].
        // These are the roots of unity used in the FwdNTT FwdT4 function
        // By creating these duplicates, we avoid extra permutations while loading the
        // roots of unity
        let w4_roots: Vec<u64> = root_of_unity_powers[n / 8..n / 4]
            .iter()
            .flat_map(|&x| std::iter::repeat(x).take(4))
            .collect();

        avx512_root_of_unity_powers.extend_from_slice(&root_of_unity_powers[0..n / 8]);
        avx512_root_of_unity_powers.extend_from_slice(&w4_roots);
        avx512_root_of_unity_powers.extend_from_slice(&w2_roots);
        avx512_root_of_unity_powers.extend_from_slice(&root_of_unity_powers[n / 2..]);

        let compute_barrett_vector = |values: &[u64], bit_shift: u32| {
            AVec::from_iter(
                64,
                values
                    .iter()
                    .map(|&value| MultiplyFactor::new(value, bit_shift, q).barrett_factor()),
            )
        };

        let precon32_root_of_unity_powers = compute_barrett_vector(&root_of_unity_powers, 32);
        let precon64_root_of_unity_powers = compute_barrett_vector(&root_of_unity_powers, 64);

        let avx512_precon52_root_of_unity_powers = if is_x86_feature_detected!("avx512ifma") {
            compute_barrett_vector(&avx512_root_of_unity_powers, 52)
        } else {
            AVec::new(0)
        };

        let (avx512_precon32_root_of_unity_powers, avx512_precon64_root_of_unity_powers) =
            if is_x86_feature_detected!("avx512dq") {
                (
                    compute_barrett_vector(&avx512_root_of_unity_powers, 32),
                    compute_barrett_vector(&avx512_root_of_unity_powers, 64),
                )
            } else {
                (AVec::new(0), AVec::new(0))
            };

        // 32-bit preconditioned inverse root of unity powers
        let precon32_inv_root_of_unity_powers =
            compute_barrett_vector(&inv_root_of_unity_powers, 32);

        // 52-bit preconditioned inverse root of unity powers
        let precon52_inv_root_of_unity_powers = if is_x86_feature_detected!("avx512ifma") {
            compute_barrett_vector(&inv_root_of_unity_powers, 52)
        } else {
            AVec::new(0)
        };

        let precon64_inv_root_of_unity_powers =
            compute_barrett_vector(&inv_root_of_unity_powers, 64);

        let inv_n = modulus.reduce_inv(n as u64);

        Ok(Self {
            n,
            q,
            log_n,
            inv_n,
            root,
            inv_root,
            root_of_unity_powers,
            precon32_root_of_unity_powers,
            precon64_root_of_unity_powers,
            avx512_root_of_unity_powers,
            avx512_precon32_root_of_unity_powers,
            avx512_precon52_root_of_unity_powers,
            avx512_precon64_root_of_unity_powers,
            precon32_inv_root_of_unity_powers,
            precon52_inv_root_of_unity_powers,
            precon64_inv_root_of_unity_powers,
            inv_root_of_unity_powers,
        })
    }

    fn poly_length(&self) -> usize {
        self.n
    }

    fn transform_inplace<
        S: primus_integer::RawData<Elem = Self::ValueT> + primus_integer::DataMut,
    >(
        &self,
        _poly: primus_poly::Polynomial<S>,
    ) -> primus_poly::NttPolynomial<S> {
        todo!()
    }

    fn inverse_transform_inplace<
        S: primus_integer::RawData<Elem = Self::ValueT> + primus_integer::DataMut,
    >(
        &self,
        _values: primus_poly::NttPolynomial<S>,
    ) -> primus_poly::Polynomial<S> {
        todo!()
    }

    fn lazy_transform_slice(&self, _poly: &mut [<Self as NttTable>::ValueT]) {
        todo!()
    }

    fn transform_slice(&self, _poly: &mut [<Self as NttTable>::ValueT]) {
        todo!()
    }

    fn lazy_inverse_transform_slice(&self, _values: &mut [<Self as NttTable>::ValueT]) {
        todo!()
    }

    fn inverse_transform_slice(&self, _values: &mut [<Self as NttTable>::ValueT]) {
        todo!()
    }

    fn transform_monomial(
        &self,
        _coeff: Self::ValueT,
        _degree: usize,
        _values: &mut [<Self as NttTable>::ValueT],
    ) {
        todo!()
    }

    fn transform_coeff_one_monomial(
        &self,
        _degree: usize,
        _values: &mut [<Self as NttTable>::ValueT],
    ) {
        todo!()
    }

    fn transform_coeff_minus_one_monomial(
        &self,
        _degree: usize,
        _values: &mut [<Self as NttTable>::ValueT],
    ) {
        todo!()
    }
}

impl HexlNttTable {
    pub fn compute_forward(
        &self,
        operand: &mut [u64],
        input_mod_factor: u64,
        output_mod_factor: u64,
    ) {
        debug_assert!(
            input_mod_factor == 1 || input_mod_factor == 2 || input_mod_factor == 4,
            "input_mod_factor must be 1, 2 or 4; got {input_mod_factor}",
        );
        debug_assert!(
            output_mod_factor == 1 || output_mod_factor == 4,
            "output_mod_factor must be 1 or 4; got {output_mod_factor}",
        );
        debug_assert_eq!(operand.len(), self.n, "operand length must be n={}", self.n);

        if is_x86_feature_detected!("avx512ifma") && self.q < MAX_FWD_IFMA_MODULUS && self.n >= 16 {
            let root_of_unity_powers = self.avx512_root_of_unity_powers();
            let precon_root_of_unity_powers = self.avx512_precon52_root_of_unity_powers();

            unsafe {
                forward_transform_to_bit_reverse_avx512::<IFMA_SHIFT_BITS>(
                    operand,
                    self.n,
                    self.q,
                    root_of_unity_powers,
                    precon_root_of_unity_powers,
                    input_mod_factor,
                    output_mod_factor,
                    0,
                    0,
                )
            };
            return;
        }

        if is_x86_feature_detected!("avx512dq") && self.n >= 16 {
            if self.q < MAX_FWD_32_MODULUS {
                let root_of_unity_powers = self.avx512_root_of_unity_powers();
                let precon_root_of_unity_powers = self.avx512_precon32_root_of_unity_powers();

                unsafe {
                    forward_transform_to_bit_reverse_avx512::<32>(
                        operand,
                        self.n,
                        self.q,
                        root_of_unity_powers,
                        precon_root_of_unity_powers,
                        input_mod_factor,
                        output_mod_factor,
                        0,
                        0,
                    )
                };
            } else {
                let root_of_unity_powers = self.avx512_root_of_unity_powers();
                let precon_root_of_unity_powers = self.avx512_precon64_root_of_unity_powers();

                unsafe {
                    forward_transform_to_bit_reverse_avx512::<DEFAULT_SHIFT_BITS>(
                        operand,
                        self.n,
                        self.q,
                        root_of_unity_powers,
                        precon_root_of_unity_powers,
                        input_mod_factor,
                        output_mod_factor,
                        0,
                        0,
                    )
                };
            }

            return;
        }

        let root_of_unity_powers = self.root_of_unity_powers();
        let precon_root_of_unity_powers = self.precon64_root_of_unity_powers();

        forward_transform_to_bit_reverse_radix2_inplace(
            operand,
            self.q,
            root_of_unity_powers,
            precon_root_of_unity_powers,
            output_mod_factor as u32,
        );
    }

    pub fn compute_inverse(
        &self,
        operand: &mut [u64],
        input_mod_factor: u64,
        output_mod_factor: u64,
    ) {
        debug_assert_eq!(operand.len(), self.n, "operand length must be n={}", self.n);
        debug_assert!(
            input_mod_factor == 1 || input_mod_factor == 2,
            "input_mod_factor must be 1 or 2; got {input_mod_factor}",
        );
        debug_assert!(
            output_mod_factor == 1 || output_mod_factor == 2,
            "output_mod_factor must be 1 or 2; got {output_mod_factor}",
        );

        if is_x86_feature_detected!("avx512ifma") && self.q < MAX_INV_IFMA_MODULUS && self.n >= 16 {
            let inv_root_of_unity_powers = self.inv_root_of_unity_powers();
            let precon_inv_root_of_unity_powers = self.precon52_inv_root_of_unity_powers();

            unsafe {
                inverse_transform_from_bit_reverse_avx512::<IFMA_SHIFT_BITS>(
                    operand,
                    self.q,
                    self.inv_n,
                    inv_root_of_unity_powers,
                    precon_inv_root_of_unity_powers,
                    input_mod_factor,
                    output_mod_factor,
                    0,
                    0,
                );
            }

            return;
        }

        if is_x86_feature_detected!("avx512dq") && self.n >= 16 {
            if self.q < MAX_INV_32_MODULUS {
                let inv_root_of_unity_powers = self.inv_root_of_unity_powers();
                let precon_inv_root_of_unity_powers = self.precon32_inv_root_of_unity_powers();
                unsafe {
                    inverse_transform_from_bit_reverse_avx512::<32>(
                        operand,
                        self.q,
                        self.inv_n,
                        inv_root_of_unity_powers,
                        precon_inv_root_of_unity_powers,
                        input_mod_factor,
                        output_mod_factor,
                        0,
                        0,
                    );
                }
            } else {
                let inv_root_of_unity_powers = self.inv_root_of_unity_powers();
                let precon_inv_root_of_unity_powers = self.precon64_inv_root_of_unity_powers();
                unsafe {
                    inverse_transform_from_bit_reverse_avx512::<DEFAULT_SHIFT_BITS>(
                        operand,
                        self.q,
                        self.inv_n,
                        inv_root_of_unity_powers,
                        precon_inv_root_of_unity_powers,
                        input_mod_factor,
                        output_mod_factor,
                        0,
                        0,
                    );
                }
            }
            return;
        }

        let _inv_root_of_unity_powers = self.inv_root_of_unity_powers();
        let _precon_inv_root_of_unity_powers = self.precon64_inv_root_of_unity_powers();

        todo!()
    }
}

#[cfg(test)]
mod tests {
    use primus_modulus::BarrettModulus;

    use super::*;

    #[test]
    fn test_hexl() {
        let _table = HexlNttTable::new(10, BarrettModulus::new(132120577)).unwrap();
    }
}
