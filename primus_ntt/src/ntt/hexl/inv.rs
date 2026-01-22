use core::arch::x86_64::*;

use primus_factor::MultiplyFactor;

use crate::ntt::hexl::internal::{check_arguments, max_inv_modulus};

use super::utils::*;

/// The Harvey butterfly: assume `X`, `Y` in `[0, 2q)`, and return `X'`, `Y'` in
/// `[0, 2q)` such that
/// `X' = X + Y (mod q)`, `Y' = W * (X - Y) (mod q)`.
///
/// # Parameters
/// - `X`: Input/output representing 8 64-bit signed integers in SIMD form.
/// - `Y`: Input/output representing 8 64-bit signed integers in SIMD form.
/// - `W`: Root of unity representing 8 64-bit signed integers in SIMD form.
/// - `W_precon`: Preconditioned `W` for `BitShift`-bit Barrett reduction.
/// - `neg_modulus`: Negative modulus, i.e. `(-q)` represented as 8 64-bit signed integers in SIMD form.
/// - `twice_modulus`: Twice the modulus, i.e. `2*q` represented as 8 64-bit signed integers in SIMD form.
/// - `input_less_than_mod`: If `true`, assumes `X, Y < q`. Otherwise assumes `X, Y < 2*q`.
///
/// # Details
/// See Algorithm 3 of https://arxiv.org/pdf/1205.2926.pdf
pub unsafe fn inv_butterfly<const BIT_SHIFT: u32, const INPUT_LESS_THAN_MOD: bool>(
    x: &mut __m512i,
    y: &mut __m512i,
    w: __m512i,
    w_precon: __m512i,
    neg_modulus: __m512i,
    twice_modulus: __m512i,
) {
    // Compute T first to allow in-place update of X
    let y_minus_2q = unsafe { _mm512_sub_epi64(*y, twice_modulus) };
    let t = unsafe { _mm512_sub_epi64(*x, y_minus_2q) };

    if INPUT_LESS_THAN_MOD {
        // No need for modulus reduction, since inputs are in [0, q)
        *x = unsafe { _mm512_add_epi64(*x, *y) };
    } else {
        // Algorithm 3 computes (X >= 2q) ? (X - 2q) : X
        // We instead compute (X - 2q >= 0) ? (X - 2q) : X
        // This allows us to use the faster _mm512_movepi64_mask rather than
        // _mm512_cmp_epu64_mask to create the mask.
        unsafe {
            *x = _mm512_add_epi64(*x, y_minus_2q);
            let sign_bits = _mm512_movepi64_mask(*x);
            *x = _mm512_mask_add_epi64(*x, sign_bits, *x, twice_modulus);
        }
    }

    if BIT_SHIFT == 32 {
        unsafe {
            let mut q = _mm512_hexl_mullo_epi_64(w_precon, t);
            q = _mm512_srli_epi64(q, 32);
            let q_p = _mm512_hexl_mullo_epi_64(q, neg_modulus);
            *y = _mm512_hexl_mullo_add_lo_epi_64(q_p, w, t);
        }
    } else if BIT_SHIFT == 52 {
        unsafe {
            let q = _mm512_hexl_mulhi_epi_52(w_precon, t);
            let q_p = _mm512_hexl_mullo_epi_52(q, neg_modulus);
            *y = _mm512_hexl_mullo_add_lo_epi_52(q_p, w, t);
        }
    } else if BIT_SHIFT == 64 {
        unsafe {
            // Perform approximate computation of Q, as described in page 7 of
            // https://arxiv.org/pdf/2003.04510.pdf
            let q = _mm512_hexl_mulhi_approx_epi_64(w_precon, t);
            let q_p = _mm512_hexl_mullo_epi_64(q, neg_modulus);
            // Compute Y in range [0, 4q)
            *y = _mm512_hexl_mullo_add_lo_epi_64(q_p, w, t);
            // Reduce Y to range [0, 2q)
            *y = _mm512_hexl_small_mod_epu64_2(*y, twice_modulus);
        }
    } else {
        debug_assert!(false, "Invalid BitShift {BIT_SHIFT}")
    }
}

pub fn inv_t1<const BIT_SHIFT: u32, const INPUT_LESS_THAN_MOD: bool>(
    x: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    w: &[u64],
    w_precon: &[u64],
) {
    // n >= 16 and n is power of 2
    unsafe {
        for ((chunk, w_chunk), w_precon_chunk) in x
            .as_chunks_unchecked_mut::<16>()
            .iter_mut()
            .zip(w.as_chunks_unchecked::<8>())
            .zip(w_precon.as_chunks_unchecked::<8>())
        {
            let (mut v_x, mut v_y) = load_inv_interleaved_t1(chunk);

            let v_w = _mm512_loadu_si512(w_chunk.as_ptr().cast());
            let v_w_precon = _mm512_loadu_si512(w_precon_chunk.as_ptr().cast());

            inv_butterfly::<BIT_SHIFT, INPUT_LESS_THAN_MOD>(
                &mut v_x,
                &mut v_y,
                v_w,
                v_w_precon,
                v_neg_modulus,
                v_twice_mod,
            );

            let v_x_pt: *mut __m512i = chunk.as_mut_ptr().cast();

            _mm512_storeu_si512(v_x_pt, v_x);
            _mm512_storeu_si512(v_x_pt.add(1), v_y);
        }
    }
}

pub fn inv_t2<const BIT_SHIFT: u32>(
    x: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    w: &[u64],
    w_precon: &[u64],
) {
    // n >= 16 and n is power of 2
    unsafe {
        for ((chunk, w_chunk), w_precon_chunk) in x
            .as_chunks_unchecked_mut::<16>()
            .iter_mut()
            .zip(w.as_chunks_unchecked::<4>())
            .zip(w_precon.as_chunks_unchecked::<4>())
        {
            let (mut v_x, mut v_y) = load_inv_interleaved_t2(chunk);

            let v_w = load_w_op_t2(w_chunk);
            let v_w_precon = load_w_op_t2(w_precon_chunk);

            inv_butterfly::<BIT_SHIFT, false>(
                &mut v_x,
                &mut v_y,
                v_w,
                v_w_precon,
                v_neg_modulus,
                v_twice_mod,
            );

            let v_x_pt: *mut __m512i = chunk.as_mut_ptr().cast();

            _mm512_storeu_si512(v_x_pt, v_x);
            _mm512_storeu_si512(v_x_pt.add(1), v_y);
        }
    }
}

pub fn inv_t4<const BIT_SHIFT: u32>(
    x: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    w: &[u64],
    w_precon: &[u64],
) {
    // n >= 16 and n is power of 2
    unsafe {
        for ((chunk, w_chunk), w_precon_chunk) in x
            .as_chunks_unchecked_mut::<16>()
            .iter_mut()
            .zip(w.as_chunks_unchecked::<2>())
            .zip(w_precon.as_chunks_unchecked::<2>())
        {
            let (mut v_x, mut v_y) = load_inv_interleaved_t4(chunk);

            let v_w = load_w_op_t4(w_chunk);
            let v_w_precon = load_w_op_t4(w_precon_chunk);

            inv_butterfly::<BIT_SHIFT, false>(
                &mut v_x,
                &mut v_y,
                v_w,
                v_w_precon,
                v_neg_modulus,
                v_twice_mod,
            );

            write_inv_interleaved_t4(v_x, v_y, chunk);
        }
    }
}

pub fn inv_t8<const BIT_SHIFT: u32>(
    operand: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    t: usize,
    w: &[u64],
    w_precon: &[u64],
) {
    let mut w_iter = w.iter().copied();
    let mut w_precon_iter = w_precon.iter().copied();

    // assume 8 | t
    for chunk in operand.chunks_exact_mut(t << 1) {
        let (x, y) = unsafe { chunk.split_at_mut_unchecked(t) };

        unsafe {
            let v_w = _mm512_set1_epi64(w_iter.next().unwrap() as i64);
            let v_w_precon = _mm512_set1_epi64(w_precon_iter.next().unwrap() as i64);

            for (x_chunk, y_chunk) in x
                .as_chunks_unchecked_mut::<8>()
                .iter_mut()
                .zip(y.as_chunks_unchecked_mut::<8>())
            {
                let mut v_x = _mm512_loadu_si512(x_chunk.as_ptr().cast());
                let mut v_y = _mm512_loadu_si512(y_chunk.as_ptr().cast());

                inv_butterfly::<BIT_SHIFT, false>(
                    &mut v_x,
                    &mut v_y,
                    v_w,
                    v_w_precon,
                    v_neg_modulus,
                    v_twice_mod,
                );

                _mm512_storeu_si512(x_chunk.as_mut_ptr().cast(), v_x);
                _mm512_storeu_si512(y_chunk.as_mut_ptr().cast(), v_y);
            }
        }
    }
}

pub unsafe fn inverse_transform_from_bit_reverse_avx512<const BIT_SHIFT: u32>(
    operand: &mut [u64],
    modulus: u64,
    inv_n: u64,
    inv_root_of_unity_powers: &[u64],
    precon_inv_root_of_unity_powers: &[u64],
    input_mod_factor: u64,
    output_mod_factor: u64,
    recursion_depth: usize,
    recursion_half: usize,
) {
    let n = operand.len();

    check_arguments(n, modulus);
    debug_assert!(
        n >= 16,
        "inverse_transform_from_bit_reverse_avx512 doesn't support small transforms. Need n >= 16, got n = {n}"
    );
    debug_assert!(
        modulus < max_inv_modulus(BIT_SHIFT),
        "modulus {modulus} too large for BitShift {BIT_SHIFT} => maximum value {}",
        max_inv_modulus(BIT_SHIFT)
    );

    debug_assert!(
        input_mod_factor == 1 || input_mod_factor == 2,
        "input_mod_factor must be 1 or 2; got {input_mod_factor}",
    );
    debug_assert!(
        output_mod_factor == 1 || output_mod_factor == 2,
        "output_mod_factor must be 1 or 2; got {output_mod_factor}",
    );

    let twice_mod = modulus << 1;
    let v_modulus = unsafe { _mm512_set1_epi64(modulus as i64) };
    let v_neg_modulus = unsafe { _mm512_set1_epi64(-(modulus as i64)) };
    let v_twice_mod = unsafe { _mm512_set1_epi64(twice_mod as i64) };

    let mut t = 1;
    let mut m = n >> 1;
    let mut w_idx = 1 + m * recursion_half;

    const BASE_NTT_SIZE: usize = 1024;

    if n <= BASE_NTT_SIZE {
        // Perform breadth-first InvNTT

        // Extract t=1, t=2, t=4 loops separately
        {
            // t = 1
            let w = &inv_root_of_unity_powers[w_idx..w_idx + m];
            let w_precon = &precon_inv_root_of_unity_powers[w_idx..w_idx + m];

            if input_mod_factor == 1 && recursion_depth == 0 {
                inv_t1::<BIT_SHIFT, true>(operand, v_neg_modulus, v_twice_mod, w, w_precon);
            } else {
                inv_t1::<BIT_SHIFT, false>(operand, v_neg_modulus, v_twice_mod, w, w_precon);
            }

            t <<= 1;
            m >>= 1;
            let mut w_idx_delta = m * ((1 << (recursion_depth + 1)) - recursion_half);
            w_idx += w_idx_delta;

            // t = 2
            let w = &inv_root_of_unity_powers[w_idx..w_idx + m];
            let w_precon = &precon_inv_root_of_unity_powers[w_idx..w_idx + m];
            inv_t2::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, w, w_precon);

            t <<= 1;
            m >>= 1;
            w_idx_delta >>= 1;
            w_idx += w_idx_delta;

            // t = 4
            let w = &inv_root_of_unity_powers[w_idx..w_idx + m];
            let w_precon = &precon_inv_root_of_unity_powers[w_idx..w_idx + m];
            inv_t4::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, w, w_precon);

            t <<= 1;
            m >>= 1;
            w_idx_delta >>= 1;
            w_idx += w_idx_delta;

            // t >= 8
            while m > 1 {
                let w = &inv_root_of_unity_powers[w_idx..w_idx + m];
                let w_precon = &precon_inv_root_of_unity_powers[w_idx..w_idx + m];
                inv_t8::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, t, w, w_precon);
                t <<= 1;
                m >>= 1;
                w_idx_delta >>= 1;
                w_idx += w_idx_delta;
            }
        }
    } else {
        unsafe {
            let (left, right) = operand.split_at_mut_unchecked(n / 2);
            inverse_transform_from_bit_reverse_avx512::<BIT_SHIFT>(
                left,
                modulus,
                inv_n,
                inv_root_of_unity_powers,
                precon_inv_root_of_unity_powers,
                input_mod_factor,
                output_mod_factor,
                recursion_depth + 1,
                recursion_half * 2,
            );
            inverse_transform_from_bit_reverse_avx512::<BIT_SHIFT>(
                right,
                modulus,
                inv_n,
                inv_root_of_unity_powers,
                precon_inv_root_of_unity_powers,
                input_mod_factor,
                output_mod_factor,
                recursion_depth + 1,
                recursion_half * 2 + 1,
            );
        }

        let mut w_idx_delta = m * ((1 << (recursion_depth + 1)) - recursion_half);
        while m > 2 {
            t <<= 1;
            m >>= 1;
            w_idx_delta >>= 1;
            w_idx += w_idx_delta;
        }
        if m == 2 {
            let w = &inv_root_of_unity_powers[w_idx..w_idx + m];
            let w_precon = &precon_inv_root_of_unity_powers[w_idx..w_idx + m];
            inv_t8::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, t, w, w_precon);
            // t <<= 1;
            // m >>= 1;
            w_idx_delta >>= 1;
            w_idx += w_idx_delta;
        }
    }

    // Final loop through data
    if recursion_depth == 0 {
        let w = inv_root_of_unity_powers[w_idx];
        let mf_inv_n = MultiplyFactor::new(inv_n, BIT_SHIFT, modulus);
        let inv_n_prime = mf_inv_n.barrett_factor();

        let inv_n_w = mf_inv_n.mul_modulo::<BIT_SHIFT>(w, modulus);
        let mf_inv_n_w = MultiplyFactor::new(inv_n_w, BIT_SHIFT, modulus);
        let inv_n_w_prime = mf_inv_n_w.barrett_factor();

        unsafe {
            let (x, y) = operand.split_at_mut_unchecked(n / 2);

            let v_inv_n = _mm512_set1_epi64(inv_n as i64);
            let v_inv_n_prime = _mm512_set1_epi64(inv_n_prime as i64);
            let v_inv_n_w = _mm512_set1_epi64(inv_n_w as i64);
            let v_inv_n_w_prime = _mm512_set1_epi64(inv_n_w_prime as i64);

            // Merge final InvNTT loop with modulus reduction baked-in
            for (x_chunk, y_chunk) in x
                .as_chunks_unchecked_mut::<8>()
                .iter_mut()
                .zip(y.as_chunks_unchecked_mut::<8>())
            {
                let mut v_x = _mm512_loadu_si512(x_chunk.as_ptr().cast());
                let mut v_y = _mm512_loadu_si512(y_chunk.as_ptr().cast());

                // Slightly different from regular InvButterfly because different W is
                // used for X and Y
                let y_minus_2q = _mm512_sub_epi64(v_y, v_twice_mod);
                let x_plus_y_mod2q = _mm512_hexl_small_add_mod_epi64(v_x, v_y, v_twice_mod);
                // T = *X + twice_mod - *Y
                let t = _mm512_sub_epi64(v_x, y_minus_2q);

                if BIT_SHIFT == 32 {
                    let mut q1 = _mm512_hexl_mullo_epi_64(v_inv_n_prime, x_plus_y_mod2q);
                    q1 = _mm512_srli_epi64::<32>(q1);
                    // X = inv_N * X_plus_Y_mod2q - Q1 * modulus;
                    let inv_n_tx = _mm512_hexl_mullo_epi_64(v_inv_n, x_plus_y_mod2q);
                    v_x = _mm512_hexl_mullo_add_lo_epi_64(inv_n_tx, q1, v_neg_modulus);

                    let mut q2 = _mm512_hexl_mullo_epi_64(v_inv_n_w_prime, t);
                    q2 = _mm512_srli_epi64::<32>(q2);

                    // Y = inv_N_W * T - Q2 * modulus;
                    let inv_n_w_t = _mm512_hexl_mullo_epi_64(v_inv_n_w, t);
                    v_y = _mm512_hexl_mullo_add_lo_epi_64(inv_n_w_t, q2, v_neg_modulus);
                } else if BIT_SHIFT == 52 {
                    let q1 = _mm512_hexl_mulhi_epi_52(v_inv_n_prime, x_plus_y_mod2q);
                    // X = inv_N * X_plus_Y_mod2q - Q1 * modulus;
                    let inv_n_tx = _mm512_hexl_mullo_epi_52(v_inv_n, x_plus_y_mod2q);
                    v_x = _mm512_hexl_mullo_add_lo_epi_52(inv_n_tx, q1, v_neg_modulus);

                    let q2 = _mm512_hexl_mulhi_epi_52(v_inv_n_w_prime, t);
                    // Y = inv_N_W * T - Q2 * modulus;
                    let inv_n_w_t = _mm512_hexl_mullo_epi_52(v_inv_n_w, t);
                    v_y = _mm512_hexl_mullo_add_lo_epi_52(inv_n_w_t, q2, v_neg_modulus);
                } else if BIT_SHIFT == 64 {
                    let q1 = _mm512_hexl_mulhi_epi_64(v_inv_n_prime, x_plus_y_mod2q);
                    // X = inv_N * X_plus_Y_mod2q - Q1 * modulus;
                    let inv_n_tx = _mm512_hexl_mullo_epi_64(v_inv_n, x_plus_y_mod2q);
                    v_x = _mm512_hexl_mullo_add_lo_epi_64(inv_n_tx, q1, v_neg_modulus);

                    let q2 = _mm512_hexl_mulhi_epi_64(v_inv_n_w_prime, t);
                    // Y = inv_N_W * T - Q2 * modulus;
                    let inv_n_w_t = _mm512_hexl_mullo_epi_64(v_inv_n_w, t);
                    v_y = _mm512_hexl_mullo_add_lo_epi_64(inv_n_w_t, q2, v_neg_modulus);
                } else {
                    debug_assert!(false);
                }

                if output_mod_factor == 1 {
                    // Modulus reduction from [0, 2q), to [0, q)
                    v_x = _mm512_hexl_small_mod_epu64_2(v_x, v_modulus);
                    v_y = _mm512_hexl_small_mod_epu64_2(v_y, v_modulus);
                }

                _mm512_storeu_si512(x_chunk.as_mut_ptr().cast(), v_x);
                _mm512_storeu_si512(y_chunk.as_mut_ptr().cast(), v_y);
            }
        }
    }
}
