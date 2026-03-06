use core::arch::x86_64::*;

use crate::ntt::hexl::internal::{check_arguments, max_fwd_modulus};

use super::utils::*;

/// The Harvey butterfly: assume `X`, `Y` in `[0, 4q)`, and return `X'`, `Y'` in
/// `[0, 4q)` such that
/// `X' = X + W*Y (mod q)`, `Y' = X - W*Y (mod q)`.
///
/// # Parameters
/// - `X`: Input/output representing 8 64-bit signed integers in SIMD form.
/// - `Y`: Input/output representing 8 64-bit signed integers in SIMD form.
/// - `W`: Root of unity represented as 8 64-bit signed integers in SIMD form.
/// - `W_precon`: Preconditioned `W` for `BitShift`-bit Barrett reduction.
/// - `neg_modulus`: Negative modulus, i.e. `(-q)` represented as 8 64-bit signed integers in SIMD form.
/// - `twice_modulus`: Twice the modulus, i.e. `2*q` represented as 8 64-bit signed integers in SIMD form.
/// - `input_less_than_mod`: If `true`, assumes `X, Y < q`. Otherwise assumes `X, Y < 4*q`.
///
/// # Details
/// See Algorithm 4 of https://arxiv.org/pdf/1205.2926.pdf
pub unsafe fn fwd_butterfly<const BIT_SHIFT: u32, const INPUT_LESS_THAN_MOD: bool>(
    x: &mut __m512i,
    y: &mut __m512i,
    w: __m512i,
    w_precon: __m512i,
    neg_modulus: __m512i,
    twice_modulus: __m512i,
) {
    if !INPUT_LESS_THAN_MOD {
        *x = unsafe { _mm512_hexl_small_mod_epu64_2(*x, twice_modulus) };
    }

    let mut t: __m512i;
    if BIT_SHIFT == 32 {
        unsafe {
            let mut q = _mm512_hexl_mullo_epi_64(w_precon, *y);
            q = _mm512_srli_epi64::<32>(q);
            let w_y = _mm512_hexl_mullo_epi_64(w, *y);
            t = _mm512_hexl_mullo_add_lo_epi_64(w_y, q, neg_modulus);
        }
    } else if BIT_SHIFT == 52 {
        unsafe {
            let q = _mm512_hexl_mulhi_epi_52(w_precon, *y);
            let w_y = _mm512_hexl_mullo_epi_52(w, *y);
            t = _mm512_hexl_mullo_add_lo_epi_52(w_y, q, neg_modulus);
        }
    } else if BIT_SHIFT == 64 {
        // Perform approximate computation of Q, as described in page 7 of
        // https://arxiv.org/pdf/2003.04510.pdf
        unsafe {
            let q = _mm512_hexl_mulhi_approx_epi_64(w_precon, *y);
            let w_y = _mm512_hexl_mullo_epi_64(w, *y);
            // Compute T in range [0, 4q)
            t = _mm512_hexl_mullo_add_lo_epi_64(w_y, q, neg_modulus);
            // Reduce T to range [0, 2q)
            t = _mm512_hexl_small_mod_epu64_2(t, twice_modulus);
        }
    } else {
        panic!("Invalid BitShift {BIT_SHIFT}");
    }

    unsafe {
        let twice_mod_minus_t = _mm512_sub_epi64(twice_modulus, t);
        *y = _mm512_add_epi64(*x, twice_mod_minus_t);
        *x = _mm512_add_epi64(*x, t);
    }
}

pub fn fwd_t1<const BIT_SHIFT: u32>(
    operand: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    w: &[u64],
    w_precon: &[u64],
) {
    unsafe {
        for ((x, v_w), v_w_precon) in operand
            .as_chunks_unchecked_mut::<16>()
            .iter_mut()
            .zip(w.as_chunks_unchecked::<8>())
            .zip(w_precon.as_chunks_unchecked::<8>())
        {
            let (mut v_x, mut v_y) = load_fwd_interleaved_t1(x);

            let v_w = _mm512_loadu_si512(v_w.as_ptr().cast());
            let v_w_precon = _mm512_loadu_si512(v_w_precon.as_ptr().cast());

            fwd_butterfly::<BIT_SHIFT, false>(
                &mut v_x,
                &mut v_y,
                v_w,
                v_w_precon,
                v_neg_modulus,
                v_twice_mod,
            );

            write_fwd_interleaved_t1(v_x, v_y, x);
        }
    }
}

pub fn fwd_t2<const BIT_SHIFT: u32>(
    operand: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    w: &[u64],
    w_precon: &[u64],
) {
    unsafe {
        for ((x, v_w), v_w_precon) in operand
            .as_chunks_unchecked_mut::<16>()
            .iter_mut()
            .zip(w.as_chunks_unchecked::<8>())
            .zip(w_precon.as_chunks_unchecked::<8>())
        {
            let (mut v_x, mut v_y) = load_fwd_interleaved_t2(x);

            let v_w = _mm512_loadu_si512(v_w.as_ptr().cast());
            let v_w_precon = _mm512_loadu_si512(v_w_precon.as_ptr().cast());

            fwd_butterfly::<BIT_SHIFT, false>(
                &mut v_x,
                &mut v_y,
                v_w,
                v_w_precon,
                v_neg_modulus,
                v_twice_mod,
            );

            let v_x_pt: *mut __m512i = x.as_mut_ptr().cast();

            _mm512_storeu_si512(v_x_pt, v_x);
            _mm512_storeu_si512(v_x_pt.add(1), v_y);
        }
    }
}

pub fn fwd_t4<const BIT_SHIFT: u32>(
    operand: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    w: &[u64],
    w_precon: &[u64],
) {
    unsafe {
        for ((x, v_w), v_w_precon) in operand
            .as_chunks_unchecked_mut::<16>()
            .iter_mut()
            .zip(w.as_chunks_unchecked::<8>())
            .zip(w_precon.as_chunks_unchecked::<8>())
        {
            let (mut v_x, mut v_y) = load_fwd_interleaved_t4(x);

            let v_w = _mm512_loadu_si512(v_w.as_ptr().cast());
            let v_w_precon = _mm512_loadu_si512(v_w_precon.as_ptr().cast());

            fwd_butterfly::<BIT_SHIFT, false>(
                &mut v_x,
                &mut v_y,
                v_w,
                v_w_precon,
                v_neg_modulus,
                v_twice_mod,
            );

            let v_x_pt: *mut __m512i = x.as_mut_ptr().cast();

            _mm512_storeu_si512(v_x_pt, v_x);
            _mm512_storeu_si512(v_x_pt.add(1), v_y);
        }
    }
}

pub fn fwd_t8_inplace<const BIT_SHIFT: u32, const INPUT_LESS_THAN_MOD: bool>(
    operand: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    t: usize,
    w: &[u64],
    w_precon: &[u64],
) {
    let mut w_iter = w.iter().copied();
    let mut w_precon_iter = w_precon.iter().copied();

    for chunk in operand.chunks_exact_mut(t << 1) {
        unsafe {
            let (x, y) = chunk.split_at_mut_unchecked(t);

            // Weights and weights' preconditions
            let v_w: __m512i = _mm512_set1_epi64(w_iter.next().unwrap() as i64);
            let v_w_precon: __m512i = _mm512_set1_epi64(w_precon_iter.next().unwrap() as i64);

            for (x_chunk, y_chunk) in x
                .as_chunks_unchecked_mut::<8>()
                .iter_mut()
                .zip(y.as_chunks_unchecked_mut::<8>())
            {
                let mut v_x = _mm512_loadu_si512(x_chunk.as_ptr().cast());
                let mut v_y = _mm512_loadu_si512(y_chunk.as_ptr().cast());

                fwd_butterfly::<BIT_SHIFT, INPUT_LESS_THAN_MOD>(
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

#[allow(clippy::too_many_arguments)]
pub unsafe fn forward_transform_to_bit_reverse_avx512<const BIT_SHIFT: u32>(
    operand: &mut [u64],
    modulus: u64,
    root_of_unity_powers: &[u64],
    precon_root_of_unity_powers: &[u64],
    input_mod_factor: u64,
    output_mod_factor: u64,
    recursion_depth: usize,
    recursion_half: usize,
) {
    let n = operand.len();

    check_arguments(n, modulus);
    debug_assert!(
        modulus < max_fwd_modulus(BIT_SHIFT),
        "modulus {modulus} too large for BitShift {BIT_SHIFT} => maximum value {}",
        max_fwd_modulus(BIT_SHIFT)
    );
    debug_assert!(
        n >= 16,
        "Don't support small transforms. Need n >= 16, got n = {n}"
    );
    debug_assert!(
        input_mod_factor == 1 || input_mod_factor == 2 || input_mod_factor == 4,
        "input_mod_factor must be 1, 2, or 4; got {input_mod_factor}"
    );
    debug_assert!(
        output_mod_factor == 1 || output_mod_factor == 4,
        "output_mod_factor must be 1 or 4; got {output_mod_factor}"
    );

    let twice_mod = modulus << 1;

    let v_modulus = unsafe { _mm512_set1_epi64(modulus as i64) };
    let v_neg_modulus = unsafe { _mm512_set1_epi64(-(modulus as i64)) };
    let v_twice_mod = unsafe { _mm512_set1_epi64(twice_mod as i64) };

    const BASE_NTT_SIZE: usize = 1024;

    if n <= BASE_NTT_SIZE {
        // Perform breadth-first NTT
        let mut t = n >> 1;
        let mut m = 1;
        let mut w_idx = (m << recursion_depth) + (recursion_half * m);

        // First iteration assumes input in [0,p)
        if m < (n >> 3) {
            let w = &root_of_unity_powers[w_idx..w_idx + m];
            let w_precon = &precon_root_of_unity_powers[w_idx..w_idx + m];

            if input_mod_factor <= 2 && recursion_depth == 0 {
                fwd_t8_inplace::<BIT_SHIFT, true>(
                    operand,
                    v_neg_modulus,
                    v_twice_mod,
                    t,
                    w,
                    w_precon,
                );
            } else {
                fwd_t8_inplace::<BIT_SHIFT, false>(
                    operand,
                    v_neg_modulus,
                    v_twice_mod,
                    t,
                    w,
                    w_precon,
                );
            }

            t >>= 1;
            m <<= 1;
            w_idx <<= 1;
        }

        while m < (n >> 3) {
            let w = &root_of_unity_powers[w_idx..w_idx + m];
            let w_precon = &precon_root_of_unity_powers[w_idx..w_idx + m];

            fwd_t8_inplace::<BIT_SHIFT, false>(operand, v_neg_modulus, v_twice_mod, t, w, w_precon);

            t >>= 1;
            m <<= 1;
            w_idx <<= 1;
        }

        // Do T=4, T=2, T=1 separately
        {
            // Correction step needed due to extra copies of roots of unity in the
            // AVX512 vectors loaded for FwdT2 and FwdT4
            let compute_new_w_idx = |idx: usize| {
                // Originally, from root of unity vector index to loop:
                // [0, N/8) => FwdT8
                // [N/8, N/4) => FwdT4
                // [N/4, N/2) => FwdT2
                // [N/2, N) => FwdT1
                // The new mapping from AVX512 root of unity vector index to loop:
                // [0, N/8) => FwdT8
                // [N/8, 5N/8) => FwdT4
                // [5N/8, 9N/8) => FwdT2
                // [9N/8, 13N/8) => FwdT1
                let n = n << recursion_depth;

                // FwdT8 range
                if idx <= n / 8 {
                    return idx;
                }
                // FwdT4 range
                if idx <= n / 4 {
                    return (idx - n / 8) * 4 + (n / 8);
                }
                // FwdT2 range
                if idx <= n / 2 {
                    return (idx - n / 4) * 2 + (5 * n / 8);
                }
                // FwdT1 range
                idx + (5 * n / 8)
            };

            let mut new_w_idx = compute_new_w_idx(w_idx);
            let mut w = &root_of_unity_powers[new_w_idx..new_w_idx + m * 4];
            let mut w_precon = &precon_root_of_unity_powers[new_w_idx..new_w_idx + m * 4];
            fwd_t4::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, w, w_precon);

            m <<= 1;
            w_idx <<= 1;
            new_w_idx = compute_new_w_idx(w_idx);
            w = &root_of_unity_powers[new_w_idx..new_w_idx + m * 2];
            w_precon = &precon_root_of_unity_powers[new_w_idx..new_w_idx + m * 2];
            fwd_t2::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, w, w_precon);

            m <<= 1;
            w_idx <<= 1;
            new_w_idx = compute_new_w_idx(w_idx);
            w = &root_of_unity_powers[new_w_idx..new_w_idx + m];
            w_precon = &precon_root_of_unity_powers[new_w_idx..new_w_idx + m];
            fwd_t1::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, w, w_precon);
        }

        if output_mod_factor == 1 {
            // n power of two at least 8 => n divisible by 8
            unsafe {
                for chunk in operand.as_chunks_unchecked_mut::<8>() {
                    let mut v_x = _mm512_loadu_si512(chunk.as_ptr().cast());

                    // Reduce from [0, 4q) to [0, q)
                    v_x = _mm512_hexl_small_mod_epu64_2(v_x, v_twice_mod);
                    v_x = _mm512_hexl_small_mod_epu64_2(v_x, v_modulus);

                    _mm512_storeu_si512(chunk.as_mut_ptr().cast(), v_x);
                }
            }
        }
    } else {
        // Perform depth-first NTT via recursive call
        let t = n >> 1;
        let w_idx = (1 << recursion_depth) + recursion_half;
        let w = &root_of_unity_powers[w_idx..w_idx + 1];
        let w_precon = &precon_root_of_unity_powers[w_idx..w_idx + 1];

        fwd_t8_inplace::<BIT_SHIFT, false>(operand, v_neg_modulus, v_twice_mod, t, w, w_precon);

        unsafe {
            let (left, right) = operand.split_at_mut_unchecked(n / 2);

            forward_transform_to_bit_reverse_avx512::<BIT_SHIFT>(
                left,
                modulus,
                root_of_unity_powers,
                precon_root_of_unity_powers,
                input_mod_factor,
                output_mod_factor,
                recursion_depth + 1,
                recursion_half * 2,
            );

            forward_transform_to_bit_reverse_avx512::<BIT_SHIFT>(
                right,
                modulus,
                root_of_unity_powers,
                precon_root_of_unity_powers,
                input_mod_factor,
                output_mod_factor,
                recursion_depth + 1,
                recursion_half * 2 + 1,
            );
        }
    }
}
