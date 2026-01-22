use core::arch::x86_64::*;

use crate::ntt::hexl::internal::{check_arguments, max_fwd_modulus};

use super::utils::*;

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

#[target_feature(enable = "avx512dq")]
pub fn fwd_butterfly_32<const INPUT_LESS_THAN_MOD: bool>(
    x: &mut __m512i,
    y: &mut __m512i,
    w: __m512i,
    w_precon: __m512i,
    neg_modulus: __m512i,
    twice_modulus: __m512i,
) {
    if !INPUT_LESS_THAN_MOD {
        *x = _mm512_hexl_small_mod_epu64_2(*x, twice_modulus);
    }

    let mut q = _mm512_hexl_mullo_epi_64(w_precon, *y);
    q = _mm512_srli_epi64::<32>(q);
    let w_y = _mm512_hexl_mullo_epi_64(w, *y);
    let t = _mm512_hexl_mullo_add_lo_epi_64(w_y, q, neg_modulus);

    let twice_mod_minus_t = _mm512_sub_epi64(twice_modulus, t);
    *y = _mm512_add_epi64(*x, twice_mod_minus_t);
    *x = _mm512_add_epi64(*x, t);
}

#[target_feature(enable = "avx512ifma")]
pub fn fwd_butterfly_52<const INPUT_LESS_THAN_MOD: bool>(
    x: &mut __m512i,
    y: &mut __m512i,
    w: __m512i,
    w_precon: __m512i,
    neg_modulus: __m512i,
    twice_modulus: __m512i,
) {
    if !INPUT_LESS_THAN_MOD {
        *x = _mm512_hexl_small_mod_epu64_2(*x, twice_modulus);
    }

    let q = _mm512_hexl_mulhi_epi_52(w_precon, *y);
    let w_y = _mm512_hexl_mullo_epi_52(w, *y);
    let t = _mm512_hexl_mullo_add_lo_epi_52(w_y, q, neg_modulus);

    let twice_mod_minus_t = _mm512_sub_epi64(twice_modulus, t);
    *y = _mm512_add_epi64(*x, twice_mod_minus_t);
    *x = _mm512_add_epi64(*x, t);
}

#[target_feature(enable = "avx512dq")]
pub fn fwd_butterfly_64<const INPUT_LESS_THAN_MOD: bool>(
    x: &mut __m512i,
    y: &mut __m512i,
    w: __m512i,
    w_precon: __m512i,
    neg_modulus: __m512i,
    twice_modulus: __m512i,
) {
    if !INPUT_LESS_THAN_MOD {
        *x = _mm512_hexl_small_mod_epu64_2(*x, twice_modulus);
    }

    // Perform approximate computation of Q, as described in page 7 of
    // https://arxiv.org/pdf/2003.04510.pdf
    let q = _mm512_hexl_mulhi_approx_epi_64(w_precon, *y);
    let w_y = _mm512_hexl_mullo_epi_64(w, *y);
    // Compute T in range [0, 4q)
    let t = _mm512_hexl_mullo_add_lo_epi_64(w_y, q, neg_modulus);
    // Reduce T to range [0, 2q)
    let t = _mm512_hexl_small_mod_epu64_2(t, twice_modulus);

    let twice_mod_minus_t = _mm512_sub_epi64(twice_modulus, t);
    *y = _mm512_add_epi64(*x, twice_mod_minus_t);
    *x = _mm512_add_epi64(*x, t);
}

pub fn fwd_t1<const BIT_SHIFT: u32>(
    operand: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    _m: usize,
    w: &[u64],
    w_precon: &[u64],
) {
    let mut v_w_pt: *const __m512i = w.as_ptr().cast();
    let mut v_w_precon_pt: *const __m512i = w_precon.as_ptr().cast();

    for x in unsafe { operand.as_chunks_unchecked_mut::<16>() } {
        unsafe {
            let (mut v_x, mut v_y) = load_fwd_interleaved_t1(x);

            let v_w = _mm512_loadu_si512(v_w_pt);
            v_w_pt = v_w_pt.add(1);
            let v_w_precon = _mm512_loadu_si512(v_w_precon_pt);
            v_w_precon_pt = v_w_precon_pt.add(1);

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
    _m: usize,
    w: &[u64],
    w_precon: &[u64],
) {
    let mut v_w_pt: *const __m512i = w.as_ptr().cast();
    let mut v_w_precon_pt: *const __m512i = w_precon.as_ptr().cast();

    for x in unsafe { operand.as_chunks_unchecked_mut::<16>() } {
        let v_x_pt: *mut __m512i = x.as_mut_ptr().cast();

        unsafe {
            let (mut v_x, mut v_y) = load_fwd_interleaved_t2(x);

            let v_w = _mm512_loadu_si512(v_w_pt);
            v_w_pt = v_w_pt.add(1);
            let v_w_precon = _mm512_loadu_si512(v_w_precon_pt);
            v_w_precon_pt = v_w_precon_pt.add(1);

            fwd_butterfly::<BIT_SHIFT, false>(
                &mut v_x,
                &mut v_y,
                v_w,
                v_w_precon,
                v_neg_modulus,
                v_twice_mod,
            );

            _mm512_storeu_si512(v_x_pt, v_x);
            _mm512_storeu_si512(v_x_pt.add(1), v_y);
        }
    }
}

pub fn fwd_t4<const BIT_SHIFT: u32>(
    operand: &mut [u64],
    v_neg_modulus: __m512i,
    v_twice_mod: __m512i,
    _m: usize,
    w: &[u64],
    w_precon: &[u64],
) {
    let mut v_w_pt: *const __m512i = w.as_ptr().cast();
    let mut v_w_precon_pt: *const __m512i = w_precon.as_ptr().cast();

    for x in unsafe { operand.as_chunks_unchecked_mut::<16>() } {
        let v_x_pt: *mut __m512i = x.as_mut_ptr().cast();

        unsafe {
            let (mut v_x, mut v_y) = load_fwd_interleaved_t4(x);

            let v_w = _mm512_loadu_si512(v_w_pt);
            v_w_pt = v_w_pt.add(1);
            let v_w_precon = _mm512_loadu_si512(v_w_precon_pt);
            v_w_precon_pt = v_w_precon_pt.add(1);

            fwd_butterfly::<BIT_SHIFT, false>(
                &mut v_x,
                &mut v_y,
                v_w,
                v_w_precon,
                v_neg_modulus,
                v_twice_mod,
            );

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
    m: usize,
    w: &[u64],
    w_precon: &[u64],
) {
    let mut j1 = 0;
    let mut w_iter = w.iter().copied();
    let mut w_precon_iter = w_precon.iter().copied();

    for _i in 0..m {
        let mut v_x_op_pt: *mut __m512i = operand[j1..].as_mut_ptr().cast();
        let mut v_y_op_pt: *mut __m512i = operand[j1 + t..].as_mut_ptr().cast();

        unsafe {
            // Weights and weights' preconditions
            let v_w: __m512i = _mm512_set1_epi64(w_iter.next().unwrap() as i64);
            let v_w_precon: __m512i = _mm512_set1_epi64(w_precon_iter.next().unwrap() as i64);

            // assume 8 | t
            for _j in 0..t / 8 {
                let mut v_x = _mm512_loadu_si512(v_x_op_pt);
                let mut v_y = _mm512_loadu_si512(v_y_op_pt);

                fwd_butterfly::<BIT_SHIFT, INPUT_LESS_THAN_MOD>(
                    &mut v_x,
                    &mut v_y,
                    v_w,
                    v_w_precon,
                    v_neg_modulus,
                    v_twice_mod,
                );

                _mm512_storeu_si512(v_x_op_pt, v_x);
                _mm512_storeu_si512(v_y_op_pt, v_y);

                v_x_op_pt = v_x_op_pt.add(1);
                v_y_op_pt = v_y_op_pt.add(1);
            }
        }
        j1 += t << 1;
    }
}

pub unsafe fn forward_transform_to_bit_reverse_avx512<const BIT_SHIFT: u32>(
    operand: &mut [u64],
    n: usize,
    modulus: u64,
    root_of_unity_powers: &[u64],
    precon_root_of_unity_powers: &[u64],
    input_mod_factor: u64,
    output_mod_factor: u64,
    recursion_depth: usize,
    recursion_half: usize,
) {
    debug_assert!(n.is_power_of_two());
    let log_n = n.trailing_zeros();
    check_arguments(log_n, modulus);

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
            let w = &root_of_unity_powers[w_idx..];
            let w_precon = &precon_root_of_unity_powers[w_idx..];

            if input_mod_factor <= 2 && recursion_depth == 0 {
                fwd_t8_inplace::<BIT_SHIFT, true>(
                    operand,
                    v_neg_modulus,
                    v_twice_mod,
                    t,
                    m,
                    w,
                    w_precon,
                );
            } else {
                fwd_t8_inplace::<BIT_SHIFT, false>(
                    operand,
                    v_neg_modulus,
                    v_twice_mod,
                    t,
                    m,
                    w,
                    w_precon,
                );
            }

            t >>= 1;
            m <<= 1;
            w_idx <<= 1;
        }

        while m < (n >> 3) {
            let w = &root_of_unity_powers[w_idx..];
            let w_precon = &precon_root_of_unity_powers[w_idx..];

            fwd_t8_inplace::<BIT_SHIFT, false>(
                operand,
                v_neg_modulus,
                v_twice_mod,
                t,
                m,
                w,
                w_precon,
            );

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
                return idx + (5 * n / 8);
            };

            let mut new_w_idx = compute_new_w_idx(w_idx);
            let mut w = &root_of_unity_powers[new_w_idx..];
            let mut w_precon = &precon_root_of_unity_powers[new_w_idx..];
            fwd_t4::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, m, w, w_precon);

            m <<= 1;
            w_idx <<= 1;
            new_w_idx = compute_new_w_idx(w_idx);
            w = &root_of_unity_powers[new_w_idx..];
            w_precon = &precon_root_of_unity_powers[new_w_idx..];
            fwd_t2::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, m, w, w_precon);

            m <<= 1;
            w_idx <<= 1;
            new_w_idx = compute_new_w_idx(w_idx);
            w = &root_of_unity_powers[new_w_idx..];
            w_precon = &precon_root_of_unity_powers[new_w_idx..];
            fwd_t1::<BIT_SHIFT>(operand, v_neg_modulus, v_twice_mod, m, w, w_precon);
        }

        if output_mod_factor == 1 {
            // n power of two at least 8 => n divisible by 8
            let mut v_x_pt: *mut __m512i = operand.as_mut_ptr().cast();
            for _i in (0..n).step_by(8) {
                let mut v_x = unsafe { _mm512_loadu_si512(v_x_pt) };

                // Reduce from [0, 4q) to [0, q)
                v_x = unsafe { _mm512_hexl_small_mod_epu64_2(v_x, v_twice_mod) };
                v_x = unsafe { _mm512_hexl_small_mod_epu64_2(v_x, v_modulus) };

                unsafe {
                    _mm512_storeu_si512(v_x_pt, v_x);
                    v_x_pt = v_x_pt.add(1);
                }
            }
        }
    } else {
        // Perform depth-first NTT via recursive call
        let t = n >> 1;
        let w_idx = (1 << recursion_depth) + recursion_half;
        let w = &root_of_unity_powers[w_idx..];
        let w_precon = &precon_root_of_unity_powers[w_idx..];

        fwd_t8_inplace::<BIT_SHIFT, false>(operand, v_neg_modulus, v_twice_mod, t, 1, w, w_precon);

        let (left, right) = operand.split_at_mut(n / 2);

        unsafe {
            forward_transform_to_bit_reverse_avx512::<BIT_SHIFT>(
                left,
                n / 2,
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
                n / 2,
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
