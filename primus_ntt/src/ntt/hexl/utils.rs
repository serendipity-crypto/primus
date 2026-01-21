use core::arch::x86_64::*;

/// Returns lower NumBits bits from a 64-bit value
#[target_feature(enable = "avx512f")]
pub fn clear_top_bits_64<const NUM_BITS: u32>(x: __m512i) -> __m512i {
    let low52b_mask = _mm512_set1_epi64(((1u64 << NUM_BITS) - 1) as i64);
    _mm512_and_epi64(x, low52b_mask)
}

#[target_feature(enable = "avx512f")]
pub fn _mm512_hexl_mulhi_epi_64(x: __m512i, y: __m512i) -> __m512i {
    // https://stackoverflow.com/questions/28807341/simd-signed-with-unsigned-multiplication-for-64-bit-64-bit-to-128-bit
    let lo_mask = _mm512_set1_epi64(0x00000000ffffffff);
    // Shuffle high bits with low bits in each 64-bit integer =>
    // x0_lo, x0_hi, x1_lo, x1_hi, x2_lo, x2_hi, ...
    let x_hi = _mm512_shuffle_epi32::<0xB1>(x);
    // y0_lo, y0_hi, y1_lo, y1_hi, y2_lo, y2_hi, ...
    let y_hi = _mm512_shuffle_epi32::<0xB1>(y);
    let z_lo_lo = _mm512_mul_epu32(x, y); // x_lo * y_lo
    let z_lo_hi = _mm512_mul_epu32(x, y_hi); // x_lo * y_hi
    let z_hi_lo = _mm512_mul_epu32(x_hi, y); // x_hi * y_lo
    let z_hi_hi = _mm512_mul_epu32(x_hi, y_hi); // x_hi * y_hi

    //                   x_hi | x_lo
    // x                 y_hi | y_lo
    // ------------------------------
    //                  [x_lo * y_lo]    // z_lo_lo
    // +           [z_lo * y_hi]         // z_lo_hi
    // +           [x_hi * y_lo]         // z_hi_lo
    // +    [x_hi * y_hi]                // z_hi_hi
    //     ^-----------^ <-- only bits needed
    //  sum_|  hi | mid | lo  |

    // Low bits of z_lo_lo are not needed
    let z_lo_lo_shift = _mm512_srli_epi64::<32>(z_lo_lo);

    //                   [x_lo  *  y_lo] // z_lo_lo
    //          + [z_lo  *  y_hi]        // z_lo_hi
    //          ------------------------
    //            |    sum_tmp   |
    //            |sum_mid|sum_lo|
    let sum_tmp = _mm512_add_epi64(z_lo_hi, z_lo_lo_shift);
    let sum_lo = _mm512_and_si512(sum_tmp, lo_mask);
    let sum_mid = _mm512_srli_epi64::<32>(sum_tmp);
    //            |       |sum_lo|
    //          + [x_hi   *  y_lo]       // z_hi_lo
    //          ------------------
    //            [   sum_mid2   ]
    let sum_mid2 = _mm512_add_epi64(z_hi_lo, sum_lo);
    let sum_mid2_hi = _mm512_srli_epi64(sum_mid2, 32);
    let sum_hi = _mm512_add_epi64(z_hi_hi, sum_mid);
    _mm512_add_epi64(sum_hi, sum_mid2_hi)
}

#[target_feature(enable = "avx512ifma")]
pub fn _mm512_hexl_mulhi_epi_52(x: __m512i, y: __m512i) -> __m512i {
    let zero = _mm512_set1_epi64(0);
    _mm512_madd52hi_epu64(zero, x, y)
}

#[target_feature(enable = "avx512f")]
pub fn _mm512_hexl_mulhi_approx_epi_64(x: __m512i, y: __m512i) -> __m512i {
    // https://stackoverflow.com/questions/28807341/simd-signed-with-unsigned-multiplication-for-64-bit-64-bit-to-128-bit
    let lo_mask = _mm512_set1_epi64(0x00000000ffffffff);
    // Shuffle high bits with low bits in each 64-bit integer =>
    // x0_lo, x0_hi, x1_lo, x1_hi, x2_lo, x2_hi, ...
    let x_hi = _mm512_shuffle_epi32::<0xB1>(x);
    // y0_lo, y0_hi, y1_lo, y1_hi, y2_lo, y2_hi, ...
    let y_hi = _mm512_shuffle_epi32::<0xB1>(y);
    let z_lo_hi = _mm512_mul_epu32(x, y_hi); // x_lo * y_hi
    let z_hi_lo = _mm512_mul_epu32(x_hi, y); // x_hi * y_lo
    let z_hi_hi = _mm512_mul_epu32(x_hi, y_hi); // x_hi * y_hi

    //                   x_hi | x_lo
    // x                 y_hi | y_lo
    // ------------------------------
    //                  [x_lo * y_lo]    // unused, resulting in approximation
    // +           [z_lo * y_hi]         // z_lo_hi
    // +           [x_hi * y_lo]         // z_hi_lo
    // +    [x_hi * y_hi]                // z_hi_hi
    //     ^-----------^ <-- only bits needed
    //  sum_|  hi | mid | lo  |

    let sum_lo = _mm512_and_si512(z_lo_hi, lo_mask);
    let sum_mid = _mm512_srli_epi64(z_lo_hi, 32);
    //            |       |sum_lo|
    //          + [x_hi   *  y_lo]       // z_hi_lo
    //          ------------------
    //            [   sum_mid2   ]
    let sum_mid2 = _mm512_add_epi64(z_hi_lo, sum_lo);
    let sum_mid2_hi = _mm512_srli_epi64(sum_mid2, 32);
    let sum_hi = _mm512_add_epi64(z_hi_hi, sum_mid);
    _mm512_add_epi64(sum_hi, sum_mid2_hi)
}

#[target_feature(enable = "avx512ifma")]
pub fn _mm512_hexl_mulhi_approx_epi_52(x: __m512i, y: __m512i) -> __m512i {
    let zero = _mm512_set1_epi64(0);
    _mm512_madd52hi_epu64(zero, x, y)
}

// pub fn _mm512_hexl_mullo_epi<const BitShift: u32>(x: __m512i, y: __m512i) -> __m512i {
//     if BitShift == 32 {
//         unimplemented!()
//     } else if BitShift == 52 {
//         _mm512_hexl_mullo_epi_52(x, y)
//     } else if BitShift == 64 {
//         _mm512_hexl_mullo_epi_64(x, y)
//     } else {
//         panic!("Unsupported BitShift value");
//     }
// }

#[target_feature(enable = "avx512dq")]
pub fn _mm512_hexl_mullo_epi_64(x: __m512i, y: __m512i) -> __m512i {
    _mm512_mullo_epi64(x, y)
}

#[target_feature(enable = "avx512ifma")]
pub fn _mm512_hexl_mullo_epi_52(x: __m512i, y: __m512i) -> __m512i {
    let zero = _mm512_set1_epi64(0);
    _mm512_madd52lo_epu64(zero, x, y)
}

#[target_feature(enable = "avx512ifma")]
pub fn _mm512_hexl_mullo_add_lo_epi_52(x: __m512i, y: __m512i, z: __m512i) -> __m512i {
    let result = _mm512_madd52lo_epu64(x, y, z);

    // Clear high 12 bits from result
    clear_top_bits_64::<52>(result)
}

#[target_feature(enable = "avx512dq")]
pub fn _mm512_hexl_mullo_add_lo_epi_64(x: __m512i, y: __m512i, z: __m512i) -> __m512i {
    let prod = _mm512_mullo_epi64(y, z);
    _mm512_add_epi64(x, prod)
}

#[target_feature(enable = "avx512f")]
pub fn _mm512_hexl_small_mod_epu64_2(x: __m512i, q: __m512i) -> __m512i {
    _mm512_min_epu64(x, _mm512_sub_epi64(x, q))
}

// #[target_feature(enable = "avx512f")]
// pub fn _mm512_hexl_small_mod_epu64<const INPUT_MOD_FACTOR: u8>(
//     mut x: __m512i,
//     q: __m512i,
//     q_times_2: Option<__m512i>,
//     q_times_4: Option<__m512i>,
// ) -> __m512i {
//     debug_assert!(
//         INPUT_MOD_FACTOR == 1
//             || INPUT_MOD_FACTOR == 2
//             || INPUT_MOD_FACTOR == 4
//             || INPUT_MOD_FACTOR == 8,
//         "InputModFactor must be 1, 2, 4, or 8"
//     );

//     if INPUT_MOD_FACTOR == 1 {
//         return x;
//     }

//     if INPUT_MOD_FACTOR == 2 {
//         return _mm512_min_epu64(x, _mm512_sub_epi64(x, q));
//     }

//     if INPUT_MOD_FACTOR == 4 {
//         debug_assert!(
//             q_times_2.is_some(),
//             "q_times_2 must be provided for INPUT_MOD_FACTOR=4"
//         );
//         let q_times_2 = q_times_2.unwrap();
//         x = _mm512_min_epu64(x, _mm512_sub_epi64(x, q_times_2));
//         return _mm512_min_epu64(x, _mm512_sub_epi64(x, q));
//     }

//     if INPUT_MOD_FACTOR == 8 {
//         debug_assert!(
//             q_times_2.is_some(),
//             "q_times_2 must be provided for INPUT_MOD_FACTOR=4"
//         );
//         debug_assert!(
//             q_times_4.is_some(),
//             "q_times_4 must be provided for INPUT_MOD_FACTOR=8"
//         );
//         let q_times_2 = q_times_2.unwrap();
//         let q_times_4 = q_times_4.unwrap();
//         x = _mm512_min_epu64(x, _mm512_sub_epi64(x, q_times_4));
//         x = _mm512_min_epu64(x, _mm512_sub_epi64(x, q_times_2));
//         return _mm512_min_epu64(x, _mm512_sub_epi64(x, q));
//     }

//     debug_assert!(false, "Invalid InputModFactor");

//     x
// }

// Returns (x + y) mod q; assumes 0 < x, y < q
#[target_feature(enable = "avx512f")]
pub fn _mm512_hexl_small_add_mod_epi64(x: __m512i, y: __m512i, q: __m512i) -> __m512i {
    _mm512_hexl_small_mod_epu64_2(_mm512_add_epi64(x, y), q)
}

// Given input: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
// Returns
// out1 =  _mm512_set_epi64(14, 6, 12, 4, 10, 2, 8, 0);
// out2 =  _mm512_set_epi64(15, 7, 13, 5, 11, 3, 9, 1);
#[target_feature(enable = "avx512f")]
pub fn load_fwd_interleaved_t1(arg: *const u64) -> (__m512i, __m512i) {
    let arg_512: *const __m512i = arg.cast();

    // 0, 1, 2, 3, 4, 5, 6, 7
    let v1 = unsafe { _mm512_loadu_si512(arg_512) };
    // 8, 9, 10, 11, 12, 13, 14, 15
    let v2 = unsafe { _mm512_loadu_si512(arg_512.add(1)) };

    let perm_idx = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1);

    // 1, 0, 3, 2, 5, 4, 7, 6
    let v1_perm = _mm512_permutexvar_epi64(perm_idx, v1);
    // 9, 8, 11, 10, 13, 12, 15, 14
    let v2_perm = _mm512_permutexvar_epi64(perm_idx, v2);

    (
        _mm512_mask_blend_epi64(0xaa, v1, v2_perm),
        _mm512_mask_blend_epi64(0xaa, v1_perm, v2),
    )
}

// Given input: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
// Returns
// out1 =  _mm512_set_epi64(14, 12, 10, 8, 6, 4, 2, 0);
// out2 =  _mm512_set_epi64(15, 13, 11, 9, 7, 5, 3, 1);
#[target_feature(enable = "avx512f")]
pub fn load_inv_interleaved_t1(arg: &[u64; 16]) -> (__m512i, __m512i) {
    let vperm_hi_idx = _mm512_set_epi64(6, 4, 2, 0, 7, 5, 3, 1);
    let vperm_lo_idx = _mm512_set_epi64(7, 5, 3, 1, 6, 4, 2, 0);
    let vperm2_idx = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);

    let arg_512: *const __m512i = arg.as_ptr().cast();

    // 7, 6, 5, 4, 3, 2, 1, 0
    let v_7to0 = unsafe { _mm512_loadu_si512(arg_512) };
    // 15, 14, 13, 12, 11, 10, 9, 8
    let v_15to8 = unsafe { _mm512_loadu_si512(arg_512.add(1)) };
    // 7, 5, 3, 1, 6, 4, 2, 0
    let perm_lo = _mm512_permutexvar_epi64(vperm_lo_idx, v_7to0);
    // 14, 12, 10, 8, 15, 13, 11, 9
    let perm_hi = _mm512_permutexvar_epi64(vperm_hi_idx, v_15to8);

    let out1 = _mm512_mask_blend_epi64(0x0f, perm_hi, perm_lo);
    let out2 = _mm512_mask_blend_epi64(0xf0, perm_hi, perm_lo);
    (out1, _mm512_permutexvar_epi64(vperm2_idx, out2))
}

// Given input: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
// Returns
// out1 =  _mm512_set_epi64(13, 12, 9, 8, 5, 4, 1, 0);
// out2 =  _mm512_set_epi64(15, 14, 11, 10, 7, 6, 3, 2)
#[target_feature(enable = "avx512f")]
pub fn load_fwd_interleaved_t2(arg: *const u64) -> (__m512i, __m512i) {
    let arg_512: *const __m512i = arg.cast();

    // 11, 10, 9, 8, 3, 2, 1, 0
    let v1 = unsafe { _mm512_loadu_si512(arg_512) };
    // 15, 14, 13, 12, 7, 6, 5, 4
    let v2 = unsafe { _mm512_loadu_si512(arg_512.add(1)) };

    let v1_perm_idx = _mm512_set_epi64(5, 4, 7, 6, 1, 0, 3, 2);

    let v1_perm = _mm512_permutexvar_epi64(v1_perm_idx, v1);
    let v2_perm = _mm512_permutexvar_epi64(v1_perm_idx, v2);

    (
        _mm512_mask_blend_epi64(0xcc, v1, v2_perm),
        _mm512_mask_blend_epi64(0xcc, v1_perm, v2),
    )
}

// Given input: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
// Returns
// out1 =  _mm512_set_epi64(14, 12, 10, 8, 6, 4, 2, 0);
// out2 =  _mm512_set_epi64(15, 13, 11, 9, 7, 5, 3, 1);
#[target_feature(enable = "avx512f")]
pub fn load_inv_interleaved_t2(arg: &[u64; 16]) -> (__m512i, __m512i) {
    let arg_512: *const __m512i = arg.as_ptr().cast();

    let v1 = unsafe { _mm512_loadu_si512(arg_512) };
    let v2 = unsafe { _mm512_loadu_si512(arg_512.add(1)) };

    let v1_perm_idx = _mm512_set_epi64(6, 7, 4, 5, 2, 3, 0, 1);

    let v1_perm = _mm512_permutexvar_epi64(v1_perm_idx, v1);
    let v2_perm = _mm512_permutexvar_epi64(v1_perm_idx, v2);

    (
        _mm512_mask_blend_epi64(0xaa, v1, v2_perm),
        _mm512_mask_blend_epi64(0xaa, v1_perm, v2),
    )
}

// Returns
// out1 =  _mm512_set_epi64(arg[11], arg[10], arg[9], arg[8],
//                           arg[3], arg[2], arg[1], arg[0]);
// out2 =  _mm512_set_epi64(arg[15], arg[14], arg[13], arg[12],
//                           arg[7], arg[6], arg[5], arg[4]);
#[target_feature(enable = "avx512f")]
pub fn load_fwd_interleaved_t4(arg: *const u64) -> (__m512i, __m512i) {
    let arg_512: *const __m512i = arg.cast();

    let vperm2_idx = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    let v_7to0 = unsafe { _mm512_loadu_si512(arg_512) };
    let v_15to8 = unsafe { _mm512_loadu_si512(arg_512.add(1)) };
    let perm_hi = _mm512_permutexvar_epi64(vperm2_idx, v_15to8);
    let out1 = _mm512_mask_blend_epi64(0x0f, perm_hi, v_7to0);
    let mut out2 = _mm512_mask_blend_epi64(0xf0, perm_hi, v_7to0);
    out2 = _mm512_permutexvar_epi64(vperm2_idx, out2);
    (out1, out2)
}

#[target_feature(enable = "avx512f")]
pub fn load_inv_interleaved_t4(arg: &[u64; 16]) -> (__m512i, __m512i) {
    let arg_512: *const __m512i = arg.as_ptr().cast();

    // 0, 1, 2, 3, 4, 5, 6, 7
    let v1 = unsafe { _mm512_loadu_si512(arg_512) };
    // 8, 9, 10, 11, 12, 13, 14, 15
    let v2 = unsafe { _mm512_loadu_si512(arg_512.add(1)) };
    let perm_idx = _mm512_set_epi64(5, 4, 7, 6, 1, 0, 3, 2);

    // 1, 0, 3, 2, 5, 4, 7, 6
    let v1_perm = _mm512_permutexvar_epi64(perm_idx, v1);
    // 9, 8, 11, 10, 13, 12, 15, 14
    let v2_perm = _mm512_permutexvar_epi64(perm_idx, v2);

    (
        _mm512_mask_blend_epi64(0xcc, v1, v2_perm),
        _mm512_mask_blend_epi64(0xcc, v1_perm, v2),
    )
}

// Given inputs
//     arg1 = _mm512_set_epi64(15, 14, 13, 12, 11, 10, 9, 8);
//     arg2 = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
// Writes out = {8,  0, 9,  1, 10, 2, 11, 3,
//               12, 4, 13, 5, 14, 6, 15, 7}
#[target_feature(enable = "avx512f")]
pub fn write_fwd_interleaved_t1(mut arg1: __m512i, mut arg2: __m512i, out: *mut __m512i) {
    let vperm2_idx = _mm512_set_epi64(3, 2, 1, 0, 7, 6, 5, 4);
    let v_x_out_idx = _mm512_set_epi64(7, 3, 6, 2, 5, 1, 4, 0);
    let v_y_out_idx = _mm512_set_epi64(3, 7, 2, 6, 1, 5, 0, 4);

    // v_Y => (4, 5, 6, 7, 0, 1, 2, 3)
    arg2 = _mm512_permutexvar_epi64(vperm2_idx, arg2);
    // 4, 5, 6, 7, 12, 13, 14, 15
    let perm_lo = _mm512_mask_blend_epi64(0x0f, arg1, arg2);

    // 8, 9, 10, 11, 0, 1, 2, 3
    let perm_hi = _mm512_mask_blend_epi64(0xf0, arg1, arg2);

    arg1 = _mm512_permutexvar_epi64(v_x_out_idx, perm_hi);
    arg2 = _mm512_permutexvar_epi64(v_y_out_idx, perm_lo);

    unsafe {
        _mm512_storeu_si512(out, arg1);
        _mm512_storeu_si512(out.add(1), arg2);
    }
}

// Given inputs
// @param arg1 = _mm512_set_epi64(15, 14, 13, 12, 11, 10, 9, 8);
// @param arg2 = _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0);
// Writes out = {8,  9,  10, 11, 0, 1, 2, 3,
//               12, 13, 14, 15, 4, 5, 6, 7}
#[target_feature(enable = "avx512f,avx")]
pub fn write_inv_interleaved_t4(arg1: __m512i, arg2: __m512i, out: *mut __m512i) {
    let x0 = _mm512_extracti64x4_epi64(arg1, 0);
    let x1 = _mm512_extracti64x4_epi64(arg1, 1);
    let y0 = _mm512_extracti64x4_epi64(arg2, 0);
    let y1 = _mm512_extracti64x4_epi64(arg2, 1);
    let out_256: *mut __m256i = out.cast();

    unsafe {
        _mm256_storeu_si256(out_256, x0);
        _mm256_storeu_si256(out_256.add(1), y0);
        _mm256_storeu_si256(out_256.add(2), x1);
        _mm256_storeu_si256(out_256.add(3), y1);
    }
}

// Returns _mm512_set_epi64(arg[3], arg[3], arg[2], arg[2],
//                          arg[1], arg[1], arg[0], arg[0]);
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx")]
pub fn load_w_op_t2(arg: &[u64; 4]) -> __m512i {
    let vperm_w_idx: __m512i = _mm512_set_epi64(3, 3, 2, 2, 1, 1, 0, 0);
    let v_w_256: __m256i = unsafe { _mm256_loadu_si256(arg.as_ptr().cast()) };
    let v_w: __m512i = _mm512_broadcast_i64x4(v_w_256);
    _mm512_permutexvar_epi64(vperm_w_idx, v_w)
}

// Returns _mm512_set_epi64(arg[1], arg[1], arg[1], arg[1],
//                          arg[0], arg[0], arg[0], arg[0]);
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512dq,avx512f,sse2")]
pub fn load_w_op_t4(arg: &[u64; 2]) -> __m512i {
    let vperm_w_idx: __m512i = _mm512_set_epi64(1, 1, 1, 1, 0, 0, 0, 0);
    let v_w_128: __m128i = unsafe { _mm_loadu_si128(arg.as_ptr().cast()) };
    let v_w: __m512i = _mm512_broadcast_i64x2(v_w_128);
    _mm512_permutexvar_epi64(vperm_w_idx, v_w)
}
