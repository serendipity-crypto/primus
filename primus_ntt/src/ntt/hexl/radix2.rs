use crate::ntt::hexl::{
    internal::check_arguments,
    number_theory::{multiply_mod_lazy, reduce_once, reduce_twice},
};

// pub fn fwd_butterfly_radix2(
//     x_r: &mut u64,
//     y_r: &mut u64,
//     x_op: u64,
//     y_op: u64,
//     w: u64,
//     w_precon: u64,
//     modulus: u64,
//     twice_modulus: u64,
// ) {
//     let tx = reduce_mod_2(x_op, twice_modulus);
//     let t = multiply_mod_lazy::<64>(y_op, w, w_precon, modulus);

//     *x_r = tx + t;
//     *y_r = tx + twice_modulus - t;
// }

/// Out-of-place Harvey butterfly: assume `X`, `Y` in `[0, 4q)`, and return
/// `X`, `Y` in `[0, 4q)` such that
/// `X = X + W*Y (mod q)`, `Y = X - W*Y (mod q)`.
///
/// # Parameters
/// - `X`: Input butterfly data.
/// - `Y`: Input butterfly data.
/// - `W`: Root of unity.
/// - `W_precon`: Preconditioned `W` for `BitShift`-bit Barrett reduction.
/// - `modulus`: Negative modulus, i.e. `(-q)` represented as 8 64-bit signed integers in SIMD form.
/// - `twice_modulus`: Twice the modulus, i.e. `2*q` represented as 8 64-bit signed integers in SIMD form.
///
/// # Details
/// See Algorithm 4 of https://arxiv.org/pdf/1205.2926.pdf
#[inline]
pub fn fwd_butterfly_radix2(
    x: &mut u64,
    y: &mut u64,
    w: u64,
    w_precon: u64,
    modulus: u64,
    twice_modulus: u64,
) {
    let tx = reduce_once(*x, twice_modulus);
    let t = multiply_mod_lazy::<64>(*y, w, w_precon, modulus);

    *x = tx + t;
    *y = tx + twice_modulus - t;
}

pub fn forward_transform_to_bit_reverse_radix2_inplace(
    operand: &mut [u64],
    modulus: u64,
    root_of_unity_powers: &[u64],
    precon_root_of_unity_powers: &[u64],
    output_mod_factor: u32,
) {
    let n = operand.len();

    check_arguments(n, modulus);

    debug_assert!(
        output_mod_factor == 1 || output_mod_factor == 4,
        "output_mod_factor must be 1 or 4; got {output_mod_factor}"
    );

    let twice_modulus = modulus << 1;

    let mut w_iter = root_of_unity_powers.iter().copied();
    let mut w_precon_iter = precon_root_of_unity_powers.iter().copied();

    w_iter.next(); // skip w[0]
    w_precon_iter.next(); // skip w_precon[0]

    let mut t = n >> 1;
    let mut m = 1;
    while m < n {
        match t {
            8 => {
                let chunks = unsafe { operand.as_chunks_unchecked_mut::<16>() };
                for chunk in chunks {
                    let w = w_iter.next().unwrap();
                    let w_precon = w_precon_iter.next().unwrap();

                    let [
                        x_r0,
                        x_r1,
                        x_r2,
                        x_r3,
                        x_r4,
                        x_r5,
                        x_r6,
                        x_r7,
                        y_r0,
                        y_r1,
                        y_r2,
                        y_r3,
                        y_r4,
                        y_r5,
                        y_r6,
                        y_r7,
                    ] = chunk;

                    fwd_butterfly_radix2(x_r0, y_r0, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r1, y_r1, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r2, y_r2, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r3, y_r3, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r4, y_r4, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r5, y_r5, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r6, y_r6, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r7, y_r7, w, w_precon, modulus, twice_modulus);
                }
            }
            4 => {
                let chunks = unsafe { operand.as_chunks_unchecked_mut::<8>() };
                for chunk in chunks {
                    let w = w_iter.next().unwrap();
                    let w_precon = w_precon_iter.next().unwrap();

                    let [x_r0, x_r1, x_r2, x_r3, y_r0, y_r1, y_r2, y_r3] = chunk;

                    fwd_butterfly_radix2(x_r0, y_r0, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r1, y_r1, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r2, y_r2, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r3, y_r3, w, w_precon, modulus, twice_modulus);
                }
            }
            2 => {
                let chunks = unsafe { operand.as_chunks_unchecked_mut::<4>() };
                for chunk in chunks {
                    let w = w_iter.next().unwrap();
                    let w_precon = w_precon_iter.next().unwrap();

                    let [x_r0, x_r1, y_r0, y_r1] = chunk;

                    fwd_butterfly_radix2(x_r0, y_r0, w, w_precon, modulus, twice_modulus);
                    fwd_butterfly_radix2(x_r1, y_r1, w, w_precon, modulus, twice_modulus);
                }
            }
            1 => {
                let chunks = unsafe { operand.as_chunks_unchecked_mut::<2>() };
                for chunk in chunks {
                    let w = w_iter.next().unwrap();
                    let w_precon = w_precon_iter.next().unwrap();

                    let [x_r, y_r] = chunk;

                    fwd_butterfly_radix2(x_r, y_r, w, w_precon, modulus, twice_modulus);
                }
            }
            _ => {
                for chunk in operand.chunks_exact_mut(t * 2) {
                    let w = w_iter.next().unwrap();
                    let w_precon = w_precon_iter.next().unwrap();

                    let (x, y) = chunk.split_at_mut(t);

                    let x_chunks = unsafe { x.as_chunks_unchecked_mut::<8>() };
                    let y_chunks = unsafe { y.as_chunks_unchecked_mut::<8>() };

                    for (x_chunk, y_chunk) in x_chunks.iter_mut().zip(y_chunks.iter_mut()) {
                        let [x_r0, x_r1, x_r2, x_r3, x_r4, x_r5, x_r6, x_r7] = x_chunk;
                        let [y_r0, y_r1, y_r2, y_r3, y_r4, y_r5, y_r6, y_r7] = y_chunk;

                        fwd_butterfly_radix2(x_r0, y_r0, w, w_precon, modulus, twice_modulus);
                        fwd_butterfly_radix2(x_r1, y_r1, w, w_precon, modulus, twice_modulus);
                        fwd_butterfly_radix2(x_r2, y_r2, w, w_precon, modulus, twice_modulus);
                        fwd_butterfly_radix2(x_r3, y_r3, w, w_precon, modulus, twice_modulus);
                        fwd_butterfly_radix2(x_r4, y_r4, w, w_precon, modulus, twice_modulus);
                        fwd_butterfly_radix2(x_r5, y_r5, w, w_precon, modulus, twice_modulus);
                        fwd_butterfly_radix2(x_r6, y_r6, w, w_precon, modulus, twice_modulus);
                        fwd_butterfly_radix2(x_r7, y_r7, w, w_precon, modulus, twice_modulus);
                    }
                }
            }
        }
        t >>= 1;
        m <<= 1;
    }

    if output_mod_factor == 1 {
        operand.iter_mut().for_each(|x| {
            *x = reduce_twice(*x, modulus, twice_modulus);
        });
    }
}
