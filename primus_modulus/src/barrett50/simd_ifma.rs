//! AVX-512-IFMA fast-path slice kernels for `Barrett50Modulus`.
//!
//! Each kernel processes 8 u64 lanes per iteration (one zmm register).
//! Tails shorter than 8 elements delegate to the inner
//! [`crate::BarrettModulus<u64>`] slice impls, which already SIMD-vectorize
//! 8-wide chunks via portable_simd and handle a per-element scalar tail.
//!
//! # Algorithm (per 8-lane chunk, both operands < m, m ∈ [2^48, 2^50))
//!
//! Notation: `µ_104 = ⌊2^104 / m⌋`, split as `µ_lo52 || µ_hi`.
//!
//! ```text
//! x_lo  = madd52lo(0, a, b)          // (a*b) mod 2^52
//! x_hi  = madd52hi(0, a, b)          // ⌊a*b / 2^52⌋, < 2^48
//! q     = x_hi * µ_hi                // SIMD mullo64, ≤ 2^52
//! q    += madd52hi(q, x_hi, µ_lo52)  // accumulate ⌊x_hi·µ_lo52 / 2^52⌋
//! q    += madd52hi(q, x_lo, µ_hi)    // accumulate ⌊x_lo·µ_hi   / 2^52⌋
//! r     = madd52lo(x_lo, q, 2^52 − m) & (2^52 − 1)
//!       // = (x_lo − q·m) mod 2^52,  r ∈ [0, 4m)
//! ```
//!
//! # Error analysis
//!
//! Let `q_true = ⌊x/m⌋`. Standard Barrett with the wider shift `2k+s` (s=4
//! for m near 2^50, larger for smaller m) tightens the error to
//! `q_Barrett ∈ {q_true−1, q_true}`. The IFMA approximation drops the
//! sub-1 fractional term `x_lo·µ_lo52/2^104` and floors two terms
//! independently, contributing at most 2 additional loss. So
//! `q_ifma ∈ [q_true − 3, q_true]` and
//! `r = x − q_ifma·m ∈ [0, 4m)`.
//!
//! Canonical reduction needs **3** `min(r, r−m)` passes; lazy form (`[0, 2m)`)
//! needs **2**. For `mul_add` family the carry `c < m` is added before
//! reduction, growing the upper bound by `m`: **4** passes canonical,
//! **3** passes lazy.

use core::arch::x86_64::*;

use primus_reduce::{
    LazyReduceMulAddSlice, LazyReduceMulSlice, ReduceAdd, ReduceDotProduct, ReduceMulAddSlice,
    ReduceMulSlice,
};

use super::Barrett50Modulus;

/// Lane count = AVX-512 zmm width / 64 bits.
const N: usize = 8;

/// Inner-chunk size for dot_product double-word accumulation.
///
/// `K · m² < 2^104` with m < 2^50 gives the safety margin: 16 widening
/// products accumulate into a `(lo52, hi52)` double-word without overflow.
const DOT_INNER_CHUNK: usize = 16;

// ===========================================================================
// Per-chunk primitives
// ===========================================================================

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
#[inline]
fn ifma_widening_mul(a: __m512i, b: __m512i) -> (__m512i, __m512i) {
    let zero = _mm512_setzero_si512();
    let lo = _mm512_madd52lo_epu64(zero, a, b);
    let hi = _mm512_madd52hi_epu64(zero, a, b);
    (lo, hi)
}

/// Compute `r = (x_lo + x_hi · 2^52) mod m`, lazily reduced to `[0, 4m)`.
///
/// Caller must guarantee `x_lo < 2^52` and `x_hi < 2^48` lane-wise (true
/// for any `a · b` with `a, b < 2^50`).
#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
#[inline]
fn ifma_lazy_reduce_4m(
    x_lo: __m512i,
    x_hi: __m512i,
    mu_lo52: __m512i,
    mu_hi: __m512i,
    neg_m: __m512i,
    pow52_mask: __m512i,
) -> __m512i {
    // q ≈ ⌊x · µ_104 / 2^104⌋
    //   = x_hi·µ_hi + ⌊x_hi·µ_lo52 / 2^52⌋ + ⌊x_lo·µ_hi / 2^52⌋ + ε
    // with ε ∈ {0, 1, 2}.
    let mut q = _mm512_mullo_epi64(x_hi, mu_hi);
    q = _mm512_madd52hi_epu64(q, x_hi, mu_lo52);
    q = _mm512_madd52hi_epu64(q, x_lo, mu_hi);

    // r = (x_lo − q·m) mod 2^52. Using madd52lo with `2^52 − m`:
    //   x_lo + low52(q · (2^52 − m))
    //     = x_lo + low52(q·2^52 − q·m)
    //     = x_lo + low52(−q·m)
    // which equals `(x_lo − q·m) mod 2^52` once we mask the carry bit.
    let raw = _mm512_madd52lo_epu64(x_lo, q, neg_m);
    _mm512_and_si512(raw, pow52_mask)
}

/// One `min(r, r − m)` round.
#[target_feature(enable = "avx512f")]
#[inline]
fn cond_sub_m(r: __m512i, m: __m512i) -> __m512i {
    _mm512_min_epu64(r, _mm512_sub_epi64(r, m))
}

/// Multiply two lanes and reduce to canonical `[0, m)`.
#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
#[inline]
fn ifma_mul_canonical(
    a: __m512i,
    b: __m512i,
    mu_lo52: __m512i,
    mu_hi: __m512i,
    neg_m: __m512i,
    pow52_mask: __m512i,
    m: __m512i,
) -> __m512i {
    let (x_lo, x_hi) = ifma_widening_mul(a, b);
    let r = ifma_lazy_reduce_4m(x_lo, x_hi, mu_lo52, mu_hi, neg_m, pow52_mask);
    // r ∈ [0, 4m) → three rounds of conditional subtract bring it to [0, m).
    let r = cond_sub_m(r, m);
    let r = cond_sub_m(r, m);
    cond_sub_m(r, m)
}

/// Multiply two lanes and reduce to lazy `[0, 2m)`.
#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
#[inline]
fn ifma_mul_lazy(
    a: __m512i,
    b: __m512i,
    mu_lo52: __m512i,
    mu_hi: __m512i,
    neg_m: __m512i,
    pow52_mask: __m512i,
    m: __m512i,
) -> __m512i {
    let (x_lo, x_hi) = ifma_widening_mul(a, b);
    let r = ifma_lazy_reduce_4m(x_lo, x_hi, mu_lo52, mu_hi, neg_m, pow52_mask);
    let r = cond_sub_m(r, m);
    cond_sub_m(r, m)
}

/// Returns 6 zmm constants used by every kernel.
#[target_feature(enable = "avx512f")]
#[inline]
fn splat_params(modulus: &Barrett50Modulus) -> (__m512i, __m512i, __m512i, __m512i, __m512i) {
    let m = _mm512_set1_epi64(modulus.value as i64);
    let mu_lo52 = _mm512_set1_epi64(modulus.mu_lo52 as i64);
    let mu_hi = _mm512_set1_epi64(modulus.mu_hi as i64);
    let neg_m = _mm512_set1_epi64(modulus.neg_m_mod_pow2_52 as i64);
    let pow52_mask = _mm512_set1_epi64(((1u64 << 52) - 1) as i64);
    (m, mu_lo52, mu_hi, neg_m, pow52_mask)
}

// ===========================================================================
// Load / store helpers
// ===========================================================================

#[inline]
unsafe fn load_u64x8(slice: &[u64]) -> __m512i {
    debug_assert!(slice.len() >= N);
    unsafe { _mm512_loadu_si512(slice.as_ptr().cast()) }
}

#[inline]
unsafe fn store_u64x8(slice: &mut [u64], v: __m512i) {
    debug_assert!(slice.len() >= N);
    unsafe { _mm512_storeu_si512(slice.as_mut_ptr().cast(), v) }
}

// ===========================================================================
// reduce_mul_slice_{assign,to}
// ===========================================================================

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn reduce_mul_slice_to(
    modulus: Barrett50Modulus,
    a: &[u64],
    b: &[u64],
    output: &mut [u64],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let chunks = a.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        let r = ifma_mul_canonical(av, bv, mu_lo52, mu_hi, neg_m, pow52_mask, m);
        unsafe { store_u64x8(&mut output[off..], r) };
    }
    modulus
        .inner
        .reduce_mul_slice_to(&a[main_len..], &b[main_len..], &mut output[main_len..]);
}

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn reduce_mul_slice_assign(modulus: Barrett50Modulus, a: &mut [u64], b: &[u64]) {
    debug_assert_eq!(a.len(), b.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let chunks = a.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        let r = ifma_mul_canonical(av, bv, mu_lo52, mu_hi, neg_m, pow52_mask, m);
        unsafe { store_u64x8(&mut a[off..], r) };
    }
    let (head, tail) = a.split_at_mut(main_len);
    let _ = head; // already written
    modulus.inner.reduce_mul_slice_assign(tail, &b[main_len..]);
}

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn lazy_reduce_mul_slice_to(
    modulus: Barrett50Modulus,
    a: &[u64],
    b: &[u64],
    output: &mut [u64],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), output.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let chunks = a.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        let r = ifma_mul_lazy(av, bv, mu_lo52, mu_hi, neg_m, pow52_mask, m);
        unsafe { store_u64x8(&mut output[off..], r) };
    }
    modulus
        .inner
        .lazy_reduce_mul_slice_to(&a[main_len..], &b[main_len..], &mut output[main_len..]);
}

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn lazy_reduce_mul_slice_assign(modulus: Barrett50Modulus, a: &mut [u64], b: &[u64]) {
    debug_assert_eq!(a.len(), b.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let chunks = a.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        let r = ifma_mul_lazy(av, bv, mu_lo52, mu_hi, neg_m, pow52_mask, m);
        unsafe { store_u64x8(&mut a[off..], r) };
    }
    let (_head, tail) = a.split_at_mut(main_len);
    modulus
        .inner
        .lazy_reduce_mul_slice_assign(tail, &b[main_len..]);
}

// ===========================================================================
// reduce_mul_add_slice_{to,assign-style}
// ===========================================================================

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn reduce_mul_add_slice_to(
    modulus: Barrett50Modulus,
    a: &[u64],
    b: &[u64],
    c: &[u64],
    output: &mut [u64],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), c.len());
    debug_assert_eq!(a.len(), output.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let chunks = a.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        let cv = unsafe { load_u64x8(&c[off..]) };
        let (x_lo, x_hi) = ifma_widening_mul(av, bv);
        let r = ifma_lazy_reduce_4m(x_lo, x_hi, mu_lo52, mu_hi, neg_m, pow52_mask);
        // r ∈ [0, 4m), c ∈ [0, m), sum ∈ [0, 5m). Four cond-subs to canonical.
        let r = _mm512_add_epi64(r, cv);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        unsafe { store_u64x8(&mut output[off..], r) };
    }
    modulus.inner.reduce_mul_add_slice_to(
        &a[main_len..],
        &b[main_len..],
        &c[main_len..],
        &mut output[main_len..],
    );
}

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn lazy_reduce_mul_add_slice_to(
    modulus: Barrett50Modulus,
    a: &[u64],
    b: &[u64],
    c: &[u64],
    output: &mut [u64],
) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), c.len());
    debug_assert_eq!(a.len(), output.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let chunks = a.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        let cv = unsafe { load_u64x8(&c[off..]) };
        let (x_lo, x_hi) = ifma_widening_mul(av, bv);
        let r = ifma_lazy_reduce_4m(x_lo, x_hi, mu_lo52, mu_hi, neg_m, pow52_mask);
        let r = _mm512_add_epi64(r, cv);
        // [0, 5m) → three cond-subs → [0, 2m).
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        unsafe { store_u64x8(&mut output[off..], r) };
    }
    modulus.inner.lazy_reduce_mul_add_slice_to(
        &a[main_len..],
        &b[main_len..],
        &c[main_len..],
        &mut output[main_len..],
    );
}

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn reduce_add_mul_slice_assign(
    modulus: Barrett50Modulus,
    acc: &mut [u64],
    a: &[u64],
    b: &[u64],
) {
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let chunks = acc.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let accv = unsafe { load_u64x8(&acc[off..]) };
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        let (x_lo, x_hi) = ifma_widening_mul(av, bv);
        let r = ifma_lazy_reduce_4m(x_lo, x_hi, mu_lo52, mu_hi, neg_m, pow52_mask);
        let r = _mm512_add_epi64(r, accv);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        unsafe { store_u64x8(&mut acc[off..], r) };
    }
    let (_head, tail) = acc.split_at_mut(main_len);
    modulus
        .inner
        .reduce_add_mul_slice_assign(tail, &a[main_len..], &b[main_len..]);
}

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn lazy_reduce_add_mul_slice_assign(
    modulus: Barrett50Modulus,
    acc: &mut [u64],
    a: &[u64],
    b: &[u64],
) {
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let chunks = acc.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let accv = unsafe { load_u64x8(&acc[off..]) };
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        let (x_lo, x_hi) = ifma_widening_mul(av, bv);
        let r = ifma_lazy_reduce_4m(x_lo, x_hi, mu_lo52, mu_hi, neg_m, pow52_mask);
        let r = _mm512_add_epi64(r, accv);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        unsafe { store_u64x8(&mut acc[off..], r) };
    }
    let (_head, tail) = acc.split_at_mut(main_len);
    modulus
        .inner
        .lazy_reduce_add_mul_slice_assign(tail, &a[main_len..], &b[main_len..]);
}

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn reduce_sub_mul_slice_assign(
    modulus: Barrett50Modulus,
    acc: &mut [u64],
    a: &[u64],
    b: &[u64],
) {
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let chunks = acc.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let accv = unsafe { load_u64x8(&acc[off..]) };
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        // Fully reduce prod to canonical [0, m) so the subsequent
        // sub_mod has well-defined semantics.
        let prod = ifma_mul_canonical(av, bv, mu_lo52, mu_hi, neg_m, pow52_mask, m);
        // (acc − prod) mod m using one cond-add of m on underflow.
        let diff = _mm512_sub_epi64(accv, prod);
        let lt_mask = _mm512_cmplt_epu64_mask(accv, prod);
        let r = _mm512_mask_add_epi64(diff, lt_mask, diff, m);
        unsafe { store_u64x8(&mut acc[off..], r) };
    }
    let (_head, tail) = acc.split_at_mut(main_len);
    modulus
        .inner
        .reduce_sub_mul_slice_assign(tail, &a[main_len..], &b[main_len..]);
}

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn lazy_reduce_sub_mul_slice_assign(
    modulus: Barrett50Modulus,
    acc: &mut [u64],
    a: &[u64],
    b: &[u64],
) {
    // Lazy version: skip reducing prod past `[0, 2m)`, then use the
    // standard 2m-correction pattern matching barrett/simd.rs.
    debug_assert_eq!(acc.len(), a.len());
    debug_assert_eq!(acc.len(), b.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);
    let two_m = _mm512_add_epi64(m, m);

    let chunks = acc.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let accv = unsafe { load_u64x8(&acc[off..]) };
        let av = unsafe { load_u64x8(&a[off..]) };
        let bv = unsafe { load_u64x8(&b[off..]) };
        let prod = ifma_mul_lazy(av, bv, mu_lo52, mu_hi, neg_m, pow52_mask, m);
        // prod ∈ [0, 2m), acc ∈ [0, m). diff_true = acc − prod ∈ (−2m, m).
        // Adding 2m on underflow gives result in [0, 2m) ⊂ [0, 2m) — lazy.
        let diff = _mm512_sub_epi64(accv, prod);
        let lt_mask = _mm512_cmplt_epu64_mask(accv, prod);
        let r = _mm512_mask_add_epi64(diff, lt_mask, diff, two_m);
        unsafe { store_u64x8(&mut acc[off..], r) };
    }
    let (_head, tail) = acc.split_at_mut(main_len);
    modulus
        .inner
        .lazy_reduce_sub_mul_slice_assign(tail, &a[main_len..], &b[main_len..]);
}

// ===========================================================================
// reduce_scalar_mul_add_slice_to
// ===========================================================================

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn reduce_scalar_mul_add_slice_to(
    modulus: Barrett50Modulus,
    scalar: u64,
    b: &[u64],
    c: &[u64],
    output: &mut [u64],
) {
    debug_assert_eq!(b.len(), c.len());
    debug_assert_eq!(b.len(), output.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);
    let sv = _mm512_set1_epi64(scalar as i64);

    let chunks = b.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let bv = unsafe { load_u64x8(&b[off..]) };
        let cv = unsafe { load_u64x8(&c[off..]) };
        let (x_lo, x_hi) = ifma_widening_mul(sv, bv);
        let r = ifma_lazy_reduce_4m(x_lo, x_hi, mu_lo52, mu_hi, neg_m, pow52_mask);
        let r = _mm512_add_epi64(r, cv);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        unsafe { store_u64x8(&mut output[off..], r) };
    }
    modulus.inner.reduce_scalar_mul_add_slice_to(
        scalar,
        &b[main_len..],
        &c[main_len..],
        &mut output[main_len..],
    );
}

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn lazy_reduce_scalar_mul_add_slice_to(
    modulus: Barrett50Modulus,
    scalar: u64,
    b: &[u64],
    c: &[u64],
    output: &mut [u64],
) {
    debug_assert_eq!(b.len(), c.len());
    debug_assert_eq!(b.len(), output.len());
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);
    let sv = _mm512_set1_epi64(scalar as i64);

    let chunks = b.len() / N;
    let main_len = chunks * N;
    for i in 0..chunks {
        let off = i * N;
        let bv = unsafe { load_u64x8(&b[off..]) };
        let cv = unsafe { load_u64x8(&c[off..]) };
        let (x_lo, x_hi) = ifma_widening_mul(sv, bv);
        let r = ifma_lazy_reduce_4m(x_lo, x_hi, mu_lo52, mu_hi, neg_m, pow52_mask);
        let r = _mm512_add_epi64(r, cv);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        unsafe { store_u64x8(&mut output[off..], r) };
    }
    modulus.inner.lazy_reduce_scalar_mul_add_slice_to(
        scalar,
        &b[main_len..],
        &c[main_len..],
        &mut output[main_len..],
    );
}

// ===========================================================================
// reduce_dot_product
//
// Outer chunk = K · N elements. Inside one outer chunk we accumulate K
// widening products into a (lo52, hi52) double-word per lane via IFMA
// `madd52lo`/`madd52hi` (the second IFMA round adds to the existing
// accumulator, so we get K_FMA per lane for free). Once the outer chunk
// is done, fold the double-word back to a single SIMD word via the same
// `ifma_lazy_reduce_4m` machinery, fully canonical-reduce to `[0, m)`,
// then cross-chunk-accumulate into a running canonical SIMD acc.
//
// Hi-limb safety: each lane-wise widening product is in `[0, m²) ⊂
// [0, 2^100)`. Accumulating K = 16 products keeps the running 104-bit
// double-word strictly below `K · m² < 2^104`, so `(hi52, lo52)` does
// not overflow. The same bound is used by the existing barrett/simd.rs
// dot_product.
// ===========================================================================

#[target_feature(enable = "avx512f,avx512dq,avx512ifma")]
pub unsafe fn reduce_dot_product(modulus: Barrett50Modulus, a: &[u64], b: &[u64]) -> u64 {
    assert_eq!(a.len(), b.len(), "reduce_dot_product: length mismatch");
    let (m, mu_lo52, mu_hi, neg_m, pow52_mask) = splat_params(&modulus);

    let outer = DOT_INNER_CHUNK * N;
    let mut total_acc = _mm512_setzero_si512();

    let outer_chunks = a.len() / outer;
    for ch in 0..outer_chunks {
        let base = ch * outer;
        let mut lo = _mm512_setzero_si512();
        let mut hi = _mm512_setzero_si512();
        for k in 0..DOT_INNER_CHUNK {
            let off = base + k * N;
            let av = unsafe { load_u64x8(&a[off..]) };
            let bv = unsafe { load_u64x8(&b[off..]) };
            // Two-stage IFMA accumulation: `madd52lo/hi` add the new
            // product limbs straight into the running accumulator.
            lo = _mm512_madd52lo_epu64(lo, av, bv);
            hi = _mm512_madd52hi_epu64(hi, av, bv);
        }
        // The lo accumulator may have spilled into bits 52+ (16 additions
        // of 52-bit values fit in 56 bits). Carry the overflow into hi
        // before reducing.
        let carry = _mm512_srli_epi64::<52>(lo);
        let lo = _mm512_and_si512(lo, pow52_mask);
        let hi = _mm512_add_epi64(hi, carry);

        let r = ifma_lazy_reduce_4m(lo, hi, mu_lo52, mu_hi, neg_m, pow52_mask);
        // r ∈ [0, 4m): canonicalize before the cross-chunk accumulator.
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        let r = cond_sub_m(r, m);
        // total_acc ∈ [0, m), r ∈ [0, m) → sum ∈ [0, 2m) → one cond-sub
        // brings it back to canonical [0, m).
        let sum = _mm512_add_epi64(total_acc, r);
        total_acc = cond_sub_m(sum, m);
    }

    // Horizontal modular sum across the 8 lanes.
    let lanes: [u64; N] = {
        let mut buf = [0u64; N];
        unsafe { _mm512_storeu_si512(buf.as_mut_ptr().cast(), total_acc) };
        buf
    };
    let mut result: u64 = 0;
    for &v in lanes.iter() {
        result = modulus.inner.reduce_add(result, v);
    }

    let tail_len = a.len() - outer_chunks * outer;
    if tail_len > 0 {
        let tail_start = outer_chunks * outer;
        let tail = modulus
            .inner
            .reduce_dot_product(&a[tail_start..], &b[tail_start..]);
        result = modulus.inner.reduce_add(result, tail);
    }
    result
}
