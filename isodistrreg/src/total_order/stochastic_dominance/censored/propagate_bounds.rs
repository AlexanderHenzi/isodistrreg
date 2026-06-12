//! Inner-loop reduction kernels for `propagate_bounds`.
//!
//! `propagate_bounds` reduces over a pair of equal-length slices to compute
//! `(max_k min(row[k], col[k]), min_k max(row[k], col[k]))`. It is the hot spot of the
//! censored-data PAVA, so each architecture uses the widest SIMD its compile-time baseline —
//! and, on x86_64, its runtime feature set — allows.
//!
//! The `Kernel` trait is the dispatch seam; `dispatch_generalized_pava` (in `fast.rs`) selects
//! an implementor at runtime by feature-detecting AVX-512F then AVX2, falling back to
//! `ScalarKernel`.
//!
//! All kernels compute the **exact** reduction over the full slices — there is no early
//! termination. Mean reduction lengths in real fits are a few dozen to a few hundred
//! elements, where mid-loop collapse checks (horizontal reduce + compare per chunk) cost more
//! than they save; dropping them also makes the result independent of chunking, so all three
//! kernels return bit-identical bounds and fit results no longer vary with the host's SIMD
//! level (collapsed-bounds midpoints used to depend on where in the row the early exit
//! fired). `kernels_agree` locks the cross-kernel equivalence in.
//!
//! ## Why three implementations
//!
//! The split mirrors what each pre-compiled binary may legally execute:
//!
//! * `Avx512Kernel` / `Avx2Kernel` — explicit intrinsics behind `#[target_feature]`. AVX-512F
//!   and AVX2 sit *above* the x86_64 baseline (SSE2), so a pre-compiled binary (CRAN, distro
//!   packages) cannot use them unconditionally; the runtime feature gate is what unlocks them.
//! * `ScalarKernel` → `propagate_bounds_kernel` — a manually unrolled scalar loop whose SIMD
//!   comes from the auto-vectorizer. SLP fuses each per-lane group into one SIMD min/max at
//!   the target's baseline width:
//!     - aarch64 (incl. Apple Silicon): NEON `fmin`/`fmax` are mandatory in ARMv8, so a
//!       pre-compiled aarch64 binary reaches 128-bit 4-wide f32 with no runtime gate.
//!     - x86_64 with neither AVX2 nor AVX-512 detected: SSE2 is baseline since 2003, so SLP
//!       gets 128-bit `vminps`/`vmaxps`, 4-wide.
//!
//! ## Why bare `min`/`max`, not `f32::min`/`max` or `core::simd::SimdFloat`
//!
//! The reduced values are survival probabilities and never NaN. The `if a < b { a } else { b }`
//! pattern lowers to a single hardware `minss`/`maxss` / `minps`/`maxps` (NaN-propagating),
//! which matches the algorithm. The IEEE-754-clean alternatives lower to `llvm.minnum` /
//! `llvm.maxnum` and emit a `vcmpunordps`-based cleanup after every op; in our benchmarks the
//! cleanup overhead exceeded the width gain.

/// Compile-time selector for the inner reduction kernel. Implementors are zero-sized; the
/// associated function is statically dispatched, so monomorphizing the algorithm tree over
/// `K: Kernel` inlines `K::apply` into every `propagate_bounds*` callsite.
pub trait Kernel {
    fn apply(row: &[f32], col: &[f32]) -> (f32, f32);
}

/// Portable scalar / SLP-vectorized fallback. Always usable.
pub struct ScalarKernel;

impl Kernel for ScalarKernel {
    #[inline(always)]
    fn apply(row: &[f32], col: &[f32]) -> (f32, f32) {
        propagate_bounds_kernel(row, col)
    }
}

#[inline(always)]
fn fast_min(a: f32, b: f32) -> f32 {
    if a < b { a } else { b }
}
#[inline(always)]
fn fast_max(a: f32, b: f32) -> f32 {
    if a > b { a } else { b }
}

/// Inner kernel of `propagate_bounds`: the exact reduction
/// `(max_k min(row[k], col[k]), min_k max(row[k], col[k]))`.
///
/// Targets baseline SIMD only — SSE2 on x86_64, NEON on aarch64. The body keeps four
/// independent `(lower, upper)` accumulator lanes; with contiguous `[f32; 4]` reads SLP fuses
/// each per-lane group into one 128-bit `vminps`/`vmaxps` (or NEON `fmin`/`fmax`), and the
/// four lanes only meet in the final horizontal reduction, keeping dependency chains short.
const LANES: usize = 4;

#[inline]
fn propagate_bounds_kernel(row: &[f32], col: &[f32]) -> (f32, f32) {
    debug_assert_eq!(row.len(), col.len());
    let n = row.len();

    let mut lo = [f32::NEG_INFINITY; LANES];
    let mut hi = [f32::INFINITY; LANES];
    let mut k = 0;

    while k + LANES <= n {
        // Borrow the 4-element groups as fixed-size array references so the compiler
        // hoists the bounds check above the body — without this, each scalar
        // `row[k + j]` re-checks.
        let row_g: &[f32; LANES] = (&row[k..k + LANES]).try_into().unwrap();
        let col_g: &[f32; LANES] = (&col[k..k + LANES]).try_into().unwrap();
        for j in 0..LANES {
            let l = fast_min(row_g[j], col_g[j]);
            let u = fast_max(row_g[j], col_g[j]);
            lo[j] = fast_max(lo[j], l);
            hi[j] = fast_min(hi[j], u);
        }
        k += LANES;
    }
    while k < n {
        let l = fast_min(row[k], col[k]);
        let u = fast_max(row[k], col[k]);
        lo[0] = fast_max(lo[0], l);
        hi[0] = fast_min(hi[0], u);
        k += 1;
    }

    let lower = fast_max(fast_max(lo[0], lo[1]), fast_max(lo[2], lo[3]));
    let upper = fast_min(fast_min(hi[0], hi[1]), fast_min(hi[2], hi[3]));
    (lower, upper)
}

/// AVX2 kernel — bare `vminps`/`vmaxps` on `ymm` (8-wide f32). Only safe to use after
/// `is_x86_feature_detected!("avx2")` confirmed it is available.
#[cfg(target_arch = "x86_64")]
pub struct Avx2Kernel;

#[cfg(target_arch = "x86_64")]
impl Kernel for Avx2Kernel {
    #[inline(always)]
    fn apply(row: &[f32], col: &[f32]) -> (f32, f32) {
        // SAFETY: only monomorphized into the algorithm tree by `dispatch_generalized_pava`
        // after `is_x86_feature_detected!("avx2")` returned true.
        unsafe { propagate_bounds_kernel_avx2(row, col) }
    }
}

/// Bare-AVX2 inner kernel — 8-wide f32, two independent accumulator pairs for ILP, masked
/// (`vmaskmovps` + blend) tail so sub-8 remainders stay vectorized.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx,avx2")]
unsafe fn propagate_bounds_kernel_avx2(row: &[f32], col: &[f32]) -> (f32, f32) {
    use core::arch::x86_64::{
        __m256, _mm_cvtss_f32, _mm_max_ps, _mm_min_ps, _mm_movehdup_ps, _mm_movehl_ps,
        _mm256_blendv_ps, _mm256_castps256_ps128, _mm256_castsi256_ps, _mm256_cmpgt_epi32,
        _mm256_extractf128_ps, _mm256_loadu_ps, _mm256_maskload_ps, _mm256_max_ps, _mm256_min_ps,
        _mm256_set1_epi32, _mm256_set1_ps, _mm256_setr_epi32,
    };

    const F32_LANES: usize = 8;

    debug_assert_eq!(row.len(), col.len());
    let n = row.len();
    let mut k = 0;

    // SAFETY: AVX/AVX2 enabled by function-level `target_feature`. Pointer adds stay
    // in-bounds because each step tests the remaining length first; the masked tail loads
    // only lanes `< n - k`.
    unsafe {
        let neg_inf = _mm256_set1_ps(f32::NEG_INFINITY);
        let pos_inf = _mm256_set1_ps(f32::INFINITY);
        let mut lo0 = neg_inf;
        let mut hi0 = pos_inf;
        let mut lo1 = neg_inf;
        let mut hi1 = pos_inf;

        let row_ptr = row.as_ptr();
        let col_ptr = col.as_ptr();

        while k + 2 * F32_LANES <= n {
            let r0 = _mm256_loadu_ps(row_ptr.add(k));
            let c0 = _mm256_loadu_ps(col_ptr.add(k));
            let r1 = _mm256_loadu_ps(row_ptr.add(k + F32_LANES));
            let c1 = _mm256_loadu_ps(col_ptr.add(k + F32_LANES));
            lo0 = _mm256_max_ps(lo0, _mm256_min_ps(r0, c0));
            hi0 = _mm256_min_ps(hi0, _mm256_max_ps(r0, c0));
            lo1 = _mm256_max_ps(lo1, _mm256_min_ps(r1, c1));
            hi1 = _mm256_min_ps(hi1, _mm256_max_ps(r1, c1));
            k += 2 * F32_LANES;
        }
        if k + F32_LANES <= n {
            let r = _mm256_loadu_ps(row_ptr.add(k));
            let c = _mm256_loadu_ps(col_ptr.add(k));
            lo0 = _mm256_max_ps(lo0, _mm256_min_ps(r, c));
            hi0 = _mm256_min_ps(hi0, _mm256_max_ps(r, c));
            k += F32_LANES;
        }
        if k < n {
            // Masked tail: load the remaining `n - k` lanes, then blend the inactive lanes
            // to (-inf, +inf) so their (min, max) pair is the reduction identity.
            let rem = _mm256_set1_epi32((n - k) as i32);
            let idx = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
            let mask = _mm256_cmpgt_epi32(rem, idx);
            let r = _mm256_maskload_ps(row_ptr.add(k), mask);
            let c = _mm256_maskload_ps(col_ptr.add(k), mask);
            let mask_ps = _mm256_castsi256_ps(mask);
            let r = _mm256_blendv_ps(neg_inf, r, mask_ps);
            let c = _mm256_blendv_ps(pos_inf, c, mask_ps);
            lo1 = _mm256_max_ps(lo1, _mm256_min_ps(r, c));
            hi1 = _mm256_min_ps(hi1, _mm256_max_ps(r, c));
        }

        let lo: __m256 = _mm256_max_ps(lo0, lo1);
        let hi: __m256 = _mm256_min_ps(hi0, hi1);

        // Horizontal reduction: 256 -> 128 -> 64 -> 32 bits via max/min + shuffles.
        let lo128 = _mm_max_ps(_mm256_castps256_ps128(lo), _mm256_extractf128_ps(lo, 1));
        let hi128 = _mm_min_ps(_mm256_castps256_ps128(hi), _mm256_extractf128_ps(hi, 1));
        let lo64 = _mm_max_ps(lo128, _mm_movehl_ps(lo128, lo128));
        let hi64 = _mm_min_ps(hi128, _mm_movehl_ps(hi128, hi128));
        let lo32 = _mm_max_ps(lo64, _mm_movehdup_ps(lo64));
        let hi32 = _mm_min_ps(hi64, _mm_movehdup_ps(hi64));
        (_mm_cvtss_f32(lo32), _mm_cvtss_f32(hi32))
    }
}

/// AVX-512F kernel — bare `vminps`/`vmaxps` on `zmm` (16-wide f32). Only safe to use after
/// `is_x86_feature_detected!("avx512f")` confirmed it is available.
#[cfg(target_arch = "x86_64")]
pub struct Avx512Kernel;

#[cfg(target_arch = "x86_64")]
impl Kernel for Avx512Kernel {
    #[inline(always)]
    fn apply(row: &[f32], col: &[f32]) -> (f32, f32) {
        // SAFETY: only monomorphized into the algorithm tree by `dispatch_generalized_pava`
        // after `is_x86_feature_detected!("avx512f")` returned true.
        unsafe { propagate_bounds_kernel_avx512(row, col) }
    }
}

/// AVX-512F inner kernel — 16-wide f32, two independent accumulator pairs for ILP, native
/// masked-load tail so sub-16 remainders stay vectorized.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn propagate_bounds_kernel_avx512(row: &[f32], col: &[f32]) -> (f32, f32) {
    use core::arch::x86_64::{
        __m512, _mm_cvtss_f32, _mm_max_ps, _mm_min_ps, _mm_movehdup_ps, _mm_movehl_ps,
        _mm512_castps512_ps128, _mm512_loadu_ps, _mm512_mask_loadu_ps, _mm512_max_ps,
        _mm512_min_ps, _mm512_set1_ps, _mm512_shuffle_f32x4,
    };

    const F32_LANES: usize = 16;

    debug_assert_eq!(row.len(), col.len());
    let n = row.len();
    let mut k = 0;

    // SAFETY: AVX-512F enforced by function-level `target_feature`. Pointer adds stay
    // in-bounds because each step tests the remaining length first; the masked tail loads
    // only lanes `< n - k`.
    unsafe {
        let neg_inf = _mm512_set1_ps(f32::NEG_INFINITY);
        let pos_inf = _mm512_set1_ps(f32::INFINITY);
        let mut lo0 = neg_inf;
        let mut hi0 = pos_inf;
        let mut lo1 = neg_inf;
        let mut hi1 = pos_inf;

        let row_ptr = row.as_ptr();
        let col_ptr = col.as_ptr();

        while k + 2 * F32_LANES <= n {
            let r0 = _mm512_loadu_ps(row_ptr.add(k));
            let c0 = _mm512_loadu_ps(col_ptr.add(k));
            let r1 = _mm512_loadu_ps(row_ptr.add(k + F32_LANES));
            let c1 = _mm512_loadu_ps(col_ptr.add(k + F32_LANES));
            lo0 = _mm512_max_ps(lo0, _mm512_min_ps(r0, c0));
            hi0 = _mm512_min_ps(hi0, _mm512_max_ps(r0, c0));
            lo1 = _mm512_max_ps(lo1, _mm512_min_ps(r1, c1));
            hi1 = _mm512_min_ps(hi1, _mm512_max_ps(r1, c1));
            k += 2 * F32_LANES;
        }
        if k + F32_LANES <= n {
            let r = _mm512_loadu_ps(row_ptr.add(k));
            let c = _mm512_loadu_ps(col_ptr.add(k));
            lo0 = _mm512_max_ps(lo0, _mm512_min_ps(r, c));
            hi0 = _mm512_min_ps(hi0, _mm512_max_ps(r, c));
            k += F32_LANES;
        }
        if k < n {
            // Masked tail: inactive row lanes read -inf and col lanes +inf, making their
            // (min, max) pair the reduction identity.
            let mask = (1u16 << (n - k)) - 1;
            let r = _mm512_mask_loadu_ps(neg_inf, mask, row_ptr.add(k));
            let c = _mm512_mask_loadu_ps(pos_inf, mask, col_ptr.add(k));
            lo1 = _mm512_max_ps(lo1, _mm512_min_ps(r, c));
            hi1 = _mm512_min_ps(hi1, _mm512_max_ps(r, c));
        }

        let lo: __m512 = _mm512_max_ps(lo0, lo1);
        let hi: __m512 = _mm512_min_ps(hi0, hi1);

        // Horizontal reduction: fold the four 128-bit blocks together (AVX-512F-only —
        // no DQ extracts), then reduce within 128 bits.
        let lo_sw = _mm512_shuffle_f32x4(lo, lo, 0b01_00_11_10);
        let hi_sw = _mm512_shuffle_f32x4(hi, hi, 0b01_00_11_10);
        let lo2 = _mm512_max_ps(lo, lo_sw);
        let hi2 = _mm512_min_ps(hi, hi_sw);
        let lo_sw2 = _mm512_shuffle_f32x4(lo2, lo2, 0b10_11_00_01);
        let hi_sw2 = _mm512_shuffle_f32x4(hi2, hi2, 0b10_11_00_01);
        let lo3 = _mm512_max_ps(lo2, lo_sw2);
        let hi3 = _mm512_min_ps(hi2, hi_sw2);
        let lo128 = _mm512_castps512_ps128(lo3);
        let hi128 = _mm512_castps512_ps128(hi3);
        let lo64 = _mm_max_ps(lo128, _mm_movehl_ps(lo128, lo128));
        let hi64 = _mm_min_ps(hi128, _mm_movehl_ps(hi128, hi128));
        let lo32 = _mm_max_ps(lo64, _mm_movehdup_ps(lo64));
        let hi32 = _mm_min_ps(hi64, _mm_movehdup_ps(hi64));
        (_mm_cvtss_f32(lo32), _mm_cvtss_f32(hi32))
    }
}

#[cfg(test)]
mod test {
    use super::*;

    /// Deterministic pseudo-random f32 in (0, 1] — the value domain of the survival kernel.
    fn lcg_f32(state: &mut u64) -> f32 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (((*state >> 33) as u32) as f32 + 1.0) / (u32::MAX as f32 + 1.0)
    }

    /// All kernels are exact reductions, so they must agree bit-for-bit with the naive
    /// definition (and therefore with each other) on every length, including empty and
    /// tail-only slices. This is what guarantees fit results do not depend on the host's
    /// SIMD level.
    #[test]
    fn kernels_agree() {
        let mut state = 0x5eed;
        for n in 0..=257 {
            let row: Vec<f32> = (0..n).map(|_| lcg_f32(&mut state)).collect();
            let col: Vec<f32> = (0..n).map(|_| lcg_f32(&mut state)).collect();

            let mut expected = (f32::NEG_INFINITY, f32::INFINITY);
            for k in 0..n {
                let l = if row[k] < col[k] { row[k] } else { col[k] };
                let u = if row[k] > col[k] { row[k] } else { col[k] };
                expected.0 = if expected.0 > l { expected.0 } else { l };
                expected.1 = if expected.1 < u { expected.1 } else { u };
            }

            assert_eq!(ScalarKernel::apply(&row, &col), expected, "scalar, n={n}");
            #[cfg(target_arch = "x86_64")]
            {
                if is_x86_feature_detected!("avx2") {
                    assert_eq!(Avx2Kernel::apply(&row, &col), expected, "avx2, n={n}");
                }
                if is_x86_feature_detected!("avx512f") {
                    assert_eq!(Avx512Kernel::apply(&row, &col), expected, "avx512, n={n}");
                }
            }
        }
    }
}
