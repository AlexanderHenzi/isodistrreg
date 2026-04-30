//! Inner-loop reduction kernels for `propagate_bounds`.
//!
//! `propagate_bounds` reduces over a pair of equal-length slices to compute
//! `(max_k min(row[k], col[k]), min_k max(row[k], col[k]))` with an early exit when the
//! running window collapses (`lower >= upper`). It is the hot spot of the censored-data PAVA,
//! so each architecture uses the widest SIMD its compile-time baseline — and, on x86_64, its
//! runtime feature set — allows.
//!
//! The `Kernel` trait is the dispatch seam; `dispatch_generalized_pava` selects an implementor
//! at runtime by feature-detecting AVX-512F then AVX2, falling back to `ScalarKernel`.
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
//!       pre-compiled aarch64 binary reaches 128-bit 2-wide f64 with no runtime gate.
//!     - x86_64 with neither AVX2 nor AVX-512 detected: SSE2 is baseline since 2003, so SLP
//!       gets 128-bit `minpd`/`maxpd`.
//!
//!   On x86_64 hardware where AVX2 or AVX-512 *is* detected, `propagate_bounds_kernel` is
//!   never reached — the explicit kernels handle it.
//!
//! ## Why bare `min`/`max`, not `f64::min`/`max` or `core::simd::SimdFloat`
//!
//! The reduced values are survival probabilities and never NaN. The `if a < b { a } else { b }`
//! pattern lowers to a single hardware `minsd`/`maxsd` / `minpd`/`maxpd` (NaN-propagating),
//! which matches the algorithm. The IEEE-754-clean alternatives lower to `llvm.minnum` /
//! `llvm.maxnum` and emit a `vcmpunordpd`-based cleanup after every op; in our benchmarks the
//! cleanup overhead exceeded the width gain (≈+10% on the indep n=1000 hot path when forcing
//! 256-bit `ymm` via `core::simd`).

/// Compile-time selector for the inner reduction kernel. Implementors are zero-sized; the
/// associated function is statically dispatched, so monomorphizing the algorithm tree over
/// `K: Kernel` inlines `K::apply` into every `propagate_bounds*` callsite.
pub trait Kernel {
    fn apply(row: &[f64], col: &[f64]) -> (f64, f64);
}

/// Portable scalar / SLP-vectorized fallback. Always usable.
pub struct ScalarKernel;

impl Kernel for ScalarKernel {
    #[inline(always)]
    fn apply(row: &[f64], col: &[f64]) -> (f64, f64) {
        propagate_bounds_kernel(row, col)
    }
}

// Survival probabilities are never NaN; the simple-comparison forms below lower to bare
// `minsd`/`maxsd` and should auto-vectorize into `vminpd`/`vmaxpd` for the chunk body.
#[inline(always)]
fn fast_min(a: f64, b: f64) -> f64 {
    if a < b { a } else { b }
}
#[inline(always)]
fn fast_max(a: f64, b: f64) -> f64 {
    if a > b { a } else { b }
}

/// Inner kernel of `propagate_bounds`: computes
/// `(max_k min(row[k], col[k]), min_k max(row[k], col[k]))` with an early exit when the
/// running bounds collapse.
///
/// Targets baseline SIMD only — SSE2 on x86_64, NEON on aarch64. The chunk body keeps four
/// independent `(lower, upper)` accumulator lanes; with contiguous `[f64; 4]` reads SLP fuses
/// each per-lane group into one 128-bit `vminpd`/`vmaxpd` (or NEON `fmin`/`fmax`). The
/// bound-collapse check happens between chunks (and per-element in the tail), so the
/// early-termination win is preserved at chunk granularity. Above-baseline x86_64 features
/// (AVX2, AVX-512F) are handled by the explicit-intrinsic kernels in this module instead.
const LANES: usize = 4;
const GROUPS_PER_CHUNK: usize = 4;
const PROPAGATE_BOUNDS_CHUNK: usize = LANES * GROUPS_PER_CHUNK;

#[inline]
fn propagate_bounds_kernel(row: &[f64], col: &[f64]) -> (f64, f64) {
    debug_assert_eq!(row.len(), col.len());
    let n = row.len();

    let mut lower = f64::NEG_INFINITY;
    let mut upper = f64::INFINITY;
    let mut k = 0;

    while k + PROPAGATE_BOUNDS_CHUNK <= n {
        // Four independent accumulator lanes. Updates within a group are mutually independent
        // and the four cross-group chains never alias, so SLP fuses each per-lane group into
        // one SIMD min/max at the target's baseline width (128-bit on SSE2 / aarch64 NEON).
        // The four chains also keep per-lane dependency height short for ILP. See module
        // docs for why we use bare `min`/`max` and why this kernel only targets baseline SIMD.
        let mut lo0 = f64::NEG_INFINITY;
        let mut lo1 = f64::NEG_INFINITY;
        let mut lo2 = f64::NEG_INFINITY;
        let mut lo3 = f64::NEG_INFINITY;
        let mut hi0 = f64::INFINITY;
        let mut hi1 = f64::INFINITY;
        let mut hi2 = f64::INFINITY;
        let mut hi3 = f64::INFINITY;

        for g in 0..GROUPS_PER_CHUNK {
            let off = k + g * LANES;
            // Borrow the 4-element groups as fixed-size array references so the compiler
            // hoists the bounds check above the body — without this, each scalar
            // `row[off + j]` re-checks.
            let row_g: &[f64; LANES] = (&row[off..off + LANES]).try_into().unwrap();
            let col_g: &[f64; LANES] = (&col[off..off + LANES]).try_into().unwrap();
            let r0 = row_g[0];
            let r1 = row_g[1];
            let r2 = row_g[2];
            let r3 = row_g[3];
            let c0 = col_g[0];
            let c1 = col_g[1];
            let c2 = col_g[2];
            let c3 = col_g[3];

            let l0 = fast_min(r0, c0);
            let l1 = fast_min(r1, c1);
            let l2 = fast_min(r2, c2);
            let l3 = fast_min(r3, c3);
            let u0 = fast_max(r0, c0);
            let u1 = fast_max(r1, c1);
            let u2 = fast_max(r2, c2);
            let u3 = fast_max(r3, c3);

            lo0 = fast_max(lo0, l0);
            lo1 = fast_max(lo1, l1);
            lo2 = fast_max(lo2, l2);
            lo3 = fast_max(lo3, l3);
            hi0 = fast_min(hi0, u0);
            hi1 = fast_min(hi1, u1);
            hi2 = fast_min(hi2, u2);
            hi3 = fast_min(hi3, u3);
        }

        // Horizontal reduction across the four lanes.
        let chunk_lower = fast_max(fast_max(lo0, lo1), fast_max(lo2, lo3));
        let chunk_upper = fast_min(fast_min(hi0, hi1), fast_min(hi2, hi3));
        lower = fast_max(lower, chunk_lower);
        upper = fast_min(upper, chunk_upper);
        k += PROPAGATE_BOUNDS_CHUNK;
        if lower >= upper {
            return (lower, upper);
        }
    }
    while k < n {
        let a = row[k];
        let b = col[k];
        let lo = fast_min(a, b);
        let hi = fast_max(a, b);
        lower = fast_max(lower, lo);
        upper = fast_min(upper, hi);
        k += 1;
        if lower >= upper {
            return (lower, upper);
        }
    }
    (lower, upper)
}

/// AVX2 kernel — bare `vminpd`/`vmaxpd` on `ymm`. Only safe to use after `is_x86_feature_detected!`
/// confirmed AVX2 is available.
#[cfg(target_arch = "x86_64")]
pub struct Avx2Kernel;

#[cfg(target_arch = "x86_64")]
impl Kernel for Avx2Kernel {
    #[inline(always)]
    fn apply(row: &[f64], col: &[f64]) -> (f64, f64) {
        // SAFETY: this impl is only ever monomorphized into the algorithm tree by
        // `dispatch_generalized_pava` after `is_x86_feature_detected!("avx2")` returned true.
        unsafe { propagate_bounds_kernel_avx2(row, col) }
    }
}

/// Bare-AVX2 inner kernel — calls `_mm256_min_pd`/`_mm256_max_pd` (`vminpd`/`vmaxpd` on `ymm`)
/// without IEEE-754 NaN cleanup. Survival probabilities are guaranteed non-NaN, so the bare
/// hardware semantics (NaN-propagating min/max) match the algorithm's needs.
///
/// SAFETY: caller must ensure the CPU supports AVX2. The runtime gate in `dispatch_kernel`
/// is the only place this is selected.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx,avx2")]
unsafe fn propagate_bounds_kernel_avx2(row: &[f64], col: &[f64]) -> (f64, f64) {
    use core::arch::x86_64::{
        _mm256_loadu_pd, _mm256_max_pd, _mm256_min_pd, _mm256_set1_pd, _mm256_storeu_pd,
    };

    debug_assert_eq!(row.len(), col.len());
    let n = row.len();

    let mut k = 0;

    // SAFETY of the intrinsic calls in the unsafe block below: the caller must ensure AVX2 is
    // available, which is expressed by the `#[target_feature(enable = "avx,avx2")]` annotation
    // and enforced by `dispatch_kernel` which only selects this function when
    // `is_x86_feature_detected!("avx2")` returned true. The pointer arithmetic via `add(off)`
    // stays in-bounds because each chunk tests `k + PROPAGATE_BOUNDS_CHUNK <= n` first, so all
    // four `4×f64` group loads at offsets `k, k+4, k+8, k+12` lie inside the slice. The store
    // targets are stack-local 4-element arrays.
    let (mut lower, mut upper) = unsafe {
        let mut lo_v = _mm256_set1_pd(f64::NEG_INFINITY);
        let mut hi_v = _mm256_set1_pd(f64::INFINITY);

        while k + PROPAGATE_BOUNDS_CHUNK <= n {
            let row_ptr = row.as_ptr();
            let col_ptr = col.as_ptr();
            for g in 0..GROUPS_PER_CHUNK {
                let off = k + g * LANES;
                let r = _mm256_loadu_pd(row_ptr.add(off));
                let c = _mm256_loadu_pd(col_ptr.add(off));
                let lo_pair = _mm256_min_pd(r, c);
                let hi_pair = _mm256_max_pd(r, c);
                lo_v = _mm256_max_pd(lo_v, lo_pair);
                hi_v = _mm256_min_pd(hi_v, hi_pair);
            }
            k += PROPAGATE_BOUNDS_CHUNK;

            // Horizontal reduce + check early termination.
            let mut lo_arr = [0.0f64; LANES];
            let mut hi_arr = [0.0f64; LANES];
            _mm256_storeu_pd(lo_arr.as_mut_ptr(), lo_v);
            _mm256_storeu_pd(hi_arr.as_mut_ptr(), hi_v);
            let cur_lo = fast_max(
                fast_max(lo_arr[0], lo_arr[1]),
                fast_max(lo_arr[2], lo_arr[3]),
            );
            let cur_hi = fast_min(
                fast_min(hi_arr[0], hi_arr[1]),
                fast_min(hi_arr[2], hi_arr[3]),
            );
            if cur_lo >= cur_hi {
                return (cur_lo, cur_hi);
            }
        }

        let mut lo_arr = [0.0f64; LANES];
        let mut hi_arr = [0.0f64; LANES];
        _mm256_storeu_pd(lo_arr.as_mut_ptr(), lo_v);
        _mm256_storeu_pd(hi_arr.as_mut_ptr(), hi_v);
        let lower = fast_max(
            fast_max(lo_arr[0], lo_arr[1]),
            fast_max(lo_arr[2], lo_arr[3]),
        );
        let upper = fast_min(
            fast_min(hi_arr[0], hi_arr[1]),
            fast_min(hi_arr[2], hi_arr[3]),
        );
        (lower, upper)
    };

    while k < n {
        let a = row[k];
        let b = col[k];
        let lo = fast_min(a, b);
        let hi = fast_max(a, b);
        lower = fast_max(lower, lo);
        upper = fast_min(upper, hi);
        k += 1;
        if lower >= upper {
            return (lower, upper);
        }
    }

    (lower, upper)
}

/// AVX-512F kernel — bare `vminpd`/`vmaxpd` on `zmm` (8-wide f64). Only safe to use after
/// `is_x86_feature_detected!("avx512f")` confirmed it is available.
#[cfg(target_arch = "x86_64")]
pub struct Avx512Kernel;

#[cfg(target_arch = "x86_64")]
impl Kernel for Avx512Kernel {
    #[inline(always)]
    fn apply(row: &[f64], col: &[f64]) -> (f64, f64) {
        // SAFETY: only monomorphized into the algorithm tree by `dispatch_generalized_pava` after
        // `is_x86_feature_detected!("avx512f")` returned true.
        unsafe { propagate_bounds_kernel_avx512(row, col) }
    }
}

/// AVX-512F kernel — same shape as the AVX2 variant but 8-wide. Each chunk holds 16 elements
/// processed as 2 groups of 8, giving 4 zmm load pairs + 4 `vminpd`/`vmaxpd` (zmm) per chunk.
/// SAFETY: caller must ensure AVX-512F is available; enforced by `dispatch_generalized_pava`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn propagate_bounds_kernel_avx512(row: &[f64], col: &[f64]) -> (f64, f64) {
    use core::arch::x86_64::{
        _mm512_loadu_pd, _mm512_max_pd, _mm512_min_pd, _mm512_set1_pd, _mm512_storeu_pd,
    };

    debug_assert_eq!(row.len(), col.len());
    let n = row.len();

    const AVX512_LANES: usize = 8;
    const AVX512_GROUPS_PER_CHUNK: usize = 2;
    const AVX512_CHUNK: usize = AVX512_LANES * AVX512_GROUPS_PER_CHUNK; // 16, matches AVX2 chunk

    #[inline(always)]
    fn reduce_max8(a: [f64; 8]) -> f64 {
        let m01 = if a[0] > a[1] { a[0] } else { a[1] };
        let m23 = if a[2] > a[3] { a[2] } else { a[3] };
        let m45 = if a[4] > a[5] { a[4] } else { a[5] };
        let m67 = if a[6] > a[7] { a[6] } else { a[7] };
        let m0123 = if m01 > m23 { m01 } else { m23 };
        let m4567 = if m45 > m67 { m45 } else { m67 };
        if m0123 > m4567 { m0123 } else { m4567 }
    }
    #[inline(always)]
    fn reduce_min8(a: [f64; 8]) -> f64 {
        let m01 = if a[0] < a[1] { a[0] } else { a[1] };
        let m23 = if a[2] < a[3] { a[2] } else { a[3] };
        let m45 = if a[4] < a[5] { a[4] } else { a[5] };
        let m67 = if a[6] < a[7] { a[6] } else { a[7] };
        let m0123 = if m01 < m23 { m01 } else { m23 };
        let m4567 = if m45 < m67 { m45 } else { m67 };
        if m0123 < m4567 { m0123 } else { m4567 }
    }

    let mut k = 0;

    // SAFETY: AVX-512F enforced by the function-level `target_feature`. Pointer adds stay
    // in-bounds because each chunk tests `k + AVX512_CHUNK <= n`. Stores target stack arrays.
    let (mut lower, mut upper) = unsafe {
        let mut lo_v = _mm512_set1_pd(f64::NEG_INFINITY);
        let mut hi_v = _mm512_set1_pd(f64::INFINITY);

        while k + AVX512_CHUNK <= n {
            let row_ptr = row.as_ptr();
            let col_ptr = col.as_ptr();
            for g in 0..AVX512_GROUPS_PER_CHUNK {
                let off = k + g * AVX512_LANES;
                let r = _mm512_loadu_pd(row_ptr.add(off));
                let c = _mm512_loadu_pd(col_ptr.add(off));
                let lo_pair = _mm512_min_pd(r, c);
                let hi_pair = _mm512_max_pd(r, c);
                lo_v = _mm512_max_pd(lo_v, lo_pair);
                hi_v = _mm512_min_pd(hi_v, hi_pair);
            }
            k += AVX512_CHUNK;

            let mut lo_arr = [0.0f64; AVX512_LANES];
            let mut hi_arr = [0.0f64; AVX512_LANES];
            _mm512_storeu_pd(lo_arr.as_mut_ptr(), lo_v);
            _mm512_storeu_pd(hi_arr.as_mut_ptr(), hi_v);
            let cur_lo = reduce_max8(lo_arr);
            let cur_hi = reduce_min8(hi_arr);
            if cur_lo >= cur_hi {
                return (cur_lo, cur_hi);
            }
        }

        let mut lo_arr = [0.0f64; AVX512_LANES];
        let mut hi_arr = [0.0f64; AVX512_LANES];
        _mm512_storeu_pd(lo_arr.as_mut_ptr(), lo_v);
        _mm512_storeu_pd(hi_arr.as_mut_ptr(), hi_v);
        (reduce_max8(lo_arr), reduce_min8(hi_arr))
    };

    while k < n {
        let a = row[k];
        let b = col[k];
        let lo = fast_min(a, b);
        let hi = fast_max(a, b);
        lower = fast_max(lower, lo);
        upper = fast_min(upper, hi);
        k += 1;
        if lower >= upper {
            return (lower, upper);
        }
    }

    (lower, upper)
}
