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
//! All kernels compute the exact reduction, with collapse checkpoints at the fixed element
//! counts 128, 256, 512, … (doubling): once the running bounds satisfy `lower >= upper` at a
//! checkpoint, the reduction returns them immediately. The crossing is permanent — `lower`
//! only grows and `upper` only shrinks — so an exit fires exactly on cells whose final
//! bounds are collapsed too; cells whose bounds stay open reduce over the full slices and
//! return the exact result. Collapsed cells return the bounds as of the checkpoint, which
//! sandwich the exact clip value just like the final pair does (`lower_k <= lower_final` and
//! `upper_k >= upper_final` cannot restore order once crossed, and in exact arithmetic the
//! bounds never cross at all — the collapse is f32 noise around bound equality).
//!
//! The checkpoint schedule is defined in *element counts*, reachable by every kernel's step
//! size, and a min/max reduction is order-independent — so all three kernels reduce exactly
//! the elements `[0, checkpoint)` when they test, make identical exit decisions, and return
//! bit-identical bounds: fit results do not vary with the host's SIMD level.
//! `kernels_agree` locks the cross-kernel equivalence in, collapsing inputs included.
//! Starting the checks at 128 keeps them free where reductions are short (discrete
//! covariates, strongly correlated data — mean lengths around a hundred); the doubling
//! bounds the check overhead at O(log) per long cell. On continuous-covariate fits with
//! weak X-Y dependence — where mean lengths reach many hundreds and ~3/4 of cell visits
//! collapse — the checkpoints skip roughly half of all reduction work.
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

use crate::total_order::stochastic_dominance::censored::structures::WALK_OBSERVED_BIT;

/// First collapse checkpoint of the bound reductions, in elements; later checkpoints double
/// (256, 512, …). Shared by all kernels — the schedule is part of the numeric contract (it
/// determines which prefix a collapsed cell's returned bounds reduce over), so every kernel
/// must test at exactly these element counts. Must be a multiple of every kernel's main-loop
/// step (32 for AVX-512, 16 for AVX2, 4 for the scalar kernel) so the loops hit each
/// checkpoint exactly.
const COLLAPSE_CHECK_START: usize = 128;

// The checkpoint tests fire on `k == next_check` equality, so a schedule entry that is not
// a multiple of some kernel's step is silently skipped by that kernel only — the kernels
// would then reduce different prefixes and return SIMD-level-dependent bounds, the exact
// failure mode this schedule exists to rule out. 32 is the largest main-loop step
// (AVX-512); the smaller steps (16 AVX2, 4 scalar) divide it, and doubling preserves
// divisibility for every later checkpoint.
const _: () = assert!(
    COLLAPSE_CHECK_START.is_multiple_of(32),
    "checkpoints must be reachable by every kernel's main-loop step",
);

/// Minimum reduction length for the checkpointed kernel mode — the caller's dispatch
/// threshold, kept next to the schedule it reasons about. Below it the straight-line
/// kernel runs: reductions shorter than `COLLAPSE_CHECK_START` cannot reach a checkpoint,
/// and for mid-length reductions the checkpoint reduces on open-bounds cells cost more
/// than the collapsed cells' exits save (their absolute savings are small); both effects
/// measured as net regressions on low-censoring and weakly correlated fits, whose
/// collapsed cells concentrate at these lengths. Long reductions keep the checkpoints:
/// one horizontal reduce is noise next to hundreds of straight-line elements, and exits
/// there skip most of the work.
pub const COLLAPSE_CHECK_MIN_LEN: usize = 384;

/// Compile-time selector for the inner reduction kernel. Implementors are zero-sized; the
/// associated functions are statically dispatched, so monomorphizing the algorithm tree over
/// `K: Kernel` inlines them into every callsite.
pub trait Kernel {
    /// Marks the semantic reference implementation ([`ScalarKernel`]). Debug builds
    /// cross-check every `apply`/`walk_scan` result of a non-reference kernel against the
    /// scalar one at the callsites — bit-equality there is what makes fit results
    /// independent of the host's SIMD level.
    #[cfg(debug_assertions)]
    const IS_REFERENCE: bool = false;

    /// The bound reduction. `CHECKED` enables the collapse checkpoints (see
    /// `COLLAPSE_CHECK_START`); without them the reduction is the exact straight-line
    /// fold over the full slices. Callers pick per cell via the `Estimates::collapsed`
    /// scheduling flag.
    fn apply<const CHECKED: bool>(row: &[f32], col: &[f32]) -> (f32, f32);

    /// The K-M walk of `update_value`: fold all items whose covariate (low 31 bits of
    /// `xs[i]`) lies in `[start, start + span]` into `(raw_value, remaining_weight)`;
    /// the high bit of `xs[i]` marks an observed (uncensored) item, `ws[i]` its weight.
    /// Items must be folded in slice order — the K-M product is order-sensitive.
    #[inline(always)]
    fn walk_scan(
        xs: &[u32],
        ws: &[f32],
        start: u32,
        span: u32,
        raw_value: f32,
        remaining_weight: f32,
    ) -> (f32, f32) {
        walk_scan_scalar(xs, ws, start, span, raw_value, remaining_weight)
    }
}

/// Scalar walk: one wrapping-subtraction range test per item (as unsigned,
/// `x - start > span` exactly when `x` is outside `[start, start + span]`); weights are
/// only read for items that pass.
#[inline]
fn walk_scan_scalar(
    xs: &[u32],
    ws: &[f32],
    start: u32,
    span: u32,
    mut raw_value: f32,
    mut remaining_weight: f32,
) -> (f32, f32) {
    debug_assert_eq!(xs.len(), ws.len());
    for (i, &x_flagged) in xs.iter().enumerate() {
        let x = x_flagged & !WALK_OBSERVED_BIT;
        if x.wrapping_sub(start) > span {
            continue;
        }
        let weight = ws[i];
        if x_flagged & WALK_OBSERVED_BIT != 0 {
            // The at-risk weight must still be strictly positive at every event fold:
            // `update_value`'s `weight_noise_floor` guard is what keeps walks away from
            // drifted-to-zero denominators, and a violation here means it under-guarded
            // (the product would blow up to inf/NaN). The factor itself may dip a few
            // ulps below 0 for a cell's last at-risk observation — that is expected
            // subtraction drift, not a bug, so only positivity is asserted.
            debug_assert!(remaining_weight > 0.0);
            raw_value *= 1.0 - weight / remaining_weight;
        }
        remaining_weight -= weight;
    }
    (raw_value, remaining_weight)
}

/// Portable scalar / SLP-vectorized fallback. Always usable.
pub struct ScalarKernel;

impl Kernel for ScalarKernel {
    #[cfg(debug_assertions)]
    const IS_REFERENCE: bool = true;

    #[inline(always)]
    fn apply<const CHECKED: bool>(row: &[f32], col: &[f32]) -> (f32, f32) {
        propagate_bounds_kernel::<CHECKED>(row, col)
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

/// Inner kernel of `propagate_bounds`: the reduction
/// `(max_k min(row[k], col[k]), min_k max(row[k], col[k]))`, exiting at the shared
/// collapse checkpoints (see `COLLAPSE_CHECK_START`).
///
/// Targets baseline SIMD only — SSE2 on x86_64, NEON on aarch64. The body keeps four
/// independent `(lower, upper)` accumulator lanes; with contiguous `[f32; 4]` reads SLP fuses
/// each per-lane group into one 128-bit `vminps`/`vmaxps` (or NEON `fmin`/`fmax`), and the
/// four lanes only meet in the final horizontal reduction, keeping dependency chains short.
const LANES: usize = 4;

#[inline]
fn propagate_bounds_kernel<const CHECKED: bool>(row: &[f32], col: &[f32]) -> (f32, f32) {
    debug_assert_eq!(row.len(), col.len());
    let n = row.len();

    let mut lo = [f32::NEG_INFINITY; LANES];
    let mut hi = [f32::INFINITY; LANES];
    let mut k = 0;
    let mut next_check = COLLAPSE_CHECK_START;

    #[inline(always)]
    fn reduce(lo: &[f32; LANES], hi: &[f32; LANES]) -> (f32, f32) {
        let lower = fast_max(fast_max(lo[0], lo[1]), fast_max(lo[2], lo[3]));
        let upper = fast_min(fast_min(hi[0], hi[1]), fast_min(hi[2], hi[3]));
        (lower, upper)
    }

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
        if CHECKED && k == next_check {
            let (lower, upper) = reduce(&lo, &hi);
            if lower >= upper {
                return (lower, upper);
            }
            next_check *= 2;
        }
    }
    while k < n {
        let l = fast_min(row[k], col[k]);
        let u = fast_max(row[k], col[k]);
        lo[0] = fast_max(lo[0], l);
        hi[0] = fast_min(hi[0], u);
        k += 1;
    }

    reduce(&lo, &hi)
}

/// AVX2 kernel — bare `vminps`/`vmaxps` on `ymm` (8-wide f32). Only safe to use after
/// `is_x86_feature_detected!("avx2")` confirmed it is available.
#[cfg(target_arch = "x86_64")]
pub struct Avx2Kernel;

#[cfg(target_arch = "x86_64")]
impl Kernel for Avx2Kernel {
    #[inline(always)]
    fn apply<const CHECKED: bool>(row: &[f32], col: &[f32]) -> (f32, f32) {
        // SAFETY: only monomorphized into the algorithm tree by `dispatch_generalized_pava`
        // after `is_x86_feature_detected!("avx2")` returned true.
        unsafe { propagate_bounds_kernel_avx2::<CHECKED>(row, col) }
    }

    #[inline(always)]
    fn walk_scan(
        xs: &[u32],
        ws: &[f32],
        start: u32,
        span: u32,
        raw_value: f32,
        remaining_weight: f32,
    ) -> (f32, f32) {
        // SAFETY: see `apply`.
        unsafe { walk_scan_avx2(xs, ws, start, span, raw_value, remaining_weight) }
    }
}

/// AVX2 walk: test 8 covariates per compare and 32 per round (four compare-masks merged
/// into one u32, one test+branch for the common no-hit case), then fold only the set mask
/// bits in index order — the K-M product is order-sensitive. The walk is usually
/// filter-dominated (a few percent of items pass), so vectorizing the range test removes
/// most of its cost.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx,avx2")]
unsafe fn walk_scan_avx2(
    xs: &[u32],
    ws: &[f32],
    start: u32,
    span: u32,
    mut raw_value: f32,
    mut remaining_weight: f32,
) -> (f32, f32) {
    use core::arch::x86_64::{
        _mm256_and_si256, _mm256_castsi256_ps, _mm256_cmpeq_epi32, _mm256_cmpgt_epi32,
        _mm256_loadu_si256, _mm256_maskload_epi32, _mm256_min_epu32, _mm256_movemask_ps,
        _mm256_set1_epi32, _mm256_setr_epi32, _mm256_sub_epi32,
    };

    debug_assert_eq!(xs.len(), ws.len());
    let n = xs.len();
    let mut k = 0;

    // SAFETY: AVX/AVX2 enabled by function-level `target_feature`; loads stay in-bounds
    // because each block tests the remaining length first, and mask-bit indices stay below
    // the block length.
    unsafe {
        let start_v = _mm256_set1_epi32(start as i32);
        let span_v = _mm256_set1_epi32(span as i32);
        let clear_v = _mm256_set1_epi32(!WALK_OBSERVED_BIT as i32);

        #[inline(always)]
        unsafe fn range_mask_8(
            xs: *const u32,
            clear_v: core::arch::x86_64::__m256i,
            start_v: core::arch::x86_64::__m256i,
            span_v: core::arch::x86_64::__m256i,
        ) -> u32 {
            // SAFETY: caller guarantees 8 readable lanes at `xs`.
            unsafe {
                let v = _mm256_loadu_si256(xs.cast());
                let x = _mm256_and_si256(v, clear_v);
                // Unsigned `x - start <= span` via min: t <= span  <=>  min(t, span) == t.
                let t = _mm256_sub_epi32(x, start_v);
                let in_range = _mm256_cmpeq_epi32(_mm256_min_epu32(t, span_v), t);
                _mm256_movemask_ps(_mm256_castsi256_ps(in_range)) as u32
            }
        }

        while k + 32 <= n {
            let p = xs.as_ptr().add(k);
            let m0 = range_mask_8(p, clear_v, start_v, span_v);
            let m1 = range_mask_8(p.add(8), clear_v, start_v, span_v);
            let m2 = range_mask_8(p.add(16), clear_v, start_v, span_v);
            let m3 = range_mask_8(p.add(24), clear_v, start_v, span_v);
            let mut mask = m0 | (m1 << 8) | (m2 << 16) | (m3 << 24);
            while mask != 0 {
                let j = mask.trailing_zeros() as usize;
                let x_flagged = *xs.get_unchecked(k + j);
                let weight = *ws.get_unchecked(k + j);
                if x_flagged & WALK_OBSERVED_BIT != 0 {
                    raw_value *= 1.0 - weight / remaining_weight;
                }
                remaining_weight -= weight;
                mask &= mask - 1;
            }
            k += 32;
        }

        while k + 8 <= n {
            let mut mask = range_mask_8(xs.as_ptr().add(k), clear_v, start_v, span_v);
            while mask != 0 {
                let j = mask.trailing_zeros() as usize;
                let x_flagged = *xs.get_unchecked(k + j);
                let weight = *ws.get_unchecked(k + j);
                if x_flagged & WALK_OBSERVED_BIT != 0 {
                    raw_value *= 1.0 - weight / remaining_weight;
                }
                remaining_weight -= weight;
                mask &= mask - 1;
            }
            k += 8;
        }

        // Masked tail (4-7 items): one masked load + compare instead of a scalar
        // per-item loop — the tail runs on every walk, and short deltas are common. The
        // lane mask must also gate the range mask: a zeroed inactive lane would pass the
        // range test whenever `start == 0`. Tails of 1-3 items stay scalar: cells that
        // are revisited on every arrival (large n over few covariate levels) see
        // 1-2-item deltas by the million, and there the masked setup costs more than
        // the items.
        if n - k >= 4 {
            let lane_idx = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
            let lmask = _mm256_cmpgt_epi32(_mm256_set1_epi32((n - k) as i32), lane_idx);
            let v = _mm256_maskload_epi32(xs.as_ptr().add(k).cast(), lmask);
            let x = _mm256_and_si256(v, clear_v);
            let t = _mm256_sub_epi32(x, start_v);
            let in_range =
                _mm256_and_si256(_mm256_cmpeq_epi32(_mm256_min_epu32(t, span_v), t), lmask);
            let mut mask = _mm256_movemask_ps(_mm256_castsi256_ps(in_range)) as u32;
            while mask != 0 {
                let j = mask.trailing_zeros() as usize;
                let x_flagged = *xs.get_unchecked(k + j);
                let weight = *ws.get_unchecked(k + j);
                if x_flagged & WALK_OBSERVED_BIT != 0 {
                    raw_value *= 1.0 - weight / remaining_weight;
                }
                remaining_weight -= weight;
                mask &= mask - 1;
            }
            k = n;
        }
    }

    walk_scan_scalar(&xs[k..], &ws[k..], start, span, raw_value, remaining_weight)
}

/// Bare-AVX2 inner kernel — 8-wide f32, two independent accumulator pairs for ILP, masked
/// (`vmaskmovps` + blend) tail so sub-8 remainders stay vectorized.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx,avx2")]
unsafe fn propagate_bounds_kernel_avx2<const CHECKED: bool>(
    row: &[f32],
    col: &[f32],
) -> (f32, f32) {
    use core::arch::x86_64::{
        __m256, _mm_cvtss_f32, _mm_max_ps, _mm_min_ps, _mm_movehdup_ps, _mm_movehl_ps,
        _mm256_blendv_ps, _mm256_castps256_ps128, _mm256_castsi256_ps, _mm256_cmpgt_epi32,
        _mm256_extractf128_ps, _mm256_loadu_ps, _mm256_maskload_ps, _mm256_max_ps, _mm256_min_ps,
        _mm256_set1_epi32, _mm256_set1_ps, _mm256_setr_epi32,
    };

    const F32_LANES: usize = 8;

    /// Horizontal reduction: 256 -> 128 -> 64 -> 32 bits via max/min + shuffles.
    #[inline]
    #[target_feature(enable = "avx,avx2")]
    unsafe fn reduce(lo: __m256, hi: __m256) -> (f32, f32) {
        let lo128 = _mm_max_ps(_mm256_castps256_ps128(lo), _mm256_extractf128_ps(lo, 1));
        let hi128 = _mm_min_ps(_mm256_castps256_ps128(hi), _mm256_extractf128_ps(hi, 1));
        let lo64 = _mm_max_ps(lo128, _mm_movehl_ps(lo128, lo128));
        let hi64 = _mm_min_ps(hi128, _mm_movehl_ps(hi128, hi128));
        let lo32 = _mm_max_ps(lo64, _mm_movehdup_ps(lo64));
        let hi32 = _mm_min_ps(hi64, _mm_movehdup_ps(hi64));
        (_mm_cvtss_f32(lo32), _mm_cvtss_f32(hi32))
    }

    debug_assert_eq!(row.len(), col.len());
    let n = row.len();
    let mut k = 0;
    let mut next_check = COLLAPSE_CHECK_START;

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
            // Collapse checkpoint (see `COLLAPSE_CHECK_START`): `k` steps by 16 from 0, so
            // it hits every multiple of 128 exactly.
            if CHECKED && k == next_check {
                let (lower, upper) = reduce(_mm256_max_ps(lo0, lo1), _mm256_min_ps(hi0, hi1));
                if lower >= upper {
                    return (lower, upper);
                }
                next_check *= 2;
            }
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

        reduce(_mm256_max_ps(lo0, lo1), _mm256_min_ps(hi0, hi1))
    }
}

/// AVX-512F kernel — bare `vminps`/`vmaxps` on `zmm` (16-wide f32). Only safe to use after
/// `is_x86_feature_detected!("avx512f")` confirmed it is available.
#[cfg(target_arch = "x86_64")]
pub struct Avx512Kernel;

#[cfg(target_arch = "x86_64")]
impl Kernel for Avx512Kernel {
    #[inline(always)]
    fn apply<const CHECKED: bool>(row: &[f32], col: &[f32]) -> (f32, f32) {
        // SAFETY: only monomorphized into the algorithm tree by `dispatch_generalized_pava`
        // after `is_x86_feature_detected!("avx512f")` returned true.
        unsafe { propagate_bounds_kernel_avx512::<CHECKED>(row, col) }
    }

    #[inline(always)]
    fn walk_scan(
        xs: &[u32],
        ws: &[f32],
        start: u32,
        span: u32,
        raw_value: f32,
        remaining_weight: f32,
    ) -> (f32, f32) {
        // SAFETY: see `apply`.
        unsafe { walk_scan_avx512(xs, ws, start, span, raw_value, remaining_weight) }
    }
}

/// AVX-512 walk: test 16 covariates per compare (native unsigned compare-to-mask) and 64
/// per round, then fold only the set mask bits in index order. See the AVX2 variant for
/// the rationale.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn walk_scan_avx512(
    xs: &[u32],
    ws: &[f32],
    start: u32,
    span: u32,
    mut raw_value: f32,
    mut remaining_weight: f32,
) -> (f32, f32) {
    use core::arch::x86_64::{
        _mm512_and_si512, _mm512_cmple_epu32_mask, _mm512_loadu_si512, _mm512_maskz_loadu_epi32,
        _mm512_set1_epi32, _mm512_sub_epi32,
    };

    debug_assert_eq!(xs.len(), ws.len());
    let n = xs.len();
    let mut k = 0;

    // SAFETY: AVX-512F enabled by function-level `target_feature`; loads stay in-bounds
    // because each block tests the remaining length first, and mask-bit indices stay below
    // the block length.
    unsafe {
        let start_v = _mm512_set1_epi32(start as i32);
        let span_v = _mm512_set1_epi32(span as i32);
        let clear_v = _mm512_set1_epi32(!WALK_OBSERVED_BIT as i32);

        #[inline(always)]
        unsafe fn range_mask_16(
            xs: *const u32,
            clear_v: core::arch::x86_64::__m512i,
            start_v: core::arch::x86_64::__m512i,
            span_v: core::arch::x86_64::__m512i,
        ) -> u16 {
            // SAFETY: caller guarantees 16 readable lanes at `xs`.
            unsafe {
                let v = _mm512_loadu_si512(xs.cast());
                let x = _mm512_and_si512(v, clear_v);
                let t = _mm512_sub_epi32(x, start_v);
                _mm512_cmple_epu32_mask(t, span_v)
            }
        }

        // 64 items per round: four compare-masks combined into one u64, one test+branch
        // for the (common) no-hit case, and hits processed in index order off the bits.
        while k + 64 <= n {
            let p = xs.as_ptr().add(k);
            let m0 = range_mask_16(p, clear_v, start_v, span_v) as u64;
            let m1 = range_mask_16(p.add(16), clear_v, start_v, span_v) as u64;
            let m2 = range_mask_16(p.add(32), clear_v, start_v, span_v) as u64;
            let m3 = range_mask_16(p.add(48), clear_v, start_v, span_v) as u64;
            let mut mask = m0 | (m1 << 16) | (m2 << 32) | (m3 << 48);
            while mask != 0 {
                let j = mask.trailing_zeros() as usize;
                let x_flagged = *xs.get_unchecked(k + j);
                let weight = *ws.get_unchecked(k + j);
                if x_flagged & WALK_OBSERVED_BIT != 0 {
                    raw_value *= 1.0 - weight / remaining_weight;
                }
                remaining_weight -= weight;
                mask &= mask - 1;
            }
            k += 64;
        }

        while k + 16 <= n {
            let mut mask = range_mask_16(xs.as_ptr().add(k), clear_v, start_v, span_v) as u32;
            while mask != 0 {
                let j = mask.trailing_zeros() as usize;
                let x_flagged = *xs.get_unchecked(k + j);
                let weight = *ws.get_unchecked(k + j);
                if x_flagged & WALK_OBSERVED_BIT != 0 {
                    raw_value *= 1.0 - weight / remaining_weight;
                }
                remaining_weight -= weight;
                mask &= mask - 1;
            }
            k += 16;
        }

        // Masked tail (4-15 items): one masked load + compare instead of a scalar
        // per-item loop — the tail runs on every walk, and short deltas are common. The
        // load mask must also gate the range mask: a zeroed inactive lane would pass the
        // range test whenever `start == 0`. Tails of 1-3 items stay scalar: cells that
        // are revisited on every arrival (large n over few covariate levels) see
        // 1-2-item deltas by the million, and there the masked setup costs more than
        // the items.
        if n - k >= 4 {
            let rem_mask = ((1u32 << (n - k)) - 1) as u16;
            let v = _mm512_maskz_loadu_epi32(rem_mask, xs.as_ptr().add(k).cast());
            let x = _mm512_and_si512(v, clear_v);
            let t = _mm512_sub_epi32(x, start_v);
            let mut mask = (_mm512_cmple_epu32_mask(t, span_v) & rem_mask) as u32;
            while mask != 0 {
                let j = mask.trailing_zeros() as usize;
                let x_flagged = *xs.get_unchecked(k + j);
                let weight = *ws.get_unchecked(k + j);
                if x_flagged & WALK_OBSERVED_BIT != 0 {
                    raw_value *= 1.0 - weight / remaining_weight;
                }
                remaining_weight -= weight;
                mask &= mask - 1;
            }
            k = n;
        }
    }

    walk_scan_scalar(&xs[k..], &ws[k..], start, span, raw_value, remaining_weight)
}

/// AVX-512F inner kernel — 16-wide f32, two independent accumulator pairs for ILP, native
/// masked-load tail so sub-16 remainders stay vectorized.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn propagate_bounds_kernel_avx512<const CHECKED: bool>(
    row: &[f32],
    col: &[f32],
) -> (f32, f32) {
    use core::arch::x86_64::{
        __m512, _mm_cvtss_f32, _mm_max_ps, _mm_min_ps, _mm_movehdup_ps, _mm_movehl_ps,
        _mm512_castps512_ps128, _mm512_loadu_ps, _mm512_mask_loadu_ps, _mm512_max_ps,
        _mm512_min_ps, _mm512_set1_ps, _mm512_shuffle_f32x4,
    };

    const F32_LANES: usize = 16;

    /// Horizontal reduction: fold the four 128-bit blocks together (AVX-512F-only —
    /// no DQ extracts), then reduce within 128 bits.
    #[inline]
    #[target_feature(enable = "avx512f")]
    unsafe fn reduce(lo: __m512, hi: __m512) -> (f32, f32) {
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

    debug_assert_eq!(row.len(), col.len());
    let n = row.len();
    let mut k = 0;
    let mut next_check = COLLAPSE_CHECK_START;

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
            // Collapse checkpoint (see `COLLAPSE_CHECK_START`): `k` steps by 32 from 0, so
            // it hits every multiple of 128 exactly.
            if CHECKED && k == next_check {
                let (lower, upper) = reduce(_mm512_max_ps(lo0, lo1), _mm512_min_ps(hi0, hi1));
                if lower >= upper {
                    return (lower, upper);
                }
                next_check *= 2;
            }
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

        reduce(_mm512_max_ps(lo0, lo1), _mm512_min_ps(hi0, hi1))
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

    /// Reference semantics of the kernels: the exact prefix reduction with an exit at the
    /// shared element-count checkpoints (128, 256, 512, …).
    fn reference_reduction(row: &[f32], col: &[f32]) -> (f32, f32) {
        let mut expected = (f32::NEG_INFINITY, f32::INFINITY);
        let mut next_check = COLLAPSE_CHECK_START;
        for k in 0..row.len() {
            let l = if row[k] < col[k] { row[k] } else { col[k] };
            let u = if row[k] > col[k] { row[k] } else { col[k] };
            expected.0 = if expected.0 > l { expected.0 } else { l };
            expected.1 = if expected.1 < u { expected.1 } else { u };
            if k + 1 == next_check && expected.0 >= expected.1 {
                return expected;
            }
            if k + 1 == next_check {
                next_check *= 2;
            }
        }
        expected
    }

    /// All kernels implement the same reduction with the same collapse-checkpoint
    /// schedule, so they must agree bit-for-bit with the naive definition (and therefore
    /// with each other) on every length — empty, tail-only, and multi-checkpoint slices,
    /// collapsing and non-collapsing inputs alike. This is what guarantees fit results do
    /// not depend on the host's SIMD level.
    #[test]
    fn kernels_agree() {
        let mut state = 0x5eed;
        // Random values in (0, 1] collapse the running bounds almost immediately, so long
        // slices exercise the checkpoint exits; the interleaved non-collapsing case keeps
        // `min(row, col) < 0.5 < max(row, col)` everywhere, so it always reduces in full.
        let lengths = (0..=257).chain([300, 511, 512, 513, 800, 1024, 1500]);
        for n in lengths {
            for collapsing in [true, false] {
                let (row, col): (Vec<f32>, Vec<f32>) = if collapsing {
                    (
                        (0..n).map(|_| lcg_f32(&mut state)).collect(),
                        (0..n).map(|_| lcg_f32(&mut state)).collect(),
                    )
                } else {
                    (
                        (0..n).map(|_| lcg_f32(&mut state) * 0.4).collect(),
                        (0..n).map(|_| 0.6 + lcg_f32(&mut state) * 0.4).collect(),
                    )
                };

                let expected_checked = reference_reduction(&row, &col);
                // The unchecked mode is the exact full reduction: the reference with the
                // checkpoints never firing.
                let mut expected_full = (f32::NEG_INFINITY, f32::INFINITY);
                for k in 0..row.len() {
                    let l = if row[k] < col[k] { row[k] } else { col[k] };
                    let u = if row[k] > col[k] { row[k] } else { col[k] };
                    expected_full.0 = if expected_full.0 > l {
                        expected_full.0
                    } else {
                        l
                    };
                    expected_full.1 = if expected_full.1 < u {
                        expected_full.1
                    } else {
                        u
                    };
                }

                assert_eq!(
                    ScalarKernel::apply::<true>(&row, &col),
                    expected_checked,
                    "scalar checked, n={n} collapsing={collapsing}"
                );
                assert_eq!(
                    ScalarKernel::apply::<false>(&row, &col),
                    expected_full,
                    "scalar full, n={n} collapsing={collapsing}"
                );
                #[cfg(target_arch = "x86_64")]
                {
                    if is_x86_feature_detected!("avx2") {
                        assert_eq!(
                            Avx2Kernel::apply::<true>(&row, &col),
                            expected_checked,
                            "avx2 checked, n={n} collapsing={collapsing}"
                        );
                        assert_eq!(
                            Avx2Kernel::apply::<false>(&row, &col),
                            expected_full,
                            "avx2 full, n={n} collapsing={collapsing}"
                        );
                    }
                    if is_x86_feature_detected!("avx512f") {
                        assert_eq!(
                            Avx512Kernel::apply::<true>(&row, &col),
                            expected_checked,
                            "avx512 checked, n={n} collapsing={collapsing}"
                        );
                        assert_eq!(
                            Avx512Kernel::apply::<false>(&row, &col),
                            expected_full,
                            "avx512 full, n={n} collapsing={collapsing}"
                        );
                    }
                }
            }
        }
    }

    /// The SIMD walk scans must fold exactly the same items in the same order as the scalar
    /// walk, so all implementations return bit-identical (raw_value, remaining_weight).
    #[test]
    fn walk_scans_agree() {
        let mut state = 0xa11c;
        for n in 0..=130 {
            // Covariate indices in a small range so the filter passes a meaningful share;
            // ~half the items observed.
            let xs: Vec<u32> = (0..n)
                .map(|_| {
                    let x = (lcg_f32(&mut state) * 40.0) as u32;
                    let observed = lcg_f32(&mut state) > 0.5;
                    x | if observed { WALK_OBSERVED_BIT } else { 0 }
                })
                .collect();
            let ws: Vec<f32> = (0..n).map(|_| lcg_f32(&mut state) + 0.5).collect();

            for (start, span) in [(0u32, 40u32), (10, 5), (35, 0), (12, 13)] {
                let expected = ScalarKernel::walk_scan(&xs, &ws, start, span, 1.0, n as f32 + 50.0);
                #[cfg(target_arch = "x86_64")]
                {
                    if is_x86_feature_detected!("avx2") {
                        assert_eq!(
                            Avx2Kernel::walk_scan(&xs, &ws, start, span, 1.0, n as f32 + 50.0),
                            expected,
                            "avx2, n={n} start={start} span={span}",
                        );
                    }
                    if is_x86_feature_detected!("avx512f") {
                        assert_eq!(
                            Avx512Kernel::walk_scan(&xs, &ws, start, span, 1.0, n as f32 + 50.0),
                            expected,
                            "avx512, n={n} start={start} span={span}",
                        );
                    }
                }
            }
        }
    }
}
