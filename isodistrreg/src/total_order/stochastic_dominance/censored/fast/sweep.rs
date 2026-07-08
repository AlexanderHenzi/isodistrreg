//! The tile sweep over one merge rectangle, and the dispatch seam that runs it inside a
//! frame with the kernel's `#[target_feature]` set enabled.
//!
//! `K::apply` and `K::walk_scan` are `#[target_feature]` functions on the AVX kernels, and
//! such functions only inline into callers that enable at least the same features. The
//! sweep loop itself is compiled for the portable baseline, so calling the kernels
//! straight from it would leave two non-inlinable calls per cell. [`SweepDispatch`]
//! therefore enters one kernel-matching `#[target_feature]` frame per merge rectangle;
//! [`sweep_rect_body`] is `#[inline(always)]`, so the whole cell loop lands inside that
//! frame and the kernel calls inline into it.

#[cfg(target_arch = "x86_64")]
use crate::total_order::stochastic_dominance::censored::propagate_bounds::{
    Avx2Kernel, Avx512Kernel,
};
use crate::total_order::stochastic_dominance::censored::propagate_bounds::{Kernel, ScalarKernel};
use crate::total_order::stochastic_dominance::censored::structures::{CompletionIndex, Estimates};
use crate::total_order::structures::CensoredContext;

/// Number of columns processed per tile in the sweep's cell-update loops. Sized so a
/// tile's columns plus their bounds/K-M state stay cache-resident across the r sweep on
/// any modern CPU (>= 512 KB L2); 64-128 measured equivalently well, smaller tiles start
/// paying per-tile overhead.
const POOL_TILE: usize = 96;

/// One merge rectangle's sweep, bundled so the dispatch trait, the `#[target_feature]`
/// wrappers, and the body all share a single signature.
pub(super) struct SweepArgs<'a, X: crate::Float, Y: crate::Float> {
    pub(super) estimates: &'a mut Estimates,
    pub(super) data_index: usize,
    /// Rows `[penultimate_start, row_end)` — already clamped by the caller's
    /// within-remainder row skip.
    pub(super) penultimate_start: usize,
    pub(super) row_end: usize,
    /// Columns `[ultimate_start, ultimate_end)`.
    pub(super) ultimate_start: usize,
    pub(super) ultimate_end: usize,
    pub(super) input: &'a CensoredContext<X, Y>,
    pub(super) epsilon: f32,
    /// Prefix maxima of the completion markers over the ultimate range:
    /// `marker_buf[s - ultimate_start]` = max of `m[ultimate_start..=s]`.
    pub(super) marker_buf: &'a [u32],
    /// Suffix maxima of the completion markers over the penultimate range:
    /// `left_marker_buf[r - penultimate_start]` = max of `m[r..penultimate_end]`.
    pub(super) left_marker_buf: &'a [u32],
}

/// Dispatch seam that runs one merge rectangle's whole tile sweep inside a frame with the
/// kernel's `#[target_feature]` set enabled, so `K::apply` and `K::walk_scan` — otherwise
/// non-inlinable `#[target_feature]` functions called twice per cell — inline into the
/// cell loop.
pub(super) trait SweepDispatch: Kernel + Sized {
    fn sweep_rect<X: crate::Float, Y: crate::Float>(args: SweepArgs<'_, X, Y>) {
        sweep_rect_body::<Self, X, Y>(args);
    }
}

impl SweepDispatch for ScalarKernel {}

/// Implement [`SweepDispatch`] for an AVX kernel: the trait method forwards into a
/// `#[target_feature]` wrapper whose only job is to give [`sweep_rect_body`] a frame the
/// kernel's functions can inline into.
#[cfg(target_arch = "x86_64")]
macro_rules! target_feature_sweep {
    ($kernel:ty, $features:literal, $wrapper:ident) => {
        impl SweepDispatch for $kernel {
            #[inline(always)]
            fn sweep_rect<X: crate::Float, Y: crate::Float>(args: SweepArgs<'_, X, Y>) {
                // SAFETY: the kernel is only monomorphized into the algorithm tree by
                // `dispatch_generalized_pava` after the matching
                // `is_x86_feature_detected!` returned true.
                unsafe { $wrapper(args) }
            }
        }

        #[target_feature(enable = $features)]
        unsafe fn $wrapper<X: crate::Float, Y: crate::Float>(args: SweepArgs<'_, X, Y>) {
            sweep_rect_body::<$kernel, X, Y>(args);
        }
    };
}

#[cfg(target_arch = "x86_64")]
target_feature_sweep!(Avx2Kernel, "avx,avx2", sweep_rect_avx2);
#[cfg(target_arch = "x86_64")]
target_feature_sweep!(Avx512Kernel, "avx512f", sweep_rect_avx512);

/// The tile sweep over one merge rectangle: rows `[penultimate_start, row_end)`, columns
/// `[ultimate_start, ultimate_end)`.
///
/// Iterating over r in the outer loop and over s in the inner loop is significantly
/// faster than the transposed order; on top of that the s range is processed in
/// tiles of `POOL_TILE` columns (tile outer, r descending inside, s ascending
/// within the tile). Cell (r, s) reads row entries (r, k<s) — written in earlier
/// tiles or earlier in this row's sweep — and column entries (k>r, s) — written at
/// larger r within the same tile — so the tiled order is a valid schedule of the
/// same DP and every cell sees identical inputs (bit-identical results). The point
/// of the tiles: wide merges revisit each column once per row, and with the full
/// span in flight the columns (plus their K-M state) fall out of L2 between visits;
/// a tile's working set stays cache-resident across the r sweep.
#[inline(always)]
fn sweep_rect_body<K: Kernel, X: crate::Float, Y: crate::Float>(args: SweepArgs<'_, X, Y>) {
    let SweepArgs {
        estimates,
        data_index,
        penultimate_start,
        row_end,
        ultimate_start,
        ultimate_end,
        input,
        epsilon,
        marker_buf,
        left_marker_buf,
    } = args;

    // One cell: kernel, then K-M update. Index and bounds flow through registers between
    // the two — bounds are never stored (only the value store remains for future readers),
    // so the consecutive-cell dependency chain skips the memory round-trip.
    #[inline(always)]
    fn cell<K: Kernel, X: crate::Float, Y: crate::Float>(
        estimates: &mut Estimates,
        data_index: usize,
        r: usize,
        s: usize,
        input: &CensoredContext<X, Y>,
        epsilon: f32,
        marker_max: u32,
    ) {
        let idx = Estimates::compute_index((r, s), estimates.len());
        // Completed cells pin to exactly 0 before `update_value` ever looks at the
        // bounds (its first branch), so the O(s - r) kernel reduction below would be
        // computed only to be discarded. Do the same stores here and skip it. The
        // check re-fires on every visit once true (monotone in `data_index`), so
        // this shortcut applies to every revisit of a completed cell.
        if CompletionIndex::completes_marker(marker_max, data_index) {
            let row_idx = Estimates::compute_row_index((r, s), estimates.len());
            estimates.cold[idx].raw_value = 0.0;
            estimates.set_value(idx, row_idx, 0.0);
            return;
        }
        let bounds = estimates.propagate_bounds_with_row::<K>(idx, r, s);
        estimates
            .update_value::<K, _, _>(data_index, idx, bounds, r, s, input, epsilon, marker_max);
    }

    let mut tile_start = ultimate_start;
    while tile_start < ultimate_end {
        let tile_end = (tile_start + POOL_TILE).min(ultimate_end);
        // Rows advance in pairs per column: cell (r0, s) is computed right before
        // (r0 - 1, s), whose column operand merely *ends* with the fresh (r0, s) value —
        // so the two kernels overlap almost entirely, doubling the independent work in
        // flight. The pair schedule is a valid order of the same DP (each cell still
        // sees row entries with smaller s and column entries with larger r already
        // computed), so results stay bit-identical.
        let mut r_high = row_end;
        while r_high >= penultimate_start + 2 {
            let r0 = r_high - 1;
            let r1 = r_high - 2;
            let left0 = left_marker_buf[r0 - penultimate_start];
            let left1 = left_marker_buf[r1 - penultimate_start];
            for s in tile_start..tile_end {
                let right = marker_buf[s - ultimate_start];
                cell::<K, _, _>(
                    estimates,
                    data_index,
                    r0,
                    s,
                    input,
                    epsilon,
                    left0.max(right),
                );
                cell::<K, _, _>(
                    estimates,
                    data_index,
                    r1,
                    s,
                    input,
                    epsilon,
                    left1.max(right),
                );
            }
            r_high -= 2;
        }
        if r_high > penultimate_start {
            let left = left_marker_buf[0];
            for s in tile_start..tile_end {
                let right = marker_buf[s - ultimate_start];
                cell::<K, _, _>(
                    estimates,
                    data_index,
                    penultimate_start,
                    s,
                    input,
                    epsilon,
                    left.max(right),
                );
            }
        }
        tile_start = tile_end;
    }
}
