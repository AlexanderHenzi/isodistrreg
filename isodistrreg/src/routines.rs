use crate::structures::{Direction, Increasing, Observation};
use std::cmp::Ordering;
use std::iter::once;
use std::mem;

/// Inner tile size for the blocked transposes. A `BLOCK × BLOCK` tile of `f32` is 9 KiB; a
/// swap holds two such tiles of scratch (18 KiB), which stays resident in the smallest
/// modern L1d (32 KiB) with headroom for the matrix rows streaming through.
///
/// 48 is chosen empirically: each `copy_from_slice` then moves a 192-byte (BLOCK·4) run,
/// long enough for the hardware prefetcher to engage on the otherwise row-scattered tiles
/// (a transpose of a multi-GB matrix is bound by scattered DRAM access). Smaller tiles
/// (32: 128-byte runs) leave throughput on the table; larger ones (64: two 16 KiB buffers
/// = a full 32 KiB L1) thrash the smallest laptop caches. The value is deliberately a
/// conservative middle, not tuned to one machine.
const BLOCK: usize = 48;

/// Outer block size for the two-level square transpose. Inner `BLOCK` tiles keep each swap
/// in L1, but on a matrix far larger than cache every tile's rows still come from distinct
/// pages, so the inner pass alone is throttled by TLB/last-level-cache misses. Grouping the
/// inner tiles into `OUTER_BLOCK × OUTER_BLOCK` regions and finishing one region (and its
/// mirror) before moving on keeps that working set — 2·`768²`·4 B ≈ 4.5 MiB for a pair —
/// resident in the last-level cache of any recent laptop, roughly halving the run time at
/// the multi-GB sizes that dominate. 768 is a deliberately conservative, machine-agnostic
/// middle: large enough to amortise the page/line misses, small enough to fit a modest L3.
const OUTER_BLOCK: usize = 768;

/// Transposes a row-major `m × n` matrix stored in `matrix`.
///
/// Both paths move each tile through a small stack buffer so every *main-memory* access
/// is contiguous (one row segment per `copy_from_slice`); the unavoidable strided access
/// of a transpose is confined to the L1-resident buffer. A naive element-wise swap streams
/// one operand down a column — a fresh cache line and TLB entry per element — which on the
/// multi-GB matrices that dominate the uncensored fit runs at a small fraction of DRAM
/// bandwidth.
///
/// The square (`m == n`) case is done in place with no extra allocation — at the sizes
/// that dominate (uncensored continuous, where `m == n` and the matrix is multiple GB),
/// a second buffer would double peak RSS. Rectangular matrices transpose into a fresh
/// buffer.
pub fn transpose<T: Copy>(matrix: &mut Vec<T>, m: usize, n: usize) {
    assert_eq!(matrix.len(), m * n);
    if m == 0 || n == 0 {
        return;
    }

    // `T: Copy` and the matrix is non-empty, so this element seeds the stack scratch
    // without needing `T: Default` or `MaybeUninit`.
    let fill = matrix[0];

    if m == n {
        transpose_square_inplace(matrix, n, fill);
        return;
    }

    // Rectangular: write into a fresh `n × m` buffer. `vec![fill; _]` memsets once (a
    // write-only pass) rather than `clone()`'s read+write of the whole matrix; every
    // element is overwritten below.
    let mut out: Vec<T> = vec![fill; m * n];
    let mut buf = [fill; BLOCK * BLOCK];
    // Output is `n × m`; tile over output rows (`ii`, source columns) and output
    // columns (`jj`, source rows).
    let mut ii = 0;
    while ii < n {
        let ie = (ii + BLOCK).min(n);
        let rci = ie - ii;
        let mut jj = 0;
        while jj < m {
            let je = (jj + BLOCK).min(m);
            let rcj = je - jj;
            // Read source tile (rows jj..je, cols ii..ie) contiguously into `buf`.
            for bj in 0..rcj {
                let base = (jj + bj) * n + ii;
                buf[bj * rci..bj * rci + rci].copy_from_slice(&matrix[base..base + rci]);
            }
            // Write its transpose into the output tile (rows ii..ie, cols jj..je),
            // one contiguous output row at a time.
            for bi in 0..rci {
                let base = (ii + bi) * m + jj;
                let dst = &mut out[base..base + rcj];
                for (bj, d) in dst.iter_mut().enumerate() {
                    *d = buf[bj * rci + bi];
                }
            }
            jj += BLOCK;
        }
        ii += BLOCK;
    }
    *matrix = out;
}

/// Transpose the diagonal `BLOCK` tile `[s, e) × [s, e)` of an `n × n` matrix in place.
fn transpose_diagonal_tile<T: Copy>(matrix: &mut [T], n: usize, s: usize, e: usize) {
    for i in s..e {
        for j in (i + 1)..e {
            matrix.swap(i * n + j, j * n + i);
        }
    }
}

/// Swap the off-diagonal `BLOCK` tile pair A = (rows `ii..ie`, cols `jj..je`) and
/// B = (rows `jj..je`, cols `ii..ie`) of an `n × n` matrix, leaving A = Bᵀ and B = Aᵀ.
/// Both tiles are staged through row-major stack buffers so every matrix read and write is
/// a contiguous row segment; the strided access of the transpose stays in the L1 buffers.
fn swap_offdiagonal_tiles<T: Copy>(
    matrix: &mut [T],
    n: usize,
    (ii, ie): (usize, usize),
    (jj, je): (usize, usize),
    buf_a: &mut [T],
    buf_b: &mut [T],
) {
    let ra = ie - ii;
    let ca = je - jj;
    for r in 0..ra {
        let base = (ii + r) * n + jj;
        buf_a[r * ca..r * ca + ca].copy_from_slice(&matrix[base..base + ca]);
    }
    for r in 0..ca {
        let base = (jj + r) * n + ii;
        buf_b[r * ra..r * ra + ra].copy_from_slice(&matrix[base..base + ra]);
    }
    // A[r][c] = Bᵀ[r][c] = buf_b[c][r]; write one contiguous row of A at a time.
    for r in 0..ra {
        let base = (ii + r) * n + jj;
        let dst = &mut matrix[base..base + ca];
        for (c, d) in dst.iter_mut().enumerate() {
            *d = buf_b[c * ra + r];
        }
    }
    // B[r][c] = Aᵀ[r][c] = buf_a[c][r].
    for r in 0..ca {
        let base = (jj + r) * n + ii;
        let dst = &mut matrix[base..base + ra];
        for (c, d) in dst.iter_mut().enumerate() {
            *d = buf_a[c * ca + r];
        }
    }
}

/// In-place transpose of an `n × n` row-major matrix, two-level blocked: `OUTER_BLOCK`
/// regions (last-level-cache locality) of `BLOCK` tiles (L1 locality). Each outer region
/// — a diagonal one, then it paired with every region to its right — is finished before
/// the next, so its pages stay cache-resident across all the inner tiles it contains.
/// Within a region, diagonal tiles transpose in place and off-diagonal tile pairs swap
/// through stack buffers (all matrix access contiguous).
fn transpose_square_inplace<T: Copy>(matrix: &mut [T], n: usize, fill: T) {
    let mut buf_a = [fill; BLOCK * BLOCK];
    let mut buf_b = [fill; BLOCK * BLOCK];

    // Inner tiles of the diagonal outer region `[o, oe) × [o, oe)`.
    let inner_diagonal_region =
        |matrix: &mut [T], o: usize, oe: usize, buf_a: &mut [T], buf_b: &mut [T]| {
            let mut ii = o;
            while ii < oe {
                let ie = (ii + BLOCK).min(oe);
                transpose_diagonal_tile(matrix, n, ii, ie);
                let mut jj = ie;
                while jj < oe {
                    let je = (jj + BLOCK).min(oe);
                    swap_offdiagonal_tiles(matrix, n, (ii, ie), (jj, je), buf_a, buf_b);
                    jj += BLOCK;
                }
                ii += BLOCK;
            }
        };

    let mut oi = 0;
    while oi < n {
        let oie = (oi + OUTER_BLOCK).min(n);
        inner_diagonal_region(matrix, oi, oie, &mut buf_a, &mut buf_b);
        // Off-diagonal outer regions to the right: all inner tile pairs between the row
        // band `[oi, oie)` and the column band `[oj, oje)`.
        let mut oj = oie;
        while oj < n {
            let oje = (oj + OUTER_BLOCK).min(n);
            let mut ii = oi;
            while ii < oie {
                let ie = (ii + BLOCK).min(oie);
                let mut jj = oj;
                while jj < oje {
                    let je = (jj + BLOCK).min(oje);
                    swap_offdiagonal_tiles(matrix, n, (ii, ie), (jj, je), &mut buf_a, &mut buf_b);
                    jj += BLOCK;
                }
                ii += BLOCK;
            }
            oj += OUTER_BLOCK;
        }
        oi += OUTER_BLOCK;
    }
}

pub fn argsort_unstable_by<D: Direction, F>(cmp: F, n: usize) -> Vec<usize>
where
    F: FnMut(usize, usize) -> Ordering,
{
    argsort_indices_unstable_by::<D, F>(cmp, (0..n).collect())
}

/// [`argsort_unstable_by`] over a caller-built index vector, which need not be
/// contiguous — e.g. pre-filtered to positive-weight observations, so dropped elements
/// are never even compared.
pub fn argsort_indices_unstable_by<D: Direction, F>(mut cmp: F, mut idx: Vec<usize>) -> Vec<usize>
where
    F: FnMut(usize, usize) -> Ordering,
{
    idx.sort_unstable_by(|&i, &j| {
        let ord = cmp(i, j);
        if D::IS_INCREASING { ord } else { ord.reverse() }
    });
    idx
}

// Compute the empirical CDF. Assumes that observations have been deduplicated so responses are
// unique. The response slot type `Y` is not used arithmetically here (the cumulative share is
// f32), so any type fits — typically the caller's `Y: Float` threshold type.
pub fn empirical_cdf<C, Y, I: Into<Observation<C, Y, (), f32>>>(
    observations: impl Iterator<Item = I>,
    total_weight: f32,
) -> Vec<f32> {
    let mut cdf: Vec<f32> = observations
        .map(Into::into)
        .scan(0.0, |acc, o| {
            *acc += o.weight / total_weight;
            // Ensure numerical imprecision won't get us out of [0, 1]
            *acc = acc.clamp(0.0, 1.0);
            Some(*acc)
        })
        .collect();
    // The final cumulative value is total_weight / total_weight — exactly 1. Pin it so
    // downstream proper-CDF checks (e.g. `prediction::mean`) see 1.0 rather than
    // accumulated f32 round-off like 0.99999994.
    if let Some(last) = cdf.last_mut() {
        *last = 1.0;
    }
    cdf
}

/// Compute the Kaplan-Meier estimator on unique responses.
///
/// Responses may repeat (tied observations, or the same threshold seen for several
/// covariates); one CDF value is emitted per unique response, reflecting every
/// observation in the tied group. Within a tied group the stream must order deaths
/// before censorings, matching the standard at-risk convention — which is the order
/// all in-tree callers produce.
pub fn kaplan_meier<C, R: Copy + PartialEq>(
    observations: impl Iterator<Item = Observation<C, R, bool, f32>>,
    total_weight: f32,
) -> Vec<f32> {
    let mut cdfs = Vec::new();
    let mut survival: f32 = 1.0;
    let mut remaining_weight = total_weight;
    let mut previous = None;
    // Survival reaches exactly 0 iff the last positive-weight observation is observed
    // (its factor is then exactly `1 - w/w`) — a purely combinatorial condition. Track
    // it so the final value can be pinned: the f32 running weight subtraction drifts,
    // leaving the product a few ulps off where mathematics says exactly 0.
    let mut last_positive_observed = false;
    for o in observations {
        // A response's value is final only once its tied group is complete.
        match previous.replace(o.y) {
            Some(value) if value != o.y => cdfs.push(survival),
            _ => {}
        }
        if o.observed {
            survival *= 1.0 - o.weight / remaining_weight;
        }
        remaining_weight -= o.weight;
        if o.weight > 0.0 {
            last_positive_observed = o.observed;
        }
    }
    if previous.is_some() {
        if last_positive_observed {
            survival = 0.0;
        }
        cdfs.push(survival);
    }
    // Clamping for numerics
    for v in &mut cdfs {
        *v = 1.0 - v.clamp(0.0, 1.0);
    }
    cdfs
}

pub fn median(elements: impl IntoIterator<Item = f64>) -> f64 {
    let allocated: Vec<_> = elements.into_iter().collect();
    assert!(!allocated.len().is_multiple_of(2));
    let median_index = allocated.len() / 2;
    select(allocated, median_index)
}

#[must_use]
pub fn select(mut elements: Vec<f64>, index: usize) -> f64 {
    // Debug-only: fires once per recursion level; in release an out-of-range index still
    // panics on the slice indexing below rather than corrupting silently.
    debug_assert!(index < elements.len());

    match elements.len() {
        0 => unreachable!(),
        1 => elements[0],
        2 => select2(elements[0], elements[1], index),
        3 => select3(elements[0], elements[1], elements[2], index),
        4 => select4(elements[0], elements[1], elements[2], elements[3], index),
        _ => {
            // recurse

            let (chunks, remainder) = elements.as_chunks::<5>();

            let size_5_groups = chunks.iter().map(|g| median5(g[0], g[1], g[2], g[3], g[4]));
            let small_medians: Vec<_> = if remainder.is_empty() {
                size_5_groups.collect()
            } else {
                let index = remainder.len() / 2;
                let last_group_value = match *remainder {
                    [a] => a,
                    [a, b] => select2(a, b, index),
                    [a, b, c] => select3(a, b, c, index),
                    [a, b, c, d] => select4(a, b, c, d, index),
                    _ => unreachable!(),
                };
                size_5_groups.chain(once(last_group_value)).collect()
            };

            let new_index = small_medians.len() / 2;
            let pivot = select(small_medians, new_index);
            let (nr_lt, nr_le) = three_way_partition(&mut elements, pivot);

            if index < nr_lt {
                elements.truncate(nr_lt);
                select(elements, index)
            } else if index < nr_le {
                pivot
            } else {
                elements.drain(..nr_le);
                select(elements, index - nr_le)
            }
        }
    }
}

/// Partitions `a` into: [ < pivot | == pivot | > pivot ]
/// Returns `(eq_start, eq_end)` so:
/// - `a[..eq_start]`         < pivot
/// - `a[eq_start..eq_end]`   == pivot
/// - `a[eq_end..]`           > pivot
pub fn three_way_partition(a: &mut [f64], pivot: f64) -> (usize, usize) {
    let mut lt = 0usize;
    let mut i = 0usize;
    let mut gt = a.len(); // exclusive upper bound

    while i < gt {
        // Compute ordering first so the borrow of `a[i]` doesn't overlap swaps.
        match a[i].total_cmp(&pivot) {
            Ordering::Less => {
                a.swap(lt, i);
                lt += 1;
                i += 1;
            }
            Ordering::Greater => {
                gt -= 1;
                a.swap(i, gt);
                // don't increment i: the swapped-in element at `i` is unclassified
            }
            Ordering::Equal => {
                i += 1;
            }
        }
    }

    (lt, gt)
}

/// Median selection network for 5 inputs, 7 compare-exchanges.
///
/// Comparator pairs (wire indices): (0,1)(2,3) (0,2)(1,3) (2,4) (1,2) (2,4)
#[inline(always)]
fn median5(mut a: f64, mut b: f64, mut c: f64, mut d: f64, mut e: f64) -> f64 {
    sort2(&mut a, &mut b);
    sort2(&mut c, &mut d);

    sort2(&mut a, &mut c);
    sort2(&mut b, &mut d);

    sort2(&mut c, &mut e);
    sort2(&mut b, &mut c);
    sort2(&mut c, &mut e);

    c // median
}

#[must_use]
#[inline(always)]
pub fn select4(mut a: f64, mut b: f64, mut c: f64, mut d: f64, index: usize) -> f64 {
    // Debug-only: these leaf selectors run once per 4-element group of every quickselect
    // (O(n) per call, O(m²) medians in the median functional); the `match` below still
    // hits `unreachable!()` on an invalid index in release.
    debug_assert!(index < 4);

    sort2(&mut a, &mut b); // a <= b
    sort2(&mut c, &mut d); // c <= d
    sort2(&mut a, &mut c); // a is global min
    sort2(&mut b, &mut d); // d is global max
    sort2(&mut b, &mut c);

    match index {
        0 => a,
        1 => b,
        2 => c,
        3 => d,
        _ => unreachable!(),
    }
}

#[must_use]
#[inline(always)]
pub fn select3(mut a: f64, mut b: f64, mut c: f64, index: usize) -> f64 {
    // Debug-only: see `select4`.
    debug_assert!(index < 3);

    // (0,1), (1,2), (0,1)
    sort2(&mut a, &mut b);
    sort2(&mut b, &mut c);
    sort2(&mut a, &mut b);

    match index {
        0 => a,
        1 => b,
        2 => c,
        _ => unreachable!(),
    }
}

#[must_use]
#[inline(always)]
pub fn select2(a: f64, b: f64, index: usize) -> f64 {
    // Debug-only: see `select4`.
    debug_assert!(index < 2);

    match index {
        0 => a.min(b),
        1 => a.max(b),
        _ => unreachable!(),
    }
}

#[inline(always)]
fn sort2(x: &mut f64, y: &mut f64) {
    if *x > *y {
        mem::swap(x, y);
    }
}

/// Performs a binary search over the integer index range `[lo, hi)`.
///
/// This function is analogous to [`slice::binary_search_by`], but it operates on
/// indices rather than elements of a slice. It repeatedly invokes the provided
/// comparator with a candidate index and uses the returned [`Ordering`] to
/// narrow the search interval until a match is found or the search space is
/// exhausted.
///
/// The comparator must define a total order consistent with the search target
/// over the range `lo..hi`:
/// - Return [`Ordering::Less`] if the value at the candidate index is strictly
///   less than the target.
/// - Return [`Ordering::Greater`] if it is strictly greater than the target.
/// - Return [`Ordering::Equal`] to indicate a match.
///
/// # Arguments
///
/// - `lo`: Lower bound (inclusive) of the search range.
/// - `hi`: Upper bound (exclusive) of the search range.
/// - `cmp`: Comparator called with a candidate index in `lo..hi`.
///
/// # Returns
///
/// - `Ok(index)` if the comparator returns `Equal` for some `index` in `lo..hi`.
/// - `Err(insertion_index)` if no match is found, where `insertion_index` is in
///   `lo..=hi` and indicates the position where a value equal to the target
///   could be inserted while preserving the order induced by `cmp`.
///
/// When multiple matching indices exist (i.e., the comparator would return
/// `Equal` at several positions), any one of them may be returned; no stability
/// is guaranteed.
///
/// # Panics
///
/// - Panics if `lo > hi`.
///
/// # Complexity
///
/// - Performs at most `O(log(hi - lo))` calls to the comparator.
///
/// # Examples
///
/// Basic usage with a conceptual sorted array:
/// ```rust
/// use isodistrreg::routines::binary_search_by_index;
///
/// // Conceptual array: [1, 3, 5, 7, 9]
/// let needle = 7;
/// let res = binary_search_by_index(0, 5, false, |i| {
///     let val = [1, 3, 5, 7, 9][i];
///     val.cmp(&needle)
/// });
/// assert_eq!(res, Ok(3));
/// ```
///
/// Getting the insertion index when not found:
/// ```rust
/// use isodistrreg::routines::binary_search_by_index;
///
/// // Conceptual array: [10, 20, 30, 40]
/// let needle = 25;
/// let res = binary_search_by_index(0, 4, false, |i| {
///     [10, 20, 30, 40][i].cmp(&needle)
/// });
/// // 25 would be inserted at position 2
/// assert_eq!(res, Err(2));
/// ```
///
/// Searching a monotonic computation without allocating:
/// ```rust
/// use isodistrreg::routines::binary_search_by_index;
///
/// // f(i) = 2*i is monotonic increasing on i >= 0
/// let needle = 14;
/// let res = binary_search_by_index(0, 10, false, |i| {
///     let v = 2 * i;
///     v.cmp(&needle)
/// });
/// assert_eq!(res, Ok(7));
/// ```
///
/// Outside the range:
/// ```rust
/// use isodistrreg::routines::binary_search_by_index;
///
/// assert_eq!(binary_search_by_index(0, 10, false, |i| (2 * i).cmp(&21)), Err(10));
/// assert_eq!(binary_search_by_index(0, 10, false, |i| (2 * i as i32).cmp(&-1_i32)), Err(0));
/// ```
///
/// Returns the left-most value:
/// ```rust
/// use isodistrreg::routines::binary_search_by_index;
///
/// assert_eq!(binary_search_by_index(0, 100, false, |i| (i / 10).cmp(&5)), Ok(50));
/// assert_eq!(binary_search_by_index(0, 100, true, |i| (i / 10).cmp(&5)), Ok(59));
/// ```
///
/// [`slice::binary_search_by`]: https://doc.rust-lang.org/std/primitive.slice.html#method.binary_search_by
/// [`Ordering`]: https://doc.rust-lang.org/std/cmp/enum.Ordering.html
/// [`Ordering::Less`]: https://doc.rust-lang.org/std/cmp/enum.Ordering.html#variant.Less
/// [`Ordering::Greater`]: https://doc.rust-lang.org/std/cmp/enum.Ordering.html#variant.Greater
/// [`Ordering::Equal`]: https://doc.rust-lang.org/std/cmp/enum.Ordering.html#variant.Equal
#[inline]
pub fn binary_search_by_index<F>(
    lo: usize,
    hi: usize,
    upper: bool,
    mut cmp: F,
) -> Result<usize, usize>
where
    F: FnMut(usize) -> Ordering,
{
    assert!(lo <= hi, "lower bound must be <= upper bound");
    let mut left = lo;
    let mut right = hi;

    if upper {
        // Find the highest index where cmp returns Equal.
        // Treat Equal like Less to push `left` as far right as possible.
        while left < right {
            let mid = left + (right - left) / 2;
            match cmp(mid) {
                Ordering::Less | Ordering::Equal => left = mid + 1,
                Ordering::Greater => right = mid,
            }
        }

        // `left` is now one past the last Equal (or the insertion point).
        // Check if the element just before `left` is Equal.
        if left > lo && cmp(left - 1) == Ordering::Equal {
            Ok(left - 1)
        } else {
            Err(left)
        }
    } else {
        // Find the lowest index where cmp returns Equal.
        // Treat Equal like Greater to pull `right` as far left as possible.
        while left < right {
            let mid = left + (right - left) / 2;
            match cmp(mid) {
                Ordering::Less => left = mid + 1,
                Ordering::Greater | Ordering::Equal => right = mid,
            }
        }

        if left < hi && cmp(left) == Ordering::Equal {
            Ok(left)
        } else {
            Err(left)
        }
    }
}

#[must_use]
pub fn lexicographic_order<T: PartialOrd>(covariates: &[T], n: usize, d: usize) -> Vec<usize> {
    assert_eq!(covariates.len(), n * d);

    let get_cdf = |i| &covariates[i * d..(i + 1) * d];
    argsort_unstable_by::<Increasing, _>(|a, b| lexicographic_cmp(get_cdf(a), get_cdf(b)), n)
}

/// Compares two vectors of equal dimension in lexicographical order.
///
/// Element-wise comparison of the two vectors from low indices to high, with -0.0 and 0.0 treated
/// as equal (when called on floats — `PartialOrd` semantics apply).
#[must_use]
#[inline]
pub fn lexicographic_cmp<T: PartialOrd>(left: &[T], right: &[T]) -> Ordering {
    for (v_l, v_r) in left.iter().zip(right) {
        if v_l < v_r {
            return Ordering::Less;
        }
        if v_l > v_r {
            return Ordering::Greater;
        }
    }
    Ordering::Equal
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::structures::{Decreasing, Increasing};

    /// One Kaplan-Meier value per unique response; within a tie, deaths apply before
    /// censorings leave the at-risk set.
    #[test]
    fn kaplan_meier_tied_responses() {
        let obs = |y: f64, observed: bool| Observation {
            x: (),
            y,
            observed,
            weight: 1.0f32,
        };

        // death@1, censor@1, death@2: S(1) = 2/3, S(2) = 2/3 * 0 -> CDF [1/3, 1].
        let cdf = kaplan_meier(
            [obs(1.0, true), obs(1.0, false), obs(2.0, true)].into_iter(),
            3.0,
        );
        assert_eq!(cdf.len(), 2);
        assert!((cdf[0] - 1.0 / 3.0).abs() < 1e-6, "{cdf:?}");
        assert!((cdf[1] - 1.0).abs() < 1e-6, "{cdf:?}");

        // Tied deaths fold into their group's value: ECDF of {1, 1, 2} = [2/3, 1].
        let cdf = kaplan_meier(
            [obs(1.0, true), obs(1.0, true), obs(2.0, true)].into_iter(),
            3.0,
        );
        assert_eq!(cdf.len(), 2);
        assert!((cdf[0] - 2.0 / 3.0).abs() < 1e-6, "{cdf:?}");
        assert!((cdf[1] - 1.0).abs() < 1e-6, "{cdf:?}");
    }

    /// The final cumulative value is total/total: exactly 1.0, independent of
    /// accumulation round-off.
    #[test]
    fn empirical_cdf_ends_at_exactly_one() {
        let weights = [1.0f32, 1.5, 0.75, 0.5, 1.5, 0.75];
        let cdf = empirical_cdf(
            weights.iter().map(|&weight| Observation {
                x: (),
                y: (),
                observed: (),
                weight,
            }),
            weights.iter().sum(),
        );
        assert_eq!(*cdf.last().unwrap(), 1.0);
        assert!(cdf.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn transpose_degenerate() {
        // Empty matrix: no-op
        let mut empty: Vec<i32> = vec![];
        transpose(&mut empty, 0, 0);
        assert!(empty.is_empty());

        // 1×N becomes N×1; the linear order is unchanged
        let mut v = vec![10, 20, 30];
        transpose(&mut v, 1, 3);
        assert_eq!(v, vec![10, 20, 30]); // now interpreted as 3×1
    }

    #[test]
    fn transpose_square() {
        let mut a = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        transpose(&mut a, 3, 3);
        assert_eq!(a, vec![1, 4, 7, 2, 5, 8, 3, 6, 9,]);
    }

    #[test]
    fn transpose_() {
        let mut a = vec![1, 2, 3, 4, 5, 6, 7, 8]; // 2×4, row-major
        transpose(&mut a, 2, 4);
        assert_eq!(
            a,
            vec![1, 5, 2, 6, 3, 7, 4, 8,] // 4×2, row-major
        );
    }

    #[test]
    fn sort_increasing() {
        let key = vec![30, 10, 20, 40];
        let index = argsort_unstable_by::<Increasing, _>(|i, j| key[i].cmp(&key[j]), 4);
        assert!(index.into_iter().map(|i| key[i]).is_sorted());
    }

    #[test]
    fn sort_decreasing() {
        let key = vec![30, 10, 20, 40];
        let index = argsort_unstable_by::<Decreasing, _>(|i, j| key[i].cmp(&key[j]), 4);
        let mut reversed = index.into_iter().map(|i| key[i]).collect::<Vec<_>>();
        reversed.reverse();
        assert!(reversed.is_sorted());
    }

    #[test]
    fn test_median() {
        assert_eq!(median([6.].into_iter()), 6.);
        assert_eq!(median([5., 6., 7.].into_iter()), 6.);
        assert_eq!(median([7., 6., 5.].into_iter()), 6.);
        assert_eq!(median([7., 6., 6.].into_iter()), 6.);
        assert_eq!(median([6., 6., 7.].into_iter()), 6.);
        assert_eq!(median([6., 7., 6.].into_iter()), 6.);
        assert_eq!(median([5., 3., 4., 3., 6.].into_iter()), 4.);
        assert_eq!(median((0..=100).map(|i| i as f64)), 50.0);
        assert_eq!(median((0..=100).rev().map(|i| i as f64)), 50.0);
    }

    #[test]
    fn test_select_small() {
        assert_eq!(select2(1.0, 2.0, 1), 2.0);
        assert_eq!(select2(2.0, 1.0, 1), 2.0);
        assert_eq!(select2(2.0, 1.0, 0), 1.0);
        assert_eq!(select3(3.0, 2.0, 1.0, 0), 1.0);
        assert_eq!(select3(3.0, 1.0, 1.0, 1), 1.0);
        assert_eq!(select4(3.0, 1.0, 1.0, 2.0, 1), 1.0);
        assert_eq!(select4(3.0, 1.0, 1.0, 2.0, 2), 2.0);
        assert_eq!(select4(3.0, 1.0, 1.0, 2.0, 3), 3.0);
    }

    #[test]
    fn test_select() {
        assert_eq!(select((0..=100).map(|i| i as f64).collect(), 17), 17.0);
        assert_eq!(select((0..=100).map(|i| i as f64).collect(), 99), 99.0);
        assert_eq!(
            select((0..=100).rev().map(|i| i as f64).collect(), 17),
            17.0
        );
        assert_eq!(
            select((0..=100).rev().map(|i| i as f64).collect(), 99),
            99.0
        );
    }

    #[test]
    fn test_zero_treatment() {
        assert_eq!(median([-1.0, 0.0, -0.0, -0.0, 0.0]), -0.0);
    }
}
