use crate::structures::{Direction, Increasing, Observation};
use std::cmp::Ordering;
use std::iter::once;
use std::mem;

/// Transposes a row-major `m × n` matrix stored in `slice`.
pub fn transpose<T: Copy>(matrix: &mut Vec<T>, m: usize, n: usize) {
    assert_eq!(matrix.len(), m * n);
    if m == 0 || n == 0 {
        return;
    }

    if m == n {
        // In-place symmetric swap across the diagonal
        for i in 0..m {
            for j in (i + 1)..n {
                matrix.swap(i * n + j, j * n + i);
            }
        }
        return;
    }

    let mut out = Vec::with_capacity(matrix.len());
    // Iterate by the new index
    for i in 0..n {
        for j in 0..m {
            out.push(matrix[j * n + i]);
        }
    }
    *matrix = out;
}

pub fn argsort_unstable_by<D: Direction, F>(mut cmp: F, n: usize) -> Vec<usize>
where
    F: FnMut(usize, usize) -> Ordering,
{
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_unstable_by(|&i, &j| {
        let ord = cmp(i, j);
        if D::IS_INCREASING { ord } else { ord.reverse() }
    });
    idx
}

// Compute the empirical CDF. Assumes that observations have been deduplicated so responses are
// unique.
pub fn empirical_cdf<C, I: Into<Observation<C, f64, ()>>>(
    observations: impl Iterator<Item = I>,
    total_weight: f64,
) -> Vec<f64> {
    observations
        .map(Into::into)
        .scan(0.0, |acc, o| {
            *acc += o.weight / total_weight;
            // Ensure numerical imprecision won't get us out of [0, 1]
            *acc = acc.clamp(0.0, 1.0);
            Some(*acc)
        })
        .collect()
}

/// Compute the Kaplan-Meier estimator on unique responses.
pub fn kaplan_meier<C, R: Copy + PartialEq>(
    observations: impl Iterator<Item = Observation<C, R, bool>>,
    total_weight: f64,
) -> Vec<f64> {
    // Observations have been deduplicated so responses are unique
    observations
        .scan(
            (1.0, total_weight, None),
            |(s, total_weight, previous), o| {
                if o.observed {
                    *s *= 1.0 - o.weight / *total_weight;
                }
                *total_weight -= o.weight;

                match previous.replace(o.y) {
                    Some(value) if value == o.y => None,
                    _ => Some(*s),
                }
            },
        )
        // Clamping for numerics
        .map(|v| 1.0 - v.clamp(0.0, 1.0))
        .collect()
}

pub fn median(elements: impl IntoIterator<Item = f64>) -> f64 {
    let allocated: Vec<_> = elements.into_iter().collect();
    assert!(!allocated.len().is_multiple_of(2));
    let median_index = allocated.len() / 2;
    select(allocated, median_index)
}

#[must_use]
pub fn select(mut elements: Vec<f64>, index: usize) -> f64 {
    assert!(index < elements.len());

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
    assert!(index < 4);

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
    assert!(index < 3);

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
    assert!(index < 2);

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
pub fn lexicographic_order(covariates: &[f64], n: usize, d: usize) -> Vec<usize> {
    assert_eq!(covariates.len(), n * d);

    let get_cdf = |i| &covariates[i * d..(i + 1) * d];
    argsort_unstable_by::<Increasing, _>(|a, b| lexicographic_cmp(get_cdf(a), get_cdf(b)), n)
}

/// Compares two vectors of equal dimension in lexicographical order.
///
/// Element-wise comparison of the two vectors from low indices to high, with -0.0 and 0.0 treated
/// as equal.
#[must_use]
#[inline]
pub fn lexicographic_cmp(left: &[f64], right: &[f64]) -> Ordering {
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
