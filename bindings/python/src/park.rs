use rayon::prelude::*;
use std::f64::consts::PI;
// ========================== Data Structures ==========================

/// Pre-computed Kaplan–Meier statistics for a single group.
struct GroupData {
    t_sorted: Vec<f64>,
    n_total: usize,
    evt_times: Vec<f64>,
    d: Vec<f64>,
    n: Vec<f64>,
    cum_log_surv: Vec<f64>,
}

/// Result returned by [`fit_park`].
#[derive(Debug, PartialEq)]
pub struct ParkFitResult {
    /// CDF grid, row-major, shape `(n_buckets, n_times)`.
    pub fit: Vec<f64>,
    /// Center value for each active bucket.
    pub bucket_centers: Vec<f64>,
    /// Evaluation time grid.
    pub times: Vec<f64>,
}

// ========================== Utility Helpers ==========================

/// `np.searchsorted(a, v, side="left")`
#[inline]
fn searchsorted_left(sorted: &[f64], value: f64) -> usize {
    sorted.partition_point(|&x| x < value)
}

/// `np.searchsorted(a, v, side="right")`
#[inline]
fn searchsorted_right(sorted: &[f64], value: f64) -> usize {
    sorted.partition_point(|&x| x <= value)
}

/// Map each element of `x` to the index of its nearest value in the
/// sorted slice `centers`.
fn closest_center(x: &[f64], centers: &[f64]) -> Vec<usize> {
    x.iter()
        .map(|&xi| {
            let pos = centers.partition_point(|&c| c < xi);
            if pos == 0 {
                0
            } else if pos >= centers.len() {
                centers.len() - 1
            } else if (xi - centers[pos - 1]).abs() <= (xi - centers[pos]).abs() {
                pos - 1
            } else {
                pos
            }
        })
        .collect()
}

/// Run-length encode a **sorted** slice into `(unique_values, counts)`.
fn unique_with_counts(sorted: &[f64]) -> (Vec<f64>, Vec<usize>) {
    if sorted.is_empty() {
        return (vec![], vec![]);
    }
    let mut vals = vec![sorted[0]];
    let mut cnts = vec![1_usize];
    for &v in &sorted[1..] {
        if v == *vals.last().unwrap() {
            *cnts.last_mut().unwrap() += 1;
        } else {
            vals.push(v);
            cnts.push(1);
        }
    }
    (vals, cnts)
}

/// `f64::max` that **propagates** NaN (matching `np.maximum`).
#[inline]
fn max_nan(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a >= b {
        a
    } else {
        b
    }
}

// ===================== Park Algorithm Internals =====================

#[inline]
fn n_at_risk(gd: &GroupData, x: f64) -> i64 {
    gd.n_total as i64 - searchsorted_left(&gd.t_sorted, x) as i64
}

#[inline]
fn m_events(gd: &GroupData, x: f64) -> usize {
    searchsorted_right(&gd.evt_times, x)
}

fn solve_k_for_q(nj: &[f64], dj: &[f64], q: f64, epsilon: f64) -> f64 {
    if q >= 0.0 {
        return f64::INFINITY;
    }

    let lb = dj
        .iter()
        .zip(nj.iter())
        .map(|(&d, &n)| d - n)
        .fold(f64::NEG_INFINITY, f64::max)
        + 1e-12;

    let f = |k: f64| -> f64 {
        nj.iter()
            .zip(dj.iter())
            .map(|(&n, &d)| (-d / (n + k)).ln_1p())
            .sum::<f64>()
            - q
    };

    let mut ub = 1.0_f64;
    while f(ub) < 0.0 {
        ub *= 2.0;
        if ub > 1e18 {
            return f64::INFINITY;
        }
    }

    let (mut lo, mut hi) = (lb, ub);
    while hi - lo >= epsilon {
        let mid = 0.5 * (lo + hi);
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn k_of_q_x(gd: &GroupData, q: f64, x: f64, epsilon: f64) -> f64 {
    let n = n_at_risk(gd, x);
    let m = m_events(gd, x);

    if n <= 0 {
        return 0.0;
    }
    if m == 0 {
        return -(n as f64);
    }

    let nj = &gd.n[..m];
    let dj = &gd.d[..m];

    let k_hat = if q >= 0.0 {
        f64::INFINITY
    } else {
        solve_k_for_q(nj, dj, q, epsilon)
    };

    (-(n as f64)).max(k_hat)
}

fn qstar_unrestricted(gd: &GroupData, x: f64) -> f64 {
    let m = m_events(gd, x);
    if m == 0 { 0.0 } else { gd.cum_log_surv[m - 1] }
}

fn block_q_hat(
    gds: &[GroupData],
    active: &[usize],
    index_start: usize,
    index_end: usize,
    x: f64,
    epsilon: f64,
) -> f64 {
    let block = &active[index_start..index_end];

    let any_events = block
        .iter()
        .any(|&i| m_events(&gds[i], x) > 0 && n_at_risk(&gds[i], x) > 0);
    if !any_events {
        return 0.0;
    }

    let sum_k = |q: f64| -> f64 {
        block
            .iter()
            .map(|&i| k_of_q_x(&gds[i], q, x, epsilon))
            .sum()
    };

    let qmin = block
        .iter()
        .map(|&i| qstar_unrestricted(&gds[i], x))
        .fold(f64::INFINITY, f64::min);

    let mut lo = if qmin.is_finite() { qmin - 1.0 } else { -1.0 };
    let mut s_lo = sum_k(lo);
    for _ in 0..80 {
        if s_lo < 0.0 {
            break;
        }
        lo *= 2.0;
        s_lo = sum_k(lo);
    }
    if s_lo >= 0.0 {
        return 0.0;
    }

    let mut hi = 0.0_f64;
    while hi - lo >= epsilon {
        let mid = 0.5 * (lo + hi);
        let s_mid = sum_k(mid);
        let s_mid = if s_mid.is_finite() { s_mid } else { 1.0 };
        if s_mid < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

// ========================== PAVA Step ===============================

struct PavaBlock {
    lo: usize,
    hi: usize,
    q: f64,
}

/// Run the pool-adjacent-violators algorithm at one time point.
/// Returns `Some(log_S_column)` when at least one group is active, else `None`.
fn process_timepoint(
    x: f64,
    gds: &[GroupData],
    n_col: &[i32],
    qstar_col: &[f64],
    epsilon: f64,
) -> Option<Vec<f64>> {
    let active: Vec<usize> = n_col
        .iter()
        .enumerate()
        .filter(|&(_, &n)| n > 0)
        .map(|(i, _)| i)
        .collect();

    if active.is_empty() {
        return None;
    }

    let n_groups = n_col.len();
    let mut log_s = vec![f64::NAN; n_groups];

    let mut blocks: Vec<PavaBlock> = Vec::new();
    for (pos, &i) in active.iter().enumerate() {
        blocks.push(PavaBlock {
            lo: pos,
            hi: pos,
            q: qstar_col[i],
        });

        while blocks.len() >= 2 {
            let len = blocks.len();
            if blocks[len - 2].q < blocks[len - 1].q {
                let b2 = blocks.pop().unwrap();
                let b1 = blocks.pop().unwrap();
                let lo = b1.lo;
                let hi = b2.hi;
                let qhat = block_q_hat(gds, &active, lo, hi + 1, x, epsilon);
                blocks.push(PavaBlock { lo, hi, q: qhat });
            } else {
                break;
            }
        }
    }

    for b in &blocks {
        for pos in b.lo..=b.hi {
            log_s[active[pos]] = b.q;
        }
    }

    Some(log_s)
}

// ====================== Group Data Preparation ======================

fn prep_group_data(tg: &[f64], eg: &[bool]) -> GroupData {
    let mut order: Vec<usize> = (0..tg.len()).collect();
    order.sort_by(|&a, &b| {
        tg[a]
            .partial_cmp(&tg[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let t_sorted: Vec<f64> = order.iter().map(|&i| tg[i]).collect();
    let e_sorted: Vec<bool> = order.iter().map(|&i| eg[i]).collect();
    let n_total = t_sorted.len();

    // sorted event times (subset of t_sorted where event occurred)
    let event_times_raw: Vec<f64> = t_sorted
        .iter()
        .zip(e_sorted.iter())
        .filter(|&(_, &e)| e)
        .map(|(&t, _)| t)
        .collect();

    let (evt_times, d_counts) = unique_with_counts(&event_times_raw);
    let d: Vec<f64> = d_counts.iter().map(|&c| c as f64).collect();

    let (n, cum_log_surv) = if evt_times.is_empty() {
        (vec![], vec![])
    } else {
        let n_at_event: Vec<f64> = evt_times
            .iter()
            .map(|&et| (n_total - searchsorted_left(&t_sorted, et)) as f64)
            .collect();

        let steps: Vec<f64> = d
            .iter()
            .zip(n_at_event.iter())
            .map(|(&dj, &nj)| (-dj / nj).ln_1p())
            .collect();

        let mut cum = Vec::with_capacity(steps.len());
        let mut acc = 0.0;
        for &s in &steps {
            acc += s;
            cum.push(acc);
        }

        (n_at_event, cum)
    };

    GroupData {
        t_sorted,
        n_total,
        evt_times,
        d,
        n,
        cum_log_surv,
    }
}

// ============================= Public API =============================

/// Pointwise constrained NPMLE of stochastically ordered survivor functions
/// (Park, Taylor & Kalbfleisch, *Biometrika* 2012) for the **simple ordering**
/// case:
///
/// \( S_1(t) \ge S_2(t) \ge \cdots \ge S_C(t) \) for all \( t \),
///
/// where groups come from mapping each `x[i]` to its closest value in
/// `centers`.
///
/// Returns a [`ParkFitResult`] whose `fit` field stores
/// \( \hat F(t \mid \text{bucket}) = 1 - \hat S(t) \) in row-major layout
/// `(n_buckets, n_times)`.
///
/// # Arguments
///
/// * `x`        – covariate value for each observation.
/// * `time`     – observed event / censoring time.
/// * `event`    – event indicator (`true` = event, `false` = censored).
/// * `centers`  – grouping centers (sorted ascending). If `None`, `unique(x)` is used.
/// * `epsilon`  – bisection tolerance (e.g. `1e-6`).
/// * `parallel` – use rayon parallelism over time points.
pub fn fit_park(
    x: &[f64],
    time: &[f64],
    event: &[bool],
    centers: Option<&[f64]>,
    epsilon: f64,
    parallel: bool,
) -> ParkFitResult {
    let n = x.len();
    assert_eq!(n, time.len());
    assert_eq!(n, event.len());
    assert!(n > 0, "no observations");

    if let Some(c) = centers {
        assert!(c.windows(2).all(|w| w[0] < w[1]));
    }

    let mut maybe_storage = None;
    let centers = centers.unwrap_or_else(|| {
        maybe_storage = Some(default_centers(x));
        maybe_storage.as_deref().unwrap()
    });

    // ---- assign groups ----
    let groups = closest_center(x, centers);

    // ---- sort by group ----
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| groups[i]);

    let sorted_groups: Vec<usize> = order.iter().map(|&i| groups[i]).collect();
    let sorted_time: Vec<f64> = order.iter().map(|&i| time[i]).collect();
    let sorted_event: Vec<bool> = order.iter().map(|&i| event[i]).collect();

    // ---- bucket boundaries ----
    let mut starts = vec![0_usize];
    for i in 1..n {
        if sorted_groups[i] != sorted_groups[i - 1] {
            starts.push(i);
        }
    }
    let mut stops: Vec<usize> = starts[1..].to_vec();
    stops.push(n);

    let buckets: Vec<usize> = starts.iter().map(|&s| sorted_groups[s]).collect();
    let n_buckets = buckets.len();

    // ---- evaluation time grid ----
    let mut time_grid: Vec<f64> = Vec::with_capacity(n);
    for i in 0..n {
        if sorted_event[i] {
            time_grid.push(sorted_time[i]);
        } else {
            time_grid.push(sorted_time[i].next_up());
        }
    }
    time_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
    time_grid.dedup();
    let n_times = time_grid.len();

    // ---- per-group data ----
    let gds: Vec<GroupData> = starts
        .iter()
        .zip(stops.iter())
        .map(|(&s, &e)| prep_group_data(&sorted_time[s..e], &sorted_event[s..e]))
        .collect();
    debug_assert_eq!(gds.len(), n_buckets);

    // ---- pre-compute N (at risk) and q* for every (group, time) cell ----
    //      row-major, shape (n_buckets, n_times)
    let mut n_matrix = vec![0_i32; n_buckets * n_times];
    let mut qstar_matrix = vec![0.0_f64; n_buckets * n_times];

    for (gi, gd) in gds.iter().enumerate() {
        let base = gi * n_times;
        for ti in 0..n_times {
            let xt = time_grid[ti];
            n_matrix[base + ti] = gd.n_total as i32 - searchsorted_left(&gd.t_sorted, xt) as i32;
            let m = searchsorted_right(&gd.evt_times, xt);
            if m > 0 {
                qstar_matrix[base + ti] = gd.cum_log_surv[m - 1];
            }
        }
    }

    // ---- main loop: PAVA at each time point ----
    let mut log_s = vec![f64::NAN; n_buckets * n_times];

    // Helper closure: extract column vectors and run the PAVA step.
    let process_tp = |ti: usize| -> (usize, Option<Vec<f64>>) {
        let n_col: Vec<i32> = (0..n_buckets)
            .map(|gi| n_matrix[gi * n_times + ti])
            .collect();
        let qstar_col: Vec<f64> = (0..n_buckets)
            .map(|gi| qstar_matrix[gi * n_times + ti])
            .collect();
        let col = process_timepoint(time_grid[ti], &gds, &n_col, &qstar_col, epsilon);
        (ti, col)
    };

    let results: Vec<(usize, Option<Vec<f64>>)> = {
        if parallel {
            (0..n_times).into_par_iter().map(process_tp).collect()
        } else {
            (0..n_times).map(process_tp).collect()
        }
    };

    for (ti, opt_col) in results {
        if let Some(col) = opt_col {
            for gi in 0..n_buckets {
                log_s[gi * n_times + ti] = col[gi];
            }
        }
    }

    // ---- post-processing ----
    // cdfs = 1 − exp(logS), clip to [0, 1]
    let mut fit: Vec<f64> = log_s.iter().map(|&v| 1.0 - v.exp()).collect();
    for v in fit.iter_mut() {
        *v = v.clamp(0.0, 1.0);
    }

    // cumulative max along time axis (NaN propagates, matching np.maximum)
    for gi in 0..n_buckets {
        let base = gi * n_times;
        for ti in 1..n_times {
            fit[base + ti] = max_nan(fit[base + ti - 1], fit[base + ti]);
        }
    }

    // mask: once CDF hits 1.0 propagate along time, then along groups
    let mut mask = vec![false; n_buckets * n_times];
    for i in 0..mask.len() {
        mask[i] = fit[i] == 1.0;
    }
    for gi in 0..n_buckets {
        let base = gi * n_times;
        for ti in 1..n_times {
            if mask[base + ti - 1] {
                mask[base + ti] = true;
            }
        }
    }
    for ti in 0..n_times {
        for gi in 1..n_buckets {
            if mask[(gi - 1) * n_times + ti] {
                mask[gi * n_times + ti] = true;
            }
        }
    }
    for i in 0..fit.len() {
        if mask[i] {
            fit[i] = 1.0;
        }
    }

    // ---- result ----
    let bucket_centers: Vec<f64> = buckets.iter().map(|&b| centers[b]).collect();

    debug_assert_eq!(bucket_centers.len(), n_buckets);
    debug_assert_eq!(time_grid.len(), n_times);

    ParkFitResult {
        fit,
        bucket_centers,
        times: time_grid,
    }
}

/// Place kernel centers at evenly-spaced quantiles of `x`,
/// with the number of buckets derived from `default_bw`.
fn default_centers(x: &[f64]) -> Vec<f64> {
    let bw = default_bw(x);
    let bin_width = (6_f64 / PI.powf(0.5)).powf(1.0 / 3.0) * bw;
    let range = x.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - x.iter().copied().fold(f64::INFINITY, f64::min);
    let bin_count = ((range / bin_width).round() as usize).clamp(1, x.len());

    // Build the odd-indexed entries of linspace(0, 1, 2*n_buckets+1),
    // i.e. indices 1, 3, 5, …, 2*n_buckets-1 of that grid.
    // These are simply (2*k + 1) / (2*n_buckets) for k in 0..n_buckets.
    let centers_quantile_space: Vec<f64> = (0..bin_count)
        .map(|k| (2 * k + 1) as f64 / (2 * bin_count) as f64)
        .collect();

    // Pre-sort once for all quantile lookups
    let mut sorted = x.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);

    centers_quantile_space
        .iter()
        .map(|&q| quantile_linear(&sorted, q))
        .collect()
}

/// Linearly-interpolated quantile on a **sorted** slice,
/// where `q` is in \[0, 1\].  Matches NumPy's default `method="linear"`.
fn quantile_linear(sorted: &[f64], q: f64) -> f64 {
    debug_assert!(!sorted.is_empty());
    let rank = q * (sorted.len() - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = lo + 1;
    let frac = rank - lo as f64;
    if hi >= sorted.len() {
        sorted[lo]
    } else {
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

/// Silverman's rule-of-thumb bandwidth with robust scale estimate.
fn default_bw(x: &[f64]) -> f64 {
    let n = x.len();
    if n < 2 {
        return 1.0;
    }

    // sample standard deviation (ddof = 1)
    let nf = n as f64;
    let mean = x.iter().sum::<f64>() / nf;
    let var = x.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
    let sd = var.sqrt();

    // IQR via linear-interpolated percentiles (numpy default)
    let mut sorted = x.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);

    let q25 = percentile_linear(&sorted, 25.0);
    let q75 = percentile_linear(&sorted, 75.0);
    let iqr = q75 - q25;

    // robust scale estimate
    let s = if iqr > 0.0 { sd.min(iqr / 1.349) } else { sd };

    // fallback if x is (nearly) constant
    let s = if !s.is_finite() || s <= 0.0 {
        f64::EPSILON + x.iter().map(|&v| v.abs()).sum::<f64>() / nf
    } else {
        s
    };

    0.9 * s * nf.powf(-0.2)
}

/// Linearly-interpolated percentile on a **sorted** slice,
/// matching NumPy's default `method="linear"`.
fn percentile_linear(sorted: &[f64], pct: f64) -> f64 {
    debug_assert!(!sorted.is_empty());
    let rank = pct / 100.0 * (sorted.len() - 1) as f64;
    let lo = rank.floor() as usize;
    let hi = lo + 1;
    let frac = rank - lo as f64;
    if hi >= sorted.len() {
        sorted[lo]
    } else {
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}
