//! Sparse Cholesky factorization of the ADMM linear system.
//!
//! The ADMM iteration for the isotonic QP repeatedly solves
//!
//! ```text
//! K x = b,    K = (1 + sigma) * I + rho * A^T A
//! ```
//!
//! where `A` is the `m x n` constraint matrix of the cover relation: every row
//! has exactly two nonzeros, sitting at the two endpoints of one cover edge. The
//! matrix is handed over in factored form,
//!
//! - `edges[r] = (i, j)` with `i != j`, the columns touched by row `r`,
//! - `coef[r] = (a, b)`, with `a` in column `edges[r].0` and `b` in column `edges[r].1`,
//!
//! and the unordered pairs `{i, j}` are unique (the edge set is a transitive
//! reduction, so there are neither duplicate cover relations nor 2-cycles).
//!
//! Because of that structure `A^T A` never has to be formed:
//!
//! ```text
//! K[c, c] = 1 + sigma + rho * sum_{r touching c} (coef of row r at c)^2
//! K[i, j] = rho * coef[r].0 * coef[r].1                for the edge r = (i, j)
//! ```
//!
//! so the upper triangle of `K` has exactly `n + m` nonzeros and the *pattern*
//! depends only on `edges`. That is what makes the split below worthwhile: the
//! symbolic analysis (fill-reducing AMD ordering + elimination tree + supernode
//! detection) runs once per fit in [`KktFactor::new`], and every threshold only
//! pays for [`KktFactor::refactor`], which is an `O(n + m)` value assembly plus
//! one numeric factorization into preallocated storage.
//!
//! `K` is symmetric positive definite with `lambda_min(K) >= 1 + sigma > 0` for
//! any `sigma >= 0`, since `A^T A` is positive semidefinite. Plain `LL^T` is
//! therefore unconditionally valid: no `LDL^T`, no pivoting, no regularization.
//!
//! # Allocation
//!
//! Everything that the numeric phase and the solve need — the `K` value arrays
//! for both precisions, the `L` value arrays for both precisions, the narrowed
//! right-hand side, and the `MemStack` scratch — is allocated in
//! [`KktFactor::new`]. [`KktFactor::refactor`] and [`KktFactor::solve_in_place`]
//! perform no heap traffic at all.
//!
//! # Determinism
//!
//! Every faer entry point is called with [`Par::Seq`]. The factorization and the
//! triangular solves are then a fixed sequence of scalar operations, so repeated
//! runs on the same input are bit-identical.

use std::cell::RefCell;

use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::linalg::cholesky::llt::factor::{LltParams, LltRegularization};
use faer::sparse::linalg::SupernodalThreshold;
use faer::sparse::linalg::cholesky::{
    CholeskySymbolicParams, LltRef, SymbolicCholesky, SymmetricOrdering,
    factorize_symbolic_cholesky,
};
use faer::sparse::{SparseColMatRef, SymbolicSparseColMatRef};
use faer::traits::ComplexField;
use faer::{Conj, MatMut, Par, Side, Spec};

/// A reusable `LL^T` factorization of `K = (1 + sigma) I + rho A^T A`.
///
/// Construct once per fit with [`KktFactor::new`], then call
/// [`KktFactor::refactor`] whenever `sigma`, `rho` or `coef` change and
/// [`KktFactor::solve_in_place`] for each right-hand side.
pub(crate) struct KktFactor {
    /// Number of columns of `A` (the dimension of `K`).
    n: usize,
    /// Number of rows of `A` (the number of cover edges).
    m: usize,

    // ---- sparsity pattern of the upper triangle of K, in CSC ----
    /// Column pointers, length `n + 1`.
    col_ptr: Vec<usize>,
    /// Row indices, length `n + m`, ascending within each column.
    row_idx: Vec<usize>,
    /// `diag_slot[c]` is the index into the value array holding `K[c, c]`.
    diag_slot: Vec<usize>,
    /// `edge_slot[r]` is the index into the value array holding the off-diagonal
    /// contributed by `edges[r]`.
    edge_slot: Vec<usize>,

    // ---- numeric buffers ----
    /// Symbolic factorization (AMD ordering); pattern-only, so shared by both
    /// precisions.
    symbolic: SymbolicCholesky<usize>,
    /// Values of the upper triangle of `K`, `f64` flavour.
    k_val_f64: Vec<f64>,
    /// Values of the upper triangle of `K`, `f32` flavour.
    k_val_f32: Vec<f32>,
    /// Values of `L`, `f64` flavour.
    l_val_f64: Vec<f64>,
    /// Values of `L`, `f32` flavour.
    l_val_f32: Vec<f32>,
    /// Narrowed right-hand side used by the `f32` solve path.
    rhs_f32: RefCell<Vec<f32>>,
    /// Scratch for faer, sized for the largest of the four operations
    /// (factorize/solve x f64/f32). `RefCell` because `solve_in_place` only
    /// takes `&self` but `MemStack::new` needs `&mut`.
    stack: RefCell<MemBuffer>,

    /// Precision of the factor currently stored in `l_val_f32` / `l_val_f64`.
    use_f32: bool,
    /// Whether [`Self::refactor`] has run at least once.
    factored: bool,
}

impl KktFactor {
    /// Symbolic phase: analyse the pattern of `K` and allocate every buffer.
    ///
    /// Runs once per fit; the pattern is a function of `edges` alone and does not
    /// change as thresholds, `sigma`, `rho` or `coef` move.
    ///
    /// # Panics
    ///
    /// Panics if any edge is a self-loop or references a column `>= n`, or if the
    /// symbolic factorization runs out of memory.
    pub(crate) fn new(n: usize, edges: &[(usize, usize)]) -> Self {
        Self::build(n, edges, SupernodalThreshold::AUTO)
    }

    /// [`Self::new`], but pinning faer's kernel choice. Test-only: the two numeric
    /// kernels are independent implementations and both have to be exercised.
    #[cfg(test)]
    pub(crate) fn new_with_threshold(
        n: usize,
        edges: &[(usize, usize)],
        threshold: SupernodalThreshold,
    ) -> Self {
        Self::build(n, edges, threshold)
    }

    fn build(n: usize, edges: &[(usize, usize)], threshold: SupernodalThreshold) -> Self {
        let m = edges.len();

        // Upper triangle in CSC: the entry for edge (a, b) lives in column
        // max(a, b) at row min(a, b); the diagonal of column c lives at row c and
        // is therefore always the *last* entry of that column.
        let mut col_ptr = vec![0usize; n + 1];
        for &(a, b) in edges {
            assert!(a != b, "cover edge must not be a self-loop");
            assert!(a < n && b < n, "cover edge endpoint out of range");
            col_ptr[a.max(b) + 1] += 1;
        }
        // One extra slot per column for the diagonal.
        for c in 0..n {
            col_ptr[c + 1] += col_ptr[c] + 1;
        }
        debug_assert_eq!(col_ptr[n], n + m);

        // Bucket the edges by their owning column, then sort each bucket by row so
        // that the CSC row indices come out ascending.
        let mut order: Vec<usize> = (0..m).collect();
        order.sort_by_key(|&r| {
            let (a, b) = edges[r];
            (a.max(b), a.min(b))
        });

        let mut row_idx = vec![0usize; n + m];
        let mut diag_slot = vec![0usize; n];
        let mut edge_slot = vec![0usize; m];
        for c in 0..n {
            diag_slot[c] = col_ptr[c + 1] - 1;
            row_idx[diag_slot[c]] = c;
        }
        let mut fill = col_ptr[..n].to_vec();
        for &r in &order {
            let (a, b) = edges[r];
            let (lo, hi) = (a.min(b), a.max(b));
            let slot = fill[hi];
            fill[hi] += 1;
            debug_assert!(slot < diag_slot[hi]);
            row_idx[slot] = lo;
            edge_slot[r] = slot;
        }
        debug_assert!((0..n).all(|c| fill[c] == diag_slot[c]));
        debug_assert!(
            (0..n).all(|c| row_idx[col_ptr[c]..col_ptr[c + 1]]
                .windows(2)
                .all(|w| w[0] < w[1])),
            "row indices must be strictly ascending within a column",
        );

        let pattern =
            SymbolicSparseColMatRef::<'_, usize>::new_checked(n, n, &col_ptr, None, &row_idx);
        let symbolic = factorize_symbolic_cholesky(
            pattern,
            Side::Upper,
            SymmetricOrdering::Amd,
            CholeskySymbolicParams {
                supernodal_flop_ratio_threshold: threshold,
                ..Default::default()
            },
        )
        .expect("symbolic Cholesky of the ADMM KKT matrix");

        let len_val = symbolic.len_val();
        let req = StackReq::any_of(&[
            symbolic
                .factorize_numeric_llt_scratch::<f64>(Par::Seq, Spec::<LltParams, f64>::default()),
            symbolic
                .factorize_numeric_llt_scratch::<f32>(Par::Seq, Spec::<LltParams, f32>::default()),
            symbolic.solve_in_place_scratch::<f64>(1, Par::Seq),
            symbolic.solve_in_place_scratch::<f32>(1, Par::Seq),
        ]);

        Self {
            n,
            m,
            col_ptr,
            row_idx,
            diag_slot,
            edge_slot,
            symbolic,
            k_val_f64: vec![0.0f64; n + m],
            k_val_f32: vec![0.0f32; n + m],
            l_val_f64: vec![0.0f64; len_val],
            l_val_f32: vec![0.0f32; len_val],
            rhs_f32: RefCell::new(vec![0.0f32; n]),
            stack: RefCell::new(MemBuffer::new(req)),
            use_f32: false,
            factored: false,
        }
    }

    /// Numeric phase: rebuild the values of `K` for the given `sigma`, `rho` and
    /// `coef`, then refactorize into the preallocated `L` storage.
    ///
    /// `O(n + m)` assembly followed by one numeric factorization over the symbolic
    /// structure computed in [`Self::new`]. `use_f32` chooses the precision the
    /// factor is stored in; [`Self::solve_in_place`] follows that choice.
    ///
    /// Allocates nothing.
    ///
    /// # Panics
    ///
    /// Panics if `edges` or `coef` do not match the length and the pairs handed to
    /// [`Self::new`] (the latter only under `debug_assertions`).
    pub(crate) fn refactor(
        &mut self,
        sigma: f64,
        rho: f64,
        edges: &[(usize, usize)],
        coef: &[(f64, f64)],
        use_f32: bool,
    ) {
        assert_eq!(
            edges.len(),
            self.m,
            "edge count changed since the symbolic phase"
        );
        assert_eq!(coef.len(), self.m, "coef must be parallel to edges");

        let Self {
            n,
            col_ptr,
            row_idx,
            diag_slot,
            edge_slot,
            symbolic,
            k_val_f64,
            k_val_f32,
            l_val_f64,
            l_val_f32,
            stack,
            ..
        } = self;

        // Assemble the upper triangle of K in f64 regardless of the requested
        // output precision: the diagonal is a sum of up to `degree` terms and is
        // worth accumulating in the wider type before narrowing once.
        let base = 1.0 + sigma;
        for c in 0..*n {
            k_val_f64[diag_slot[c]] = base;
        }
        for (r, (&(a, b), &(ca, cb))) in edges.iter().zip(coef).enumerate() {
            debug_assert_eq!(row_idx[edge_slot[r]], a.min(b));
            k_val_f64[diag_slot[a]] += rho * ca * ca;
            k_val_f64[diag_slot[b]] += rho * cb * cb;
            k_val_f64[edge_slot[r]] = rho * ca * cb;
        }

        let mut mem = stack.borrow_mut();
        let mem_stack = MemStack::new(&mut mem);
        if use_f32 {
            for (dst, src) in k_val_f32.iter_mut().zip(k_val_f64.iter()) {
                *dst = *src as f32;
            }
            factorize(symbolic, l_val_f32, col_ptr, row_idx, k_val_f32, mem_stack);
        } else {
            factorize(symbolic, l_val_f64, col_ptr, row_idx, k_val_f64, mem_stack);
        }
        drop(mem);

        self.use_f32 = use_f32;
        self.factored = true;
    }

    /// Solve `K x = rhs`, overwriting `rhs` with `x`.
    ///
    /// `rhs` is `f64` on entry and on exit; when the stored factor is `f32` the
    /// narrowing and widening happen inside.
    ///
    /// Allocates nothing.
    ///
    /// # Panics
    ///
    /// Panics if `rhs.len() != n`, or if [`Self::refactor`] has not run yet.
    pub(crate) fn solve_in_place(&self, rhs: &mut [f64]) {
        assert_eq!(rhs.len(), self.n, "right-hand side has the wrong length");
        assert!(self.factored, "solve_in_place before refactor");

        let mut mem = self.stack.borrow_mut();
        let mem_stack = MemStack::new(&mut mem);
        if self.use_f32 {
            let mut narrowed = self.rhs_f32.borrow_mut();
            for (dst, src) in narrowed.iter_mut().zip(rhs.iter()) {
                *dst = *src as f32;
            }
            solve(&self.symbolic, &self.l_val_f32, &mut narrowed, mem_stack);
            for (dst, src) in rhs.iter_mut().zip(narrowed.iter()) {
                *dst = f64::from(*src);
            }
        } else {
            solve(&self.symbolic, &self.l_val_f64, rhs, mem_stack);
        }
    }

    /// Dimension of `K`.
    #[cfg(test)]
    pub(crate) fn dim(&self) -> usize {
        self.n
    }

    /// Which of faer's two numeric kernels the symbolic analysis picked. faer
    /// chooses automatically from the predicted flop ratio, and the two kernels
    /// are separate implementations, so the tests check both.
    #[cfg(test)]
    pub(crate) fn is_supernodal(&self) -> bool {
        matches!(
            self.symbolic.raw(),
            faer::sparse::linalg::cholesky::SymbolicCholeskyRaw::Supernodal(_)
        )
    }
}

/// Numeric `LL^T` over the shared symbolic structure, for either precision.
///
/// Every temporary faer needs comes out of `mem_stack`; nothing here touches the
/// heap.
fn factorize<T: ComplexField>(
    symbolic: &SymbolicCholesky<usize>,
    l_val: &mut [T],
    col_ptr: &[usize],
    row_idx: &[usize],
    k_val: &[T],
    mem_stack: &mut MemStack,
) {
    let n = col_ptr.len() - 1;
    let pattern = SymbolicSparseColMatRef::<'_, usize>::new_checked(n, n, col_ptr, None, row_idx);
    let k = SparseColMatRef::<'_, usize, T>::new(pattern, k_val);
    symbolic
        .factorize_numeric_llt(
            l_val,
            k,
            Side::Upper,
            LltRegularization::default(),
            Par::Seq,
            mem_stack,
            Spec::<LltParams, T>::default(),
        )
        .expect("K = (1 + sigma) I + rho A^T A is positive definite by construction");
}

/// Triangular solves against a stored factor, for either precision.
fn solve<T: ComplexField>(
    symbolic: &SymbolicCholesky<usize>,
    l_val: &[T],
    rhs: &mut [T],
    mem_stack: &mut MemStack,
) {
    let n = symbolic.nrows();
    let llt = LltRef::<'_, usize, T>::new(symbolic, l_val);
    let rhs = MatMut::from_column_major_slice_mut(rhs, n, 1);
    llt.solve_in_place_with_conj(Conj::No, rhs, Par::Seq, mem_stack);
}

#[cfg(test)]
mod test {
    use super::*;

    // ---------------------------------------------------------------------
    // Independent dense reference.
    //
    // Deliberately built the long way round: form A explicitly, multiply out
    // A^T A, add (1 + sigma) I, and solve by Gaussian elimination with partial
    // pivoting. Nothing here shares code with the sparse path.
    // ---------------------------------------------------------------------

    /// Dense `K = (1 + sigma) I + rho A^T A`, built from an explicit dense `A`.
    fn dense_kkt(
        n: usize,
        sigma: f64,
        rho: f64,
        edges: &[(usize, usize)],
        coef: &[(f64, f64)],
    ) -> Vec<Vec<f64>> {
        let m = edges.len();
        let mut a = vec![vec![0.0f64; n]; m];
        for r in 0..m {
            a[r][edges[r].0] = coef[r].0;
            a[r][edges[r].1] = coef[r].1;
        }
        let mut k = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                let mut acc = 0.0;
                for row in a.iter() {
                    acc += row[i] * row[j];
                }
                k[i][j] = rho * acc;
            }
            k[i][i] += 1.0 + sigma;
        }
        k
    }

    /// Gaussian elimination with partial pivoting.
    fn dense_solve(k: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
        let n = b.len();
        let mut a: Vec<Vec<f64>> = k.to_vec();
        let mut x = b.to_vec();
        for col in 0..n {
            let pivot = (col..n)
                .max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))
                .unwrap();
            a.swap(col, pivot);
            x.swap(col, pivot);
            let d = a[col][col];
            assert!(d.abs() > 0.0, "dense reference hit a zero pivot");
            let pivot_row = a[col].clone();
            let pivot_rhs = x[col];
            for row in col + 1..n {
                let f = a[row][col] / d;
                if f == 0.0 {
                    continue;
                }
                for (dst, src) in a[row].iter_mut().zip(&pivot_row).skip(col) {
                    *dst -= f * src;
                }
                x[row] -= f * pivot_rhs;
            }
        }
        for col in (0..n).rev() {
            let mut acc = x[col];
            for c in col + 1..n {
                acc -= a[col][c] * x[c];
            }
            x[col] = acc / a[col][col];
        }
        x
    }

    fn dense_matvec(k: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
        k.iter()
            .map(|row| row.iter().zip(x).map(|(a, b)| a * b).sum())
            .collect()
    }

    fn max_norm(v: &[f64]) -> f64 {
        v.iter().fold(0.0f64, |acc, x| acc.max(x.abs()))
    }

    fn rel_err(got: &[f64], want: &[f64]) -> f64 {
        let diff: Vec<f64> = got.iter().zip(want).map(|(a, b)| a - b).collect();
        max_norm(&diff) / max_norm(want).max(f64::MIN_POSITIVE)
    }

    /// Largest eigenvalue of a dense SPD matrix, by power iteration.
    fn lambda_max(k: &[Vec<f64>]) -> f64 {
        let n = k.len();
        let mut v = vec![1.0f64; n];
        for (i, entry) in v.iter_mut().enumerate() {
            *entry = 1.0 + ((i as f64) * 0.37) % 1.0;
        }
        let mut lambda = 0.0;
        for _ in 0..4000 {
            let w = dense_matvec(k, &v);
            let norm = max_norm(&w);
            lambda = v.iter().zip(&w).map(|(a, b)| a * b).sum::<f64>()
                / v.iter().map(|a| a * a).sum::<f64>();
            for (dst, src) in v.iter_mut().zip(&w) {
                *dst = src / norm;
            }
        }
        lambda
    }

    /// Deterministic 64-bit LCG, so the "random" tests are reproducible.
    struct Lcg(u64);

    impl Lcg {
        fn next_u64(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            self.0
        }

        /// Uniform in `[-1, 1)`.
        fn signed_unit(&mut self) -> f64 {
            2.0 * ((self.next_u64() >> 11) as f64 / (1u64 << 53) as f64) - 1.0
        }

        fn below(&mut self, bound: usize) -> usize {
            (self.next_u64() >> 33) as usize % bound
        }
    }

    /// A connected graph on `n` nodes: a random spanning tree plus `extra` chords,
    /// with unordered pairs unique (as a transitive reduction guarantees).
    fn random_connected(n: usize, extra: usize, rng: &mut Lcg) -> Vec<(usize, usize)> {
        let mut seen = vec![false; n * n];
        let mut edges = Vec::new();
        for child in 1..n {
            let parent = rng.below(child);
            seen[parent * n + child] = true;
            seen[child * n + parent] = true;
            // Orientation deliberately mixed: the pair is unordered, but the API
            // takes an ordered tuple and must handle both.
            if rng.below(2) == 0 {
                edges.push((parent, child));
            } else {
                edges.push((child, parent));
            }
        }
        let mut added = 0;
        let mut guard = 0;
        while added < extra && guard < 100 * extra + 100 {
            guard += 1;
            let a = rng.below(n);
            let b = rng.below(n);
            if a == b || seen[a * n + b] {
                continue;
            }
            seen[a * n + b] = true;
            seen[b * n + a] = true;
            edges.push((a, b));
            added += 1;
        }
        edges
    }

    /// Solve through `KktFactor`, and through the dense reference, and return
    /// both plus three error measures.
    struct Checked {
        sparse: Vec<f64>,
        dense: Vec<f64>,
        /// `||x_sparse - x_dense||_inf / ||x_dense||_inf`. Bounded by
        /// `~ cond(K) * eps`, so it grows with `rho`.
        forward: f64,
        /// `||K x - b||_inf / ||b||_inf`. Also `rho`-sensitive, because `||K||`
        /// grows with `rho` while `||b||` does not.
        residual: f64,
        /// `||K x - b||_inf / (||K||_inf ||x||_inf + ||b||_inf)`, the normwise
        /// relative backward error. Governed by the working precision alone, so
        /// this is the measure that stays flat across condition numbers.
        backward: f64,
    }

    fn check(
        factor: &KktFactor,
        n: usize,
        sigma: f64,
        rho: f64,
        edges: &[(usize, usize)],
        coef: &[(f64, f64)],
        b: &[f64],
    ) -> Checked {
        let k = dense_kkt(n, sigma, rho, edges, coef);
        let dense = dense_solve(&k, b);
        let mut sparse = b.to_vec();
        factor.solve_in_place(&mut sparse);
        let residual: Vec<f64> = dense_matvec(&k, &sparse)
            .iter()
            .zip(b)
            .map(|(a, b)| a - b)
            .collect();
        let k_norm = k
            .iter()
            .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
            .fold(0.0f64, f64::max);
        let scale = k_norm * max_norm(&sparse) + max_norm(b);
        Checked {
            forward: rel_err(&sparse, &dense),
            residual: max_norm(&residual) / max_norm(b).max(f64::MIN_POSITIVE),
            backward: max_norm(&residual) / scale.max(f64::MIN_POSITIVE),
            sparse,
            dense,
        }
    }

    // ---------------------------------------------------------------------
    // Tests
    // ---------------------------------------------------------------------

    #[test]
    fn single_node_no_edges() {
        let edges: Vec<(usize, usize)> = Vec::new();
        let coef: Vec<(f64, f64)> = Vec::new();
        let mut factor = KktFactor::new(1, &edges);
        assert_eq!(factor.dim(), 1);

        // K = [1 + sigma], so x = b / (1 + sigma).
        for &sigma in &[0.0, 0.25, 7.0] {
            factor.refactor(sigma, 3.0, &edges, &coef, false);
            let mut x = vec![5.0];
            factor.solve_in_place(&mut x);
            let want = 5.0 / (1.0 + sigma);
            assert!(
                (x[0] - want).abs() <= 1e-15 * want.abs(),
                "sigma = {sigma}: got {}, want {want}",
                x[0]
            );
        }
    }

    #[test]
    fn two_nodes_one_edge_hand_derived() {
        // A = [1, -1], sigma = 0, rho = 1
        //   =>  K = [[2, -1], [-1, 2]],  det K = 3,  K^-1 = 1/3 [[2, 1], [1, 2]]
        // b = (1, 0)  =>  x = (2/3, 1/3)
        let edges = [(0usize, 1usize)];
        let coef = [(1.0f64, -1.0f64)];
        let mut factor = KktFactor::new(2, &edges);
        factor.refactor(0.0, 1.0, &edges, &coef, false);

        let mut x = vec![1.0, 0.0];
        factor.solve_in_place(&mut x);
        assert!(
            (x[0] - 2.0 / 3.0).abs() < 1e-15 && (x[1] - 1.0 / 3.0).abs() < 1e-15,
            "got {x:?}, want [2/3, 1/3]"
        );

        // b = (1, 1) is in the null space of A, so x = b / (1 + sigma) = b.
        let mut x = vec![1.0, 1.0];
        factor.solve_in_place(&mut x);
        assert!(
            (x[0] - 1.0).abs() < 1e-15 && (x[1] - 1.0).abs() < 1e-15,
            "got {x:?}, want [1, 1]"
        );

        // Reversed edge orientation must give the identical matrix.
        let edges_rev = [(1usize, 0usize)];
        let coef_rev = [(-1.0f64, 1.0f64)];
        let mut rev = KktFactor::new(2, &edges_rev);
        rev.refactor(0.0, 1.0, &edges_rev, &coef_rev, false);
        let mut y = vec![1.0, 0.0];
        rev.solve_in_place(&mut y);
        assert!((y[0] - 2.0 / 3.0).abs() < 1e-15 && (y[1] - 1.0 / 3.0).abs() < 1e-15);

        // sigma != 0: K = [[1 + s + r, -r], [-r, 1 + s + r]] with s = 0.5, r = 2.
        //   d = 3.5, off = -2, det = 3.5^2 - 4 = 8.25
        //   b = (1, 0) => x = (3.5, 2) / 8.25
        factor.refactor(0.5, 2.0, &edges, &coef, false);
        let mut x = vec![1.0, 0.0];
        factor.solve_in_place(&mut x);
        let det = 3.5f64 * 3.5 - 4.0;
        assert!(
            (x[0] - 3.5 / det).abs() < 1e-14 && (x[1] - 2.0 / det).abs() < 1e-14,
            "got {x:?}"
        );
    }

    #[test]
    fn chain_of_five() {
        let n = 5;
        let edges = [(0usize, 1usize), (1, 2), (2, 3), (3, 4)];
        // Unequal coefficients: energy-coordinate scaling makes the two nonzeros
        // of a row differ in magnitude, not just in sign.
        let coef = [(1.0f64, -1.0f64), (2.0, -0.5), (0.25, -3.0), (-1.5, 0.75)];
        let mut factor = KktFactor::new(n, &edges);
        let b = [1.0, -2.0, 0.5, 3.0, -0.25];

        for &(sigma, rho) in &[(0.0, 1.0), (0.0, 1e-3), (0.3, 0.1), (2.0, 100.0)] {
            factor.refactor(sigma, rho, &edges, &coef, false);
            let out = check(&factor, n, sigma, rho, &edges, &coef, &b);
            assert!(
                out.forward < 1e-11,
                "sigma = {sigma}, rho = {rho}: forward {:.3e}, sparse {:?}, dense {:?}",
                out.forward,
                out.sparse,
                out.dense
            );
            assert!(out.residual < 1e-12, "residual {:.3e}", out.residual);
        }
    }

    #[test]
    fn random_connected_graph_f64() {
        let n = 40;
        let mut rng = Lcg(0x5eed_1234_abcd_0001);
        let edges = random_connected(n, 25, &mut rng);
        // Unequal, sign-opposed coefficients drawn per row.
        let coef: Vec<(f64, f64)> = (0..edges.len())
            .map(|_| {
                let a = 0.25 + rng.signed_unit().abs() * 3.0;
                let b = 0.25 + rng.signed_unit().abs() * 3.0;
                (a, -b)
            })
            .collect();
        let b: Vec<f64> = (0..n).map(|_| rng.signed_unit()).collect();

        let mut factor = KktFactor::new(n, &edges);
        for &rho in &[1e-3f64, 0.1, 1.0, 100.0] {
            for &sigma in &[0.0f64, 1.0] {
                factor.refactor(sigma, rho, &edges, &coef, false);
                let out = check(&factor, n, sigma, rho, &edges, &coef, &b);
                assert!(
                    out.forward < 1e-9,
                    "rho = {rho}, sigma = {sigma}: forward {:.3e}",
                    out.forward
                );
                assert!(
                    out.residual < 1e-12,
                    "rho = {rho}, sigma = {sigma}: residual {:.3e}",
                    out.residual
                );
            }
        }
    }

    #[test]
    fn refactor_is_reusable() {
        let n = 30;
        let mut rng = Lcg(0x5eed_1234_abcd_0002);
        let edges = random_connected(n, 18, &mut rng);
        let coef: Vec<(f64, f64)> = (0..edges.len())
            .map(|_| {
                (
                    0.5 + rng.signed_unit().abs(),
                    -(0.5 + rng.signed_unit().abs()),
                )
            })
            .collect();
        let b: Vec<f64> = (0..n).map(|_| rng.signed_unit()).collect();

        let mut reused = KktFactor::new(n, &edges);
        let rhos = [1e-3f64, 1.0, 1e4, 0.1, 1e-3, 250.0];

        // Also assert that our own buffers are never reallocated: same pointer and
        // same length before and after every refactor / solve.
        let fingerprint = |f: &KktFactor| {
            (
                f.k_val_f64.as_ptr(),
                f.k_val_f64.len(),
                f.k_val_f32.as_ptr(),
                f.k_val_f32.len(),
                f.l_val_f64.as_ptr(),
                f.l_val_f64.len(),
                f.l_val_f32.as_ptr(),
                f.l_val_f32.len(),
                f.col_ptr.as_ptr(),
                f.row_idx.as_ptr(),
                f.rhs_f32.borrow().as_ptr(),
                f.rhs_f32.borrow().len(),
            )
        };
        let before = fingerprint(&reused);

        for (i, &rho) in rhos.iter().enumerate() {
            let sigma = 0.1 * i as f64;
            reused.refactor(sigma, rho, &edges, &coef, false);
            let mut from_reused = b.clone();
            reused.solve_in_place(&mut from_reused);

            // A brand-new factor, symbolic phase and all, must agree bit for bit:
            // the symbolic structure genuinely is a function of `edges` alone.
            let mut fresh = KktFactor::new(n, &edges);
            fresh.refactor(sigma, rho, &edges, &coef, false);
            let mut from_fresh = b.clone();
            fresh.solve_in_place(&mut from_fresh);
            assert_eq!(
                from_reused, from_fresh,
                "rho = {rho}: reused factor differs from a fresh one"
            );

            // ... and both must match the dense reference.
            let out = check(&reused, n, sigma, rho, &edges, &coef, &b);
            assert!(
                out.forward < 1e-8,
                "rho = {rho}: forward {:.3e}",
                out.forward
            );
        }
        assert_eq!(before, fingerprint(&reused), "a buffer was reallocated");

        // Switching precision back and forth must keep working.
        reused.refactor(0.0, 1.0, &edges, &coef, true);
        let mut x32 = b.clone();
        reused.solve_in_place(&mut x32);
        reused.refactor(0.0, 1.0, &edges, &coef, false);
        let mut x64 = b.clone();
        reused.solve_in_place(&mut x64);
        assert!(
            rel_err(&x32, &x64) < 1e-4,
            "f32 <-> f64 round trip diverged"
        );
        assert_eq!(before, fingerprint(&reused), "a buffer was reallocated");
    }

    #[test]
    fn unequal_coefficients_scaled_rows() {
        // Row scaling: A' = D A for diagonal D changes A^T A to A^T D^2 A, which
        // is exactly `rho -> rho * d_r^2` per row. Check the assembly handles the
        // resulting wildly unequal magnitudes.
        let n = 12;
        let mut rng = Lcg(0x5eed_1234_abcd_0003);
        let edges = random_connected(n, 6, &mut rng);
        let coef: Vec<(f64, f64)> = (0..edges.len())
            .map(|r| {
                let scale = 10f64.powi((r % 7) as i32 - 3);
                (scale, -scale * 0.75)
            })
            .collect();
        let b: Vec<f64> = (0..n).map(|_| rng.signed_unit()).collect();

        let mut factor = KktFactor::new(n, &edges);
        for &rho in &[1e-2f64, 1.0] {
            factor.refactor(0.0, rho, &edges, &coef, false);
            let out = check(&factor, n, 0.0, rho, &edges, &coef, &b);
            // `||K||_inf` reaches ~1e6 here, so judge the solve by its backward
            // error rather than by a residual scaled against `||b||_inf = O(1)`.
            assert!(
                out.backward < 1e-14,
                "rho = {rho}: backward {:.3e} (residual {:.3e})",
                out.backward,
                out.residual
            );
        }
    }

    /// faer has two independent numeric kernels and picks between them from a
    /// predicted flop ratio. Cover parts of the space a real fit could land in by
    /// pinning each kernel and checking they agree with each other and with the
    /// dense reference.
    #[test]
    fn both_faer_kernels_agree() {
        let n = 50;
        let mut rng = Lcg(0x5eed_1234_abcd_0005);
        // Dense-ish graph so the supernodal kernel actually forms real supernodes.
        let mut edges = Vec::new();
        let mut coef = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                edges.push((i, j));
                coef.push((1.0 + 0.5 * ((i + j) % 3) as f64, -1.0));
            }
        }
        let b: Vec<f64> = (0..n).map(|_| rng.signed_unit()).collect();

        let mut simplicial =
            KktFactor::new_with_threshold(n, &edges, SupernodalThreshold::FORCE_SIMPLICIAL);
        let mut supernodal =
            KktFactor::new_with_threshold(n, &edges, SupernodalThreshold::FORCE_SUPERNODAL);
        assert!(!simplicial.is_supernodal());
        assert!(supernodal.is_supernodal());
        // What the shipping configuration actually picks for this shape.
        println!(
            "  n = {n}, m = {}: auto-selected kernel is {}",
            edges.len(),
            if KktFactor::new(n, &edges).is_supernodal() {
                "supernodal"
            } else {
                "simplicial"
            }
        );

        for &rho in &[1e-3f64, 1.0, 10.0] {
            for &use_f32 in &[false, true] {
                simplicial.refactor(0.0, rho, &edges, &coef, use_f32);
                supernodal.refactor(0.0, rho, &edges, &coef, use_f32);

                let out_simp = check(&simplicial, n, 0.0, rho, &edges, &coef, &b);
                let out_super = check(&supernodal, n, 0.0, rho, &edges, &coef, &b);

                let tol = if use_f32 { 1e-4 } else { 1e-9 };
                assert!(
                    rel_err(&out_simp.sparse, &out_super.sparse) < tol,
                    "rho = {rho}, f32 = {use_f32}: kernels disagree"
                );
                assert!(
                    out_simp.backward < if use_f32 { 1e-5 } else { 1e-13 },
                    "rho = {rho}, f32 = {use_f32}: simplicial bwd {:.3e}",
                    out_simp.backward
                );
                assert!(
                    out_super.backward < if use_f32 { 1e-5 } else { 1e-13 },
                    "rho = {rho}, f32 = {use_f32}: supernodal bwd {:.3e}",
                    out_super.backward
                );
            }
        }
    }

    #[test]
    fn dense_reference_is_independent() {
        // Guard the reference itself: for the 2x2 hand case it must reproduce the
        // closed form.
        let k = dense_kkt(2, 0.0, 1.0, &[(0, 1)], &[(1.0, -1.0)]);
        assert_eq!(k, vec![vec![2.0, -1.0], vec![-1.0, 2.0]]);
        let x = dense_solve(&k, &[1.0, 0.0]);
        assert!((x[0] - 2.0 / 3.0).abs() < 1e-15 && (x[1] - 1.0 / 3.0).abs() < 1e-15);
    }

    /// The shipping decision: how much accuracy does the `f32` factor cost?
    ///
    /// `A` is the +-1 incidence matrix of a connected graph, so `A^T A` is a graph
    /// Laplacian: its null space is spanned by the all-ones vector and
    /// `lambda_min(A^T A) = 0`. With `sigma = 0` that pins
    /// `lambda_min(K) = 1` exactly and `cond(K) = 1 + rho * lambda_max(A^T A)`,
    /// so `rho` is a dial for the condition number.
    #[test]
    fn f32_versus_f64_accuracy() {
        let n = 64;
        let mut rng = Lcg(0x5eed_1234_abcd_0004);
        let edges = random_connected(n, 40, &mut rng);
        let coef: Vec<(f64, f64)> = (0..edges.len()).map(|_| (1.0, -1.0)).collect();
        let b: Vec<f64> = (0..n).map(|_| rng.signed_unit()).collect();

        // lambda_max(A^T A) from the rho = 1, sigma = 0 matrix: K = I + A^T A.
        let lam_max = lambda_max(&dense_kkt(n, 0.0, 1.0, &edges, &coef)) - 1.0;

        let mut factor = KktFactor::new(n, &edges);
        println!(
            "\n  n = {n}, m = {}, lambda_max(A^T A) = {lam_max:.6}",
            edges.len()
        );
        println!(
            "  {:>9}  {:>11}  {:>11}  {:>11}  {:>11}  {:>11}  {:>11}  {:>11}",
            "cond(K)", "rho", "f64 fwd", "f32 fwd", "f64 resid", "f32 resid", "f64 bwd", "f32 bwd"
        );

        for &target in &[1e2f64, 1e4, 1e6] {
            let rho = (target - 1.0) / lam_max;
            let cond = 1.0 + rho * lam_max;

            factor.refactor(0.0, rho, &edges, &coef, false);
            let out64 = check(&factor, n, 0.0, rho, &edges, &coef, &b);

            factor.refactor(0.0, rho, &edges, &coef, true);
            let out32 = check(&factor, n, 0.0, rho, &edges, &coef, &b);

            println!(
                "  {cond:>9.3e}  {rho:>11.4e}  {:>11.3e}  {:>11.3e}  {:>11.3e}  {:>11.3e}  {:>11.3e}  {:>11.3e}",
                out64.forward,
                out32.forward,
                out64.residual,
                out32.residual,
                out64.backward,
                out32.backward,
            );

            // f64: forward error tracks cond * 2^-53, with slack.
            assert!(
                out64.forward < 100.0 * cond * f64::EPSILON,
                "cond {cond:.2e}: f64 forward error {:.3e} exceeds the model",
                out64.forward
            );
            // f32: forward error tracks cond * 2^-24, with slack. This is the
            // number that decides whether the f32 factor is shippable.
            assert!(
                out32.forward < 100.0 * cond * f32::EPSILON as f64,
                "cond {cond:.2e}: f32 forward error {:.3e} exceeds the model",
                out32.forward
            );
            // The backward error is what a stable factorization actually promises,
            // and unlike the two columns above it is flat in cond: it tracks the
            // working precision times a small growth factor.
            assert!(
                out32.backward < 100.0 * f32::EPSILON as f64,
                "cond {cond:.2e}: f32 backward error {:.3e}",
                out32.backward
            );
            assert!(
                out64.backward < 100.0 * f64::EPSILON,
                "cond {cond:.2e}: f64 backward error {:.3e}",
                out64.backward
            );
        }
    }
}
