//! The constraint matrix `A` of the partial-order quadratic program.

/// `A` stored by row, as the cover edges plus one coefficient pair per row.
///
/// Every row of `A` encodes one cover edge `(i, j)` and carries exactly two nonzeros,
/// of opposite sign, in columns `i` and `j`. A general sparse format would spend an
/// index array recovering what the edge list already says, so the solver works
/// directly from the edge list: `coef[r].0` is the entry in column `edges[r].0` and
/// `coef[r].1` the entry in column `edges[r].1`.
///
/// The edges come from the transitive reduction, so the unordered pairs `{i, j}` are
/// distinct across rows and `i != j`.
pub(crate) struct Constraints<'a> {
    pub(crate) edges: &'a [(usize, usize)],
    pub(crate) coef: &'a [(f64, f64)],
    pub(crate) n: usize,
}

impl Constraints<'_> {
    pub(crate) fn m(&self) -> usize {
        self.edges.len()
    }

    /// `out <- A u`.
    pub(crate) fn mul(&self, u: &[f64], out: &mut [f64]) {
        debug_assert_eq!(u.len(), self.n);
        debug_assert_eq!(out.len(), self.m());
        for ((out_r, &(i, j)), &(a, b)) in out.iter_mut().zip(self.edges).zip(self.coef) {
            *out_r = a * u[i] + b * u[j];
        }
    }

    /// `out <- A^T y`.
    pub(crate) fn mul_transpose(&self, y: &[f64], out: &mut [f64]) {
        debug_assert_eq!(y.len(), self.m());
        debug_assert_eq!(out.len(), self.n);
        out.fill(0.0);
        for ((&(i, j), &(a, b)), &y_r) in self.edges.iter().zip(self.coef).zip(y) {
            out[i] += a * y_r;
            out[j] += b * y_r;
        }
    }

    /// A rigorous upper bound on `lambda_max(A^T A)`, via `||A||_2^2 <= ||A||_1 ||A||_inf`.
    ///
    /// The ADMM system is `K = (1 + sigma) I + rho A^T A`, and `P = I` puts a floor of
    /// `1 + sigma` under its spectrum, so `cond(K) <= 1 + rho * lambda_max(A^T A) / (1 + sigma)`.
    /// That is what decides whether a given `rho` may be factorized in f32
    /// (see `kkt::KktFactor`), so the bound has to be an upper bound, not an estimate.
    ///
    /// `scratch` is an `n`-length buffer; its contents on entry are ignored.
    pub(crate) fn spectral_norm_squared_bound(&self, scratch: &mut [f64]) -> f64 {
        debug_assert_eq!(scratch.len(), self.n);
        scratch.fill(0.0);
        let mut norm_inf = 0.0f64;
        for (&(i, j), &(a, b)) in self.edges.iter().zip(self.coef) {
            norm_inf = norm_inf.max(a.abs() + b.abs());
            scratch[i] += a.abs();
            scratch[j] += b.abs();
        }
        let norm_1 = scratch.iter().copied().fold(0.0f64, f64::max);
        norm_1 * norm_inf
    }
}

#[cfg(test)]
mod test {
    use super::Constraints;

    /// Two nodes, one edge `x_1 - x_0 >= 0`.
    #[test]
    fn single_edge_products() {
        let edges = [(0usize, 1usize)];
        let coef = [(-1.0, 1.0)];
        let a = Constraints {
            edges: &edges,
            coef: &coef,
            n: 2,
        };

        let mut au = [0.0; 1];
        a.mul(&[3.0, 5.0], &mut au);
        assert_eq!(au, [2.0]);

        let mut aty = [0.0; 2];
        a.mul_transpose(&[7.0], &mut aty);
        assert_eq!(aty, [-7.0, 7.0]);
    }

    /// `A^T y` must accumulate across every row touching a column, not overwrite.
    #[test]
    fn transpose_accumulates_over_shared_columns() {
        // A star: node 0 below nodes 1, 2, 3.
        let edges = [(0usize, 1usize), (0, 2), (0, 3)];
        let coef = [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)];
        let a = Constraints {
            edges: &edges,
            coef: &coef,
            n: 4,
        };

        let mut aty = [0.0; 4];
        a.mul_transpose(&[1.0, 2.0, 4.0], &mut aty);
        // Column 0 sees -1 from each of the three rows.
        assert_eq!(aty, [-7.0, 1.0, 2.0, 4.0]);
    }

    /// `<Au, y> == <u, A^T y>` for arbitrary coefficients -- the adjoint identity that
    /// makes the two routines a matched pair. Exact in binary floating point here
    /// because every value is a small dyadic rational.
    #[test]
    fn products_are_adjoint() {
        let edges = [(0usize, 2usize), (1, 2), (0, 1)];
        let coef = [(-0.5, 1.0), (-0.25, 1.0), (-1.0, 0.5)];
        let a = Constraints {
            edges: &edges,
            coef: &coef,
            n: 3,
        };

        let u = [1.5, -2.0, 0.75];
        let y = [2.0, -1.0, 0.5];

        let mut au = [0.0; 3];
        a.mul(&u, &mut au);
        let mut aty = [0.0; 3];
        a.mul_transpose(&y, &mut aty);

        let lhs: f64 = au.iter().zip(&y).map(|(p, q)| p * q).sum();
        let rhs: f64 = u.iter().zip(&aty).map(|(p, q)| p * q).sum();
        assert_eq!(lhs, rhs);
    }

    /// The f32 gate depends on this being a genuine *upper* bound on the top
    /// eigenvalue, so check it against a power iteration on `A^T A`.
    #[test]
    fn spectral_bound_dominates_power_iteration() {
        // A chain of 6 nodes with assorted coefficients.
        let edges: Vec<(usize, usize)> = (0..5).map(|i| (i, i + 1)).collect();
        let coef: Vec<(f64, f64)> = [
            (-1.0, 1.0),
            (-0.5, 1.0),
            (-1.0, 0.25),
            (-0.75, 1.0),
            (-1.0, 1.0),
        ]
        .to_vec();
        let n = 6;
        let a = Constraints {
            edges: &edges,
            coef: &coef,
            n,
        };

        let mut scratch = vec![0.0; n];
        let bound = a.spectral_norm_squared_bound(&mut scratch);

        // Power iteration on A^T A.
        let mut v = vec![1.0; n];
        let mut au = vec![0.0; edges.len()];
        let mut atav = vec![0.0; n];
        let mut lambda = 0.0;
        for _ in 0..2000 {
            a.mul(&v, &mut au);
            a.mul_transpose(&au, &mut atav);
            lambda = atav.iter().zip(&v).map(|(p, q)| p * q).sum::<f64>()
                / v.iter().map(|p| p * p).sum::<f64>();
            let norm = atav.iter().map(|p| p * p).sum::<f64>().sqrt();
            if norm == 0.0 {
                break;
            }
            for (dst, src) in v.iter_mut().zip(&atav) {
                *dst = src / norm;
            }
        }
        assert!(
            bound >= lambda,
            "bound {bound} must dominate the true lambda_max {lambda}"
        );
        // And it should not be uselessly loose.
        assert!(
            bound <= 8.0 * lambda,
            "bound {bound} vs lambda_max {lambda}"
        );
    }
}
