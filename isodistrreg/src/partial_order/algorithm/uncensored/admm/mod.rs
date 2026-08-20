//! The quadratic-program solver for the partial-order estimator.
//!
//! Solves `min 0.5 u^T u + q^T u  s.t.  A u >= 0` by operator splitting (Stellato et al.,
//! 2020), specialized to the structure this estimator produces: `P = I` in energy
//! coordinates, and every constraint row carrying exactly two nonzeros of opposite sign at
//! the endpoints of a cover edge.
//!
//! That specialization is what makes an in-tree solver worthwhile rather than merely
//! possible:
//!
//!  * eliminating the slack leaves the normal equations `K = (1 + sigma) I + rho A^T A`,
//!    whose sparsity pattern is exactly the cover graph plus a diagonal. It is assembled in
//!    `O(n + m)` straight from the edge list, so `A^T A` is never formed.
//!  * `P = I` floors the spectrum at `1 + sigma`, making `K` unconditionally positive
//!    definite -- plain Cholesky, no pivoting, no regularization -- and giving
//!    `cond(K) = 1 + rho lambda_max(A^T A)` in closed form, which is what lets the
//!    factorization decide its own working precision.
//!  * the pattern never changes across a threshold sequence, so the symbolic factorization
//!    runs once per fit and later thresholds refresh values only.
pub mod constraints;
pub mod hildreth;
pub mod kkt;
pub mod polish;
pub mod solver;
#[cfg(test)]
mod test;
