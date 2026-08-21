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
//!
//! The deeper specialization is combinatorial: the program is isotonic regression on a
//! DAG, so its optimum is piecewise constant on a partition into pooled components, exact
//! once the active edge set is known. [`polish`] computes that exact minimizer for a
//! candidate active set in `O(n + m)`, peels the candidate's exact multipliers off its
//! spanning forest, and the duality gap of the pair is a rigorous optimality certificate.
//! [`solver`] wraps the splitting iteration in that machinery: a *predictor* tries the
//! previous threshold's active set before touching the factorization -- consecutive
//! thresholds usually share their structure, so most solves finish with zero iterations --
//! and a *corrector* retries at failing residual checks with exponential backoff, cutting
//! off the splitting iteration's slow tail the moment the structure settles. Negative
//! peeled multipliers drive classical one-row active-set exchanges between attempts.
//! ADMM remains the fallback that handles genuinely moved structure.
pub mod constraints;
pub mod hildreth;
pub mod kkt;
pub mod polish;
pub mod solver;
#[cfg(test)]
mod test;
