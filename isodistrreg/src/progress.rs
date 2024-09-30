//! Progress reporting for fits.
//!
//! Implementations are passed as a generic parameter to the inner algorithms so that the no-op
//! variant ([`NoProgress`]) is fully monomorphized away, leaving zero runtime cost when progress
//! tracking is disabled.

/// Tracks the progress of a fit. Implementors must be thread-safe so that progress can be
/// reported from parallel subagging workers.
pub trait ProgressTracker: Send + Sync {
    /// Set the total number of work units expected for the fit. Called once by the subagging
    /// wrapper before any inner fits run.
    fn set_total(&self, total: usize);

    /// Mark one work unit (typically: one threshold) as completed.
    fn increment(&self);

    fn finish(&self);
}

/// Zero-cost no-op tracker. Use this when progress reporting is disabled — the trait methods
/// are inlined to nothing, so the compiler can elide every progress-related instruction in the
/// monomorphized fit.
pub struct NoProgress;

impl ProgressTracker for NoProgress {
    #[inline(always)]
    fn set_total(&self, _total: usize) {}

    #[inline(always)]
    fn increment(&self) {}

    #[inline(always)]
    fn finish(&self) {}
}
