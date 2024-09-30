use crate::functionals::{ClippingWrapper, KaplanMeier};
use crate::partial_order::routines::{compute_transitive_closure, derive_transitive_reduction};
use crate::partial_order::{
    AlgorithmOutput, BitSet, ExtendedAlgorithmContext, OrderingInfo, QualityIndicators, functionals,
};
use crate::progress::ProgressTracker;
use crate::structures::{Direction, Observation};

#[must_use]
pub fn algorithm<D: Direction>(
    context: &ExtendedAlgorithmContext,
    progress: &dyn ProgressTracker,
) -> AlgorithmOutput {
    if context.n() == 0 {
        return AlgorithmOutput {
            cdfs: Vec::with_capacity(0),
            ordering_info: OrderingInfo::empty(),
            quality_indicators: QualityIndicators {
                precision: 0.0,
                convergence_fraction: 0.0,
            },
        };
    }

    match (context.n() - 1) / BitSet::<1>::capacity() {
        0..1 => algorithm_inner::<D, 1>(context, progress),
        1..2 => algorithm_inner::<D, 2>(context, progress),
        2..4 => algorithm_inner::<D, 4>(context, progress),
        4..8 => algorithm_inner::<D, 8>(context, progress),
        8..16 => algorithm_inner::<D, 16>(context, progress),
        16..32 => algorithm_inner::<D, 32>(context, progress),
        32..64 => algorithm_inner::<D, 64>(context, progress),
        64..128 => algorithm_inner::<D, 128>(context, progress),
        _ => unimplemented!(
            "data set is too large, largest supported is n = {}",
            BitSet::<128>::capacity(),
        ),
    }
}

fn algorithm_inner<D: Direction, const B: usize>(
    context: &ExtendedAlgorithmContext,
    progress: &dyn ProgressTracker,
) -> AlgorithmOutput {
    // Build comparable pairs and reduce to cover edges
    let constraint_edges = derive_transitive_reduction(
        &context.x_unique,
        context.n_covariate(),
        context.covariate_dimension(),
    );

    progress.set_total(context.n_threshold());

    let (successors, predecessors, topological_order, reverse_topological_order) =
        compute_transitive_closure::<B>(context.n_covariate(), &constraint_edges);

    // Build reverse map from covariate to observation index
    let mut data_start_indices = {
        let mut data_counts = vec![0; context.n_covariate()];
        for &i in &context.x {
            data_counts[reverse_topological_order[i]] += 1;
        }
        for i in 1..context.n_covariate() {
            data_counts[i] += data_counts[i - 1];
        }
        data_counts
    };
    let mut data_indices = vec![0; context.n()];
    for (j, &i) in context.x.iter().enumerate() {
        let topological_index = reverse_topological_order[i];
        data_start_indices[topological_index] -= 1;
        data_indices[data_start_indices[topological_index]] = j;
    }

    let indexer = |i| {
        let end = if i < context.n_covariate() - 1 {
            data_start_indices[i + 1]
        } else {
            context.n()
        };
        data_indices[data_start_indices[i]..end]
            .iter()
            .map(|&j| Observation {
                x: (),
                y: context.y[j],
                observed: context.y_observed[j],
                weight: context.weight[j],
            })
    };

    // Inclusive closures for convenience
    let successors_inclusive: Vec<_> = successors
        .into_iter()
        .enumerate()
        .map(|(i, set)| set.with_singleton(i))
        .collect();

    let mut cdfs = Vec::with_capacity(context.n_threshold() * context.n_covariate());

    let mut data_index = 0;
    loop {
        let current_threshold = context.y[data_index];
        while data_index < context.n()
            && (!context.y_observed[data_index] || context.y[data_index] == current_threshold)
        {
            data_index += 1;
        }

        let inner = KaplanMeier::new(current_threshold);
        let functional = ClippingWrapper::new(inner);

        let solution = functionals::algorithm_pre_sorted_inner::<D::REVERSE, _, _, _, _>(
            context.n_covariate(),
            &topological_order,
            &successors_inclusive,
            &predecessors,
            &indexer,
            &functional,
        );

        cdfs.extend(solution);
        progress.increment();

        if data_index == context.n() {
            break;
        }
    }

    AlgorithmOutput {
        cdfs,
        ordering_info: OrderingInfo::from_edges(constraint_edges, context.n_covariate()),
        quality_indicators: QualityIndicators {
            precision: 0.0,
            convergence_fraction: 0.0,
        },
    }
}
