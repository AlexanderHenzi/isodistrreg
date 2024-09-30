use crate::error::Error;
use crate::partial_order::prediction::Interpolation;
use crate::partial_order::preprocessing::{
    preprocess_censored, preprocess_uncensored, unify_group_orders, validate,
};
use crate::partial_order::{censored, uncensored};
use crate::routines::{empirical_cdf, kaplan_meier, transpose};
use crate::structures::{Observation, StochasticOrder};
use crate::{
    CovariateInterpolator, Decreasing, Increasing, IntoCdfIterator,
    IsotonicDistributionalRegressionFit, ProgressTracker,
};
use itertools::Itertools;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fmt::Display;
use std::mem;
use std::str::FromStr;

/// A computed IDR solution that can be used to make distributional predictions.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Fit {
    /// Whether the fit is increasing (i.e., nondecreasing) w.r.t. the covariate ordering.
    pub increasing: bool,

    /// For each covariate, the estimated distribution.
    ///
    /// These values will be interpolated when predicting a CDF for an unseen value.
    ///
    /// Represents a covariate-major matrix with in each minor index a CDF value for the
    /// corresponding threshold.
    pub cdfs: Vec<f64>,
    /// Global average CDF, disregarding completely the covariates.
    ///
    /// Used as a prediction for covariates that are incomparable to any of the others.
    pub global_cdf: Vec<f64>,

    /// Which dimensions of the covariate are equivalent?
    ///
    /// Used to transform the newly provided covariate to a quantity of the same dimension with
    /// which ordering can be determined through component-wise comparisons only as a pre-processing
    /// step.
    pub covariate_groups: CovariateGroups,
    /// Unique covariate values.
    ///
    /// Are used to determine the order constraints against unseen data. These have already been
    /// "unified" in terms of the ordering relation on the different groups by `unify_group_orders`,
    /// which other covariates need to be converted by also before prediction.
    ///
    /// Represents a covariate-major matrix.
    pub covariates: Vec<f64>,
    /// Indices of which covariates are below which others.
    ///
    /// Together with the raw covariates, these are used to place unseen covariates in the partial
    /// order graph.
    pub ordering_info: OrderingInfo,

    /// Unique thresholds sorted from low to high.
    pub thresholds: Vec<f64>,
    /// How well the numerical solve worked.
    pub quality_indicators: QualityIndicators,
}

#[derive(Clone)]
pub struct Config {
    /// Settings passed to OSQP.
    pub osqp_settings: osqp::Settings,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            osqp_settings: osqp::Settings::default()
                .verbose(false)
                .eps_abs(1e-5)
                .eps_rel(1e-5)
                .max_iter(10_000),
        }
    }
}

/// Data structures needed in the `find_neighbors` method.
///
/// Initialize with a small size and hope it doesn't grow to the number of covariates.
pub struct PredictionWorkspace {
    pub stack: Vec<NodeID>,
    pub visited: HashSet<NodeID>,
    pub found: HashSet<NodeID>,
}
impl PredictionWorkspace {
    #[must_use]
    pub fn new() -> Self {
        Self {
            stack: Vec::new(),
            visited: HashSet::new(),
            found: HashSet::new(),
        }
    }
}
impl Default for PredictionWorkspace {
    fn default() -> Self {
        Self::new()
    }
}

impl IsotonicDistributionalRegressionFit for Fit {
    type Covariate<'a> = &'a [f64];
    type CovariateOrder = CovariateGroups;
    type Config = Config;

    /// Compute the multivariate IDR.
    ///
    /// # Arguments
    ///
    /// The number of observations is `n` and their dimensionality is `p`. The number of groups is `g`.
    ///
    /// - `covariates`: A `n x p` matrix containing a row for each observation.
    /// - `responses`: A length-`n` vector.
    /// - `weights`: A length-`n` vector with positive numbers.
    /// - `groups`: A length-`p` vector with the group index (`0`, `1`, ..., `g - 1`) of each covariate
    ///   index.
    /// - `orders`: A length-`g` vector of the orders of the groups.
    ///
    /// # Details
    ///
    /// This implements the multivariate branch of the R code for `stoch == "sd"`:
    /// - Transform X "by groups" according to their `orders` (only `ComponentWise` and
    ///   `StochasticDominance` are required here; `IncreasingConvexOrder` is also supported for
    ///   completeness).
    /// - Aggregate duplicate covariate rows to unique rows; keep, for each unique row j,
    ///   all responses `y_j` and the sum of weights `w_j`.
    /// - Build a partial order using the componentwise order on the transformed (prepared) X
    ///   and reduce it to cover relations (a transitive reduction) to keep the constraint
    ///   set small.
    /// - For each threshold z in `unique(responses)` except the last, solve the QP
    ///   min 1/2 x' P x + q' x s.t. A x >= 0
    ///   where `P = diag(w_j)` and `q_j = -w_j * F_j(z)` with `F_j(z)` the empirical CDF of the
    ///   responses attached to row j at threshold z. Warm start and cost updates are used
    ///   to accelerate repeated solves.
    /// - Clamp to [0,1], record whether the solver hit the iteration limit, compute a diagnostic
    ///   "precision" (maximal downward step across thresholds), run a final PAVA along thresholds
    ///   to ensure monotonicity in z, and append a trailing column of ones.
    ///
    /// # Returns
    ///
    /// Returns `N x K_plus_1` CDFs (`K_plus_1` = number of unique thresholds), `precision`,
    /// and the fraction of thresholds with `MaxIterationsReached`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// // Two exchangeable covariates in one group using SD order (sort rows decreasing).
    /// // We create four observations, two unique covariate rows after SD-preparation.
    /// use isodistrreg::partial_order::{CovariateGroups, Config, Fit, PartialOrder};
    /// use isodistrreg::{Increasing, IsotonicDistributionalRegressionFit, NoProgress, StochasticOrder};
    ///
    /// let covariates = vec![
    ///     0.2, 0.6, // row 0
    ///     0.6, 0.2, // row 1  -> same multiset as row 0, SD-sorted both become [0.6,0.2]
    ///     0.9, 0.8, // row 2
    ///     0.8, 0.9, // row 3  -> same multiset as row 2, SD-sorted both become [0.9,0.8]
    /// ];
    /// let responses = vec![0.0, 1.0, 0.5, 0.6];
    /// let weights   = vec![1.0; 4];
    /// // both predictors in the same group
    /// let covariate_order = CovariateGroups::parse([("sd", [0, 1])], 2).unwrap();
    ///
    /// let fit = Fit::fit(
    ///     &covariates,
    ///     &responses,
    ///     None,
    ///     Some(weights).as_deref(),
    ///     covariate_order,
    ///     StochasticOrder::StochasticDominance,
    ///     false,
    ///     Config::default(),
    ///     &NoProgress,
    /// ).unwrap();
    ///
    /// // Two unique covariate rows remain; CDF columns equals unique thresholds count.
    /// assert_eq!(fit.covariates.len(), 2 * 2);
    /// assert_eq!(fit.cdfs.len(), 2 * 4);
    ///
    /// let expected = [
    ///     0.5, 0.5, 0.75, 1.0,
    ///     0.0, 0.5, 0.75, 1.0,
    /// ];
    /// for (realized, expected) in fit.cdfs.into_iter().zip(expected.into_iter()) {
    ///     assert!((realized - expected).abs() < 1e-3);
    /// }
    /// ```
    fn fit(
        x: &[f64],
        y: &[f64],
        y_observed: Option<&[bool]>,
        sample_weight: Option<&[f64]>,
        covariate_order: Self::CovariateOrder,
        response_order: StochasticOrder,
        decreasing: bool,
        config: Self::Config,
        progress: &dyn ProgressTracker,
    ) -> Result<Self, Error> {
        let n = validate(x, y, y_observed, sample_weight, &covariate_order)?;

        let mut weight_allocation = None;
        let weight_to_use = sample_weight.unwrap_or_else(|| {
            weight_allocation = Some(vec![1.0; n]);
            weight_allocation.as_deref().unwrap()
        });

        let uncensored_case = |covariate_order| -> Self {
            let mut algorithm_context =
                preprocess_uncensored(x, y, weight_to_use, &covariate_order);
            let algo_result = match response_order {
                StochasticOrder::StochasticDominance => match decreasing {
                    false => uncensored::<Increasing, false>(&algorithm_context, config, progress),
                    true => uncensored::<Decreasing, false>(&algorithm_context, config, progress),
                },
                StochasticOrder::HazardRateOrder => match decreasing {
                    false => uncensored::<Increasing, true>(&algorithm_context, config, progress),
                    true => uncensored::<Decreasing, true>(&algorithm_context, config, progress),
                },
            };
            let covariates = mem::take(&mut algorithm_context.x);
            let global_cdf = empirical_cdf(
                algorithm_context
                    .y
                    .iter()
                    .zip(algorithm_context.weight.iter())
                    .chunk_by(|&(&response, _)| response)
                    .into_iter()
                    .map(|(response, group)| Observation {
                        x: (),
                        y: response,
                        observed: (),
                        weight: group.into_iter().map(|(_, w)| w).sum::<f64>(),
                    }),
                algorithm_context.weight.iter().sum(),
            );
            let output = Self {
                increasing: !decreasing,
                cdfs: algo_result.cdfs,
                global_cdf,

                covariate_groups: covariate_order,
                covariates,
                ordering_info: algo_result.ordering_info,

                thresholds: algorithm_context.thresholds,
                quality_indicators: algo_result.quality_indicators,
            };
            debug_assert!({
                output.assert_consistent();
                true
            });
            output
        };

        match y_observed {
            None => Ok(uncensored_case(covariate_order)),
            Some(indicators) if indicators.iter().all(|&b| !b) => {
                Ok(uncensored_case(covariate_order))
            }
            Some(indicators) => match response_order {
                StochasticOrder::StochasticDominance => {
                    let algorithm_context =
                        preprocess_censored(x, y, indicators, weight_to_use, &covariate_order);
                    if indicators.iter().all(|&b| b) {
                        let empty = Fit {
                            increasing: !decreasing,
                            cdfs: Vec::with_capacity(0),
                            global_cdf: Vec::with_capacity(0),

                            covariate_groups: covariate_order,
                            covariates: Vec::with_capacity(0),
                            ordering_info: OrderingInfo::empty(),

                            thresholds: Vec::with_capacity(0),
                            quality_indicators: QualityIndicators {
                                precision: 0.0,
                                convergence_fraction: 0.0,
                            },
                        };
                        Ok(empty)
                    } else {
                        let mut result = match decreasing {
                            false => censored::<Increasing>(&algorithm_context, progress),
                            true => censored::<Decreasing>(&algorithm_context, progress),
                        };
                        transpose(
                            &mut result.cdfs,
                            algorithm_context.n_threshold(),
                            algorithm_context.n_covariate(),
                        );
                        let global_cdf = kaplan_meier(
                            (0..algorithm_context.n()).map(|i| Observation {
                                x: (),
                                y: algorithm_context.y[i],
                                observed: algorithm_context.y_observed[i],
                                weight: algorithm_context.weight[i],
                            }),
                            algorithm_context.weight.iter().sum(),
                        );

                        let output = Fit {
                            increasing: !decreasing,
                            cdfs: result.cdfs,
                            global_cdf,

                            covariate_groups: covariate_order,
                            covariates: algorithm_context.x_unique,
                            ordering_info: result.ordering_info,

                            thresholds: algorithm_context.thresholds,
                            quality_indicators: result.quality_indicators,
                        };
                        debug_assert!({
                            output.assert_consistent();
                            true
                        });
                        Ok(output)
                    }
                }
                StochasticOrder::HazardRateOrder => Err(Error::NotImplemented(
                    "censoring is not supported with partially ordered covariates under the hazard order constraint",
                )),
            },
        }
    }

    fn interpolate_covariate<'a>(
        &'a self,
        covariate: Self::Covariate<'_>,
    ) -> impl CovariateInterpolator + IntoCdfIterator + 'a {
        let target: Vec<_> = unify_group_orders(covariate, &self.covariate_groups).collect();
        Interpolation::new(
            &target,
            &self.covariates,
            &self.ordering_info,
            self.increasing,
            (&self.cdfs, &self.global_cdf),
            &mut PredictionWorkspace::new(),
        )
    }

    fn thresholds(&self) -> &[f64] {
        &self.thresholds
    }

    fn assert_consistent(&self) {
        assert!(!self.thresholds.is_empty());
        assert!(self.thresholds.windows(2).all(|w| w[0] < w[1]));
        let n_threshold = self.thresholds.len();
        let n_covariate = self.cdfs.len() / n_threshold;

        assert_eq!(self.cdfs.len(), n_covariate * n_threshold);
        assert!(self.cdfs.iter().all(|v| (0.0..=1.0).contains(v)));
        assert!(
            self.cdfs
                .chunks_exact(self.thresholds.len())
                .all(|cdf| cdf.is_sorted())
        );

        assert!(self.covariate_groups.is_consistent());
        let covariate_dimension = self.covariate_groups.dimension;
        assert_eq!(self.covariates.len(), n_covariate * covariate_dimension);

        self.ordering_info.assert_consistent();
    }
}

impl Fit {
    pub fn interpolate_covariate_with_workspace(
        &self,
        covariate: <Self as IsotonicDistributionalRegressionFit>::Covariate<'_>,
        workspace: &mut PredictionWorkspace,
    ) -> Interpolation<'_> {
        let target: Vec<_> = unify_group_orders(covariate, &self.covariate_groups).collect();
        Interpolation::new(
            &target,
            &self.covariates,
            &self.ordering_info,
            self.increasing,
            (&self.cdfs, &self.global_cdf),
            workspace,
        )
    }
}

/// Describes which groups of columns (i.e., dimensions of our covariate) are interchangeable and
/// how they should be compared.
///
/// For each group, the order to be used and the included columns.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CovariateGroups {
    /// Dimensions can be grouped and an order applied to them together.
    pub groups: Vec<(PartialOrder, Vec<usize>)>,
    /// Dimensions that are not grouped
    pub ungrouped: Vec<usize>,
    /// Dimensionality of the covariate (size of all groups plus the number of ungrouped).
    pub dimension: usize,
}
impl CovariateGroups {
    #[must_use]
    pub fn empty(dimension: usize) -> Self {
        Self {
            groups: vec![],
            ungrouped: (0..dimension).collect(),
            dimension,
        }
    }
    pub fn parse(
        input: impl IntoIterator<Item = (impl AsRef<str>, impl IntoIterator<Item = usize>)>,
        covariate_dimension: usize,
    ) -> Result<Self, String> {
        let mut seen = vec![false; covariate_dimension];

        let mut groups = Vec::new();
        let mut ungrouped = Vec::new();
        for (name, members) in input {
            let partial_order = name.as_ref().parse()?;
            let mut new_members = Vec::new();
            for member in members {
                if member >= covariate_dimension {
                    return Err(format!(
                        "column index {member} above covariate dimension {covariate_dimension}"
                    ));
                }
                if seen[member] {
                    return Err(format!("column index {member} occurs at least twice"));
                }
                seen[member] = true;
                new_members.push(member);
            }
            if new_members.len() > 1 {
                groups.push((partial_order, new_members));
            } else {
                ungrouped.extend(new_members);
            }
        }

        ungrouped.extend(
            seen.into_iter()
                .enumerate()
                .filter(|(_, seen)| !seen)
                .map(|(i, _)| i),
        );

        Ok(Self {
            groups,
            ungrouped,
            dimension: covariate_dimension,
        })
    }
    #[must_use]
    pub fn is_consistent(&self) -> bool {
        let n_columns =
            self.ungrouped.len() + self.groups.iter().map(|(_, cs)| cs.len()).sum::<usize>();
        let mut seen = vec![false; n_columns];
        for &member in &self.ungrouped {
            if member >= n_columns || seen[member] {
                return false;
            }
            seen[member] = true;
        }
        for (_, members) in &self.groups {
            for &member in members {
                if member >= n_columns || seen[member] {
                    return false;
                }
                seen[member] = true;
            }
        }

        // is implied by pigeon-hole principle
        if !seen.into_iter().all(|b| b) {
            return false;
        }

        if self.dimension == 0 {
            return false;
        }

        true
    }
}

/// A family of partial orders for comparing vectors.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PartialOrder {
    /// Product order on vectors: \(x \preceq y \iff \forall i,\; x_i \le y_i\).
    ComponentWise,
    /// First-order stochastic dominance (FSD):
    /// \(F_X(t) \ge F_Y(t)\ \forall t\), or equivalently
    /// \( \mathbb{E}[\varphi(X)] \le \mathbb{E}[\varphi(Y)] \) for all increasing \( \varphi \).
    StochasticDominance,
    /// Increasing convex order (ICX):
    /// \( \mathbb{E}[\varphi(X)] \le \mathbb{E}[\varphi(Y)] \) for all increasing convex \( \varphi \);
    /// equivalently \( \mathbb{E}[(X - t)_+] \le \mathbb{E}[(Y - t)_+] \ \forall t \).
    IncreasingConvexOrder,
}
impl Display for PartialOrder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let str = match self {
            PartialOrder::ComponentWise => "comp",
            PartialOrder::StochasticDominance => "sd",
            PartialOrder::IncreasingConvexOrder => "icx",
        };
        f.write_str(str)
    }
}
impl FromStr for PartialOrder {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "comp" => Ok(Self::ComponentWise),
            "sd" => Ok(Self::StochasticDominance),
            "icx" => Ok(Self::IncreasingConvexOrder),
            unknown => Err(format!(
                "unknown PartialOrder identifier \"{unknown}\", I know \"comp\", \"sd\" and \"icx\"",
            )),
        }
    }
}

/// Pre-processing results that the algorithm needs
pub struct AlgorithmContext {
    /// Unique covariates in a flattened covariate-major matrix
    pub x: Vec<f64>,
    /// Total weight of all observations at the covariate, same length as `covariates` field
    pub x_weight: Vec<f64>,

    /// Responses sorted increasing, not deduplicated
    pub y: Vec<f64>,
    /// Weights in the same order as the responses
    pub weight: Vec<f64>,
    /// Covariate index of this observation
    pub x_indices: Vec<usize>,

    /// Unique response values (must have an uncensored observation)
    pub thresholds: Vec<f64>,
}

impl AlgorithmContext {
    #[must_use]
    pub fn n(&self) -> usize {
        self.y.len()
    }
    #[must_use]
    pub fn n_covariate(&self) -> usize {
        self.x_weight.len()
    }
    #[must_use]
    pub fn n_threshold(&self) -> usize {
        self.thresholds.len()
    }
    #[must_use]
    pub fn dimension(&self) -> usize {
        debug_assert_eq!(self.x.len() % self.n_covariate(), 0);

        self.x.len() / self.n_covariate()
    }
    pub fn validate(&self) {
        assert_eq!(self.x.len(), self.n_covariate());
        let covariate_dimension = self.x.len() / self.y.len();
        assert!(self.x.chunks_exact(covariate_dimension).is_sorted());
        assert_eq!(self.x_weight.len(), self.n_covariate());

        assert_eq!(self.y.len(), self.n());
        assert_eq!(self.weight.len(), self.n());
        assert_eq!(self.x_indices.len(), self.n());

        assert!(self.y.is_sorted());

        let mut seen = vec![false; self.n_covariate()];
        for &index in &self.x_indices {
            assert!(index < self.n_covariate());
            seen[index] = true;
        }
        assert!(seen.into_iter().all(|b| b));
    }
}

pub struct ExtendedAlgorithmContext {
    // Covariate index of this observation
    pub x: Vec<usize>,
    // Responses sorted increasing, deduplicated
    pub y: Vec<usize>,
    // Censoring indicators in the same order as the responses, true is observed, false is
    // right-censored.
    pub y_observed: Vec<bool>,
    // Weights in the same order as the responses
    pub weight: Vec<f64>,

    // Unique covariates in a flattened covariate-major matrix
    pub x_unique: Vec<f64>,
    // Total weight of all observations at the covariate, one for each unique covariate
    pub x_weight: Vec<f64>,

    /// Only thresholds that have at least one uncensored observation
    pub thresholds: Vec<f64>,
}

impl ExtendedAlgorithmContext {
    #[must_use]
    pub fn n(&self) -> usize {
        self.y.len()
    }
    #[must_use]
    pub fn n_covariate(&self) -> usize {
        self.x_weight.len()
    }
    #[must_use]
    pub fn n_threshold(&self) -> usize {
        self.thresholds.len()
    }
    #[must_use]
    pub fn covariate_dimension(&self) -> usize {
        debug_assert_eq!(self.x_unique.len() % self.x_weight.len(), 0);

        self.x_unique.len() / self.x_weight.len()
    }
    pub fn validate(&self) {
        assert_eq!(self.y.len(), self.n());
        assert_eq!(self.weight.len(), self.n());
        assert_eq!(self.x.len(), self.n());

        assert!(self.y.is_sorted());

        let mut seen = vec![false; self.n_covariate()];
        for &index in &self.x {
            assert!(index < self.n_covariate());
            seen[index] = true;
        }
        assert!(seen.into_iter().all(|b| b));
    }
}

pub struct AlgorithmOutput {
    pub cdfs: Vec<f64>,
    pub ordering_info: OrderingInfo,
    pub quality_indicators: QualityIndicators,
}

pub type NodeID = usize;

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Eq, PartialEq)]
pub struct OrderingInfo {
    /// Size of the graph (the number of covariates)
    pub n: usize,
    /// For each node, its children in the Hasse diagram
    pub smaller: Csr,
    /// For each node, its parents in the Hasse diagram
    pub larger: Csr,

    /// Minimal elements
    pub min: Vec<NodeID>,
    /// Maximal elements
    pub max: Vec<NodeID>,
}

impl OrderingInfo {
    #[must_use]
    pub fn empty() -> Self {
        Self {
            n: 0,
            smaller: Csr {
                offsets: vec![],
                nodes: vec![],
            },
            larger: Csr {
                offsets: vec![],
                nodes: vec![],
            },
            min: vec![],
            max: vec![],
        }
    }
    #[must_use]
    pub fn from_edges(mut edges: Vec<(NodeID, NodeID)>, n: usize) -> Self {
        edges.sort_unstable_by_key(|&(small, _)| small);
        debug_assert!(edges.windows(2).all(|w| w[0] != w[1]));
        let larger = Csr::from_sorted(&edges, |(small, large)| (small, large), n);
        edges.sort_unstable_by_key(|&(_, large)| large);
        let smaller = Csr::from_sorted(&edges, |(small, large)| (large, small), n);

        let min = (0..n).filter(|&i| smaller.get(i).is_empty()).collect();
        let max = (0..n).filter(|&i| larger.get(i).is_empty()).collect();

        Self {
            n,
            smaller,
            larger,
            min,
            max,
        }
    }
    pub fn assert_consistent(&self) {
        assert!(self.min.iter().all(|&i| i < self.n));
        assert!(self.max.iter().all(|&i| i < self.n));

        assert_eq!(
            self.min.iter().collect::<HashSet<_>>().len(),
            self.min.len(),
        );
        assert_eq!(
            self.max.iter().collect::<HashSet<_>>().len(),
            self.max.len(),
        );
    }
}

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Eq, PartialEq)]
pub struct Csr {
    // Nodes for node `i` are in `offsets[i]..offsets[i + 1]` (could be an empty range)
    offsets: Vec<usize>,
    // IDs of other nodes
    nodes: Vec<NodeID>,
}

impl Csr {
    pub fn from_sorted<F: Fn((NodeID, NodeID)) -> (NodeID, NodeID)>(
        edges: &[(NodeID, NodeID)],
        indexer: F,
        n: usize,
    ) -> Self {
        debug_assert!(edges.iter().map(|&pair| indexer(pair).0).is_sorted());

        let mut offsets = Vec::with_capacity(n);
        let mut nodes = Vec::with_capacity(edges.len());

        let mut edges_index = 0;
        for small in 0..n {
            offsets.push(nodes.len());
            while edges_index < edges.len() && indexer(edges[edges_index]).0 == small {
                nodes.push(indexer(edges[edges_index]).1);
                edges_index += 1;
            }
        }

        Csr { offsets, nodes }
    }
    #[must_use]
    pub fn get(&self, index: usize) -> &[NodeID] {
        debug_assert!(index < self.offsets.len());

        if index < self.offsets.len() - 1 {
            &self.nodes[self.offsets[index]..self.offsets[index + 1]]
        } else {
            debug_assert_eq!(index, self.offsets.len() - 1);
            &self.nodes[self.offsets[self.offsets.len() - 1]..self.nodes.len()]
        }
    }
}

#[derive(Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct QualityIndicators {
    pub precision: f64,
    pub convergence_fraction: f64,
}
