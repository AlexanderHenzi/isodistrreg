use crate::structures::sealed::Sealed;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::str::FromStr;

#[derive(Copy, Clone, Eq, PartialEq)]
pub enum StochasticOrder {
    StochasticDominance,
    HazardRateOrder,
}
impl FromStr for StochasticOrder {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "sd" => Ok(Self::StochasticDominance),
            "hazard" => Ok(Self::HazardRateOrder),
            unknown => Err(format!(
                "unknown StochasticOrder identifier \"{unknown}\", I know \"sd\" and \"hro\"",
            )),
        }
    }
}

#[allow(private_bounds)]
pub trait Direction: Sealed + Default {
    /// The left value is not allowed to have this ordering w.r.t. the right value.
    ///
    /// E.g., if we have `FORBIDDEN_ORDERING == Greater`, then we work toward a state with
    /// `left <= right`, i.e., non-decreasing.
    const FORBIDDEN_ORDERING: Ordering;
    const IS_INCREASING: bool;
    type REVERSE: Direction;
}

/// Increasing from left to right (non-decreasing, to be precise).
#[derive(Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Increasing;
impl Direction for Increasing {
    const FORBIDDEN_ORDERING: Ordering = Ordering::Greater;
    const IS_INCREASING: bool = true;
    type REVERSE = Decreasing;
}
impl Sealed for Increasing {}
impl Default for Increasing {
    fn default() -> Self {
        Self
    }
}

/// Decreasing from left to right (non-increasing, to be precise).
#[derive(Copy, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Decreasing;
impl Direction for Decreasing {
    const FORBIDDEN_ORDERING: Ordering = Ordering::Less;
    const IS_INCREASING: bool = false;
    type REVERSE = Increasing;
}
impl Sealed for Decreasing {}
impl Default for Decreasing {
    fn default() -> Self {
        Self
    }
}

mod sealed {
    pub(super) trait Sealed {}
}

pub trait ExecutionMode {
    const PARALLEL: bool;
}
pub struct Serial;
impl ExecutionMode for Serial {
    const PARALLEL: bool = false;
}
pub struct Parallel;
impl ExecutionMode for Parallel {
    const PARALLEL: bool = cfg!(feature = "parallel");
}

/// The data for which we compute the algorithm, closely related to the data provided by the caller
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Observation<X, Y, S> {
    /// Index in `0..n_covariate` or original data value
    pub x: X,
    /// Index in `0..n_threshold` or original data value
    pub y: Y,
    /// Can indicate censoring for survival applications.
    ///
    /// If not used and `S = ()`, the outcome is implicitly always observed. If `S = bool`, then
    /// `true` indicates observed and `false` indicates right-censored.
    pub observed: S,
    /// Non-negative and finite weight
    pub weight: f64,
}

impl<T: Into<f64>> From<T> for Observation<(), f64, ()> {
    fn from(y: T) -> Self {
        Observation {
            x: (),
            y: y.into(),
            observed: (),
            weight: 1.0,
        }
    }
}

impl<T: Into<f64>, U: Into<f64>> From<(T, U)> for Observation<f64, f64, ()> {
    fn from((x, y): (T, U)) -> Self {
        Observation {
            x: x.into(),
            y: y.into(),
            observed: (),
            weight: 1.0,
        }
    }
}

impl<T: Into<f64>, U: Into<f64>, V: Into<f64>> From<(T, U, V)> for Observation<f64, f64, ()> {
    fn from((x, y, weight): (T, U, V)) -> Self {
        Observation {
            x: x.into(),
            y: y.into(),
            observed: (),
            weight: weight.into(),
        }
    }
}

impl<T: Into<f64>, U: Into<f64>, V: Into<bool>, W: Into<f64>> From<(T, U, V, W)>
    for Observation<f64, f64, bool>
{
    fn from((x, y, censored, weight): (T, U, V, W)) -> Self {
        Observation {
            x: x.into(),
            y: y.into(),
            observed: censored.into(),
            weight: weight.into(),
        }
    }
}

impl<T: Into<f64>, U: Into<f64>, V: Into<bool>> From<(T, U, V)> for Observation<f64, f64, bool> {
    fn from((x, y, observed): (T, U, V)) -> Self {
        Observation {
            x: x.into(),
            y: y.into(),
            observed: observed.into(),
            weight: 1.0,
        }
    }
}
