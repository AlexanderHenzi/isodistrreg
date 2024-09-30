use std::fmt::Display;

#[derive(Debug)]
pub enum Error {
    /// Config data could not be read.
    ConfigParseError(&'static str),
    /// Arguments are for a total order but this is a partial ordered context, or vice versa.
    CovariateDimensionMismatch {
        shape: Vec<usize>,
        message: &'static str,
    },
    /// The specified covariate order doesn't make sense (partial orders only).
    CovariateOrderInconsistency,
    /// No data provided for fit.
    Empty,
    /// Input to the fit algorithm has incompatible shapes.
    IncompatibleShapes {
        covariate_len: usize,
        response_len: usize,
        weight_len: Option<usize>,
        y_observed_len: Option<usize>,
    },
    /// Some floats were NAN or infinite.
    NonFiniteFloats,
    /// Some weights were strictly below 0.0.
    NegativeWeights,
    /// The requested functionality is not yet implemented.
    NotImplemented(&'static str),
    /// Operation doesn't make sense on this order
    OrderMismatch,
    /// Could not broadcast arguments passed to a prediction method.
    ShapeMismatch {
        covariate_shape: Vec<usize>,
        other_shape: Vec<usize>,
    },
    /// Parameters provided about subagging are inconsistent.
    SubaggingParameterInconsistency(&'static str),
}

impl Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use Error::*;
        match self {
            ConfigParseError(message) => f.write_str(message),
            CovariateDimensionMismatch { shape, message } => {
                write!(f, "covariates had shape {shape:?}, but {message}")
            }
            CovariateOrderInconsistency => write!(f, "covariate order inconsistent"),
            Empty => write!(f, "all arguments were of length 0"),
            IncompatibleShapes {
                covariate_len,
                response_len,
                weight_len: weights_len,
                y_observed_len: censoring_len,
            } => {
                write!(
                    f,
                    "incompatible shapes; covariate {covariate_len}, response {response_len}",
                )?;
                if let Some(len) = weights_len {
                    write!(f, ", weights {len}")?;
                }
                if let Some(len) = censoring_len {
                    write!(f, ", censoring {len}")?;
                }
                Ok(())
            }
            NegativeWeights => write!(f, "at least one weight value was negative"),
            NonFiniteFloats => write!(f, "some float values were infinite, nan or subnormal"),
            NotImplemented(message) => f.write_str(message),
            OrderMismatch => write!(
                f,
                "operation is not supported, you might be mixing up partial and total orders"
            ),
            ShapeMismatch {
                covariate_shape,
                other_shape,
            } => {
                write!(
                    f,
                    "arguments of shape {covariate_shape:?} (after possible removal of covariate dimension) and {other_shape:?} can't be broadcast together",
                )
            }
            SubaggingParameterInconsistency(message) => f.write_str(message),
        }
    }
}
