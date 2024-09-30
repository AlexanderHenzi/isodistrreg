use crate::error::Error;

pub fn validate<'a>(
    covariates: impl Iterator<Item = &'a [f64]>,
    responses: &[f64],
    censoring: Option<&[bool]>,
    weights: Option<&[f64]>,
) -> Result<usize, Error> {
    let n = {
        let mut i = 0;
        for vs in covariates {
            if vs.iter().any(|v| !v.is_finite()) {
                return Err(Error::NonFiniteFloats);
            }
            i += 1;
        }
        i
    };

    let censoring_len = censoring.map(|slice| slice.len());
    let weights_len = weights.map(|slice| slice.len());
    let lengths_equal = responses.len() == n
        && weights_len.is_none_or(|l| l == n)
        && censoring_len.is_none_or(|l| l == n);
    if !lengths_equal {
        Err(Error::IncompatibleShapes {
            covariate_len: n,
            response_len: responses.len(),
            weight_len: weights_len,
            y_observed_len: censoring_len,
        })
    } else if n == 0 {
        Err(Error::Empty)
    } else if responses.iter().any(|f| !f.is_finite())
        | weights.is_some_and(|slice| slice.iter().any(|f| !f.is_finite()))
    {
        Err(Error::NonFiniteFloats)
    } else if weights.is_some_and(|slice| slice.iter().any(|&w| w < 0.0)) {
        Err(Error::NegativeWeights)
    } else {
        Ok(n)
    }
}
