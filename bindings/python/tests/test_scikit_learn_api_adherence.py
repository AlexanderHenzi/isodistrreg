import pytest

pytest.importorskip("sklearn")

from sklearn.utils.estimator_checks import parametrize_with_checks

from isodistrreg import IsotonicDistributionalRegressor


def expected_failed_checks(estimator):
    if estimator.subsamples is None:
        return {}
    # The same reason scikit-learn's own bootstrap ensembles give: random
    # subsamples of the weighted data are not the same draws as random
    # subsamples of the data with rows repeated.
    reason = "sample_weight is not equivalent to repeating rows when subsampling"
    return {
        "check_sample_weight_equivalence_on_dense_data": reason,
        "check_sample_weight_equivalence_on_sparse_data": reason,
    }


@parametrize_with_checks(
    [
        IsotonicDistributionalRegressor(),
        IsotonicDistributionalRegressor(decreasing=True, response_order="hazard"),
        IsotonicDistributionalRegressor(subsamples=3, random_state=0),
    ],
    expected_failed_checks=expected_failed_checks,
)
def test_api_adherence(estimator, check):
    check(estimator)
