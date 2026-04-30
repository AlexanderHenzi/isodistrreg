from sklearn.utils.estimator_checks import check_estimator

from isodistrreg import IsotonicDistributionalRegressor


def test_api_adherence():
    assert check_estimator(IsotonicDistributionalRegressor())
