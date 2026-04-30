import numpy as np
import pytest

from isodistrreg import IDR


@pytest.mark.parametrize(
    "X, y, y_observed, expected",
    [
        ([[1], [2], [3]], [4, 3, 2], [False, True, True], [[1 / 3, 2 / 3]] * 3),
        ([[1], [2], [3]], [4, 3, 2], [False] * 3, [] * 3),
    ],
)
def test_isotonic_distributional_regression(X, y, y_observed, expected):
    fit = IDR(y, X, y_observed)
    result = fit.cdf(X)
    assert result.shape == np.array(expected).reshape(len(X), -1).shape
    assert np.allclose(result, expected)
    assert np.allclose(fit.X, np.unique(np.array(X).flat).reshape(-1, 1))
    assert np.allclose(fit.thresholds, np.unique(np.array(y)[np.array(y_observed)]))


def test_outside_range():
    fit = IDR(np.arange(10, 30)[::-1], np.arange(20), [True] * 19 + [False])
    result = fit.cdf_at([5.0, 15.0], [0.0, 40.0])
    assert np.allclose(result, [[0.0, 1.0]] * 2)

    fit = IDR(
        np.arange(10, 30)[::-1], np.arange(20), [True] * 19 + [False], decreasing=True
    )
    result = fit.cdf_at([5.0, 15.0], [0.0, 40.0])
    assert np.allclose(result, [[0.0, 1.0]] * 2)


def test_outside_range_subagging():
    fit = IDR(
        np.arange(10, 30)[::-1],
        np.arange(20),
        [True] * 19 + [False],
        subsamples=5,
        subsample_size=0.5,
    )
    result = fit.cdf_at([5.0, 15.0], [0.0, 40.0])
    assert np.allclose(result, [[0.0, 1.0]] * 2)
