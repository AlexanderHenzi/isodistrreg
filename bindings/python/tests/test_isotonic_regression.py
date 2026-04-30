import numpy as np
import pytest

from isodistrreg import isotonic_regression


@pytest.mark.parametrize(
    "X, y, w, expected",
    [
        ([1, 2, 3], [4, 3, 2], None, [3, 3, 3]),
        ([1, 2, 3], [4, 3, 2], [2] * 3, [3, 3, 3]),
    ],
)
def test_isotonic_regression(X, y, w, expected):
    result = isotonic_regression(y, X, w)
    assert result.shape == np.array(expected).shape
    assert np.allclose(result, expected)


def test_broadcasting():
    # 1-dim index covariate
    assert isotonic_regression(np.arange(10), sample_weight=1.0).shape == (10,)
    assert isotonic_regression(np.arange(10), sample_weight=np.arange(10)).shape == (
        10,
    )
    # waiting for https://github.com/PyO3/rust-numpy/pull/496 assert isotonic_regression(np.arange(10), sample_weight=np.arange(10).reshape(-1, 1)).shape == (10, 1)
    assert isotonic_regression(
        np.arange(10), sample_weight=np.arange(10).reshape(1, -1)
    ).shape == (1, 10)
    assert isotonic_regression(
        np.arange(12).reshape(3, 4), sample_weight=np.arange(12).reshape(1, 1, 3, 4)
    ).shape == (1, 1, 3, 4)

    # 1-dim explicit unique covariate
    assert isotonic_regression(
        np.arange(10), np.arange(10), sample_weight=1.0
    ).shape == (10,)
    assert isotonic_regression(
        np.arange(10), np.arange(10), sample_weight=np.arange(10)
    ).shape == (10,)
    # waiting for https://github.com/PyO3/rust-numpy/pull/496 assert isotonic_regression(np.arange(10), np.arange(10), sample_weight=np.arange(10).reshape(-1, 1)).shape == (10, 1)
    assert isotonic_regression(
        np.arange(10), np.arange(10), sample_weight=np.arange(10).reshape(1, -1)
    ).shape == (1, 10)
    assert isotonic_regression(
        np.arange(12).reshape(3, 4),
        np.arange(12).reshape(3, 4),
        sample_weight=np.arange(12).reshape(1, 1, 3, 4),
    ).shape == (1, 1, 3, 4)

    # 1-dim explicit non-unique covariate
    assert isotonic_regression(np.arange(10), np.ones(10), sample_weight=1.0).shape == (
        1,
    )
    assert isotonic_regression(
        np.arange(10), np.ones(10), sample_weight=np.arange(10)
    ).shape == (1,)
    # waiting for https://github.com/PyO3/rust-numpy/pull/496 assert isotonic_regression(np.arange(10), np.ones(10), sample_weight=np.arange(10).reshape(-1, 1)).shape == (10, 1)
    assert isotonic_regression(
        np.arange(10), np.ones(10), sample_weight=np.arange(10).reshape(1, -1)
    ).shape == (1, 1)
    assert isotonic_regression(
        np.arange(12).reshape(3, 4),
        np.ones(12).reshape(3, 4),
        sample_weight=np.arange(12).reshape(1, 1, 3, 4),
    ).shape == (1, 1, 3, 1)

    with pytest.raises(Exception):
        isotonic_regression(
            np.arange(12).reshape(3, 4),
            np.stack([np.arange(4), np.ones(4), np.arange(4)]),
            sample_weight=np.arange(12).reshape(1, 1, 3, 4),
        )

    # multidim explicit unique covariate
    covariate = np.tile(np.arange(10)[:, np.newaxis], (1, 5))
    assert isotonic_regression(np.arange(10), covariate, sample_weight=1.0).shape == (
        10,
    )
    assert isotonic_regression(
        np.arange(10), covariate, sample_weight=np.arange(10) + 1
    ).shape == (10,)
    # waiting for https://github.com/PyO3/rust-numpy/pull/496 assert isotonic_regression(np.arange(10), covariate, sample_weight=np.arange(10).reshape(-1, 1)).shape == (10, 1)
    assert isotonic_regression(
        np.arange(10), covariate, sample_weight=np.arange(10).reshape(1, -1) + 1
    ).shape == (1, 10)

    # multidim explicit non-unique covariate
    covariate = np.ones((10, 2))
    assert isotonic_regression(np.arange(10), covariate, sample_weight=1.0).shape == (
        1,
    )
    assert isotonic_regression(
        np.arange(10), covariate, sample_weight=np.arange(10) + 1
    ).shape == (1,)
    # waiting for https://github.com/PyO3/rust-numpy/pull/496 assert isotonic_regression(np.arange(10), covariate, sample_weight=np.arange(10).reshape(-1, 1) + 1).shape == (1, 1)
    assert isotonic_regression(
        np.arange(10), covariate, sample_weight=np.arange(10).reshape(1, -1) + 1
    ).shape == (1, 1)

    with pytest.raises(Exception):
        isotonic_regression(
            np.arange(12).reshape(3, 4),
            np.tile(np.stack([np.arange(4), np.ones(4), np.arange(4)]), 2),
            sample_weight=np.arange(12).reshape(1, 1, 3, 4),
        )

    # general partial orders
    assert np.allclose(isotonic_regression([1, 2], constraints=[(1, 0)]), [1.5, 1.5])
    assert np.allclose(
        isotonic_regression([[1, 2]], constraints=[(1, 0)]), [[1.5, 1.5]]
    )
    assert np.allclose(
        isotonic_regression([1, 2], sample_weight=[[[2.0]]], constraints=[(1, 0)]),
        [[[1.5, 1.5]]],
    )
