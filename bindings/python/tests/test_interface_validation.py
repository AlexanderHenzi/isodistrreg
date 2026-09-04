"""Regression tests for input validation at the Python binding boundary.

Each case used to raise a Rust panic (pyo3_runtime.PanicException) or silently
return wrong results; they must all be clean ValueErrors / correct values now.
"""

import numpy as np
import pytest

from isodistrreg import IDR, fit_park, isotonic_regression, kaplan_meier


@pytest.fixture
def fit():
    return IDR([10.0, 20.0, 15.0, 30.0], [1.0, 2.0, 3.0, 4.0])


class TestPredictionNanValidation:
    def test_predict_nan_covariate(self, fit):
        with pytest.raises(ValueError, match="NaN"):
            fit.predict([np.nan])

    def test_cdf_nan_covariate(self, fit):
        with pytest.raises(ValueError, match="NaN"):
            fit.cdf([np.nan])

    def test_cdf_at_nan_response(self, fit):
        with pytest.raises(ValueError, match="NaN"):
            fit.cdf_at([2.0], [np.nan])

    def test_quantile_nan_covariate(self, fit):
        with pytest.raises(ValueError, match="NaN"):
            fit.quantile([np.nan], [0.5])

    def test_multivariate_nan_covariate(self):
        fit2 = IDR([1.0, 2.0, 3.0], np.array([[1.0, 0.0], [2.0, 1.0], [3.0, 2.0]]))
        with pytest.raises(ValueError, match="NaN"):
            fit2.cdf(np.array([[np.nan, 1.0]]))


class TestCdfGrid:
    def test_empty_covariates(self, fit):
        assert fit.cdf_grid(np.array([]), np.array([15.0])).shape == (0, 1)

    def test_nan_rejected(self, fit):
        with pytest.raises(ValueError, match="NaN"):
            fit.cdf_grid([np.nan], [15.0])

    def test_unsorted_rejected(self, fit):
        with pytest.raises(ValueError, match="sorted"):
            fit.cdf_grid([1.0, 3.0], [30.0, 10.0, 20.0])
        with pytest.raises(ValueError, match="sorted"):
            fit.cdf_grid([3.0, 1.0], [10.0, 20.0])

    def test_matches_cdf_at(self, fit):
        cov = np.array([1.0, 2.5, 4.0])
        thr = np.array([10.0, 20.0, 30.0])
        grid = fit.cdf_grid(cov, thr)
        reference = fit.cdf_at(cov[:, None], thr)
        np.testing.assert_allclose(grid, reference)

    def test_unsqueezed_univariate_fit(self):
        fit2 = IDR([1.0, 2.0, 3.0], np.array([[1.0], [2.0], [3.0]]))
        assert fit2.cdf_grid([1.0, 2.0], [1.0, 2.0]).shape == (2, 2)


class TestAllCensoredFit:
    def test_cdf_at_and_grid_are_zero(self):
        fit = IDR([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0], y_observed=[False] * 4)
        np.testing.assert_array_equal(fit.cdf_at([2.0], [15.0]), [0.0])
        np.testing.assert_array_equal(
            fit.cdf_grid([1.0, 2.0], [1.5, 2.5]), np.zeros((2, 2))
        )


class TestSettingsValidation:
    @pytest.mark.parametrize(
        "settings",
        [{"max_iter": 0}, {"eps_abs": -1.0}, {"eps_rel": -1.0}, {"eps_abs": np.nan}],
    )
    def test_invalid_solver_settings(self, settings):
        X2 = np.array([[1.0, 0.0], [2.0, 1.0], [3.0, 2.0]])
        with pytest.raises(ValueError):
            IDR([1.0, 2.0, 3.0], X2, settings={"solver_settings": settings})


class TestFromCdfsValidation:
    def test_nan_rejected(self):
        with pytest.raises(ValueError, match="finite"):
            IDR.from_cdfs(
                np.array([[np.nan, 1.0]], dtype=np.float32), [1.0], [1.0, 2.0]
            )

    def test_out_of_range_rejected(self):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            IDR.from_cdfs(np.array([[1.5, -0.5]], dtype=np.float32), [1.0], [1.0, 2.0])

    def test_non_monotone_rejected(self):
        with pytest.raises(ValueError, match="non-decreasing"):
            IDR.from_cdfs(np.array([[0.9, 0.1]], dtype=np.float32), [1.0], [1.0, 2.0])

    def test_valid_round_trip(self):
        fit = IDR.from_cdfs(np.array([[0.5, 1.0]], dtype=np.float32), [1.0], [1.0, 2.0])
        np.testing.assert_array_equal(fit.quantile([1.0], [0.5]), [1.0])


class TestIsotonicRegressionValidation:
    def test_nan_response(self):
        with pytest.raises(ValueError, match="NaN"):
            isotonic_regression([3.0, np.nan, 2.0])

    def test_nan_covariate(self):
        with pytest.raises(ValueError, match="NaN"):
            isotonic_regression([3.0, 1.0, 2.0], X=[1.0, np.nan, 3.0])

    @pytest.mark.parametrize(
        "weights", [[1.0, -1.0, 1.0], [1.0, np.nan, 1.0], [1.0, np.inf, 1.0]]
    )
    def test_bad_weights(self, weights):
        with pytest.raises(ValueError, match="finite non-negative"):
            isotonic_regression([3.0, 1.0, 2.0], sample_weight=weights)

    def test_all_zero_weights(self):
        with pytest.raises(ValueError, match="positive"):
            isotonic_regression([3.0, 1.0, 2.0], sample_weight=[0.0, 0.0, 0.0])

    def test_constraint_out_of_bounds(self):
        with pytest.raises(ValueError, match="out of bounds"):
            isotonic_regression([1.0, 2.0], constraints=[(5, 0)])

    def test_constraint_self_loop(self):
        with pytest.raises(ValueError, match="self-loop"):
            isotonic_regression([1.0, 2.0], constraints=[(0, 0)])

    def test_constraint_cycle(self):
        with pytest.raises(ValueError, match="cycle"):
            isotonic_regression([1.0, 2.0], constraints=[(0, 1), (1, 0)])

    def test_valid_constraint_still_works(self):
        np.testing.assert_array_equal(
            isotonic_regression([1.0, 2.0], constraints=[(1, 0)]), [1.5, 1.5]
        )


class TestKaplanMeierValidation:
    def test_nan_times(self):
        with pytest.raises(ValueError, match="NaN"):
            kaplan_meier([1.0, np.nan, 3.0], [1, 1, 1])

    def test_valid_inputs_still_work(self):
        times, survival = kaplan_meier([1.0, 2.0], [True, False])
        np.testing.assert_array_equal(times, [1.0])
        np.testing.assert_array_equal(survival, [0.5])


# Every entry point taking an event indicator shares one validator: a forced cast
# to bool would read any nonzero value as an event, so a 1 = censored / 2 = event
# status coding would fit as if nothing were censored.
INDICATOR_ENTRY_POINTS = {
    "IDR": lambda t, d: IDR(t, np.arange(len(t), dtype=float), y_observed=d),
    "kaplan_meier": lambda t, d: kaplan_meier(t, d),
    "fit_park": lambda t, d: fit_park(np.arange(len(t), dtype=float), t, d),
}


@pytest.mark.parametrize(
    "call", INDICATOR_ENTRY_POINTS.values(), ids=INDICATOR_ENTRY_POINTS
)
class TestEventIndicatorValidation:
    times = [1.0, 2.0, 3.0, 4.0]

    @pytest.mark.parametrize(
        "indicator",
        [
            [1, 0, 2, 1],
            [1.0, 0.0, np.nan, 1.0],
            np.array([1, 0, -1, 1], dtype=np.int8),
            np.array([1.0, 0.0, 0.5, 1.0], dtype=np.float32),
        ],
        ids=["two", "nan", "negative", "half"],
    )
    def test_non_binary_values_are_rejected(self, call, indicator):
        with pytest.raises(ValueError, match="0/1 or boolean"):
            call(self.times, indicator)

    @pytest.mark.parametrize(
        "indicator", [["yes", "no", "yes", "yes"], [1 + 0j, 0j, 1 + 0j, 1 + 0j]]
    )
    def test_non_numeric_dtypes_are_rejected(self, call, indicator):
        with pytest.raises(ValueError, match="boolean or 0/1"):
            call(self.times, indicator)

    @pytest.mark.parametrize(
        "indicator",
        [
            [True, False, True, True],
            [1, 0, 1, 1],
            [1.0, 0.0, 1.0, 1.0],
            np.array([1, 0, 1, 1], dtype=np.uint8),
            np.array([1, 0, 1, 1], dtype=np.float32),
        ],
        ids=["bool", "int", "float", "uint8", "float32"],
    )
    def test_binary_values_are_accepted(self, call, indicator):
        call(self.times, indicator)


TEMPORAL_ENTRY_POINTS = {
    "IDR y": lambda v: IDR(v, np.arange(len(v), dtype=float)),
    "IDR X": lambda v: IDR(np.arange(len(v), dtype=float), v),
    "IDR sample_weight": lambda v: IDR(
        np.arange(len(v), dtype=float), np.arange(len(v), dtype=float), sample_weight=v
    ),
    "cdf_at": lambda v: IDR([1.0, 2.0], [1.0, 2.0]).cdf_at([1.0] * len(v), v),
    "cdf_grid": lambda v: IDR([1.0, 2.0], [1.0, 2.0]).cdf_grid([1.0], v),
    "from_cdfs": lambda v: IDR.from_cdfs(np.eye(len(v), dtype=np.float32), v, v),
    "kaplan_meier": lambda v: kaplan_meier(v, [True] * len(v)),
    "fit_park": lambda v: fit_park(np.arange(len(v), dtype=float), v, [True] * len(v)),
    "isotonic_regression": lambda v: isotonic_regression(v),
}


@pytest.mark.parametrize(
    "call", TEMPORAL_ENTRY_POINTS.values(), ids=TEMPORAL_ENTRY_POINTS
)
@pytest.mark.parametrize(
    "values",
    [
        np.array([1, 2, 3], dtype="timedelta64[D]"),
        np.array([1, 2, 3], dtype="timedelta64[ns]"),
        np.array(["2020-01-01", "2020-01-02", "2020-01-03"], dtype="datetime64[D]"),
    ],
    ids=["timedelta-days", "timedelta-ns", "datetime"],
)
def test_temporal_dtypes_are_rejected_with_a_hint(call, values):
    # A forced cast would silently turn these into raw unit counts: [ns] values
    # near 1.6e18, and [h] values 24 times the [D] values of the same durations.
    with pytest.raises(ValueError, match=r"np\.timedelta64\(1"):
        call(values)


def test_temporal_pandas_series_is_rejected_before_conversion():
    pd = pytest.importorskip("pandas")
    durations = pd.Series(np.array([1, 2, 3], dtype="timedelta64[D]"))
    with pytest.raises(ValueError, match="timedelta64"):
        IDR(durations, [1.0, 2.0, 3.0])


def test_structured_y_with_temporal_time_field_is_rejected_with_a_hint():
    y = np.array([(1, True), (2, False)], dtype=[("time", "m8[D]"), ("event", "?")])
    with pytest.raises(ValueError, match="timedelta64 times are not supported"):
        IDR(y, [1.0, 2.0])


class TestSubsampleSizeDefaults:
    def _data(self):
        rng = np.random.default_rng(42)
        x = rng.uniform(size=50)
        return x + rng.normal(scale=0.3, size=50), x

    def test_subagging_defaults_to_half_the_data(self):
        y, x = self._data()
        default = IDR(y, x, subsamples=5, seed=1)
        explicit = IDR(y, x, subsamples=5, subsample_size=0.5, seed=1)
        np.testing.assert_array_equal(default.cdf(x), explicit.cdf(x))
        # Half-sample subagging must differ from the plain fit (the old default
        # silently reproduced it).
        plain = IDR(y, x)
        assert not np.array_equal(default.cdf(x), plain.cdf(x))

    def test_bootstrap_defaults_to_full_size(self):
        y, x = self._data()
        default = IDR(y, x, subsamples=5, replace=True, seed=1)
        explicit = IDR(y, x, subsamples=5, subsample_size=1.0, replace=True, seed=1)
        np.testing.assert_array_equal(default.cdf(x), explicit.cdf(x))

    def test_full_size_without_replacement_rejected(self):
        y, x = self._data()
        with pytest.raises(ValueError, match="plain fit"):
            IDR(y, x, subsamples=5, subsample_size=1.0)
        with pytest.raises(ValueError, match="plain fit"):
            IDR(y, x, subsamples=5, subsample_size=50)
        # A single full-size subsample is the plain fit and stays allowed.
        single = IDR(y, x, subsamples=1, subsample_size=1.0)
        plain = IDR(y, x)
        np.testing.assert_array_equal(single.cdf(x), plain.cdf(x))
