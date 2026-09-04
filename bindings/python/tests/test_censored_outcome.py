"""Right-censored outcomes supplied as a structured array, or as its field views.

A structured array's field views (``y["time"]``) stride by the record size, so
their elements are not aligned to their own dtype. The binding used to refuse
such input outright; it now copies it into an aligned buffer. The structured
array itself is also accepted in place of the ``(y, y_observed)`` pair.
"""

import numpy as np
import pytest

from isodistrreg import IDR, fit_park, isotonic_regression, kaplan_meier
from isodistrreg._core import split_censored_outcome

FIELD_ORDERS = {
    "time-first": [("time", "f8"), ("event", "?")],
    "event-first": [("event", "?"), ("time", "f8")],  # scikit-survival's layout
}


@pytest.fixture
def sample():
    rng = np.random.default_rng(0)
    n = 80
    x = rng.normal(size=n)
    time = np.exp(x + rng.normal(size=n))
    event = rng.random(n) < 0.7
    return x, time, event


def structured(time, event, dtype):
    y = np.empty(len(time), dtype=dtype)
    y["time"] = time
    y["event"] = event
    return y


def small_multivariate(sample):
    """A censored 2-D design small enough for the censored partial-order solver,
    whose cost grows exponentially with the width of the covariate poset."""
    x, time, event = sample
    n = 12
    return np.column_stack([x[:n], x[:n] ** 2]), time[:n], event[:n]


def unaligned_copy(values):
    """A C-contiguous array whose data pointer is one byte off alignment."""
    values = np.ascontiguousarray(values)
    buffer = np.zeros(values.nbytes + 1, dtype=np.uint8)
    view = np.frombuffer(buffer, dtype=values.dtype, offset=1, count=values.size)
    view = view.reshape(values.shape)
    np.copyto(view, values)
    assert not view.flags.aligned
    return view


def test_field_views_are_unaligned(sample):
    x, time, event = sample
    y = structured(time, event, FIELD_ORDERS["time-first"])
    assert not y["time"].flags.aligned
    assert y["time"].strides == (9,)


class TestIDR:
    def test_field_views_match_aligned_arrays(self, sample):
        x, time, event = sample
        y = structured(time, event, FIELD_ORDERS["time-first"])
        reference = IDR(time, x, y_observed=event)
        fit = IDR(y["time"], x, y_observed=y["event"])
        np.testing.assert_array_equal(fit.thresholds, reference.thresholds)
        np.testing.assert_array_equal(fit.cdf(x), reference.cdf(x))

    @pytest.mark.parametrize("dtype", FIELD_ORDERS.values(), ids=FIELD_ORDERS)
    def test_structured_y_matches_separate_arguments(self, sample, dtype):
        x, time, event = sample
        reference = IDR(time, x, y_observed=event)
        fit = IDR(structured(time, event, dtype), x)
        np.testing.assert_array_equal(fit.thresholds, reference.thresholds)
        np.testing.assert_array_equal(fit.cdf(x), reference.cdf(x))

    def test_structured_y_with_partial_order(self, sample):
        X, time, event = small_multivariate(sample)
        reference = IDR(time, X, y_observed=event)
        fit = IDR(structured(time, event, FIELD_ORDERS["event-first"]), X)
        np.testing.assert_array_equal(fit.cdf(X), reference.cdf(X))

    @pytest.mark.parametrize(
        "time_dtype, threshold_dtype",
        [("f4", np.float32), ("f8", np.float64), ("i8", np.float64)],
    )
    def test_time_field_dtype_sets_threshold_dtype(
        self, sample, time_dtype, threshold_dtype
    ):
        x, time, event = sample
        y = structured(time, event, [("time", time_dtype), ("event", "?")])
        assert IDR(y, x).thresholds.dtype == threshold_dtype

    def test_prediction_accepts_unaligned_input(self, sample):
        x, time, event = sample
        fit = IDR(time, x, y_observed=event)
        rec = np.empty(len(x), dtype=[("x", "f8"), ("t", "f8"), ("q", "f8")])
        rec["x"] = x
        rec["t"] = time
        rec["q"] = np.linspace(0.0, 1.0, len(x))
        np.testing.assert_array_equal(fit.predict(rec["x"]), fit.predict(x))
        np.testing.assert_array_equal(fit.cdf(rec["x"]), fit.cdf(x))
        np.testing.assert_array_equal(
            fit.cdf_at(rec["x"], rec["t"]), fit.cdf_at(x, time)
        )
        np.testing.assert_array_equal(
            fit.quantile(rec["x"], rec["q"]), fit.quantile(x, rec["q"].copy())
        )
        np.testing.assert_array_equal(
            fit.cdf_grid(np.sort(rec["x"]), np.sort(rec["t"])),
            fit.cdf_grid(np.sort(x), np.sort(time)),
        )

    def test_unaligned_data_pointer(self, sample):
        # Contiguous strides but an off-by-one data pointer: the other way an
        # array can fail alignment.
        x, time, event = sample
        reference = IDR(time, x, y_observed=event)
        fit = IDR(unaligned_copy(time), unaligned_copy(x), y_observed=event)
        np.testing.assert_array_equal(fit.cdf(x), reference.cdf(x))
        np.testing.assert_array_equal(
            fit.cdf(unaligned_copy(x[:, None])[:, 0]), reference.cdf(x)
        )

    def test_unaligned_covariates_and_weights(self, sample):
        x, time, event = sample
        rec = np.empty(len(x), dtype=[("x", "f8"), ("w", "f8"), ("flag", "?")])
        rec["x"] = x
        rec["w"] = np.linspace(1.0, 2.0, len(x))
        reference = IDR(time, x, y_observed=event, sample_weight=rec["w"].copy())
        fit = IDR(time, rec["x"], y_observed=event, sample_weight=rec["w"])
        np.testing.assert_array_equal(fit.cdf(x), reference.cdf(x))

    def test_unaligned_two_dimensional_covariates(self, sample):
        X, time, event = small_multivariate(sample)
        reference = IDR(time, X, y_observed=event)
        fit = IDR(time, unaligned_copy(X), y_observed=event)
        np.testing.assert_array_equal(fit.cdf(X), reference.cdf(X))

    def test_from_cdfs_accepts_unaligned_input(self, sample):
        x, time, event = sample
        fit = IDR(time, x, y_observed=event)
        rebuilt = IDR.from_cdfs(
            unaligned_copy(fit.cdf(fit.X)),
            unaligned_copy(fit.X),
            unaligned_copy(fit.thresholds),
        )
        np.testing.assert_array_equal(rebuilt.cdf(x), fit.cdf(x))

    def test_structured_y_rejects_y_observed(self, sample):
        x, time, event = sample
        y = structured(time, event, FIELD_ORDERS["time-first"])
        with pytest.raises(ValueError, match="y_observed must be omitted"):
            IDR(y, x, y_observed=event)

    @pytest.mark.parametrize(
        "dtype, message",
        [
            ([("time", "f8"), ("event", "i1")], "both.*numeric|store the event"),
            ([("a", "?"), ("b", "?")], "a: bool and b: bool"),
            ([("label", "U3"), ("event", "?")], "label: <U3"),
            ([("time", "f8")], "found 1 field"),
            ([("time", "f8"), ("event", "?"), ("id", "i8")], "found 3 field"),
        ],
    )
    def test_unusable_record_layouts_are_rejected(self, dtype, message):
        y = np.zeros(4, dtype=dtype)
        with pytest.raises(ValueError, match=message):
            IDR(y, [1.0, 2.0, 3.0, 4.0])


class TestKaplanMeier:
    def test_field_views_and_structured_y(self, sample):
        _, time, event = sample
        y = structured(time, event, FIELD_ORDERS["event-first"])
        times, survival = kaplan_meier(time, event)
        for args in [(y["time"], y["event"]), (y,)]:
            got_times, got_survival = kaplan_meier(*args)
            np.testing.assert_array_equal(got_times, times)
            np.testing.assert_array_equal(got_survival, survival)

    def test_structured_y_with_weights(self, sample):
        _, time, event = sample
        y = structured(time, event, FIELD_ORDERS["time-first"])
        weight = np.linspace(1.0, 2.0, len(time))
        times, survival = kaplan_meier(time, event, weight)
        got_times, got_survival = kaplan_meier(y, weight=weight)
        np.testing.assert_array_equal(got_times, times)
        np.testing.assert_array_equal(got_survival, survival)

    def test_integer_time_field_keeps_its_dtype(self, sample):
        _, time, event = sample
        y = structured(np.round(time * 10), event, [("time", "i8"), ("event", "?")])
        times, _ = kaplan_meier(y)
        assert times.dtype == np.int64

    def test_y_observed_required_for_plain_y(self, sample):
        _, time, _ = sample
        with pytest.raises(ValueError, match="y_observed is required"):
            kaplan_meier(time)

    def test_structured_y_rejects_y_observed(self, sample):
        _, time, event = sample
        y = structured(time, event, FIELD_ORDERS["time-first"])
        with pytest.raises(ValueError, match="y_observed must be omitted"):
            kaplan_meier(y, event)


class TestFitPark:
    def test_structured_y_and_field_views(self, sample):
        x, time, event = sample
        y = structured(time, event, FIELD_ORDERS["time-first"])
        centers = np.linspace(-2.0, 2.0, 10)
        reference = fit_park(x, time, event, centers=centers)
        for args in [(x, y["time"], y["event"]), (x, y)]:
            for got, expected in zip(fit_park(*args, centers=centers), reference):
                np.testing.assert_array_equal(got, expected)

    def test_y_observed_required_for_plain_y(self, sample):
        x, time, _ = sample
        with pytest.raises(ValueError, match="y_observed is required"):
            fit_park(x, time)


def test_isotonic_regression_accepts_unaligned_input(sample):
    x, time, _ = sample
    rec = np.empty(len(x), dtype=[("x", "f8"), ("y", "f8"), ("w", "f8")])
    rec["x"] = x
    rec["y"] = time
    rec["w"] = np.linspace(1.0, 2.0, len(x))
    np.testing.assert_array_equal(
        isotonic_regression(rec["y"], rec["x"], sample_weight=rec["w"]),
        isotonic_regression(time, x, sample_weight=rec["w"].copy()),
    )


def test_split_censored_outcome_passes_plain_inputs_through(sample):
    _, time, event = sample
    assert split_censored_outcome(time) == (time, None)
    assert split_censored_outcome(time, event) == (time, event)
    as_list = [1.0, 2.0]
    assert split_censored_outcome(as_list, None) == (as_list, None)


def test_split_censored_outcome_returns_field_views(sample):
    _, time, event = sample
    y = structured(time, event, FIELD_ORDERS["event-first"])
    got_time, got_event = split_censored_outcome(y)
    np.testing.assert_array_equal(got_time, time)
    np.testing.assert_array_equal(got_event, event)
    assert got_event.dtype == np.bool_


class TestEstimator:
    sklearn = pytest.importorskip("sklearn")

    @pytest.mark.parametrize("dtype", FIELD_ORDERS.values(), ids=FIELD_ORDERS)
    def test_structured_y_matches_separate_arguments(self, sample, dtype):
        from isodistrreg import IsotonicDistributionalRegressor

        x, time, event = sample
        X = x[:, None]
        reference = IsotonicDistributionalRegressor().fit(X, time, y_observed=event)
        fit = IsotonicDistributionalRegressor().fit(X, structured(time, event, dtype))
        np.testing.assert_array_equal(fit.thresholds_, reference.thresholds_)
        np.testing.assert_array_equal(fit.cdf(X), reference.cdf(X))

    def test_field_views(self, sample):
        from isodistrreg import IsotonicDistributionalRegressor

        x, time, event = sample
        X = x[:, None]
        y = structured(time, event, FIELD_ORDERS["time-first"])
        reference = IsotonicDistributionalRegressor().fit(X, time, y_observed=event)
        fit = IsotonicDistributionalRegressor().fit(X, y["time"], y_observed=y["event"])
        np.testing.assert_array_equal(fit.cdf(X), reference.cdf(X))

    def test_structured_y_rejects_y_observed(self, sample):
        from isodistrreg import IsotonicDistributionalRegressor

        x, time, event = sample
        y = structured(time, event, FIELD_ORDERS["time-first"])
        with pytest.raises(ValueError, match="y_observed must be omitted"):
            IsotonicDistributionalRegressor().fit(x[:, None], y, y_observed=event)
