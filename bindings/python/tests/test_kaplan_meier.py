import numpy as np
import pandas as pd
import pytest

from isodistrreg import kaplan_meier


TIME_VALUES = [1, 3, 5, 5, 7, 9, 11]
EVENT_VALUES = [True, False, True, True, True, False, True]

EXPECTED_TIMES = np.array([1, 5, 7, 11])
EXPECTED_SURVIVAL = np.array([6 / 7, (6 / 7) * (3 / 5), (6 / 7) * (3 / 5) * (2 / 3), 0.0])


@pytest.mark.parametrize(
    "dtype",
    [np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64],
)
def test_preserves_integer_dtype(dtype):
    times, survival = kaplan_meier(np.asarray(TIME_VALUES, dtype=dtype), EVENT_VALUES)
    assert times.dtype == dtype
    assert survival.dtype == np.float64
    np.testing.assert_array_equal(times, EXPECTED_TIMES.astype(dtype))
    np.testing.assert_allclose(survival, EXPECTED_SURVIVAL)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_preserves_float_dtype(dtype):
    times, survival = kaplan_meier(np.asarray(TIME_VALUES, dtype=dtype), EVENT_VALUES)
    assert times.dtype == dtype
    assert survival.dtype == np.float64
    np.testing.assert_allclose(times, EXPECTED_TIMES.astype(dtype))
    np.testing.assert_allclose(survival, EXPECTED_SURVIVAL)


def test_accepts_pandas_series_int64():
    series = pd.Series(TIME_VALUES, dtype="int64")
    times, survival = kaplan_meier(series, EVENT_VALUES)
    assert times.dtype == np.int64
    np.testing.assert_array_equal(times, EXPECTED_TIMES)
    np.testing.assert_allclose(survival, EXPECTED_SURVIVAL)


def test_accepts_pandas_dt_days_and_values_agree():
    train = pd.DataFrame(
        {
            "time": pd.to_timedelta(TIME_VALUES, unit="D").astype("timedelta64[s]"),
            "event": EVENT_VALUES,
        }
    )
    t1, s1 = kaplan_meier(train["time"].dt.days, train["event"])
    t2, s2 = kaplan_meier(train["time"].dt.days.values, train["event"])
    assert t1.dtype == np.int64
    assert t2.dtype == np.int64
    np.testing.assert_array_equal(t1, t2)
    np.testing.assert_allclose(s1, s2)


def test_python_list_default_int():
    times, _ = kaplan_meier(TIME_VALUES, EVENT_VALUES)
    assert times.dtype == np.asarray(TIME_VALUES).dtype


def test_weights_match_unit_weight_unweighted():
    weights = np.ones(len(TIME_VALUES), dtype=np.float64)
    t_unw, s_unw = kaplan_meier(TIME_VALUES, EVENT_VALUES)
    t_w, s_w = kaplan_meier(TIME_VALUES, EVENT_VALUES, weights)
    np.testing.assert_array_equal(t_unw, t_w)
    np.testing.assert_allclose(s_unw, s_w)


def test_weights_accept_int_dtype():
    weights = np.ones(len(TIME_VALUES), dtype=np.int32)
    t, s = kaplan_meier(TIME_VALUES, EVENT_VALUES, weights)
    np.testing.assert_array_equal(t, EXPECTED_TIMES)
    np.testing.assert_allclose(s, EXPECTED_SURVIVAL)


def test_events_accept_int_dtype():
    events = np.asarray(EVENT_VALUES, dtype=np.int32)
    t, s = kaplan_meier(TIME_VALUES, events)
    np.testing.assert_array_equal(t, EXPECTED_TIMES)
    np.testing.assert_allclose(s, EXPECTED_SURVIVAL)


def test_rejects_mismatched_lengths():
    with pytest.raises(ValueError, match="same length"):
        kaplan_meier([1, 2, 3], [True, False])


def test_rejects_negative_weights():
    with pytest.raises(ValueError, match="nonnegative"):
        kaplan_meier([1.0, 2.0], [True, False], [-1.0, 1.0])


def test_rejects_2d_input():
    with pytest.raises(ValueError, match="1-dimensional"):
        kaplan_meier(np.zeros((2, 2)), [True, False, True, False])


def test_rejects_complex_dtype():
    with pytest.raises(ValueError, match="unsupported dtype"):
        kaplan_meier(np.array([1 + 0j, 2 + 0j]), [True, False])
