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


def test_quantile_levels_outside_range_raise():
    fit = IDR([1.0, 2.0, 3.0], [[1], [2], [3]])
    # Valid levels (boundaries included) pass.
    fit.quantile([[2]], [0.0, 0.5, 1.0])
    for bad in (-0.1, 1.5, np.nan):
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            fit.quantile([[2]], [0.5, bad])


def test_quantile_levels_accept_integer_arrays():
    fit = IDR([1.0, 2.0, 3.0], [[1], [2], [3]])
    # Integer levels used to be refused with a bare type error; they are now
    # coerced like every other array argument.
    np.testing.assert_array_equal(
        fit.quantile([[2]], np.array([0, 1])), fit.quantile([[2]], [0.0, 1.0])
    )


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


def test_quantile_monotonicity():
    rng = np.random.default_rng(11)
    n = 120
    x = rng.uniform(0, 10, n)
    y = rng.gamma(shape=2.0, scale=(x + 1) / 4)
    raw = IDR(y, x)

    xs = np.linspace(0.2, 9.8, 300)
    q90 = np.array([[0.9]])
    result = raw.quantile(xs[:, None], q90)[:, 0]

    assert np.all(np.diff(result) >= 0.0)


# `.X` and `.thresholds` round-trip the user's input dtype:
#   - np.float64 in  → np.float64 out
#   - np.float32 in  → np.float32 out
#   - everything else (ints, lists, Python scalars) → np.float64 out
#
# Internal storage is f64 so close-but-distinct f64 inputs don't collide on a f32 cast;
# for f32 inputs the f64 internal values originated from an `as f64` upcast so the f32
# we expose is bit-exact in the original f32 values.


def test_input_dtype_f64_preserved():
    X = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    y = np.array([4.0, 3.0, 2.0], dtype=np.float64)
    fit = IDR(y, X)
    assert fit.X.dtype == np.float64
    assert fit.thresholds.dtype == np.float64


def test_input_dtype_f32_preserved():
    X = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    y = np.array([4.0, 3.0, 2.0], dtype=np.float32)
    fit = IDR(y, X)
    assert fit.X.dtype == np.float32
    assert fit.thresholds.dtype == np.float32


def test_input_dtype_int_becomes_f64():
    # Integer inputs default to f64 (NumPy's default for arithmetic).
    X = np.array([1, 2, 3], dtype=np.int64)
    y = np.array([4, 3, 2], dtype=np.int64)
    fit = IDR(y, X)
    assert fit.X.dtype == np.float64
    assert fit.thresholds.dtype == np.float64


def test_input_dtype_pythonlist_becomes_f64():
    # Plain Python lists also default to f64.
    fit = IDR([4.0, 3.0, 2.0], [[1.0], [2.0], [3.0]])
    assert fit.X.dtype == np.float64
    assert fit.thresholds.dtype == np.float64


def test_input_dtype_mixed_is_tracked_per_array():
    # `.X` and `.thresholds` track independently.
    X = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    y = np.array([4.0, 3.0, 2.0], dtype=np.float64)
    fit = IDR(y, X)
    assert fit.X.dtype == np.float32
    assert fit.thresholds.dtype == np.float64


def test_input_dtype_f32_round_trip_is_exact():
    # An f32 → f64 → f32 round-trip preserves the bit pattern. The covariate values
    # we get back must equal the original f32 array exactly (not just approximately).
    X = np.array([0.1, 0.2, 0.3], dtype=np.float32)
    y = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    fit = IDR(y, X)
    np.testing.assert_array_equal(fit.X, np.sort(np.unique(X)))


def test_input_dtype_2d_X_preserved():
    # Multivariate covariates: the dtype propagation should work on 2-D arrays too.
    X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], dtype=np.float32)
    y = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    fit = IDR(y, X)
    assert fit.X.dtype == np.float32
    assert fit.thresholds.dtype == np.float32


def test_input_dtype_survives_pickle():
    # Round-tripping through pickle preserves the dtype tag.
    import pickle

    X = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    y = np.array([4.0, 3.0, 2.0], dtype=np.float64)
    fit = IDR(y, X)
    restored = pickle.loads(pickle.dumps(fit))
    assert restored.X.dtype == np.float32
    assert restored.thresholds.dtype == np.float64


def test_from_cdfs_preserves_input_dtype():
    cdfs = np.array([[0.25, 0.5, 1.0], [0.2, 0.4, 1.0]], dtype=np.float32)
    X_f32 = np.array([1.0, 2.0], dtype=np.float32)
    y_f64 = np.array([10.0, 20.0, 30.0], dtype=np.float64)
    fit = IDR.from_cdfs(cdfs, X_f32, y_f64)
    assert fit.X.dtype == np.float32
    assert fit.thresholds.dtype == np.float64


def test_from_cdfs_silently_narrows_f64_cdfs_to_f32():
    # CDFs are stored as f32 internally regardless of how the fit was produced.
    # Passing an f64 CDF array to `from_cdfs` is allowed but silently narrowed:
    # the documented behaviour is "no precision loss vs. a fresh fit (storage is
    # f32 either way), but precision is lost relative to the user's f64 input."
    #
    # We pick CDF values that are exactly representable in both f32 and f64 so the
    # f64→f32 narrowing is bit-exact. The same values fitted via the f32 entry
    # point must produce a bit-identical CDF — this is the "no precision loss vs.
    # a fresh fit" guarantee.
    cdfs_f64 = np.array([[0.25, 0.5, 1.0], [0.125, 0.5, 1.0]], dtype=np.float64)
    cdfs_f32 = cdfs_f64.astype(np.float32)
    X = np.array([1.0, 2.0], dtype=np.float64)
    y = np.array([10.0, 20.0, 30.0], dtype=np.float64)

    fit_from_f64 = IDR.from_cdfs(cdfs_f64, X, y)
    fit_from_f32 = IDR.from_cdfs(cdfs_f32, X, y)

    # The fit succeeded with f64 input (no exception, despite the dtype mismatch
    # against the declared `PyArrayLike2<f32, AllowTypeChange>` parameter).
    assert fit_from_f64 is not None

    # Observable evidence of f32 internal storage: predicted CDF values come out
    # as np.float32, and the f64-input fit is bit-identical to the f32-input fit.
    eval_X = np.array([1.0, 2.0], dtype=np.float64)
    cdf_from_f64 = fit_from_f64.cdf(eval_X)
    cdf_from_f32 = fit_from_f32.cdf(eval_X)
    assert cdf_from_f64.dtype == np.float32
    assert cdf_from_f32.dtype == np.float32
    np.testing.assert_array_equal(cdf_from_f64, cdf_from_f32)

    # And the CDF values are the original (f32-narrowed) ones.
    np.testing.assert_array_equal(cdf_from_f64, cdfs_f32)


# Tests for the typed-FitImpl rewrite: storage now matches the user's input dtype
# directly (no internal f64 widening), so f32 inputs stay in f32 storage. Each of the
# four (X, Y) dtype combos is exercised end-to-end.


# Parametrize over all four (X, Y) combinations. The same data fitted in each combo
# should produce equivalent CDFs (to ~ulp precision), confirming that the algorithm
# behaves consistently regardless of which combo the user picks.
_DTYPES = [np.float32, np.float64]


@pytest.mark.parametrize("x_dtype", _DTYPES)
@pytest.mark.parametrize("y_dtype", _DTYPES)
def test_all_four_combinations_fit(x_dtype, y_dtype):
    X = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=x_dtype)
    y = np.array([10.0, 30.0, 20.0, 50.0, 40.0], dtype=y_dtype)
    fit = IDR(y, X)
    assert fit.X.dtype == x_dtype
    assert fit.thresholds.dtype == y_dtype
    # CDF output is always f32, independent of storage precision.
    assert fit.cdf(X).dtype == np.float32


@pytest.mark.parametrize("x_dtype", _DTYPES)
@pytest.mark.parametrize("y_dtype", _DTYPES)
def test_all_four_combinations_produce_equal_cdfs(x_dtype, y_dtype):
    # Integer-valued data round-trips losslessly through f32; CDFs should match
    # the f64-storage baseline exactly.
    X_f64 = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y_f64 = np.array([10.0, 30.0, 20.0, 50.0, 40.0])
    baseline = IDR(y_f64, X_f64).cdf(X_f64)

    X = X_f64.astype(x_dtype)
    y = y_f64.astype(y_dtype)
    fit = IDR(y, X)
    np.testing.assert_array_equal(fit.cdf(X), baseline)


@pytest.mark.parametrize("x_dtype", _DTYPES)
@pytest.mark.parametrize("y_dtype", _DTYPES)
def test_partial_order_all_four_combinations(x_dtype, y_dtype):
    # The multivariate (partial-order) path solves in f64. The wrapper upcasts on
    # entry to the solver and narrows on exit; storage stays at
    # the user's chosen precision.
    X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]], dtype=x_dtype)
    y = np.array([1.0, 2.0, 3.0, 4.0], dtype=y_dtype)
    fit = IDR(y, X)
    assert fit.X.dtype == x_dtype
    assert fit.thresholds.dtype == y_dtype
    assert fit.X.shape == (4, 2)


@pytest.mark.parametrize("x_dtype", _DTYPES)
@pytest.mark.parametrize("y_dtype", _DTYPES)
def test_all_four_combinations_predict(x_dtype, y_dtype):
    X = np.array([1.0, 2.0, 3.0], dtype=x_dtype)
    y = np.array([10.0, 20.0, 30.0], dtype=y_dtype)
    fit = IDR(y, X)
    # predict returns the model's response dtype (Y). Internal accumulation
    # stays in f64 so the f32 result is still accurate to f64 round-off.
    means = fit.predict(X.astype(np.float64))
    assert means.dtype == y_dtype
    np.testing.assert_allclose(means, [10.0, 20.0, 30.0], atol=1e-6)


@pytest.mark.parametrize("x_dtype", _DTYPES)
@pytest.mark.parametrize("y_dtype", _DTYPES)
def test_all_four_combinations_cdf_at(x_dtype, y_dtype):
    X = np.array([1.0, 2.0, 3.0], dtype=x_dtype)
    y = np.array([10.0, 20.0, 30.0], dtype=y_dtype)
    fit = IDR(y, X)
    # cdf_at takes f64 thresholds at the API boundary; output is f32.
    result = fit.cdf_at(X.astype(np.float64), np.array([15.0]))
    assert result.dtype == np.float32


@pytest.mark.parametrize("x_dtype", _DTYPES)
@pytest.mark.parametrize("y_dtype", _DTYPES)
def test_all_four_combinations_quantile(x_dtype, y_dtype):
    X = np.array([1.0, 2.0, 3.0], dtype=x_dtype)
    y = np.array([10.0, 20.0, 30.0], dtype=y_dtype)
    fit = IDR(y, X)
    q = fit.quantile(X.astype(np.float64), np.array([0.5]))
    assert q.dtype == np.float64
    # NumPy broadcasting: (3,) covariate vs (1,) probability collapses to shape (3,)
    assert q.shape == (3,)


@pytest.mark.parametrize("x_dtype", _DTYPES)
@pytest.mark.parametrize("y_dtype", _DTYPES)
def test_pickle_round_trip_all_four_combinations(x_dtype, y_dtype):
    import pickle

    X = np.array([1.0, 2.0, 3.0], dtype=x_dtype)
    y = np.array([4.0, 3.0, 2.0], dtype=y_dtype)
    fit = IDR(y, X)
    restored = pickle.loads(pickle.dumps(fit))
    assert restored.X.dtype == x_dtype
    assert restored.thresholds.dtype == y_dtype
    np.testing.assert_array_equal(restored.cdf(X), fit.cdf(X))


@pytest.mark.parametrize("x_dtype", _DTYPES)
@pytest.mark.parametrize("y_dtype", _DTYPES)
def test_repr_contains_dtype_info(x_dtype, y_dtype):
    # repr and str should expose the storage precision so users can tell at a glance
    # which (X, Y) combo a fit is using.
    X = np.array([1.0, 2.0, 3.0], dtype=x_dtype)
    y = np.array([10.0, 20.0, 30.0], dtype=y_dtype)
    fit = IDR(y, X)
    x_name = "float32" if x_dtype == np.float32 else "float64"
    y_name = "float32" if y_dtype == np.float32 else "float64"
    rep = repr(fit)
    assert f"X={x_name}" in rep
    assert f"y={y_name}" in rep
    s = str(fit)
    assert x_name in s
    assert y_name in s


def test_x_view_is_zero_copy_and_readonly():
    # After the typed-FitImpl rewrite, .X is a direct view into the typed storage —
    # not a freshly allocated cast. The array should be read-only (writes refused)
    # and share memory with subsequent reads.
    X = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    y = np.array([4.0, 3.0, 2.0], dtype=np.float32)
    fit = IDR(y, X)
    view = fit.X
    assert view.dtype == np.float32
    assert not view.flags.writeable
    # The view's data pointer is stable across re-reads (zero-copy).
    view2 = fit.X
    assert view.ctypes.data == view2.ctypes.data


def test_thresholds_view_is_zero_copy_and_readonly():
    X = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    y = np.array([4.0, 3.0, 2.0], dtype=np.float32)
    fit = IDR(y, X)
    view = fit.thresholds
    assert view.dtype == np.float32
    assert not view.flags.writeable
    view2 = fit.thresholds
    assert view.ctypes.data == view2.ctypes.data


def test_f32_round_trip_is_bit_exact_for_subnormal_input():
    # Internal storage now matches the user's input dtype, so any representable f32
    # round-trips bit-exactly through .X and .thresholds without going through f64.
    X = np.array([1e-30, 1.5e-30, 2e-30, 1e-20, 1.5e-20], dtype=np.float32)
    y = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
    fit = IDR(y, X)
    np.testing.assert_array_equal(fit.X, np.sort(np.unique(X)))
    np.testing.assert_array_equal(fit.thresholds, np.sort(np.unique(y)))


@pytest.mark.parametrize(
    "y_dtype, x_dtype",
    [
        (np.float64, np.float32),
        (np.float32, np.float64),
    ],
)
def test_mixed_dtypes_are_independent(y_dtype, x_dtype):
    # X and Y precision are independent — each tracked through its own type
    # parameter in the rust enum. This test makes sure the dispatch picks the
    # right combo (not just a "majority vote").
    X = np.array([1.0, 2.0, 3.0], dtype=x_dtype)
    y = np.array([10.0, 20.0, 30.0], dtype=y_dtype)
    fit = IDR(y, X)
    assert fit.X.dtype == x_dtype
    assert fit.thresholds.dtype == y_dtype


@pytest.mark.parametrize("x_dtype", _DTYPES)
@pytest.mark.parametrize("y_dtype", _DTYPES)
@pytest.mark.parametrize("w_dtype", _DTYPES)
def test_sample_weight_dtype_is_independent(x_dtype, y_dtype, w_dtype):
    # W is now a separate type parameter on the rust `fit` method — passing
    # f32 or f64 weights works for any (X, Y) storage combo. Weights don't
    # leak into the stored types; they only affect the algorithm's internal
    # numerics (narrowed to f32 for total-order, f64 for partial-order).
    X = np.array([1.0, 2.0, 3.0], dtype=x_dtype)
    y = np.array([10.0, 20.0, 30.0], dtype=y_dtype)
    weight = np.array([1.0, 2.0, 3.0], dtype=w_dtype)
    fit = IDR(y, X, sample_weight=weight)
    assert fit.X.dtype == x_dtype
    assert fit.thresholds.dtype == y_dtype


@pytest.mark.parametrize("w_dtype", _DTYPES)
def test_f32_and_f64_weights_produce_equal_cdfs(w_dtype):
    # When the weight values are exactly representable in both precisions, the
    # algorithm must produce identical CDFs regardless of the input weight dtype.
    X = np.array([1.0, 2.0, 3.0, 4.0])
    y = np.array([10.0, 20.0, 30.0, 40.0])
    weight_ref = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
    weight_test = np.array([1.0, 2.0, 3.0, 4.0], dtype=w_dtype)

    fit_ref = IDR(y, X, sample_weight=weight_ref)
    fit_test = IDR(y, X, sample_weight=weight_test)

    np.testing.assert_array_equal(fit_ref.cdf(X), fit_test.cdf(X))


def test_partial_order_accepts_f32_weights():
    # The partial-order solver runs in f64 internally; f32 weights are narrowed
    # at preprocessing. Exercises the partial-order path to make sure W=f32 doesn't
    # break the Vec<f64> Context layout.
    X = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    y = np.array([10.0, 20.0, 30.0])
    weight = np.array([1.0, 1.0, 1.0], dtype=np.float32)
    fit = IDR(y, X, sample_weight=weight)
    means = fit.predict(X.astype(np.float64))
    # Solver convergence tolerance ~1e-5, hence the loose atol.
    np.testing.assert_allclose(means, [10.0, 20.0, 30.0], atol=1e-4)


# Backward-compat: older versions accepted a top-level `epsilon` key for the
# total-order Config. That field is gone and silently ignoring it would let
# users upgrade with no signal that their tuning has stopped taking effect.
# These tests pin down that unknown top-level keys are rejected explicitly.


def test_settings_rejects_epsilon_with_migration_hint():
    X = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 3.0, 2.0])
    with pytest.raises(ValueError, match="epsilon"):
        IDR(y, X, settings={"epsilon": 1e-9})


def test_settings_rejects_unknown_top_level_key():
    X = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 3.0, 2.0])
    with pytest.raises(ValueError, match="unknown settings key"):
        IDR(y, X, settings={"unknown_key": 42})


def test_settings_solver_settings_accepted():
    # Sanity check that the known-key path didn't regress. Use a multivariate X
    # so the partial-order path actually runs and the max_iter knob is exercised.
    X = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    y = np.array([10.0, 20.0, 30.0])
    fit = IDR(y, X, settings={"solver_settings": {"max_iter": 100}})
    assert fit.X.shape == (3, 2)


def test_settings_osqp_settings_names_its_replacement():
    # The key was renamed, not dropped. Silently ignoring the old spelling would
    # leave a caller's tuning quietly unapplied, so it must name the new key.
    X = np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
    y = np.array([10.0, 20.0, 30.0])
    with pytest.raises(ValueError, match="solver_settings"):
        IDR(y, X, settings={"osqp_settings": {"max_iter": 100}})
