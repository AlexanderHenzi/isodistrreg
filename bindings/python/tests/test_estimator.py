"""Behaviour of the scikit-learn estimator beyond the generic API checks.

The generic checks live in ``test_scikit_learn_api_adherence.py``; these tests
pin down what the estimator does with its own parameters and methods.
"""

import subprocess
import sys

import numpy as np
import pytest

pytest.importorskip("sklearn")

import sklearn
from sklearn.base import clone
from sklearn.exceptions import NotFittedError
from sklearn.model_selection import GridSearchCV, KFold, cross_validate
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils._param_validation import InvalidParameterError

from isodistrreg import IDR, IsotonicDistributionalRegressor


@pytest.fixture
def data():
    rng = np.random.default_rng(0)
    X = rng.uniform(size=(40, 1))
    return X, X[:, 0] + rng.uniform(size=40)


@pytest.fixture
def data2d():
    rng = np.random.default_rng(1)
    X = rng.uniform(size=(40, 2))
    return X, X.sum(1) + rng.uniform(size=40)


class TestParameters:
    def test_hyperparameters_live_in_init(self):
        est = IsotonicDistributionalRegressor(decreasing=True, subsamples=3)
        params = est.get_params()
        assert params["decreasing"] is True
        assert params["subsamples"] == 3
        assert clone(est).get_params() == params
        assert est.set_params(decreasing=False).decreasing is False

    def test_grid_search_over_hyperparameters(self, data):
        X, y = data
        search = GridSearchCV(
            IsotonicDistributionalRegressor(),
            {"decreasing": [False, True], "response_order": ["sd", "hazard"]},
            cv=3,
        ).fit(X, y)
        assert search.best_params_["decreasing"] is False

    @pytest.mark.parametrize(
        "params",
        [
            {"decreasing": 1.0},
            {"decreasing": "yes"},
            {"response_order": "up"},
            {"subsamples": 0},
            {"subsample_size": 0.0},
            {"subsample_size": "half"},
            {"settings": []},
            {"covariate_order": "comp"},
            {"random_state": -1},
            {"n_jobs": 1.5},
        ],
    )
    def test_invalid_parameters_are_rejected_at_fit(self, data, params):
        X, y = data
        est = IsotonicDistributionalRegressor(**params)
        with pytest.raises(InvalidParameterError, match=next(iter(params))):
            est.fit(X, y)

    def test_random_state_reproduces_subagging(self, data):
        X, y = data
        fits = [
            IsotonicDistributionalRegressor(
                subsamples=5, subsample_size=0.5, random_state=1
            ).fit(X, y)
            for _ in range(2)
        ]
        np.testing.assert_array_equal(fits[0].predict(X), fits[1].predict(X))
        other = IsotonicDistributionalRegressor(
            subsamples=5, subsample_size=0.5, random_state=2
        ).fit(X, y)
        assert not np.array_equal(fits[0].predict(X), other.predict(X))

    def test_random_state_instance_is_consumed(self, data):
        X, y = data
        rng = np.random.RandomState(0)
        est = IsotonicDistributionalRegressor(subsamples=5, random_state=rng)
        first = est.fit(X, y).predict(X)
        second = est.fit(X, y).predict(X)
        assert not np.array_equal(first, second)

    def test_plain_fit_leaves_global_rng_alone(self, data):
        X, y = data
        np.random.seed(0)
        before = np.random.get_state()[1]
        IsotonicDistributionalRegressor().fit(X, y)
        np.testing.assert_array_equal(np.random.get_state()[1], before)

    def test_n_jobs_all_processors(self, data):
        X, y = data
        est = IsotonicDistributionalRegressor(subsamples=4, n_jobs=-1, random_state=0)
        assert est.fit(X, y).predict(X[:2]).shape == (2,)


class TestFit:
    def test_sets_fitted_attributes(self, data):
        X, y = data
        est = IsotonicDistributionalRegressor()
        with pytest.raises(NotFittedError):
            est.predict(X)
        est.fit(X, y)
        assert est.n_features_in_ == 1
        assert est.X_unique_.shape == (len(np.unique(X[:, 0])), 1)
        np.testing.assert_array_equal(est.thresholds_, np.unique(y))

    def test_sample_weight_and_y_observed(self, data):
        X, y = data
        est = IsotonicDistributionalRegressor().fit(
            X, y, sample_weight=np.ones(40), y_observed=np.ones(40, dtype=bool)
        )
        assert est.predict(X[:2]).shape == (2,)
        with pytest.raises(ValueError, match="Negative values"):
            IsotonicDistributionalRegressor().fit(X, y, sample_weight=-np.ones(40))
        with pytest.raises(ValueError, match="0/1"):
            IsotonicDistributionalRegressor().fit(X, y, y_observed=2 * np.ones(40))

    @pytest.mark.parametrize("dtype", ["timedelta64[D]", "datetime64[D]"])
    def test_temporal_input_is_rejected_with_the_core_hint(self, data, dtype):
        X, y = data
        temporal = (y * 100).astype("int64").astype(dtype)
        with pytest.raises(ValueError, match=r"np\.timedelta64\(1"):
            IsotonicDistributionalRegressor().fit(X, temporal)
        with pytest.raises(ValueError, match=r"np\.timedelta64\(1"):
            IsotonicDistributionalRegressor().fit(temporal.reshape(-1, 1), y)

    @pytest.mark.parametrize(
        "indicator", [np.ones(40) + 1, np.array(["yes"] * 40), np.array([1.5] * 40)]
    )
    def test_event_indicator_values_are_checked_by_the_core(self, data, indicator):
        X, y = data
        with pytest.raises(ValueError, match="0/1"):
            IsotonicDistributionalRegressor().fit(X, y, y_observed=indicator)

    def test_non_numeric_y_is_rejected(self, data):
        X, _ = data
        with pytest.raises(ValueError, match="could not convert|numeric"):
            IsotonicDistributionalRegressor().fit(X, np.array(["a"] * 40))

    def test_sample_weight_reaches_fit_under_metadata_routing(self, data):
        X, y = data
        w = np.arange(40) % 2 == 0  # drop every other row from the fit
        cv = KFold(3)
        with sklearn.config_context(enable_metadata_routing=True):
            est = (
                IsotonicDistributionalRegressor()
                .set_fit_request(sample_weight=True)
                .set_score_request(sample_weight=False)
            )
            result = cross_validate(
                est, X, y, cv=cv, params={"sample_weight": w}, return_estimator=True
            )
        for fitted, (train, _) in zip(result["estimator"], cv.split(X)):
            manual = IsotonicDistributionalRegressor().fit(
                X[train], y[train], sample_weight=w[train]
            )
            np.testing.assert_array_equal(fitted.predict(X), manual.predict(X))

    def test_covariate_order_on_a_single_feature(self, data):
        X, y = data
        est = IsotonicDistributionalRegressor(covariate_order=[("comp", [0])])
        assert est.fit(X, y).predict(X[:2]).shape == (2,)

    def test_numpy_integer_covariate_order(self):
        rng = np.random.default_rng(1)
        X = rng.uniform(size=(20, 3))
        est = IsotonicDistributionalRegressor(
            covariate_order=[("comp", list(np.array([0, 1, 2])))]
        ).fit(X, X.sum(1))
        assert est.predict(X[:2]).shape == (2,)

    @pytest.mark.parametrize(
        "order, match",
        [
            ([("comp", [0], "extra")], "pair"),
            ([(1, [0])], "kind"),
            ([("comp", [0.0])], "ints"),
            ([("comp", [True])], "ints"),
            ([("lex", [0, 1])], "covariate groups"),
            ([("comp", [0, 5])], "covariate groups"),
            ([("comp", [0, 0])], "covariate groups"),
        ],
    )
    def test_invalid_covariate_order(self, data2d, order, match):
        X, y = data2d
        with pytest.raises(ValueError, match=match):
            IsotonicDistributionalRegressor(covariate_order=order).fit(X, y)


class TestPredictionMethods:
    def test_cdf_matches_core_grid(self, data):
        X, y = data
        est = IsotonicDistributionalRegressor().fit(X, y)
        out = est.cdf(X[:5])
        assert out.shape == (5, len(est.thresholds_))
        np.testing.assert_array_equal(out, IDR(y, X).cdf(X[:5]))

    def test_cdf_at_evaluates_every_row_at_every_threshold(self, data):
        X, y = data
        est = IsotonicDistributionalRegressor().fit(X, y)
        thresholds = [0.2, 0.5, 0.8]
        out = est.cdf_at(X[:5], thresholds)
        assert out.shape == (5, 3)
        assert np.all(np.diff(out, axis=1) >= 0)
        for j, t in enumerate(thresholds):
            np.testing.assert_array_equal(out[:, j], est.cdf_at(X[:5], t))
        assert est.cdf_at(X[:5], 0.5).shape == (5,)

    def test_quantile_evaluates_every_row_at_every_level(self, data):
        X, y = data
        est = IsotonicDistributionalRegressor().fit(X, y)
        levels = [0.1, 0.5, 0.9]
        # n_samples == n_levels used to be silently paired elementwise.
        out = est.quantile(X[:3], levels)
        assert out.shape == (3, 3)
        assert np.all(np.diff(out, axis=1) >= 0)
        for j, q in enumerate(levels):
            np.testing.assert_array_equal(out[:, j], est.quantile(X[:3], q))
        assert est.quantile(X[:5], levels).shape == (5, 3)
        assert est.quantile(X[:5], 0.5).shape == (5,)
        np.testing.assert_array_equal(
            est.quantile(X[:5], 0.5), IDR(y, X).quantile(X[:5], 0.5)
        )

    def test_quantile_upper(self, data):
        X, y = data
        est = IsotonicDistributionalRegressor().fit(X, y)
        assert np.all(est.quantile(X, 0.5, upper=True) >= est.quantile(X, 0.5))

    def test_quantile_multivariate(self, data2d):
        X, y = data2d
        est = IsotonicDistributionalRegressor().fit(X, y)
        assert est.quantile(X[:4], [0.25, 0.75]).shape == (4, 2)

    def test_grid_arguments_are_validated(self, data):
        X, y = data
        est = IsotonicDistributionalRegressor().fit(X, y)
        with pytest.raises(ValueError, match="between 0 and 1"):
            est.quantile(X[:2], [0.5, 1.5])
        with pytest.raises(ValueError, match="1-D"):
            est.quantile(X[:2], [[0.5]])
        with pytest.raises(ValueError, match="1-D"):
            est.cdf_at(X[:2], [[0.5]])

    def test_pipeline_with_dataframe(self, data):
        pd = pytest.importorskip("pandas")
        X, y = data
        frame = pd.DataFrame(X, columns=["x"])
        pipe = make_pipeline(StandardScaler(), IsotonicDistributionalRegressor())
        assert pipe.fit(frame, y).predict(frame.iloc[:3]).shape == (3,)
        est = IsotonicDistributionalRegressor().fit(frame, y)
        np.testing.assert_array_equal(est.feature_names_in_, ["x"])
        with pytest.raises(ValueError, match="feature names"):
            est.cdf_at(frame.rename(columns={"x": "z"}), 0.5)


def test_missing_scikit_learn_is_reported_by_name():
    script = (
        "import sys; sys.modules['sklearn'] = None\n"
        "import isodistrreg\n"
        "try:\n"
        "    from isodistrreg import IsotonicDistributionalRegressor\n"
        "except ImportError as e:\n"
        "    print(e)\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=True
    )
    assert "isodistrreg[scikit-learn]" in out.stdout
