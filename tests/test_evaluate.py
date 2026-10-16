"""Grouped folds, nested selection, metrics, the corrected t-test, and the 2022 merge bug."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.datasets import make_classification

import codon_utils as cu


@pytest.mark.parametrize("col", ["group_species", "group_genus"])
def test_outer_folds_never_share_a_group(clean, col):
    df, _ = clean
    y, g = df.kingdom.to_numpy(), df[col].to_numpy()
    for train, test in cu.splitter("species", 5, seed=0).split(np.zeros(len(y)), y, g):
        assert not set(g[train]) & set(g[test])
        assert set(y[test]) == set(y)  # every class in every test fold


class Recorder(ClassifierMixin, BaseEstimator):
    """Predicts the class priors and records the row ids (column 0) of every fit."""

    fits = []

    def __init__(self, c=1):
        self.c = c

    def fit(self, X, y):
        Recorder.fits.append(set(X[:, 0].astype(int)))
        self.classes_ = np.unique(y)
        self.prior_ = np.bincount(y) / len(y)
        return self

    def predict_proba(self, X):
        return np.tile(self.prior_, (len(X), 1))

    def predict(self, X):
        return self.classes_[self.predict_proba(X).argmax(1)]


def toy(n=300):
    X, y = make_classification(n, 6, n_informative=4, n_classes=3, random_state=0)
    return np.column_stack([np.arange(n), X]), y, np.repeat(np.arange(n // 3), 3)


def test_tuning_never_sees_the_outer_test_fold():
    X, y, groups = toy()
    Recorder.fits = []
    _, _, oof = cu.nested_cv(X, y, groups, {"rec": (Recorder(), {"c": [1, 2]})}, n_splits=3, repeats=2)
    # per outer fold: 2 settings x 3 inner folds, then 1 refit = 7 fits
    assert len(Recorder.fits) == 2 * 3 * 7
    test_sets = [set(f.row) for _, f in oof.groupby(["repeat", "fold"])]
    for k, test in enumerate(test_sets):
        for fitted_rows in Recorder.fits[7 * k : 7 * (k + 1)]:
            assert not fitted_rows & test


def test_every_row_predicted_once_per_repeat_on_reproducible_folds():
    X, y, groups = toy()
    labels = np.array(["a", "b", "c"])[y]
    models = {"majority": cu.make_models()["majority"]}
    scores, _, oof = cu.nested_cv(X, labels, groups, models, n_splits=3, repeats=2)
    assert len(scores) == 6 and scores.roc_auc_ovr_macro.eq(0.5).all()
    for _, r in oof.groupby("repeat"):
        assert sorted(r.row) == list(range(len(y)))
    assert set(oof.y_true) == {"a", "b", "c"}
    again, _, _ = cu.nested_cv(X, labels, groups, models, n_splits=3, repeats=2)
    assert again.n_test.tolist() == scores.n_test.tolist()  # same seed, same folds


def test_auc_uses_probabilities_not_hard_labels():
    y = np.array([0, 0, 1, 1])
    proba = np.array([[0.6, 0.4], [0.7, 0.3], [0.45, 0.55], [0.55, 0.45]])
    s = cu.score_fold(y, proba.argmax(1), proba, labels=np.array([0, 1]))
    assert s["roc_auc_ovr_macro"] == 1.0  # the ranking is perfect...
    assert s["accuracy"] == 0.75  # ...although one hard label is wrong


def test_corrected_ttest_formula_and_ci():
    rng = np.random.default_rng(1)
    a = 0.90 + rng.normal(0, 0.01, 15)
    b = a - 0.004 + rng.normal(0, 0.005, 15)
    d = a - b
    mean, t, p = cu.corrected_ttest(a, b, n_train=8000, n_test=2000)
    assert t == pytest.approx(d.mean() / np.sqrt((1 / 15 + 2000 / 8000) * d.var(ddof=1)))
    assert p == pytest.approx(2 * stats.t.sf(abs(t), df=14))
    assert p > stats.ttest_rel(a, b).pvalue  # the correction is conservative
    lo, hi = cu.corrected_ci(a, b, 8000, 2000)
    assert lo < mean < hi
    assert (lo > 0 or hi < 0) == (p < 0.05)  # the CI excludes 0 exactly when p < 0.05


def test_index_merge_after_dropping_rows_misaligns_metrics():
    """The 2022 results-table bug (notebook 05, cells 33-52), in miniature."""
    models = ["kNN", "kNN", "SHAP", "XGB", "SHAP", "XGB-Hyperopt", "GNB"]
    results = pd.DataFrame(
        {"model_type": models, "assessment": [{"test_f1": i / 10, "fitted": m} for i, m in enumerate(models)]}
    )
    results = results[results.model_type != "SHAP"]  # index labels 0, 1, 3, 5, 6 remain
    assess = pd.DataFrame(list(results.assessment))  # fresh labels 0..4: the bug
    buggy = pd.merge(results, assess, left_index=True, right_index=True)
    assert len(buggy) == 3  # XGB-Hyperopt and GNB silently dropped
    assert buggy.set_index("model_type").loc["XGB", "fitted"] == "XGB-Hyperopt"  # wrong metrics
    # safe pattern: build the metrics on the same index, and make pandas check the join
    assess = pd.DataFrame(list(results.assessment), index=results.index)
    fixed = results.join(assess, validate="one_to_one")
    assert len(fixed) == len(results) and (fixed.model_type == fixed.fitted).all()


def test_colsample_bytree_zero_grows_single_feature_trees():
    """ERRATUM item 2: the 2022 search passed int(colsample_bytree), which is always 0."""
    from xgboost import XGBClassifier

    X, y = make_classification(400, 10, n_informative=6, n_classes=3, random_state=0)

    def max_features_per_tree(colsample):
        model = XGBClassifier(n_estimators=20, max_depth=4, colsample_bytree=colsample, verbosity=0).fit(X, y)
        trees = model.get_booster().trees_to_dataframe()
        return trees[trees.Feature != "Leaf"].groupby("Tree").Feature.nunique().max()

    assert int(0.61) == 0  # what int() did to every value drawn from U(0.5, 1)
    assert max_features_per_tree(0) == 1  # accepted silently: one feature per tree
    assert max_features_per_tree(0.8) > 1
