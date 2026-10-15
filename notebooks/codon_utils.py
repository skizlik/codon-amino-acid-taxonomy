"""Shared, tested helpers for the codon-usage notebooks.

Data: read the raw UCI file and repair its two comma-shifted rows; codebooks;
the explicit class filter; amino-acid features; cross-validation groups.
Evaluation: models, nested grouped cross-validation, the Nadeau-Bengio
corrected resampled t-test, and summaries.
"""
import csv
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import clone
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
from xgboost import XGBClassifier

DATA = Path(__file__).resolve().parent.parent / "data" / "raw" / "codon_usage.csv"

# Column order of the UCI file.
CODONS = [
    "UUU", "UUC", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "AUU", "AUC", "AUA", "AUG",
    "GUU", "GUC", "GUA", "GUG", "GCU", "GCC", "GCA", "GCG", "CCU", "CCC", "CCA", "CCG",
    "UGG", "GGU", "GGC", "GGA", "GGG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC", "ACU",
    "ACC", "ACA", "ACG", "UAU", "UAC", "CAA", "CAG", "AAU", "AAC", "UGU", "UGC", "CAU",
    "CAC", "AAA", "AAG", "CGU", "CGC", "CGA", "CGG", "AGA", "AGG", "GAU", "GAC", "GAA",
    "GAG", "UAA", "UAG", "UGA",
]

# Synonymous-codon groups under the STANDARD genetic code (NCBI table 1).
# AUG appears once (methionine): the 2022 'start' column duplicated it exactly,
# because aggregate frequencies cannot separate initiator from internal AUGs.
# Animal mitochondria (about 21% of rows) use other codes (UGA = Trp, AUA = Met,
# AGA/AGG = stop or Ser), so for those rows these are codon groups, not
# amino-acid composition. Choosing the right code per row would require the
# compartment and lineage -- the targets -- so the standard code is the only
# label-free choice.
AMINO_CODONS = {
    "alanine": ["GCU", "GCC", "GCA", "GCG"],
    "arginine": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "asparagine": ["AAU", "AAC"],
    "aspartic_acid": ["GAU", "GAC"],
    "cysteine": ["UGU", "UGC"],
    "glutamine": ["CAA", "CAG"],
    "glutamic_acid": ["GAA", "GAG"],
    "glycine": ["GGU", "GGC", "GGA", "GGG"],
    "histidine": ["CAU", "CAC"],
    "isoleucine": ["AUU", "AUC", "AUA"],
    "leucine": ["CUU", "CUC", "CUA", "CUG", "UUA", "UUG"],
    "lysine": ["AAA", "AAG"],
    "methionine": ["AUG"],
    "phenylalanine": ["UUU", "UUC"],
    "proline": ["CCU", "CCC", "CCA", "CCG"],
    "serine": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "threonine": ["ACU", "ACC", "ACA", "ACG"],
    "tryptophan": ["UGG"],
    "tyrosine": ["UAU", "UAC"],
    "valine": ["GUU", "GUC", "GUA", "GUG"],
    "stop": ["UAA", "UAG", "UGA"],
}
AMINOS = list(AMINO_CODONS)
assert sorted(c for cs in AMINO_CODONS.values() for c in cs) == sorted(CODONS)

# UCI codebook (archive.ics.uci.edu/dataset/577/codon+usage).
KINGDOM_CODES = {
    "arc": "archaea", "bct": "bacteria", "phg": "bacteriophage", "plm": "plasmid",
    "pln": "plant", "inv": "invertebrate", "vrt": "vertebrate", "mam": "mammal",
    "rod": "rodent", "pri": "primate", "vrl": "virus",
}
DNATYPE_CODES = {
    0: "genomic", 1: "mitochondrial", 2: "chloroplast", 3: "cyanelle", 4: "plastid",
    5: "nucleomorph", 6: "secondary_endosymbiont", 7: "chromoplast", 8: "leucoplast",
    9: "unknown",  # "NA" in the codebook: keep it as a code, never as a missing value
    10: "proplastid", 11: "apicoplast", 12: "kinetoplast",
}


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def repair_row(fields):
    """Realign a row whose SpeciesName contains unquoted commas.

    The name spills into the codon columns, so the last codon value(s) are
    absent from the file. Frequencies are counts / Ncodons (5 dp), so the
    missing counts are Ncodons minus the counts present. Returns the repaired
    fields and the number of codons recovered; raises ValueError if the
    missing values cannot be identified.
    """
    k = 5
    while not all(_is_number(x) for x in fields[k:]):
        k += 1
    values = np.array(fields[k:], dtype=float)
    n_missing = len(fields) - 5 - len(values)
    ncodons = int(fields[3])
    counts = values * ncodons
    if np.abs(counts - counts.round()).max() > 0.05:
        raise ValueError(f"row {fields[2]}: frequencies are not counts / Ncodons")
    residual = ncodons - int(counts.round().sum())
    if n_missing == 1:
        recovered = [residual / ncodons]
    elif residual == 0:
        recovered = [0.0] * n_missing  # non-negative counts that sum to 0
    else:
        raise ValueError(f"row {fields[2]}: {n_missing} missing codons, residual {residual}")
    name = ",".join(fields[4:k])
    return fields[:4] + [name] + fields[k:] + [f"{v:.5f}" for v in recovered], n_missing


def read_raw(path=DATA):
    """Read codon_usage.csv, repairing comma-shifted rows.

    Adds `repaired` (bool) and the decoded labels `kingdom` and `dna_type`.
    DNAtype stays an integer code (code 9 is not turned into a missing value).
    """
    with open(path, newline="") as fh:
        header, *body = list(csv.reader(fh))
    rows, repaired = [], []
    for fields in body:
        if len(fields) != len(header):
            raise ValueError(f"{len(fields)} fields, expected {len(header)}")
        is_bad = not all(_is_number(x) for x in fields[5:])
        rows.append(repair_row(fields)[0] if is_bad else fields)
        repaired.append(is_bad)
    df = pd.DataFrame(rows, columns=header)
    df[CODONS] = df[CODONS].astype(float)
    df[["DNAtype", "SpeciesID", "Ncodons"]] = df[["DNAtype", "SpeciesID", "Ncodons"]].astype(int)
    df["repaired"] = repaired
    df["kingdom"] = df.Kingdom.map(KINGDOM_CODES)
    df["dna_type"] = df.DNAtype.map(DNATYPE_CODES)
    if df[["kingdom", "dna_type"]].isna().any().any():
        raise ValueError("a Kingdom or DNAtype code is not in the codebook")
    return df


def filter_classes(df):
    """The project's class filter, step by step. Returns (rows kept, row-count log)."""
    steps = [("raw", len(df))]
    df = df[df.kingdom != "plasmid"]
    steps.append(("drop Kingdom 'plm' (plasmid)", len(df)))
    df = df[df.dna_type.isin(["genomic", "mitochondrial", "chloroplast"])]
    steps.append(("keep DNAtype in {0 genomic, 1 mitochondrial, 2 chloroplast}", len(df)))
    # 'chloroplast Toxoplasma gondii', filed under inv: an apicoplast, not a chloroplast
    df = df[~((df.kingdom == "invertebrate") & (df.dna_type == "chloroplast"))]
    steps.append(("drop inv + chloroplast (Toxoplasma apicoplast) row", len(df)))
    log = pd.DataFrame(steps, columns=["step", "rows"])
    log["removed"] = (-log.rows.diff()).fillna(0).astype(int)
    return df.copy(), log


def amino_features(df):
    """Sum codon frequencies within each amino acid (and stop): 21 columns."""
    return pd.DataFrame({a: df[cs].sum(axis=1) for a, cs in AMINO_CODONS.items()}, index=df.index)


_SKIP = {"mitochondrion", "chloroplast", "x", "candidatus", "uncultured", "unidentified", "[candida]"}


def name_token(name):
    """First informative word of SpeciesName: the genus for most cellular
    organisms, often a host word ('human', 'tomato') for viruses."""
    words = [w for w in str(name).replace("[", " [").split() if w.lower() not in _SKIP]
    return words[0].lower() if words else str(name).lower()


def add_groups(df):
    """group_species = SpeciesID (the NCBI taxid); group_genus = name token,
    made constant within SpeciesID so that genus groups contain whole species."""
    df = df.copy()
    df["group_species"] = df.SpeciesID
    df["group_genus"] = df.SpeciesName.map(name_token).groupby(df.SpeciesID).transform("first")
    return df


def load_project(path=DATA):
    """read_raw -> filter_classes -> add_groups. Returns (modelling frame, row-count log)."""
    df, log = filter_classes(read_raw(path))
    return add_groups(df), log


# --- Evaluation ---------------------------------------------------------------

METRICS = ["accuracy", "balanced_accuracy", "f1_macro", "roc_auc_ovr_macro"]


def make_models(seed=0, n_jobs=2):
    """name -> (estimator, grid searched by inner CV). kNN and naive Bayes are
    tuned; RF, XGBoost and MLP use fixed, conventional settings."""
    return {
        "majority": (DummyClassifier(strategy="prior"), {}),
        "gnb": (GaussianNB(), {"var_smoothing": [1e-9, 1e-6, 1e-3]}),
        "knn": (make_pipeline(StandardScaler(), KNeighborsClassifier(weights="distance")),
                {"kneighborsclassifier__n_neighbors": [1, 3, 7, 15]}),
        "rf": (RandomForestClassifier(n_estimators=200, max_features="sqrt", n_jobs=n_jobs,
                                      random_state=seed), {}),
        "xgb": (XGBClassifier(n_estimators=150, learning_rate=0.2, max_depth=6, max_bin=64,
                              subsample=0.8, colsample_bytree=0.8, tree_method="hist",
                              n_jobs=n_jobs, random_state=seed, verbosity=0), {}),
        "mlp": (make_pipeline(StandardScaler(),
                              MLPClassifier(hidden_layer_sizes=(128,), alpha=1e-3, early_stopping=True,
                                            max_iter=300, n_iter_no_change=10, random_state=seed)), {}),
    }


def score_fold(y_true, y_pred, proba, labels):
    """Four metrics for one test fold; ROC-AUC (one-vs-rest, macro) from probabilities."""
    try:
        auc = (roc_auc_score(y_true, proba[:, 1]) if len(labels) == 2 else
               roc_auc_score(y_true, proba, multi_class="ovr", average="macro", labels=labels))
    except ValueError:  # a class is absent from this test fold
        auc = float("nan")
    return {"accuracy": accuracy_score(y_true, y_pred),
            "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
            "f1_macro": f1_score(y_true, y_pred, average="macro"),
            "roc_auc_ovr_macro": auc}


def splitter(scheme, n_splits, seed):
    """StratifiedKFold for 'ungrouped' (leakage comparison only); otherwise
    StratifiedGroupKFold, so no group is in both training and test folds."""
    if scheme == "ungrouped":
        return StratifiedKFold(n_splits, shuffle=True, random_state=seed)
    return StratifiedGroupKFold(n_splits, shuffle=True, random_state=seed)


def nested_cv(X, y, groups, models, scheme="species", n_splits=5, repeats=3,
              inner_splits=3, seed=0, tag=None):
    """Repeated outer CV with tuning by inner grouped CV on the outer-training part only.

    Repeat r uses outer seed `seed + r` and inner seed `seed + 100 + r`, so all
    models and feature sets see identical folds (paired comparisons).
    Returns (fold scores, chosen parameters, out-of-fold predictions).
    """
    X = np.asarray(X, dtype=float)
    enc = LabelEncoder().fit(np.asarray(y))
    y = enc.transform(np.asarray(y))
    labels = np.arange(len(enc.classes_))
    groups = np.asarray(groups)
    tag = tag or {}
    scores, params, oof = [], [], []
    for r in range(repeats):
        outer = splitter(scheme, n_splits, seed + r)
        split_args = (X, y) if scheme == "ungrouped" else (X, y, groups)
        for f, (tr, te) in enumerate(outer.split(*split_args)):
            if scheme != "ungrouped":
                assert not set(groups[tr]) & set(groups[te]), "a group is in both train and test"
            for name, (est, grid) in models.items():
                t0 = time.perf_counter()
                with warnings.catch_warnings():
                    # MLP early stopping can end at max_iter on a few inner folds
                    warnings.simplefilter("ignore", ConvergenceWarning)
                    if grid:
                        inner = splitter(scheme, inner_splits, seed + 100 + r)
                        fit_kw = {} if scheme == "ungrouped" else {"groups": groups[tr]}
                        gs = GridSearchCV(clone(est), grid, scoring="f1_macro", cv=inner, n_jobs=1)
                        model, best = gs.fit(X[tr], y[tr], **fit_kw).best_estimator_, gs.best_params_
                    else:
                        model, best = clone(est).fit(X[tr], y[tr]), {}
                    proba, pred = model.predict_proba(X[te]), model.predict(X[te]).astype(int)
                ids = {**tag, "scheme": scheme, "model": name, "repeat": r, "fold": f}
                scores.append({**ids, "n_train": len(tr), "n_test": len(te),
                               **score_fold(y[te], pred, proba, labels),
                               "seconds": time.perf_counter() - t0})
                params.append({**ids, "params": repr(best)})
                oof.append(pd.DataFrame({**ids, "row": te, "y_true": enc.classes_[y[te]],
                                         "y_pred": enc.classes_[pred]}))
    return pd.DataFrame(scores), pd.DataFrame(params), pd.concat(oof, ignore_index=True)


def corrected_ttest(a, b, n_train, n_test):
    """Nadeau & Bengio (2003) corrected resampled t-test for J = k x r paired fold scores.

    The variance of the mean difference is var/J for independent folds; training
    sets overlap, so it is inflated to (1/J + n_test/n_train) * var.
    Returns (mean difference, t, two-sided p).
    """
    d = np.asarray(a, float) - np.asarray(b, float)
    J, var = d.size, d.var(ddof=1)
    if var == 0:
        return d.mean(), float("nan"), float("nan")
    t = d.mean() / np.sqrt((1 / J + n_test / n_train) * var)
    return d.mean(), t, 2 * stats.t.sf(abs(t), df=J - 1)


def corrected_ci(a, b, n_train, n_test, level=0.95):
    """Confidence interval for mean(a - b) with the same variance correction."""
    d = np.asarray(a, float) - np.asarray(b, float)
    se = np.sqrt((1 / d.size + n_test / n_train) * d.var(ddof=1))
    h = stats.t.ppf(0.5 + level / 2, d.size - 1) * se
    return d.mean() - h, d.mean() + h


def summarize(scores, by=("scheme", "task", "features", "model")):
    """Mean and SD (over folds) of each metric."""
    g = scores.groupby(list(by), sort=False)[METRICS]
    out = pd.concat([g.mean().add_suffix("_mean"), g.std(ddof=1).add_suffix("_sd")], axis=1)
    out["n_folds"] = g.size()
    return out[[f"{m}_{s}" for m in METRICS for s in ("mean", "sd")] + ["n_folds"]].reset_index()
