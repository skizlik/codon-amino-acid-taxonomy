# Erratum: results published 2022–2026

The model-comparison results originally reported in this repository (README tables,
`data/rawbest.csv`, the 2022 report and slides; preserved at tag
[`v0-2022`](https://github.com/skizlik/codon-amino-acid-taxonomy/tree/v0-2022)) contain errors.
They are withdrawn. A corrected analysis, produced by tested code in this repository, will
replace them in the [README](README.md).

## What was wrong

1. **Mislabelled XGBoost results.** When the results table was flattened (notebook 05,
   cells 33–52), metrics were merged on index labels after rows had been dropped. As a result,
   rows labelled "XGB" carried the XGB-Hyperopt metrics and rows labelled "XGB-Hyperopt"
   carried the Gaussian naive Bayes metrics. This is why the 2022 report found the metrics of
   XGB-Hyperopt and naive Bayes "exactly identical". The XGBoost-default metrics are not
   recoverable. The published "XGBoost" row also paired XGB-Hyperopt scores with the
   XGB-default runtime.
2. **Hyperparameter search defect.** The XGBoost search cast `colsample_bytree` to an integer,
   which set it to 0 on every trial; XGBoost then grows each tree on a single feature. One
   sampled parameter (`reg_lambda`) was never passed, and the ranges for `gamma` and
   `reg_alpha` excluded their defaults. The "XGBoost defaults" recorded for comparison used
   `max_depth = 3`; the library default is 6.
3. **Evaluation design.** k-NN and random forest settings were chosen on the test set;
   k-NN and random forest were fitted to one-hot targets (multilabel) while the other models
   used class labels, and their ROC-AUC was computed from hard predictions. There was no
   cross-validation and no grouping of related organisms.
4. **Data.** Two records whose species names contain commas were misparsed, shifting all of
   their codon values. An amino-acid feature (`start`) duplicated methionine. Two records with
   DNA-type code 9 were silently lost as missing values, so 65 rows were removed in total, not
   62 as the report states.

## How the join error mapped rows to metrics

`results` carried index labels 0–111. Notebook 05 cell 33 dropped the 12 SHAP rows
(labels 84–87, 92–95, 104–107); cell 42 built the metrics frame with a fresh 0–99 index;
cell 52 merged the two on index labels.

| rows (labels) | model shown | metrics actually attached |
|---|---|---|
| 0–39 | kNN | kNN |
| 40–83 | RFC | RFC |
| 88–91 | XGB | **XGB-Hyperopt** |
| 96–99 | XGB-Hyperopt | **GNB** |
| 100–103, 108–111 | MLP-Hyperopt, GNB | dropped by the inner join (re-added correctly later) |
