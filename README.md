# Codon Usage Analysis: Predicting Taxonomic Groups and DNA Types

> **Erratum (2026):** The model-comparison numbers originally published here were affected by a
> results-table join error and two evaluation defects. They are withdrawn pending a corrected
> analysis. See [ERRATUM.md](ERRATUM.md).

A machine learning analysis of codon usage patterns to predict taxonomic classification and DNA type using multiple algorithms including k-Nearest Neighbors, Random Forest, XGBoost, Multi-Layer Perceptron, and Naive Bayes.

## Project Overview

This project replicates and extends the findings of Khomtchouk (2020) on using codon frequency analysis for taxonomic prediction. The study demonstrates that both taxonomic group and DNA type can be predicted with high accuracy from simple frequency analysis of codons in DNA samples.

**Key Innovation**: This analysis extends beyond codon frequencies to examine whether amino acid frequencies (derived from the codons) can provide comparable predictive power, potentially offering computational advantages in certain applications.

## Dataset

- **Source**: UCI Machine Learning Repository - Codon Usage Dataset
- **Size**: 13,028 observations with 69 attributes
- **Features**: Relative frequencies of all 64 codons plus organism metadata
- **Targets**: Taxonomic group ('kingdom') and DNA type classification

### Data Preprocessing

- Addressed class imbalance by removing extremely sparse categories
- Created amino acid frequency features from codon data
- Applied 85/15 train-test splits with separate validation sets
- One-hot encoding for categorical variables

The original 2022 project report is `docs/codon-amino-acid-taxonomy-report.pdf`. Its model-comparison results are withdrawn (see [ERRATUM.md](ERRATUM.md)).

## Installation and Setup

```bash
git clone https://github.com/skizlik/codon-amino-acid-taxonomy.git
cd codon-amino-acid-taxonomy
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

The 2022 notebooks predate pandas 2.0 and XGBoost 2.0 and do not run unchanged on current versions.

## Key Findings

Withdrawn. The model-performance table and the codon-versus-amino-acid comparison previously shown here
were affected by the defects described in [ERRATUM.md](ERRATUM.md). A corrected analysis
(grouped, repeated cross-validation with nested tuning) will replace this section.

## Technical Approach

### Machine Learning Pipeline
1. **Exploratory Data Analysis**: Distribution analysis, correlation matrices, PCA
2. **Feature Engineering**: Conversion of codon frequencies to amino acid frequencies
3. **Model Training**: Five different algorithms with hyperparameter optimization
4. **Evaluation**: Single train/validation/test split with F1-score, accuracy, and ROC-AUC metrics
5. **Feature Importance**: SHAP analysis and Gini impurity rankings

### Hyperparameter Optimization
- Manual tuning for k-NN and Random Forest (settings were chosen on the test split; see ERRATUM.md)
- Automated optimization using Hyperopt for XGBoost and MLP

## Notebooks Description

1. **Notebook 1**: Data loading and initial exploration
2. **Notebook 2**: Exploratory data analysis and visualization
3. **Notebook 3**: Data preprocessing and feature engineering
4. **Notebook 4**: Machine learning model training and evaluation
5. **Notebook 5**: Feature importance analysis and results visualization

## Key Dependencies

- pandas: Data manipulation and analysis
- numpy: Numerical computing
- scikit-learn: Machine learning algorithms
- xgboost: Gradient boosting implementation
- hyperopt: Hyperparameter optimization
- shap: Model interpretability
- matplotlib/seaborn: Data visualization

## Results and Implications

### Biological Insights
- Certain codons (particularly CUA and UGA) show high predictive importance across models
![SHAP analysis of codon and amino acid frequencies to predict DNA.  Most predictive codons were UGA, CUA, and UUA.  Most predictive amino acids were aspartic acid, leucine, and isoleucine. Results were for MLP classifier.](assets/SHAP-MLP-DNA.png)
- Taxonomic relationships are reflected in codon usage patterns
- DNA type classification achieves higher accuracy than taxonomic classification

### Methodological Contributions
- Examines amino acid frequencies as an alternative to codon frequencies
![Boxplots of relative frequency of leucine by taxonomic group, all DNA types pooled. The vertebrate-invertebrate difference is mostly a compartment effect: 78% of vertebrate rows are mitochondrial, and mitochondrial rows have higher leucine frequencies.](assets/leucine_boxplot.png)

## Future Work

- Corrected re-analysis (the anomalous XGBoost results are explained in ERRATUM.md)
- Expand analysis to include sparse taxonomic categories
- Implement ensemble methods combining top-performing algorithms
- Investigate correlation patterns between related amino acids and codons

## References

Khomtchouk, B.B. (2020). Codon Usage Bias Levels Predict Taxonomic Identity and Genetic Composition. BioRxiv. doi: 10.1101/2020.10.26.356295

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

Stephen Kizlik · stephen.kizlik@gmail.com

Project link: https://github.com/skizlik/codon-amino-acid-taxonomy
