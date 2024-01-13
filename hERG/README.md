# hERG Performance Analysis

Explore the hERG models in the Ersilia Model Hub performance through a comprehensive analysis using an unbiased dataset. The test dataset has been obtained from the [NCATS repository](https://github.com/ncats/herg-ml/tree/master/data/train_valid). Molecules in the test set that existed in the training sets of the available models have been removed.

## Data folder

- **Combined Training Datasets:** Find the combined training datasets of hERG models in the hub in the 'data' folder, named `combined_dataset_eos.csv`.

- **Unbiased Test Dataset:** Access the unbiased dataset with a 10 μM cut-off in the 'data' folder, named `test_dataset.csv`.

## Notebooks folder

### 1. Comparing Datasets
- File: `compare_datasets.ipynb`
- Description: This notebook calculates and prints the overlap percentage of Ersilia hERG model datasets.

### 2. Comparing EOS and Test Data
- File: `compare_eos_and_test.ipynb`
- Description: This notebook compares the unbiased dataset with Ersilia hERG training datasets.

### 3. Model Evaluation
- File: `evaluate_model.ipynb`
- Description: Assess individual model performance. This notebook provides accuracy, precision, recall, and F1 score, along with an AUROC curve.

### 4. Scatter Plots and Correlation
- File: `scatter_plots.ipynb`
- Description: Explore the relationships between model performances through scatter plots. This notebook also generates a correlation matrix for a deeper understanding.

Dive into the `#hERG_Datasets` folder for a comprehensive understanding of hERG model performances in the Ersilia model hub.

