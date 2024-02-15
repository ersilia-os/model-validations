# hERG Performance Analysis
## Overview
The hERG Performance Analysis repository is a dedicated to validating the various hERG models hosted in the Ersilia model hub. This README provides detailed information on the models included in the validation, the data acquisition process, data pre-processing techniques, and the construction of an unbiased tests dataset for evaluation.

## Models Used in Validation
The following ersilia hERG models were included in the validation process:

eos2ta5 - [Repository Link](https://github.com/ersilia-os/eos2ta5)

eos4tcc - [Repository Link](https://github.com/ersilia-os/eos4tcc)

eos30f3 - [Repository Link](https://github.com/ersilia-os/eos30f3)

eos30gr - [Repository Link](https://github.com/ersilia-os/eos30gr)



## Data Acquisition

### Training dataset
The datasets in the model_datasets folder comprise the training data utilized for each hERG model, sourced directly from the respective model's codebase.

For the eos2ta5 model, the training dataset was retrieved from the source code available [here](https://github.com/Abdulk084/CardioTox/tree/master/data), incorporating files  `external_test_set_neg.csv`, `external_test_set_pos.csv`, `external_test_set_new.csv`, and `train_validation_cardio_tox_data.tar.xz`.

Regarding the eos4tcc model, its training dataset consists of files obtained from the finetuning directory within the source code accessible [here](https://github.com/GIST-CSBL/BayeshERG/tree/main/data/Finetuning).

The eos30f3 model was trained using a dataset introduced by Cai et al. in their work published in J Chem Inf Model, 2019. This dataset comprises 7889 molecules with various cutoffs for hERG blocking activity, available [here](https://github.com/AI-amateur/DMPNN-hERG/blob/main/Figure2/Cai_TableS3_fixed.csv)

Lastly, for the eos30gr model, the training dataset was directly provided within the source code, accessible [here](https://github.com/ChengF-Lab/deephERG/blob/master/Table%20S6.xlsx)

These datasets undergo preprocessing where we clean up the original files and add the InChiKey of the smiles if not available. After preprocessing, each dataset is stored under data/model_datasets/{model_name}_processed.csv. These processed files are then concatenated and duplicates are removed based on Inchikeys to create our training_set.csv

### Test dataset
The test dataset was obtained from the [NCATS repository](https://github.com/ncats/herg-ml/tree/master/data/train_valid). Molecules already present in the training_set.csv were eliminated. This is to ensure fair and unbiased evaluation of the hERG models included in the validation process.


## Repository Structure
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
