# 📊 GRaph Adaptive horseshoe for Compositional REgression (GRACE)

This repository provides the implementation for the **Graph-Adaptive Horseshoe** model, a Bayesian variable selection method tailored for **high-dimensional compositional data** such as microbiome relative abundances.

## 🔍 Overview

This method extends traditional log-contrast regression by:

- Enforcing the **sum-to-zero constraint** through a novel parameterization;
- Learning the **feature similarity graph** directly from data, rather than relying on fixed structures like phylogenetic trees;
- Applying a **graph-structured horseshoe prior** for **sparse and smooth shrinkage**;
- Using a custom **Gibbs sampler with elliptical slice sampling** for scalable inference.

The model improves prediction and feature selection by identifying **outcome-relevant clusters** of covariates. We validate the method through simulations and apply it to a real oral microbiome dataset to uncover associations with **insulin resistance**.

## 📂 Contents

- `code/` – Core R implementation
- `data/` – Simulated and example microbiome datasets
- `scripts/` – Analysis and plotting scripts
- `results/` – Output from simulation and real data analysis

## 🧪 Data Simulation

- **`Simulate_Pure_Communities.R`**  
  Simulates data from a regression model where the **true coefficients form two distinct clusters**:
  - One group of features with **positive effects**,
  - Another group with **negative effects**.

This setup is designed to test the model’s ability to recover outcome-relevant clusters and validate the accuracy of graph-based shrinkage and feature selection.

## 🧾 Code Files

- **`Functions_Graph_Horseshoe.R`**  
  Contains all functions required to execute the **compositional variable-selection algorithm** in the **absence of prior knowledge**. This version **learns the graph of coefficient scales from the data**, adapting to unknown relationships between predictors and their effects on the outcome.

- **`Graph-Horseshoe-Informed-Priors.R`**  
  Contains the corresponding functions for the **informed version** of the compositional variable-selection algorithm. This version assumes **prior knowledge of the relationships among predictors** (e.g., from phylogeny or co-abundance) and **centers the prior distribution of edge weights** on those known relationships while still allowing data-adaptive refinement.

## 📜 Analysis Scripts

- **`Run-Sampler-Graph-Results.R`**  
  Generates the **inferred graph structure** from the fitted model. If a ground-truth graph is available, the script also computes the **Adjusted Rand Index (ARI)** to assess clustering accuracy. Additionally, it includes evaluation metrics such as:
  - Number of **false positives** and **false negatives** (edges),
  - **AUROC**, **AUPR**, **F1 score**, and **Matthews Correlation Coefficient (MCC)**.

- **`Run-Sampler-Prediction-Results.R`**  
  Computes **prediction metrics** based on the fitted model:
  - **Root Mean Squared Error (RMSE)** on a test set,
  - **L2 loss** between the true and estimated regression coefficients,
  - Number of **false positives** and **false negatives** in coefficient selection.

These scripts are intended for post-processing MCMC outputs and evaluating model performance under both simulated and real data settings.


## 📖 Citation

If you use this code, please cite:


## 📬 Contact

For questions or collaborations, reach out to [Satabdi Saha](mailto:ssaha1@mdanderson.org) 






