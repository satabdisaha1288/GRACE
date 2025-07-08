# 📊 Graph-Adaptive Horseshoe for Compositional Regression

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

## 📖 Citation

If you use this code, please cite:


## 📬 Contact

For questions or collaborations, reach out to [Satabdi Saha](mailto:ssaha1@mdanderson.org) 


Simulate_Pure_Communities.R - This file is used to simulate data from a regression model with two distinct clusters of coefficients,one with positive 
effects and the other with negative effects.

Functions_Graph_Horseshoe.R - This file contains all the functions required to execute the compositional variable‐selection algorithm, which infers 
the weight graph among coefficient scales. It’s designed to run when no prior knowledge of predictor relationships or their effects on the outcome 
is available.

Graph-Horseshoe-Informed-Priors.R - This file contains all the functions required to execute the compositional variable‐selection algorithm, which infers 
the weight graph among coefficient scales. It’s designed to run when prior knowledge of predictor relationships or knowledge about their effects on the 
outcome is available.

Run-Sampler-Graph-Results.R - This file is used to generate the estimated graph output from the algorithm. It computes the inferred graph structure and,
if the true graph is available, calculates the Adjusted Rand Index to assess clustering accuracy. Additionally, it includes functions to compute the
number of false positive and false negative edges, as well as evaluation metrics such as AUROC, AUPR, F1 score, and Matthews Correlation Coefficient (MCC).

Run-Sampler-Prediction-Results.R - This file is used to generate the prediction outputs from the algorithm. It calculates the prediction Root Mean Squared 
Error (RMSE) using a test dataset and computes the L2 loss between the true and estimated regression coefficients. It also evaluates model performance by 
reporting the number of false positive and false negative coefficient selections.
