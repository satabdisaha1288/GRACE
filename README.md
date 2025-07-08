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
