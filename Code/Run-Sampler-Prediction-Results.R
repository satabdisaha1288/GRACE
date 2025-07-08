# Load the data
# p = 100
#load("/Users/satab/Documents/Horseshoe/Simulations/Data/Pure-Communities_n200_p100_R2_0.89_simdata.RData")
dat <- data.frame( y =sim200$y, X = sim200$X_clr)
n <- nrow(dat)
q    <- n/2
p<- ncol(sim200$X_clr)
W_true <- sim200$W_xy
W_x <-sim200$W_x
dat1 <- dat[1:q, ]
dat2 <- dat[(q+1):n, ]

#Load the main function file
source("Functions_Graph_Horseshoe.R") # No graph related informative priors
source("Graph-Horseshoe-Informed-Priors.R")

# Define number of chains
n_chains <- 2
# Set different seeds for reproducibility
set.seed(629)  # Set global seed #123
seeds <- sample(1:1e6, n_chains)
# Run four chains
gibbs_samples_list <- lapply(1:n_chains, function(chain_id) {
  set.seed(seeds[chain_id])
  run_gibbs_sampler(X = as.matrix(X), y = y, n_iter = 20000, c = 100)
})


########################################## Run Parallel ########################
# parse the array ID as your chain index
# library(parallel)
# gibbs_samples_list <- mclapply(
#   1:n_chains,
#   function(chain_id) {
#     set.seed(seeds[chain_id])
#     run_gibbs_sampler(X, y, n_iter=20000, c=100)
#   },
#   mc.cores = n_chains   # uses 4 cores on one node
# )

##############################################################################

# Extract Samples
chain_mats <- lapply(gibbs_samples_list , function(x)  lapply(x$samples, function(x) x$beta))     # each is [iter × (p+1)]
all_samps  <- do.call(cbind, chain_mats)               # ((4*n_iter) × (p+1))


burn <- 10000

# 1) drop the first burn iterations in each chain
chain_post <- lapply(chain_mats, function(chain) chain[-seq_len(burn)])

# 2) how many post-burn draws per chain and total chains
n_post    <- length(chain_post[[1]])  # should be 2000
n_chains  <- length(chain_post)       # should be 4
p         <- nrow(chain_post[[1]][[1]])  # your dimension, e.g. 100


# 3) build an array [p × n_post × n_chains]
beta_array <- array(NA_real_, dim = c(p, n_post, n_chains))
for (i in seq_len(n_chains)) {
  # each chain_post[[i]] is a list of n_post matrices (p×1)
  beta_array[ , , i] <-
    do.call(cbind, lapply(chain_post[[i]], function(mat) as.numeric(mat)))
}

# 4) flatten across chains to get [p × (n_post * n_chains)]
beta_mat <- do.call(cbind, lapply(seq_len(n_chains),
                                  function(i) beta_array[ , , i]))


# ──────────────────────────────────────────────────────────────────────
# 4. Posterior mean on the *scaled* covariate scale
# ──────────────────────────────────────────────────────────────────────
beta_hat_scaled <- rowMeans(beta_mat)                # length-p vectorburn <- 1000
sgn       <- sign(beta_hat_scaled)

X_raw_test   <- as.matrix(dat2[ , -1])
X_mu_test    <- colMeans(X_raw_test)
X_sd_test    <- apply(X_raw_test, 2, sd)
X_std_test   <- scale(X_raw_test, center = X_mu_test, scale = X_sd_test)
X_test = X_std_test

# response: centre only
y_raw_test   <- dat2[ , 1]
y_mu_test    <- mean(y_raw_test)
y_ctr_test   <- y_raw_test - y_mu_test
y_test = y_ctr_test

# 2. back-transform the slopes: β_orig = β_std / sd(X)
beta_mat_raw <- sweep(beta_mat, 1, X_sd_test, "/")
beta_hat_raw  <- rowMeans(beta_mat_raw)
signal_pos<-c(18:23, 41:45)
beta_hat_raw[signal_pos]

library(BhGLM)
pred_HS = apply(beta_mat, 2, function(x) X_std_test %*% x)
pred_HS_mean <- rowMeans(pred_HS)
#pred_HS_mean_unctr <- pred_HS_mean + y_mu_test
#(a4_HS = measure.glm(y_raw_test, pred_HS_mean_unctr, family="gaussian"))
(a4_HS = measure.glm(y_ctr_test, pred_HS_mean, family="gaussian"))
RMSE = sqrt(a4_HS)
RMSE # Reported Root Mean squared error


mcmc_beta <- as.mcmc(t(beta_mat))
mcmc_beta <- as.mcmc(t(beta_mat_raw))
beta_quantile<-t(apply(mcmc_beta, 2, function(x) quantile(x, probs = c(0.025, 0.50, 0.975))))
beta_quantile[which(abs(beta_quantile[,2])>0.1),]

# Calculate raw beta estimates
beta_mcmc_raw_quantile <- t(apply(beta_mat_raw, 1, function(x) quantile(x, probs = c(0.025,0.05, 0.50, 0.95, 0.975))))
beta_mcmc_raw_quantile[signal_pos,]
beta <- sim200$beta_true
# lower‐ and upper‐CI bounds
ci_lower <- beta_mcmc_raw_quantile[, 1]
ci_upper <- beta_mcmc_raw_quantile[, 5]  # note: use the 97.5% quantile here

#ci_lower <- beta_quantile[, 1]
#ci_upper <- beta_quantile[, 3]  # note: use the 97.5% quantile here

#small safeguard for zero values
epsilon = 0.05


# prediction: CI excludes zero ⇒ we “detect” a non‐zero effect
predicted_nonzero <- (ci_lower >  epsilon) |
  (ci_upper < -epsilon)


#
# truth: which betas are truly non‐zero
actual_nonzero    <- sim200$beta_true != 0

# compute confusion‐matrix entries
TP <- sum(predicted_nonzero &  actual_nonzero)
FP <- sum(predicted_nonzero & !actual_nonzero)
FN <- sum(!predicted_nonzero &  actual_nonzero)
TN <- sum(!predicted_nonzero & !actual_nonzero)

# display as a matrix
conf_mat <- matrix(
  c(TP, FP,
    FN, TN),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    "Predicted" = c("NonZero","Zero"),
    "Actual"    = c("NonZero","Zero")
  )
)
print(conf_mat)

# Calculate kappas as a variable selection metric

iter2kappa <- function(it, n) {
  exp2psi <- exp(2 * it$psi)      # length-p
  1 / (1 + n / it$sigma2 * it$zeta2 * exp2psi)
}

kappa_chain <- lapply(gibbs_samples_list, function(ch) {
  post <- ch$samples[-seq_len(burn)]
  t(vapply(post, iter2kappa, numeric(p), n = n))    # matrix draws × p
})
kappa_all <- do.call(rbind, kappa_chain)            # stack all chains
kappa_hat <- apply(kappa_all, 2, median)

# Select variables for which kappa hat is less than 0.5
ind_selected <- which(kappa_hat<0.5)
print(ind_selected)
kappa_hat[which(kappa_hat<0.5)]

#Kappa based prediction
predicted_nonzero <- kappa_hat < 0.5
# truth: which betas are truly non‐zero
actual_nonzero    <- sim200$beta_true != 0

# compute confusion‐matrix entries
TP <- sum(predicted_nonzero &  actual_nonzero)
FP <- sum(predicted_nonzero & !actual_nonzero)
FN <- sum(!predicted_nonzero &  actual_nonzero)
TN <- sum(!predicted_nonzero & !actual_nonzero)

# display as a matrix
conf_mat <- matrix(
  c(TP, FP,
    FN, TN),
  nrow = 2, byrow = TRUE,
  dimnames = list(
    "Predicted" = c("NonZero","Zero"),
    "Actual"    = c("NonZero","Zero")
  )
)

print(conf_mat) # This is the confusion matrix based on the variable selection metric kappa

beta_hat_raw*(1-kappa_hat) # this is the amount of shrinkage achieved
beta_Adjusted<-beta_hat_raw*(1-kappa_hat)

# Compute accuracy of the estimated betas
rmse <- sqrt(mean((beta_mcmc_raw_quantile[,3] - sim200$beta_true)^2))
print(rmse)



