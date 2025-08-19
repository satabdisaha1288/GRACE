############################################################
##  Parameters  ------------------------------
burn        <- 10000         # iterations to discard in every chain
K           <- 50          # target number of edges in the final graph
eps_small   <- 1e-3        # “practically non-zero” threshold for PIP
pip_cut     <- 0.50        # median-probability rule
############################################################

library(Matrix)
library(igraph)
library(tidygraph)
library(ggraph)
library(dplyr)
library(mclust)    # for adjustedRandIndex
library(aricode)   # for NMI

## 1) Gather post-burn-in samples ------------------------------------------
all_samples <- unlist(
  lapply(gibbs_samples_list, \(ch) ch$samples[-seq_len(burn)]),
  recursive = FALSE
)
omega_list <- lapply(all_samples, `[[`, "omega")        # list of raw weights

## Concatenate all packed vectors into one long vector -------------------
w_all <- unlist(omega_list, use.names = FALSE)   # same as do.call(c, omega_list)


# build the row–col index matrix for the upper triangle:
edge_idx <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)

w_mat <-  do.call(cbind, omega_list)
wbar      <- rowMeans(w_mat)                              # posterior mean per edge
PIP       <- rowMeans(abs(w_mat) > eps_small)             # inclusion probability

## ---- 0.  edge index you already built ---------------------------------
## edge_idx   :  N_edge × 2 matrix of (row, col) pairs for upper-triangle
## PIP        :  length N_edge vector   =  Pr( |ω_{jk}| > eps_small | data )

cluster1<- 18:23
cluster2<- 40:45

## ---- 1.  logical mask for “block edges” -------------------------------
in_block <- (edge_idx[,1] %in% cluster1 & edge_idx[,2] %in% cluster1) |
  (edge_idx[,1] %in% cluster2 & edge_idx[,2] %in% cluster2)

## ---- 2.  pull out the PIP’s you care about ----------------------------
PIP_block <- PIP[ in_block ]      # within-cluster edges
PIP_cross <- PIP[!in_block ]      # all the rest (cross-block or noise)

## ---- 3.  quick summaries ----------------------------------------------
cat(sprintf("Block edges      : %d\n", length(PIP_block)))
cat(sprintf("  mean PIP       : %.3f\n", mean(PIP_block)))
cat(sprintf("  median PIP     : %.3f\n", median(PIP_block)))
cat(sprintf("Cross-block edges: %d\n", length(PIP_cross)))
cat(sprintf("  median PIP       : %.3f\n", median(PIP_cross)))

cut <- 0.50                         # pick your tolerance
keep_block  <- sum(PIP_block  > cut) # true positives
keep_cross  <- sum(PIP_cross  > cut) # false positives

cat(sprintf("TP kept  (block) : %d / 30\n", keep_block))
cat(sprintf("FP added (cross) : %d / 4920\n", keep_cross))

FN <- sum(PIP_block  <= cut)
# true negatives  = cross edges you correctly drop (<= cut)
TN <- sum(PIP_cross  <= cut)

## 3) Uncertainty-aware edge selection -------------------------------------
keep <- PIP >= pip_cut                      # median-probability model

# Keep top 50 edges only for sparsity, true graph has 30 edges
if (sum(keep) > K) {                        # trim to exactly K if needed
  cand        <- which(keep)
  top_idx     <- cand[order(PIP[cand], decreasing = TRUE)[1:K]]
  keep        <- logical(m);  keep[top_idx] <- TRUE
}

###################################### Build estimated graph $###############################


## Sparse adjacency matrix

#Estimated Graph
i <- edge_idx[keep, 1]
j <- edge_idx[keep, 2]
x <- wbar[keep]

W_hat <- sparseMatrix(i = i, j = j, x = x, dims = c(p, p))
W_hat <- W_hat + t(W_hat)                  # symmetrise once

colnames(W_hat) <- rownames(W_hat) <- as.character(seq_len(ncol(W_hat)))

# True graph
thr <- 0.5                     # threshold
W_true <- sim200$W_xy               # start with a copy
# 1. zero any entry ≤ thr
W_true[ W_true <= thr ] <- 0
# 2. remove self-loops
diag(W_true) <- 0


# Convert adjacency matrices to sparse matrices
W_hat<- as(W_hat, "dgCMatrix")
W_true <- as(W_true, "dgCMatrix")
#  Build igraph objects (unweighted or weighted) —
g1 <- graph_from_adjacency_matrix(W_true, mode="undirected", weighted=TRUE)
g2 <- graph_from_adjacency_matrix(as.matrix(W_hat), mode="undirected", weighted=TRUE)

# — 3. Compare with ARI and NMI —
ari <- adjustedRandIndex(comm1, comm2)
nmi <- NMI(comm1, comm2)
ari # Reported ARI value

########################################################################################


total_block  <- length(PIP_block)     # 30 in your simulation
total_cross  <- length(PIP_cross)     # 4 920

calc_counts <- function(cut) {
  # raw counts (still integers)
  TP <- sum(PIP_block  > cut)
  FN <- total_block  - TP
  FP <- sum(PIP_cross > cut)
  TN <- total_cross  - FP

  # coerce to double to avoid integer overflow
  TPd <- as.numeric(TP)
  FNd <- as.numeric(FN)
  FPd <- as.numeric(FP)
  TNd <- as.numeric(TN)

  # safe denominator in doubles
  denom <- sqrt(
    (TPd + FPd) *
      (TPd + FNd) *
      (TNd + FPd) *
      (TNd + FNd)
  )

  MCC <- ifelse(denom == 0, NA, (TPd * TNd - FPd * FNd) / denom)
  Precision = TPd / (TPd + FPd)
  Recall    = TPd / total_block
  F1 <- ifelse(Precision + Recall==0, 0, 2*Precision*Recall/(Precision+Recall))

  data.frame(
    cut       = cut,
    TP        = TP,
    FN        = FN,
    FP        = FP,
    TN        = TN,
    Precision = Precision,
    Recall    = Recall,
    MCC       = MCC,
    F1        = F1
  )
}

# then:
cuts <- seq(0.35, 0.75, by = 0.05)
df<-do.call(rbind, lapply(cuts, calc_counts))

FP = df$FP[4]
FN = df$FN[4]
TP = df$TP[4]
TN = df$TN[4]

# Report false positives and false negatives
print(FP)
print(FN)

#AUROC and AUPR calculations

# 1a. build labels and scores
labels <- c( rep(1, length(PIP_block)),    # positives = block edges
             rep(0, length(PIP_cross)) )   # negatives = cross-block edges
scores <- c( PIP_block, PIP_cross )

# 1b. compute ROC & AUC
library(pROC)
roc_obj <- roc(labels, scores, quiet = TRUE)
print(auc(roc_obj))
plot(roc_obj, main = sprintf("ROC curve (AUC = %.3f)", auc(roc_obj)))
aroc <- roc_obj$auc


library(PRROC)
pr <- pr.curve(
  scores.class0 = PIP_block,   # positives (class 0 in PRROC’s terminology)
  scores.class1 = PIP_cross,   # negatives (class 1)
  curve         = TRUE         # return full (recall, precision) path
)
# area under the PR curve
ap <- pr$auc.integral          # a.k.a. average precision
print(ap)

# to plot
plot(
  pr$curve[,1], pr$curve[,2], type = "l", lwd = 2,
  xlab = "Recall", ylab = "Precision",
  main = sprintf("PR curve  (AP = %.3f)", ap)
)

