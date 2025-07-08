# ======================================================================
#  TWO-GRAPH  simulation with INNER / OUTER exponential decay
#  – W_xy  drives ψ → λ → β → y               (X-y graph)
#  – W_x   drives Σ_X, Bray distances, etc.   (X-X graph)
# ======================================================================
library(MASS)     # mvrnorm
library(Matrix)   # nearPD
library(gtools)   # rdirichlet
set.seed(786)

# ---------- GLOBAL CONSTANTS ------------------------------------------
p_max      <- 600
signal_pos <- c(18:23, 40:45)        # functional clusters
func_blocks <- list(18:23, 40:45)    # for W_xy
pred_blocks <- list(1:10, 30:40, 60:70, 120:180, 240:300)  # for W_x
# -- decay parameters (your original values) ---------------------------
decay_in  <- 0.05     # strong ties  (slow decay)  inside clusters
decay_out <- 2        # weak  ties   (fast decay)  elsewhere
edge_scale <- 1

# ---------- HELPERS ----------------------------------------------------
row_normalise <- function(M) {
  rs <- rowSums(M);  M[rs > 0, ] <- sweep(M[rs > 0,,drop=FALSE], 1, rs[rs>0], "/");  M
}

############################### Load function hub_sprinkle ##########################

hub_sprinkle <- function(W,
                         hub_frac     = 0.10,
                         add_per_hub  = 25,
                         drop_prop    = 0.10,
                         jitter_sdlog = 0,
                         do_normalise = TRUE,
                         seed         = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # 1.  build graph --------------------------------------------------
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)

  # 2.  randomly drop edges -----------------------------------------
  if (drop_prop > 0 && igraph::ecount(g) > 0) {
    n_drop <- round(drop_prop * igraph::ecount(g))
    n_drop <- min(n_drop, igraph::ecount(g))         # safety
    if (n_drop > 0)
      g <- igraph::delete_edges(g, sample(igraph::E(g), n_drop))
  }

  # 3.  choose hubs --------------------------------------------------
  p    <- nrow(W)
  n_h  <- max(1, round(hub_frac * p))
  hubs <- sample(p, n_h)

  # 4.  add edges hub-by-hub ----------------------------------------
  old_w <- if (igraph::ecount(g) > 0) igraph::E(g)$weight else 1
  for (h in hubs) {
    ## recompute complement each time -------------------------------
    g_comp   <- igraph::complementer(g, loops = FALSE)
    cand     <- igraph::neighbors(g_comp, h, mode = "out") |> as.integer()
    if (length(cand) == 0) next

    k <- min(add_per_hub, length(cand))
    tgt <- sample(cand, k, replace = FALSE)
    new_w <- sample(old_w, k, replace = TRUE)
    g <- igraph::add_edges(g, c(rbind(rep(h, k), tgt)),
                           attr = list(weight = new_w))
    old_w <- c(old_w, new_w)                       # enlarge pool
  }

  # 5.  optional jitter ---------------------------------------------
  if (jitter_sdlog > 0 && igraph::ecount(g) > 0) {
    w  <- igraph::E(g)$weight
    nn <- length(w)
    idx <- sample(nn, round(0.5 * nn))
    w[idx] <- w[idx] * exp(stats::rnorm(length(idx), 0, jitter_sdlog))
    igraph::E(g)$weight <- w
  }

  # 6.  back to matrix + normalise ----------------------------------
  W2 <- igraph::as_adjacency_matrix(g, attr = "weight",
                                    sparse = FALSE) |> as.matrix()
  diag(W2) <- 0
  if (do_normalise)
    W2 <- row_normalise(W2)

  W2
}
############################################################################

# helper 1: functional (block-exp) kernel
make_Wxy <- function(p, blocks, decay_in = 0.05, decay_out = 2) {
  D  <- abs(outer(1:p, 1:p, "-"))
  W  <- exp(-decay_out * D)                 # fast-decay background
  for (B in blocks) {
    W[B, B] <- exp(-decay_in * D[B, B])     # slow-decay inside block
  }
  diag(W) <- 0
  W[W < 1e-6] <- 0
  W                                             # raw, symmetric
}

make_Wx_exp <- function(p, k = 3, w1 = 0.40, decay_rate = 0.5) {
  D <- abs(outer(1:p, 1:p, "-"))
  W <- w1 * decay_rate^(D - 1)      # exponential of distance
  W[D == 0 | D > k] <- 0            # zero diagonal and far lags
  W
}

# edge_key <- function(M) {
#   ## string “i-j” for every upper-triangular, non-zero edge
#   idx <- which(M > 0 & upper.tri(M), arr.ind = TRUE)
#   paste(idx[,1], idx[,2], sep = "-")
# }


# ---------- BUILD FULL-SIZE GRAPHS -------------------------------------
W_xy_full <- make_Wxy(p_max, func_blocks)                # strong kernel
W_x_orig  <- make_Wx_exp(p_max, k = 3, w1 = 0.45, decay_rate = 0.8)         # weak kernel


################################### Low overlap block #############################

# ----- remove all edges that coincide with W_xy -------------------
## 3. prune all true edges (drop TP)
true_mask  <- (W_xy_full > 0) & upper.tri(W_xy_full)
W_x_pruned <- W_x_orig
W_x_pruned[ true_mask | t(true_mask) ] <- 0

## 4. re-add a fraction of true edges (controlled TP)
overlap_frac <- 0.20                              # e.g. 20% TP
idx_true     <- which(true_mask, arr.ind=TRUE)
n_true       <- nrow(idx_true)
n_keep       <- round(overlap_frac * n_true)
keep_rows    <- sample(n_true, n_keep)
add_true_idx <- idx_true[keep_rows, , drop=FALSE]

add_mask_true <- matrix(FALSE, p_max, p_max)
add_mask_true[ add_true_idx ] <- TRUE
add_mask_true <- add_mask_true | t(add_mask_true)

W_x_step4 <- W_x_pruned
W_x_step4[ add_mask_true ] <- W_xy_full[ add_mask_true ]


## 5. add random noise edges (controlled FP)
noise_frac <- 0.10                               # e.g. add 10% noise
zero_mask  <- (W_xy_full == 0) & upper.tri(W_xy_full)
idx_zero   <- which(zero_mask, arr.ind=TRUE)
n_zero     <- nrow(idx_zero)
n_noise    <- round(noise_frac * n_zero)
noise_rows <- sample(n_zero, n_noise)
add_noise_idx <- idx_zero[ noise_rows, , drop=FALSE ]

add_mask_noise <- matrix(FALSE, p_max, p_max)
add_mask_noise[ add_noise_idx ] <- TRUE
add_mask_noise <- add_mask_noise | t(add_mask_noise)

# assign random small weights to those noise edges
# here we draw Uniform(0, max(W_x_orig)) but you can choose
W_x_full <- W_x_step4
max_w    <- max(W_x_orig)
W_x_full[ add_mask_noise ] <- runif(n_noise, 0, max_w)[ rep(1, each=1) ]
# (the rep here ensures correct length‐matching; you can also vectorize)

W_x_hub_sprinkle <- hub_sprinkle(
  W_xy_full,
  hub_frac     = 0.40,     # more hubs
  add_per_hub  = 40,       # many more noise edges
  drop_prop    = 0.50,     # drop half of the true edges
  jitter_sdlog = 1.2,
  do_normalise = TRUE,
  seed         = 123
)

# ---------- OPTIONAL: inject partial overlap ---------------------------
# overlap_frac <- 0.10      # 0 = disjoint, 1 = identical
# if (overlap_frac > 0) {
#   mask <- (runif(p_max^2) < overlap_frac) &
#     (row(W_x_full) < col(W_x_full)) & (W_xy_full > 0)
#   mask <- mask | t(mask)                  # symmetrise
#   W_x_full[mask] <- W_xy_full[mask]       # copy edge weights
#   #W_x_full <- row_normalise(W_x_full)
# }

# ---------- DRAW ψ ONCE (depends only on W_xy_full) --------------------
phi <- 0.6;  tau2 <- 0.30;  eps <- 1e-6
Q_xy <- diag(rowSums(W_xy_full)) - phi * W_xy_full + eps * diag(p_max)
psi_full <- as.vector(mvrnorm(1, rep(0, p_max), tau2 * solve(Q_xy)))
lambda_full <- exp(psi_full)

# ---------- FIXED MAGNITUDES FOR 12 SIGNAL βs --------------------------
b1 <-  0.5 * c( 2.18, 1.99, 2.12, 1.81, 1.86, 2.34)   # positive cluster
b2 <- -0.5 * c( 2.41, 1.65, 2.51, 1.95, 1.93, 1.85)   # negative cluster

# ---------- MAIN GENERATOR ---------------------------------------------
simulate_once <- function(n = 800, p = 500,
                          sigma_y = 2,
                          seed = NULL,
                          overlap = c("high", "low", "hub-sprinkle")) {

  overlap <- match.arg(overlap)

  if (!is.null(seed)) set.seed(seed)
  if (p > p_max)
    stop(sprintf("p cannot exceed p_max = %d", p_max))

  ## 1. β (fixed pattern, centred, rounded) ----------------------------
  beta      <- numeric(p)
  beta_idx  <- signal_pos[signal_pos <= p]

  beta[beta_idx[1:6]]  <- b1 * lambda_full[beta_idx[1:6]]
  beta[beta_idx[7:12]] <- b2 * lambda_full[beta_idx[7:12]]

  beta[beta_idx] <- beta[beta_idx] - mean(beta[beta_idx])
  beta[beta_idx] <- round(beta[beta_idx], 2)         # round *after* centring
  beta[tail(beta_idx, 1)] <- beta[tail(beta_idx, 1)] -
    sum(beta[beta_idx])     # enforce exact zero-sum

  ## 2. Graphs ---------------------------------------------------------
  W_xy <- W_xy_full[1:p, 1:p]
  W_x <- if (overlap == "high") {
    W_x_orig[1:p, 1:p]
  } else if (overlap == "low") {
    W_x_full[1:p, 1:p]
  } else {
    W_x_hub_sprinkle[1:p, 1:p]
  }

  ## 3. Σ_X = AR(1) + predictor boost ---------------------------------
  rho_ar <- 0.2
  Dp     <- abs(outer(1:p, 1:p, "-"))
  boost  <- 0.15

  Sigma <- rho_ar ^ Dp + boost * (W_x > 0)
  diag(Sigma) <- 1
  Sigma <- as.matrix(nearPD(Sigma, corr = TRUE)$mat)

  ## 4. Generate compositional X --------------------------------------
  X_raw <- mvrnorm(n, rep(0, p), Sigma)
  X_pos <- exp(sweep(X_raw, 2, colMeans(X_raw), "+"))
  X_cmp <- X_pos / rowSums(X_pos)
  X_clr <- log(X_cmp)

  ## 5. Response -------------------------------------------------------
  y <- drop(X_clr %*% beta) + rnorm(n, 0, sigma_y)

  ## 6. Return ---------------------------------------------------------
  list(
    n          = n,
    p          = p,
    y          = y,
    X_rel      = X_cmp,
    X_clr      = X_clr,
    beta_true  = beta,
    sigma_y    = sigma_y,
    rho_ar     = rho_ar,
    boost      = boost,
    phi        = phi,
    tau2       = tau2,
    psi_true   = psi_full[1:p],
    W_xy       = W_xy,
    W_x        = W_x
  )
}


# ---------- EXAMPLE RUN -----------------------------------------------
sim200 <- simulate_once(n = 400, p = 100, seed = 1, overlap = "low")
