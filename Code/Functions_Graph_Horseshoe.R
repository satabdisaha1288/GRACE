

# Example usage:
library(MASS)
library(BayesLogit)
library(GIGrvg)
library(MCMCpack)
library(SeBR)
# Example usage:
#
# X <- as.matrix(dat1[,-1])
# #X <- scale(X)
# p <- ncol(X)
# n <- nrow(X)
# y <- dat1[,1]
# #y<- scale(y)
# #c = 100*sqrt(p)
# c = 100
# #c = 300
# XtX  <- crossprod(X)
# Xt_y <- crossprod(X, y)
# #W3 <- rewire_W(W_true, frac = 0.3)
#W_obs<- W_true
# W_pos<- pmax(W_obs, 1e-6)

################################
#For simulated data
# predictors: centre & scale
X_raw   <- as.matrix(dat1[ , -1])
X_mu    <- colMeans(X_raw)
X_sd    <- apply(X_raw, 2, sd)
X_std   <- scale(X_raw, center = X_mu, scale = X_sd)
X = X_std

# response: centre only
y_raw   <- dat1[ , 1]
y_mu    <- mean(y_raw)
y_ctr   <- y_raw - y_mu
y = y_ctr

p <- ncol(X)
n <- nrow(X)
c = 100

XtX  <- crossprod(X)
Xt_y <- crossprod(X, y)
W_obs<- sim200$W_x
W_pos<- pmax(W_obs, 1e-6)
###########################################################

as_sq <- function(A, p) if (is.matrix(A)) A else matrix(A, p, p)

safe_draw_beta <- function(Vinv, Xt_y_over_sigma2, min_ridge = 1e-6) {
  p   <- nrow(Vinv)
  rid <- min_ridge
  repeat {
    ok <- TRUE
    R  <- tryCatch(chol(Vinv + diag(rid, p)), error = function(e) { ok <<- FALSE })
    if (ok) break
    rid <- rid * 10                      # inflate ridge and retry
  }

  ## (a) posterior mean  μ = Vinv⁻¹ Xᵀy / σ²  using R
  mu <- backsolve(R, backsolve(R, Xt_y_over_sigma2, transpose = TRUE))

  ## (b) draw β = μ + R⁻¹ z     (z ∼ N(0,Iₚ))
  z  <- rnorm(p)
  beta <- mu + backsolve(R, z, transpose = TRUE)
  beta
}


# -----------------------------------------------------------
#  chol_pd()  –  safe Cholesky factor
# -----------------------------------------------------------
# M            : a symmetric numeric matrix
# init_ridge   : starting value for the diagonal ridge
#
# value        : the upper-triangular Cholesky factor R
#                such that  R' R = M + δ I   with δ ≥ init_ridge
# -----------------------------------------------------------
chol_pd <- function(M, init_ridge = 1e-6) {
  p   <- nrow(M)
  rid <- init_ridge
  repeat {
    ## try the Cholesky; if it fails, inflate the ridge and retry
    R <- try(chol(M + diag(rid, p)), silent = TRUE)
    if (!inherits(R, "try-error"))
      return(R)          # success: return the factor
    rid <- rid * 10      # otherwise add 10× more ridge and loop
  }
}

# Initialization
initialize_sampler <- function(p, n, X, y) {
  beta <- rnorm(p)
  eta <- rnorm(p)
  psi <- rnorm(p)
  sigma2 <- 29
  a_sigma <- 2.2 
  b_sigma <- 2.9 #2.55
  m0 = 0.10*p
  par_ratio <- (m0/(p - m0))*( 1/sqrt(n))    
  #zeta2 <- 0.1 #Simulation
  #zeta2 <- 1 #10 # 1 #0.1
  nu <- 0.1
  tau2 <- rep(1, p) #Main one
  omega <- W_pos #+ 0.01                # as you had
  #omega <- W_true
  omega[omega == 0] <- 1e-3             # add tiny mass to empty edges
  #omega_raw <- omega
  ########################Block for real data#####################
  #zeta2 <- 0.5   #Main one                   # start larger
  zeta2 <- 1
  nu_shape <- 0.5
  nu_rate  <-  0.5#0.04 #5      #Main one             # mean E[zeta2] ≈ 0.1
  #nu_rate  <- 1                   # E[ζ²] now ≳1
  a_tau <- 2.3 
  b_tau <- 1                     # heavier tails for each tau_j^2
  c <- 100                          # weaker zero‑sum penalty
  ########################Block for real data#####################

  ## ------------------  Hyper‑priors  (high trust) -------------------------
  # mask_edge <- W_obs > 0                      # keep zero elsewhere
  # alpha_mat            <- matrix(4, p, p)     # stronger shape
  # q_mat                <- matrix(5, p, p)
  # alpha_mat[!mask_edge] <- 4                  # safe, but irrelevant
  # q_mat    [!mask_edge] <- 5
  #
  # kappa_mat <- ((q_mat - 1) / alpha_mat) * W_obs
  # kappa_mat <- pmax(kappa_mat, 1e-6)          # numerical guard
  # diag(alpha_mat) <- diag(q_mat) <- diag(kappa_mat) <- 0

  ############### 1.  Base (conservative) prior for every edge ######################
  # #Widen the Gamma tails
  # alpha_mat <- matrix(0.7, p, p); diag(alpha_mat) <- 0 #Main one
  # #alpha_mat <- matrix(3, p, p); diag(alpha_mat) <- 0 #Main one (for stronger confidence on the true graph)
  # #alpha_mat <-  hyper$alpha ; diag(alpha_mat) <- 0 #Main one
  # q_mat     <- matrix(2.5, p, p); diag(q_mat) <- 2.1 #Main one
  # #q_mat     <- matrix(1.5, p, p); diag(q_mat) <- 1.1 #Main one
  # #q_mat     <- matrix(1, p, p); diag(q_mat) <- 1 #Main one
  # #q_mat     <- hyper$q
  # #kappa_mat <- matrix(400,p,p); diag(kappa_mat)<- 1
  #
  # #W_obs_pos  <- pmax(W_obs, 0)          # or abs(W_obs)
  # kappa0    <- 0.05
  # kappa_mat  <- kappa0 + ((q_mat - 1) / alpha_mat) * W_obs # Main one
  # #kappa_mat  <- ((q_mat - 1) / alpha_mat)  # Main one
  # diag(kappa_mat) <- 1
  ##################################################################################
  ########################.   Weakly informative prior ###########################

  alpha0 <- 1                       # shape  – weak tails
  q0     <- 2.1                       # q > 2 so variance exists
  kappa0 <- W_obs                # centre on observed |weight|

  alpha_mat <- matrix(alpha0, p, p);  diag(alpha_mat) <- 0
  q_mat     <- matrix(q0   , p, p);   diag(q_mat)     <- q0
  kappa_mat <- ((q_mat - 1) / alpha_mat) * kappa0
  diag(kappa_mat) <- 1

  ################## High confidence Edges #########################
  ## ---------------------------------------------------------------------
  ## 0.  Identify the edges you “trust”
  ## ---------------------------------------------------------------------
  # in_block <- (outer(seq_len(p), seq_len(p), Vectorize(function(i, j)
  #   (i %in% cluster1 && j %in% cluster1) ||
  #     (i %in% cluster2 && j %in% cluster2))))   # p × p logical
  # in_block[ lower.tri(in_block, diag = TRUE) ] <- FALSE   # keep upper-tri only
  #
  # cross_block <- !in_block & upper.tri(in_block)          # complement

  ## ---------------------------------------------------------------------
  ## ---------------------------------------------------------------------
  ## 2.  Boost the slab *inside* trusted blocks  -------------------------
  ##    (pick ONE of the patterns below — or combine modestly)
  ## ---------------------------------------------------------------------

  ## ------------------------------------------------------------------
  ## 2.  Boost *within-block* ─ higher α, lower κ
  ## ------------------------------------------------------------------
  # alpha_mat[in_block]  <- 8                   # ↑ shape  ⇒ stronger Bayes factor
  # kappa_mat[in_block] <- 0.1 * kappa_mat[in_block]  # ↓ κ ⇒ larger mean
  #
  # ## 3.  (Optional) Harden the *cross* edges
  # alpha_mat[cross_block] <- 1.5               # smaller shape
  # kappa_mat[cross_block] <- 2 * kappa_mat[cross_block]  # higher κ ⇒ lower mean
  # diag(alpha_mat) <- 0
  ############################################################

  #Fix nu to obtain horseshoe like shrinkage on the global scale zeta
  #Simulation
  #nu_shape <- 0.5
  #nu_rate  <-  50#10
  #For cooedinate wise
  #step_var <- rep(0.5, p) # MALA step-size variance for each ψ component (tune these values as needed)
  #For vectorized
  step_var <- 1e-3
  #kappa_mat <- ((q_mat - 1) / alpha_mat) * omega
  rho   <- matrix(1.0, p, p)
  rho_curr <- rho
  omega_raw = W_pos
  c <- 100 # Define the constant c
  T <- rbind(diag(p), c * matrix(1, nrow = 1, ncol = p))
  list(beta = beta, eta = eta, psi = psi, sigma2 = sigma2, a_sigma = a_sigma, b_sigma = b_sigma,
       zeta2 = zeta2, nu = nu, tau2 = tau2, a_tau = a_tau, b_tau = b_tau, omega = omega,
       alpha_mat = alpha_mat, kappa_mat  = kappa_mat ,q_mat = q_mat, T = T,
       step_var = step_var, rho = rho, nu_shape = nu_shape,
       nu_rate  = nu_rate, omega_raw = omega_raw, par_ratio = par_ratio)
}

pack_omega   <- \(W)  as.numeric(W[upper.tri(W, diag = FALSE)])
unpack_omega <- \(v,p){
  M <- matrix(0, p, p)
  M[upper.tri(M, diag = FALSE)] <- v
  M + t(M)
}


make_U_ops <- function(p, c) {
  A <- diag(p) + c^2 * tcrossprod(rep(1, p))   # <- keep only this
  U        <- chol(A)
  U_inv    <- backsolve(U, diag(p))
  list(apply_U    = function(b)  U %*% b,
       apply_Uinv = function(b)  U_inv %*% b,
       U = U, U_inv = U_inv,
       v = rep(1, p),            # <- un-normalised one-vector
       c2 = c^2)
}



Uop <- make_U_ops(p, c)        # c = 100 still fine



safe_mvn_prec <- function(Vinv, b, min_ridge = 1e-6) {
  R  <- chol_pd(Vinv, min_ridge)
  mu <- backsolve(R, backsolve(R, b, transpose = TRUE))
  mu + backsolve(R, rnorm(length(b)), transpose = TRUE)
}

## ---- fast β-draw with cached Cholesky --------------------------
make_beta_drawer_p <- function(XtX, Xt_y, p) {
  R_prev <- NULL          # cached factor
  Vinv_prev <- NULL
  function(sigma2, zeta2, psi, Uop) {

    p     <- length(psi)                # length of the parameter vector
    #lambda_vec <- as.numeric(exp(-psi)) #/ (sigma2 * zeta2))
    lambda_vec <- as.numeric( exp(-2 * psi) )     # plain numeric, length p
    λinv       <- diag(lambda_vec,  p,  p)        # p × p diagonal *matrix*
    vcol       <- matrix(Uop$v, ncol = 1)         # p × 1 column *matrix*
    Av <- λinv %*% vcol                           # now conformable
    # p × 1
    prior <- λinv +
      (Uop$c2 / (1 + Uop$c2 * sum(vcol * Av))) * tcrossprod(Av)
    Vinv  <- XtX/sigma2 + prior/(sigma2*zeta2)



    if (is.null(R_prev) || !identical(dim(Vinv), dim(Vinv_prev)) ||
        max(abs(Vinv - Vinv_prev)) > 1e-12) {
      R_prev   <<- chol_pd(Vinv)       # expensive, done rarely
      Vinv_prev<<- Vinv
    }
    safe_mvn_prec(Vinv_prev, Xt_y/sigma2, min_ridge = 1e-8)
  }
}



## 0.  Helper: numerically safe Cholesky with automatic jitter
chol_spd <- function(A, jitter = 1e-8, max_tries = 6) {
  for (k in 0:max_tries) {
    R <- try(chol(A + diag(jitter * 10^k, nrow(A))), silent = TRUE)
    if (!inherits(R, "try-error")) return(R)
  }
  stop("chol_spd(): matrix not positive‑definite even after jitter")
}

## 1.  Helper: guard against Inf / NA precisions
safe_scales <- function(sigma2, zeta2, psi,
                        floor = 1e-12, ceil = 1e+12) {
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- floor
  if (!is.finite(zeta2)  || zeta2  <= 0) zeta2  <- floor
  psi[!is.finite(psi)] <- 0
  lambda_inv <- pmin(pmax(exp(-2*psi) / (sigma2 * zeta2), floor), ceil)
  lambda_inv
}

## 4.  n × n Woodbury drawer  (for  p > n)

eta_draw <- make_beta_drawer_p(XtX, Xt_y, p)
make_beta_drawer_n <- function(eta_draw, Uop){
  function(sigma2, zeta2, psi){
    eta <- eta_draw(sigma2, zeta2, psi, Uop)      # draw in η‑space
    backsolve(Uop$U, eta, transpose = FALSE)      # β = U⁻¹ η
  }
}



make_beta_drawer_woodbury <- function(X, y) {
  n <- nrow(X);  p <- ncol(X)

  function(sigma2, zeta2, psi, Uop = NULL) {
    ## ---- 1. Pre‑compute pieces -------------------------------------------
    Xs <- X / sqrt(sigma2)        # n×p   (standardised design  √σ² in denom)
    ys <- y / sqrt(sigma2)        # n‑vec (standardised response)

    lambda  <- exp(-2*psi)          # p‑vec  λ_j
    #D_inv   <- (1 / zeta2) / lambda       # diag entries of  D^{-1}
    D_diag     <- zeta2 / lambda
    #D_diag  <- 1 / D_inv                   # (= zeta2 * λ)  diag entries of D

    ## ---- 2. Build the small n×n matrix  M = I + X D X' -------------------
    XD      <- sweep(Xs, 2, D_diag, "*")   # n×p  =  X D
    M       <- diag(n) + XD %*% t(Xs)      # n×n

    ## ---- 3. Posterior mean  μ -------------------------------------------
    b_star  <- crossprod(Xs, ys)           # p‑vec  X'y / σ²
    w_bar   <- solve(M, XD %*% b_star)     # n‑vec
    mu      <- D_diag * (b_star - crossprod(Xs, w_bar))

    ## ---- 4. Draw η  ~  N(0, Q^{-1})  -------------------------------------
    z1  <- rnorm(p, sd = sqrt(D_diag))     # p‑vec  ~ N(0,D)
    rhs <- Xs %*% z1 + rnorm(n)            # n‑vec  (Rue 2001)
    w2  <- solve(M, rhs)                   # n‑vec
    eta <- D_diag * (z1 - crossprod(Xs, w2))

    ## ---- 5. One posterior draw -------------------------------------------
    mu + eta
  }
}




## 5.  Initialise the correct β‑drawer once, before the Gibbs loop
# init_beta_drawer <- function(X, y) {
#   n <- nrow(X); p <- ncol(X)
#   if (p <= n) {
#     XtX  <- crossprod(X)
#     Xt_y <- crossprod(X, y)
#     make_beta_drawer_p(XtX, Xt_y, p)
#   } else {
#     make_beta_drawer_n(X, y)
#   }
# }


#beta_draw <- make_beta_drawer_n(X,y)


###########################
# 1. Define the Functions #
###########################

# Gradient of the log-posterior for ψ
# -------------------------------------------------------------
# Gradient of  log p(ψ | rest)
# -------------------------------------------------------------
# Arguments
#   psi     : numeric length‑p vector  ψ
#   omega   : p × p weight matrix  ω_{jk}
#   beta    : numeric length‑p vector  β
#   sigma2  : noise variance  σ²
#   zeta2   : global scale     ζ²
#   Uop     : list returned by make_U_ops(), must contain  apply_U()
#   tau2    : either a scalar  τ²  or length‑p vector  τ_j²
#
# Value     : numeric length‑p gradient  ∂/∂ψ_j  log p
# -------------------------------------------------------------
gradient_log_posterior <- function(psi, omega, beta,
                                   sigma2, zeta2, Uop, tau2)
{
  p      <- length(psi)

  ## 1.  Compute  U β  once (vector length p)
  Ubeta  <- Uop$apply_U(beta)

  ## 2.  Likelihood–part term  −½ (Uβ)_j² exp(−2ψ_j)/(σ² ζ²)
  term1  <- (Ubeta^2) * exp(-2 * psi) / (sigma2 * zeta2)   # length p

  ## 3.  Prior CAR term  −½/τ_j² ∑_k ω_{jk}(ψ_j−ψ_k)²
  diff_mat   <- sweep(matrix(psi, p, p), 2, psi, "-")      # ψ_j − ψ_k
  term2_vec  <- (omega * diff_mat) / tau2                  # divides each row by τ_j²
  term2      <- rowSums(term2_vec)                         # length p

  ## 4.  Linear term  −∑_j ψ_j   →   derivative  −1
  ## ------------------------------------------------
  #- term1 - term2 - 1
  + term1 - term2 - 1
  #- term1 - term2                 # drop the "−1"

}


# Log-posterior function for ψ.
# Here we assume the log posterior is (up to additive constants)
#
#    log f(ψ) ∝ -0.5 * sum_j (Uβ)_j^2 exp(-2 ψ[j])/(sigma2*zeta2)
#             -0.5 * sum_j [1/tau2[j]] * sum_{k in Neighbors(j)} ω[j,k]*(ψ[j]-ψ[k])^2
#             - sum_j ψ[j]
#
# -------------------------------------------------------------
# log_posterior(ψ)   —  up to an additive constant
# -------------------------------------------------------------
#   psi     : numeric length‑p vector   ψ
#   omega   : p × p weight matrix       ω_{jk}
#   beta    : numeric length‑p vector   β
#   sigma2  : noise variance            σ²
#   zeta2   : global scale              ζ²
#   Uop     : list from make_U_ops(), must contain apply_U()
#   tau2    : either a scalar τ² or a length‑p vector τ_j²
# -------------------------------------------------------------
log_posterior <- function(psi, omega, beta,
                          sigma2, zeta2, Uop, tau2)
{
  p <- length(psi)

  ## 1.  Likelihood part   −½ Σ_j (Uβ)_j² exp(−2ψ_j) / (σ² ζ²)
  Ubeta    <- Uop$apply_U(beta)                           # length p
  log_like <- -0.5 * sum((Ubeta^2) * exp(-2 * psi)) / (sigma2 * zeta2)

  ## 2.  CAR prior part    −½ Σ_{j<k} ω_{jk}(ψ_j−ψ_k)² / τ_j²
  diff_mat <- sweep(matrix(psi, p, p), 2, psi, "-")       # ψ_j − ψ_k
  if (length(tau2) == 1L)          # scalar τ²
    tau_mat <- tau2
  else                              # vector τ_j² : replicate row‑wise
    tau_mat <- matrix(tau2, p, p)

  log_prior <- -0.5 * sum(omega * diff_mat^2 / tau_mat)

  ## 3.  Linear term       − Σ_j ψ_j
  log_linear <- - sum(psi)

  log_like + log_prior + log_linear
  # log_like + log_prior
}

safe_dnorm <- function(x, mean, sd) {
  if (!is.finite(mean) || !is.finite(sd) || sd <= 0)
    return(-Inf)
  val <- dnorm(x, mean, sd, log = TRUE)
  if (!is.finite(val)) -Inf else val                    # <- never NaN
}



#################################################################
##  metric–pre-conditioned vector-MALA                         ##
#################################################################
## Arguments as before:
##   psi      – current p-vector
##   omega    – p × p neighbour matrix
##   beta     – current β
##   sigma2   – noise variance
##   zeta2    – global scale
##   Uop      – list with Uop$apply_U
##   tau2     – scalar or length-p vector τ²
##   step_var – scalar or length-p      ε²  ( proposal variance )
##
## Returns
##   list(psi = new ψ, accepted = TRUE / FALSE)
#################################################################

# ---- (1)  diagonal metric  -----------------------------------
make_metric <- function(psi, omega, beta, sigma2, zeta2, Uop, tau2) {
  ## Fisher-info–like diagonal:   (Uβ)_j² e^{−2ψ_j}/(σ² ζ²)
  Ub  <- Uop$apply_U(beta)
  #term1 <- (Ub^2) * exp(-2 * psi) / (sigma2 * zeta2)
  term1 <- 2*(Ub^2) * exp(-2 * psi) / (sigma2 * zeta2)
  ## CAR prior diagonal: Σ_k ω_{jk}/τ²
  omega_row <- rowSums(omega)                 # length p
  term2 <- omega_row / tau2                  #  τ² is scalar or vector

  #term1 + term2 + 1                          # +1 from linear term
  term1 + term2
}

# ======================================================
#  Helper:  quadratic form  Q_β  =  β' U Λ^{-1} U' β
#           uses the operator  Uop$apply_U()
# ======================================================
quad_form <- function(beta, Uop, psi) {
  Ubeta  <- Uop$apply_U(beta)            # U β   (length p)
  lambda <- exp(psi)
  sum((Ubeta^2) / lambda^2)             # Σ (Uβ)_j² / λ_j²
}

# ======================================================
#  σ² | y, β           —  Inverse‑Gamma
# ======================================================
sample_sigma2 <- function(y, X, beta,
                          a_sigma, b_sigma,
                          Uop, psi, zeta2) {

  RSS     <- sum((y - X %*% beta)^2)
  Q_beta  <- quad_form(beta, Uop, psi)

  shape   <- a_sigma + (length(y) + ncol(X)) / 2
  rate    <- b_sigma + 0.5 * (RSS + Q_beta / zeta2)   # divide by ζ²
  1 / rgamma(1, shape = shape, rate = rate)
}

# ======================================================
#  ζ² | β, σ²          —  Inverse‑Gamma
# ======================================================
sample_zeta2 <- function(beta, sigma2,
                         Uop, psi, nu) {
  #maybe fix nu = 1
  Q_beta <- quad_form(beta, Uop, psi)

  shape  <- 0.5 + length(beta) / 2
  #shape  <- 0.5 + m0 / 2
  rate   <- Q_beta / (2 * sigma2) + 1 / nu
  1 / rgamma(1, shape = shape, rate = rate)
}

# ======================================================
#  ν   | ζ²            —  Inverse‑Gamma 
# ======================================================
sample_nu <- function(zeta2, par_ratio, sigma2) {
  shape <- 0.5
  rate  <- (1/ (par_ratio^2 * sigma2)) + 1 / zeta2
  1 / rgamma(1, shape = shape, rate = rate)
}

# ======================================================
#  τ²  | ψ, ω          —  Inverse‑Gamma  (global τ²)
# ======================================================

sample_tau2 <- function(a_tau, b_tau, omega, psi) {
  p     <- length(psi)

  # total unique‐edge weight
  sum_w <- sum(omega[upper.tri(omega)])       # now ≈ number of non-isolated nodes / 2

  # quadratic form ψᵀQψ via pairwise diffs
  diff2   <- outer(psi, psi, "-")^2
  sum_q   <- sum(omega[upper.tri(omega)] * diff2[upper.tri(diff2)])

  # shape grows like p/2, not edges/2
  shape <- a_tau + 0.5 * sum_w

  # rate is the usual ½ ψᵀQψ + b_tau
  rate  <- b_tau + 0.5 * sum_q

  1 / rgamma(1, shape, rate)
}


# ------------------------------------------------------
#  Example Gibbs Sampler for omega_{jk} and rho_{jk}
# ------------------------------------------------------
## helper: safe Gamma draw  (vectorised)
safe_rgamma <- function(idx, shape_vec) {
  if (!length(idx)) return(invisible())
  ok      <- is.finite(shape_vec) & (shape_vec > 0)
  n_ok    <- sum(ok)
  if (n_ok)
    rho_curr[idx][ok] <<- rgamma(n_ok, shape = shape_vec[ok], rate = 1)
  # all bad shapes are left unchanged
}



# Gibbs update for omega_{jk} and rho_{jk}
# --------------------------------------------------------------
# alpha_mat : p × p matrix of α_{jk} (shape for ω)
# kappa_mat : p × p matrix of κ_{jk} (coupling constant)
# q_mat     : p × p matrix of q_{jk} (shape parameter for ρ)
# psi       : numeric length-p vector ψ
# tau2      : scalar τ² or length-p vector τ_j²
# omega_curr: current p × p symmetric matrix ω
# rho_curr  : current p × p symmetric matrix ρ
# Returns   : list(omega = updated ω, rho = updated ρ)
# --------------------------------------------------------------
gibbs_sampler_omega_rho <- function(alpha_mat, kappa_mat, q_mat,
                                    psi, tau2,
                                    omega_curr, rho_curr) {
  # Flatten psi to vector
  psi_vec <- as.numeric(psi)
  p       <- length(psi_vec)

  # Build full τ matrix
  if (length(tau2) == 1L) {
    tau_mat <- matrix(tau2, p, p)
  } else if (length(tau2) == p) {
    tau_mat <- matrix(tau2, nrow = p, ncol = p, byrow = TRUE)
  } else {
    stop("tau2 must be scalar or length-p vector")
  }
  tau_mat <- pmax(tau_mat, 1e-8)

  # Squared differences (ψ_j - ψ_k)^2
  Psi_mat   <- matrix(psi_vec, p, p)
  psi_diff2 <- (Psi_mat - t(Psi_mat))^2

  # Lambda and Psi matrices for rho
  lambda_mat <- q_mat - alpha_mat
  psi_mat    <- 2 * kappa_mat * omega_curr

  # Off-diagonal mask
  mask_off <- upper.tri(alpha_mat)

  # Rho update masks
  idx_gam  <- mask_off & (psi_mat == 0) & (lambda_mat > 0)
  idx_igam <- mask_off & (psi_mat == 0) & (lambda_mat < 0)
  idx_gig  <- mask_off & (psi_mat > 0) & is.finite(psi_mat) &
    is.finite(lambda_mat) & (lambda_mat != 0)

  # ρ | ... updates
  if (any(idx_gam)) {
    rho_curr[idx_gam] <- rgamma(sum(idx_gam),
                                shape = pmax(lambda_mat[idx_gam], 1e-8),
                                rate  = 1)
  }
  if (any(idx_igam)) {
    rho_curr[idx_igam] <- 1 / rgamma(sum(idx_igam),
                                     shape = pmax(-lambda_mat[idx_igam], 1e-8),
                                     rate  = 1)
  }
  if (any(idx_gig)) {
    rho_curr[idx_gig] <- rgig(sum(idx_gig),
                              lambda = lambda_mat[idx_gig],
                              chi    = rep(2, sum(idx_gig)),
                              psi    = psi_mat[idx_gig])
  }
  # Mirror rho and zero diagonal
  rho_curr[lower.tri(rho_curr)] <- t(rho_curr)[lower.tri(rho_curr)]
  diag(rho_curr) <- 0

  # Omega update: only upper triangle
  rate_mat <- kappa_mat / rho_curr + psi_diff2 / (2 * tau_mat)
  rate_mat <- pmax(rate_mat, 1e-6)

  idx_up <- upper.tri(alpha_mat)
  if (any(idx_up)) {
    omega_curr[idx_up] <- rgamma(sum(idx_up),
                                 shape = pmax(alpha_mat[idx_up], 1e-3),
                                 rate  = rate_mat[idx_up])
  }
  # Mirror omega and zero diagonal
  omega_curr[lower.tri(omega_curr)] <- t(omega_curr)[lower.tri(omega_curr)]
  diag(omega_curr) <- 0

  list(omega = omega_curr, rho = rho_curr)
}

#beta_draw <- make_beta_drawer(XtX, Xt_y, ncol(X))

run_gibbs_sampler <- function(X, y, n_iter = 10000, c = 100) {

  X        <- as.matrix(X)
  p        <- ncol(X)
  n        <- nrow(X)

  ## -- 0.  initial state  -----------------------------------------
  params   <- initialize_sampler(p, n, X, y)
  # --- once, after initialize_sampler() ---------------------------

  ## -- 1.  build operator for U and U^{-1} ------------------------
  Uop      <- make_U_ops(p, c)          # returns apply_U / apply_Uinv

  samples  <- vector("list", n_iter)
  acc_counts   <- numeric(p)   # accepted moves per coord
  prop_counts  <- numeric(p)   # proposed moves per coord

  acc_total  <- 0L    # how many sweeps accepted
  prop_total <- 0L    # how many sweeps proposed
  beta_draw <- make_beta_drawer_n(eta_draw, Uop)
  ## -- 2.  Gibbs loop  -------------------------------------------
  for (iter in seq_len(n_iter)) {

    ## ---- β | rest  ---------------------------------------------
    if (iter %% 100 == 0) {
      cat("Iteration =", iter, "\n")
    }

    #beta_draw <- make_beta_drawer_woodbury(X,y)
    params$beta <- beta_draw(params$sigma2, params$zeta2,
                             params$psi)
    params$beta  <- params$beta  - mean(params$beta )   # after each update


    ## ---- ψ | rest   (one full MALA sweep) ----------------------
    # res <- mala_sample_psi(
    #   psi      = params$psi,
    #   omega    = params$omega,
    #   beta     = params$beta,
    #   sigma2   = params$sigma2,
    #   zeta2    = params$zeta2,
    #   Uop      = Uop,
    #   tau2     = params$tau2,
    #   step_var = params$step_var)
    #
    # params$psi <- res$psi


    # params$psi <- update_psi_slice(psi      = params$psi,
    #                         omega    = params$omega,
    #                         beta     = params$beta,
    #                         sigma2   = params$sigma2,
    #                         zeta2    = params$zeta2,
    #                         Uop      = Uop,
    #                         tau2     = params$tau2)
    # params$psi <- update_psi_elliptical (psi = params$psi,
    #                        omega =params$omega,
    #                        beta = params$beta,
    #                        sigma = params$sigma2,
    #                        zeta2 = params$zeta2,
    #                        Uop = Uop,
    #                        tau2 = params$tau2)
    slice_freq = 100
    if (iter %% slice_freq == 0) {
      params$psi <- update_psi_slice(
        psi    = params$psi,
        omega  = params$omega,
        beta   = params$beta,
        sigma2 = params$sigma2,
        zeta2  = params$zeta2,
        Uop    = Uop,
        tau2   = params$tau2
      )
    } else {
      params$psi <- update_psi_elliptical(
        psi        = params$psi,
        omega      = params$omega,
        beta       = params$beta,
        sigma2     = params$sigma2,
        zeta2      = params$zeta2,
        Uop        = Uop,
        tau2       = params$tau2,
        init_ridge = 1e-6
      )
    }



    # acc_total  <- acc_total  + as.integer(res$accepted)
    # prop_total <- prop_total + 1L

    # burn_in <- 2000
    # if (iter <= burn_in) {
    #   target <- 0.57;  gamma <- 0.02
    #   params$step_var <- params$step_var *
    #     exp(gamma * (as.integer(res$accepted) - target))
    # }
    ## ---- σ² | rest  --------------------------------------------
    params$sigma2 <- sample_sigma2(
      y, X, params$beta,
      a_sigma = params$a_sigma,
      b_sigma = params$b_sigma,
      Uop     = Uop,
      psi     = params$psi,
      zeta2   = params$zeta2)

    ## ---- ζ² | rest  --------------------------------------------
    # params$zeta2 <- sample_zeta2(
    #   beta   = params$beta,
    #   sigma2 = params$sigma2,
    #   Uop    = Uop,
    #   psi    = params$psi,
    #   nu     = params$nu,
    #   m0     = 12)
    params$zeta2 <- params$zeta2
    params$zeta2 <- sample_zeta2(
      beta   = params$beta,
      sigma2 = params$sigma2,
      Uop    = Uop,
      psi    = params$psi,
      nu    = params$nu)


    ## ---- ν  | ζ²  ---------------------------------------------
    params$nu <- sample_nu(params$zeta2, params$par_ratio, params$sigma2)
    #params$nu<- min( params$nu, 10)

    ## ---- ω, ρ | rest  ------------------------------------------
    omega_rho <- gibbs_sampler_omega_rho(
      alpha_mat = params$alpha_mat,
      kappa_mat = params$kappa_mat,
      q_mat     = params$q_mat,
      psi       = params$psi,
      tau2      = params$tau2,
      omega_curr = params$omega,
      rho_curr  = params$rho)

    omega_raw <- omega_rho$omega         # symmetric, shrinkage-updated
    diag(omega_raw) <- 0                 # just to be sure
    # --- 2. store raw version for later diagnostics / plotting -----------------
    params$omega_raw <- omega_raw        
    #params$omega <- omega_rho$omega   # keep for NEXT iteration
    params$rho       <- omega_rho$rho
    # Row-standardise only the nonzero rows
    #W0 <- params$omega

    diag(omega_raw) <- 0
    rs <- rowSums(omega_raw)
    omega_raw[rs > 0, ] <- sweep(
      omega_raw[rs > 0, , drop = FALSE],  # only those rows with any edges
      1,
      rs[rs > 0],
      "/"
    )
    
    params$omega_deg <- rs
    params$omega <-omega_raw
    # now recalibrate the τ² hyperprior
    sum_w     <- sum(omega_raw[upper.tri(omega_raw)])
    params$a_tau <- 1
    params$b_tau <- 0.5 * sum_w    # = number_of_connected_nodes/2


    ## ---- τ² | rest  (global scalar here) -----------------------
    params$tau2 <- sample_tau2(
      a_tau = params$a_tau,
      b_tau = params$b_tau,
      omega = params$omega,
      psi   = params$psi)
    #params$tau2 = 1
    ## ---- save current state ------------------------------------
    samples[[iter]] <- list(beta = params$beta, eta = params$eta, psi = params$psi,
                            sigma2 = params$sigma2, zeta2 = params$zeta2, nu =params$nu,
                            rho = params$rho,
                            omega_raw = pack_omega(params$omega_raw),
                            tau2 = params$tau2, omega = pack_omega(params$omega), omega_deg = params$omega_deg)
  }
  #acc_rate_overall <- acc_total / prop_total

  list(samples = samples,
       acc_total = acc_total,
       prop_total = prop_total)#,
  #acc_rate_overall = acc_rate_overall)

  # return list of length n_iter
}


#Elliptical slice
update_psi_elliptical <- function(psi, omega, beta,
                                  sigma2, zeta2, Uop, tau2,
                                  init_ridge = 1e-6) {
  p <- length(psi)
  ## 1) build intrinsic CAR precision Q_int = D - Wscaled
  if (length(tau2) == 1L) {
    Wscaled <- omega / tau2
  } else {
    Wscaled <- sweep(omega, 1, tau2, "/")
  }
  Q_int <- diag(rowSums(Wscaled)) - Wscaled
  ## 2) get a PD Cholesky via your helper
  R <- chol_pd(Q_int, init_ridge = init_ridge)  # will add ridge as needed

  ## 3) draw the auxiliary nu ~ N(0, Q^{-1})
  z  <- rnorm(p)
  nu <- backsolve(R, z)  # solves t(R) x = z

  ## 4) define the non-Gaussian part of the log-target
  logL <- function(ψ) {
    Ub <- Uop$apply_U(beta)
    ll1 <- -0.5 * sum((Ub^2) * exp(-2 * ψ)) / (sigma2 * zeta2)
    ll2 <- - sum(ψ)
    ll1 + ll2
  }

  ## 5) slice threshold
  logy <- logL(psi) + log(runif(1))

  ## 6) bracket for θ
  theta     <- runif(1, 0, 2*pi)
  th_lo <- theta - 2*pi
  th_hi <- theta

  ## 7) propose on the ellipse until accepted
  repeat {
    psi_prop <- psi * cos(theta) + nu * sin(theta)
    if (logL(psi_prop) > logy) return(psi_prop)
    # shrink bracket
    if (theta < 0) th_lo <- theta else th_hi <- theta
    theta <- runif(1, th_lo, th_hi)
  }
}

# one pass of coordinate-wise slice for ψ
# make sure SeBR is installed and loaded
#install.packages("SeBR")
library(SeBR)

update_psi_slice <- function(psi, omega, beta,
                             sigma2, zeta2, Uop, tau2,
                             w = 1, m = 10) {
  p <- length(psi)
  for (j in seq_len(p)) {
    # 1) define log-posterior as a function of a single coordinate xj
    logf_j <- function(xj) {
      tmp       <- psi
      tmp[j]    <- xj
      log_posterior(tmp,
                    omega, beta,
                    sigma2, zeta2,
                    Uop, tau2)
    }
    # 2) call SeBR:::uni.slice with the correct argument name 'g'
    psi[j] <- SeBR:::uni.slice(
      x0    = psi[j],
      g     = logf_j,
      w     = w,
      m     = m,
      lower = -Inf,
      upper = +Inf
    )
  }
  psi
}
