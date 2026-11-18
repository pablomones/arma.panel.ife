################################################################################
##  QMLE and GMM for dynamic panels with interactive effects and ARMA errors.
##
##  OPTIMIZATION PARAMETER LAYOUT:
##    psi = [AR reflections (or phi), MA reflections (or theta),
##           sigma2_{p+1-q}..sigma2_T,          (T-p+q entries, heterosk.)
##           vech(L_Psi),                        (r*(r+1)/2 entries, Psi=L_Psi L_Psi')
##           vec(A2), A2 = F2 L_Psi,             ((T-r)*r entries)
##           vech(L_eta),                         (p*(p+1)/2 entries)
##           sigma_eta_s_eps_k cross-terms ]     (K(p,q) entries)
##  Estimates and inference are reported in the paper's parameterization,
##  which decodes dynamic reflections into physical phi/theta, replaces
##  vech(L_Psi) by vech(Psi), and replaces the internal Schur-factor
##  coordinates vech(L_eta) by the Sigma_eta block of Pi_(p,q): all initial
##  variances first, then covariances (s,t), s<t, in lexicographic (s,t) order.
##  The eta/innovation cross-covariances follow in their original paper order.
##  Inference uses physical dynamics/Pi and stable triangular factor coordinates.
##  The delta method maps the sandwich covariance to the reported parameters.
##
##  CONSTRAINTS (box constraints plus internal feasibility parameterization):
##    sigma2_t          >= 1e-4
##    diag(L_Psi)       >= 1e-3  (scalar Psi >= 1e-6)
##    diag(L_eta)       >= 1e-3
##    stationary=TRUE, invertible=TRUE: reflections in [-(1-1e-4),1-1e-4]
##    FALSE: corresponding physical coefficients in [-bound,bound], default bound=5
##    off-diagonal L_Psi, F2, cross-terms: unconstrained
##
##  GMM with regressors uses numDeriv for moment Jacobians.
################################################################################

if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("Package 'numDeriv' is required. Install with install.packages('numDeriv').")
}
library(numDeriv)

# ============================================================================ #
#  SECTION 1. MATRIX CONSTRUCTION HELPERS
# ============================================================================ #

# ---------------------------------------------------------------------------- #
#  1.1  B_(p)  and  Gamma_(p) 
# ---------------------------------------------------------------------------- #

build_B_p <- function(phi, big_T) {
  p <- length(phi)
  if (p == 0L) return(diag(big_T))   # pure MA: B_(0) = I_T
  stopifnot(big_T > p)

  T_minus_p <- big_T - p
  B <- diag(big_T)   # start with I_T; top-left p×p block stays as I_p

  # ---- C_(p): (T-p)×(T-p) lower-triangular Toeplitz  ----------------------
  # Diagonal = 1 (already set by diag(big_T)).
  # k-th subdiagonal = -phi_k for k = 1..p.
  for (k in seq_len(p)) {
    if (k < T_minus_p) {
      rows_C <- (p + k + 1L):big_T          # rows of B corresponding to C
      cols_C <- (p + 1L):(big_T - k)        # matching cols
      B[cbind(rows_C, cols_C)] <- -phi[k]
    }
  }

  # ---- A_(p): (T-p)×p block, bottom-left of B  ----------------------------
  # Row k of A_(p) (k=1..min(p, T-p)):
  #   nonzero in columns k..p with values (-phi_p, -phi_{p-1}, ..., -phi_k)
  #   (reverse-descending phi, prepended with zeros for cols 1..k-1).
  for (k in seq_len(min(p, T_minus_p))) {
    B[p + k, k:p] <- -rev(phi[k:p])
  }

  return(B)
}

# Build Gamma_(p) = B_(p)^{-1}  [Paper eq. 11]
build_Gamma_p <- function(phi, big_T) {
  B <- build_B_p(phi, big_T)
  solve(B)
}

# ---------------------------------------------------------------------------- #
#  1.2  Theta_(q) MA mapping matrix  [Paper eq. 18]
# ---------------------------------------------------------------------------- #

build_Theta_q <- function(theta, big_T) {
  q     <- length(theta)
  Theta <- matrix(0, nrow = big_T, ncol = big_T + q)
  coefs <- c(rev(theta), 1)   # (theta_q, ..., theta_1, theta_0=1), length q+1
  for (t in seq_len(big_T)) {
    Theta[t, t:(t + q)] <- coefs
  }
  return(Theta)
}

# ---------------------------------------------------------------------------- #
#  1.3  K(p,q): number of cross-covariance nuisance parameters  [eq. 16]
# ---------------------------------------------------------------------------- #

compute_K_pq <- function(p, q) {
  if (p == 0L || q == 0L) return(0L)
  if (q <= p) {
    return(as.integer(q * (q + 1L) / 2L))
  } else {
    return(as.integer(p * (q - p) + p * (p + 1L) / 2L))
  }
}

.chol_diag_positions <- function(r) {
  if (r == 0L) return(integer(0))
  idx <- matrix(0L, r, r)
  idx[lower.tri(idx, diag = TRUE)] <- seq_len(r * (r + 1L) / 2L)
  as.integer(diag(idx))
}

.Psi_from_lower <- function(Psi_lower, r) {
  if (r == 0L) return(matrix(0, 0, 0))
  L_Psi <- matrix(0, r, r)
  L_Psi[lower.tri(L_Psi, diag = TRUE)] <- Psi_lower
  L_Psi %*% t(L_Psi)
}

# ---------------------------------------------------------------------------- #
#  1.4  Sigma_u  [Paper eq. 18]
# ---------------------------------------------------------------------------- #

build_Sigma_u <- function(theta, sigma2_vec, p, big_T, heteroskedastic) {
  q         <- length(theta)
  Tm        <- big_T - p          # number of non-projected rows
  if (Tm == 0L) return(matrix(0, 0, 0))

  n_identifiable <- big_T - p + q   # = T-p+q

  if (!heteroskedastic) {
    # Homoskedastic: sigma^2_t = sigma^2 for all t.
    sigma2_scalar <- sigma2_vec[1L]
    # Sigma_u = sigma^2 * Theta_sub %*% t(Theta_sub)  (Theta_sub has unit entries)
    # Equivalently build via the formula directly.
    Theta_full <- build_Theta_q(theta, big_T)
    Theta_sub  <- Theta_full[(p + 1L):big_T, , drop = FALSE]  # Tm × (T+q)
    return(sigma2_scalar * Theta_sub %*% t(Theta_sub))
  }

  # Heteroskedastic: sigma2_vec has length T-p+q.
  if (length(sigma2_vec) == 1L) {
    sigma2_vec <- rep(sigma2_vec, n_identifiable)
  }
  if (length(sigma2_vec) != n_identifiable) {
    stop(sprintf(
      paste0("build_Sigma_u: sigma2_vec must have length T-p+q = %d ",
             "(got %d). Supply sigma^2 for epochs p+1-q=%d through T=%d."),
      n_identifiable, length(sigma2_vec), p + 1L - q, big_T))
  }

  # Build D: (T+q)×(T+q) diagonal, D[k] = sigma^2_{epoch k-q}.
  # D[1..p] = 0 (unidentifiable pre-sample epochs).
  # D[p+1 .. T+q] = sigma2_vec[1 .. T-p+q].
  d_vec  <- c(rep(0, p), sigma2_vec)   # length T+q
  D      <- diag(d_vec, nrow = big_T + q)

  Theta_full <- build_Theta_q(theta, big_T)
  Theta_sub  <- Theta_full[(p + 1L):big_T, , drop = FALSE]   # Tm × (T+q)

  Sigma_u <- Theta_sub %*% D %*% t(Theta_sub)
  return(Sigma_u)
}

# ---------------------------------------------------------------------------- #
#  1.5  Unpack Pi_(p,q) nuisance vector  [eq. 17]
# ---------------------------------------------------------------------------- #

unpack_Pi <- function(Pi_pq, p, q) {
  if (p == 0L) {
    return(list(Sigma_eta = matrix(0, 0, 0),
                cross_s = integer(0),
                cross_k = integer(0),
                cross_v = numeric(0)))
  }

  pos <- 1L

  # (a) Variances of projected errors sigma_eta_s^2, s = 1..p
  var_eta <- Pi_pq[pos:(pos + p - 1L)];  pos <- pos + p

  # (b) Off-diagonal covariances of projected errors
  n_offdiag <- p * (p - 1L) / 2L
  cov_eta   <- if (n_offdiag > 0L) Pi_pq[pos:(pos + n_offdiag - 1L)] else numeric(0)
  pos       <- pos + n_offdiag

  # Reconstruct Sigma_eta_(p) (symmetric p×p)
  Sigma_eta <- matrix(0, p, p)
  diag(Sigma_eta) <- var_eta
  if (n_offdiag > 0L) {
    idx_off <- 0L
    for (ii in seq_len(p - 1L)) {
      for (jj in (ii + 1L):p) {
        idx_off <- idx_off + 1L
        Sigma_eta[ii, jj] <- cov_eta[idx_off]
        Sigma_eta[jj, ii] <- cov_eta[idx_off]
      }
    }
  }

  # (c) Cross-covariance parameters sigma_{eta_s, eps_k}
  #     s from max(1, p+1-q) to p;  k from p+1-q to s  (for each s).
  K_pq    <- compute_K_pq(p, q)
  cross_s <- integer(0)
  cross_k <- integer(0)
  cross_v <- numeric(0)

  if (K_pq > 0L) {
    s_min <- max(1L, p + 1L - q)
    k_min <- p + 1L - q              # same lower bound for k for all s
    for (s in s_min:p) {
      for (k in k_min:s) {
        cross_s <- c(cross_s, s)
        cross_k <- c(cross_k, k)
        cross_v <- c(cross_v, Pi_pq[pos])
        pos     <- pos + 1L
      }
    }
  }

  list(Sigma_eta = Sigma_eta,
       cross_s   = cross_s,
       cross_k   = cross_k,
       cross_v   = cross_v)
}

# Convert between the paper's Pi_(p,q) coordinates and optimization
# coordinates that impose the exact Schur-complement feasibility condition.
# Internally,
#   Sigma_eta = A D_0^{-1} A' + L_eta L_eta',
# with a positive diagonal in L_eta.  The cross-covariances in A remain in
# their original paper ordering.
.build_A_pq <- function(cross_v, p, q) {
  A <- matrix(0, p, q)
  if (p == 0L || q == 0L || length(cross_v) == 0L) return(A)

  pos <- 1L
  first_epoch <- p + 1L - q
  for (s in max(1L, first_epoch):p) {
    for (k in first_epoch:s) {
      A[s, k - first_epoch + 1L] <- cross_v[pos]
      pos <- pos + 1L
    }
  }
  A
}

.early_innovation_variances <- function(sigma2, p, q, big_T,
                                        heteroskedastic = TRUE) {
  if (q == 0L) return(numeric(0))
  n_eps <- big_T - p + q
  sigma2_eps <- if (heteroskedastic) {
    if (length(sigma2) != n_eps) stop("Incorrect innovation-variance length.")
    sigma2
  } else {
    if (length(sigma2) < 1L) stop("Missing innovation variance.")
    rep(sigma2[1L], n_eps)
  }
  sigma2_eps[seq_len(q)]
}

.pack_Sigma_eta <- function(Sigma_eta) {
  p <- nrow(Sigma_eta)
  if (p == 0L) return(numeric(0))
  out <- diag(Sigma_eta)
  if (p > 1L) {
    for (ii in seq_len(p - 1L)) {
      for (jj in (ii + 1L):p) out <- c(out, Sigma_eta[ii, jj])
    }
  }
  out
}

.Pi_internal_to_paper <- function(Pi_internal, sigma2, p, q, big_T,
                                  heteroskedastic = TRUE) {
  if (p == 0L) return(numeric(0))
  n_eta <- p * (p + 1L) / 2L
  K_pq <- compute_K_pq(p, q)
  if (length(Pi_internal) != n_eta + K_pq)
    stop("Incorrect internal Pi_(p,q) length.")

  eta_lower <- Pi_internal[seq_len(n_eta)]
  cross_v <- if (K_pq > 0L) Pi_internal[n_eta + seq_len(K_pq)] else numeric(0)
  L_eta <- matrix(0, p, p)
  L_eta[lower.tri(L_eta, diag = TRUE)] <- eta_lower
  if (any(!is.finite(L_eta)) || any(diag(L_eta) <= 0))
    stop("The diagonal of L_eta must be positive.")

  Sigma_eta <- L_eta %*% t(L_eta)
  if (q > 0L) {
    v0 <- .early_innovation_variances(
      sigma2, p, q, big_T, heteroskedastic
    )
    if (any(!is.finite(v0)) || any(v0 <= 0))
      stop("Innovation variances must be positive.")
    A <- .build_A_pq(cross_v, p, q)
    Sigma_eta <- Sigma_eta + tcrossprod(sweep(A, 2L, sqrt(v0), "/"))
  }

  c(.pack_Sigma_eta(Sigma_eta), cross_v)
}

.Pi_paper_to_internal <- function(Pi_pq, sigma2, p, q, big_T,
                                  heteroskedastic = TRUE) {
  if (p == 0L) return(numeric(0))
  nu <- unpack_Pi(Pi_pq, p, q)
  Schur <- nu$Sigma_eta
  if (q > 0L) {
    v0 <- .early_innovation_variances(
      sigma2, p, q, big_T, heteroskedastic
    )
    if (any(!is.finite(v0)) || any(v0 <= 0))
      stop("Innovation variances must be positive.")
    A <- .build_A_pq(nu$cross_v, p, q)
    Schur <- Schur - tcrossprod(sweep(A, 2L, sqrt(v0), "/"))
  }
  Schur <- (Schur + t(Schur)) / 2
  L_eta <- t(chol(Schur))
  c(L_eta[lower.tri(L_eta, diag = TRUE)], nu$cross_v)
}

# ---------------------------------------------------------------------------- #
#  1.6  Sigma_xi_(p,q)  [Paper eq. 13]
# ---------------------------------------------------------------------------- #

build_Sigma_xi <- function(phi, theta, sigma2_vec, Pi_pq, big_T,
                            heteroskedastic = TRUE) {
  p <- length(phi)
  q <- length(theta)

  # Pure MA: no projection; Sigma_xi = Sigma_u^+
  if (p == 0L) {
    return(build_Sigma_u(theta, sigma2_vec, 0L, big_T, heteroskedastic))
  }

  stopifnot(big_T > p)
  Tm <- big_T - p   # number of non-projected periods

  # ---- Bottom-right block: Sigma_u^+_(p,q)  [eq. 18] ----------------------
  Sigma_u_block <- build_Sigma_u(theta, sigma2_vec, p, big_T, heteroskedastic)

  # ---- Unpack nuisance -------------------------------------------------------
  nu        <- unpack_Pi(Pi_pq, p, q)
  Sigma_eta <- nu$Sigma_eta   # p×p

  # ---- Cross-covariance block  [eq. 15] ------------------------------------
  # Sigma_{u,eta}_(p,q) is p×Tm; entry (s, t-p) for s=1..p, t=p+1..T.
  #
  # Cov(eta_{is}, u_{it}) = sum_{l=0}^{q-(t-s)} theta_{t-s+l} * sigma_{eta_s, eps_{s-l}}
  #
  # We iterate over stored (s, k) pairs.  Each pair provides sigma_{eta_s, eps_k}.
  # Given k, l = s - k >= 0, so theta index = t - s + l = t - s + (s - k) = t - k.
  # The sum over l includes term l = s - k as long as 0 <= l <= q - (t-s),
  # i.e., 0 <= s - k <= q - (t - s), i.e., t - s <= q - (s - k) = q - l.
  # This is already guaranteed by the sparsity condition t <= s + q.
  # So: Cov(eta_{is}, u_{it}) += theta_{t-k} * sigma_{eta_s, eps_k}
  #     whenever p+1 <= t <= s+q  AND  0 <= t-k <= q (i.e., theta_{t-k} exists).
  theta_full <- c(1, theta)   # theta_0=1, theta_1..theta_q; 1-indexed as theta_full[m+1]=theta_m
  Sigma_cross <- matrix(0, nrow = p, ncol = Tm)  # p × (T-p)

  if (length(nu$cross_s) > 0L) {
    for (idx in seq_along(nu$cross_s)) {
      s   <- nu$cross_s[idx]   # s in 1..p
      k   <- nu$cross_k[idx]   # k in p+1-q .. s
      val <- nu$cross_v[idx]   # sigma_{eta_s, eps_k}

      # t ranges over p+1 .. s+q (condition from eq. 15)
      t_max <- min(big_T, s + q)
      if (p + 1L > t_max) next   # no valid t

      for (t in (p + 1L):t_max) {
        m <- t - k   # index into theta_full: theta_{t-k}
        if (m < 0L || m > q) next   # theta_m does not exist (out of [0,q])
        Sigma_cross[s, t - p] <- Sigma_cross[s, t - p] + theta_full[m + 1L] * val
      }
    }
  }

  # ---- Assemble Sigma_xi_(p,q)  [eq. 13] -----------------------------------
  Sigma_xi <- matrix(0, big_T, big_T)
  Sigma_xi[1:p,              1:p             ] <- Sigma_eta
  Sigma_xi[1:p,              (p+1L):big_T    ] <- Sigma_cross
  Sigma_xi[(p+1L):big_T,     1:p             ] <- t(Sigma_cross)
  Sigma_xi[(p+1L):big_T,     (p+1L):big_T   ] <- Sigma_u_block

  return(Sigma_xi)
}

# ---------------------------------------------------------------------------- #
#  1.7  Full model covariance Sigma_(p,q,r)(vartheta)  [Paper eq. 19]
# ---------------------------------------------------------------------------- #
build_Sigma_model <- function(vartheta, p, q, r, big_T, heteroskedastic = TRUE) {

  phi   <- if (p > 0L) vartheta$phi   else numeric(0)
  theta <- if (q > 0L) vartheta$theta else numeric(0)
  stopifnot(length(phi) == p, length(theta) == q)
  if (length(vartheta$F2) != (big_T - r) * r)
    stop("F2 has the wrong length for the supplied r and T.")
  if (length(vartheta$Psi_lower) != r * (r + 1L) / 2L)
    stop("Psi_lower has the wrong length for the supplied r.")

  # Innovation variances
  sigma2_vec <- vartheta$sigma2
  if (!heteroskedastic) sigma2_vec <- sigma2_vec[1L]

  # Factor path F_(r): T×r with top r×r block = I_r (Assumption A1).
  # NOTE on r=0: diag(0) in R returns a 0×0 matrix (NOT numeric(0)), so
  # rbind(diag(0), matrix(numeric(0), big_T, 0)) produces a T×0 matrix
  # correctly.  F_r %*% Psi_r %*% t(F_r) is then a T×T zero matrix as
  # required.  The r=0 branch sets factor_part = 0 directly for clarity.
  if (r > 0L) {
    F2  <- matrix(vartheta$F2, nrow = big_T - r, ncol = r)
    F_r <- rbind(diag(r), F2)    # T×r
  } else {
    F_r <- matrix(0, big_T, 0)   # T×0 — no factor contribution
  }

  # Factor covariance Psi_(r): recover r×r PD matrix from lower Cholesky.
  # r=0: Psi_lower is numeric(0); no matrix needed; factor_part = 0 below.
  if (r > 0L) {
    L_Psi <- matrix(0, r, r)
    L_Psi[lower.tri(L_Psi, diag = TRUE)] <- vartheta$Psi_lower
    # Enforce PD: diagonal entries of L must be positive
    # (enforced via box constraint in optimizer; guard here just in case)
    if (any(diag(L_Psi) <= 0)) return(NULL)
    Psi_r <- .Psi_from_lower(vartheta$Psi_lower, r)
  } else {
    Psi_r <- matrix(0, 0, 0)
  }

  # AR filter matrices
  Gamma_p <- build_Gamma_p(phi, big_T)   # T×T

  # Projected-error covariance Sigma_xi_(p,q)
  Sigma_xi <- build_Sigma_xi(phi, theta, sigma2_vec, vartheta$Pi_pq,
                              big_T, heteroskedastic)

  # Factor contribution
  factor_part <- if (r > 0L) F_r %*% Psi_r %*% t(F_r) else matrix(0, big_T, big_T)

  # Full covariance (eq. 19)
  inner <- factor_part + Sigma_xi
  Sigma <- Gamma_p %*% inner %*% t(Gamma_p)
  Sigma <- (Sigma + t(Sigma)) / 2   # symmetrize (numerical noise)

  return(Sigma)
}

# ============================================================================ #
#  SECTION 2.  PARAMETER COUNTING AND PACKING
# ============================================================================ #

count_params <- function(p, q, r, big_T, heteroskedastic = TRUE) {
  n_phi   <- p
  n_theta <- q
  n_sigma <- if (heteroskedastic) (big_T - p + q) else 1L
  n_Psi   <- r * (r + 1L) / 2L
  n_F2    <- (big_T - r) * r
  n_eta   <- p * (p + 1L) / 2L
  n_cross <- compute_K_pq(p, q)
  as.integer(n_phi + n_theta + n_sigma + n_Psi + n_F2 + n_eta + n_cross)
}

.seq_names <- function(prefix, n) {
  if (n > 0L) paste0(prefix, seq_len(n)) else character(0)
}

.parameter_names <- function(p, q, r, big_T, heteroskedastic, internal,
                                  stationary = TRUE, invertible = TRUE) {
  n_eta <- p * (p + 1L) / 2L
  K_pq <- compute_K_pq(p, q)
  sigma_names <- if (heteroskedastic)
    paste0("sigma2_", seq.int(p + 1L - q, big_T)) else "sigma2"
  c(.seq_names(if (internal && stationary) "kappa_ar_" else "phi_", p),
    .seq_names(if (internal && invertible) "kappa_ma_" else "theta_", q), sigma_names,
    .seq_names(if (internal) "L_Psi_" else "Psi_", r * (r + 1L) / 2L),
    .seq_names(if (internal) "A2_" else "F2_", (big_T - r) * r),
    .seq_names(if (internal) "L_eta_" else "Sigma_eta_", n_eta),
    .seq_names("eta_eps_", K_pq))
}

# Internal optimization-vector layout:
#
#  Slot          Size           Constraint
#  AR block      p              reflections in [-(1-1e-4),1-1e-4] if stationary; otherwise physical phi
#  MA block      q              reflections in [-(1-1e-4),1-1e-4] if invertible; otherwise physical theta
#  sigma2        T-p+q (hetero) >= 1e-4
#                1 (homo)       >= 1e-4
#  Psi_lower     r*(r+1)/2      diagonal entries >= 1e-3; off-diag unconstrained
#  A2=F2 L_Psi   (T-r)*r        unconstrained
#  Pi_internal:
#    vech(L_eta) p*(p+1)/2      diagonal entries >= 1e-3; off-diag unconstrained
#    cross       K(p,q)         unconstrained

#' Pack structural parameters into optimization vector psi
pack_params <- function(vartheta, p, q, r, big_T, heteroskedastic = TRUE,
                                  stationary = TRUE, invertible = TRUE) {
  phi   <- if (p > 0L) vartheta$phi   else numeric(0)
  theta <- if (q > 0L) vartheta$theta else numeric(0)

  sigma2    <- if (heteroskedastic) vartheta$sigma2 else vartheta$sigma2[1L]
  Psi_lower <- if (r > 0L) vartheta$Psi_lower else numeric(0)
  if (r > 0L) {
    L_Psi <- matrix(0, r, r)
    L_Psi[lower.tri(L_Psi, diag = TRUE)] <- Psi_lower
    A2 <- as.vector(matrix(vartheta$F2, big_T - r, r) %*% L_Psi)
  } else {
    A2 <- numeric(0)
  }

  Pi_internal <- .Pi_paper_to_internal(
    vartheta$Pi_pq, sigma2, p, q, big_T, heteroskedastic
  )

  dynamic <- .encode_dynamics(phi, theta, stationary, invertible)
  c(dynamic$ar, dynamic$ma, sigma2, Psi_lower, A2, Pi_internal)
}

#' Unpack optimization vector psi into structural parameters
unpack_params <- function(psi, p, q, r, big_T, heteroskedastic = TRUE,
                                  stationary = TRUE, invertible = TRUE) {
  pos <- 1L

  phi   <- if (p > 0L) { v <- psi[pos:(pos+p-1L)]; pos <- pos+p; v } else { numeric(0) }
  theta <- if (q > 0L) { v <- psi[pos:(pos+q-1L)]; pos <- pos+q; v } else { numeric(0) }

  dynamic <- .decode_dynamics(phi, theta, stationary, invertible)
  phi <- dynamic$phi
  theta <- dynamic$theta

  n_sigma <- if (heteroskedastic) (big_T - p + q) else 1L
  sigma2  <- psi[pos:(pos + n_sigma - 1L)];  pos <- pos + n_sigma

  n_Psi     <- r * (r + 1L) / 2L
  Psi_lower <- if (r > 0L) { v <- psi[pos:(pos+n_Psi-1L)]; pos <- pos+n_Psi; v } else { numeric(0); pos <- pos + n_Psi; numeric(0) }

  n_F2 <- (big_T - r) * r
  A2   <- if (n_F2 > 0L) { v <- psi[pos:(pos+n_F2-1L)]; pos <- pos+n_F2; v } else { pos <- pos + n_F2; numeric(0) }
  if (r > 0L) {
    L_Psi <- matrix(0, r, r)
    L_Psi[lower.tri(L_Psi, diag = TRUE)] <- Psi_lower
    F2 <- as.vector(matrix(A2, big_T - r, r) %*% solve(L_Psi))
  } else {
    F2 <- numeric(0)
  }

  n_eta   <- p * (p + 1L) / 2L
  K_pq    <- compute_K_pq(p, q)
  n_Pi    <- n_eta + K_pq
  Pi_internal <- if (n_Pi > 0L) psi[pos:(pos + n_Pi - 1L)] else numeric(0)
  Pi_pq <- .Pi_internal_to_paper(
    Pi_internal, sigma2, p, q, big_T, heteroskedastic
  )

  list(phi       = phi,
       theta     = theta,
       sigma2    = sigma2,
       Psi_lower = Psi_lower,
       F2        = F2,
       A2        = A2,
       Pi_pq     = Pi_pq,
       Pi_internal = Pi_internal)
}

# Build Sigma directly from the internal optimization coordinates.  This is
# algebraically identical to build_Sigma_model(unpack_params(...)), but avoids
# recovering F2 and the paper-coordinate Pi_(p,q) at every objective evaluation.
.build_Sigma_internal <- function(psi, p, q, r, big_T,
                                  heteroskedastic = TRUE,
                                  stationary = TRUE, invertible = TRUE) {
  if (length(psi) != count_params(p, q, r, big_T, heteroskedastic))
    return(NULL)

  pos <- 1L
  phi <- if (p > 0L) {
    v <- psi[pos:(pos + p - 1L)]; pos <- pos + p; v
  } else numeric(0)
  theta <- if (q > 0L) {
    v <- psi[pos:(pos + q - 1L)]; pos <- pos + q; v
  } else numeric(0)

  dynamic <- .decode_dynamics(phi, theta, stationary, invertible)
  phi <- dynamic$phi
  theta <- dynamic$theta

  n_sigma <- if (heteroskedastic) big_T - p + q else 1L
  sigma2 <- psi[pos:(pos + n_sigma - 1L)]; pos <- pos + n_sigma

  n_Psi <- r * (r + 1L) / 2L
  Psi_lower <- if (n_Psi > 0L) {
    v <- psi[pos:(pos + n_Psi - 1L)]; pos <- pos + n_Psi; v
  } else numeric(0)

  n_A2 <- (big_T - r) * r
  A2 <- if (n_A2 > 0L) {
    v <- psi[pos:(pos + n_A2 - 1L)]; pos <- pos + n_A2; v
  } else numeric(0)

  n_eta <- p * (p + 1L) / 2L
  eta_lower <- if (n_eta > 0L) {
    v <- psi[pos:(pos + n_eta - 1L)]; pos <- pos + n_eta; v
  } else numeric(0)

  K_pq <- compute_K_pq(p, q)
  cross_v <- if (K_pq > 0L) psi[pos:(pos + K_pq - 1L)] else numeric(0)

  if (any(!is.finite(psi)) || any(sigma2 <= 0)) return(NULL)

  # F Psi F' = G G', where G = (L_Psi', A2')'.
  if (r > 0L) {
    L_Psi <- matrix(0, r, r)
    L_Psi[lower.tri(L_Psi, diag = TRUE)] <- Psi_lower
    if (any(diag(L_Psi) <= 0)) return(NULL)
    G <- rbind(L_Psi, matrix(A2, big_T - r, r))
    factor_part <- tcrossprod(G)
  } else {
    factor_part <- matrix(0, big_T, big_T)
  }

  # MA mapping for the innovations that enter u_(p+1),...,u_T.
  n_eps <- big_T - p + q
  sigma2_eps <- if (heteroskedastic) sigma2 else rep(sigma2[1L], n_eps)
  Theta_full <- build_Theta_q(theta, big_T)
  Theta_pq <- Theta_full[(p + 1L):big_T,
                          (p + 1L):(big_T + q), drop = FALSE]
  weighted_Theta <- sweep(Theta_pq, 2L, sqrt(sigma2_eps), "*")
  Sigma_u <- tcrossprod(weighted_Theta)

  if (p > 0L) {
    L_eta <- matrix(0, p, p)
    L_eta[lower.tri(L_eta, diag = TRUE)] <- eta_lower
    if (any(diag(L_eta) <= 0)) return(NULL)

    A <- .build_A_pq(cross_v, p, q)
    Sigma_eta <- tcrossprod(L_eta)
    if (q > 0L) {
      v0 <- sigma2_eps[seq_len(q)]
      Sigma_eta <- Sigma_eta +
        tcrossprod(sweep(A, 2L, sqrt(v0), "/"))
      Sigma_cross <- A %*%
        t(Theta_pq[, seq_len(q), drop = FALSE])
    } else {
      Sigma_cross <- matrix(0, p, big_T - p)
    }

    Sigma_xi <- rbind(
      cbind(Sigma_eta, Sigma_cross),
      cbind(t(Sigma_cross), Sigma_u)
    )
  } else {
    Sigma_xi <- Sigma_u
  }

  Gamma_p <- build_Gamma_p(phi, big_T)
  Sigma <- Gamma_p %*% (factor_part + Sigma_xi) %*% t(Gamma_p)
  (Sigma + t(Sigma)) / 2
}

.validate_coefficient_bound <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value <= 0)
    stop(name, " must be a finite positive scalar.", call. = FALSE)
  as.numeric(value)
}


# Strict reflection domain shared by every constrained estimator and start.
.dynamic_reflection_limit <- 1 - 1e-4

# Validate scalar finite-difference and initialization optimizer controls.
.validate_positive_optimizer_control <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value <= 0)
    stop(name, " must be a finite positive scalar.", call. = FALSE)
  as.numeric(value)
}

# Validate requested dynamic domains independently of covariance feasibility.
.validate_dynamic_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value))
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  value
}

# Inverse Schur recursion; the optimization margin is validated separately.
.ar_to_pacf <- function(phi, tol = 1e-8) {
  phi <- as.numeric(phi)
  if (any(!is.finite(phi))) stop("AR coefficients must be finite.")
  p <- length(phi)
  if (p == 0L) return(numeric(0))
  original <- phi
  kappa <- numeric(p)
  for (m in p:1L) {
    k <- phi[m]
    if (!is.finite(k) || abs(k) >= 1)
      stop("Coefficients must define a strictly stationary/invertible polynomial.",
           call. = FALSE)
    kappa[m] <- k
    if (m > 1L) {
      previous <- phi[seq_len(m - 1L)]
      phi <- numeric(m - 1L)
      for (j in seq_len(ceiling((m - 1L) / 2))) {
        other <- m - j
        if (j == other) {
          phi[j] <- previous[j] / (1 - k)
        } else {
          phi[other] <- (previous[other] + k * previous[j]) /
            ((1 - k) * (1 + k))
          phi[j] <- previous[j] + k * phi[other]
        }
      }
    }
  }
  if (max(abs(.pacf_to_ar(kappa) - original)) >
      8 * tol * max(1, max(abs(original))))
    stop("The inverse reflection mapping failed its reconstruction check.",
         call. = FALSE)
  kappa
}


# Near multiple reflection bounds, physical coefficients lose digits in the
# inverse recursion. A bounded reconstruction is accepted only at roundoff-level
# physical error; it does not project arbitrary supplied coefficients into the domain.
.ar_to_bounded_pacf <- function(phi) {
  phi <- as.numeric(phi)
  p <- length(phi)
  if (p == 0L) return(numeric(0))
  if (any(!is.finite(phi))) stop("AR/MA coefficients must be finite.")
  limit <- .dynamic_reflection_limit
  tolerance <- 128 * .Machine$double.eps * p * max(1, abs(phi))
  clip <- function(k) pmax(-limit, pmin(limit, k))
  error <- function(k) max(abs(.pacf_to_ar(k) - phi))
  raw <- tryCatch(.ar_to_pacf(phi), error = function(e) NULL)
  if (!is.null(raw)) {
    candidate <- clip(raw)
    if (error(candidate) <= tolerance) return(candidate)
  }
  if (abs(phi[p]) > limit + 32 * .Machine$double.eps)
    stop("Constrained reflection coefficients must lie in [-(1-1e-4), 1-1e-4].",
         call. = FALSE)
  last <- clip(phi[p])
  # At fixed kappa_p, orders two and three have segment/quadrilateral images.
  # Their edges provide an exact bounded recovery for ill-conditioned inverses.
  if (p <= 3L) {
    if (p == 1L) {
      if (error(last) <= tolerance) return(last)
    } else {
      if (p == 2L) {
        stable <- c(phi[1L] / (1 - last), last)
      } else {
        a2 <- (phi[2L] + last * phi[1L]) / ((1 - last) * (1 + last))
        a1 <- phi[1L] + last * a2
        stable <- c(a1 / (1 - a2), a2, last)
      }
      if (all(is.finite(stable)) && all(abs(stable) <= limit) &&
          error(stable) <= tolerance) return(stable)
      vertices <- as.matrix(expand.grid(rep(list(c(-limit, limit)), p - 1L)))
      vertices <- cbind(vertices, last)
      for (ii in seq_len(nrow(vertices))) for (j in seq_len(p - 1L)) {
        lo <- hi <- vertices[ii, ]
        hi[j] <- -lo[j]
        a <- .pacf_to_ar(lo)
        direction <- .pacf_to_ar(hi) - a
        denominator <- sum(direction^2)
        fraction <- if (denominator > 0)
          max(0, min(1, sum(direction * (phi - a)) / denominator)) else 0
        candidate <- clip(lo + fraction * (hi - lo))
        if (error(candidate) <= tolerance) return(candidate)
      }
    }
  } else {
    # General fallback: bounded coordinate reconstruction, with the highest
    # reflection fixed by the final physical coefficient.
    candidate <- if (is.null(raw)) rep(0, p) else clip(raw)
    candidate[p] <- last
    for (iteration in seq_len(200L)) {
      for (j in rev(seq_len(p - 1L))) {
        lo <- hi <- candidate
        lo[j] <- -limit; hi[j] <- limit
        a <- .pacf_to_ar(lo)
        direction <- .pacf_to_ar(hi) - a
        denominator <- sum(direction^2)
        if (denominator > 0) {
          fraction <- max(0, min(1, sum(direction * (phi - a)) / denominator))
          candidate[j] <- clip(-limit + 2 * limit * fraction)
        }
      }
      if (error(candidate) <= tolerance) return(candidate)
    }
    # Extremely ill-conditioned corner polynomials can lose the low-order
    # reflection signs. Check the finite corner set for modest model orders.
    if (p <= 12L) {
      corners <- as.matrix(expand.grid(rep(list(c(-limit, limit)), p - 1L)))
      for (ii in seq_len(nrow(corners))) {
        candidate <- c(corners[ii, ], last)
        if (error(candidate) <= tolerance) return(candidate)
      }
    }
  }
  stop("Coefficients are outside the constrained reflection domain or cannot be reconstructed at numerical precision.",
       call. = FALSE)
}

.encode_dynamics <- function(phi, theta, stationary, invertible) {
  stationary <- .validate_dynamic_flag(stationary, "stationary")
  invertible <- .validate_dynamic_flag(invertible, "invertible")
  ar <- if (stationary) .ar_to_bounded_pacf(phi) else as.numeric(phi)
  ma <- if (invertible) .ar_to_bounded_pacf(-theta) else as.numeric(theta)
  limit <- .dynamic_reflection_limit
  roundoff <- 32 * .Machine$double.eps
  if ((stationary && any(abs(ar) > limit + roundoff)) ||
      (invertible && any(abs(ma) > limit + roundoff)))
    stop("Constrained reflection coefficients must lie in [-(1-1e-4), 1-1e-4].",
         call. = FALSE)
  if (stationary) ar <- pmax(-limit, pmin(limit, ar))
  if (invertible) ma <- pmax(-limit, pmin(limit, ma))
  list(ar = ar, ma = ma)
}

.decode_dynamics <- function(ar, ma, stationary, invertible) {
  stationary <- .validate_dynamic_flag(stationary, "stationary")
  invertible <- .validate_dynamic_flag(invertible, "invertible")
  limit <- .dynamic_reflection_limit
  if ((stationary && any(!is.finite(ar) | abs(ar) > limit)) ||
      (invertible && any(!is.finite(ma) | abs(ma) > limit)))
    stop("Reflection coefficients must lie in [-(1-1e-4), 1-1e-4].",
         call. = FALSE)
  list(phi = if (stationary) .pacf_to_ar(ar) else as.numeric(ar),
       theta = if (invertible) -.pacf_to_ar(ma) else as.numeric(ma))
}

make_box_constraints <- function(p, q, r, big_T, heteroskedastic = TRUE,
                                 phi_bound = 5, theta_bound = 5,
                                  stationary = TRUE, invertible = TRUE) {
  stationary <- .validate_dynamic_flag(stationary, "stationary")
  invertible <- .validate_dynamic_flag(invertible, "invertible")
  phi_bound <- .validate_coefficient_bound(phi_bound, "phi_bound")
  theta_bound <- .validate_coefficient_bound(theta_bound, "theta_bound")
  limit <- .dynamic_reflection_limit
  n_sigma   <- if (heteroskedastic) (big_T - p + q) else 1L
  n_Psi     <- r * (r + 1L) / 2L
  diag_pos  <- .chol_diag_positions(r)
  n_F2      <- (big_T - r) * r
  n_eta     <- p * (p + 1L) / 2L
  eta_diag  <- .chol_diag_positions(p)
  K_pq      <- compute_K_pq(p, q)
  n_Pi      <- n_eta + K_pq

  lower <- c(
    rep(if (stationary) -limit else -phi_bound, p), # AR reflections or phi
    rep(if (invertible) -limit else -theta_bound, q), # MA reflections or theta
    rep(1e-4,  n_sigma),     # sigma2_t variance floor
    {                        # Psi_lower: Cholesky diagonal floor
      lb_Psi <- rep(-Inf, n_Psi)
      lb_Psi[diag_pos] <- 1e-3
      lb_Psi
    },
    rep(-Inf,  n_F2),        # A2 = F2 L_Psi: unconstrained
    {                        # L_eta: Cholesky diagonal floor
      lb_Pi <- rep(-Inf, n_Pi)
      if (n_eta > 0L) lb_Pi[eta_diag] <- 1e-3
      lb_Pi
    }
  )

  upper <- c(
    rep(if (stationary) limit else phi_bound, p), # strict root domain
    rep(if (invertible) limit else theta_bound, q), # strict root domain
    rep(Inf,  n_sigma),      # sigma2_t
    rep(Inf,  n_Psi),        # Psi_lower
    rep(Inf,  n_F2),         # A2
    rep(Inf,  n_Pi)          # Pi_pq
  )

  list(lower = lower, upper = upper,
       phi_bound = phi_bound, theta_bound = theta_bound,
       stationary = stationary, invertible = invertible)
}

.bound_hit_diagnostics <- function(psi, box, tol = 1e-8) {
  lower_hits <- upper_hits <- rep(FALSE, length(psi))
  finite_lower <- is.finite(box$lower)
  finite_upper <- is.finite(box$upper)
  lower_hits[finite_lower] <-
    psi[finite_lower] - box$lower[finite_lower] <=
    tol * pmax(1, abs(psi[finite_lower]), abs(box$lower[finite_lower]))
  upper_hits[finite_upper] <-
    box$upper[finite_upper] - psi[finite_upper] <=
    tol * pmax(1, abs(psi[finite_upper]), abs(box$upper[finite_upper]))
  list(lower = names(psi)[lower_hits],
       upper = names(psi)[upper_hits],
       any = any(lower_hits | upper_hits))
}

# ============================================================================ #
#  SECTION 3.  STARTING VALUE CONSTRUCTION
# ============================================================================ #

#' Construct sensible starting values for optimization
#'
#' Conditional strategy:
#'   phi <- supplied value, or pooled OLS (stabilized only if stationary=TRUE);
#'   F, Psi and Sigma_xi <- the sparse-plus-low-rank decomposition of
#'             B(phi) S_N B(phi)';
#'   theta, sigma2 and Pi_(p,q) <- recovered from the identified blocks of
#'             Sigma_xi.  No independent theta starting value is used.
#'
#' @param Y           N×T data matrix.
#' @param p,q,r       Model orders.
#' @param big_T       T.
#' @param heteroskedastic Logical.
#' @param phi_start   Optional AR coefficient vector. If NULL, use pooled OLS.
#' @param theta_bound Physical MA limit for start recovery when invertible=FALSE (default 5).
#' @return Named list compatible with pack_params.

# ---------------------------------------------------------------------------- #
#  Helper: check the covariance restrictions in Assumption 8
# ---------------------------------------------------------------------------- #

.is_uniformly_pd <- function(M, min_eig = 1e-8) {
  if (is.null(M) || !is.matrix(M) || nrow(M) != ncol(M)) return(FALSE)
  if (nrow(M) == 0L) return(TRUE)
  if (any(!is.finite(M))) return(FALSE)
  M <- (M + t(M)) / 2
  !is.null(tryCatch(chol(M - diag(min_eig, nrow(M))),
                    error = function(e) NULL))
}

.build_primitive_covariance <- function(vartheta, p, q, big_T,
                                        heteroskedastic = TRUE) {
  n_eps <- big_T - p + q
  n_Pi  <- p * (p + 1L) / 2L + compute_K_pq(p, q)
  if (n_eps <= 0L || length(vartheta$Pi_pq) != n_Pi) return(NULL)

  sigma2_eps <- if (heteroskedastic) {
    if (length(vartheta$sigma2) != n_eps) return(NULL)
    vartheta$sigma2
  } else {
    if (length(vartheta$sigma2) < 1L) return(NULL)
    rep(vartheta$sigma2[1L], n_eps)
  }
  if (any(!is.finite(sigma2_eps))) return(NULL)

  nu <- tryCatch(unpack_Pi(vartheta$Pi_pq, p, q), error = function(e) NULL)
  if (is.null(nu)) return(NULL)

  if (any(!is.finite(nu$cross_v))) return(NULL)
  A <- .build_A_pq(nu$cross_v, p, q)
  v0 <- if (q > 0L) sigma2_eps[seq_len(q)] else numeric(0)
  Schur <- nu$Sigma_eta
  if (p > 0L && q > 0L) {
    if (any(v0 <= 0)) return(NULL)
    Schur <- Schur - tcrossprod(sweep(A, 2L, sqrt(v0), "/"))
  }

  list(Schur = Schur, Sigma_eta = nu$Sigma_eta,
       A = A, v0 = v0, sigma2_eps = sigma2_eps)
}

.check_pd_components <- function(vartheta, p, q, r, big_T,
                                 heteroskedastic = TRUE, min_eig = 1e-8,
                                 Sigma_override = NULL) {
  finite_parameters <- all(is.finite(c(vartheta$phi, vartheta$theta,
                                        vartheta$sigma2, vartheta$Psi_lower,
                                        vartheta$F2, vartheta$Pi_pq)))

  primitive <- tryCatch(
    .build_primitive_covariance(vartheta, p, q, big_T, heteroskedastic),
    error = function(e) NULL
  )
  valid_variances <- !is.null(primitive) &&
    all(primitive$sigma2_eps > min_eig)
  pd_eta <- p == 0L || (!is.null(primitive) &&
    .is_uniformly_pd(primitive$Sigma_eta, 0))

  if (r == 0L) {
    pd_Psi <- TRUE
  } else {
    diag_pos <- .chol_diag_positions(r)
    Psi <- tryCatch(.Psi_from_lower(vartheta$Psi_lower, r),
                    error = function(e) NULL)
    pd_Psi <- length(vartheta$Psi_lower) == r * (r + 1L) / 2L &&
      all(vartheta$Psi_lower[diag_pos] > 0) &&
      .is_uniformly_pd(Psi, min_eig)
  }

  Sigma_u <- tryCatch(
    build_Sigma_u(vartheta$theta, vartheta$sigma2, p, big_T,
                       heteroskedastic),
    error = function(e) NULL
  )
  pd_u <- !is.null(Sigma_u) && .is_uniformly_pd(Sigma_u, 0)

  Sigma_xi <- tryCatch(
    build_Sigma_xi(vartheta$phi, vartheta$theta, vartheta$sigma2,
                   vartheta$Pi_pq, big_T, heteroskedastic),
    error = function(e) NULL
  )
  pd_xi <- !is.null(Sigma_xi) && .is_uniformly_pd(Sigma_xi, min_eig)

  Sigma <- if (!is.null(Sigma_override)) Sigma_override else tryCatch(
    build_Sigma_model(vartheta, p, q, r, big_T, heteroskedastic),
    error = function(e) NULL
  )
  pd_Sigma <- !is.null(Sigma) && .is_uniformly_pd(Sigma, min_eig)

  all_valid <- finite_parameters && valid_variances && pd_Psi &&
    pd_eta && pd_u && pd_xi && pd_Sigma

  list(innovation_variances = valid_variances,
       Psi = pd_Psi,
       Sigma_eta = pd_eta,
       Sigma_u = pd_u,
       Sigma_xi = pd_xi,
       Sigma = pd_Sigma,
       all_valid = all_valid,
       all_pd = all_valid)
}

# Make the paper-coordinate Pi_(p,q) satisfy primitive covariance feasibility.
# This helper is needed by the starting-value construction below.
make_Pi_admissible <- function(Pi_pq, p, min_var = 1e-8, q = 0L,
                               sigma2 = NULL, big_T = NULL,
                               heteroskedastic = TRUE) {
  if (p == 0L || length(Pi_pq) == 0L) return(Pi_pq)
  expected <- p * (p + 1L) / 2L + compute_K_pq(p, q)
  if (length(Pi_pq) != expected)
    stop("Pi_pq has the wrong length for the supplied p and q.")
  if (q > 0L && (is.null(sigma2) || is.null(big_T)))
    stop("For q>0, supply sigma2 and big_T to check primitive feasibility.")

  nu <- unpack_Pi(Pi_pq, p, q)
  Schur <- nu$Sigma_eta
  if (q > 0L) {
    v0 <- .early_innovation_variances(
      sigma2, p, q, big_T, heteroskedastic
    )
    if (any(!is.finite(v0)) || any(v0 <= 0))
      stop("Innovation variances must be positive.")
    A <- .build_A_pq(nu$cross_v, p, q)
    Schur <- Schur - tcrossprod(sweep(A, 2L, sqrt(v0), "/"))
  }
  Schur <- (Schur + t(Schur)) / 2
  if (any(!is.finite(Schur))) stop("Pi_pq contains nonfinite values.")
  min_eig <- min(eigen(Schur, symmetric = TRUE, only.values = TRUE)$values)
  shift <- max(0, 2 * min_var - min_eig)
  if (shift > 0) nu$Sigma_eta <- nu$Sigma_eta + diag(shift, p)
  c(.pack_Sigma_eta(nu$Sigma_eta), nu$cross_v)
}


.pooled_ols_phi <- function(Y, p, big_T) {
  if (p == 0L) return(numeric(0))
  N <- nrow(Y)
  lhs <- as.vector(Y[, (p + 1L):big_T, drop = FALSE])
  rhs <- matrix(0, N * (big_T - p), p)
  for (s in seq_len(p))
    rhs[, s] <- as.vector(Y[, (p + 1L - s):(big_T - s), drop = FALSE])
  phi <- as.numeric(lm.fit(rhs, lhs)$coefficients)
  phi[!is.finite(phi)] <- 0
  phi
}

.is_stationary_start <- function(phi, tol = 1e-6) {
  if (length(phi) == 0L) return(TRUE)
  roots <- polyroot(c(1, -phi))
  length(roots) == 0L || all(is.finite(roots) & Mod(roots) > 1 + tol)
}

.stabilize_phi_start <- function(phi, tol = 1e-6) {
  if (length(phi) == 0L || .is_stationary_start(phi, tol)) return(phi)
  scale <- 1
  while (scale > 1e-12 && !.is_stationary_start(scale * phi, tol))
    scale <- 0.95 * scale
  if (scale <= 1e-12) return(rep(0, length(phi)))

  lower <- scale
  upper <- min(1, scale / 0.95)
  for (iter in seq_len(60L)) {
    middle <- (lower + upper) / 2
    if (.is_stationary_start(middle * phi, tol)) lower <- middle else upper <- middle
  }
  lower * (1 - 1e-8) * phi
}

.pacf_to_ar <- function(kappa) {
  p <- length(kappa)
  if (p == 0L) return(numeric(0))
  phi <- kappa[1L]
  if (p > 1L) for (m in 2L:p) {
    old <- phi
    phi <- c(old - kappa[m] * rev(old), kappa[m])
  }
  as.numeric(phi)
}


# Project generated OLS/spectral starts in reflection space, then rebuild their
# physical coefficients before recovering the conditional nuisance parameters.
.project_ar_start <- function(phi) {
  phi <- .stabilize_phi_start(phi)
  kappa <- tryCatch(.ar_to_pacf(phi),
                    error = function(e) .ar_to_bounded_pacf(phi))
  limit <- .dynamic_reflection_limit
  .pacf_to_ar(pmax(-limit, pmin(limit, kappa)))
}

.stationary_phi_design <- function(n, p, kappa_limit = 0.99) {
  if (n <= 0L || p == 0L) return(matrix(numeric(0), n, p))
  U <- matrix(0, n, p)
  for (j in seq_len(p)) U[, j] <- (sample.int(n) - runif(n)) / n
  K <- kappa_limit * (2 * U - 1)
  out <- matrix(0, n, p)
  for (i in seq_len(n)) out[i, ] <- .pacf_to_ar(K[i, ])
  out
}

# Numerical optimization controls shared by all estimators:
# n_starts = 20, main factr = 1e4, ndeps = 1e-5 (in optimization units),
# start_factr = 1e4 for conditional MA initialization. The finite-difference
# step is used in main L-BFGS-B, MA starts, and factor-completion BFGS.
# nlminb fallback tolerances and inference derivative settings are separate.
.fit_heteroskedastic_ma_start <- function(noise, p, q, big_T, theta,
                                          min_var, theta_bound = 5,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  tail_index <- (p + 1L):big_T
  target <- noise[tail_index, tail_index, drop = FALSE]
  band <- lower.tri(target, diag = TRUE) &
    row(target) - col(target) <= q
  target_band <- target[band]
  variance_floor <- max(1e-6, min_var,
                        1e-4 * median(abs(diag(target))))
  n_sigma <- big_T - p + q
  spectral <- tryCatch(
    .recover_ma_homoskedastic(noise, p, q, big_T, min_var),
    error = function(e) NULL
  )
  if (!is.null(spectral)) theta <- spectral$theta
  theta <- if (invertible) -.project_ar_start(-theta) else
    pmax(-theta_bound, pmin(theta_bound, theta))
  initial <- if (invertible) .ar_to_pacf(-theta) else theta
  dynamic_bound <- if (invertible) .dynamic_reflection_limit else theta_bound

  profile_sigma2 <- function(dynamic_candidate, return_parameters = FALSE) {
    theta_candidate <- if (invertible) -.pacf_to_ar(dynamic_candidate) else
      dynamic_candidate
    X <- matrix(0, length(target_band), n_sigma)
    for (k in seq_len(n_sigma)) {
      unit <- numeric(n_sigma); unit[k] <- 1
      X[, k] <- build_Sigma_u(
        theta_candidate, unit, p, big_T, TRUE
      )[band]
    }
    sigma2 <- tryCatch(
      as.numeric(qr.solve(X, target_band, tol = 1e-10)),
      error = function(e) NULL
    )
    if (is.null(sigma2) || any(!is.finite(sigma2)))
      return(if (return_parameters) NULL else 1e15)
    sigma2 <- pmax(sigma2, variance_floor)
    value <- mean((target_band - X %*% sigma2)^2)
    if (return_parameters)
      list(theta = theta_candidate, sigma2 = sigma2, value = value) else value
  }
  fit <- optim(
    initial, profile_sigma2, method = "L-BFGS-B",
    lower = rep(-dynamic_bound, q), upper = rep(dynamic_bound, q),
    control = list(maxit = 1000L, factr = start_factr,
                   ndeps = rep(ndeps, q))
  )
  result <- profile_sigma2(fit$par, return_parameters = TRUE)
  if (is.null(result) || !is.finite(result$value) ||
      any(!is.finite(c(result$theta, result$sigma2))) ||
      any(result$sigma2 <= 0) || abs(result$theta[q]) <= 1e-8)
    stop("The finite-sample MA covariance projection failed.")
  result[c("theta", "sigma2")]
}

.recover_ma_heteroskedastic <- function(noise, p, q, big_T,
                                        min_var = 1e-8,
                                        theta_bound = 5,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  a <- p + 1L
  if (q == 0L) {
    sigma2 <- diag(noise)[a:big_T]
    if (any(!is.finite(sigma2)) || any(sigma2 <= min_var))
      stop("The recovered innovation variances are not positive.")
    return(list(theta = numeric(0), sigma2 = sigma2))
  }

  m_q <- as.integer(ceiling(q / 2))
  if (big_T < p + q + m_q + 2L)
    stop("T is too small for the conditional MA recovery.")

  g <- rep(NA_real_, big_T)
  t_outer <- a:(big_T - q)
  g[t_outer] <- noise[cbind(t_outer, t_outer + q)]
  theta_low <- numeric(m_q)
  rho <- numeric(m_q)

  for (k in seq_len(m_q)) {
    tt <- (a + k):(big_T - q)
    R_kt <- noise[cbind(tt, tt + q - k)]
    if (k > 1L) for (j in seq_len(k - 1L))
      R_kt <- R_kt - theta_low[j] * rho[k - j] * g[tt - j]
    X <- cbind(g[tt], g[tt - k])
    if (any(!is.finite(c(X, R_kt))) || qr(X, tol = 1e-10)$rank < 2L)
      stop("The conditional MA recovery system is singular.")
    coef_k <- as.numeric(qr.solve(X, R_kt, tol = 1e-10))
    rho[k] <- coef_k[1L]
    theta_low[k] <- coef_k[2L]
  }

  base_index <- q - m_q
  base_theta <- if (base_index == 0L) 1 else theta_low[base_index]
  if (!is.finite(rho[m_q]) || abs(rho[m_q]) <= 1e-8)
    stop("The recovered MA scale is numerically unidentified.")
  theta <- numeric(q)
  theta[seq_len(m_q)] <- theta_low
  theta[q] <- base_theta / rho[m_q]
  for (k in seq_len(m_q)) {
    reflected <- q - k
    if (reflected > m_q && reflected < q)
      theta[reflected] <- rho[k] * theta[q]
  }
  if (any(!is.finite(theta)) || any(abs(theta) <= 1e-8))
    stop("The recovered MA coefficients are invalid or numerically zero.")

  .fit_heteroskedastic_ma_start(
    noise, p, q, big_T, theta, min_var, theta_bound
  ,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      )
}

.recover_ma_homoskedastic <- function(noise, p, q, big_T,
                                      min_var = 1e-8) {
  a <- p + 1L
  gamma <- vapply(0:q, function(h) {
    tt <- a:(big_T - h)
    mean(noise[cbind(tt, tt + h)])
  }, numeric(1L))
  if (any(!is.finite(gamma)) || gamma[1L] <= min_var)
    stop("The recovered MA autocovariances are invalid.")
  if (q == 0L) return(list(theta = numeric(0), sigma2 = gamma[1L]))

  exponents <- 0:(2L * q)
  spectral_poly <- gamma[abs(exponents - q) + 1L]
  roots <- polyroot(spectral_poly)
  if (length(roots) != 2L * q || any(!is.finite(roots)))
    stop("The homoskedastic MA spectral factorization failed.")
  ord <- order(Mod(roots), decreasing = TRUE)
  selected <- roots[ord[seq_len(q)]]
  if (any(Mod(selected) <= 1 + 1e-8))
    stop("The homoskedastic MA covariance has no separated invertible factor.")

  coefs <- 1 + 0i
  for (root in selected) coefs <- c(coefs, 0) + c(0, -coefs / root)
  if (max(abs(Im(coefs))) > 1e-6 * max(1, max(abs(Re(coefs)))))
    stop("The homoskedastic MA factor is not numerically real.")
  theta <- Re(coefs[-1L])
  sigma2 <- gamma[1L] / sum(c(1, theta)^2)
  if (any(!is.finite(theta)) || !is.finite(sigma2) || sigma2 <= min_var)
    stop("The recovered homoskedastic MA parameters are invalid.")
  list(theta = theta, sigma2 = sigma2)
}

.recover_ma_start <- function(noise, p, q, big_T, heteroskedastic,
                              min_var = 1e-8, theta_bound = 5,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  if (heteroskedastic)
    .recover_ma_heteroskedastic(noise, p, q, big_T, min_var,
                                theta_bound,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      ) else
    .recover_ma_homoskedastic(noise, p, q, big_T, min_var)
}

make_start_values <- function(Y, p, q, r, big_T, heteroskedastic = TRUE,
                              phi_start = NULL, theta_bound = 5,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  theta_bound <- .validate_coefficient_bound(theta_bound, "theta_bound")
  stationary <- .validate_dynamic_flag(stationary, "stationary")
  invertible <- .validate_dynamic_flag(invertible, "invertible")
  N <- nrow(Y)
  phi0 <- if (is.null(phi_start)) {
    pooled <- .pooled_ols_phi(Y, p, big_T)
    if (stationary) .project_ar_start(pooled) else pooled
  } else {
    as.numeric(phi_start)
  }
  if (length(phi0) != p || any(!is.finite(phi0)))
    stop("phi_start must be a finite vector of length p.")

  if (stationary) .encode_dynamics(phi0, numeric(0), TRUE, FALSE)

  Y_c <- sweep(Y, 2L, colMeans(Y))
  start <- tryCatch({
    S_N <- crossprod(Y_c) / N
    B0 <- build_B_p(phi0, big_T)
    inner <- B0 %*% S_N %*% t(B0)
    inner <- (inner + t(inner)) / 2

    if (r > 0L) {
      ee <- eigen(inner, symmetric = TRUE)
      noise_level <- if (r < big_T)
        median(pmax(ee$values[(r + 1L):big_T], 0)) else 0
      A <- ee$vectors[, seq_len(r), drop = FALSE] %*%
        diag(sqrt(pmax(ee$values[seq_len(r)] - noise_level, 1e-8)), r)
      H <- A[seq_len(r), , drop = FALSE]
      if (qr(H)$rank < r) stop("singular PCA normalization")
      A <- A %*% solve(H, t(chol(H %*% t(H))))

      pure_factor <- matrix(FALSE, big_T, big_T)
      if (big_T > 1L) for (ii in 2L:big_T) for (jj in seq_len(ii - 1L)) {
        pure_factor[ii, jj] <-
          (ii > p && jj > p && ii - jj > q) ||
          (ii > p && jj <= p && ii > jj + q)
      }
      if (any(pure_factor)) {
        target <- inner[pure_factor]
        completion <- tryCatch(
          optim(
            as.vector(A),
            function(a) {
              AA <- tcrossprod(matrix(a, big_T, r))
              mean((target - AA[pure_factor])^2)
            },
            method = "BFGS", control = list(maxit = 1000L, reltol = 1e-10,
                                             ndeps = rep(ndeps, length(A)))
          ),
          error = function(e) NULL
        )
        if (!is.null(completion) && is.finite(completion$value)) {
          A_try <- matrix(completion$par, big_T, r)
          if (qr(A_try[seq_len(r), , drop = FALSE])$rank == r) A <- A_try
        }
      }
      H <- A[seq_len(r), , drop = FALSE]
      L_Psi <- t(chol(H %*% t(H)))
      A <- A %*% solve(H, L_Psi)
      F2 <- if (big_T > r)
        as.vector(A[(r + 1L):big_T, , drop = FALSE] %*% solve(L_Psi)) else
        numeric(0)
      factor_part <- tcrossprod(A)
      Psi_lower <- L_Psi[lower.tri(L_Psi, diag = TRUE)]
    } else {
      factor_part <- matrix(0, big_T, big_T)
      Psi_lower <- numeric(0)
      F2 <- numeric(0)
    }

    noise <- (inner - factor_part + t(inner - factor_part)) / 2
    min_var <- max(1e-8, 1e-10 * mean(diag(inner)))
    ma <- .recover_ma_start(
      noise, p, q, big_T, heteroskedastic, min_var, theta_bound
    ,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      )
    theta0 <- ma$theta
    if (invertible) theta0 <- -.project_ar_start(-theta0)
    sigma2 <- ma$sigma2

    K_pq <- compute_K_pq(p, q)
    if (p > 0L) {
      eta_start_floor <- max(min_var, 1e-6)
      Sigma_eta <- noise[seq_len(p), seq_len(p), drop = FALSE]
      Sigma_eta <- (Sigma_eta + t(Sigma_eta)) / 2
      min_eta <- min(eigen(Sigma_eta, symmetric = TRUE,
                           only.values = TRUE)$values)
      if (min_eta < 2 * eta_start_floor)
        Sigma_eta <- Sigma_eta + diag(2 * eta_start_floor - min_eta, p)

      cross <- numeric(K_pq)
      if (K_pq > 0L) {
        X <- matrix(0, p * (big_T - p), K_pq)
        for (k in seq_len(K_pq)) {
          unit <- numeric(K_pq); unit[k] <- 1
          Pi_unit <- c(rep(1, p), rep(0, p * (p - 1L) / 2L), unit)
          Sx <- build_Sigma_xi(
            phi0, theta0, sigma2, Pi_unit, big_T, heteroskedastic
          )
          X[, k] <- as.vector(
            Sx[seq_len(p), (p + 1L):big_T, drop = FALSE]
          )
        }
        if (qr(X, tol = 1e-10)$rank < K_pq)
          stop("The projected-error bridge system is singular.")
        cross <- as.numeric(qr.solve(
          X,
          as.vector(noise[seq_len(p), (p + 1L):big_T, drop = FALSE]),
          tol = 1e-10
        ))
        if (any(!is.finite(cross)))
          stop("The projected-error bridge recovery failed.")
      }
      Pi_pq <- make_Pi_admissible(
        c(.pack_Sigma_eta(Sigma_eta), cross), p,
        min_var = eta_start_floor,
        q = q, sigma2 = sigma2, big_T = big_T,
        heteroskedastic = heteroskedastic
      )
    } else {
      Pi_pq <- numeric(0)
    }

    list(phi = phi0, theta = theta0, sigma2 = sigma2,
         Psi_lower = Psi_lower, F2 = F2, Pi_pq = Pi_pq)
  }, error = function(e) {
    stop("Unable to recover a feasible start conditional on phi: ",
         conditionMessage(e), call. = FALSE)
  })
  start
}

.validate_n_starts <- function(n_starts) {
  if (length(n_starts) != 1L || !is.finite(n_starts) ||
      n_starts < 1L || n_starts != floor(n_starts)) {
    stop("n_starts must be an integer greater than or equal to 1.")
  }
  as.integer(n_starts)
}

.validate_inputs <- function(Y, p, q, r, heteroskedastic,
                             Pi_only = FALSE, warn_T = TRUE) {
  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) < 1L || ncol(Y) < 1L ||
      any(!is.finite(Y)))
    stop("Y must be a nonempty finite numeric matrix without missing values.")
  orders <- list(p = p, q = q, r = r)
  invalid <- vapply(orders, function(x)
    length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x < 0 || x != floor(x), logical(1))
  if (any(invalid)) stop("p, q, and r must be nonnegative integers.")

  big_T <- ncol(Y)
  if (p >= big_T) stop("p must be strictly smaller than T.")
  if (r >= big_T) stop("r must be strictly smaller than T.")
  minimum_T <- 2L * (r + p + q + 1L) - min(p, q)
  if (warn_T && big_T < minimum_T)
    warning(sprintf("T=%d may be insufficient for (p,q,r)=(%d,%d,%d); need T>=%d.",
                    big_T, as.integer(p), as.integer(q), as.integer(r),
                    as.integer(minimum_T)), call. = FALSE)

  d <- if (Pi_only) p * (p + 1L) / 2L + compute_K_pq(p, q) else
    count_params(p, q, r, big_T, heteroskedastic)
  d <- as.integer(d)
  dot_T <- as.integer(big_T * (big_T + 1L) / 2L)
  if (dot_T < d)
    stop(sprintf("The model has %d covariance moments but %d parameters.", dot_T, d))
  list(N = nrow(Y), big_T = big_T, dot_T = as.integer(dot_T))
}

.draw_full_start <- function(psi_start, box, p, q, r, big_T,
                             heteroskedastic, randomize = TRUE, spread = 1) {
  if (!randomize) return(pmax(box$lower, pmin(box$upper, psi_start)))

  n_phi <- p
  n_theta <- q
  n_sigma <- if (heteroskedastic) big_T - p + q else 1L
  n_Psi <- r * (r + 1L) / 2L
  n_F2 <- (big_T - r) * r
  n_eta <- p * (p + 1L) / 2L
  K_pq <- compute_K_pq(p, q)
  diag_pos_Psi <- .chol_diag_positions(r)
  diag_pos_eta <- .chol_diag_positions(p)

  psi_init <- psi_start
  pos <- 1L
  if (n_phi > 0L) {
    psi_init[pos:(pos + n_phi - 1L)] <-
      psi_start[pos:(pos + n_phi - 1L)] + rnorm(n_phi, 0, 0.2 * spread)
  }
  pos <- pos + n_phi
  if (n_theta > 0L) {
    psi_init[pos:(pos + n_theta - 1L)] <-
      psi_start[pos:(pos + n_theta - 1L)] + rnorm(n_theta, 0, 0.2 * spread)
  }
  pos <- pos + n_theta

  sigma_pos <- pos:(pos + n_sigma - 1L)
  psi_init[sigma_pos] <- pmax(
    box$lower[sigma_pos],
    psi_start[sigma_pos] *
      exp(rnorm(n_sigma, 0, 0.5 * spread))
  )
  pos <- pos + n_sigma

  if (n_Psi > 0L) {
    block <- psi_start[pos:(pos + n_Psi - 1L)]
    diag_pos_full <- pos + diag_pos_Psi - 1L
    block[diag_pos_Psi] <- pmax(
      box$lower[diag_pos_full],
      psi_start[diag_pos_full] *
        exp(rnorm(length(diag_pos_Psi), 0, 0.4 * spread))
    )
    off_diag <- setdiff(seq_len(n_Psi), diag_pos_Psi)
    if (length(off_diag) > 0L) {
      block[off_diag] <- psi_start[pos + off_diag - 1L] +
        rnorm(length(off_diag), 0, 0.2 * spread)
    }
    psi_init[pos:(pos + n_Psi - 1L)] <- block
  }
  pos <- pos + n_Psi

  if (n_F2 > 0L) {
    psi_init[pos:(pos + n_F2 - 1L)] <-
      psi_start[pos:(pos + n_F2 - 1L)] +
        rnorm(n_F2, 0, 0.5 * spread)
  }
  pos <- pos + n_F2

  if (n_eta + K_pq > 0L) {
    block <- psi_start[pos:(pos + n_eta + K_pq - 1L)]
    if (length(diag_pos_eta) > 0L) {
      diag_pos_full <- pos + diag_pos_eta - 1L
      block[diag_pos_eta] <- pmax(
        box$lower[diag_pos_full],
        psi_start[diag_pos_full] *
          exp(rnorm(length(diag_pos_eta), 0, 0.4 * spread))
      )
    }
    remaining <- setdiff(seq_len(n_eta + K_pq), diag_pos_eta)
    if (length(remaining) > 0L) {
      block[remaining] <- psi_start[pos + remaining - 1L] +
        rnorm(length(remaining), 0, 0.2 * spread)
    }
    psi_init[pos:(pos + n_eta + K_pq - 1L)] <- block
  }

  pmax(box$lower, pmin(box$upper, psi_init))
}

.valid_initial_vector <- function(psi, box, p, q, r, big_T,
                                  heteroskedastic) {
  d <- count_params(p, q, r, big_T, heteroskedastic)
  if (is.null(psi) || length(psi) != d || any(!is.finite(psi)) ||
      any(psi < box$lower) || any(psi > box$upper)) return(FALSE)
  Sigma <- tryCatch(
    .build_Sigma_internal(psi, p, q, r, big_T, heteroskedastic,
      stationary = box$stationary, invertible = box$invertible),
    error = function(e) NULL
  )
  !is.null(Sigma) && all(is.finite(Sigma)) &&
    !is.null(tryCatch(chol(Sigma), error = function(e) NULL))
}

# Deterministic last-resort start used only after the ordinary optimization
# attempts for a start slot have failed.  In paper coordinates it sets
# phi=theta=F2=0, sigma2_t=1, Psi=I_r, Sigma_eta=I_p, and all cross
# covariances to zero.  The internal coordinates therefore contain identity
# Cholesky factors for Psi and the Schur complement of Sigma_eta.
.make_stable_start <- function(box, p, q, r, big_T, heteroskedastic) {
  n_sigma <- if (heteroskedastic) big_T - p + q else 1L
  L_Psi <- if (r > 0L) diag(r) else matrix(0, 0L, 0L)
  Psi_lower <- if (r > 0L) L_Psi[lower.tri(L_Psi, diag = TRUE)] else
    numeric(0)
  n_F2 <- (big_T - r) * r

  L_eta <- if (p > 0L) diag(p) else matrix(0, 0L, 0L)
  eta_lower <- if (p > 0L) L_eta[lower.tri(L_eta, diag = TRUE)] else
    numeric(0)
  cross <- numeric(compute_K_pq(p, q))

  psi <- c(
    rep(0, p),
    rep(0, q),
    rep(1, n_sigma),
    Psi_lower,
    rep(0, n_F2),
    eta_lower,
    cross
  )
  psi <- pmax(box$lower, pmin(box$upper, psi))
  if (!.valid_initial_vector(
    psi, box, p, q, r, big_T, heteroskedastic
  )) {
    stop("Internal error: the deterministic stable start is infeasible.",
         call. = FALSE)
  }
  list(psi = psi, phi = rep(0, p), type = "stable_deterministic")
}

.conditional_start_from_phi <- function(Y, phi, box, p, q, r, big_T,
                                        heteroskedastic,
                                  ndeps = 1e-5, start_factr = 1e4) {
  if (p > 0L && !box$stationary)
    phi <- pmax(-box$phi_bound, pmin(box$phi_bound, phi))
  sv <- tryCatch(
    make_start_values(
      Y, p, q, r, big_T, heteroskedastic, phi_start = phi,
      theta_bound = box$theta_bound
    ,
      stationary = box$stationary, invertible = box$invertible,
        ndeps = ndeps, start_factr = start_factr
      ),
    error = function(e) NULL
  )
  if (is.null(sv)) return(NULL)
  psi <- tryCatch(
    pack_params(sv, p, q, r, big_T, heteroskedastic,
      stationary = box$stationary, invertible = box$invertible),
    error = function(e) NULL
  )
  if (is.null(psi) || any(!is.finite(psi))) return(NULL)
  psi <- pmax(box$lower, pmin(box$upper, psi))
  if (!.valid_initial_vector(
    psi, box, p, q, r, big_T, heteroskedastic
  )) return(NULL)
  psi
}

.draw_random_conditional_start <- function(Y, box, p, q, r, big_T,
                                           heteroskedastic,
                                           max_attempts = 1000L,
                                  ndeps = 1e-5, start_factr = 1e4) {
  if (p == 0L) {
    psi <- .conditional_start_from_phi(
      Y, numeric(0), box, p, q, r, big_T, heteroskedastic,
      ndeps = ndeps, start_factr = start_factr
    )
    if (is.null(psi)) stop("Unable to construct the unique p=0 conditional start.")
    return(list(psi = psi, phi = numeric(0), attempts = 1L))
  }
  for (attempt in seq_len(max_attempts)) {
    phi <- rnorm(p, mean = 0, sd = 0.5)
    if (box$stationary && !.is_stationary_start(phi)) next
    psi <- .conditional_start_from_phi(
      Y, phi, box, p, q, r, big_T, heteroskedastic,
      ndeps = ndeps, start_factr = start_factr
    )
    if (!is.null(psi)) return(list(psi = psi, phi = phi, attempts = attempt))
  }
  stop(sprintf(
    "Unable to recover a valid conditional start after %d random phi draws.",
    max_attempts
  ))
}

.draw_valid_full_start <- function(anchor, box, p, q, r, big_T,
                                   heteroskedastic,
                                   max_attempts = 1000L) {
  for (attempt in seq_len(max_attempts)) {
    psi <- .draw_full_start(
      anchor, box, p, q, r, big_T, heteroskedastic,
      randomize = TRUE, spread = 3
    )
    ar <- if (p > 0L) psi[seq_len(p)] else numeric(0)
    phi <- if (box$stationary) .pacf_to_ar(ar) else ar
    if (.valid_initial_vector(
      psi, box, p, q, r, big_T, heteroskedastic
    )) return(list(psi = psi, phi = phi, attempts = attempt))
  }
  stop(sprintf(
    "Unable to draw a valid full-parameter start after %d attempts.",
    max_attempts
  ))
}

# Start 1 uses phi_start (or pooled OLS); other nonmultiples of five use a
# space-filling stationary-AR design; starts 5, 10, ... perturb every block.
# A failed conditional reconstruction is replaced by iid N(0, 0.5^2) AR draws.
.make_initial_starts <- function(Y, box, p, q, r, big_T,
                                 heteroskedastic, n_starts,
                                 phi_start = NULL,
                                 max_attempts = 1000L,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  if (!is.null(phi_start)) {
    phi_start <- as.numeric(phi_start)
    if (length(phi_start) != p || any(!is.finite(phi_start)))
      stop("phi_start must be a finite vector of length p.")
    if (box$stationary) .encode_dynamics(phi_start, numeric(0), TRUE, FALSE)
    if (!box$stationary && any(abs(phi_start) > box$phi_bound))
      stop("phi_start exceeds phi_bound.")
  }

  starts <- vector("list", n_starts)
  types <- character(n_starts)
  initial_phi <- vector("list", n_starts)

  primary_phi <- if (is.null(phi_start)) {
    pooled <- .pooled_ols_phi(Y, p, big_T)
    if (box$stationary) .project_ar_start(pooled) else
      pmax(-box$phi_bound, pmin(box$phi_bound, pooled))
  } else {
    phi_start
  }
  primary <- .conditional_start_from_phi(
    Y, primary_phi, box, p, q, r, big_T, heteroskedastic,
    ndeps = ndeps, start_factr = start_factr
  )
  if (is.null(primary)) {
    replacement <- tryCatch(
      .draw_random_conditional_start(
        Y, box, p, q, r, big_T, heteroskedastic, max_attempts,
        ndeps = ndeps, start_factr = start_factr
      ),
      error = function(e) NULL
    )
    if (is.null(replacement)) {
      replacement <- .make_stable_start(
        box, p, q, r, big_T, heteroskedastic
      )
    }
    primary <- replacement$psi
    primary_phi <- replacement$phi
    types[1L] <- if (is.null(replacement$type)) "conditional_random" else
      replacement$type
  } else {
    types[1L] <- if (is.null(phi_start)) "conditional_ols" else
      "conditional_phi_start"
  }
  starts[[1L]] <- primary
  initial_phi[[1L]] <- if (p > 0L) {
        ar <- primary[seq_len(p)]
        if (box$stationary) .pacf_to_ar(ar) else ar
      } else numeric(0)
  if (n_starts == 1L) {
    attr(starts, "start_type") <- types
    attr(starts, "initial_phi") <- initial_phi
    return(starts)
  }

  structured_slots <- which(seq_len(n_starts) > 1L &
    seq_len(n_starts) %% 5L != 0L)
  if (p == 0L) structured_slots <- integer(0)
  phi_grid <- .stationary_phi_design(length(structured_slots), p)
  grid_pos <- 0L

  for (s in 2L:n_starts) {
    full_start <- s %% 5L == 0L || p == 0L
    if (full_start) {
      draw <- tryCatch(
        .draw_valid_full_start(
          primary, box, p, q, r, big_T, heteroskedastic, max_attempts
        ),
        error = function(e) NULL
      )
      if (is.null(draw)) {
        draw <- .make_stable_start(box, p, q, r, big_T, heteroskedastic)
      }
      starts[[s]] <- draw$psi
      initial_phi[[s]] <- if (p > 0L) {
        ar <- draw$psi[seq_len(p)]
        if (box$stationary) .pacf_to_ar(ar) else ar
      } else numeric(0)
      types[s] <- if (is.null(draw$type)) "full_random" else draw$type
    } else {
      grid_pos <- grid_pos + 1L
      phi <- phi_grid[grid_pos, ]
      psi <- .conditional_start_from_phi(
        Y, phi, box, p, q, r, big_T, heteroskedastic,
        ndeps = ndeps, start_factr = start_factr
      )
      if (is.null(psi)) {
        draw <- tryCatch(
          .draw_random_conditional_start(
            Y, box, p, q, r, big_T, heteroskedastic, max_attempts,
            ndeps = ndeps, start_factr = start_factr
          ),
          error = function(e) NULL
        )
        if (is.null(draw)) {
          draw <- .make_stable_start(box, p, q, r, big_T, heteroskedastic)
        }
        psi <- draw$psi
        phi <- draw$phi
        types[s] <- if (is.null(draw$type)) "conditional_random" else draw$type
      } else {
        types[s] <- "conditional_grid"
      }
      starts[[s]] <- psi
      initial_phi[[s]] <- if (p > 0L) {
        ar <- psi[seq_len(p)]
        if (box$stationary) .pacf_to_ar(ar) else ar
      } else numeric(0)
    }
  }
  attr(starts, "start_type") <- types
  attr(starts, "initial_phi") <- initial_phi
  starts
}

.draw_replacement_start <- function(start, anchor, Y, box, p, q, r, big_T,
                                    heteroskedastic,
                                    max_attempts = 1000L,
                                  ndeps = 1e-5, start_factr = 1e4) {
  if (start %% 5L == 0L || p == 0L) {
    draw <- .draw_valid_full_start(
      anchor, box, p, q, r, big_T, heteroskedastic, max_attempts
    )
    draw$type <- "full_random"
  } else {
    draw <- .draw_random_conditional_start(
      Y, box, p, q, r, big_T, heteroskedastic, max_attempts,
      ndeps = ndeps, start_factr = start_factr
    )
    draw$type <- "conditional_random"
  }
  draw
}

.run_box_optim <- function(par, fn, lower, upper, maxit, factr,
                                  ndeps = 1e-5) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  res <- tryCatch(
    optim(par, fn, method = "L-BFGS-B", lower = lower, upper = upper,
          control = list(maxit = maxit, factr = factr,
                          ndeps = rep(ndeps, length(par)))),
    error = function(e) list(value = Inf, convergence = 999L, par = par,
                              message = conditionMessage(e))
  )
  if (!is.finite(res$value) || res$convergence != 0L) {
    nb <- tryCatch(
      nlminb(res$par, fn, lower = lower, upper = upper,
             control = list(iter.max = maxit, eval.max = 5L * maxit)),
      error = function(e) NULL
    )
    if (!is.null(nb) && is.finite(nb$objective) && nb$convergence == 0L)
      res <- list(par = nb$par, value = nb$objective, convergence = 0L,
                  message = nb$message)
  }
  res
}


# ============================================================================ #
#  SECTION 5.  QMLE OBJECTIVE  [Paper eq. 20]
# ============================================================================ #

# Helper: vech (lower-triangular stack, column-major)
vech <- function(M) M[lower.tri(M, diag = TRUE)]

.to_paper_parameters <- function(psi, p, q, r, big_T,
                                 heteroskedastic = TRUE,
                                  stationary = TRUE, invertible = TRUE) {
  vartheta <- unpack_params(psi, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
  Psi <- .Psi_from_lower(vartheta$Psi_lower, r)
  c(vartheta$phi, vartheta$theta, vartheta$sigma2, vech(Psi),
    vartheta$F2, vartheta$Pi_pq)
}

.se_from_variance <- function(V, N, roundoff_tol = 1e-8) {
  v <- diag(V) / N
  finite_v <- is.finite(v)
  finite_V <- is.finite(V)
  variance_scale <- if (any(finite_V)) {
    max(.Machine$double.xmin, max(abs(V[finite_V])) / N)
  } else {
    1
  }
  roundoff_negative <- finite_v & v < 0 &
    abs(v) <= roundoff_tol * variance_scale
  v[roundoff_negative] <- 0
  sqrt(ifelse(is.finite(v) & v >= 0, v, NA_real_))
}

.set_empty_inference <- function(result) {
  d0 <- length(result$vartheta_hat)
  pnames <- names(result$vartheta_hat)
  result$se <- setNames(rep(NA_real_, d0), pnames)
  result$ci <- matrix(NA_real_, 2, d0,
                      dimnames = list(c("lower_95", "upper_95"), pnames))
  result$V_hat <- matrix(NA_real_, d0, d0, dimnames = list(pnames, pnames))
  result
}

.failed_optimizer_result <- function(estimator, p, q, r, big_T, N,
                                     heteroskedastic, all_runs, obj_val,
                                  stationary = TRUE, invertible = TRUE,
                                     phi_bound = 5, theta_bound = 5) {
  d0 <- count_params(p, q, r, big_T, heteroskedastic)
  internal_names <- .parameter_names(
    p, q, r, big_T, heteroskedastic, internal = TRUE
  ,
      stationary = stationary, invertible = invertible)
  paper_names <- .parameter_names(
    p, q, r, big_T, heteroskedastic, internal = FALSE
  ,
      stationary = stationary, invertible = invertible)
  n_sigma <- if (heteroskedastic) big_T - p + q else 1L
  n_Psi <- r * (r + 1L) / 2L
  n_F2 <- (big_T - r) * r
  n_Pi <- p * (p + 1L) / 2L + compute_K_pq(p, q)

  psi_hat <- setNames(rep(NA_real_, d0), internal_names)
  vartheta_hat <- setNames(rep(NA_real_, d0), paper_names)
  coef <- list(
    phi = rep(NA_real_, p),
    theta = rep(NA_real_, q),
    sigma2 = rep(NA_real_, n_sigma),
    Psi_lower = rep(NA_real_, n_Psi),
    Psi = matrix(NA_real_, r, r),
    F2 = rep(NA_real_, n_F2),
    A2 = rep(NA_real_, n_F2),
    Pi_pq = rep(NA_real_, n_Pi),
    Pi_internal = rep(NA_real_, n_Pi)
  )
  pd_components <- list(
    innovation_variances = FALSE,
    Psi = FALSE,
    Sigma_eta = FALSE,
    Sigma_u = FALSE,
    Sigma_xi = FALSE,
    Sigma = FALSE,
    all_valid = FALSE,
    all_pd = FALSE
  )
  successful <- sum(vapply(
    all_runs,
    function(x) is.list(x) && isTRUE(x$success),
    logical(1L)
  ))
  result <- list(
    coef = coef,
    psi_hat = psi_hat,
    vartheta_hat = vartheta_hat,
    estimator = estimator,
    obj_val = obj_val,
    optimization_obj = Inf,
    converged = FALSE,
    convergence_code = 999L,
    pd_components = pd_components,
    Sigma_hat = matrix(NA_real_, big_T, big_T),
    all_runs = all_runs,
    successful_starts = successful,
    failed_starts = length(all_runs) - successful,
    failure_reason = "No start produced a valid converged optimization.",
    hansen_test = NULL,
    p = p, q = q, r = r, big_T = big_T, N = N
  )
  result <- .set_dynamic_domain(result, stationary, invertible, phi_bound, theta_bound)
  class(result) <- c("arma_panel_fit", "list")
  .set_empty_inference(result)
}

.set_reported_inference <- function(result, V_hat, N) {
  pnames <- names(result$vartheta_hat)
  V_hat <- (V_hat + t(V_hat)) / 2
  dimnames(V_hat) <- list(pnames, pnames)
  se <- setNames(.se_from_variance(V_hat, N), pnames)
  ci <- rbind(result$vartheta_hat - 1.96 * se,
              result$vartheta_hat + 1.96 * se)
  rownames(ci) <- c("lower_95", "upper_95")

  result$V_hat <- V_hat
  result$se <- se
  result$ci <- ci
  result
}


# Record dynamic coordinates separately from paper-coordinate estimates/inference.
.set_dynamic_domain <- function(result, stationary, invertible,
                                phi_bound = 5, theta_bound = 5) {
  result$stationary <- stationary
  result$invertible <- invertible
  result$phi_bound <- phi_bound
  result$theta_bound <- theta_bound
  result$reflection_limit <- .dynamic_reflection_limit
  result$parameterization$version <- "dynamic-reflections-v1"
  result$parameterization$optimization_dynamic <- c(
    AR = if (stationary) "reflection coefficients [-(1-1e-4),1-1e-4]" else "physical phi",
    MA = if (invertible) "reflection coefficients [-(1-1e-4),1-1e-4]" else "physical theta"
  )
  ar <- if (result$p > 0L) unname(result$psi_hat[seq_len(result$p)]) else numeric(0)
  ma <- if (result$q > 0L) unname(result$psi_hat[result$p + seq_len(result$q)]) else numeric(0)
  ar_boundary <- if (stationary) {
    if (any(!is.finite(ar))) NA else any(abs(ar) >= .dynamic_reflection_limit - 1e-8)
  } else FALSE
  ma_boundary <- if (invertible) {
    if (any(!is.finite(ma))) NA else any(abs(ma) >= .dynamic_reflection_limit - 1e-8)
  } else FALSE
  min_root <- function(coefficients, sign) {
    if (any(!is.finite(coefficients))) return(NA_real_)
    roots <- polyroot(c(1, sign * coefficients))
    if (length(roots) == 0L) Inf else min(Mod(roots))
  }
  result$dynamic_diagnostics <- list(
    kappa_ar = if (stationary) ar else NULL,
    kappa_ma = if (invertible) ma else NULL,
    ar_reflection_boundary = ar_boundary,
    ma_reflection_boundary = ma_boundary,
    reflection_boundary = ar_boundary | ma_boundary,
    ar_min_root_modulus = min_root(result$coef$phi, -1),
    ma_min_root_modulus = min_root(result$coef$theta, 1),
    boundary_inference_note = if (isTRUE(ar_boundary | ma_boundary))
      "A reflection bound is active or near-active; ordinary interior Wald inference is not assured." else NULL
  )
  result
}

.warn_dynamic_roots <- function(vartheta, estimator, tol = 1e-8) {
  if (length(vartheta$phi) > 0L && all(is.finite(vartheta$phi))) {
    ar_roots <- polyroot(c(1, -vartheta$phi))
    if (length(ar_roots) > 0L && any(Mod(ar_roots) <= 1 + tol)) {
      warning(sprintf(
        "The selected %s estimate has an AR root on, near, or inside the unit circle",
        estimator
      ), call. = FALSE)
    }
  }
  if (length(vartheta$theta) > 0L && all(is.finite(vartheta$theta))) {
    ma_roots <- polyroot(c(1, vartheta$theta))
    if (length(ma_roots) > 0L && any(Mod(ma_roots) <= 1 + tol)) {
      warning(sprintf(
        "The selected %s estimate has an MA root on, near, or inside the unit circle; ordinary interior Wald inference may not apply.",
        estimator
      ), call. = FALSE)
    }
  }
  invisible(NULL)
}

# Helper: evaluate concentrated quasi-log-likelihood at psi
# Returns the NEGATIVE QLL (for minimization by optim).
# Returns 1e15 (large penalty) if Sigma is not PD or if any error occurs.
.qll_neg <- function(psi, p, q, r, big_T, S_N, heteroskedastic,
                                  stationary = TRUE, invertible = TRUE) {
  Sigma <- tryCatch(
    .build_Sigma_internal(psi, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible),
    error = function(e) NULL
  )
  if (is.null(Sigma)) return(1e15)

  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if (is.null(R)) return(1e15)

  # QLL (eq. 20, concentrated out delta): -ln|Sigma| - tr[S_N Sigma^{-1}]
  logdet <- 2 * sum(log(diag(R)))
  Sigma_inv_S <- backsolve(R, forwardsolve(t(R), S_N))
  objective <- logdet + sum(diag(Sigma_inv_S))
  if (!is.finite(objective)) return(1e15)
  objective
}

# ============================================================================ #
#  SECTION 6.  QMLE ESTIMATOR  [Paper Section 3.1]
# ============================================================================ #

# QMLE uses analytical covariance and mean derivatives.
# Derivatives and sandwiches use stable factor coordinates. Reporting maps
# covariance back to paper parameters; optimization remains unchanged.
# hessian_type selects the QMLE sandwich bread.  With "both", the observed
# Hessian remains primary and the information-Hessian inference is also returned.
# stationary/invertible default TRUE and use reflection boxes [-(1-1e-4),1-1e-4].
# phi_bound/theta_bound apply to physical coefficients only when the corresponding flag is FALSE.
# [-bound, bound] for AR and MA coefficients, respectively.
.estimate_qmle_no_x <- function(Y, p, q, r,
                          heteroskedastic = TRUE,
                          n_starts        = 20L,
                          maxit           = 3000L,
                          factr           = 1e4,
                          seed            = NULL,
                          verbose         = FALSE,
                          compute_se      = TRUE,
                          warn_dynamic    = TRUE,
                          phi_start       = NULL,
                          hessian_type    = c("obs", "info", "both"),
                          phi_bound       = 5,
                          theta_bound     = 5,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  if (!is.null(seed)) set.seed(seed)
  hessian_type <- match.arg(hessian_type)
  n_starts <- .validate_n_starts(n_starts)
  dims <- .validate_inputs(Y, p, q, r, heteroskedastic)
  N <- dims$N; big_T <- dims$big_T

  Y_bar <- colMeans(Y);  Y_c <- sweep(Y, 2L, Y_bar);  S_N <- crossprod(Y_c) / N

  d0        <- count_params(p, q, r, big_T, heteroskedastic)
  box       <- make_box_constraints(
    p, q, r, big_T, heteroskedastic,
    phi_bound = phi_bound, theta_bound = theta_bound
  ,
      stationary = stationary, invertible = invertible)
  initial_starts <- .make_initial_starts(
    Y, box, p, q, r, big_T, heteroskedastic, n_starts,
    phi_start = phi_start,
    ndeps = ndeps, start_factr = start_factr
  )
  start_types <- attr(initial_starts, "start_type")
  start_phis <- attr(initial_starts, "initial_phi")
  start_anchor <- initial_starts[[1L]]

  best_val <- Inf
  best_par <- NULL
  best_res <- list(convergence = 999L)
  all_runs <- vector("list", n_starts)
  max_start_attempts <- 25L
  stable_start <- .make_stable_start(
    box, p, q, r, big_T, heteroskedastic
  )

  for (start in seq_len(n_starts)) {
    start_success <- FALSE
    last_message <- NA_character_
    for (attempt in seq_len(max_start_attempts)) {
      replacement_error <- NULL
      psi_init <- if (attempt == 1L) {
        current_type <- start_types[start]
        current_phi <- start_phis[[start]]
        initial_starts[[start]]
      } else if (attempt == max_start_attempts) {
        current_type <- stable_start$type
        current_phi <- stable_start$phi
        stable_start$psi
      } else {
        replacement <- tryCatch(
          .draw_replacement_start(
            start, start_anchor, Y, box, p, q, r, big_T,
            heteroskedastic,
            ndeps = ndeps, start_factr = start_factr
          ),
          error = function(e) {
            replacement_error <<- conditionMessage(e)
            NULL
          }
        )
        if (is.null(replacement)) {
          current_type <- "replacement_generation_failed"
          current_phi <- rep(NA_real_, p)
          NULL
        } else {
          current_type <- replacement$type
          current_phi <- replacement$phi
          replacement$psi
        }
      }

      res <- if (is.null(psi_init)) {
        list(
          value = Inf, convergence = 999L, par = rep(NA_real_, d0),
          message = replacement_error
        )
      } else {
        .run_box_optim(
          psi_init,
          function(psi) .qll_neg(
            psi, p, q, r, big_T, S_N, heteroskedastic
          ,
      stationary = stationary, invertible = invertible),
          box$lower, box$upper, maxit, factr,
          ndeps = ndeps
        )
      }

      run_ok <- is.finite(res$value) && res$convergence == 0L &&
        res$value < 1e15
      if (!run_ok) {
        last_message <- if (!is.null(res$message) &&
                            length(res$message) == 1L &&
                            !is.na(res$message) && nzchar(res$message)) {
          as.character(res$message)
        } else {
          "nonfinite, penalized, or nonconverged optimization"
        }
        if (verbose) cat(sprintf(
          "  QMLE start %2d/%2d attempt %d failed.\n",
          start, n_starts, attempt
        ))
        next
      }
      start_success <- TRUE
      break
    }

    if (!start_success) {
      all_runs[[start]] <- list(
        obj = Inf,
        phi = rep(NA_real_, p),
        theta = rep(NA_real_, q),
        conv = 999L,
        attempts = max_start_attempts,
        start_type = current_type,
        initial_phi = current_phi,
        success = FALSE,
        message = last_message
      )
      warning(sprintf(
        paste0("QMLE start %d failed after %d attempts; ",
               "continuing with the remaining starts."),
        start, max_start_attempts
      ), call. = FALSE)
      next
    }

    fitted_dynamic <- .decode_dynamics(
      if (p > 0L) res$par[seq_len(p)] else numeric(0),
      if (q > 0L) res$par[p + seq_len(q)] else numeric(0),
      stationary, invertible
    )
    phi_r <- round(fitted_dynamic$phi, 4)
    theta_r <- round(fitted_dynamic$theta, 4)
    if (verbose)
      cat(sprintf("  QMLE start %2d/%2d  obj=%11.4f  conv=%d  phi=[%s]  theta=[%s]\n",
                  start, n_starts, -res$value, res$convergence,
                  paste(sprintf("%.4f", phi_r), collapse=", "),
                  paste(sprintf("%.4f", theta_r), collapse=", ")))
    all_runs[[start]] <- list(obj=res$value, phi=phi_r, theta=theta_r,
                              conv=res$convergence, attempts=attempt,
                              start_type=current_type,
                              initial_phi=current_phi,
                              success=TRUE,
                              message=res$message)

    if (res$value < best_val) {
      best_val <- res$value;  best_par <- res$par;  best_res <- res
    }
  }

  if (is.null(best_par)) {
    result <- .failed_optimizer_result(
      "QMLE", p, q, r, big_T, N, heteroskedastic, all_runs,
      obj_val = -Inf
    ,
      stationary = stationary, invertible = invertible,
      phi_bound = phi_bound, theta_bound = theta_bound)
    result$hessian_type <- hessian_type
    result$S_N <- S_N
    result$pd_tolerance <- 0
    result$Sigma_hat_construction <- "unavailable: all starts failed"
    result$parameterization <- list(
      optimization = "vech(L_Psi), vec(A2=F2 L_Psi), Schur-factor vech(L_eta)",
      reported = "vech(Psi); Pi_(p,q): initial variances, then lexicographic (s,t), s<t, covariances, then eta/innovation cross-covariances"
    )
    result <- .set_dynamic_domain(result, stationary, invertible, phi_bound, theta_bound)
    return(result)
  }
  psi_hat <- best_par
  vartheta <- unpack_params(psi_hat, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
  Sigma_hat <- tryCatch(
    .build_Sigma_internal(psi_hat, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible),
    error = function(e) NULL
  )
  if (is.null(Sigma_hat) || any(!is.finite(Sigma_hat)))
    stop("Internal error: the selected QMLE solution has a nonfinite covariance matrix.")
  pd_fitted <- .check_pd_components(vartheta, p, q, r, big_T,
                                    heteroskedastic, min_eig = 0,
                                    Sigma_override = Sigma_hat)
  if (!isTRUE(pd_fitted$all_valid)) {
    component_names <- c("innovation_variances", "Psi",
                         "Sigma_eta", "Sigma_u", "Sigma_xi", "Sigma")
    failed_components <- component_names[
      !vapply(pd_fitted[component_names], isTRUE, logical(1L))
    ]
    warning(paste0("The selected QMLE solution violates the covariance ",
                   "restrictions. Failed component(s): ",
                   paste(failed_components, collapse = ", "), "."))
  }

  vartheta$Psi <- .Psi_from_lower(vartheta$Psi_lower, r)
  vartheta_hat <- .to_paper_parameters(
    psi_hat, p, q, r, big_T, heteroskedastic
  ,
      stationary = stationary, invertible = invertible)
  names(psi_hat) <- .parameter_names(
    p, q, r, big_T, heteroskedastic, internal = TRUE
  ,
      stationary = stationary, invertible = invertible)
  names(vartheta_hat) <- .parameter_names(
    p, q, r, big_T, heteroskedastic, internal = FALSE
  ,
      stationary = stationary, invertible = invertible)
  result <- list(
    coef=vartheta, psi_hat=psi_hat, vartheta_hat=vartheta_hat,
    estimator="QMLE", obj_val=-best_val,
    converged=(best_res$convergence == 0L && isTRUE(pd_fitted$all_valid)),
    convergence_code=best_res$convergence,
    pd_components=pd_fitted, Sigma_hat=Sigma_hat,
    pd_tolerance=0,
    Sigma_hat_construction="internal optimization coordinates",
    all_runs=all_runs,
    successful_starts=sum(vapply(all_runs, function(x) isTRUE(x$success), logical(1L))),
    failed_starts=sum(vapply(all_runs, function(x) !isTRUE(x$success), logical(1L))),
    parameterization=list(
      optimization="vech(L_Psi), vec(A2=F2 L_Psi), Schur-factor vech(L_eta)",
      reported="vech(Psi); Pi_(p,q): initial variances, then lexicographic (s,t), s<t, covariances, then eta/innovation cross-covariances"
    ),
    hessian_type=hessian_type,
    bound_hits=.bound_hit_diagnostics(psi_hat, box),
    S_N=S_N, p=p, q=q, r=r, big_T=big_T, N=N
  )
  result <- .set_dynamic_domain(result, stationary, invertible, phi_bound, theta_bound)
  class(result) <- c("arma_panel_fit", "list")

  if (!compute_se) {
    result <- .set_empty_inference(result)
    if (warn_dynamic) .warn_dynamic_roots(vartheta, "QMLE")
    result <- .set_dynamic_domain(result, stationary, invertible, phi_bound, theta_bound)
    return(result)
  }

  result <- .stable_qmle_inference(result,
    list(N=N,big_T=big_T,Y_c=Y_c),p,q,r,heteroskedastic,
    hessian_type, with_x = FALSE)
  if (warn_dynamic) .warn_dynamic_roots(vartheta,"QMLE")
  result
}

# ============================================================================ #
#  SECTION 7.  GMM MOMENTS AND OBJECTIVE  [Paper Section 3.2]
# ============================================================================ #

.validate_W <- function(W, dot_T) {
  if (!is.matrix(W) || !is.numeric(W) ||
      any(dim(W) != c(dot_T, dot_T)))
    stop("W must be a numeric dot_T by dot_T matrix.")
  if (any(!is.finite(W))) stop("W must contain only finite values.")

  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(W)))
  if (max(abs(W - t(W))) > tolerance)
    stop("W must be symmetric.")
  W <- (W + t(W)) / 2
  if (is.null(tryCatch(chol(W), error = function(e) NULL)))
    stop("W must be positive definite.")
  W
}

# Compute the feasible individual moments and the paper's data-based estimator
# of Upsilon.  The columns of m_i follow the same order as vech().
.gmm_data_moments <- function(Y, compute_Upsilon = TRUE) {
  N <- nrow(Y)
  big_T <- ncol(Y)
  if (N < 2L)
    stop("GMM requires N >= 2 for the finite-sample correction.")
  c_N <- (N - 1) / N
  Y_bar <- colMeans(Y)
  Y_c <- sweep(Y, 2L, Y_bar)
  S_N <- crossprod(Y_c) / N
  m_N <- vech(S_N)

  Upsilon_hat <- NULL
  if (compute_Upsilon) {
    index <- which(
      lower.tri(matrix(FALSE, big_T, big_T), diag = TRUE),
      arr.ind = TRUE
    )
    m_i <- Y_c[, index[, "row"], drop = FALSE] *
      Y_c[, index[, "col"], drop = FALSE]
    centered_m_i <- sweep(m_i, 2L, m_N)
    Upsilon_hat <- crossprod(centered_m_i) / (N * c_N^2)
    Upsilon_hat <- (Upsilon_hat + t(Upsilon_hat)) / 2
  }

  list(Y_bar = Y_bar, Y_c = Y_c, S_N = S_N, m_N = m_N, c_N = c_N,
       Upsilon_hat = Upsilon_hat)
}

# m(psi) = vech(Sigma(psi)).  Sigma is built directly from the internal
# optimization coordinates, as in the QMLE objective above.
.sigma_vec <- function(psi, p, q, r, big_T, heteroskedastic,
                                  stationary = TRUE, invertible = TRUE) {
  dot_T <- as.integer(big_T * (big_T + 1L) / 2L)
  tryCatch({
    Sigma <- .build_Sigma_internal(
      psi, p, q, r, big_T, heteroskedastic
    ,
      stationary = stationary, invertible = invertible)
    if (is.null(Sigma)) return(rep(NA_real_, dot_T))
    value <- vech(Sigma)
    if (any(!is.finite(value))) rep(NA_real_, dot_T) else value
  }, error = function(e) rep(NA_real_, dot_T))
}

# Normalized GMM criterion Q_N(psi) = g_N(psi)' W g_N(psi).
.gmm_obj <- function(psi, m_N, c_N, W, p, q, r, big_T, heteroskedastic,
                                  stationary = TRUE, invertible = TRUE) {
  model_moments <- .sigma_vec(
    psi, p, q, r, big_T, heteroskedastic
  ,
      stationary = stationary, invertible = invertible)
  if (any(!is.finite(model_moments))) return(1e15)
  g_N <- m_N / c_N - model_moments
  objective <- as.numeric(crossprod(g_N, W %*% g_N))
  if (!is.finite(objective)) 1e15 else objective
}

# ============================================================================ #
#  SECTION 8.  GMM WITH AN ARBITRARY WEIGHTING MATRIX  [Paper eq. 25]
# ============================================================================ #

#' GMM estimator for an arbitrary positive definite weighting matrix
#'
#' @param Y               N by T data matrix.
#' @param p,q,r           Model orders.
#' @param W               Positive definite weighting matrix. Default: identity.
#' @param heteroskedastic Logical.
#' @param n_starts        Number of optimization starts.
#' @param phi_start       Optional AR vector for the first conditional start.
#' @param phi_bound      Physical AR coefficient limit when stationary=FALSE (default 5).
#' @param theta_bound    Physical MA coefficient limit when invertible=FALSE (default 5).
#' @param compute_se      If FALSE, skip inference calculations.
#' @param efficient_weight Internal indicator used by estimate_gmm_efficient().
#' @return Named list containing the estimates and GMM inference.
.estimate_gmm_no_x <- function(Y, p, q, r,
                         W                = NULL,
                         heteroskedastic  = TRUE,
                         n_starts         = 20L,
                         maxit            = 3000L,
                         factr            = 1e4,
                         seed             = NULL,
                         verbose          = FALSE,
                         compute_se       = TRUE,
                         warn_dynamic     = TRUE,
                         efficient_weight = FALSE,
                         .warn_T          = TRUE,
                         .gmm_data        = NULL,
                         phi_start       = NULL,
                         phi_bound       = 5,
                         theta_bound     = 5,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  if (!is.null(seed)) set.seed(seed)
  n_starts <- .validate_n_starts(n_starts)
  dims <- .validate_inputs(
    Y, p, q, r, heteroskedastic, warn_T = .warn_T
  )
  N <- dims$N
  big_T <- dims$big_T
  dot_T <- dims$dot_T

  gmm_data <- if (is.null(.gmm_data))
    .gmm_data_moments(Y, compute_Upsilon = compute_se) else .gmm_data
  S_N <- gmm_data$S_N
  m_N <- gmm_data$m_N
  c_N <- gmm_data$c_N
  Upsilon_hat <- gmm_data$Upsilon_hat

  if (is.null(W)) W <- diag(dot_T)
  W <- .validate_W(W, dot_T)

  d0 <- count_params(p, q, r, big_T, heteroskedastic)
  box <- make_box_constraints(
    p, q, r, big_T, heteroskedastic,
    phi_bound = phi_bound, theta_bound = theta_bound
  ,
      stationary = stationary, invertible = invertible)
  initial_starts <- .make_initial_starts(
    Y, box, p, q, r, big_T, heteroskedastic, n_starts,
    phi_start = phi_start,
    ndeps = ndeps, start_factr = start_factr
  )
  start_types <- attr(initial_starts, "start_type")
  start_phis <- attr(initial_starts, "initial_phi")
  start_anchor <- initial_starts[[1L]]

  best_val <- Inf
  best_par <- NULL
  best_res <- list(convergence = 999L)
  all_runs <- vector("list", n_starts)
  max_start_attempts <- 25L
  stable_start <- .make_stable_start(
    box, p, q, r, big_T, heteroskedastic
  )

  for (start in seq_len(n_starts)) {
    start_success <- FALSE
    last_message <- NA_character_
    for (attempt in seq_len(max_start_attempts)) {
      replacement_error <- NULL
      psi_init <- if (attempt == 1L) {
        current_type <- start_types[start]
        current_phi <- start_phis[[start]]
        initial_starts[[start]]
      } else if (attempt == max_start_attempts) {
        current_type <- stable_start$type
        current_phi <- stable_start$phi
        stable_start$psi
      } else {
        replacement <- tryCatch(
          .draw_replacement_start(
            start, start_anchor, Y, box, p, q, r, big_T,
            heteroskedastic,
            ndeps = ndeps, start_factr = start_factr
          ),
          error = function(e) {
            replacement_error <<- conditionMessage(e)
            NULL
          }
        )
        if (is.null(replacement)) {
          current_type <- "replacement_generation_failed"
          current_phi <- rep(NA_real_, p)
          NULL
        } else {
          current_type <- replacement$type
          current_phi <- replacement$phi
          replacement$psi
        }
      }

      res <- if (is.null(psi_init)) {
        list(
          value = Inf, convergence = 999L, par = rep(NA_real_, d0),
          message = replacement_error
        )
      } else {
        .run_box_optim(
          psi_init,
          function(psi) .gmm_obj(
            psi, m_N, c_N, W, p, q, r, big_T, heteroskedastic
          ,
      stationary = stationary, invertible = invertible),
          box$lower, box$upper, maxit, factr,
          ndeps = ndeps
        )
      }

      run_ok <- is.finite(res$value) && res$convergence == 0L &&
        res$value < 1e15
      if (!run_ok) {
        last_message <- if (!is.null(res$message) &&
                            length(res$message) == 1L &&
                            !is.na(res$message) && nzchar(res$message)) {
          as.character(res$message)
        } else {
          "nonfinite, penalized, or nonconverged optimization"
        }
        if (verbose) cat(sprintf(
          "  GMM start %2d/%2d attempt %d failed.\n",
          start, n_starts, attempt
        ))
        next
      }
      start_success <- TRUE
      break
    }

    if (!start_success) {
      all_runs[[start]] <- list(
        obj = Inf,
        phi = rep(NA_real_, p),
        theta = rep(NA_real_, q),
        conv = 999L,
        attempts = max_start_attempts,
        start_type = current_type,
        initial_phi = current_phi,
        success = FALSE,
        message = last_message
      )
      warning(sprintf(
        paste0("GMM start %d failed after %d attempts; ",
               "continuing with the remaining starts."),
        start, max_start_attempts
      ), call. = FALSE)
      next
    }

    fitted_dynamic <- .decode_dynamics(
      if (p > 0L) res$par[seq_len(p)] else numeric(0),
      if (q > 0L) res$par[p + seq_len(q)] else numeric(0),
      stationary, invertible
    )
    phi_r <- round(fitted_dynamic$phi, 4)
    theta_r <- round(fitted_dynamic$theta, 4)
    if (verbose) cat(sprintf(
      paste0("  GMM  start %2d/%2d  obj=%11.4f  conv=%d  ",
             "phi=[%s]  theta=[%s]\n"),
      start, n_starts, res$value, res$convergence,
      paste(sprintf("%.4f", phi_r), collapse = ", "),
      paste(sprintf("%.4f", theta_r), collapse = ", ")
    ))
    all_runs[[start]] <- list(
      obj = res$value, phi = phi_r, theta = theta_r,
      conv = res$convergence, attempts = attempt,
      start_type = current_type, initial_phi = current_phi,
      success = TRUE, message = res$message
    )

    if (res$value < best_val) {
      best_val <- res$value
      best_par <- res$par
      best_res <- res
    }
  }

  if (is.null(best_par)) {
    estimator_name <- if (isTRUE(efficient_weight)) "Efficient GMM" else "GMM"
    result <- .failed_optimizer_result(
      estimator_name, p, q, r, big_T, N, heteroskedastic, all_runs,
      obj_val = Inf
    ,
      stationary = stationary, invertible = invertible,
      phi_bound = phi_bound, theta_bound = theta_bound)
    result$parameterization <- list(
      optimization = paste0(
        "vech(L_Psi), vec(A2=F2 L_Psi), ",
        "Schur-factor vech(L_eta)"
      ),
      reported = "vech(Psi); Pi_(p,q): initial variances, then lexicographic (s,t), s<t, covariances, then eta/innovation cross-covariances"
    )
    result$S_N <- S_N
    result$W <- W
    result$Upsilon_hat <- if (is.null(Upsilon_hat))
      matrix(NA_real_, dot_T, dot_T) else Upsilon_hat
    result$G <- matrix(
      NA_real_, dot_T, d0,
      dimnames = list(NULL, names(result$vartheta_hat))
    )
    result$rank_G <- NA_integer_
    result$rcond_GWG <- NA_real_
    result$J_stat <- NA_real_
    result$df_overid <- NA_integer_
    result$J_pval <- NA_real_
    result <- .set_dynamic_domain(result, stationary, invertible, phi_bound, theta_bound)
    return(result)
  }

  psi_hat <- best_par
  vartheta <- unpack_params(
    psi_hat, p, q, r, big_T, heteroskedastic
  ,
      stationary = stationary, invertible = invertible)
  pd_fitted <- .check_pd_components(
    vartheta, p, q, r, big_T, heteroskedastic,
    min_eig = 0
  )
  if (!isTRUE(pd_fitted$all_valid)) {
    component_names <- c("innovation_variances", "Psi",
                         "Sigma_eta", "Sigma_u", "Sigma_xi", "Sigma")
    failed_components <- component_names[
      !vapply(pd_fitted[component_names], isTRUE, logical(1L))
    ]
    warning(paste0("The selected GMM solution violates the covariance ",
                   "restrictions. Failed component(s): ",
                   paste(failed_components, collapse = ", "), "."))
  }

  vartheta$Psi <- .Psi_from_lower(vartheta$Psi_lower, r)
  vartheta_hat <- .to_paper_parameters(
    psi_hat, p, q, r, big_T, heteroskedastic
  ,
      stationary = stationary, invertible = invertible)
  names(psi_hat) <- .parameter_names(
    p, q, r, big_T, heteroskedastic, internal = TRUE
  ,
      stationary = stationary, invertible = invertible)
  names(vartheta_hat) <- .parameter_names(
    p, q, r, big_T, heteroskedastic, internal = FALSE
  ,
      stationary = stationary, invertible = invertible)
  Sigma_hat <- .build_Sigma_internal(
    psi_hat, p, q, r, big_T, heteroskedastic,
    stationary=stationary,invertible=invertible
  )

  result <- list(
    coef = vartheta,
    psi_hat = psi_hat,
    vartheta_hat = vartheta_hat,
    estimator = if (isTRUE(efficient_weight)) "Efficient GMM" else "GMM",
    obj_val = best_val,
    converged = best_res$convergence == 0L &&
      isTRUE(pd_fitted$all_valid),
    convergence_code = best_res$convergence,
    pd_components = pd_fitted,
    Sigma_hat = Sigma_hat,
    all_runs = all_runs,
    successful_starts = sum(vapply(
      all_runs, function(x) isTRUE(x$success), logical(1L)
    )),
    failed_starts = sum(vapply(
      all_runs, function(x) !isTRUE(x$success), logical(1L)
    )),
    parameterization = list(
      optimization = paste0(
        "vech(L_Psi), vec(A2=F2 L_Psi), ",
        "Schur-factor vech(L_eta)"
      ),
      reported = "vech(Psi); Pi_(p,q): initial variances, then lexicographic (s,t), s<t, covariances, then eta/innovation cross-covariances"
    ),
    bound_hits = .bound_hit_diagnostics(psi_hat, box),
    S_N = S_N,
    W = W,
    Upsilon_hat = if (is.null(Upsilon_hat))
      matrix(NA_real_, dot_T, dot_T) else Upsilon_hat,
    hansen_test = NULL,
    p = p, q = q, r = r, big_T = big_T, N = N
  )
  result <- .set_dynamic_domain(result, stationary, invertible, phi_bound, theta_bound)
  class(result) <- c("arma_panel_fit", "list")

  if (!compute_se) {
    result <- .set_empty_inference(result)
    result$G <- matrix(
      NA_real_, dot_T, d0,
      dimnames = list(NULL, names(vartheta_hat))
    )
    result$rank_G <- NA_integer_
    result$rcond_GWG <- NA_real_
    result$J_stat <- NA_real_
    result$df_overid <- NA_integer_
    result$J_pval <- NA_real_
    result$rank_weighted_G <- NA_integer_
    result$weighted_G_singular_values <- rep(NA_real_, d0)
    result$inference_diagnostics <- list(first_failure = "not_computed")
    if (warn_dynamic) .warn_dynamic_roots(vartheta, result$estimator)
    result <- .set_dynamic_domain(result, stationary, invertible, phi_bound, theta_bound)
    return(result)
  }

  result <- .stable_gmm_inference(result,
    list(N=N,big_T=big_T,Y=Y,Y_c=sweep(Y,2,colMeans(Y)),S_N=S_N),
    p, q, r, heteroskedastic, with_x = FALSE)
  if (warn_dynamic) .warn_dynamic_roots(vartheta,result$estimator)
  result
}

# ============================================================================ #
#  SECTION 9.  EFFICIENT GMM WITH DATA-BASED UPSILON  [Paper eq. 27]
# ============================================================================ #

#' Efficient GMM without a preliminary GMM estimation
#'
#' Upsilon_hat is computed once from the centered individual second moments,
#' W_opt = Upsilon_hat^{-1} is formed once, and estimate_gmm() is run once.
#' @param phi_start Optional AR vector for the first conditional start.
#' @param phi_bound Physical AR coefficient limit when stationary=FALSE (default 5).
#' @param theta_bound Physical MA coefficient limit when invertible=FALSE (default 5).
.estimate_gmm_efficient_no_x <- function(Y, p, q, r,
                                   heteroskedastic = TRUE,
                                   n_starts        = 20L,
                                   maxit           = 3000L,
                                   factr           = 1e4,
                                   seed            = NULL,
                                   verbose         = FALSE,
                                   compute_se      = TRUE,
                                   warn_dynamic    = TRUE,
                                   phi_start       = NULL,
                                   phi_bound       = 5,
                                   theta_bound     = 5,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  phi_bound <- .validate_coefficient_bound(phi_bound, "phi_bound")
  theta_bound <- .validate_coefficient_bound(theta_bound, "theta_bound")
  n_starts <- .validate_n_starts(n_starts)
  dims <- .validate_inputs(Y, p, q, r, heteroskedastic)
  N <- dims$N
  big_T <- dims$big_T
  dot_T <- dims$dot_T

  minimum_N <- if (big_T == 1L) 3L else dot_T + 1L
  if (N < minimum_N) {
    stop(sprintf(
      paste0("Efficient GMM requires N >= %d when T=%d (dot_T=%d); ",
             "received N=%d. Otherwise Upsilon_hat cannot be positive definite."),
      minimum_N, big_T, dot_T, N
    ))
  }

  gmm_data <- .gmm_data_moments(Y, compute_Upsilon = TRUE)
  Upsilon_hat <- gmm_data$Upsilon_hat
  R_Upsilon <- tryCatch(
    chol(Upsilon_hat),
    error = function(e) NULL
  )
  if (is.null(R_Upsilon)) {
    rank_Upsilon <- tryCatch(qr(Upsilon_hat)$rank,
                             error = function(e) NA_integer_)
    stop(sprintf(
      paste0("Efficient GMM requires the data-based Upsilon_hat to be ",
             "positive definite (estimated rank %s of %d)."),
      ifelse(is.na(rank_Upsilon), "NA", rank_Upsilon), dot_T
    ))
  }
  W_opt <- chol2inv(R_Upsilon)

  if (verbose)
    cat("=== Efficient GMM with data-based Upsilon_hat ===\n")
  result <- estimate_gmm(
    Y = Y, p = p, q = q, r = r,
    W = W_opt,
    heteroskedastic = heteroskedastic,
    n_starts = n_starts,
    phi_start = phi_start,
    maxit = maxit,
    factr = factr,
    seed = seed,
    verbose = verbose,
    compute_se = compute_se,
    phi_bound = phi_bound,
    theta_bound = theta_bound,
    warn_dynamic = FALSE,
    efficient_weight = TRUE,
    .warn_T = FALSE,
    .gmm_data = gmm_data
  ,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      )
  result$W_opt <- W_opt
  result$Upsilon_hat <- Upsilon_hat
  if (warn_dynamic)
    .warn_dynamic_roots(result$coef, "efficient GMM")
  result
}

# ============================================================================ #
#  SECTION 10.  CONCISE RESULTS SUMMARY
# ============================================================================ #

.summary_coefficient_table <- function(object) {
  estimate <- if (!is.null(object$parameter_hat)) object$parameter_hat else object$vartheta_hat
  parameter_names <- names(estimate)
  if (is.null(parameter_names))
    parameter_names <- paste0("parameter_", seq_along(estimate))

  standard_error <- object$se
  if (is.null(standard_error) || length(standard_error) != length(estimate) ||
      all(!is.finite(standard_error))) {
    table <- cbind(Estimate = estimate)
  } else {
    z_value <- ifelse(
      is.finite(standard_error) & standard_error > 0,
      estimate / standard_error,
      NA_real_
    )
    p_value <- ifelse(
      is.finite(z_value),
      2 * pnorm(abs(z_value), lower.tail = FALSE),
      NA_real_
    )
    table <- cbind(
      Estimate = estimate,
      `Std. Error` = standard_error,
      `z value` = z_value,
      `Pr(>|z|)` = p_value
    )
  }
  rownames(table) <- parameter_names
  table
}

summary.arma_panel_fit <- function(object, ...) {
  all_coefficients <- .summary_coefficient_table(object)
  main_coefficients <- grepl("^(phi|theta|beta)_[0-9]+$",
                             rownames(all_coefficients))
  result <- list(
    estimator = object$estimator,
    p = object$p,
    q = object$q,
    r = object$r,
    N = object$N,
    big_T = object$big_T,
    converged = isTRUE(object$converged),
    convergence_code = object$convergence_code,
    feasible = isTRUE(object$pd_components$all_valid),
    objective = object$obj_val,
    coefficients = all_coefficients[main_coefficients, , drop = FALSE],
    all_coefficients = all_coefficients,
    hansen_test = object$hansen_test
  )
  class(result) <- c("summary.arma_panel_fit", "list")
  result
}

print.summary.arma_panel_fit <- function(x,
                                         digits = max(3L, getOption("digits") - 3L),
                                         ...) {
  cat(sprintf("%s results\n", x$estimator))
  cat(sprintf("Model: (p, q, r) = (%d, %d, %d)\n", x$p, x$q, x$r))
  cat(sprintf("Observations: N = %d, T = %d\n", x$N, x$big_T))
  cat(sprintf(
    "Converged: %s (code %d); covariance feasible: %s\n",
    if (x$converged) "yes" else "no",
    x$convergence_code,
    if (x$feasible) "yes" else "no"
  ))
  objective_label <- if (identical(x$estimator, "QMLE"))
    "Concentrated QMLE criterion" else "Normalized GMM criterion"
  cat(sprintf("%s: %.*g\n\n", objective_label, digits, x$objective))

  cat("Coefficients:\n")
  if (nrow(x$coefficients) == 0L) {
    cat("  No phi, theta, or beta coefficients in this specification.\n")
  } else if (ncol(x$coefficients) == 4L) {
    printCoefmat(
      x$coefficients,
      digits = digits,
      signif.stars = getOption("show.signif.stars"),
      na.print = "NA",
      P.values = TRUE,
      has.Pvalue = TRUE
    )
  } else {
    print(x$coefficients, digits = digits)
  }

  if (grepl("GMM", x$estimator, fixed = TRUE)) {
    cat("\nHansen's J test of overidentifying restrictions:\n")
    if (identical(x$estimator, "Efficient GMM")) {
      if (!is.null(x$hansen_test)) {
        cat(sprintf(
          "  J = %.*g, df = %d, p-value = %s\n",
          digits,
          unname(x$hansen_test$statistic),
          unname(x$hansen_test$parameter),
          format.pval(x$hansen_test$p.value, digits = digits)
        ))
      } else {
        cat("  Not available: the model is not overidentified or inference failed.\n")
      }
    } else {
      cat("  Not reported: Hansen's reference test requires efficient weighting.\n")
    }
  }

  invisible(x)
}

print.arma_panel_fit <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}

# ============================================================================ #
#  SECTION 11. QMLE WITH STRICTLY EXOGENOUS REGRESSORS [Paper Section 5.1]
# ============================================================================ #
# Supply X as an N x T matrix for one regressor or an N x T x k array.
# W columns are time-major: (x_i1', ..., x_iT').  Initial c_s coefficients
# use this COMPLETE history, and beta enters only the equations dated t>p.
# Xi ordering: beta, c_1,...,c_p (each length kT), vec(A_mundlak).
# vartheta retains the original paper and internal covariance coordinates.
# Only delta is profiled; vartheta and Xi are optimized simultaneously.
#
# Additional output:
#   Xi_hat: beta, c_s, vec(A_mundlak); delta_hat: recovered T-vector.
#   parameter_hat: reported (vartheta, Xi, delta); optimization_hat: (psi, Xi).
#   coef$c is p x kT; coef$A_mundlak is r x kT (distinct from internal A2).
#   fitted, residuals: N x T; S_N: residual cross-product divided by N.
#   se, ci, V_hat refer to parameter_hat; V_hat is the sqrt(N) asymptotic
#   covariance, as in the original code. se_vartheta/se_Xi/se_delta are slices.
# Rank-deficient unrestricted c/A histories are rejected rather than silently
# dropping columns or imposing restrictions absent from the paper.

.qmle_x_data <- function(Y, X) {
  if (!is.matrix(Y) || !is.numeric(Y) || any(!is.finite(Y)) ||
      nrow(Y) < 1L || ncol(Y) < 1L)
    stop("Y must be a nonempty finite numeric N x T matrix.")
  N <- nrow(Y)
  big_T <- ncol(Y)
  dx <- dim(X)
  if (!is.numeric(X) || is.null(dx) || !(length(dx) %in% c(2L, 3L)))
    stop("X must be a numeric N x T matrix or N x T x k array.")
  if (dx[1L] != N || dx[2L] != big_T)
    stop("The individual and date dimensions of X must match Y (N x T).")
  k <- if (length(dx) == 2L) 1L else dx[3L]
  if (k < 1L || any(!is.finite(X)))
    stop("X must contain at least one regressor and no missing or nonfinite values.")
  W <- matrix(0, N, big_T * k)
  for (tt in seq_len(big_T)) {
    cols <- (tt - 1L) * k + seq_len(k)
    W[, cols] <- if (length(dx) == 2L) X[, tt] else
      matrix(X[, tt, , drop = FALSE], N, k)
  }
  Y_bar <- colMeans(Y)
  W_bar <- colMeans(W)
  Y_c <- sweep(Y, 2L, Y_bar, "-")
  W_c <- sweep(W, 2L, W_bar, "-")
  # Cross-sectional QR compression makes conditional-GLS starts independent
  # of N in memory cost.  It leaves the least-squares coefficient solution
  # unchanged, including when W_c has deficient rank.
  qr_W <- qr(W_c, LAPACK = TRUE)
  R_W <- qr.R(qr_W)[, order(qr_W$pivot), drop = FALSE]
  Y_Q <- qr.qty(qr_W, Y_c)[seq_len(nrow(R_W)), , drop = FALSE]
  list(Y = Y, W = W, N = N, big_T = big_T, k = k,
       Y_bar = Y_bar, W_bar = W_bar, Y_c = Y_c, W_c = W_c,
       R_W = R_W, Y_Q = Y_Q)
}

.qmle_x_count_Xi <- function(p, r, k, big_T) {
  as.integer(k + (p + r) * k * big_T)
}

.qmle_x_Xi_names <- function(p, r, k, big_T) {
  L <- k * big_T
  c(paste0("beta_", seq_len(k)),
    if (p > 0L) unlist(lapply(seq_len(p), function(s)
      paste0("c_", s, "_w", seq_len(L))), use.names = FALSE) else character(0),
    if (r > 0L) unlist(lapply(seq_len(L), function(j)
      paste0("A_mundlak_", seq_len(r), "_w", j)), use.names = FALSE) else character(0))
}

.qmle_x_unpack_Xi <- function(Xi_vec, p, r, k, big_T) {
  if (length(Xi_vec) != .qmle_x_count_Xi(p, r, k, big_T))
    stop("Incorrect Xi-vector length.")
  L <- k * big_T
  pos <- k
  c_init <- if (p > 0L) matrix(Xi_vec[pos + seq_len(p * L)],
                              p, L, byrow = TRUE) else matrix(0, 0L, L)
  pos <- pos + p * L
  A_mundlak <- if (r > 0L) matrix(Xi_vec[pos + seq_len(r * L)], r, L) else
    matrix(0, 0L, L)
  list(beta = unname(Xi_vec[seq_len(k)]), c = c_init, A_mundlak = A_mundlak)
}

.qmle_x_C <- function(vartheta, Xi, p, r, k, big_T) {
  C <- matrix(0, big_T, k * big_T)
  if (p > 0L) C[seq_len(p), ] <- Xi$c
  for (tt in seq.int(p + 1L, big_T))
    C[tt, (tt - 1L) * k + seq_len(k)] <- Xi$beta
  if (r > 0L) {
    F <- rbind(diag(r), matrix(vartheta$F2, big_T - r, r))
    C <- C + F %*% Xi$A_mundlak
  }
  C
}

.qmle_x_neg <- function(optimization_vec, data, p, q, r, heteroskedastic,
                                  stationary = TRUE, invertible = TRUE) {
  d_v <- count_params(p, q, r, data$big_T, heteroskedastic)
  n_Xi <- .qmle_x_count_Xi(p, r, data$k, data$big_T)
  if (length(optimization_vec) != d_v + n_Xi ||
      any(!is.finite(optimization_vec))) return(1e15)
  tryCatch({
    psi <- optimization_vec[seq_len(d_v)]
    vartheta <- unpack_params(psi, p, q, r, data$big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
    Xi <- .qmle_x_unpack_Xi(optimization_vec[d_v + seq_len(n_Xi)],
                            p, r, data$k, data$big_T)
    C <- .qmle_x_C(vartheta, Xi, p, r, data$k, data$big_T)
    Gamma <- build_Gamma_p(vartheta$phi, data$big_T)
    Z <- data$Y_c - data$W_c %*% t(Gamma %*% C)
    .qll_neg(psi, p, q, r, data$big_T, crossprod(Z) / data$N,
             heteroskedastic,
      stationary = stationary, invertible = invertible)
  }, error = function(e) 1e15)
}

# Conditional GLS is used for initialization and rank diagnostics ONLY.
# All Xi coordinates remain in the simultaneous numerical optimization.
.qmle_x_gls_start <- function(psi, data, p, q, r, heteroskedastic,
                                  stationary = TRUE, invertible = TRUE) {
  vartheta <- unpack_params(psi, p, q, r, data$big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
  Sigma <- .build_Sigma_internal(psi, p, q, r, data$big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
  R <- chol(Sigma)
  Gamma <- build_Gamma_p(vartheta$phi, data$big_T)
  n_Xi <- .qmle_x_count_Xi(p, r, data$k, data$big_T)
  design <- matrix(0, nrow(data$R_W) * data$big_T, n_Xi)
  for (a in seq_len(n_Xi)) {
    unit <- numeric(n_Xi)
    unit[a] <- 1
    C_a <- .qmle_x_C(vartheta,
                      .qmle_x_unpack_Xi(unit, p, r, data$k, data$big_T),
                      p, r, data$k, data$big_T)
    mu_a <- data$R_W %*% t(Gamma %*% C_a)
    design[, a] <- as.vector(forwardsolve(t(R), t(mu_a)))
  }
  scales <- sqrt(colSums(design^2))
  if (any(!is.finite(scales)) || any(scales <= 0))
    stop("The centered mean design has an unidentified Xi coefficient.")
  normalized <- sweep(design, 2L, scales, "/")
  qr_design <- qr(normalized, tol = 1e-10)
  if (qr_design$rank < n_Xi)
    stop("The centered mean design is not full column rank for Xi.")
  response <- as.vector(forwardsolve(t(R), t(data$Y_Q)))
  Xi_start <- as.numeric(qr.coef(qr_design, response)) / scales
  if (any(!is.finite(Xi_start))) stop("Unable to construct a finite GLS Xi start.")
  list(Xi = Xi_start, rank = qr_design$rank,
       rcond = tryCatch(rcond(crossprod(normalized)), error = function(e) NA_real_))
}

.qmle_x_validate_design <- function(data, p, r) {
  if (p + r > 0L) {
    L <- data$k * data$big_T
    scales <- sqrt(colSums(data$W_c^2))
    rank_W <- if (any(scales <= 0)) 0L else
      qr(sweep(data$W_c, 2L, scales, "/"), tol = 1e-10)$rank
    if (rank_W < L)
      stop(paste0("The unrestricted c_s and A specification requires the centered ",
                  "complete regressor history W to have full column rank. ",
                  "Time-invariant, duplicated, or common regressors require ",
                  "explicit restrictions or a reduced parameterization."))
    if (data$big_T - p <= r)
      stop("Identifying beta separately from unrestricted A requires T-p > r.")
  } else {
    design <- do.call(rbind, lapply(seq_len(data$big_T), function(tt)
      data$W_c[, (tt - 1L) * data$k + seq_len(data$k), drop = FALSE]))
    scales <- sqrt(colSums(design^2))
    if (any(scales <= 0) ||
        qr(sweep(design, 2L, scales, "/"), tol = 1e-10)$rank < data$k)
      stop("The centered contemporaneous regressor design is not full column rank for beta.")
  }
  invisible(TRUE)
}

.estimate_qmle_with_x <- function(Y, p, q, r, heteroskedastic = TRUE,
                                  n_starts = 20L, maxit = 3000L, factr = 1e4,
                                  seed = NULL, verbose = FALSE, compute_se = TRUE,
                                  warn_dynamic = TRUE, phi_start = NULL,
                                  hessian_type = c("obs", "info", "both"),
                                  phi_bound = 5, theta_bound = 5, X,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  if (!is.null(seed)) set.seed(seed)
  hessian_type <- match.arg(hessian_type)
  n_starts <- .validate_n_starts(n_starts)
  dims <- .validate_inputs(Y, p, q, r, heteroskedastic)
  data <- .qmle_x_data(Y, X)
  .qmle_x_validate_design(data, p, r)
  N <- dims$N
  big_T <- dims$big_T
  k <- data$k
  d_v <- count_params(p, q, r, big_T, heteroskedastic)
  n_Xi <- .qmle_x_count_Xi(p, r, k, big_T)
  d_opt <- d_v + n_Xi
  box <- make_box_constraints(p, q, r, big_T, heteroskedastic,
                              phi_bound = phi_bound, theta_bound = theta_bound,
      stationary = stationary, invertible = invertible)
  joint_box <- box
  joint_box$lower <- c(box$lower, rep(-Inf, n_Xi))
  joint_box$upper <- c(box$upper, rep(Inf, n_Xi))

  # Regressor-adjusted covariance starts. Residualization is only a starting
  # heuristic; the criterion always uses the original Y and W.
  qr_W <- qr(data$W_c, tol = 1e-10)
  Y_start <- qr.resid(qr_W, data$Y_c)
  if (any(!is.finite(Y_start)) || sum(Y_start^2) <= .Machine$double.eps * sum(data$Y_c^2))
    Y_start <- data$Y_c
  starts <- .make_initial_starts(Y_start, box, p, q, r, big_T,
                                 heteroskedastic, n_starts, phi_start = phi_start,
                                   ndeps = ndeps, start_factr = start_factr
                                 )
  start_types <- attr(starts, "start_type")
  start_phis <- attr(starts, "initial_phi")
  anchor <- starts[[1L]]
  stable <- .make_stable_start(box, p, q, r, big_T, heteroskedastic)
  best_val <- Inf
  best_par <- NULL
  best_res <- list(convergence = 999L)
  all_runs <- vector("list", n_starts)
  max_start_attempts <- 25L

  for (start in seq_len(n_starts)) {
    success <- FALSE
    last_message <- NA_character_
    for (attempt in seq_len(max_start_attempts)) {
      replacement_error <- NULL
      if (attempt == 1L) {
        psi_init <- starts[[start]]
        current_type <- start_types[start]
        current_phi <- start_phis[[start]]
      } else if (attempt == max_start_attempts) {
        psi_init <- stable$psi
        current_type <- stable$type
        current_phi <- stable$phi
      } else {
        replacement <- tryCatch(.draw_replacement_start(
          start, anchor, Y_start, box, p, q, r, big_T, heteroskedastic,
          ndeps = ndeps, start_factr = start_factr
        ), error = function(e) {
          replacement_error <<- conditionMessage(e)
          NULL
        })
        psi_init <- if (is.null(replacement)) NULL else replacement$psi
        current_type <- if (is.null(replacement)) "replacement_generation_failed" else replacement$type
        current_phi <- if (is.null(replacement)) rep(NA_real_, p) else replacement$phi
      }
      Xi_init <- if (is.null(psi_init)) NULL else tryCatch(
        .qmle_x_gls_start(psi_init, data, p, q, r, heteroskedastic,
      stationary = stationary, invertible = invertible)$Xi,
        error = function(e) numeric(n_Xi)
      )
      res <- if (is.null(psi_init)) {
        list(value = Inf, convergence = 999L, par = rep(NA_real_, d_opt),
             message = replacement_error)
      } else .run_box_optim(
        c(psi_init, Xi_init),
        function(v) .qmle_x_neg(v, data, p, q, r, heteroskedastic,
      stationary = stationary, invertible = invertible),
        joint_box$lower, joint_box$upper, maxit, factr,
        ndeps = ndeps
      )
      success <- is.finite(res$value) && res$convergence == 0L && res$value < 1e15
      if (!success) {
        last_message <- if (!is.null(res$message) && length(res$message) == 1L &&
                            !is.na(res$message) && nzchar(res$message)) res$message else
          "nonfinite, penalized, or nonconverged optimization"
        if (verbose) cat(sprintf("  QMLE with X start %2d/%2d attempt %d failed.\n",
                                  start, n_starts, attempt))
        next
      }
      break
    }
    fitted_dynamic <- if (success) .decode_dynamics(
      if (p > 0L) res$par[seq_len(p)] else numeric(0),
      if (q > 0L) res$par[p + seq_len(q)] else numeric(0),
      stationary, invertible
    ) else list(phi = rep(NA_real_, p), theta = rep(NA_real_, q))
    phi_r <- round(fitted_dynamic$phi, 4)
    theta_r <- round(fitted_dynamic$theta, 4)
    all_runs[[start]] <- list(obj = if (success) res$value else Inf,
                             phi = phi_r, theta = theta_r,
                             conv = if (success) res$convergence else 999L,
                             attempts = attempt, start_type = current_type,
                             initial_phi = current_phi, success = success,
                             message = if (success) res$message else last_message)
    if (!success) {
      warning(sprintf("QMLE with X start %d failed after %d attempts; continuing with remaining starts.",
                       start, max_start_attempts), call. = FALSE)
      next
    }
    if (verbose) cat(sprintf("  QMLE with X start %2d/%2d obj=%11.4f conv=%d\n",
                              start, n_starts, -res$value, res$convergence))
    if (res$value < best_val) {
      best_val <- res$value
      best_par <- res$par
      best_res <- res
    }
  }

  Xi_names <- .qmle_x_Xi_names(p, r, k, big_T)
  delta_names <- paste0("delta_", seq_len(big_T))
  if (is.null(best_par)) {
    result <- .failed_optimizer_result("QMLE", p, q, r, big_T, N,
                                        heteroskedastic, all_runs, -Inf,
      stationary = stationary, invertible = invertible,
      phi_bound = phi_bound, theta_bound = theta_bound)
    result$Xi_hat <- setNames(rep(NA_real_, n_Xi), Xi_names)
    result$delta_hat <- setNames(rep(NA_real_, big_T), delta_names)
    result$parameter_hat <- c(result$vartheta_hat, result$Xi_hat, result$delta_hat)
    result$optimization_hat <- c(result$psi_hat, result$Xi_hat)
    result$coef$beta <- rep(NA_real_, k)
    result$coef$c <- matrix(NA_real_, p, k * big_T)
    result$coef$A_mundlak <- matrix(NA_real_, r, k * big_T)
    result$coef$delta <- result$delta_hat
    result$fitted <- result$residuals <- matrix(NA_real_, N, big_T)
    result$S_N <- matrix(NA_real_, big_T, big_T)
    result$C_hat <- result$K_hat <- matrix(NA_real_, big_T, k * big_T)
    result$k <- k
    result$has_regressors <- TRUE
    result$hessian_type <- hessian_type
    result$inference_diagnostics <- list(first_failure = "optimization_failed")
    return(.qmle_x_empty_inference(result))
  }
  psi_hat <- best_par[seq_len(d_v)]
  vartheta <- unpack_params(psi_hat, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
  vartheta$Psi <- .Psi_from_lower(vartheta$Psi_lower, r)
  Xi_hat <- setNames(best_par[d_v + seq_len(n_Xi)], Xi_names)
  Xi <- .qmle_x_unpack_Xi(Xi_hat, p, r, k, big_T)
  C_hat <- .qmle_x_C(vartheta, Xi, p, r, k, big_T)
  Gamma_hat <- build_Gamma_p(vartheta$phi, big_T)
  K_hat <- Gamma_hat %*% C_hat
  delta_hat <- setNames(drop(build_B_p(vartheta$phi, big_T) %*% data$Y_bar -
                               C_hat %*% data$W_bar), delta_names)
  residuals <- data$Y_c - data$W_c %*% t(K_hat)
  Sigma_hat <- .build_Sigma_internal(psi_hat, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
  pd_fitted <- .check_pd_components(vartheta, p, q, r, big_T,
                                    heteroskedastic, min_eig = 0,
                                    Sigma_override = Sigma_hat)
  mean_rank <- tryCatch(.qmle_x_gls_start(psi_hat, data, p, q, r, heteroskedastic,
      stationary = stationary, invertible = invertible),
                        error = function(e) NULL)
  if (is.null(mean_rank)) warning("The fitted centered mean design is rank deficient; Xi is not identified.",
                                  call. = FALSE)
  if (!isTRUE(pd_fitted$all_valid)) warning("The selected regressor QMLE solution violates covariance restrictions.",
                                           call. = FALSE)
  vartheta_hat <- setNames(.to_paper_parameters(psi_hat, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible),
                           .parameter_names(p, q, r, big_T, heteroskedastic, FALSE,
      stationary = stationary, invertible = invertible))
  names(psi_hat) <- .parameter_names(p, q, r, big_T, heteroskedastic, TRUE,
      stationary = stationary, invertible = invertible)
  vartheta$beta <- setNames(Xi$beta, paste0("beta_", seq_len(k)))
  vartheta$c <- Xi$c
  vartheta$A_mundlak <- Xi$A_mundlak
  vartheta$delta <- delta_hat
  optimization_hat <- c(psi_hat, Xi_hat)
  result <- list(
    coef = vartheta, psi_hat = psi_hat, vartheta_hat = vartheta_hat,
    Xi_hat = Xi_hat, delta_hat = delta_hat,
    parameter_hat = c(vartheta_hat, Xi_hat, delta_hat),
    optimization_hat = optimization_hat,
    estimator = "QMLE", obj_val = -best_val,
    converged = best_res$convergence == 0L && isTRUE(pd_fitted$all_valid) && !is.null(mean_rank),
    convergence_code = best_res$convergence, pd_components = pd_fitted,
    Sigma_hat = Sigma_hat, pd_tolerance = 0,
    Sigma_hat_construction = "internal optimization coordinates",
    C_hat = C_hat, K_hat = K_hat, residuals = residuals, fitted = Y - residuals,
    S_N = crossprod(residuals) / N, all_runs = all_runs,
    successful_starts = sum(vapply(all_runs, function(x) isTRUE(x$success), logical(1L))),
    failed_starts = sum(vapply(all_runs, function(x) !isTRUE(x$success), logical(1L))),
    hessian_type = hessian_type,
    bound_hits = .bound_hit_diagnostics(optimization_hat, joint_box),
    mean_design_rank = if (is.null(mean_rank)) NA_integer_ else mean_rank$rank,
    mean_design_rcond = if (is.null(mean_rank)) NA_real_ else mean_rank$rcond,
    parameterization = list(
      optimization = "dynamic reflection/physical blocks and covariance nuisance coordinates plus Xi; delta concentrated",
      reported = "paper coordinates (vartheta, Xi, delta)"
    ),
    has_regressors = TRUE, k = k, p = p, q = q, r = r, big_T = big_T, N = N
  )
  result <- .set_dynamic_domain(result, stationary, invertible, phi_bound, theta_bound)
  class(result) <- c("arma_panel_fit", "list")
  if (compute_se) result <- .stable_qmle_inference(
    result, data, p, q, r, heteroskedastic, hessian_type, with_x = TRUE
  ) else {
    result <- .qmle_x_empty_inference(result)
    result$inference_diagnostics <- list(first_failure = "not_computed")
  }
  if (warn_dynamic) .warn_dynamic_roots(vartheta, "QMLE")
  result
}

# stationary/invertible: TRUE uses reflection bounds [-(1-1e-4),1-1e-4].
# FALSE uses physical coefficients bounded by phi_bound/theta_bound (default 5).
# See dynamic_diagnostics for active reflection bounds and root diagnostics.
estimate_qmle <- function(Y, p, q, r,
                          heteroskedastic = TRUE,
                          n_starts = 20L,
                          maxit = 3000L,
                          factr = 1e4,
                          seed = NULL,
                          verbose = FALSE,
                          compute_se = TRUE,
                          warn_dynamic = TRUE,
                          phi_start = NULL,
                          hessian_type = c("obs", "info", "both"),
                          phi_bound = 5,
                          theta_bound = 5,
                          X = NULL,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  stationary <- .validate_dynamic_flag(stationary, "stationary")
  invertible <- .validate_dynamic_flag(invertible, "invertible")
  # Dispatch without changing the data or RNG; both branches share dynamic domains.
  if (is.null(X)) return(.estimate_qmle_no_x(
    Y = Y, p = p, q = q, r = r, heteroskedastic = heteroskedastic,
    n_starts = n_starts, maxit = maxit, factr = factr, seed = seed,
    verbose = verbose, compute_se = compute_se, warn_dynamic = warn_dynamic,
    phi_start = phi_start, hessian_type = hessian_type,
    phi_bound = phi_bound, theta_bound = theta_bound
  ,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      ))
  .estimate_qmle_with_x(
    Y = Y, p = p, q = q, r = r, heteroskedastic = heteroskedastic,
    n_starts = n_starts, maxit = maxit, factr = factr, seed = seed,
    verbose = verbose, compute_se = compute_se, warn_dynamic = warn_dynamic,
    phi_start = phi_start, hessian_type = hessian_type,
    phi_bound = phi_bound,
    theta_bound = theta_bound, X = X
  ,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      )
}

# Optimization profiles delta. Joint inference includes delta and uses stable
# factor coordinates before mapping the covariance to (vartheta, Xi, delta).

.qmle_x_inference_slices <- function(result) {
  n_v <- length(result$vartheta_hat)
  n_Xi <- length(result$Xi_hat)
  idx_v <- seq_len(n_v)
  idx_Xi <- n_v + seq_len(n_Xi)
  idx_delta <- n_v + n_Xi + seq_len(length(result$delta_hat))
  for (block in c("vartheta", "Xi", "delta")) {
    idx <- switch(block, vartheta = idx_v, Xi = idx_Xi, delta = idx_delta)
    result[[paste0("se_", block)]] <- result$se[idx]
    result[[paste0("ci_", block)]] <- result$ci[, idx, drop = FALSE]
    result[[paste0("V_hat_", block)]] <- result$V_hat[idx, idx, drop = FALSE]
  }
  result
}

.qmle_x_empty_inference <- function(result) {
  tmp <- .set_empty_inference(list(vartheta_hat = result$parameter_hat))
  for (nm in c("se", "ci", "V_hat")) result[[nm]] <- tmp[[nm]]
  .qmle_x_inference_slices(result)
}

.qmle_x_report_inference <- function(result, V_hat, N) {
  tmp <- .set_reported_inference(list(vartheta_hat = result$parameter_hat), V_hat, N)
  for (nm in c("se", "ci", "V_hat")) result[[nm]] <- tmp[[nm]]
  .qmle_x_inference_slices(result)
}

# Scale the coordinates for GMM numerical Jacobians. Retry with smaller steps
# if a derivative stencil crosses the admissible covariance domain.
.gmm_x_numerical_jacobian <- function(func, x, expected_dim) {
  parameter_scale <- pmax(1, abs(x))
  centered_func <- function(u) func(x + parameter_scale * u)
  for (step_scale in 0.25^(0:12)) {
    derivative_args <- list(eps = 1e-4 * step_scale, d = 1e-4,
                            zero.tol = 1e-12, r = 4L, v = 2)
    out <- tryCatch(
      numDeriv::jacobian(centered_func, rep(0, length(x)),
                         method = "Richardson", method.args = derivative_args),
      error = function(e) NULL
    )
    if (!is.null(out) && is.matrix(out) &&
        all(dim(out) == expected_dim) && all(is.finite(out))) {
      out <- sweep(out, 2L, parameter_scale, "/")
      return(list(value = out, step_scale = step_scale, success = TRUE))
    }
  }
  list(value = matrix(NA_real_, expected_dim[1L], expected_dim[2L]),
       step_scale = NA_real_, success = FALSE)
}

# ============================================================================ #
#  SECTION 12. CENTERED GMM WITH REGRESSORS [Paper Section 5.1]
# ============================================================================ #
# Supply X as in estimate_qmle().
#
# Both blocks use c_N=(N-1)/N:
#   g = [vech(Z_c'Z_c/(N*c_N)-Sigma); vec(Z_c'W_c/(N*c_N))],
#   Z_c=Y_c-W_c (Gamma C)', m=T(T+1)/2+k*T^2.
# This defines centered GMM and recovers delta from the sample mean relation;
# it is not exact profiling of an arbitrary full weighted uncentered criterion.
# Xi and vartheta are optimized simultaneously.  Only the intercept is removed.
# The current summary methods display phi/theta/beta and retain all_coefficients.
#
# Efficient weighting: one unrestricted reduced-form OLS regression supplies
# pilot residuals, then the entire moment covariance (including cross blocks)
# is inverted once. No preliminary constrained GMM or QMLE is performed.
# The actual weight remains fixed during all optimizer starts. Final residual
# moments supply robust inference, and recovered delta has joint influence SEs.
# Upsilon_weight_hat stores the pilot covariance; Upsilon_hat stores the final
# covariance used for the sandwich. They need not agree in finite samples.
# Redundant histories are not silently dropped: efficient weighting requires
# independent complete-history instruments and positive definite Upsilon.

.gmm_x_data <- function(Y, X) {
  data <- .qmle_x_data(Y, X)
  if (data$N < 2L) stop("GMM requires N >= 2 for the finite-sample correction.")
  data$c_N <- (data$N - 1) / data$N
  data$dot_T <- as.integer(data$big_T * (data$big_T + 1L) / 2L)
  data$m_count <- as.integer(data$dot_T + data$big_T^2 * data$k)
  data$moment_index <- which(lower.tri(matrix(FALSE, data$big_T, data$big_T),
                                       diag = TRUE), arr.ind = TRUE)
  data$moment_names <- c(
    paste0("cov_", data$moment_index[, "row"], "_", data$moment_index[, "col"]),
    unlist(lapply(seq_len(data$k * data$big_T), function(j)
      paste0("orthog_z", seq_len(data$big_T), "_w", j)), use.names = FALSE)
  )
  data$S_yy <- crossprod(data$Y_c) / data$N
  data$S_yw <- crossprod(data$Y_c, data$W_c) / data$N
  data$S_ww <- crossprod(data$W_c) / data$N
  data
}

.gmm_x_individual <- function(Z, data) {
  if (!is.matrix(Z) || any(dim(Z) != c(data$N, data$big_T)) || any(!is.finite(Z)))
    stop("Z must be a finite N x T residual matrix.")
  idx <- data$moment_index
  covariance_rows <- Z[, idx[, "row"], drop = FALSE] *
    Z[, idx[, "col"], drop = FALSE]
  cross_rows <- do.call(cbind, lapply(seq_len(ncol(data$W_c)), function(j)
    sweep(Z, 1L, data$W_c[, j], "*")))
  out <- cbind(covariance_rows, cross_rows) / data$c_N
  colnames(out) <- data$moment_names
  out
}

.gmm_x_moment_covariance <- function(Z, data) {
  rows <- .gmm_x_individual(Z, data)
  rows <- sweep(rows, 2L, colMeans(rows), "-")
  out <- crossprod(rows) / data$N
  (out + t(out)) / 2
}

.gmm_x_weight_data <- function(data) {
  if (data$N <= data$m_count)
    stop(sprintf(paste0("Efficient GMM with regressors requires N > m=%d ",
                        "(received N=%d) for empirical moment covariance inversion."),
                  data$m_count, data$N))
  scales <- sqrt(colSums(data$W_c^2))
  rank_W <- if (any(scales <= 0)) 0L else
    qr(sweep(data$W_c, 2L, scales, "/"), tol = 1e-10)$rank
  if (rank_W < ncol(data$W_c))
    stop(paste0("Efficient GMM requires independent centered complete-history ",
                "regressors. Redundant instruments require an explicit reduced basis."))
  pilot_residuals <- qr.resid(qr(data$W_c, tol = 1e-10), data$Y_c)
  Upsilon_hat <- .gmm_x_moment_covariance(pilot_residuals, data)
  R <- tryCatch(chol(Upsilon_hat), error = function(e) NULL)
  if (is.null(R)) {
    estimated_rank <- tryCatch(qr(Upsilon_hat, tol = 1e-10)$rank,
                               error = function(e) NA_integer_)
    stop(sprintf(paste0("The reduced-form OLS moment covariance is not positive ",
                        "definite (rank %s of %d); efficient weighting is unavailable."),
                  estimated_rank, data$m_count))
  }
  list(Upsilon_hat = Upsilon_hat, W_opt = chol2inv(R),
       pilot_residuals = pilot_residuals)
}

.gmm_x_validate_weight <- function(W, m) {
  if (!is.matrix(W) || !is.numeric(W) || any(dim(W) != c(m, m)))
    stop(sprintf("With regressors, W must be a numeric %d x %d matrix for the full augmented moment vector.", m, m))
  .validate_W(W, m)
}

.gmm_x_moments_from_model <- function(vartheta, Sigma, data, p, r) {
  if (is.null(Sigma) || any(!is.finite(Sigma)) ||
      is.null(tryCatch(chol(Sigma), error = function(e) NULL)))
    return(rep(NA_real_, data$m_count))
  C <- .qmle_x_C(vartheta, vartheta$Xi, p, r, data$k, data$big_T)
  K <- build_Gamma_p(vartheta$phi, data$big_T) %*% C
  # Algebraically identical to residual cross-products, but independent of N
  # in evaluation cost. The original data moments remain fixed, not residuals.
  KS_yw <- K %*% t(data$S_yw)
  S_z <- data$S_yy - KS_yw - t(KS_yw) + K %*% data$S_ww %*% t(K)
  S_z <- (S_z + t(S_z)) / 2
  S_zw <- data$S_yw - K %*% data$S_ww
  c(vech(S_z / data$c_N - Sigma), as.vector(S_zw / data$c_N))
}

.gmm_x_moments_internal <- function(optimization_vec, data, p, q, r,
                                    heteroskedastic,
                                  stationary = TRUE, invertible = TRUE) {
  d_v <- count_params(p, q, r, data$big_T, heteroskedastic)
  n_Xi <- .qmle_x_count_Xi(p, r, data$k, data$big_T)
  if (length(optimization_vec) != d_v + n_Xi || any(!is.finite(optimization_vec)))
    return(rep(NA_real_, data$m_count))
  tryCatch({
    psi <- optimization_vec[seq_len(d_v)]
    vartheta <- unpack_params(psi, p, q, r, data$big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
    vartheta$Xi <- .qmle_x_unpack_Xi(optimization_vec[d_v + seq_len(n_Xi)],
                                    p, r, data$k, data$big_T)
    Sigma <- .build_Sigma_internal(psi, p, q, r, data$big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
    .gmm_x_moments_from_model(vartheta, Sigma, data, p, r)
  }, error = function(e) rep(NA_real_, data$m_count))
}

.gmm_x_obj <- function(optimization_vec, data, W, p, q, r, heteroskedastic,
                                  stationary = TRUE, invertible = TRUE) {
  g <- .gmm_x_moments_internal(optimization_vec, data, p, q, r, heteroskedastic,
      stationary = stationary, invertible = invertible)
  if (any(!is.finite(g))) return(1e15)
  objective <- as.numeric(crossprod(g, W %*% g))
  if (is.finite(objective)) objective else 1e15
}

.estimate_gmm_with_x <- function(Y, p, q, r, W = NULL,
                                 heteroskedastic = TRUE, n_starts = 20L,
                                 maxit = 3000L, factr = 1e4, seed = NULL,
                                 verbose = FALSE, compute_se = TRUE,
                                 warn_dynamic = TRUE, efficient_weight = FALSE,
                                 .warn_T = TRUE, .gmm_data = NULL,
                                 phi_start = NULL,
                                 phi_bound = 5, theta_bound = 5, X,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  if (!is.null(seed)) set.seed(seed)
  if (!is.logical(efficient_weight) || length(efficient_weight) != 1L || is.na(efficient_weight))
    stop("efficient_weight must be TRUE or FALSE.")
  n_starts <- .validate_n_starts(n_starts)
  dims <- .validate_inputs(Y, p, q, r, heteroskedastic, warn_T = .warn_T)
  data <- if (is.null(.gmm_data)) .gmm_x_data(Y, X) else .gmm_data
  if (is.null(data$m_count) || data$N != nrow(Y) || data$big_T != ncol(Y))
    stop(".gmm_data must be the internally prepared regressor GMM data.")
  .qmle_x_validate_design(data, p, r)
  N <- dims$N
  big_T <- dims$big_T
  k <- data$k
  m <- data$m_count
  d_v <- count_params(p, q, r, big_T, heteroskedastic)
  n_Xi <- .qmle_x_count_Xi(p, r, k, big_T)
  d_opt <- d_v + n_Xi
  if (m < d_opt) stop("There are fewer augmented moments than optimized parameters.")
  if (is.null(W)) W <- diag(m)
  W <- .gmm_x_validate_weight(W, m)
  Upsilon_weight_hat <- data$Upsilon_weight_hat
  if (efficient_weight) {
    if (is.null(Upsilon_weight_hat) ||
        is.null(tryCatch(chol(Upsilon_weight_hat), error = function(e) NULL)) ||
        !isTRUE(all.equal(W, chol2inv(chol(Upsilon_weight_hat)), tolerance = 1e-8,
                         check.attributes = FALSE)))
      stop("Use estimate_gmm_efficient(..., X=X) to obtain a verified efficient regressor weight.")
  }
  box <- make_box_constraints(p, q, r, big_T, heteroskedastic,
                              phi_bound = phi_bound, theta_bound = theta_bound,
      stationary = stationary, invertible = invertible)
  joint_box <- box
  joint_box$lower <- c(box$lower, rep(-Inf, n_Xi))
  joint_box$upper <- c(box$upper, rep(Inf, n_Xi))
  Y_start <- if (!is.null(data$pilot_residuals)) data$pilot_residuals else
    qr.resid(qr(data$W_c, tol = 1e-10), data$Y_c)
  if (any(!is.finite(Y_start)) ||
      sum(Y_start^2) <= .Machine$double.eps * sum(data$Y_c^2)) Y_start <- data$Y_c
  starts <- .make_initial_starts(Y_start, box, p, q, r, big_T, heteroskedastic,
                                 n_starts, phi_start = phi_start,
                                   ndeps = ndeps, start_factr = start_factr
                                 )
  start_types <- attr(starts, "start_type")
  start_phis <- attr(starts, "initial_phi")
  anchor <- starts[[1L]]
  stable <- .make_stable_start(box, p, q, r, big_T, heteroskedastic)
  best_val <- Inf
  best_par <- NULL
  best_res <- list(convergence = 999L)
  all_runs <- vector("list", n_starts)
  max_start_attempts <- 25L
  for (start in seq_len(n_starts)) {
    success <- FALSE
    last_message <- NA_character_
    for (attempt in seq_len(max_start_attempts)) {
      replacement_error <- NULL
      if (attempt == 1L) {
        psi_init <- starts[[start]]
        current_type <- start_types[start]
        current_phi <- start_phis[[start]]
      } else if (attempt == max_start_attempts) {
        psi_init <- stable$psi
        current_type <- stable$type
        current_phi <- stable$phi
      } else {
        replacement <- tryCatch(.draw_replacement_start(
          start, anchor, Y_start, box, p, q, r, big_T, heteroskedastic,
          ndeps = ndeps, start_factr = start_factr
        ), error = function(e) { replacement_error <<- conditionMessage(e); NULL })
        psi_init <- if (is.null(replacement)) NULL else replacement$psi
        current_type <- if (is.null(replacement)) "replacement_generation_failed" else replacement$type
        current_phi <- if (is.null(replacement)) rep(NA_real_, p) else replacement$phi
      }
      Xi_init <- if (is.null(psi_init)) NULL else tryCatch(
        .qmle_x_gls_start(psi_init, data, p, q, r, heteroskedastic,
      stationary = stationary, invertible = invertible)$Xi,
        error = function(e) numeric(n_Xi)
      )
      res <- if (is.null(psi_init)) list(
        value = Inf, convergence = 999L, par = rep(NA_real_, d_opt),
        message = replacement_error
      ) else .run_box_optim(
        c(psi_init, Xi_init), function(v) .gmm_x_obj(v, data, W, p, q, r, heteroskedastic,
      stationary = stationary, invertible = invertible),
        joint_box$lower, joint_box$upper, maxit, factr,
        ndeps = ndeps
      )
      success <- is.finite(res$value) && res$convergence == 0L && res$value < 1e15
      if (!success) {
        last_message <- if (!is.null(res$message) && length(res$message) == 1L &&
                            !is.na(res$message) && nzchar(res$message)) res$message else
          "nonfinite, penalized, or nonconverged optimization"
        if (verbose) cat(sprintf("  GMM with X start %2d/%2d attempt %d failed.\n",
                                  start, n_starts, attempt))
        next
      }
      break
    }
    fitted_dynamic <- if (success) .decode_dynamics(
      if (p > 0L) res$par[seq_len(p)] else numeric(0),
      if (q > 0L) res$par[p + seq_len(q)] else numeric(0),
      stationary, invertible
    ) else list(phi = rep(NA_real_, p), theta = rep(NA_real_, q))
    phi_r <- round(fitted_dynamic$phi, 4)
    theta_r <- round(fitted_dynamic$theta, 4)
    all_runs[[start]] <- list(obj = if (success) res$value else Inf,
                             phi = phi_r, theta = theta_r,
                             conv = if (success) res$convergence else 999L,
                             attempts = attempt, start_type = current_type,
                             initial_phi = current_phi, success = success,
                             message = if (success) res$message else last_message)
    if (!success) {
      warning(sprintf("GMM with X start %d failed after %d attempts; continuing with remaining starts.",
                       start, max_start_attempts), call. = FALSE)
      next
    }
    if (verbose) cat(sprintf("  GMM with X start %2d/%2d obj=%11.4f conv=%d\n",
                              start, n_starts, res$value, res$convergence))
    if (res$value < best_val) { best_val <- res$value; best_par <- res$par; best_res <- res }
  }
  estimator <- if (efficient_weight) "Efficient GMM" else "GMM"
  Xi_names <- .qmle_x_Xi_names(p, r, k, big_T)
  delta_names <- paste0("delta_", seq_len(big_T))
  if (is.null(best_par)) {
    result <- .failed_optimizer_result(estimator, p, q, r, big_T, N,
                                        heteroskedastic, all_runs, Inf,
      stationary = stationary, invertible = invertible,
      phi_bound = phi_bound, theta_bound = theta_bound)
    result$Xi_hat <- setNames(rep(NA_real_, n_Xi), Xi_names)
    result$delta_hat <- setNames(rep(NA_real_, big_T), delta_names)
    result$parameter_hat <- c(result$vartheta_hat, result$Xi_hat, result$delta_hat)
    result$optimization_hat <- c(result$psi_hat, result$Xi_hat)
    result$coef$beta <- rep(NA_real_, k)
    result$coef$c <- matrix(NA_real_, p, k * big_T)
    result$coef$A_mundlak <- matrix(NA_real_, r, k * big_T)
    result$coef$delta <- result$delta_hat
    result$fitted <- result$residuals <- matrix(NA_real_, N, big_T)
    result$S_N <- matrix(NA_real_, big_T, big_T)
    result$C_hat <- result$K_hat <- matrix(NA_real_, big_T, k * big_T)
    result$moment_hat <- setNames(rep(NA_real_, m), data$moment_names)
    result$W <- W
    result$Upsilon_weight_hat <- Upsilon_weight_hat
    result$efficient_weight <- efficient_weight
    result$has_regressors <- TRUE
    result$k <- k
    result$inference_diagnostics <- list(first_failure = "optimization_failed")
    return(.gmm_x_empty_inference(result, data))
  }
  psi_hat <- best_par[seq_len(d_v)]
  vartheta <- unpack_params(psi_hat, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
  vartheta$Psi <- .Psi_from_lower(vartheta$Psi_lower, r)
  Xi_hat <- setNames(best_par[d_v + seq_len(n_Xi)], Xi_names)
  Xi <- .qmle_x_unpack_Xi(Xi_hat, p, r, k, big_T)
  C_hat <- .qmle_x_C(vartheta, Xi, p, r, k, big_T)
  K_hat <- build_Gamma_p(vartheta$phi, big_T) %*% C_hat
  delta_hat <- setNames(drop(build_B_p(vartheta$phi, big_T) %*% data$Y_bar -
                               C_hat %*% data$W_bar), delta_names)
  residuals <- data$Y_c - data$W_c %*% t(K_hat)
  Sigma_hat <- .build_Sigma_internal(psi_hat, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible)
  pd_fitted <- .check_pd_components(vartheta, p, q, r, big_T, heteroskedastic,
                                    min_eig = 0, Sigma_override = Sigma_hat)
  mean_rank <- tryCatch(.qmle_x_gls_start(psi_hat, data, p, q, r, heteroskedastic,
      stationary = stationary, invertible = invertible),
                        error = function(e) NULL)
  if (is.null(mean_rank)) warning("The fitted centered mean design is rank deficient; Xi is not identified.",
                                  call. = FALSE)
  if (!isTRUE(pd_fitted$all_valid)) warning("The selected regressor GMM solution violates covariance restrictions.",
                                           call. = FALSE)
  vartheta_hat <- setNames(.to_paper_parameters(psi_hat, p, q, r, big_T, heteroskedastic,
      stationary = stationary, invertible = invertible),
                           .parameter_names(p, q, r, big_T, heteroskedastic, FALSE,
      stationary = stationary, invertible = invertible))
  names(psi_hat) <- .parameter_names(p, q, r, big_T, heteroskedastic, TRUE,
      stationary = stationary, invertible = invertible)
  vartheta$beta <- setNames(Xi$beta, paste0("beta_", seq_len(k)))
  vartheta$c <- Xi$c
  vartheta$A_mundlak <- Xi$A_mundlak
  vartheta$delta <- delta_hat
  optimization_hat <- c(psi_hat, Xi_hat)
  result <- list(
    coef = vartheta, psi_hat = psi_hat, vartheta_hat = vartheta_hat,
    Xi_hat = Xi_hat, delta_hat = delta_hat,
    parameter_hat = c(vartheta_hat, Xi_hat, delta_hat),
    optimization_hat = optimization_hat, estimator = estimator, obj_val = best_val,
    converged = best_res$convergence == 0L && isTRUE(pd_fitted$all_valid) && !is.null(mean_rank),
    convergence_code = best_res$convergence, pd_components = pd_fitted,
    Sigma_hat = Sigma_hat, pd_tolerance = 0,
    Sigma_hat_construction = "internal optimization coordinates",
    C_hat = C_hat, K_hat = K_hat, residuals = residuals, fitted = Y - residuals,
    S_N = crossprod(residuals) / N,
    moment_hat = setNames(.gmm_x_moments_internal(best_par, data, p, q, r, heteroskedastic,
      stationary = stationary, invertible = invertible),
                           data$moment_names),
    W = W, Upsilon_weight_hat = Upsilon_weight_hat, efficient_weight = efficient_weight,
    all_runs = all_runs,
    successful_starts = sum(vapply(all_runs, function(x) isTRUE(x$success), logical(1L))),
    failed_starts = sum(vapply(all_runs, function(x) !isTRUE(x$success), logical(1L))),
    bound_hits = .bound_hit_diagnostics(optimization_hat, joint_box),
    mean_design_rank = if (is.null(mean_rank)) NA_integer_ else mean_rank$rank,
    mean_design_rcond = if (is.null(mean_rank)) NA_real_ else mean_rank$rcond,
    parameterization = list(optimization = "dynamic reflection/physical blocks and covariance nuisance coordinates plus Xi; centered mean relation",
                             reported = "paper coordinates (vartheta, Xi, delta)"),
    has_regressors = TRUE, k = k, p = p, q = q, r = r, big_T = big_T, N = N
  )
  result <- .set_dynamic_domain(result, stationary, invertible, phi_bound, theta_bound)
  class(result) <- c("arma_panel_fit", "list")
  if (compute_se) result <- .gmm_x_inference(result, data, p, q, r, heteroskedastic) else {
    result <- .gmm_x_empty_inference(result, data)
    result$inference_diagnostics <- list(first_failure = "not_computed")
  }
  if (warn_dynamic) .warn_dynamic_roots(vartheta, estimator)
  result
}

# stationary/invertible: TRUE uses reflection bounds [-(1-1e-4),1-1e-4].
# FALSE uses physical coefficients bounded by phi_bound/theta_bound (default 5).
# See dynamic_diagnostics for active reflection bounds and root diagnostics.
estimate_gmm <- function(Y, p, q, r,
                         W = NULL,
                         heteroskedastic = TRUE,
                         n_starts = 20L,
                         maxit = 3000L,
                         factr = 1e4,
                         seed = NULL,
                         verbose = FALSE,
                         compute_se = TRUE,
                         warn_dynamic = TRUE,
                         efficient_weight = FALSE,
                         .warn_T = TRUE,
                         .gmm_data = NULL,
                         phi_start = NULL,
                         phi_bound = 5,
                         theta_bound = 5,
                         X = NULL,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  stationary <- .validate_dynamic_flag(stationary, "stationary")
  invertible <- .validate_dynamic_flag(invertible, "invertible")
  if (is.null(X)) return(.estimate_gmm_no_x(
    Y = Y, p = p, q = q, r = r, W = W, heteroskedastic = heteroskedastic,
    n_starts = n_starts, maxit = maxit, factr = factr, seed = seed,
    verbose = verbose, compute_se = compute_se, warn_dynamic = warn_dynamic,
    efficient_weight = efficient_weight, .warn_T = .warn_T, .gmm_data = .gmm_data,
    phi_start = phi_start,
    phi_bound = phi_bound, theta_bound = theta_bound
  ,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      ))
  .estimate_gmm_with_x(
    Y = Y, p = p, q = q, r = r, W = W, heteroskedastic = heteroskedastic,
    n_starts = n_starts, maxit = maxit, factr = factr, seed = seed,
    verbose = verbose, compute_se = compute_se, warn_dynamic = warn_dynamic,
    efficient_weight = efficient_weight, .warn_T = .warn_T, .gmm_data = .gmm_data,
    phi_start = phi_start,
    phi_bound = phi_bound, theta_bound = theta_bound, X = X
  ,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      )
}

# stationary/invertible: TRUE uses reflection bounds [-(1-1e-4),1-1e-4].
# FALSE uses physical coefficients bounded by phi_bound/theta_bound (default 5).
# See dynamic_diagnostics for active reflection bounds and root diagnostics.
estimate_gmm_efficient <- function(Y, p, q, r,
                                   heteroskedastic = TRUE,
                                   n_starts = 20L,
                                   maxit = 3000L,
                                   factr = 1e4,
                                   seed = NULL,
                                   verbose = FALSE,
                                   compute_se = TRUE,
                                   warn_dynamic = TRUE,
                                   phi_start = NULL,
                                   phi_bound = 5,
                                   theta_bound = 5,
                                   X = NULL,
                                  stationary = TRUE, invertible = TRUE,
                                  ndeps = 1e-5, start_factr = 1e4) {
  ndeps <- .validate_positive_optimizer_control(ndeps, "ndeps")
  start_factr <- .validate_positive_optimizer_control(start_factr, "start_factr")
  stationary <- .validate_dynamic_flag(stationary, "stationary")
  invertible <- .validate_dynamic_flag(invertible, "invertible")
  if (is.null(X)) return(.estimate_gmm_efficient_no_x(
    Y = Y, p = p, q = q, r = r, heteroskedastic = heteroskedastic,
    n_starts = n_starts, maxit = maxit, factr = factr, seed = seed,
    verbose = verbose, compute_se = compute_se, warn_dynamic = warn_dynamic,
    phi_start = phi_start,
    phi_bound = phi_bound, theta_bound = theta_bound
  ,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      ))
  phi_bound <- .validate_coefficient_bound(phi_bound, "phi_bound")
  theta_bound <- .validate_coefficient_bound(theta_bound, "theta_bound")
  n_starts <- .validate_n_starts(n_starts)
  .validate_inputs(Y, p, q, r, heteroskedastic)
  data <- .gmm_x_data(Y, X)
  .qmle_x_validate_design(data, p, r)
  weight_data <- .gmm_x_weight_data(data)
  data$Upsilon_weight_hat <- weight_data$Upsilon_hat
  data$pilot_residuals <- weight_data$pilot_residuals
  if (verbose) cat("=== Efficient centered GMM with OLS-based augmented weight ===\n")
  result <- .estimate_gmm_with_x(
    Y = Y, p = p, q = q, r = r, W = weight_data$W_opt,
    heteroskedastic = heteroskedastic, n_starts = n_starts,
    maxit = maxit, factr = factr, seed = seed, verbose = verbose,
    compute_se = compute_se, warn_dynamic = FALSE, efficient_weight = TRUE,
    .warn_T = FALSE, .gmm_data = data, phi_start = phi_start,
    phi_bound = phi_bound,
    theta_bound = theta_bound, X = X
  ,
      stationary = stationary, invertible = invertible,
        ndeps = ndeps, start_factr = start_factr
      )
  result$W_opt <- weight_data$W_opt
  result$Upsilon_weight_hat <- weight_data$Upsilon_hat
  result$pilot_residuals <- weight_data$pilot_residuals
  if (warn_dynamic) .warn_dynamic_roots(result$coef, "Efficient GMM")
  result
}

# Inference for GMM with covariates and an unrestricted recovered intercept.
# The sample moments are centered covariance/orthogonality moments.  The full
# reported variance includes the sampling variation of Y_bar and W_bar.

.gmm_x_empty_inference <- function(result, data) {
  result <- .qmle_x_empty_inference(result)
  outer_names <- names(c(result$vartheta_hat, result$Xi_hat))
  d <- length(outer_names)
  m <- data$m_count
  moment_names <- if (!is.null(data$moment_names)) data$moment_names else
    paste0("moment_", seq_len(m))
  result$G <- matrix(NA_real_, m, d,
                    dimnames = list(moment_names, outer_names))
  result$moment_jacobian <- result$G
  result$Upsilon_hat <- matrix(NA_real_, m, m,
                              dimnames = list(moment_names, moment_names))
  result$H_hat <- result$J_hat <- matrix(NA_real_, d, d,
                                        dimnames = list(outer_names, outer_names))
  result$A_hat <- matrix(NA_real_, d, m,
                        dimnames = list(outer_names, moment_names))
  result$D_delta_eta <- matrix(NA_real_, data$big_T, d,
                              dimnames = list(names(result$delta_hat), outer_names))
  result$rank_G <- result$rank_weighted_G <- NA_integer_
  result$weighted_G_singular_values <- rep(NA_real_, d)
  result$rcond_GWG <- NA_real_
  result$J_stat <- NA_real_
  result$df_overid <- NA_integer_
  result$J_pval <- NA_real_
  result$hansen_test <- NULL
  result
}

.gmm_x_inference <- function(result, data, p, q, r, heteroskedastic) {
  .stable_gmm_inference(result, data, p, q, r, heteroskedastic, with_x = TRUE)
}


# ============================================================================ #
# EXACT ANALYTICAL QMLE DERIVATIVES AND INFERENCE
# ============================================================================ #

# A diagonal congruence changes units only: H^{-1}=D (D H D)^{-1} D.
# No eigenvalue clipping, ridge, pseudoinverse or numerical derivative fallback.
# Preserve existing policy: regressor bread and information bread require PD;
# covariance-only observed bread accepts a finite invertible indefinite matrix.
.qmle_analytic_bread <- function(H, d, require_pd = FALSE) {
  finite <- is.matrix(H) && identical(dim(H), c(d, d)) && all(is.finite(H))
  H <- if (finite) (H + t(H)) / 2 else matrix(NA_real_, d, d)
  safe <- function(expr, default) tryCatch(expr, error = function(e) default)
  eigenvalues <- if (finite) safe(eigen(H, symmetric = TRUE, only.values = TRUE)$values,
                                  rep(NA_real_, d)) else rep(NA_real_, d)
  rank <- if (finite) safe(qr(H, tol = 1e-10)$rank, NA_integer_) else NA_integer_
  rc <- if (finite) safe(rcond(H), NA_real_) else NA_real_
  scale <- rep(1, d); Hs <- H
  if (finite) {
    magnitude <- abs(diag(H))
    bad <- !is.finite(magnitude) | magnitude == 0
    if (any(bad)) magnitude[bad] <- apply(abs(H), 1L, max)[bad]
    magnitude[magnitude == 0] <- 1
    scale <- 1 / sqrt(magnitude)
    Hs <- sweep(sweep(H, 1L, scale, "*"), 2L, scale, "*")
  }
  scaled_finite <- finite && all(is.finite(Hs))
  Rs <- if (scaled_finite) safe(chol(Hs), NULL) else NULL
  pd <- !is.null(Rs) && all(is.finite(Rs)) && all(diag(Rs) > 0)
  if (!pd) Rs <- NULL
  inverse_scaled <- if (pd) safe(chol2inv(Rs), NULL) else
    if (scaled_finite && !require_pd) safe(solve(Hs), NULL) else NULL
  inverse <- if (!is.null(inverse_scaled))
    sweep(sweep(inverse_scaled, 1L, scale, "*"), 2L, scale, "*") else NULL
  ok <- !is.null(inverse) && all(is.finite(inverse))
  solve_rhs <- if (ok) function(rhs) {
    rhs_scaled <- sweep(rhs, 1L, scale, "*")
    answer <- if (pd) backsolve(Rs, forwardsolve(t(Rs), rhs_scaled)) else
      solve(Hs, rhs_scaled)
    sweep(answer, 1L, scale, "*")
  } else NULL
  list(H = H, finite = finite, rank = rank, rcond = rc,
       eigenvalues = eigenvalues, positive_definite = pd,
       inverse = if (ok) (inverse + t(inverse)) / 2 else matrix(NA_real_, d, d),
       inverse_method = if (ok) if (pd) "scaled_chol" else "scaled_solve" else "failed",
       scaled_finite = scaled_finite,
       scaled_rcond = if (scaled_finite) safe(rcond(Hs), NA_real_) else NA_real_,
       scaled_rank = if (scaled_finite) safe(qr(Hs, tol = 1e-10)$rank, NA_integer_) else NA_integer_,
       scale = scale, solve_rhs = solve_rhs,
       inverse_backward_error = if (ok) norm(H %*% inverse - diag(d), "I") /
         (norm(H, "I") * norm(inverse, "I") + 1) else NA_real_)
}

.qmle_analytic_sandwich <- function(bread, scores, N) {
  d <- ncol(scores)
  if (is.null(bread$solve_rhs) || any(!is.finite(scores))) return(matrix(NA_real_, d, d))
  influence <- tryCatch(bread$solve_rhs(t(scores)), error = function(e) NULL)
  if (is.null(influence) || any(!is.finite(influence))) return(matrix(NA_real_, d, d))
  # Nonnegative diagonal by construction; reporting subsequently divides by N.
  tcrossprod(influence) / N
}




# Stable optimizer seed for one model; independent of candidate enumeration.
model_optimizer_seed <- function(seed, p, q, r) {
  max_integer <- .Machine$integer.max
  valid_integer <- function(x) is.numeric(x) && is.null(dim(x)) &&
    length(x) == 1L && is.finite(x) && x == floor(x) &&
    x >= 0 && x <= max_integer
  if (!valid_integer(seed))
    stop("seed must be one integer in 0..2147483647.", call. = FALSE)
  orders <- list(p = p, q = q, r = r)
  if (!all(vapply(orders, valid_integer, logical(1L))))
    stop("p, q and r must be nonnegative integers no larger than 2147483647.", call. = FALSE)
  max_exact <- 2^53 - 1
  pair <- function(a, b) {
    a <- as.double(a); b <- as.double(b)
    if (a > max_exact - b)
      stop("Model orders exceed the exact seed-pairing range.", call. = FALSE)
    total <- a + b
    # Divide an even factor first, avoiding an inexact intermediate product.
    left <- if (total %% 2 == 0) total / 2 else total
    right <- if (total %% 2 == 0) total + 1 else (total + 1) / 2
    if (right > 0 && left > (max_exact - b) / right)
      stop("Model orders exceed the exact seed-pairing range.", call. = FALSE)
    left * right + b
  }
  code <- pair(pair(p, q), r) + 1
  if (!is.finite(code) || code > floor((max_exact - as.double(seed)) / 100000))
    stop("Model orders exceed the exact seed-offset range.", call. = FALSE)
  as.integer((as.double(seed) + 100000 * code) %% 2147483647)
}


# Stable inference module: stable_derivatives.R
# Analytical covariance and mean derivatives in stable factor coordinates.
# The coordinate order is physical phi/theta, innovation variances,
# vech(L), vec(A2), and physical Pi; A = rbind(L, A2).
# For regressors, the former Mundlak block is R = L^{-1} M, so F M = A R.
# These functions do not decode optimizer reflections/Schur factors or alter
# the optimization problem. All nuisance coordinates are independent here.

.stable_covariance <- function(vec, p, q, r, big_T,
                               heteroskedastic = TRUE, second = TRUE) {
  d <- count_params(p, q, r, big_T, heteroskedastic)
  if (big_T <= p || r > big_T || p < 0L || q < 0L || r < 0L)
    stop("Invalid dimensions for stable covariance derivatives.")
  if (length(vec) < d || any(!is.finite(vec[seq_len(d)])))
    stop("Incorrect or nonfinite stable covariance parameter vector.")
  x <- vec[seq_len(d)]
  T <- as.integer(big_T)
  zero <- matrix(0, T, T)
  blanks <- function(nr, nc) replicate(d, matrix(0, nr, nc), simplify = FALSE)
  type <- rep("", d)
  pos <- 0L

  B <- diag(T)
  B1 <- blanks(T, T)
  phi_index <- seq_len(p)
  rows <- seq.int(p + 1L, T)
  if (p > 0L) for (j in seq_len(p)) {
    a <- pos <- pos + 1L
    type[a] <- "phi"
    B[cbind(rows, rows - j)] <- -x[a]
    B1[[a]][cbind(rows, rows - j)] <- -1
  }

  # Innovation epochs index columns 1-q,...,T, including the q presample slots.
  M <- matrix(0, T - p, T + q)
  M1 <- blanks(T - p, T + q)
  theta_index <- p + seq_len(q)
  if (q > 0L) for (j in seq_len(q)) {
    a <- pos <- pos + 1L
    type[a] <- "theta"
    ij <- cbind(rows - p, rows - j + q)
    M[ij] <- x[a]
    M1[[a]][ij] <- 1
  }
  M[cbind(rows - p, rows + q)] <- 1

  ns <- if (heteroskedastic) T - p + q else 1L
  sigma_index <- pos + seq_len(ns)
  D <- matrix(0, T + q, T + q)
  D1 <- blanks(T + q, T + q)
  if (heteroskedastic) for (j in seq_len(ns)) {
    a <- pos <- pos + 1L
    type[a] <- "sigma"
    D[p + j, p + j] <- x[a]
    D1[[a]][p + j, p + j] <- 1
  } else {
    a <- pos <- pos + 1L
    type[a] <- "sigma"
    ix <- seq.int(p + 1L, T + q)
    D[cbind(ix, ix)] <- x[a]
    D1[[a]][cbind(ix, ix)] <- 1
  }

  # A is affine in vech(L), vec(A2); its first r rows are lower triangular.
  A <- matrix(0, T, r)
  A_index <- matrix(0L, T, r)
  A_row <- A_factor <- integer(d)
  if (r > 0L) {
    for (j in seq_len(r)) for (i in seq.int(j, r)) {
      a <- pos <- pos + 1L
      type[a] <- "factor"
      A[i, j] <- x[a]
      A_index[i, j] <- a
      A_row[a] <- i
      A_factor[a] <- j
    }
    if (T > r) for (j in seq_len(r)) for (i in seq.int(r + 1L, T)) {
      a <- pos <- pos + 1L
      type[a] <- "factor"
      A[i, j] <- x[a]
      A_index[i, j] <- a
      A_row[a] <- i
      A_factor[a] <- j
    }
  }

  # Pi: initial variances, lexicographic initial covariances, then bridges.
  pi_start <- pos
  E <- matrix(0, p, p)
  E1 <- blanks(p, p)
  Q <- matrix(0, p, T + q)
  Q1 <- blanks(p, T + q)
  if (p > 0L) {
    for (i in seq_len(p)) {
      a <- pos <- pos + 1L
      type[a] <- "eta"
      E[i, i] <- x[a]
      E1[[a]][i, i] <- 1
    }
    if (p > 1L) for (i in seq_len(p - 1L)) for (j in seq.int(i + 1L, p)) {
      a <- pos <- pos + 1L
      type[a] <- "eta"
      E[i, j] <- E[j, i] <- x[a]
      E1[[a]][i, j] <- E1[[a]][j, i] <- 1
    }
    if (q > 0L) for (s in seq.int(max(1L, p + 1L - q), p))
      for (k in seq.int(p + 1L - q, s)) {
        a <- pos <- pos + 1L
        type[a] <- "bridge"
        Q[s, k + q] <- x[a]
        Q1[[a]][s, k + q] <- 1
      }
  }
  if (pos != d || any(!nzchar(type)))
    stop("Stable covariance parameter-order mismatch.")

  assemble <- function(top, bridge, bottom) {
    z <- zero
    if (p > 0L) {
      z[seq_len(p), seq_len(p)] <- top
      z[seq_len(p), rows] <- bridge
      z[rows, seq_len(p)] <- t(bridge)
    }
    z[rows, rows] <- bottom
    z
  }
  C <- tcrossprod(A) + assemble(E, Q %*% t(M), M %*% D %*% t(M))
  C1 <- lapply(seq_len(d), function(a) {
    if (type[a] == "phi") return(zero)
    if (type[a] == "factor") {
      z <- zero
      z[A_row[a], ] <- A[, A_factor[a]]
      z[, A_row[a]] <- z[, A_row[a]] + A[, A_factor[a]]
      return(z)
    }
    assemble(E1[[a]], Q1[[a]] %*% t(M) + Q %*% t(M1[[a]]),
      M1[[a]] %*% D %*% t(M) + M %*% D1[[a]] %*% t(M) +
        M %*% D %*% t(M1[[a]]))
  })

  C2 <- function(a, b) {
    z <- zero
    if (type[a] == "factor" && A_factor[a] == A_factor[b]) {
      z[A_row[a], A_row[b]] <- z[A_row[a], A_row[b]] + 1
      z[A_row[b], A_row[a]] <- z[A_row[b], A_row[a]] + 1
      return(z)
    }
    # The idiosyncratic covariance is quadratic in theta, and bilinear in
    # theta x sigma2 and theta x initial-innovation bridges.
    if (!any(c(type[a], type[b]) == "theta")) return(z)
    ma <- M1[[a]] %*% D1[[b]] %*% t(M) +
      M1[[a]] %*% D %*% t(M1[[b]]) + M1[[b]] %*% D1[[a]] %*% t(M) +
      M %*% D1[[a]] %*% t(M1[[b]]) + M1[[b]] %*% D %*% t(M1[[a]]) +
      M %*% D1[[b]] %*% t(M1[[a]])
    bridge <- Q1[[a]] %*% t(M1[[b]]) + Q1[[b]] %*% t(M1[[a]])
    assemble(matrix(0, p, p), bridge, ma)
  }

  # solve() also supports complex-step verification; forwardsolve() drops
  # imaginary perturbations. Real production evaluations use the triangular solve.
  Gamma <- if (is.complex(B)) solve(B) else forwardsolve(B, diag(T))
  Gamma_first <- lapply(B1, function(z) -Gamma %*% z %*% Gamma)
  valid_index <- function(a) length(a) == 1L && is.numeric(a) &&
    is.finite(a) && a >= 1L && a <= d && a == floor(a)
  Gamma_second <- function(a, b) {
    if (!valid_index(a) || !valid_index(b))
      stop("Invalid stable covariance derivative index.")
    if (a > p || b > p) return(zero)
    Gamma %*% B1[[a]] %*% Gamma %*% B1[[b]] %*% Gamma +
      Gamma %*% B1[[b]] %*% Gamma %*% B1[[a]] %*% Gamma
  }
  Sigma <- Gamma %*% C %*% t(Gamma)
  first <- lapply(seq_len(d), function(a) {
    if (a <= p) {
      z <- Gamma_first[[a]] %*% C %*% t(Gamma)
      z + t(z)
    } else Gamma %*% C1[[a]] %*% t(Gamma)
  })
  S2 <- function(a, b) {
    if (!valid_index(a) || !valid_index(b))
      stop("Invalid stable covariance derivative index.")
    if (a <= p && b <= p) {
      gab <- Gamma_second(a, b)
      z <- gab %*% C %*% t(Gamma)
      ab <- Gamma_first[[a]] %*% C %*% t(Gamma_first[[b]])
      return(z + t(z) + ab + t(ab))
    }
    if (a <= p || b <= p) {
      j <- if (a <= p) a else b
      k <- if (a <= p) b else a
      z <- Gamma_first[[j]] %*% C1[[k]] %*% t(Gamma)
      return(z + t(z))
    }
    Gamma %*% C2(a, b) %*% t(Gamma)
  }
  list(Sigma = Sigma, first = first, second = if (second) S2 else NULL,
       Gamma = Gamma, Gamma_first = Gamma_first, Gamma_second = Gamma_second,
       A = A, L = A[seq_len(r), , drop = FALSE], A_index = A_index,
       indices = list(phi = phi_index, theta = theta_index,
                      sigma2 = sigma_index, factor = which(type == "factor"),
                      Pi = if (pos > pi_start) seq.int(pi_start + 1L, pos) else integer(0)))
}


.stable_mean <- function(vec, p, q, r, big_T, k,
                         heteroskedastic = TRUE, covariance = NULL,
                         second = TRUE) {
  d_v <- count_params(p, q, r, big_T, heteroskedastic)
  T <- as.integer(big_T)
  nw <- k * T
  n_Xi <- k + (p + r) * nw
  d <- d_v + n_Xi + T
  if (T <= p || r > T || k < 1L)
    stop("Invalid dimensions for stable regressor mean derivatives.")
  if (length(vec) != d || any(!is.finite(vec)))
    stop("Incorrect or nonfinite stable joint parameter vector.")
  if (is.null(covariance)) covariance <- .stable_covariance(
    vec, p, q, r, T, heteroskedastic, second = FALSE)
  Gamma <- covariance$Gamma
  Gamma_first <- covariance$Gamma_first
  Gamma_second <- covariance$Gamma_second
  A <- covariance$A
  A_index <- covariance$A_index
  if (!identical(dim(Gamma), c(T, T)) || length(Gamma_first) < p ||
      !is.function(Gamma_second) || !identical(dim(A), c(T, as.integer(r))) ||
      !identical(dim(A_index), c(T, as.integer(r))))
    stop("Incomplete stable covariance derivatives supplied to mean derivatives.")

  zero_T <- matrix(0, T, T)
  zero_map <- matrix(0, T, 1L + nw)
  beta_index <- d_v + seq_len(k)
  c_index <- if (p > 0L)
    matrix(d_v + k + seq_len(p * nw), p, nw, byrow = TRUE) else
      matrix(integer(0), 0L, nw)
  R_index <- if (r > 0L)
    matrix(d_v + k + p * nw + seq_len(r * nw), r, nw) else
      matrix(integer(0), 0L, nw)
  delta_index <- d_v + n_Xi + seq_len(T)
  R <- if (r > 0L) matrix(vec[R_index], r, nw) else matrix(0, 0L, nw)
  M <- zero_map
  M[, 1L] <- vec[delta_index]
  if (p > 0L) M[seq_len(p), -1L] <- matrix(vec[c_index], p, nw)
  for (tt in seq.int(p + 1L, T))
    M[tt, 1L + (tt - 1L) * k + seq_len(k)] <- vec[beta_index]
  if (r > 0L) M[, -1L] <- M[, -1L, drop = FALSE] + A %*% R

  M_first <- replicate(d, zero_map, simplify = FALSE)
  for (j in seq_len(k)) {
    E <- zero_map
    dates <- seq.int(p + 1L, T)
    E[cbind(dates, 1L + (dates - 1L) * k + j)] <- 1
    M_first[[beta_index[j]]] <- E
  }
  if (p > 0L) for (s in seq_len(p)) for (w in seq_len(nw)) {
    E <- zero_map
    E[s, 1L + w] <- 1
    M_first[[c_index[s, w]]] <- E
  }
  A_row <- A_factor <- R_column <- R_factor <- integer(d)
  if (r > 0L) {
    for (j in seq_len(r)) {
      for (w in seq_len(nw)) {
        E <- zero_map
        E[, 1L + w] <- A[, j]
        M_first[[R_index[j, w]]] <- E
        R_column[R_index[j, w]] <- 1L + w
        R_factor[R_index[j, w]] <- j
      }
      for (tt in seq_len(T)) if (A_index[tt, j] > 0L) {
        a <- A_index[tt, j]
        E <- zero_map
        E[tt, -1L] <- R[j, ]
        M_first[[a]] <- E
        A_row[a] <- tt
        A_factor[a] <- j
      }
    }
  }
  for (tt in seq_len(T)) {
    E <- zero_map
    E[tt, 1L] <- 1
    M_first[[delta_index[tt]]] <- E
  }
  M_second <- function(a, b) {
    E <- zero_map
    if (A_factor[a] > 0L && A_factor[a] == R_factor[b])
      E[A_row[a], R_column[b]] <- 1
    if (A_factor[b] > 0L && A_factor[b] == R_factor[a])
      E[A_row[b], R_column[a]] <- 1
    E
  }
  first <- lapply(seq_len(d), function(a) {
    if (a <= p) Gamma_first[[a]] %*% M else Gamma %*% M_first[[a]]
  })
  valid_index <- function(a) length(a) == 1L && is.numeric(a) &&
    is.finite(a) && a >= 1L && a <= d && a == floor(a)
  S2 <- function(a, b) {
    if (!valid_index(a) || !valid_index(b))
      stop("Invalid stable mean derivative index.")
    if (a <= p && b <= p) return(Gamma_second(a, b) %*% M)
    if (a <= p || b <= p) {
      j <- if (a <= p) a else b
      h <- if (a <= p) b else a
      return(Gamma_first[[j]] %*% M_first[[h]])
    }
    Gamma %*% M_second(a, b)
  }
  list(map = Gamma %*% M, first = first, second = if (second) S2 else NULL,
       Gamma = Gamma, d_v = d_v, n_Xi = n_Xi, A = A, R = R,
       indices = list(beta = beta_index, c = c_index, R = R_index,
                      delta = delta_index, factor = A_index))
}



# Stable inference module: stable_reporting_qmle.R
# Inference-only stable factor chart; point objectives and optimizer coordinates
# are not changed. Physical dynamics and Pi remain independent coordinates.
.stable_layout <- function(p, q, r, big_T, heteroskedastic = TRUE, k = 0L) {
  ns <- if (heteroskedastic) big_T - p + q else 1L
  pre <- p + q + ns; nl <- r * (r + 1L) / 2L
  dv <- count_params(p, q, r, big_T, heteroskedastic)
  np <- if (k > 0L) .qmle_x_count_Xi(p, r, k, big_T) else 0L
  list(pre = pre, nl = nl, dv = dv, nXi = np,
       L = if (nl) pre + seq_len(nl) else integer(),
       A2 = if (r && big_T > r) pre + nl + seq_len((big_T-r)*r) else integer(),
       M = if (k > 0L && r > 0L) dv+k+p*k*big_T+seq_len(r*k*big_T) else integer(),
       full = dv+np+if (k > 0L) big_T else 0L)
}

.stable_inference_parameters <- function(result, with_x = FALSE) {
  p <- result$p; q <- result$q; r <- result$r; T <- result$big_T
  hetero <- if (is.null(result$heteroskedastic)) TRUE else result$heteroskedastic
  k <- if (with_x) result$k else 0L
  ix <- .stable_layout(p,q,r,T,hetero,k)
  paper <- if (with_x) result$parameter_hat else result$vartheta_hat
  if (length(paper) != ix$full || any(!is.finite(paper)))
    stop("Missing or nonfinite fitted parameters for stable inference.")
  ans <- as.numeric(paper); names(ans) <- names(paper)
  L <- matrix(0,r,r)
  internal <- result$psi_hat
  if (r > 0L) {
    if (length(internal) >= ix$dv && all(is.finite(internal[c(ix$L,ix$A2)]))) {
      ans[c(ix$L,ix$A2)] <- internal[c(ix$L,ix$A2)]
      L[lower.tri(L,diag=TRUE)] <- ans[ix$L]
    } else {
      # Legacy saved parameter vectors may not retain optimizer coordinates.
      P <- matrix(0,r,r); P[lower.tri(P,diag=TRUE)] <- paper[ix$L]
      P <- P+t(P)-diag(diag(P),nrow=r)
      L <- t(chol(P))
      F2 <- matrix(paper[ix$A2],T-r,r)
      ans[ix$L] <- L[lower.tri(L,diag=TRUE)]
      ans[ix$A2] <- as.vector(F2%*%L)
    }
    if (any(diag(L) <= 0)) stop("Nonpositive triangular factor diagonal.")
    ij <- which(lower.tri(L,diag=TRUE),arr.ind=TRUE)
    names(ans)[ix$L] <- paste0("L_Psi_",ij[,1],"_",ij[,2])
    if(length(ix$A2)) names(ans)[ix$A2] <- unlist(lapply(seq_len(r),function(j)paste0("A_",seq.int(r+1L,T),"_",j)))
    if (with_x) {
      M <- matrix(paper[ix$M],r,k*T)
      ans[ix$M] <- as.vector(forwardsolve(L,M))
      names(ans)[ix$M] <- unlist(lapply(seq_len(k*T),function(j)paste0("R_mundlak_",seq_len(r),"_",j)))
    }
  }
  ans
}

.stable_reporting_map <- function(vec, p, q, r, big_T,
                                  heteroskedastic = TRUE, k = 0L) {
  ix <- .stable_layout(p,q,r,big_T,heteroskedastic,k)
  if (length(vec) != ix$full) stop("Incorrect stable reporting dimension.")
  value <- as.numeric(vec); J <- diag(ix$full)
  if (r > 0L) {
    L <- matrix(0,r,r); L[lower.tri(L,diag=TRUE)] <- vec[ix$L]
    if (any(!is.finite(L)) || any(diag(L) <= 0)) stop("Invalid reporting factor.")
    A2 <- matrix(vec[ix$A2],big_T-r,r)
    right_solve <- function(B) t(backsolve(t(L),t(B)))
    F2 <- right_solve(A2); P <- tcrossprod(L)
    value[ix$L] <- P[lower.tri(P,diag=TRUE)]
    value[ix$A2] <- as.vector(F2)
    J[c(ix$L,ix$A2),] <- 0
    ij <- which(lower.tri(L,diag=TRUE),arr.ind=TRUE)
    for (j in seq_len(ix$nl)) {
      E <- matrix(0,r,r); E[ij[j,1],ij[j,2]] <- 1
      dP <- E%*%t(L)+L%*%t(E)
      J[ix$L,ix$L[j]] <- dP[lower.tri(P,diag=TRUE)]
      J[ix$A2,ix$L[j]] <- as.vector(-right_solve(F2%*%E))
    }
    if (length(ix$A2)) for (j in seq_along(ix$A2)) {
      E <- matrix(0,big_T-r,r); E[j] <- 1
      J[ix$A2,ix$A2[j]] <- as.vector(right_solve(E))
    }
    if (k > 0L) {
      R <- matrix(vec[ix$M],r,k*big_T); value[ix$M] <- as.vector(L%*%R)
      J[ix$M,] <- 0
      for (j in seq_len(ix$nl)) {
        E <- matrix(0,r,r); E[ij[j,1],ij[j,2]] <- 1
        J[ix$M,ix$L[j]] <- as.vector(E%*%R)
      }
      for (j in seq_along(ix$M)) {
        E <- matrix(0,r,k*big_T); E[j] <- 1
        J[ix$M,ix$M[j]] <- as.vector(L%*%E)
      }
    }
  }
  list(value = value, jacobian = J)
}

.stable_report_inference <- function(result, V, N, with_x = FALSE) {
  vec <- result$inference_parameter_hat
  k <- if (with_x) result$k else 0L
  pp <- if (with_x) result$parameter_hat else result$vartheta_hat
  mp <- tryCatch(.stable_reporting_map(vec,result$p,result$q,result$r,result$big_T,
    result$heteroskedastic,k),error=function(e)NULL)
  Vr <- matrix(NA_real_,length(pp),length(pp))
  rows_ok <- rep(FALSE,length(pp))
  if (!is.null(mp) && all(is.finite(V))) {
    rows_ok <- apply(mp$jacobian,1,function(v)all(is.finite(v)))
    R <- mp$jacobian[rows_ok,,drop=FALSE]
    Vr[rows_ok,rows_ok] <- R%*%V%*%t(R)
  }
  result <- if (with_x) .qmle_x_report_inference(result, Vr, N) else
    .set_reported_inference(result, Vr, N)
  dimnames(V) <- list(names(vec),names(vec))
  result$V_inference <- (V + t(V)) / 2
  result$reporting_jacobian <- if (is.null(mp)) NULL else mp$jacobian
  if (!is.null(result$reporting_jacobian))
    dimnames(result$reporting_jacobian) <- list(names(pp),names(vec))
  if (is.null(result$inference_diagnostics)) result$inference_diagnostics <- list()
  result$inference_diagnostics$reporting_coordinates <- if(with_x)"paper (vartheta, Xi, delta)"else"paper vartheta"
  result$inference_diagnostics$reporting_jacobian_finite <- !is.null(mp)&&all(is.finite(mp$jacobian))
  result$inference_diagnostics$paper_standard_errors_finite <- all(is.finite(result$se))
  result$inference_diagnostics$V_inference_finite <- all(is.finite(V))
  result$inference_diagnostics$V_hat_finite <- all(is.finite(result$V_hat))
  target <- grep("^(phi_|theta_|beta_)",names(result$se))
  result$inference_diagnostics$target_standard_errors_finite <- all(is.finite(result$se[target]))
  if(identical(result$inference_diagnostics$first_failure,"none") &&
     !result$inference_diagnostics$target_standard_errors_finite)
    result$inference_diagnostics$first_failure <- "reported_target_variance"
  result$inference_diagnostics$nuisance_reporting_failure <-
    result$inference_diagnostics$target_standard_errors_finite && !all(is.finite(result$se))
  result
}

.stable_qmle_components <- function(parameters,p,q,r,big_T,heteroskedastic,
                                    Y,W=NULL,Sigma_hat=NULL,hessian_type="obs") {
  N<-nrow(Y);d<-length(parameters);dv<-count_params(p,q,r,big_T,heteroskedastic)
  with_x<-!is.null(W);obs<-hessian_type%in%c("obs","both");info<-hessian_type%in%c("info","both")
  cv<-.stable_covariance(parameters,p,q,r,big_T,heteroskedastic,second=obs)
  mn<-if(with_x).stable_mean(parameters,p,q,r,big_T,ncol(W)/big_T,heteroskedastic,cv,
    second=obs)else NULL
  Sigma<-cv$Sigma;R<-chol(Sigma);K<-chol2inv(R)
  reconstruction<-if(is.null(Sigma_hat))0 else max(abs(Sigma-Sigma_hat))/max(1,max(abs(Sigma_hat)))
  if(with_x){
    aug<-cbind(1,W);res<-Y-aug%*%t(mn$map);qw<-qr(aug,LAPACK=TRUE)
    RW<-qr.R(qw)[,order(qw$pivot),drop=FALSE]/sqrt(N)
    means<-vapply(mn$first,function(M)as.vector(forwardsolve(t(R),M)%*%t(RW)),numeric(big_T*nrow(RW)))
    mean_info<-crossprod(means);resK<-res%*%K;QZ<-t(crossprod(aug,resK)/N)
    mixed<-lapply(cv$first,function(M)K%*%M%*%QZ)
  }else{res<-Y;resK<-res%*%K;mean_info<-matrix(0,d,d)}
  S<-crossprod(res)/N;KS<-K%*%S%*%K;A<-K-KS;scores<-matrix(0,N,d)
  invR<-backsolve(R,diag(big_T));cols<-matrix(0,big_T^2,d)
  for(a in seq_len(d)){
    if(with_x)scores[,a]<-rowSums((aug%*%t(mn$first[[a]]))*resK)
    if(a<=dv){Sa<-cv$first[[a]];scores[,a]<-scores[,a]+.5*(rowSums((resK%*%Sa)*resK)-sum(K*Sa))
      if(info)cols[,a]<-as.vector(crossprod(invR,Sa%*%invR))}
  }
  Hi<-if(info)mean_info+crossprod(cols)/2 else NULL;Ho<-if(obs)mean_info else NULL
  if(obs){
    for(b in seq_len(dv)){
      Sb<-cv$first[[b]];Ab<--K%*%Sb%*%K+K%*%Sb%*%KS+KS%*%Sb%*%K
      for(a in seq_len(b)){
        h<-.5*(sum(Ab*cv$first[[a]])+sum(A*cv$second(a,b)))
        Ho[a,b]<-Ho[a,b]+h;if(a!=b)Ho[b,a]<-Ho[b,a]+h
      }
    }
    if(with_x)for(b in seq_len(d))for(a in seq_len(b)){
      h<--sum(mn$second(a,b)*QZ)
      if(b<=dv)h<-h+sum(mn$first[[a]]*mixed[[b]])
      if(a<=dv)h<-h+sum(mn$first[[b]]*mixed[[a]])
      Ho[a,b]<-Ho[a,b]+h;if(a!=b)Ho[b,a]<-Ho[b,a]+h
    }
  }
  if(obs)Ho<-(Ho+t(Ho))/2
  if(any(!is.finite(scores))||(obs&&any(!is.finite(Ho)))||(info&&any(!is.finite(Hi))))stop("Nonfinite stable likelihood derivatives.")
  list(H_obs=Ho,H_info=Hi,J=crossprod(scores)/N,scores=scores,
    Sigma_jac=do.call(cbind,lapply(cv$first,as.vector)),
    mean_jac=if(with_x)do.call(cbind,lapply(mn$first,as.vector))else NULL,
    covariance_reconstruction_relative_error=reconstruction)
}

.stable_qmle_inference <- function(result,data,p,q,r,heteroskedastic,hessian_type,
                                   with_x=FALSE){
  result$p<-p;result$q<-q;result$r<-r;result$big_T<-data$big_T;result$N<-data$N
  result$heteroskedastic<-heteroskedastic;if(with_x)result$k<-data$k
  parameters<-.stable_inference_parameters(result,with_x);d<-length(parameters);N<-data$N
  obs<-hessian_type%in%c("obs","both");info<-hessian_type%in%c("info","both")
  primary<-if(hessian_type=="info")"info"else"obs";err<-NULL
  z<-tryCatch(.stable_qmle_components(parameters,p,q,r,data$big_T,heteroskedastic,
      if(with_x)data$Y else data$Y_c,if(with_x)data$W else NULL,result$Sigma_hat,hessian_type),
    error=function(e){err<<-conditionMessage(e);NULL})
  if(is.null(z))z<-list(H_obs=matrix(NA_real_,d,d),H_info=matrix(NA_real_,d,d),J=matrix(NA_real_,d,d),
    scores=matrix(NA_real_,N,d),Sigma_jac=NULL,mean_jac=NULL,covariance_reconstruction_relative_error=NA_real_)
  ob<-if(obs).qmle_analytic_bread(z$H_obs,d,require_pd=with_x)else NULL
  ib<-if(info).qmle_analytic_bread(z$H_info,d,require_pd=TRUE)else NULL
  bread<-if(primary=="info")ib else ob;V<-.qmle_analytic_sandwich(bread,z$scores,N)
  failure<-if(!is.null(err))"analytical_derivatives"else if(!bread$finite)"H_hat"else if(!bread$scaled_finite)"H_scaling"else
    if((with_x||primary=="info")&&!bread$positive_definite)"H_hat_not_positive_definite"else
      if(bread$inverse_method=="failed")"H_inverse"else if(any(!is.finite(z$J)))"J_hat"else if(any(!is.finite(V)))"V_hat"else"none"
  named<-function(M){dimnames(M)<-list(names(parameters),names(parameters));M}
  result$inference_parameter_hat<-parameters
  result$parameterization$inference <- "physical dynamics/Pi and stable triangular A factor coordinates"
  result$H_inference<-result$H_hat<-named(bread$H);result$J_inference<-result$J_hat<-named(z$J)
  if(obs)result$H_obs<-named(ob$H);if(info)result$H_info<-named(ib$H)
  result$inference_diagnostics<-list(first_failure=failure,hessian_type=hessian_type,primary_hessian=primary,
    analytical_error=err,
    inference_coordinates="physical dynamics/variances/Pi; triangular A factor chart; joint mean R when X",
    H_hat_finite=bread$finite,rank_H_hat=bread$rank,rcond_H_hat=bread$rcond,H_hat_positive_definite=bread$positive_definite,
    H_inverse_method=bread$inverse_method,H_eigenvalues=bread$eigenvalues,H_scaled_finite=bread$scaled_finite,
    H_scaled_rcond=bread$scaled_rcond,H_scaled_rank=bread$scaled_rank,H_diagonal_scale=bread$scale,
    H_inverse_backward_error=bread$inverse_backward_error,
    H_obs_finite=obs&&!is.null(ob)&&ob$finite,H_obs_positive_definite=obs&&!is.null(ob)&&ob$positive_definite,
    H_info_finite=info&&!is.null(ib)&&ib$finite,H_info_positive_definite=info&&!is.null(ib)&&ib$positive_definite,
    H_info_eigenvalues=if(info)ib$eigenvalues else NULL,
    Sigma_jac_finite=!is.null(z$Sigma_jac)&&all(is.finite(z$Sigma_jac)),
    rank_Sigma_jac=if(!is.null(z$Sigma_jac))qr(z$Sigma_jac,tol=1e-10)$rank else NA_integer_,
    mean_jac_finite=if(with_x)!is.null(z$mean_jac)&&all(is.finite(z$mean_jac))else NULL,
    rank_mean_jac=if(with_x&&!is.null(z$mean_jac))qr(z$mean_jac,tol=1e-10)$rank else NULL,
    Sigma_hat_cholesky_ok=is.null(err),J_hat_finite=all(is.finite(z$J)),V_hat_finite=all(is.finite(V)),
    score_column_means=setNames(colMeans(z$scores),names(parameters)),max_abs_mean_score=max(abs(colMeans(z$scores))),
    covariance_reconstruction_relative_error=z$covariance_reconstruction_relative_error)
  result<-.stable_report_inference(result,V,N,with_x)
  if(hessian_type=="both"){
    result$V_hat_obs<-result$V_hat;result$se_obs<-result$se;result$ci_obs<-result$ci
    vi<-.qmle_analytic_sandwich(ib,z$scores,N)
    tmp<-.stable_report_inference(result,vi,N,with_x)
    result$V_hat_info<-tmp$V_hat;result$se_info<-tmp$se;result$ci_info<-tmp$ci
    result$inference_diagnostics$V_hat_info_finite<-all(is.finite(vi))
    result$inference_diagnostics$information_standard_errors_finite<-all(is.finite(tmp$se))
  }
  result
}



# Stable inference module: stable_gmm.R
# GMM post-estimation inference in stable factor coordinates.
# Point objectives, optimizer calls, weighting matrices and moment definitions
# remain in their existing routines. No singular-value truncation is performed:
# all parameter directions must pass the existing full-rank threshold.

.stable_gmm_weighted_inverse <- function(G, W) {
  d <- ncol(G)
  R_W <- tryCatch(chol(W), error = function(e) NULL)
  weighted <- if (!is.null(R_W)) tryCatch(R_W %*% G, error = function(e) NULL) else NULL
  decomp <- if (!is.null(weighted) && nrow(weighted) >= d &&
                 all(is.finite(weighted))) tryCatch(
    svd(weighted, nu = d, nv = d), error = function(e) NULL) else NULL
  singular <- if (is.null(decomp)) rep(NA_real_, d) else decomp$d
  rank <- if (all(is.finite(singular))) {
    cutoff <- max(dim(weighted)) * .Machine$double.eps * max(singular)
    sum(singular > cutoff)
  } else NA_integer_
  complete <- !is.na(rank) && rank == d && !is.null(decomp)
  inverse <- matrix(NA_real_, d, d)
  influence <- matrix(NA_real_, d, nrow(G))
  if (complete) {
    inverse_columns <- sweep(decomp$v, 2L, singular, "/")
    inverse <- tcrossprod(inverse_columns)
    influence <- inverse_columns %*% t(decomp$u) %*% R_W
    complete <- all(is.finite(inverse)) && all(is.finite(influence))
  }
  H <- if (!is.null(R_W) && all(is.finite(G))) crossprod(G, W %*% G) else
    matrix(NA_real_, d, d)
  list(full_rank = complete, rank = rank, singular_values = singular,
       inverse = inverse, influence = influence, H = H,
       rcond_H = if (all(is.finite(H))) tryCatch(rcond(H), error = function(e) NA_real_) else NA_real_,
       W_positive_definite = !is.null(R_W))
}

.stable_gmm_no_x_inference <- function(result, data, p, q, r,
                                      heteroskedastic) {
  N <- data$N; big_T <- data$big_T
  parameters <- .stable_inference_parameters(result, with_x = FALSE)
  d <- length(parameters); m <- as.integer(big_T * (big_T + 1L) / 2L)
  pnames <- names(parameters)
  derivative_error <- NULL
  cv <- tryCatch(.stable_covariance(parameters, p, q, r, big_T,
    heteroskedastic, second = FALSE), error = function(e) {
      derivative_error <<- conditionMessage(e); NULL
    })
  G <- if (!is.null(cv)) vapply(cv$first, vech, numeric(m)) else
    matrix(NA_real_, m, d)
  dimnames(G) <- list(NULL, pnames)
  valid_G <- all(dim(G) == c(m, d)) && all(is.finite(G))
  rank_G <- if (valid_G) tryCatch(qr(G, tol = 1e-16)$rank,
                                 error = function(e) NA_integer_) else NA_integer_
  Upsilon <- result$Upsilon_hat
  if ((!is.matrix(Upsilon) || any(dim(Upsilon) != c(m,m)) ||
       any(!is.finite(Upsilon))) && is.matrix(data$Y)) {
    Upsilon <- .gmm_data_moments(data$Y, compute_Upsilon = TRUE)$Upsilon_hat
    result$Upsilon_hat <- Upsilon
  }
  valid_Upsilon <- is.matrix(Upsilon) && all(dim(Upsilon) == c(m, m)) &&
    all(is.finite(Upsilon))
  W <- result$W
  valid_W <- is.matrix(W) && all(dim(W) == c(m, m)) && all(is.finite(W))
  Upsilon_rcond <- if (valid_Upsilon) tryCatch(rcond(Upsilon), error = function(e) NA_real_) else NA_real_
  solve <- if (valid_G && valid_W) .stable_gmm_weighted_inverse(G, W) else NULL
  H <- if (!is.null(solve)) solve$H else matrix(NA_real_, d, d)
  J <- if (valid_G && valid_W && valid_Upsilon)
    crossprod(G, W %*% Upsilon %*% W %*% G) else matrix(NA_real_, d, d)
  V <- matrix(NA_real_, d, d)
  if (!is.null(solve) && solve$full_rank && valid_Upsilon) {
    V <- if (isTRUE(result$efficient_weight)) solve$inverse else
      solve$influence %*% Upsilon %*% t(solve$influence)
    V <- (V + t(V)) / 2
  }
  for(nm in c("H", "J", "V")) {
    z <- get(nm); dimnames(z) <- list(pnames, pnames); assign(nm, z)
  }
  regular <- !is.null(solve) && solve$full_rank && !is.na(rank_G) &&
    rank_G == d && all(is.finite(V))
  J_stat <- if (isTRUE(result$efficient_weight) && regular) N * result$obj_val else NA_real_
  df_overid <- if (regular) as.integer(m - d) else NA_integer_
  J_pval <- if (is.finite(J_stat) && !is.na(df_overid) && df_overid > 0L)
    pchisq(J_stat, df_overid, lower.tail = FALSE) else NA_real_
  result$inference_parameter_hat <- parameters
  result$G_inference <- G
  result$H_inference <- H
  result$J_inference <- J
  result$V_inference <- V
  result$Sigma_inference <- if (!is.null(cv)) cv$Sigma else NULL
  result$G <- G
  result$rank_G <- rank_G
  result$rank_weighted_G <- if (!is.null(solve)) solve$rank else NA_integer_
  result$weighted_G_singular_values <- if (!is.null(solve)) solve$singular_values else rep(NA_real_, d)
  result$rcond_GWG <- if (!is.null(solve)) solve$rcond_H else NA_real_
  first_failure <- if (!valid_G) "G" else if (!valid_Upsilon) "Upsilon_hat" else
    if (is.null(solve) || !solve$W_positive_definite) "W_positive_definite" else
      if (!solve$full_rank) "weighted_G_inverse" else if (!all(is.finite(V))) "V_hat" else "none"
  result$inference_diagnostics <- list(
    first_failure = first_failure,
    inference_coordinates = "physical dynamics/variances/Pi; stable factor (L_Psi,A2)",
    reported_coordinates = "paper vartheta", derivative_method = "analytical",
    analytical_error = derivative_error, G_finite = valid_G,
    Upsilon_hat_finite = valid_Upsilon, rcond_Upsilon_hat = Upsilon_rcond,
    rank_G = rank_G, rank_weighted_G = result$rank_weighted_G,
    rcond_GWG = result$rcond_GWG,
    covariance_reconstruction_relative_error = if (!is.null(cv) &&
      is.matrix(result$Sigma_hat) && all(is.finite(result$Sigma_hat)))
        max(abs(cv$Sigma-result$Sigma_hat))/(1+max(abs(result$Sigma_hat))) else NA_real_)
  # Reporting appends Jacobian/target/nuisance diagnostics; do not overwrite
  # them afterward or turn nuisance-only reporting trouble into target failure.
  result <- .stable_report_inference(result, V, N, with_x = FALSE)
  result$J_stat <- J_stat; result$df_overid <- df_overid; result$J_pval <- J_pval
  result$hansen_test <- NULL
  if (isTRUE(result$efficient_weight) && is.finite(J_stat) &&
      !is.na(df_overid) && df_overid > 0L && is.finite(J_pval)) {
    result$hansen_test <- structure(list(statistic = setNames(J_stat, "J"),
      parameter = setNames(df_overid, "df"), p.value = J_pval,
      method = "Hansen's J test of overidentifying restrictions",
      data.name = "efficient GMM moment conditions"), class = "htest")
  }
  result
}

# The Xi coefficient called A_mundlak by the paper parser contains R=L^-1 M
# in this internal inference vector; no point-estimation Xi is overwritten.
.stable_gmm_x_model <- function(outer, data, p, q, r, heteroskedastic) {
  T <- data$big_T; d_v <- count_params(p, q, r, T, heteroskedastic)
  n_Xi <- .qmle_x_count_Xi(p, r, data$k, T)
  if (length(outer) != d_v+n_Xi || any(!is.finite(outer)))
    stop("Incorrect or nonfinite stable GMM regressor parameter vector.")
  cv <- .stable_covariance(outer[seq_len(d_v)], p, q, r, T,
                           heteroskedastic, second = FALSE)
  if (is.null(tryCatch(chol(cv$Sigma), error = function(e) NULL)))
    stop("Stable GMM covariance is not positive definite.")
  Xi <- .qmle_x_unpack_Xi(outer[d_v+seq_len(n_Xi)], p, r, data$k, T)
  C <- matrix(0, T, data$k*T)
  if (p > 0L) C[seq_len(p),] <- Xi$c
  for (tt in seq.int(p+1L,T)) C[tt,(tt-1L)*data$k+seq_len(data$k)] <- Xi$beta
  if (r > 0L) C <- C + cv$A %*% Xi$A_mundlak
  list(Sigma = cv$Sigma, C = C, K = cv$Gamma %*% C,
       phi = if (p > 0L) outer[seq_len(p)] else numeric(0))
}

.stable_gmm_x_moments <- function(outer, data, p, q, r, heteroskedastic) {
  tryCatch({
    z <- .stable_gmm_x_model(outer, data, p, q, r, heteroskedastic)
    KS_yw <- z$K %*% t(data$S_yw)
    S_z <- data$S_yy-KS_yw-t(KS_yw)+z$K%*%data$S_ww%*%t(z$K)
    S_z <- (S_z+t(S_z))/2
    S_zw <- data$S_yw-z$K%*%data$S_ww
    c(vech(S_z/data$c_N-z$Sigma),as.vector(S_zw/data$c_N))
  }, error = function(e) rep(NA_real_,data$m_count))
}

.stable_gmm_x_delta <- function(outer, data, p, q, r, heteroskedastic) {
  tryCatch({
    z <- .stable_gmm_x_model(outer,data,p,q,r,heteroskedastic)
    drop(build_B_p(z$phi,data$big_T)%*%data$Y_bar-z$C%*%data$W_bar)
  },error=function(e)rep(NA_real_,data$big_T))
}

.stable_gmm_x_inference <- function(result, data, p, q, r,
                                   heteroskedastic) {
  N <- data$N; T <- data$big_T; m <- data$m_count
  full_hat <- .stable_inference_parameters(result,with_x=TRUE)
  outer <- head(full_hat,length(full_hat)-T)
  d <- length(outer);d_full <- length(full_hat)
  outer_names <- names(outer);full_names <- names(full_hat)
  moment_names <- data$moment_names
  derivative <- .gmm_x_numerical_jacobian(function(v)
    .stable_gmm_x_moments(v,data,p,q,r,heteroskedastic),outer,c(m,d))
  delta_derivative <- .gmm_x_numerical_jacobian(function(v)
    .stable_gmm_x_delta(v,data,p,q,r,heteroskedastic),outer,c(T,d))
  D <- derivative$value;D_delta <- delta_derivative$value
  dimnames(D) <- list(moment_names,outer_names)
  dimnames(D_delta) <- list(names(result$delta_hat),outer_names)
  n_Xi <- length(result$Xi_hat)
  mean_rank <- result$mean_design_rank
  if(is.null(mean_rank))mean_rank <- tryCatch({
    mn <- .stable_mean(full_hat,p,q,r,T,data$k,heteroskedastic,second=FALSE)
    dv <- count_params(p,q,r,T,heteroskedastic)
    design <- vapply(dv+seq_len(n_Xi),function(j)
      as.vector(data$W_c%*%t(mn$first[[j]][,-1L,drop=FALSE])),numeric(N*T))
    scales <- sqrt(colSums(design^2))
    if(any(scales<=0))0L else qr(sweep(design,2,scales,"/"),tol=1e-10)$rank
  },error=function(e)NA_integer_)
  mean_valid <- length(mean_rank)==1L&&is.finite(mean_rank)&&mean_rank==n_Xi
  covariance_valid <- is.null(result$pd_components)||isTRUE(result$pd_components$all_valid)
  products <- tryCatch(.gmm_x_individual(result$residuals,data),error=function(e)NULL)
  products_valid <- is.matrix(products)&&all(dim(products)==c(N,m))&&all(is.finite(products))
  Upsilon <- matrix(NA_real_,m,m);u <- NULL
  if(products_valid){
    u <- sweep(products,2L,colMeans(products),"-")
    Upsilon <- crossprod(u)/N;Upsilon <- (Upsilon+t(Upsilon))/2
  }
  dimnames(Upsilon)<-list(moment_names,moment_names)
  Upsilon_valid <- all(is.finite(Upsilon))
  rank_Upsilon <- if(Upsilon_valid)tryCatch(qr(Upsilon,tol=1e-10)$rank,error=function(e)NA_integer_)else NA_integer_
  rcond_Upsilon <- if(Upsilon_valid)tryCatch(rcond(Upsilon),error=function(e)NA_real_)else NA_real_
  Upsilon_chol <- if(Upsilon_valid)tryCatch(chol(Upsilon),error=function(e)NULL)else NULL
  W <- result$W;W_valid <- is.matrix(W)&&all(dim(W)==c(m,m))&&all(is.finite(W))
  solve <- if(derivative$success&&W_valid).stable_gmm_weighted_inverse(D,W)else NULL
  rank_D <- if(derivative$success)tryCatch(qr(D,tol=1e-10)$rank,error=function(e)NA_integer_)else NA_integer_
  H <- if(!is.null(solve))solve$H else matrix(NA_real_,d,d)
  J <- if(derivative$success&&W_valid&&Upsilon_valid)
    crossprod(D,W%*%Upsilon%*%W%*%D)else matrix(NA_real_,d,d)
  dimnames(H)<-dimnames(J)<-list(outer_names,outer_names)
  A <- if(!is.null(solve))solve$influence else matrix(NA_real_,d,m)
  dimnames(A)<-list(outer_names,moment_names)
  V <- matrix(NA_real_,d_full,d_full,dimnames=list(full_names,full_names))
  influence_means<-setNames(rep(NA_real_,d_full),full_names)
  can_infer<-derivative$success&&delta_derivative$success&&mean_valid&&covariance_valid&&
    products_valid&&Upsilon_valid&&!is.null(solve)&&solve$full_rank&&all(is.finite(A))
  if(can_infer){
    IF_outer <- -u%*%t(A)
    phi <- if(p>0L)outer[seq_len(p)]else numeric(0)
    IF_delta <- result$residuals%*%t(build_B_p(phi,T))+IF_outer%*%t(D_delta)
    IF_full <- cbind(IF_outer,IF_delta)
    influence_means <- setNames(colMeans(IF_full),full_names)
    IF_full <- sweep(IF_full,2L,influence_means,"-")
    V <- crossprod(IF_full)/N;dimnames(V)<-list(full_names,full_names)
  }
  first_failure<-if(!derivative$success)"moment_jacobian"else
    if(!delta_derivative$success)"delta_jacobian"else if(!mean_valid)"mean_design_rank"else
      if(!covariance_valid)"covariance_restrictions"else if(!products_valid||!Upsilon_valid)"Upsilon_hat"else
        if(is.null(solve)||!solve$W_positive_definite)"W_positive_definite"else
          if(!solve$full_rank)"weighted_G_inverse"else if(!all(is.finite(V)))"V_hat"else "none"
  pilot<-result$Upsilon_weight_hat
  pilot_chol<-if(is.matrix(pilot)&&all(dim(pilot)==c(m,m))&&all(is.finite(pilot)))
    tryCatch(chol(pilot),error=function(e)NULL)else NULL
  J_eligible<-isTRUE(result$efficient_weight)&&!is.null(pilot_chol)&&can_infer&&
    is.finite(result$obj_val)&&result$obj_val>=0
  J_stat<-if(J_eligible)N*result$obj_val else NA_real_
  df_overid<-if(can_infer)as.integer(m-d)else NA_integer_
  J_pval<-if(is.finite(J_stat)&&!is.na(df_overid)&&df_overid>0L)
    pchisq(J_stat,df_overid,lower.tail=FALSE)else NA_real_
  result$inference_parameter_hat<-full_hat
  result$G_inference<-result$G<-result$moment_jacobian<-D
  result$Upsilon_hat<-Upsilon
  result$A_hat<-A;result$D_delta_eta<-D_delta
  result$H_inference<-result$H_hat<-H
  result$J_inference<-result$J_hat<-J
  result$V_inference<-V
  result$rank_G<-rank_D
  result$rank_weighted_G<-if(!is.null(solve))solve$rank else NA_integer_
  result$weighted_G_singular_values<-if(!is.null(solve))solve$singular_values else rep(NA_real_,d)
  result$rcond_GWG<-if(!is.null(solve))solve$rcond_H else NA_real_
  result$J_stat<-J_stat;result$df_overid<-df_overid;result$J_pval<-J_pval;result$hansen_test<-NULL
  if (J_eligible && !is.na(df_overid) && df_overid > 0L && is.finite(J_pval)) {
    result$hansen_test <- structure(list(statistic = setNames(J_stat, "J"),
      parameter = setNames(df_overid, "df"), p.value = J_pval,
      method = "Hansen's J test of overidentifying restrictions",
      data.name = "efficient GMM moment conditions"), class = "htest")
  }
  result$inference_diagnostics<-list(first_failure=first_failure,
    inference_coordinates="physical dynamics/variances/Pi; stable factors and joint (beta,c,R,delta)",
    reported_coordinates="paper (vartheta,Xi,delta)",
    jacobian_definition="derivative of centered sample moments in stable factor coordinates",
    derivative_method="Richardson",G_finite=derivative$success,
    delta_jacobian_finite=delta_derivative$success,mean_design_valid=mean_valid,
    mean_design_rank=mean_rank,covariance_restrictions_valid=covariance_valid,
    rank_G=rank_D,rank_weighted_G=result$rank_weighted_G,rcond_GWG=result$rcond_GWG,
    W_positive_definite=!is.null(solve)&&solve$W_positive_definite,
    Upsilon_hat_finite=Upsilon_valid,rank_Upsilon_hat=rank_Upsilon,
    rcond_Upsilon_hat=rcond_Upsilon,Upsilon_hat_positive_definite=!is.null(Upsilon_chol),
    pilot_Upsilon_positive_definite=!is.null(pilot_chol),efficient_weight=isTRUE(result$efficient_weight),
    hansen_test_eligible=J_eligible,V_hat_finite=all(is.finite(V)),
    influence_column_means=influence_means,
    derivative_step_scales=c(moment_jacobian=derivative$step_scale,delta_jacobian=delta_derivative$step_scale))
  .stable_report_inference(result,V,N,with_x=TRUE)
}

.stable_gmm_inference <- function(result,data,p,q,r,heteroskedastic,
                                  with_x=FALSE) {
  result$p <- p;result$q <- q;result$r <- r
  result$big_T <- data$big_T;result$N <- data$N
  result$heteroskedastic <- heteroskedastic
  result$has_regressors <- with_x
  if(is.null(result$efficient_weight))
    result$efficient_weight <- identical(result$estimator,"Efficient GMM")
  if(with_x)result$k <- data$k
  ans <- if(with_x).stable_gmm_x_inference(result,data,p,q,r,heteroskedastic)else
    .stable_gmm_no_x_inference(result,data,p,q,r,heteroskedastic)
  ans$parameterization$inference <- "physical dynamics/Pi and stable triangular A factor coordinates"
  dd <- ans$inference_diagnostics
  dd$reporting_coordinates <- if(with_x)"paper (vartheta, Xi, delta)"else"paper vartheta"
  dd$reporting_jacobian_finite <- !is.null(ans$reporting_jacobian)&&all(is.finite(ans$reporting_jacobian))
  dd$paper_standard_errors_finite <- all(is.finite(ans$se))
  target <- grep("^(phi_|theta_|beta_)",names(ans$se))
  dd$target_standard_errors_finite <- all(is.finite(ans$se[target]))
  dd$nuisance_reporting_failure <- dd$target_standard_errors_finite&&!dd$paper_standard_errors_finite
  if(identical(dd$first_failure,"variance_diagonal")&&dd$nuisance_reporting_failure)dd$first_failure<-"none"
  if(identical(dd$first_failure,"none")&&!dd$target_standard_errors_finite)dd$first_failure<-"reported_target_variance"
  ans$inference_diagnostics <- dd
  ans
}

