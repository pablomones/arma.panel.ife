################################################################################
##  Model_selection_function.R
##
##  INFORMATION CRITERION FOR JOINT ORDER SELECTION  [Paper Section 4.2]
##
##  QMLE selection minimizes the fitted negative normalized quasi-log-likelihood
##  plus kappa * d(p, q, r) * log(N) / N [Paper eqs. (18)-(19)].
##  GMM selection uses the minimized GMM objective with the same penalty.
##
##  Default candidates satisfy 2(p + q + r + 1) - min(p, q) <= T.
##  With the default heteroskedasticity setting, count_params() gives
##    d(p, q, r) = p + q + (T - p + q) + r(r+1)/2 + (T-r)r + p(p+1)/2 + K(p,q).
##
##  Dependencies: Estimation_functions.R  (must be sourced before this file).
################################################################################

# ============================================================================ #
#  SECTION 1.  HELPERS
# ============================================================================ #

# Returns TRUE when (p, q, r) satisfies Assumption A6 for a panel with big_T
# time periods: 2(p + q + r + 1) - min(p, q) <= T.
.is_feasible <- function(p, q, r, big_T) {
  2L * (p + q + r + 1L) - min(p, q) <= big_T
}

# Enumerate every triple (p, q, r) with p in p_range[1]:p_range[2],
# q in q_range[1]:q_range[2], r in r_range[1]:r_range[2]
# that satisfies Assumption A6.
# Returns a list of integer vectors named c(p = ., q = ., r = .).
.enumerate_candidates <- function(p_range, q_range, r_range, big_T) {
  out <- list()
  for (p in p_range[1L]:p_range[2L]) {
    for (q in q_range[1L]:q_range[2L]) {
      for (r in r_range[1L]:r_range[2L]) {
        if (.is_feasible(p, q, r, big_T))
          out[[length(out) + 1L]] <- c(p = p, q = q, r = r)
      }
    }
  }
  out
}

# ============================================================================ #
#  SECTION 2.  MAIN MODEL SELECTION FUNCTION
# ============================================================================ #

#' Joint model order selection using the GMM criterion
#'
#' For each feasible candidate triple (p, q, r) the GMM estimator is run and
#' the BIC-type information criterion
#'
#'   IC_N(p, q, r) = Q_hat_{N(p,q,r)} + kappa * d(p, q, r) * log(N) / N
#'
#' is evaluated.  The triple minimising IC_N is returned as the selected order.
#' A custom candidates list can restrict the order comparison.
#'
#' @param Y               N x T data matrix (N rows = individuals, T cols = time).
#' @param p_range         Integer vector c(min, max) giving the AR order search
#'                        range.  Use c(k, k) to fix the AR order to k exactly.
#'                        NULL (default) searches over all feasible AR orders.
#' @param q_range         Integer vector c(min, max) for the MA order search range.
#'                        NULL (default) searches over all feasible MA orders.
#' @param r_range         Integer vector c(min, max) for the factor count search
#'                        range.  NULL (default) searches over all feasible counts.
#' @param candidates      Optional list of integer triples c(p = ., q = ., r = .)
#'                        specifying a custom candidate set.  When supplied,
#'                        p_range, q_range, and r_range are ignored entirely.
#'                        Infeasible candidates are warned about but still fitted.
#' @param kappa           Penalty parameter: a finite positive scalar or numeric
#'                        vector (all > 0).  Default: 1.  When a vector is
#'                        supplied each value yields its own selected order; the
#'                        return value is a named list keyed by "kappa=<value>"
#'                        rather than a single result list.  GMM is estimated
#'                        only once per candidate regardless of length(kappa).
#' @param heteroskedastic Logical.  TRUE (default) allows time-varying innovation
#'                        variances, matching the paper's baseline specification.
#' @param W               Positive definite dot_T x dot_T weighting matrix passed
#'                        to estimate_gmm(). NULL (default) uses the identity
#'                        matrix, matching estimate_gmm(). One W is shared by all candidates.
#' @param n_starts        Number of optimisation restarts per candidate. Default 15.
#' @param maxit           Maximum optimiser iterations.  Default 3000.
#' @param factr           L-BFGS-B tolerance factor for optim(). Default 1e6.
#' @param seed            Root optimizer seed. Each model receives its own seed
#'                        from model_optimizer_seed(seed,p,q,r). NULL draws one root seed.
#' @param ndeps           Numerical-gradient step passed to estimation; default 1e-5.
#' @param start_factr     MA-initialization tolerance passed to estimation; default 1e4.
#' @param verbose         If TRUE (default), print candidate and optimization-start progress.
#' @param stationary,invertible Logical; both default TRUE. TRUE constrains the
#'   corresponding reflection coefficients to [-0.9999,0.9999], using a
#'   fixed margin of 1e-4 from +/-1. FALSE uses physical coefficient bounds.
#' @param phi_bound,theta_bound Positive physical coefficient limits when the
#'   respective constraint is FALSE; both default 5. Ignored for constrained blocks.
#' @return When kappa is a scalar: a named list with elements
#'   p_hat, q_hat, r_hat, IC_min, Q_hat_min, d_selected, kappa, N, big_T, results
#'   (results is a data.frame with columns p, q, r, d, Q_hat, IC, converged).
#'   fit_diagnostics records one row per candidate; weight_diagnostics records
#'   how the shared weighting matrix was prepared. Failed, non-converged, or
#'   covariance-invalid fits have Q_hat = IC = Inf and are excluded; the search
#'   stops only if no eligible fit remains. Diagnostics retain the raw objective.
#'   When kappa is a vector: a named list with one entry per kappa value (names
#'   "kappa=<value>"), each entry having the same structure as above.
select_model_order_gmm <- function(
    Y,
    p_range         = NULL,
    q_range         = NULL,
    r_range         = NULL,
    candidates      = NULL,
    kappa           = 1,
    heteroskedastic = TRUE,
    W               = NULL,
    n_starts        = 15L,
    maxit           = 3000L,
    factr           = 1e6,
    seed            = NULL,
    verbose         = TRUE,
    stationary      = TRUE,
    invertible      = TRUE,
    phi_bound       = 5,
    theta_bound     = 5,
    ndeps           = 1e-5,
    start_factr     = 1e4
) {
  optimizer_seed_base <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
  invisible(model_optimizer_seed(optimizer_seed_base, 0L, 0L, 0L))
  optimizer_seed_base <- as.integer(optimizer_seed_base)

  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) < 1L || ncol(Y) < 1L ||
      any(!is.finite(Y)))
    stop("Y must be a nonempty finite numeric matrix without missing values.")

  if (!is.logical(stationary) || length(stationary) != 1L || is.na(stationary))
    stop("stationary must be TRUE or FALSE.")
  if (!is.logical(invertible) || length(invertible) != 1L || is.na(invertible))
    stop("invertible must be TRUE or FALSE.")
  phi_bound <- .validate_coefficient_bound(phi_bound, "phi_bound")
  theta_bound <- .validate_coefficient_bound(theta_bound, "theta_bound")

  if (!is.numeric(kappa) || !is.null(dim(kappa)) ||
      length(kappa) == 0L || any(!is.finite(kappa)) || any(kappa <= 0))
    stop("kappa must contain finite positive value(s).")
  kappa <- as.numeric(kappa)
  scalar_kappa <- (length(kappa) == 1L)

  N     <- nrow(Y)
  big_T <- ncol(Y)

  # ---- Build the candidate list --------------------------------------------
  if (!is.null(candidates)) {
    # candidates dominates p_range / q_range / r_range entirely
    if (!is.list(candidates) || length(candidates) == 0L)
      stop("candidates must be a non-empty list of integer triples c(p,q,r).")
    for (i in seq_along(candidates)) {
      cand <- candidates[[i]]
      if (!is.numeric(cand) || length(cand) != 3L ||
          any(!is.finite(cand)) || any(cand < 0) ||
          any(cand != floor(cand)) || any(cand > .Machine$integer.max))
        stop(sprintf("candidates[[%d]] must be a non-negative integer triple.", i))
      if (!is.null(names(cand))) {
        if (anyDuplicated(names(cand)) ||
            !setequal(names(cand), c("p", "q", "r")))
          stop(sprintf("candidates[[%d]] must have names p, q, r or no names.", i))
        cand <- cand[c("p", "q", "r")]
      }
      candidates[[i]] <- setNames(as.integer(cand), c("p", "q", "r"))
      cand <- candidates[[i]]
      if (!.is_feasible(cand[1L], cand[2L], cand[3L], big_T))
        warning(sprintf(
          "candidates[[%d]] = (%d,%d,%d) violates Assumption A6 for T=%d; there is no theoretical guarantee that fitting it will work.",
          i, cand[1L], cand[2L], cand[3L], big_T
        ), call. = FALSE)
    }
    cand_list <- lapply(candidates, function(cand) {
      x <- as.integer(cand)
      if (is.null(names(x))) names(x) <- c("p", "q", "r")
      x
    })
  } else {
    # Validate a supplied range vector; returns a checked integer c(min, max).
    .check_range <- function(rng, name) {
      if (!is.numeric(rng) || length(rng) != 2L ||
          any(!is.finite(rng)) || any(rng != floor(rng)) ||
          any(rng > .Machine$integer.max) ||
          rng[1L] < 0 || rng[2L] < rng[1L])
        stop(name, " must be a length-2 integer vector c(min, max) with 0 <= min <= max.")
      as.integer(rng)
    }

    # NULL => search over all feasible orders: upper bound floor((T-2)/2).
    max_order <- max(0L, as.integer(floor((big_T - 2L) / 2L)))

    p_range <- if (is.null(p_range)) c(0L, max_order) else .check_range(p_range, "p_range")
    q_range <- if (is.null(q_range)) c(0L, max_order) else .check_range(q_range, "q_range")
    r_range <- if (is.null(r_range)) c(0L, max_order) else .check_range(r_range, "r_range")

    cand_list <- .enumerate_candidates(p_range, q_range, r_range, big_T)
  }

  n_cand <- length(cand_list)
  if (n_cand == 0L)
    stop("No feasible candidate orders found.  T may be too small relative ",
         "to the supplied p_range / q_range / r_range bounds.")

  # Warn about large candidate sets before starting the (potentially long) search.
  if (n_cand >= 1000L) {
    warning(
      "The number of candidate models is at least 1000, selecting the model ",
      "order can be extremely time consuming.",
      call. = FALSE, immediate. = TRUE
    )
  } else if (n_cand > 100L) {
    warning(
      "The number of candidate models exceeds 100, selecting the model order ",
      "can be very time consuming.",
      call. = FALSE, immediate. = TRUE
    )
  }

  if (verbose)
    cat(sprintf(
      "Model selection: %d candidate(s), N = %d, T = %d, kappa = %s.\n",
      n_cand, N, big_T,
      paste(sprintf("%.4g", kappa), collapse = ", ")
    ))

  # ---- Prepare one weighting matrix and one set of sample moments ----------
  dot_T <- as.integer(big_T * (big_T + 1L) / 2L)
  identity_default_W <- is.null(W)
  if (identity_default_W) W <- diag(dot_T)
  W <- .validate_W(W, dot_T)
  gmm_data <- .gmm_data_moments(Y, compute_Upsilon = FALSE)
  weight_diagnostics <- list(
    weight_mode = if (identity_default_W) "identity" else "supplied",
    moment_count = dot_T,
    Upsilon_rcond = NA_real_,
    W_rcond = as.numeric(rcond(W)),
    Upsilon_chol_min_diag = NA_real_
  )

  # ---- Pre-allocate result vectors ----------------------------------------
  # Failed fits receive an infinite criterion and cannot be selected.
  Q_hat_vec     <- rep(Inf, n_cand)
  raw_objective_vec <- rep(NA_real_, n_cand)
  eligible_vec <- rep(FALSE, n_cand)
  d_vec         <- rep(NA_integer_, n_cand)
  converged_vec <- rep(FALSE,      n_cand)
  convergence_code_vec <- rep(NA_integer_, n_cand)
  pd_valid_vec <- rep(NA, n_cand)
  successful_starts_vec <- rep(NA_integer_, n_cand)
  failed_starts_vec <- rep(NA_integer_, n_cand)
  elapsed_sec_vec <- rep(NA_real_, n_cand)
  warnings_vec <- rep("", n_cand)
  bound_hit_any_vec <- rep(NA, n_cand)
  bound_hit_lower_vec <- rep(NA_character_, n_cand)
  bound_hit_upper_vec <- rep(NA_character_, n_cand)
  optimizer_seed_vec <- rep(NA_integer_, n_cand)

  # ---- Loop over all candidate models (GMM estimated once per candidate) --
  for (idx in seq_len(n_cand)) {
    p <- as.integer(cand_list[[idx]]["p"])
    q <- as.integer(cand_list[[idx]]["q"])
    r <- as.integer(cand_list[[idx]]["r"])
    optimizer_seed_vec[idx] <- model_optimizer_seed(optimizer_seed_base, p, q, r)

    d_pqr      <- count_params(p, q, r, big_T, heteroskedastic)
    d_vec[idx] <- d_pqr

    if (verbose)
      cat(sprintf("  (%d/%d) (p=%d, q=%d, r=%d)  d=%d  ...",
                  idx, n_cand, p, q, r, d_pqr))

    fit_warnings <- character()
    start_time <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      withCallingHandlers(
        estimate_gmm(
          Y               = Y,
          p               = p,
          q               = q,
          r               = r,
          W               = W,
          heteroskedastic = heteroskedastic,
          n_starts        = n_starts,
          maxit           = maxit,
          factr           = factr,
          ndeps           = ndeps,
          start_factr     = start_factr,
          seed            = optimizer_seed_vec[idx],
          compute_se      = FALSE,
          warn_dynamic    = FALSE,
          verbose         = verbose,
          efficient_weight = FALSE,
          .warn_T         = FALSE,  # A6 has already been applied
          .gmm_data       = gmm_data,
          stationary      = stationary,
          invertible      = invertible,
          phi_bound       = phi_bound,
          theta_bound     = theta_bound
        ),
        warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    elapsed_sec_vec[idx] <- as.numeric(proc.time()[["elapsed"]] - start_time)
    if (inherits(fit, "error")) {
      fit_warnings <- c(fit_warnings, paste0("ERROR: ", conditionMessage(fit)))
      if (verbose) cat(sprintf("  FAILED (%s)\n", conditionMessage(fit)))
    }
    warnings_vec[idx] <- paste(unique(fit_warnings), collapse = " | ")
    has_fit <- is.list(fit) && !inherits(fit, "error")
    if (has_fit) {
      if (length(fit$obj_val) == 1L && is.numeric(fit$obj_val))
        raw_objective_vec[idx] <- fit$obj_val
      converged_vec[idx] <- isTRUE(fit$converged)
      if (length(fit$convergence_code) == 1L)
        convergence_code_vec[idx] <- as.integer(fit$convergence_code)
      pd_valid_vec[idx] <- isTRUE(fit$pd_components$all_valid)
      if (length(fit$successful_starts) == 1L)
        successful_starts_vec[idx] <- as.integer(fit$successful_starts)
      if (length(fit$failed_starts) == 1L)
        failed_starts_vec[idx] <- as.integer(fit$failed_starts)
      if (!is.null(fit$bound_hits)) {
        bound_hit_any_vec[idx] <- isTRUE(fit$bound_hits$any)
        bound_hit_lower_vec[idx] <- paste(fit$bound_hits$lower, collapse = ";")
        bound_hit_upper_vec[idx] <- paste(fit$bound_hits$upper, collapse = ";")
      }
    }
    valid_objective <- has_fit && length(fit$obj_val) == 1L &&
      is.numeric(fit$obj_val) && is.finite(fit$obj_val) && fit$obj_val < 1e15
    valid_fit <- valid_objective && isTRUE(fit$converged) &&
      isTRUE(fit$pd_components$all_valid)
    eligible_vec[idx] <- valid_fit

    if (valid_fit) {
      Q_hat_vec[idx] <- fit$obj_val
      if (verbose) {
        if (scalar_kappa) {
          ic_val <- Q_hat_vec[idx] + kappa * d_pqr * log(N) / N
          cat(sprintf("  Q_hat = %10.6f   IC = %10.6f\n", Q_hat_vec[idx], ic_val))
        } else {
          cat(sprintf("  Q_hat = %10.6f\n", Q_hat_vec[idx]))
        }
      }
    } else if (has_fit) {
      reason <- if (!valid_objective) "non-finite or penalized objective" else
        "non-converged fit or invalid covariance components"
      warnings_vec[idx] <- paste(c(fit_warnings, paste0("EXCLUDED: ", reason)),
                                collapse = " | ")
      if (verbose) cat(sprintf("  FAILED (%s)\n", reason))
    }
  }

  # ---- Base results table (kappa-independent columns) ---------------------
  base_df <- data.frame(
    p         = vapply(cand_list, function(x) as.integer(x["p"]), integer(1L)),
    q         = vapply(cand_list, function(x) as.integer(x["q"]), integer(1L)),
    r         = vapply(cand_list, function(x) as.integer(x["r"]), integer(1L)),
    d         = d_vec,
    Q_hat     = Q_hat_vec,
    converged = converged_vec,
    stringsAsFactors = FALSE
  )
  fit_diagnostics <- data.frame(
    p = base_df$p, q = base_df$q, r = base_df$r,
    optimizer_seed = optimizer_seed_vec,
    d = base_df$d, Q_hat = base_df$Q_hat,
    converged = base_df$converged,
    eligible = eligible_vec,
    raw_objective = raw_objective_vec,
    convergence_code = convergence_code_vec,
    pd_valid = pd_valid_vec,
    successful_starts = successful_starts_vec,
    failed_starts = failed_starts_vec,
    elapsed_sec = elapsed_sec_vec,
    warnings = warnings_vec,
    bound_hit_any = bound_hit_any_vec,
    bound_hit_lower = bound_hit_lower_vec,
    bound_hit_upper = bound_hit_upper_vec,
    stringsAsFactors = FALSE
  )

  # ---- Build the result for a single kappa value --------------------------
  .build_one <- function(k) {
    ic_vec <- Q_hat_vec + k * d_vec * log(N) / N
    if (!any(is.finite(ic_vec)))
      stop("All candidate GMM estimations failed or returned non-finite IC values.")
    best_idx <- which.min(ic_vec)
    p_hat    <- as.integer(cand_list[[best_idx]]["p"])
    q_hat    <- as.integer(cand_list[[best_idx]]["q"])
    r_hat    <- as.integer(cand_list[[best_idx]]["r"])
    if (verbose) {
      if (scalar_kappa) {
        cat(sprintf(
          "\nSelected: (p=%d, q=%d, r=%d)  IC = %.6f  Q_hat = %.6f  d = %d\n",
          p_hat, q_hat, r_hat,
          ic_vec[best_idx], Q_hat_vec[best_idx], d_vec[best_idx]
        ))
      } else {
        cat(sprintf(
          "\nkappa = %.4g  ->  Selected: (p=%d, q=%d, r=%d)  IC = %.6f  Q_hat = %.6f  d = %d\n",
          k, p_hat, q_hat, r_hat,
          ic_vec[best_idx], Q_hat_vec[best_idx], d_vec[best_idx]
        ))
      }
    }
    res_df    <- base_df
    res_df$IC <- ic_vec
    res_df    <- res_df[, c("p", "q", "r", "d", "Q_hat", "IC", "converged")]
    list(
      p_hat      = p_hat,
      q_hat      = q_hat,
      r_hat      = r_hat,
      IC_min     = ic_vec[best_idx],
      Q_hat_min  = Q_hat_vec[best_idx],
      d_selected = d_vec[best_idx],
      kappa      = k,
      stationary = stationary,
      invertible = invertible,
      phi_bound  = phi_bound,
      theta_bound = theta_bound,
      optimizer_seed_base = optimizer_seed_base,
      N          = N,
      big_T      = big_T,
      results    = res_df,
      fit_diagnostics = fit_diagnostics,
      weight_diagnostics = weight_diagnostics
    )
  }

  # ---- Return: scalar kappa -> same format as before; vector -> named list -
  if (scalar_kappa) {
    .build_one(kappa)
  } else {
    out        <- lapply(kappa, .build_one)
    names(out) <- paste0("kappa=", kappa)
    out
  }
}

# ============================================================================ #
#  SECTION 3.  QMLE MODEL SELECTION
# ============================================================================ #

#' Joint model order selection using QMLE [Paper eqs. (18)-(19)]
#'
#' Fits estimate_qmle() once per candidate, then minimizes, for each kappa,
#'   IC_N = neg_loglik + kappa * d(p,q,r) * log(N) / N.
#' estimate_qmle()$obj_val is the maximized normalized quasi-log-likelihood:
#'   loglik = -log|Sigma_hat| - tr(S_N Sigma_hat^{-1}),
#' where S_N uses divisor N. Hence neg_loglik = -fit$obj_val. This equals
#' -2/N times the conventional concentrated Gaussian log likelihood up to a
#' constant common to all candidates; no further N or 1/2 scaling is applied.
#'
#' @param Y N x T finite numeric matrix, with individuals in rows.
#' @param p_range,q_range,r_range Inclusive integer c(min,max) order ranges.
#'   NULL searches all feasible orders. candidates overrides all three ranges.
#' @param candidates Optional list of nonnegative integer (p,q,r) triples.
#'   Explicit infeasible triples are warned about but still attempted.
#' @param kappa Finite positive scalar or vector; default 1. Each candidate is
#'   fitted once irrespective of length(kappa).
#' @param heteroskedastic Allow time-varying innovation variances; default TRUE.
#' @param n_starts,maxit,factr Optimization controls; defaults 15, 3000, 1e6.
#' @param seed Root optimizer seed; model_optimizer_seed gives each candidate
#'   an order-independent seed. NULL draws one root seed.
#' @param ndeps,start_factr Estimation controls; defaults 1e-5 and 1e4.
#' @param verbose Print candidate and optimization-start progress and each selected order; default TRUE.
#' @param stationary,invertible Logical; both default TRUE. TRUE constrains the
#'   corresponding reflection coefficients to [-0.9999,0.9999], using a
#'   fixed margin of 1e-4 from +/-1. FALSE uses physical coefficient bounds.
#' @param phi_bound,theta_bound Positive physical coefficient limits when the
#'   respective constraint is FALSE; both default 5. Ignored for constrained blocks.
#' @return For scalar kappa, a list containing p_hat, q_hat, r_hat, IC_min,
#'   neg_loglik_min, loglik_hat, d_selected, kappa, N, big_T, results, and
#'   fit_diagnostics. results has p, q, r, d, loglik, neg_loglik, IC, converged.
#'   For vector kappa, a named list ("kappa=<value>") of such result lists.
#'   Failed/non-converged/covariance-invalid fits have neg_loglik = IC = Inf.
#'   They remain in the diagnostics but cannot be selected. If all fits fail,
#'   the function stops. No standard errors or GMM weight are computed.
select_model_order_qmle <- function(
    Y,
    p_range         = NULL,
    q_range         = NULL,
    r_range         = NULL,
    candidates      = NULL,
    kappa           = 1,
    heteroskedastic = TRUE,
    n_starts        = 15L,
    maxit           = 3000L,
    factr           = 1e6,
    seed            = NULL,
    verbose         = TRUE,
    stationary      = TRUE,
    invertible      = TRUE,
    phi_bound       = 5,
    theta_bound     = 5,
    ndeps           = 1e-5,
    start_factr     = 1e4
) {
  optimizer_seed_base <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
  invisible(model_optimizer_seed(optimizer_seed_base, 0L, 0L, 0L))
  optimizer_seed_base <- as.integer(optimizer_seed_base)

  if (!is.matrix(Y) || !is.numeric(Y) || nrow(Y) < 1L || ncol(Y) < 1L ||
      any(!is.finite(Y)))
    stop("Y must be a nonempty finite numeric matrix without missing values.")

  if (!is.logical(stationary) || length(stationary) != 1L || is.na(stationary))
    stop("stationary must be TRUE or FALSE.")
  if (!is.logical(invertible) || length(invertible) != 1L || is.na(invertible))
    stop("invertible must be TRUE or FALSE.")
  phi_bound <- .validate_coefficient_bound(phi_bound, "phi_bound")
  theta_bound <- .validate_coefficient_bound(theta_bound, "theta_bound")

  if (!is.numeric(kappa) || !is.null(dim(kappa)) ||
      length(kappa) == 0L || any(!is.finite(kappa)) || any(kappa <= 0))
    stop("kappa must contain finite positive value(s).")
  kappa <- as.numeric(kappa)
  scalar_kappa <- (length(kappa) == 1L)

  N     <- nrow(Y)
  big_T <- ncol(Y)

  # ---- Build the candidate list --------------------------------------------
  if (!is.null(candidates)) {
    # candidates dominates p_range / q_range / r_range entirely
    if (!is.list(candidates) || length(candidates) == 0L)
      stop("candidates must be a non-empty list of integer triples c(p,q,r).")
    for (i in seq_along(candidates)) {
      cand <- candidates[[i]]
      if (!is.numeric(cand) || length(cand) != 3L ||
          any(!is.finite(cand)) || any(cand < 0) ||
          any(cand != floor(cand)) || any(cand > .Machine$integer.max))
        stop(sprintf("candidates[[%d]] must be a non-negative integer triple.", i))
      if (!is.null(names(cand))) {
        if (anyDuplicated(names(cand)) ||
            !setequal(names(cand), c("p", "q", "r")))
          stop(sprintf("candidates[[%d]] must have names p, q, r or no names.", i))
        cand <- cand[c("p", "q", "r")]
      }
      candidates[[i]] <- setNames(as.integer(cand), c("p", "q", "r"))
      cand <- candidates[[i]]
      if (!.is_feasible(cand[1L], cand[2L], cand[3L], big_T))
        warning(sprintf(
          "candidates[[%d]] = (%d,%d,%d) violates Assumption A6 for T=%d; there is no theoretical guarantee that fitting it will work.",
          i, cand[1L], cand[2L], cand[3L], big_T
        ), call. = FALSE)
    }
    cand_list <- lapply(candidates, function(cand) {
      x <- as.integer(cand)
      if (is.null(names(x))) names(x) <- c("p", "q", "r")
      x
    })
  } else {
    # Validate a supplied range vector; returns a checked integer c(min, max).
    .check_range <- function(rng, name) {
      if (!is.numeric(rng) || length(rng) != 2L ||
          any(!is.finite(rng)) || any(rng != floor(rng)) ||
          any(rng > .Machine$integer.max) ||
          rng[1L] < 0 || rng[2L] < rng[1L])
        stop(name, " must be a length-2 integer vector c(min, max) with 0 <= min <= max.")
      as.integer(rng)
    }

    # NULL => search over all feasible orders: upper bound floor((T-2)/2).
    max_order <- max(0L, as.integer(floor((big_T - 2L) / 2L)))

    p_range <- if (is.null(p_range)) c(0L, max_order) else .check_range(p_range, "p_range")
    q_range <- if (is.null(q_range)) c(0L, max_order) else .check_range(q_range, "q_range")
    r_range <- if (is.null(r_range)) c(0L, max_order) else .check_range(r_range, "r_range")

    cand_list <- .enumerate_candidates(p_range, q_range, r_range, big_T)
  }

  n_cand <- length(cand_list)
  if (n_cand == 0L)
    stop("No feasible candidate orders found.  T may be too small relative ",
         "to the supplied p_range / q_range / r_range bounds.")

  # Warn about large candidate sets before starting the (potentially long) search.
  if (n_cand >= 1000L) {
    warning(
      "The number of candidate models is at least 1000, selecting the model ",
      "order can be extremely time consuming.",
      call. = FALSE, immediate. = TRUE
    )
  } else if (n_cand > 100L) {
    warning(
      "The number of candidate models exceeds 100, selecting the model order ",
      "can be very time consuming.",
      call. = FALSE, immediate. = TRUE
    )
  }

  if (verbose)
    cat(sprintf(
      "Model selection: %d candidate(s), N = %d, T = %d, kappa = %s.\n",
      n_cand, N, big_T,
      paste(sprintf("%.4g", kappa), collapse = ", ")
    ))

  # ---- Pre-allocate result vectors ----------------------------------------
  # Failed fits receive an infinite criterion and cannot be selected.
  neg_loglik_vec     <- rep(Inf, n_cand)
  raw_objective_vec <- rep(NA_real_, n_cand)
  eligible_vec <- rep(FALSE, n_cand)
  d_vec         <- rep(NA_integer_, n_cand)
  converged_vec <- rep(FALSE,      n_cand)
  convergence_code_vec <- rep(NA_integer_, n_cand)
  pd_valid_vec <- rep(NA, n_cand)
  successful_starts_vec <- rep(NA_integer_, n_cand)
  failed_starts_vec <- rep(NA_integer_, n_cand)
  elapsed_sec_vec <- rep(NA_real_, n_cand)
  warnings_vec <- rep("", n_cand)
  bound_hit_any_vec <- rep(NA, n_cand)
  bound_hit_lower_vec <- rep(NA_character_, n_cand)
  bound_hit_upper_vec <- rep(NA_character_, n_cand)
  optimizer_seed_vec <- rep(NA_integer_, n_cand)

  # ---- Loop over all candidate models (QMLE estimated once per candidate) --
  for (idx in seq_len(n_cand)) {
    p <- as.integer(cand_list[[idx]]["p"])
    q <- as.integer(cand_list[[idx]]["q"])
    r <- as.integer(cand_list[[idx]]["r"])
    optimizer_seed_vec[idx] <- model_optimizer_seed(optimizer_seed_base, p, q, r)

    d_pqr      <- count_params(p, q, r, big_T, heteroskedastic)
    d_vec[idx] <- d_pqr

    if (verbose)
      cat(sprintf("  (%d/%d) (p=%d, q=%d, r=%d)  d=%d  ...",
                  idx, n_cand, p, q, r, d_pqr))

    fit_warnings <- character()
    start_time <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      withCallingHandlers(
        estimate_qmle(
          Y               = Y,
          p               = p,
          q               = q,
          r               = r,
          heteroskedastic = heteroskedastic,
          n_starts        = n_starts,
          maxit           = maxit,
          factr           = factr,
          ndeps           = ndeps,
          start_factr     = start_factr,
          seed            = optimizer_seed_vec[idx],
          compute_se      = FALSE,
          warn_dynamic    = FALSE,
          verbose         = verbose,
          stationary      = stationary,
          invertible      = invertible,
          phi_bound       = phi_bound,
          theta_bound     = theta_bound
        ),
        warning = function(w) {
          fit_warnings <<- c(fit_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    elapsed_sec_vec[idx] <- as.numeric(proc.time()[["elapsed"]] - start_time)
    if (inherits(fit, "error")) {
      fit_warnings <- c(fit_warnings, paste0("ERROR: ", conditionMessage(fit)))
      if (verbose) cat(sprintf("  FAILED (%s)\n", conditionMessage(fit)))
    }
    warnings_vec[idx] <- paste(unique(fit_warnings), collapse = " | ")
    has_fit <- is.list(fit) && !inherits(fit, "error")
    if (has_fit) {
      if (length(fit$obj_val) == 1L && is.numeric(fit$obj_val))
        raw_objective_vec[idx] <- fit$obj_val
      converged_vec[idx] <- isTRUE(fit$converged)
      if (length(fit$convergence_code) == 1L)
        convergence_code_vec[idx] <- as.integer(fit$convergence_code)
      pd_valid_vec[idx] <- isTRUE(fit$pd_components$all_valid)
      if (length(fit$successful_starts) == 1L)
        successful_starts_vec[idx] <- as.integer(fit$successful_starts)
      if (length(fit$failed_starts) == 1L)
        failed_starts_vec[idx] <- as.integer(fit$failed_starts)
      if (!is.null(fit$bound_hits)) {
        bound_hit_any_vec[idx] <- isTRUE(fit$bound_hits$any)
        bound_hit_lower_vec[idx] <- paste(fit$bound_hits$lower, collapse = ";")
        bound_hit_upper_vec[idx] <- paste(fit$bound_hits$upper, collapse = ";")
      }
    }
    valid_objective <- has_fit && length(fit$obj_val) == 1L &&
      is.numeric(fit$obj_val) && is.finite(fit$obj_val) && (-fit$obj_val) < 1e15
    valid_fit <- valid_objective && isTRUE(fit$converged) &&
      isTRUE(fit$pd_components$all_valid)
    eligible_vec[idx] <- valid_fit

    if (valid_fit) {
      neg_loglik_vec[idx] <- -fit$obj_val
      if (verbose) {
        if (scalar_kappa) {
          ic_val <- neg_loglik_vec[idx] + kappa * d_pqr * log(N) / N
          cat(sprintf("  neg_loglik = %10.6f   IC = %10.6f\n", neg_loglik_vec[idx], ic_val))
        } else {
          cat(sprintf("  neg_loglik = %10.6f\n", neg_loglik_vec[idx]))
        }
      }
    } else if (has_fit) {
      reason <- if (!valid_objective) "non-finite or penalized objective" else
        "non-converged fit or invalid covariance components"
      warnings_vec[idx] <- paste(c(fit_warnings, paste0("EXCLUDED: ", reason)),
                                collapse = " | ")
      if (verbose) cat(sprintf("  FAILED (%s)\n", reason))
    }
  }

  # ---- Base results table (kappa-independent columns) ---------------------
  base_df <- data.frame(
    p         = vapply(cand_list, function(x) as.integer(x["p"]), integer(1L)),
    q         = vapply(cand_list, function(x) as.integer(x["q"]), integer(1L)),
    r         = vapply(cand_list, function(x) as.integer(x["r"]), integer(1L)),
    d         = d_vec,
    loglik    = -neg_loglik_vec,
    neg_loglik = neg_loglik_vec,
    converged = converged_vec,
    stringsAsFactors = FALSE
  )
  fit_diagnostics <- data.frame(
    p = base_df$p, q = base_df$q, r = base_df$r,
    optimizer_seed = optimizer_seed_vec,
    d = base_df$d, neg_loglik = base_df$neg_loglik,
    converged = base_df$converged,
    eligible = eligible_vec,
    raw_objective = raw_objective_vec,
    convergence_code = convergence_code_vec,
    pd_valid = pd_valid_vec,
    successful_starts = successful_starts_vec,
    failed_starts = failed_starts_vec,
    elapsed_sec = elapsed_sec_vec,
    warnings = warnings_vec,
    bound_hit_any = bound_hit_any_vec,
    bound_hit_lower = bound_hit_lower_vec,
    bound_hit_upper = bound_hit_upper_vec,
    stringsAsFactors = FALSE
  )

  # ---- Build the result for a single kappa value --------------------------
  .build_one <- function(k) {
    ic_vec <- neg_loglik_vec + k * d_vec * log(N) / N
    if (!any(is.finite(ic_vec)))
      stop("All candidate QMLE estimations failed or returned non-finite IC values.")
    best_idx <- which.min(ic_vec)
    p_hat    <- as.integer(cand_list[[best_idx]]["p"])
    q_hat    <- as.integer(cand_list[[best_idx]]["q"])
    r_hat    <- as.integer(cand_list[[best_idx]]["r"])
    if (verbose) {
      if (scalar_kappa) {
        cat(sprintf(
          "\nSelected: (p=%d, q=%d, r=%d)  IC = %.6f  neg_loglik = %.6f  d = %d\n",
          p_hat, q_hat, r_hat,
          ic_vec[best_idx], neg_loglik_vec[best_idx], d_vec[best_idx]
        ))
      } else {
        cat(sprintf(
          "\nkappa = %.4g  ->  Selected: (p=%d, q=%d, r=%d)  IC = %.6f  neg_loglik = %.6f  d = %d\n",
          k, p_hat, q_hat, r_hat,
          ic_vec[best_idx], neg_loglik_vec[best_idx], d_vec[best_idx]
        ))
      }
    }
    res_df    <- base_df
    res_df$IC <- ic_vec
    res_df    <- res_df[, c("p", "q", "r", "d", "loglik", "neg_loglik", "IC", "converged")]
    list(
      p_hat      = p_hat,
      q_hat      = q_hat,
      r_hat      = r_hat,
      IC_min     = ic_vec[best_idx],
      neg_loglik_min = neg_loglik_vec[best_idx],
      loglik_hat = -neg_loglik_vec[best_idx],
      d_selected = d_vec[best_idx],
      kappa      = k,
      stationary = stationary,
      invertible = invertible,
      phi_bound  = phi_bound,
      theta_bound = theta_bound,
      optimizer_seed_base = optimizer_seed_base,
      N          = N,
      big_T      = big_T,
      results    = res_df,
      fit_diagnostics = fit_diagnostics
    )
  }

  # ---- Return: scalar kappa -> same format as before; vector -> named list -
  if (scalar_kappa) {
    .build_one(kappa)
  } else {
    out        <- lapply(kappa, .build_one)
    names(out) <- paste0("kappa=", kappa)
    out
  }
}
