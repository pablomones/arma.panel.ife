\name{select_model_order_qmle}
\alias{select_model_order_qmle}
\title{Joint Model Order Selection by a QMLE Information Criterion}
\description{Jointly select AR order, MA order, and factor count using a penalized fitted loss.}
\usage{
select_model_order_qmle(Y, p_range = NULL, q_range = NULL, r_range = NULL,
  candidates = NULL, kappa = 1, heteroskedastic = TRUE, n_starts = 20L, maxit
  = 3000L, factr = 1e4, seed = NULL, verbose = FALSE, stationary = TRUE,
  invertible = TRUE, phi_bound = 5, theta_bound = 5, ndeps = 1e-5, start_factr
  = 1e4)
}
\arguments{
\item{Y}{Required finite numeric N by T matrix: individuals in rows, dates in columns. No missing values.}
\item{p_range}{Inclusive integer pair c(min, max) for AR orders; NULL searches feasible orders.}
\item{q_range}{Inclusive integer pair for MA orders; NULL searches feasible orders.}
\item{r_range}{Inclusive integer pair for factor counts; NULL searches feasible orders.}
\item{candidates}{Optional nonempty list of nonnegative integer triples \code{c(p, q, r)}. Named triples must have exactly the names \code{p}, \code{q}, and \code{r}. Overrides the three ranges. Explicit infeasible triples generate warnings but are still attempted.}
\item{kappa}{Positive penalty scalar or vector. Each candidate is fitted once for all penalties.}
\item{heteroskedastic}{Logical. If true, estimate \eqn{T-p+q} innovation variances; otherwise estimate one common innovation variance. Default is true.}
\item{n_starts}{Positive integer number of optimization starts.}
\item{maxit}{Maximum optimizer iterations per run.}
\item{factr}{L-BFGS-B convergence tolerance factor.}
\item{seed}{Optional random seed. Selection derives a seed per candidate. Calls may change the random-number state.}
\item{verbose}{Print progress.}
\item{stationary}{Constrain AR reflection coefficients to [-0.9999, 0.9999].}
\item{invertible}{Constrain MA reflection coefficients to [-0.9999, 0.9999].}
\item{phi_bound}{Positive AR coefficient bound when stationary is FALSE.}
\item{theta_bound}{Positive MA coefficient bound when invertible is FALSE.}
\item{ndeps}{Positive numerical-gradient step used in optimization.}
\item{start_factr}{Positive optimizer tolerance factor for MA starting values.}
}
\details{Only Y is required; restrict the search with candidates or order ranges.
Regressors are not supported. Automatic candidates satisfy
\eqn{T \ge 2(p+q+r+1)-\min(p,q)}. IC adds \eqn{\kappa d\log(N)/N}
to the fitted loss, where d is the covariance parameter count.
Failed, nonconverged, or covariance-invalid fits receive infinite IC.
If all candidates fail, the function stops. Standard errors are not computed.
The loss is the negative normalized quasi-log-likelihood, \code{-fit$obj_val}.}
\value{For scalar kappa, a list with \code{p_hat}, \code{q_hat}, \code{r_hat},
\code{IC_min}, \code{d_selected}, \code{results} (orders, parameter counts, loss,
IC, and convergence), and \code{fit_diagnostics} (eligibility, seeds, optimizer
codes, covariance validity, warnings, bounds, and timings). For vector kappa,
a named list of these results, one per penalty.
 Additional fields are \code{neg_loglik_min} and \code{loglik_hat}.}

\examples{
set.seed(47)
Y <- matrix(rnorm(400), nrow = 100, ncol = 4)
orders <- list(c(0, 0, 0), c(1, 0, 0))  # Each triple is (p, q, r)
selected <- select_model_order_qmle(Y, candidates = orders,
    kappa = 1, n_starts = 1, seed = 123)
c(p = selected$p_hat, q = selected$q_hat, r = selected$r_hat)
}
