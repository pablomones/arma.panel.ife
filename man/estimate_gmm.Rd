\name{estimate_gmm}
\alias{estimate_gmm}
\title{GMM for ARMA Panels with a Supplied or Identity Weight}
\description{Estimate ARMA panel models with interactive effects, optionally with strictly exogenous regressors.}
\usage{
estimate_gmm(Y, p, q, r, W = NULL, heteroskedastic = TRUE, n_starts = 20L,
  maxit = 3000L, factr = 1e4, seed = NULL, verbose = FALSE, compute_se = TRUE,
  warn_dynamic = TRUE, efficient_weight = FALSE, .warn_T = TRUE, .gmm_data =
  NULL, phi_start = NULL, phi_bound = 5, theta_bound = 5, X = NULL, stationary
  = TRUE, invertible = TRUE, ndeps = 1e-5, start_factr = 1e4)
}
\arguments{
\item{Y}{Required finite numeric N by T matrix: individuals in rows, dates in columns. No missing values.}
\item{p}{Required nonnegative integer AR order, smaller than T.}
\item{q}{Required nonnegative integer MA order.}
\item{r}{Required nonnegative integer factor count, smaller than T.}
\item{W}{Optional finite symmetric positive-definite weighting matrix. Without regressors its dimension is \eqn{m \times m}, where \eqn{m=T(T+1)/2}. With regressors, \eqn{m=T(T+1)/2+kT^2}. \code{NULL} uses the identity matrix. In model selection one weight is shared by all candidates.}
\item{heteroskedastic}{Logical. If true, estimate \eqn{T-p+q} innovation variances; otherwise estimate one common innovation variance. Default is true.}
\item{n_starts}{Positive integer number of optimization starts.}
\item{maxit}{Maximum optimizer iterations per run.}
\item{factr}{L-BFGS-B convergence tolerance factor.}
\item{seed}{Optional random seed. Selection derives a seed per candidate. Calls may change the random-number state.}
\item{verbose}{Print progress.}
\item{compute_se}{Compute standard errors and confidence intervals; otherwise return missing inference values.}
\item{warn_dynamic}{Warn about fitted AR/MA roots.}
\item{efficient_weight}{Internal efficient-weight indicator; normally leave FALSE. Use estimate_gmm_efficient for automatic weighting.}
\item{.warn_T}{Internal time-dimension warning switch; normally leave TRUE.}
\item{.gmm_data}{Internal precomputed-data input; normally leave NULL.}
\item{phi_start}{Optional initial AR coefficient vector of length p.}
\item{phi_bound}{Positive AR coefficient bound when stationary is FALSE.}
\item{theta_bound}{Positive MA coefficient bound when invertible is FALSE.}
\item{X}{Optional strictly exogenous regressors: a finite numeric \eqn{N \times T} matrix for one regressor, or an \eqn{N \times T \times k} array for \eqn{k} regressors. Individual and time dimensions must match \code{Y}. \code{NULL} selects estimation without regressors.}
\item{stationary}{Constrain AR reflection coefficients to [-0.9999, 0.9999].}
\item{invertible}{Constrain MA reflection coefficients to [-0.9999, 0.9999].}
\item{ndeps}{Positive numerical-gradient step used in optimization.}
\item{start_factr}{Positive optimizer tolerance factor for MA starting values.}
}
\details{
For \eqn{t=p+1,\ldots,T}, the model is
\deqn{y_{it}=\sum_{s=1}^{p}\phi_s y_{i,t-s}+x_{it}'\beta+
\delta_t+\lambda_i'f_t+u_{it},\qquad
u_{it}=\varepsilon_{it}+\sum_{\ell=1}^{q}\theta_\ell\varepsilon_{i,t-\ell}.}{
y_it = sum_s phi_s y_i,t-s + x_it' beta + delta_t + lambda_i' f_t + u_it,
where u_it = epsilon_it + sum_l theta_l epsilon_i,t-l.}
The regressor term is omitted when \code{X = NULL}. The MA coefficients use
this plus-sign convention: the MA polynomial is
\eqn{1+\theta_1 z+\cdots+\theta_q z^q}{1 + theta_1 z + ... + theta_q z^q}, while the AR polynomial is
\eqn{1-\phi_1 z-\cdots-\phi_p z^p}{1 - phi_1 z - ... - phi_p z^p}.

The first \eqn{p} columns of \code{Y} contain the initial observations.
These observations are retained in estimation; the autoregressive equation is
imposed from period \eqn{p+1} onward. Without regressors, they are represented as
\deqn{y_{is}=\delta_s^*+\lambda_i'f_s^*+\eta_{is},\qquad s=1,\ldots,p.}{
y_is = delta_s^* + lambda_i' f_s^* + eta_is, s = 1, ..., p.}
With regressors, their representation additionally includes projections on the
complete observed regressor history. \code{Pi_pq} contains the initial-residual
covariance parameters and their permitted covariances with early innovations.
No presample observations of \eqn{y_{it}}{y_it} need to be supplied.

The factor matrix in the stacked system is normalized as
\deqn{F=(I_r',F_2')',\qquad F_2\in\mathbb{R}^{(T-r)\times r}.}{
F = rbind(I_r, F2), where F2 is a (T-r) by r matrix.}
Its first \eqn{p} rows correspond to the initial-observation projections, and
its remaining rows correspond to the factor terms in the dynamic equations.
\code{F2} stores the unrestricted part of this normalized matrix. When
\eqn{r=0}, the factor component is absent.

Only Y, p, q, and r are required; other arguments have defaults.
The sufficient time-dimension condition is \eqn{T \ge 2(p+q+r+1)-\min(p,q)}.
With X and p+r > 0, centered complete histories must have full column rank
and T-p > r. Missing data are not supported.
Check convergence, covariance validity, and inference diagnostics separately:
a converged fit can still have unavailable standard errors. All-start failure
can return missing estimates with \code{converged = FALSE}.
}
\value{A list of class \code{arma_panel_fit}, with:
\describe{
\item{\code{coef}}{A list of parameter components, including:
\describe{
\item{\code{phi}}{Vector of \eqn{p} autoregressive coefficients, ordered by lag.}
\item{\code{theta}}{Vector of \eqn{q} moving-average coefficients, ordered by
lag, using the plus-sign convention in the model equation.}
\item{\code{sigma2}}{With \code{heteroskedastic = TRUE}, a vector of
\eqn{T-p+q} innovation variances for dates \eqn{p+1-q,\ldots,T}, in
chronological order. Otherwise, one common innovation variance.}
\item{\code{F2}}{Vector of length \eqn{(T-r)r}, containing the columns of the
\eqn{(T-r)\times r}{(T-r) by r} matrix \eqn{F_2} consecutively. The complete normalized
factor matrix is \eqn{F=(I_r',F_2')'}. Empty when \eqn{r=0}.}
\item{\code{Psi_lower}, \code{Psi}}{\code{Psi_lower} is a vector of length
\eqn{r(r+1)/2}, containing the lower-triangular entries, including the diagonal,
of the Cholesky factor \eqn{L_\Psi}, read column by column. The loading
covariance matrix used in the covariance specification is
\eqn{\Psi=L_\Psi L_\Psi'}, returned as the \eqn{r\times r}{r by r} matrix \code{Psi}.
\code{Psi_lower} is empty when \eqn{r=0}.}
\item{\code{Pi_pq}}{Vector containing, in order:
\enumerate{
\item Initial-residual variances \eqn{\operatorname{Var}(\eta_{is})}{Var(eta_is)},
for \eqn{s=1,\ldots,p}.
\item Initial-residual covariances
\eqn{\operatorname{Cov}(\eta_{is},\eta_{it})}{Cov(eta_is, eta_it)}, for \eqn{1\le s<t\le p}{1 <= s < t <= p}, ordered by
increasing \eqn{s}, then increasing \eqn{t}.
\item Covariances \eqn{\operatorname{Cov}(\eta_{is},\varepsilon_{ik})}{Cov(eta_is, epsilon_ik)},
for \eqn{s=\max(1,p+1-q),\ldots,p}{s = max(1,p+1-q), ..., p} and, within each \eqn{s},
\eqn{k=p+1-q,\ldots,s}, ordered by increasing \eqn{s}, then increasing \eqn{k}.
}
Its length is \eqn{p(p+1)/2+K(p,q)}, where
\eqn{K(p,q)=h(2q-h+1)/2} and \eqn{h=\min(p,q)}{h = min(p,q)}. The vector is empty when
\eqn{p=0}; its third block is empty when \eqn{q=0}.
}
}
}
\item{\code{vartheta_hat}}{Named reported parameter vector, ordered as
\deqn{\vartheta=(\phi',\theta',(\sigma^2)',\operatorname{vech}(\Psi)',
\operatorname{vec}(F_2)',\Pi_{p,q}')'.}{
vartheta = (phi, theta, sigma2, vech(Psi), vec(F2), Pi_pq).}
Here, \eqn{\operatorname{vech}}{vech} reads the lower triangle, including the diagonal,
column by column, and \eqn{\operatorname{vec}}{vec} reads all columns consecutively.
Its covariance block contains entries of \eqn{\Psi}, whereas
\code{coef$Psi_lower} contains entries of its Cholesky factor.
}
\item{\code{psi_hat}}{Internal optimization coordinates, distinct from reported parameters.}
\item{\code{se}, \code{ci}, \code{V_hat}}{Standard errors, a two-row matrix of
normal-approximation 95 percent confidence intervals, and the estimated
asymptotic covariance of the square-root-\eqn{N} scaled estimation error.
The estimator covariance is \code{V_hat / N}, and
\code{se = sqrt(diag(V_hat) / N)}, subject to numerical validity checks.
Rows of \code{ci} are named \code{lower_95} and \code{upper_95}; the limits
are the estimates plus or minus \eqn{1.96} standard errors. Names and ordering
follow \code{vartheta_hat} without regressors and \code{parameter_hat} with
regressors, including both dimensions of \code{V_hat}. Standard errors named
\code{Psi_...} refer to entries of \eqn{\Psi}.
}
\item{\code{converged}, \code{convergence_code}}{Convergence/feasibility indicator and optimizer code.}
\item{\code{pd_components}, \code{bound_hits}, \code{dynamic_diagnostics}, \code{inference_diagnostics}}{Covariance validity, active bounds, dynamic roots, and inference diagnostics.}
\item{\code{all_runs}, \code{successful_starts}, \code{failed_starts}}{Per-start diagnostics and counts.}
\item{\code{Sigma_hat}, \code{S_N}}{Fitted covariance and sample/residual cross-product divided by N.}
}
With X, additional fields include \code{Xi_hat} (slopes and history projections),
\code{delta_hat} (time intercepts), \code{parameter_hat} (full reported vector),
\code{fitted}, and \code{residuals} (N by T). Inference refers to
\code{parameter_hat = c(vartheta_hat, Xi_hat, delta_hat)} with X and
\code{vartheta_hat} otherwise.
\code{obj_val} is the minimized GMM criterion. \code{W} stores the weight;
\code{Upsilon_hat}, \code{rank_G}, and \code{rcond_GWG} describe moment covariance
and inference conditioning. Hansen-test fields are available only when efficient
weighting and inference permit an overidentification test.
}

\examples{
set.seed(47)
Y <- matrix(rnorm(400), nrow = 100, ncol = 4)
fit <- estimate_gmm(Y, p = 0, q = 0, r = 0,
    n_starts = 1, seed = 123, compute_se = FALSE)
fit$coef$sigma2  # Estimated variance for each period
}
