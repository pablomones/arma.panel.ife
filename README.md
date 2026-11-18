# arma.panel.ife

`arma.panel.ife` provides QMLE, GMM, and efficient GMM estimation and
inference for short ARMA panels with interactive effects, together
with joint selection of the ARMA orders and number of factors.

This package implements the methods developed in *Likelihood and
Order Selection for ARMA Dynamics in Panel Models with Interactive
Effects* by Jushan Bai and Pablo Mones (2026).

Version 0.1.0.
## Installation

Install the latest version from GitHub:

```r
install.packages("remotes")
remotes::install_github("pablomones/arma.panel.ife")
library(arma.panel.ife)
```

The required dependencies, including `numDeriv`, are installed automatically
if needed.

To download the source code without installing it, visit the
[GitHub repository](https://github.com/pablomones/arma.panel.ife) and choose
**Code > Download ZIP**, then extract the downloaded archive.

## Functions

- `estimate_qmle`: quasi-maximum likelihood.
- `estimate_gmm`: GMM with a supplied or identity weight.
- `estimate_gmm_efficient`: GMM with a data-based efficient weight.
- `select_model_order_qmle`: QMLE-based order selection.
- `select_model_order_gmm`: GMM-based order selection.

The estimators accept optional strictly exogenous regressors. Each function has an R help page.

License: GPL-3 (GNU General Public License version 3).
