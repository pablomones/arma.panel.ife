# arma.panel.ife

QMLE, GMM, efficient GMM, and model-order selection for ARMA panels with
interactive effects. Version 0.1.0.

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

The estimators accept optional strictly exogenous regressors. Order selection
does not accept regressors. Each function has an R help page.

License: GPL-3 (GNU General Public License version 3).
