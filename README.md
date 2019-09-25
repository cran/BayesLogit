
[![Build
Status](https://travis-ci.org/tmsalab/BayesLogit.svg)](https://travis-ci.org/tmsalab/BayesLogit)
[![Package-License](http://img.shields.io/badge/license-GPL%20\(%3E=3\)-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
[![CRAN Version
Badge](http://www.r-pkg.org/badges/version/BayesLogit)](https://cran.r-project.org/package=BayesLogit)
[![CRAN
Status](https://cranchecks.info/badges/worst/BayesLogit)](https://cran.r-project.org/web/checks/check_results_BayesLogit.html)
[![RStudio CRAN Mirror’s Monthly
Downloads](http://cranlogs.r-pkg.org/badges/BayesLogit?color=brightgreen)](http://www.r-pkg.org/pkg/BayesLogit)
[![RStudio CRAN Mirror’s Total
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/BayesLogit?color=brightgreen)](http://www.r-pkg.org/pkg/BayesLogit)

# `BayesLogit` R package

Perform posterior simulation for binomial and multinomial logistic
regression using the Polya-Gamma latent variable technique. This method
is fully automatic, exact, and fast. A routine to efficiently sample
from the Polya-Gamma class of distributions is included.

## Installation

You can install `BayesLogit` from CRAN using:

``` r
install.packages("BayesLogit")
```

Or, you can be on the cutting-edge development version on GitHub using:

``` r
if(!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("tmsalab/BayesLogit")
```

## Usage

To use the `BayesLogit` package, load it into *R* using:

``` r
library("BayesLogit")
```

## Authors

Nicholas G. Polson, James G. Scott, and Jesse Windle

Fixes to relist on CRAN: James Joseph Balamuta

## License

GPL (\>= 3)
