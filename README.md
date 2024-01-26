
<!-- README.md is generated from README.Rmd. Please edit this file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/yboulag/cTOST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yboulag/cTOST/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# `cTOST` Overview

This repository holds the `cTOST` R package. This package contains the
function `tost` which provides an assessment of equivalence in the
univariate framework based on the state-of-the-art Two One-Sided Tests
(TOST). In addition, the package contains the functions `atost` and
`dtost`, two corrective procedures applied to the TOST in the univariate
framework in order to ensure the preservation of the Type I error rate
at the desired nominal level and a uniform increase in power. These two
functions output an assessment of equivalence in the univariate
framework after their respective corrections is applied. More details
can be found in Boulaguiem et al. (2023) that you can access via this
[link](https://onlinelibrary.wiley.com/doi/10.1002/sim.9993).

# Install Instructions

The `cTOST` package is available on GitHub at the moment. It is subject
to ongoing updates that may lead to stability issues.

In order to install the package, it is required to pre-install the
`devtools` dependency. Run the following command if you do not have it
already installed:

``` r
install.packages("devtools")
```

The package is then installed with the following command:

``` r
devtools::install_github("yboulag/cTOST")
```

Note that Windows users are assumed that have Rtools installed (if this
is not the case, please visit this
[link](https://cran.r-project.org/bin/windows/Rtools/).

# How to cite

``` r
@Manual{boulaguiem2023ctost,
  title = {cTOST: Finite Sample Correction of The TOST in The Univariate Framework},
  author = {Younes Boulaguiem and Stéphane Guerrier and Dominique-Laurent Couturier},
  year = {2023},
  note = {R package version 1.0.0},
  url = {https://github.com/yboulag/cTOST},
}
```

# License

The license this source code is released under is the GNU AFFERO GENERAL
PUBLIC LICENSE (AGPL) v3.0. Please see the LICENSE file for full text.
Otherwise, please consult
[GNU](https://www.gnu.org/licenses/agpl-3.0.en.html) which will provide
a synopsis of the restrictions placed upon the code.

# References

Boulaguiem, Younes, Julie Quartier, Maria Lapteva, Yogeshvar N Kalia,
Maria-Pia Victoria-Feser, Stéphane Guerrier, and Dominique-Laurent
Couturier. 2023. “Finite Sample Adjustments for Average Equivalence
Testing”, *Statistics in Medicine*,
<https://onlinelibrary.wiley.com/doi/10.1002/sim.9993>.
