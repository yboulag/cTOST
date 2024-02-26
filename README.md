
<!-- README.md is generated from README.Rmd. Please edit this file -->

# `cTOST` Overview <a href="https://yboulag.github.io/cTOST/"><img src="man/figures/hex-cTOST.png" alt="" align="right" height="138" width="125" /></a>

<!-- badges: start -->

[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--02--26-green.svg)](https://github.com/yboulag/cTOST)
[![R-CMD-check](https://github.com/yboulag/cTOST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yboulag/cTOST/actions/workflows/R-CMD-check.yaml)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/cTOST)](https://www.r-pkg.org/pkg/cTOST)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/cTOST)](https://www.r-pkg.org/pkg/cTOST)
<!-- badges: end -->

## Overview

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
[link](https://www.biorxiv.org/content/10.1101/2023.03.11.532179v3).

## Install Instructions

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

## How to cite

    @Manual{boulaguiem2023ctost,
      title = {cTOST: Finite Sample Correction of The TOST in The Univariate Framework},
      author = {Younes Boulaguiem and Stéphane Guerrier and Dominique-Laurent Couturier},
      year = {2023},
      note = {R package version 1.0.0},
      url = {https://github.com/yboulag/cTOST},
    }

## License

The license this source code is released under is the GNU AFFERO GENERAL
PUBLIC LICENSE (AGPL) v3.0. Please see the LICENSE file for full text.
Otherwise, please consult
[GNU](https://www.gnu.org/licenses/agpl-3.0.en.html) which will provide
a synopsis of the restrictions placed upon the code.

## References

Boulaguiem, Y., Quartier, J., Lapteva, M., Kalia, Y. N., Victoria-Feser,
M. P., Guerrier, S., & Couturier, D. L., “*Finite Sample Adjustments for
Average Equivalence Testing*”, Statistics in Medicine, 2023,
<https://doi.org/10.1002/sim.9993>.
