
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview <a href='https://astamm.github.io/nevada/'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![check-standard](https://github.com/astamm/nevada/workflows/R-CMD-check/badge.svg)](https://github.com/astamm/nevada/actions)
[![test-coverage](https://github.com/astamm/nevada/workflows/test-coverage/badge.svg)](https://github.com/astamm/nevada/actions)
[![Codecov test
coverage](https://codecov.io/gh/astamm/nevada/branch/master/graph/badge.svg)](https://codecov.io/gh/astamm/nevada?branch=master)
[![pkgdown](https://github.com/astamm/nevada/workflows/pkgdown/badge.svg)](https://github.com/astamm/nevada/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/nevada)](https://CRAN.R-project.org/package=nevada)
<!-- badges: end -->

The package `nevada` (NEtwork-VAlued Data Analysis) is an R package for
the statistical analysis of network-valued datasets. In this setting, a
sample is made of statistical units that are networks themselves. The
package provides a set of matrix representations for networks so that
network-valued data can be transformed into matrix-valued data.
Subsequently, a number of distances between matrices is provided as well
to quantify how far two networks are from each other and a number of
distance-based statistics is proposed for testing equality in
distribution between samples of networks using exact permutation testing
procedures. The implementation is largely made in C++ and the matrix of
inter- and intra-sample distances is pre-computed, which alleviates the
computational burden often associated with permutation tests. In
details:

-   the `repr_*()` functions return the chosen matrix representation of
    the input graph,
-   the `dist_*()` functions return the chosen distance between two
    networks,
-   the `stat_*()` functions return the value of the chosen test
    statistic,
-   the `test2_global()` function returns the p-value of a permutation
    test in which the null hypothesis is that the two samples come from
    the same distribution of networks,
-   the `power2()` function returns a Monte-Carlo estimate of the power
    of the test in some specific scenarios.

See the vignette *NEtwork-VAlued Data Analysis* for the details of each
function.

## Installation

You can install the latest stable version of **nevada** on CRAN with:

``` r
install.packages("nevada")
```

Or you can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("astamm/nevada")
```

## Usage

**Example 1**

In this first example, we compare two populations of networks generated
according to two different models (Watts-Strogatz and Barabasi), using
the modularity matrix representation of networks, the Hamming distance
to compare single networks and the average statistic to summarize
information and perform the permutation test.

``` r
n <- 10L
x <- nevada::nvd("smallworld", n)
y <- nevada::nvd("pa", n)
t1 <- nevada::test2_global(x, y, representation = "modularity")
t1$pvalue
[1] 0.0009936031
```

**Example 2**

In this second example, we compare two populations of networks generated
according to the sane model (Watts-Strogatz), using once again the
modularity matrix representation of networks, the Hamming distance to
compare single networks and the average statistic to summarize
information and perform the permutation test.

``` r
n <- 10L
x <- nevada::nvd("smallworld", n)
y <- nevada::nvd("smallworld", n)
t2 <- nevada::test2_global(x, y, representation = "modularity")
t2$pvalue
[1] 0.1368577
```
