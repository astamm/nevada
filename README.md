
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/astamm/nevada.svg?branch=master)](https://travis-ci.org/astamm/nevada)

Overview of the `nevada` package
--------------------------------

The package `nevada` (NEtwork-VAlued Data Analysis) contains tools for the statistical analysis of network-valued data. In this framework, the statistical unit is a network and the package deals with samples of networks. With nevada it is possible to test in a non-parametric manner if the two samples came from the same population via two test statistics based on three matrix representations for a network and four different types of distances:

-   the `repr_*()` functions return the selected matrix representation of the input graph
-   the `dist_*()` functions return the chosen distance between two networks
-   the `stat_*()` functions return the value of the preferred test statistic
-   the `test_twosample()` function returns the p-value of the resulting permutation test.
-   the `power_twosample()` function returns a Monte-Carlo estimate of the power of the test in some specific scenarios.

See the vignette *NEtwork-VAlued Data Analysis* for the details of each function.

Installation
------------

You can install `nevada` from github with:

``` r
# install.packages("devtools")
devtools::install_github("ilovato/nevada")
```

It relies on the `igraph` package. If you encounter bugs or for questions and comments, please contact the maintainer of the package.

Usage
-----

**Example 1**

In this first example, we compare two populations of networks generated according to two different models (Watts-Strogatz and Barabasi), using the modularity matrix representation of networks, the Hamming distance to compare single networks and the average statistic to summarize information and perform the permutation test.

``` r
n <- 10L
x <- nevada::nvd("smallworld", n)
y <- nevada::nvd("pa", n)
test1 <- nevada::test_twosample(x, y, "modularity")
test1$pvalue
#> [1] 0.0009936031
```

**Example 2**

In this second example, we compare two populations of networks generated according to the sane model (Watts-Strogatz), using once again the modularity matrix representation of networks, the Hamming distance to compare single networks and the average statistic to summarize information and perform the permutation test.

``` r
n <- 10L
x <- nevada::nvd("smallworld", n)
y <- nevada::nvd("smallworld", n)
test2 <- nevada::test_twosample(x, y, "modularity")
test2$pvalue
#> [1] 0.415579
```
