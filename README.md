
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

The package [**nevada**](https://astamm.github.io/nevada/)
(NEtwork-VAlued Data Analysis) is an R package for the statistical
analysis of network-valued data. In this setting, a sample is made of
statistical units that are networks themselves. The package provides a
set of matrix representations for networks so that network-valued data
can be transformed into matrix-valued data. Subsequently, a number of
distances between matrices is provided as well to quantify how far two
networks are from each other and several test statistics are proposed
for testing equality in distribution between samples of networks using
exact permutation testing procedures. The permutation scheme is carried
out by the [**flipr**](https://astamm.github.io/flipr/) package which
also provides a number of test statistics based on inter-point distances
that play nicely with network-valued data. The implementation is largely
made in C++ and the matrix of inter- and intra-sample distances is
pre-computed, which alleviates the computational burden often associated
with permutation tests.

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

### Example 1

In this first example, we compare two populations of networks generated
according to two different models (Watts-Strogatz and Barabasi), using
the adjacency matrix representation of networks, the Frobenius distance
to compare single networks and the combination of Student-like and
Fisher-like statistics based on inter-point distances to summarize
information and perform the permutation test.

``` r
set.seed(123)
n <- 10L
x <- nevada::nvd("smallworld", n)
y <- nevada::nvd("pa", n)
```

By default the `nvd()` constructor generates networks with 25 nodes. One
can wonder whether there is a difference between the distributions that
generated these two samples (which there is given the models that we
used). The `test2_global()` function provides an answer to this
question:

``` r
t1_global <- nevada::test2_global(x, y, seed = 1234)
t1_global$pvalue
[1] 0.0009962984
```

The p-value is very small, leading to the conclusion that we should
reject the null hypothesis of equal distributions.

Although this is a fake example, we could create a partition to try to
localize differences along this partition:

``` r
partition <- as.integer(c(1:5, each = 5))
```

The `test2_local()` function provides an answer to this question:

``` r
t1_local <- nevada::test2_local(x, y, partition, seed = 1234)
t1_local
$intra
# A tibble: 5 × 3
  E     pvalue truncated
  <chr>  <dbl> <lgl>    
1 P1     0.549 TRUE     
2 P2     0.549 TRUE     
3 P3     0.549 TRUE     
4 P4     0.549 TRUE     
5 P5     0.549 TRUE     

$inter
# A tibble: 10 × 4
   E1    E2      pvalue truncated
   <chr> <chr>    <dbl> <lgl>    
 1 P1    P2    0.000996 FALSE    
 2 P1    P3    0.000996 FALSE    
 3 P1    P4    0.000996 FALSE    
 4 P1    P5    0.000996 FALSE    
 5 P2    P3    0.000996 FALSE    
 6 P2    P4    0.000996 FALSE    
 7 P2    P5    0.000996 FALSE    
 8 P3    P4    0.000996 FALSE    
 9 P3    P5    0.000996 FALSE    
10 P4    P5    0.000996 FALSE    
```

### Example 2

In this second example, we compare two populations of networks generated
according to the same model (Watts-Strogatz), using the adjacency matrix
representation of networks, the Frobenius distance to compare single
networks and the combination of Student-like and Fisher-like statistics
based on inter-point distances to summarize information and perform the
permutation test.

``` r
n <- 10L
x <- nevada::nvd("smallworld", n)
y <- nevada::nvd("smallworld", n)
```

One can wonder whether there is a difference between the distributions
that generated these two samples (which there is given the models that
we used). The `test2_global()` function provides an answer to this
question:

``` r
t2 <- nevada::test2_global(x, y, seed = 1234)
t2$pvalue
[1] 0.2257715
```

The p-value is larger than 5% or even 10%, leading us to failing to
reject the null hypothesis of equal distributions at these significance
thresholds.
