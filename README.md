
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview of the `nevada` package
--------------------------------

The package `nevada` (NEtwork-VAlued Data Analysis) contains tools for the statistical analysis of network-valued data. In this framework, the statistical unit is a network and the package deals with samples of networks. With nevada it is possible to test in a non-parametric manner if the two samples came from the same population via two test statistics based on three matrix representations for a network and four different types of distances:

-   the `get-representation` functions return the selected matrix representation of the input graph
-   the `get-distance` functions return the chosen distance between two networks
-   the `get-statistic` functions return the value of the preferred test statistic
-   the `twosample_test()` function returns the p-value of the resulting permutation test.
-   the `twosample_power()` function returns a Monte-Carlo estimate of the power of the test in some specific scenarios.

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
x <- replicate(n, igraph::sample_smallworld(dim = 1, size = 25, nei = 3, p = 0.05), simplify = FALSE)
y <- replicate(n, igraph::sample_pa(n = 25, power = 2, m = 3, directed = FALSE), simplify = FALSE)

test1 <- nevada::twosample_test(x, y, "modularity")
#>  - P-value resolution: 0.001
#>  - Computing approximate p-value using 1000 random permutations.
#>  - P-value will not drop below 5.41254411223451e-06 on average over repeated Monte-Carlo estimates.
test1$pvalue
#> [1] 0.0009936031
```

**Example 2**

In this second example, we compare two populations of networks generated according to the sane model (Watts-Strogatz), using once again the modularity matrix representation of networks, the Hamming distance to compare single networks and the average statistic to summarize information and perform the permutation test.

``` r
n <- 10L
x <- replicate(n, igraph::sample_smallworld(dim = 1, size = 25, nei = 3, p = 0.05), simplify = FALSE)
y <- replicate(n, igraph::sample_smallworld(dim = 1, size = 25, nei = 3, p = 0.05), simplify = FALSE)

test2 <- nevada::twosample_test(x, y, "modularity")
#>  - P-value resolution: 0.001
#>  - Computing approximate p-value using 1000 random permutations.
#>  - P-value will not drop below 5.41254411223451e-06 on average over repeated Monte-Carlo estimates.
test2$pvalue
#> [1] 0.5114831
```
