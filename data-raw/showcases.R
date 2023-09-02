# Interesting simulation showing that not all repr/distances "see" the grouping
# structure
library(nevada)
n <- 100
gnp_params <- list(p = 1/3)
k_regular_params <- list(k = 8L)
x <- nvd(model = "gnp", n = n, model_params = gnp_params)
y <- nvd(model = "k_regular", n = n, model_params = k_regular_params)
z <- as_nvd(c(x, y))
mb <- c(rep(1, length(x)), rep(2, length(y)))
plot_mds  <- ggplot2::autoplot(z, memberships = mb, method = "mds")
plot_tsne <- ggplot2::autoplot(z, memberships = mb, method = "tsne")
plot_umap <- ggplot2::autoplot(z, memberships = mb, method = "umap")
