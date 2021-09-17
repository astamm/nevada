n <- 10L

# Two different models for the two populations
x <- nvd("pa", n)
y <- nvd("pa", n)

t1 <- test2_global(x, y, stats = "student_euclidean", seed = 1234)
t1$pvalue
t2 <- test2_global(x, y, stats = "welch_euclidean", seed = 1234)
t2$pvalue
t3 <- test2_global(x, y, stats = "original_edge_count", seed = 1234)
t3$pvalue
t4 <- test2_global(x, y, stats = "generalized_edge_count", seed = 1234)
t4$pvalue
t5 <- test2_global(x, y, stats = "weighted_edge_count", seed = 1234)
t5$pvalue
t6 <- test2_global(x, y, stats = "flipr:student_ip", seed = 1234)
t6$pvalue
t7 <- test2_global(x, y, stats = "flipr:fisher_ip", seed = 1234)
t7$pvalue
t8 <- test2_global(x, y, stats = "flipr:bg_ip", seed = 1234)
t8$pvalue
t9 <- test2_global(x, y, stats = "flipr:energy_ip", seed = 1234)
t9$pvalue
t10 <- test2_global(x, y, stats = "flipr:cq_ip", seed = 1234)
t10$pvalue

test2_global(x, y, stats = c("student_euclidean", "welch_euclidean"), seed = 1234)$pvalue
test2_global(x, y, stats = c("student_euclidean", "generalized_edge_count"), seed = 1234)$pvalue
test2_global(x, y, stats = c("flipr:bg_ip", "flipr:cq_ip"), seed = 1234)$pvalue
test2_global(x, y, stats = c("original_edge_count", "generalized_edge_count"), seed = 1234)$pvalue
test2_global(x, y, stats = c("flipr:bg_ip", "generalized_edge_count"), seed = 1234)$pvalue
