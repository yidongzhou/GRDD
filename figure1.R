# Figure 1 — Simulation results: (a) relative bias of network GRDD; (b) empirical rejection under H0 / power.
#
# Inputs:  output/rdata/estimation.RData (`results`), output/rdata/inference.RData (`reject_array`).
#          Shipped files reproduce paper plots without re-running simulations.
# Outputs: output/figures/figure1a_simulation_estimation.pdf, figure1b_simulation_inference.pdf
#          Console: summary tables (lm for bias scaling; rejection rates vs delta).
#
# Seeds: fixed inside simulation scripts (not used here — only loads saved MC output).
#
# Usage (from repository root):  Rscript figure1.R

source(file.path("R", "paths.R"))
library(ggplot2)

load(grdd_path("output", "rdata", "estimation.RData"))
sample_sizes <- as.integer(gsub("^n=", "", names(results)))
n_reps <- length(results[[1L]])

plot_data <- data.frame(
  n = rep(sample_sizes, each = n_reps),
  bias = unlist(results)
)
p_sim <- ggplot(plot_data, aes(x = factor(n), y = bias, fill = factor(n))) +
  geom_boxplot(alpha = 0.5) +
  labs(x = "Sample size", y = "Relative bias") +
  theme_minimal() +
  theme(text = element_text(size = 20)) +
  guides(fill = "none")

ggsave(grdd_path("output", "figures", "figure1a_simulation_estimation.pdf"), p_sim, width = 8, height = 6)

ab <- sapply(results, mean)
fitn <- lm(log(ab) ~ log(sample_sizes))
cat("\n--- simulation_estimation: lm(log(mean bias) ~ log(n)) ---\n")
print(summary(fitn))

load(grdd_path("output", "rdata", "inference.RData"))
alpha <- 0.05
B_boot <- 1000L
n_reps <- dim(reject_array)[3L]
sample_sizes <- as.integer(dimnames(reject_array)$n)
delta_grid <- as.numeric(dimnames(reject_array)$delta)

empirical_reject <- apply(reject_array, c(1, 2), mean, na.rm = TRUE)
se <- sqrt(empirical_reject * (1 - empirical_reject) / n_reps)
n_missing <- apply(is.na(reject_array), c(1, 2), sum)
summary_tbl <- expand.grid(n = sample_sizes, delta = delta_grid)
summary_tbl$prop_reject <- as.vector(empirical_reject)
summary_tbl$se <- as.vector(se)
summary_tbl$n_reps <- n_reps
summary_tbl$n_missing <- as.vector(n_missing)

cat("\n--- simulation_inference: network bootstrap test ---\n")
cat("alpha = ", alpha, ", B = ", B_boot, "\n", sep = "")
print(summary_tbl, row.names = FALSE)

plot_data <- summary_tbl
plot_data$n <- factor(plot_data$n)
p_test <- ggplot(plot_data, aes(x = delta, y = prop_reject, color = n, linetype = n, group = n)) +
  geom_hline(yintercept = alpha, color = "gray50", linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = expression(delta), y = "Rejection probability", color = "Sample size", linetype = "Sample size") +
  theme_minimal() +
  theme(text = element_text(size = 20), legend.position = "top")

ggsave(grdd_path("output", "figures", "figure1b_simulation_inference.pdf"), p_test, width = 8, height = 6)
