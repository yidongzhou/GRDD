# Figure 2 — UK election application: binned vote shares and local compositional fits vs. lagged margin.
#
# Inputs:  scripts/application_uk_election.R (loads CSVs; loads output/rdata/uk.RData when present
#          so grdd() is skipped; still runs inference, rdrobust, and plot prep).
# Outputs: output/figures/figure2_uk_election.pdf
#
# Usage (from repository root):  Rscript figure2.R

library(ggplot2)

source(file.path("scripts", "application_uk_election.R"))

p <- ggplot() +
  geom_line(
    data = data.frame(
      x = rep(c(seq(-0.15, 0, by = 0.01), seq(0, 0.15, by = 0.01)), each = 3),
      y = c(unlist(resuk0$predict), unlist(resuk1$predict)),
      trt = rep(c("Control", "Treat"), each = 16 * 3),
      p = rep(c("Conservative", "Labour", "Liberal"), 16 * 2)
    ),
    linewidth = 1,
    aes(x = x, y = y, group = trt, color = trt)
  ) +
  geom_point(data = avg_points, aes(x = x, y = y), alpha = 0.3, shape = 1) +
  facet_wrap(vars(p)) +
  labs(x = "Conservative margin", y = "Vote share") +
  theme_minimal() +
  theme(text = element_text(size = 15), legend.position = "none")

ggsave(file.path("output", "figures", "figure2_uk_election.pdf"), p, width = 10, height = 5)
