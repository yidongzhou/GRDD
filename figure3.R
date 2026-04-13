# Figure 3 — Taipei supplement: estimated mean CO curves by hour (near vs. far from highway).
#
# Inputs:  output/rdata/taipei_metro.RData (from scripts/supplement_taipei_metro.R if missing).
# Outputs: output/figures/figure3_taipei_metro.pdf
#
# Usage (from repository root):  Rscript figure3.R
# Heavy step: run scripts/supplement_taipei_metro.R first if .RData is absent.

source(file.path("R", "paths.R"))
library(ggplot2)

rda <- grdd_path("output", "rdata", "taipei_metro.RData")
if (!file.exists(rda) || Sys.getenv("GRDD_FORCE_REFIT_TAIPEI", "") == "1") {
  message("Running supplement_taipei_metro.R …")
  source(grdd_path("scripts", "supplement_taipei_metro.R"))
}
load(rda, envir = .GlobalEnv)

p <- ggplot(
  data = data.frame(
    hour = rep(1:24 - 0.5, 4),
    value = c(
      resmt_highway$tau$left, resmt_highway$tau$right,
      resmt_nonhighway$tau$left, resmt_nonhighway$tau$right
    ),
    group = factor(
      rep(rep(c("Before policy", "After policy"), each = 24), 2),
      levels = c("Before policy", "After policy")
    ),
    highway = factor(
      rep(c("Near highway", "Far from highway"), each = 48),
      levels = c("Near highway", "Far from highway")
    )
  ),
  aes(x = hour, y = value, color = group)
) +
  geom_smooth(se = FALSE, span = 0.2, linewidth = 1) +
  labs(x = "Hour of day", y = "CO Concentration", color = NULL) +
  facet_wrap(vars(highway)) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  theme(text = element_text(size = 15), legend.position = "top")

ggsave(grdd_path("output", "figures", "figure3_taipei_metro.pdf"), p, width = 10, height = 5)
