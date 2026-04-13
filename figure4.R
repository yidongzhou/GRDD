# Figure 4 — Taipei supplement: interactive 3D surfaces of predicted CO (near vs. far stations).
#
# Inputs:  output/rdata/taipei_metro.RData (same as figure3.R).
# Outputs: output/figures/co_surface_near_highway.html, co_surface_far_highway.html
#
# Usage (from repository root):  Rscript figure4.R

rda <- file.path("output", "rdata", "taipei_metro.RData")
if (!file.exists(rda)) {
  message("Running supplement_taipei_metro.R …")
  source(file.path("scripts", "supplement_taipei_metro.R"))
}
load(rda, envir = .GlobalEnv)

htmlwidgets::saveWidget(p_near, file.path("output", "figures", "co_surface_near_highway.html"))
htmlwidgets::saveWidget(p_far, file.path("output", "figures", "co_surface_far_highway.html"))
