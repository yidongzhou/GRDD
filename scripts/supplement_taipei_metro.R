# Supplement: Taipei air quality (Chen & Whalley AEJ applied paper data).
# Functional GRDD on smoothed hourly CO curves around a policy date (RDD cutoff).
#
# Inputs:  data/all_full_report_stations_94_07_hourly.dta (~200 MB; not in git — download from
#          https://doi.org/10.3886/E114778V1 , path data_programs/data/taipei/ in the archive).
# Outputs: output/rdata/taipei_metro.RData with fitted objects and plotly figures for supplement plots.
# Runtime: long (smoothing + GRDD + inference + local fits + 3D plotly).
#
# Paper: Figures 3–4 — run figure3.R / figure4.R after this file, or rely on saved .RData.
#
# Usage (from repository root):  Rscript scripts/supplement_taipei_metro.R

source(file.path("R", "paths.R"))
source(grdd_path("R", "kerFctn.R"))
source(grdd_path("R", "local_linear.R"))
source(grdd_path("R", "lfr_com.R"))
source(grdd_path("R", "lfr_fun.R"))
source(grdd_path("R", "lfr_mea.R"))
source(grdd_path("R", "lfr_net.R"))
source(grdd_path("R", "lfr_spd.R"))
source(grdd_path("R", "lfr_euc.R"))
source(grdd_path("R", "lcm.R"))
source(grdd_path("R", "grdd.R"))
source(grdd_path("R", "grdd_inference.R"))

library(haven)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(plotly)
library(lubridate)
library(fdapace)

dta_path <- grdd_path("data", "all_full_report_stations_94_07_hourly.dta")
if (!file.exists(dta_path)) {
  stop(
    "Missing ", dta_path, ".\n",
    "Download Chen & Whalley replication from https://doi.org/10.3886/E114778V1 ",
    "and copy data_programs/data/taipei/all_full_report_stations_94_07_hourly.dta to data/.",
    call. = FALSE
  )
}

taipei <- read_dta(dta_path)

cutoff_date <- as.Date("1996-03-28")
start_date <- as.Date("1995-03-28")
end_date <- as.Date("1997-03-28")

smooth_station_day_lwls <- function(xin, yin, xout = 1:24, bw = 1, kernel = "gauss") {
  if (length(xin) < 4 || length(unique(xin)) < 2) {
    return(rep(NA, length(xout)))
  }
  ord <- order(xin)
  xin <- xin[ord]
  yin <- yin[ord]
  tryCatch(
    Lwls1D(
      bw = bw,
      kernel_type = kernel,
      xin = xin,
      yin = yin,
      xout = xout
    ),
    error = function(e) rep(NA, length(xout))
  )
}

taipei_filtered <- taipei %>%
  mutate(date = as.Date(Date, origin = "1960-01-01")) %>%
  filter(date >= start_date & date <= end_date) %>%
  mutate(group = ifelse(station_highr == 1, "Near highway", "Far from highway")) %>%
  select(group, station, date, Hour, CO)

station_day_list <- taipei_filtered %>%
  group_by(group, station, date) %>%
  group_split()

smoothed_curves <- map(station_day_list, function(df) {
  xin <- df$Hour
  yin <- df$CO
  smooth_station_day_lwls(xin, yin)
})

smoothed_list <- map2(
  station_day_list,
  smoothed_curves,
  function(meta_df, smoothed_vec) {
    tibble(
      group = unique(meta_df$group),
      station = unique(meta_df$station),
      date = unique(meta_df$date),
      smoothed = list(smoothed_vec)
    )
  }
)

curve_df <- bind_rows(smoothed_list) %>%
  unnest_wider(smoothed, names_sep = "") %>%
  rename_with(~ paste0("h", 1:24), starts_with("smoothed"))

group_avg_curves <- curve_df %>%
  group_by(group, date) %>%
  summarise(across(starts_with("h"), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(running_var = as.numeric(date - as.Date("1996-03-28")))

Ymt_highway <- group_avg_curves %>%
  filter(group == "Near highway") %>%
  select(starts_with("h")) %>%
  asplit(1)

Rmt_highway <- group_avg_curves %>%
  filter(group == "Near highway") %>%
  pull(running_var)

Ymt_nonhighway <- group_avg_curves %>%
  filter(group == "Far from highway") %>%
  select(starts_with("h")) %>%
  asplit(1)

Rmt_nonhighway <- group_avg_curves %>%
  filter(group == "Far from highway") %>%
  pull(running_var)

resmt_highway <- grdd(
  y = Ymt_highway, x = Rmt_highway, cutoff = 0,
  optns = list(type = "function", supp = 1:24 - 0.5, bw_optns = list(bw_max = 60))
)
resmt_nonhighway <- grdd(
  y = Ymt_nonhighway, x = Rmt_nonhighway, cutoff = 0,
  optns = list(type = "function", supp = 1:24 - 0.5, bw_optns = list(bw_max = 60))
)

distance(resmt_highway$tau$right, resmt_highway$tau$left, optns = list(type = "function", supp = 1:24 - 0.5))
distance(resmt_nonhighway$tau$right, resmt_nonhighway$tau$left, optns = list(type = "function", supp = 1:24 - 0.5))

set.seed(1)
pv_highway <- grdd_inference(resmt_highway, alpha = 0.05, B = 1000, seed = 1)
pv_nonhighway <- grdd_inference(resmt_nonhighway, alpha = 0.05, B = 1000, seed = 1)
print(pv_highway)
print(pv_nonhighway)

res_highway_before <- lfr_fun(
  y = Ymt_highway[Rmt_highway < 0],
  x = Rmt_highway[Rmt_highway < 0],
  xOut = seq(-100, 0, by = 2),
  optns = list(bw = resmt_highway$optns$bw, kernel = "triangular")
)

res_highway_after <- lfr_fun(
  y = Ymt_highway[Rmt_highway >= 0],
  x = Rmt_highway[Rmt_highway >= 0],
  xOut = seq(0, 100, by = 2),
  optns = list(bw = resmt_highway$optns$bw, kernel = "triangular")
)

res_nonhighway_before <- lfr_fun(
  y = Ymt_nonhighway[Rmt_nonhighway < 0],
  x = Rmt_nonhighway[Rmt_nonhighway < 0],
  xOut = seq(-100, 0, by = 2),
  optns = list(bw = resmt_nonhighway$optns$bw, kernel = "triangular")
)

res_nonhighway_after <- lfr_fun(
  y = Ymt_nonhighway[Rmt_nonhighway >= 0],
  x = Rmt_nonhighway[Rmt_nonhighway >= 0],
  xOut = seq(0, 100, by = 2),
  optns = list(bw = resmt_nonhighway$optns$bw, kernel = "triangular")
)

df_highway <- data.frame(
  hour = rep(1:24, 102) - 0.5,
  date = rep(c(seq(-100, 0, by = 2), seq(0, 100, by = 2)) + cutoff_date, each = 24),
  value = c(unlist(res_highway_before$predict), unlist(res_highway_after$predict)),
  id = rep(1:102, each = 24),
  group = rep(c("Before policy", "After policy"), each = 24 * 51),
  type = "Near highway"
)

df_nonhighway <- data.frame(
  hour = rep(1:24, 102) - 0.5,
  date = rep(c(seq(-100, 0, by = 2), seq(0, 100, by = 2)) + cutoff_date, each = 24),
  value = c(unlist(res_nonhighway_before$predict), unlist(res_nonhighway_after$predict)),
  id = rep(1:102, each = 24),
  group = rep(c("Before policy", "After policy"), each = 24 * 51),
  type = "Far from highway"
)

zlim <- range(bind_rows(df_highway, df_nonhighway)$value, na.rm = TRUE)

p_near <- plot_ly() %>%
  add_surface(
    x = rep(1:24, 102) - 0.5,
    y = seq(-100, 0, by = 2) + cutoff_date,
    z = matrix(df_highway$value[df_highway$group == "Before policy"], ncol = 24, byrow = TRUE),
    showscale = FALSE,
    colorscale = list(list(0, "#F8766D"), list(1, "#F8766D")),
    name = "Before Policy"
  ) %>%
  add_surface(
    x = rep(1:24, 102) - 0.5,
    y = seq(0, 100, by = 2) + cutoff_date,
    z = matrix(df_highway$value[df_highway$group == "After policy"], ncol = 24, byrow = TRUE),
    showscale = FALSE,
    colorscale = list(list(0, "#00BFC4"), list(1, "#00BFC4")),
    name = "After Policy"
  ) %>%
  add_surface(
    x = rep(1:24, 2) - 0.5,
    y = rep(cutoff_date, 2),
    z = matrix(c(
      df_highway$value[df_highway$group == "Before policy" & df_highway$date == cutoff_date],
      df_highway$value[df_highway$group == "After policy" & df_highway$date == cutoff_date]
    ), nrow = 2, byrow = TRUE),
    showscale = FALSE,
    colorscale = list(list(0, "#404040"), list(1, "#404040")),
    name = "Treatment Boundary"
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "Hour of Day"),
      yaxis = list(title = ""),
      zaxis = list(title = "CO Concentration", range = zlim),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    ),
    title = list(text = "CO Concentration Surface for Near Highway Stations", y = 0.95)
  )

p_far <- plot_ly() %>%
  add_surface(
    x = rep(1:24, 102) - 0.5,
    y = seq(-100, 0, by = 2) + cutoff_date,
    z = matrix(df_nonhighway$value[df_nonhighway$group == "Before policy"], ncol = 24, byrow = TRUE),
    showscale = FALSE,
    colorscale = list(list(0, "#F8766D"), list(1, "#F8766D")),
    name = "Before Policy"
  ) %>%
  add_surface(
    x = rep(1:24, 102) - 0.5,
    y = seq(0, 100, by = 2) + cutoff_date,
    z = matrix(df_nonhighway$value[df_nonhighway$group == "After policy"], ncol = 24, byrow = TRUE),
    showscale = FALSE,
    colorscale = list(list(0, "#00BFC4"), list(1, "#00BFC4")),
    name = "After Policy"
  ) %>%
  add_surface(
    x = rep(1:24, 2) - 0.5,
    y = rep(cutoff_date, 2),
    z = matrix(c(
      df_nonhighway$value[df_nonhighway$group == "Before policy" & df_nonhighway$date == cutoff_date],
      df_nonhighway$value[df_nonhighway$group == "After policy" & df_nonhighway$date == cutoff_date]
    ), nrow = 2, byrow = TRUE),
    showscale = FALSE,
    colorscale = list(list(0, "#404040"), list(1, "#404040")),
    name = "Treatment Boundary"
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "Hour of Day"),
      yaxis = list(title = ""),
      zaxis = list(title = "CO Concentration", range = zlim),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    ),
    title = list(text = "CO Concentration Surface for Far from Highway Stations", y = 0.95)
  )

out_rda <- grdd_path("output", "rdata", "taipei_metro.RData")
save(
  resmt_highway,
  resmt_nonhighway,
  df_highway,
  df_nonhighway,
  cutoff_date,
  zlim,
  p_near,
  p_far,
  Ymt_highway,
  Rmt_highway,
  Ymt_nonhighway,
  Rmt_nonhighway,
  pv_highway,
  pv_nonhighway,
  file = out_rda
)
message("Saved: ", out_rda)
