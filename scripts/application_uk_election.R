# UK constituency application: compositional vote shares vs. lagged Conservative margin (RDD at 0).
# Compares GRDD (composition) to separate rdrobust runs per party.
#
# Inputs:  data/incumbency_advantage_* CSV files (Eggers & Spirling replication).
# Outputs: Printed tables; writes output/rdata/uk.RData with resuk, resuk0, resuk1, avg_points
#          (for fast figure2.R). If that file exists, resuk is loaded from disk.
#
# Usage (from repository root):  source("scripts/application_uk_election.R")

source(file.path("R", "kerFctn.R"))
source(file.path("R", "ll_weights.R"))
source(file.path("R", "lfr_com.R"))
source(file.path("R", "lfr_fun.R"))
source(file.path("R", "lfr_mea.R"))
source(file.path("R", "lfr_net.R"))
source(file.path("R", "lfr_spd.R"))
source(file.path("R", "lfr_euc.R"))
source(file.path("R", "lcm.R"))
source(file.path("R", "grdd.R"))
source(file.path("R", "grdd_inference.R"))

library(dplyr)
library(readr)
library(lubridate)
library(ggplot2)
library(rdrobust)
library(tidyr)

data_dir <- file.path("data")
con_df <- read.csv(file.path(data_dir, "incumbency_advantage_con_based_all_20150122_including_unionists.csv"))
lab_df <- read.csv(file.path(data_dir, "incumbency_advantage_lab_based_all_20150122.csv"))
lib_df <- read.csv(file.path(data_dir, "incumbency_advantage_lib_based_all_20150122.csv"))

con_votes <- con_df %>%
  select(year = `e.year`, month = `e.month`, cid = `c.id`, con_votes = ref_votes)

lab_votes <- lab_df %>%
  select(year = `e.year`, month = `e.month`, cid = `c.id`, lab_votes = ref_votes)

lib_votes <- lib_df %>%
  select(year = `e.year`, month = `e.month`, cid = `c.id`, lib_votes = ref_votes)

votes_all <- con_votes %>%
  inner_join(lab_votes, by = c("year", "month", "cid")) %>%
  inner_join(lib_votes, by = c("year", "month", "cid")) %>%
  filter(con_votes > 0 & lab_votes > 0 & lib_votes > 0) %>%
  mutate(
    total_votes = con_votes + lab_votes + lib_votes,
    election_date = make_date(year, month, 1),
    con_share = con_votes / total_votes,
    lab_share = lab_votes / total_votes,
    lib_share = lib_votes / total_votes,
    top_rival_votes = pmax(lab_votes, lib_votes),
    margin_pct = (con_votes - top_rival_votes) / total_votes
  )

votes_all <- votes_all %>%
  arrange(cid, election_date) %>%
  group_by(cid) %>%
  mutate(prev_margin = lag(margin_pct)) %>%
  ungroup()

final_data <- votes_all %>%
  filter(!is.na(prev_margin)) %>%
  select(R = prev_margin, Conservative = con_share, Labour = lab_share, Liberal = lib_share)

Ruk <- final_data$R
Yuk <- lapply(seq_len(nrow(final_data)), function(i) {
  Yi <- as.numeric(final_data[i, c("Conservative", "Labour", "Liberal")])
  names(Yi) <- c("Conservative", "Labour", "Liberal")
  Yi
})

Y_mat <- do.call(rbind, Yuk)
colnames(Y_mat) <- c("Conservative", "Labour", "Liberal")

comp_tol <- 1e-10
bad_comp <- vapply(Yuk, function(v) abs(sum(v) - 1) > comp_tol, logical(1))
if (any(bad_comp)) {
  warning("Some compositions do not sum to 1; check data merge.")
}

uk_rda <- file.path("output", "rdata", "uk.RData")

if (file.exists(uk_rda)) {
  message("Loading cached objects from ", uk_rda)
  load(uk_rda)
}
if (!exists("resuk")) {
  message("Fitting GRDD (composition) …")
  resuk <- grdd(y = Yuk, x = Ruk, cutoff = 0, optns = list(type = "composition"))
}

distance(resuk$tau$right, resuk$tau$left, optns = list(type = "composition"))
set.seed(1)
pvuk <- grdd_inference(resuk, alpha = 0.05, B = 5000, seed = 1)
print(pvuk)

grdd_jump <- as.numeric(resuk$tau$right - resuk$tau$left)
names(grdd_jump) <- colnames(Y_mat)

rd_extract <- function(rd, i = 2L) {
  coef_v <- as.numeric(rd$coef)
  se_v <- as.numeric(rd$se)
  pv_v <- as.numeric(rd$pv)
  ci_mat <- as.matrix(rd$ci)
  c(
    estimate = coef_v[i],
    se = se_v[i],
    pv = pv_v[i],
    ci_l = ci_mat[i, 1L],
    ci_u = ci_mat[i, 2L]
  )
}

party_names <- colnames(Y_mat)
rd_scalar <- vector("list", length(party_names))
names(rd_scalar) <- party_names

for (pn in party_names) {
  yj <- Y_mat[, pn]
  rd_scalar[[pn]] <- rdrobust(yj, Ruk, c = 0)
}

uk_rd_compare <- data.frame(
  party = party_names,
  grdd_tau_left = resuk$tau$left[party_names],
  grdd_tau_right = resuk$tau$right[party_names],
  grdd_tau = grdd_jump[party_names],
  grdd_bw = rep(resuk$optns$bw, length(party_names)),
  rd_tau_bc_left = vapply(rd_scalar, function(rd) as.numeric(rd$tau_bc)[1], numeric(1)),
  rd_tau_bc_right = vapply(rd_scalar, function(rd) as.numeric(rd$tau_bc)[2], numeric(1)),
  rd_h_left = vapply(rd_scalar, function(rd) as.numeric(rd$bws["h", "left"]), numeric(1)),
  rd_h_right = vapply(rd_scalar, function(rd) as.numeric(rd$bws["h", "right"]), numeric(1)),
  rd_b_left = vapply(rd_scalar, function(rd) as.numeric(rd$bws["b", "left"]), numeric(1)),
  rd_b_right = vapply(rd_scalar, function(rd) as.numeric(rd$bws["b", "right"]), numeric(1)),
  rd_tau = vapply(rd_scalar, function(rd) rd_extract(rd, 2L)["estimate"], numeric(1)),
  rd_se = vapply(rd_scalar, function(rd) rd_extract(rd, 2L)["se"], numeric(1)),
  rd_pv = vapply(rd_scalar, function(rd) rd_extract(rd, 2L)["pv"], numeric(1)),
  rd_ci_l = vapply(rd_scalar, function(rd) rd_extract(rd, 2L)["ci_l"], numeric(1)),
  rd_ci_u = vapply(rd_scalar, function(rd) rd_extract(rd, 2L)["ci_u"], numeric(1)),
  stringsAsFactors = FALSE
)
print(uk_rd_compare)

resuk0 <- lfr_com(
  y = Yuk[Ruk < 0], x = Ruk[Ruk < 0], xOut = seq(-0.15, 0, by = 0.01),
  optns = list(bw = resuk$optns$bw, kernel = "triangular")
)
resuk1 <- lfr_com(
  y = Yuk[Ruk >= 0], x = Ruk[Ruk >= 0], xOut = seq(0, 0.15, by = 0.01),
  optns = list(bw = resuk$optns$bw, kernel = "triangular")
)

bin_df <- data.frame(R = Ruk, con = Y_mat[, 1], lab = Y_mat[, 2], lib = Y_mat[, 3])
bin_df <- bin_df %>% filter(R >= -0.15, R <= 0.15)
bin_df$bin <- cut(bin_df$R, breaks = seq(-0.15, 0.15, by = 0.001), include.lowest = TRUE)

avg_points <- bin_df %>%
  group_by(bin) %>%
  summarise(
    x = mean(R),
    con = mean(con),
    lab = mean(lab),
    lib = mean(lib),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c("con", "lab", "lib"), names_to = "p", values_to = "y") %>%
  mutate(p = recode(p, con = "Conservative", lab = "Labour", lib = "Liberal"))

save(resuk, resuk0, resuk1, avg_points, file = uk_rda)
message("Saved ", uk_rda, " (resuk + plot inputs for figure2.R).")
