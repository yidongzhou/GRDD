# Network bias simulation (Monte Carlo relative bias of GRDD for network outcomes).
# Replicable RNG via base_seed; parallel over replications (doSNOW).
#
# Inputs:  none (DGP defined below).
# Outputs: output/rdata/estimation.RData containing object `results` (named list by sample size).
# Runtime: hours on a multi-core machine (500 reps x 4 sample sizes; parallel workers).
#
# Paper: Figure 1(a) — run scripts/figure1.R after this (or use shipped estimation.RData).
#
# Usage (from repository root):  Rscript scripts/simulation_estimation.R

source(file.path("R", "kerFctn.R"))
source(file.path("R", "local_linear.R"))
source(file.path("R", "lfr_com.R"))
source(file.path("R", "lfr_fun.R"))
source(file.path("R", "lfr_mea.R"))
source(file.path("R", "lfr_net.R"))
source(file.path("R", "lfr_spd.R"))
source(file.path("R", "lfr_euc.R"))
source(file.path("R", "lcm.R"))
source(file.path("R", "grdd.R"))

library(truncnorm)
library(foreach)
library(doSNOW)
library(parallel)

# DGP constants (keep in sync with scripts/simulation_inference.R)
P_BLOCK_SBM <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2)
EPS_UNIF_MAX <- 0.1

tau_ij <- function(block_i, block_j) {
  if (block_i != block_j) {
    return(0)
  }
  if (block_i == 1L) {
    return(1)
  }
  -1
}

generate_network <- function(r, n_nodes = 10) {
  n_blocks <- 2
  block_size <- n_nodes / n_blocks
  p_block <- P_BLOCK_SBM
  block_assignments <- rep(1:n_blocks, each = block_size)

  adj <- matrix(0, n_nodes, n_nodes)
  for (i in 1:(n_nodes - 1)) {
    for (j in (i + 1):n_nodes) {
      block_i <- block_assignments[i]
      block_j <- block_assignments[j]
      if (rbinom(1, 1, p_block[block_i, block_j])) {
        eps <- runif(1, min = 0, max = EPS_UNIF_MAX)
        base <- cos(pi / 2 * r)
        tau <- tau_ij(block_i, block_j)
        weight <- base + ifelse(r < 0, 0, tau) + eps
        adj[i, j] <- adj[j, i] <- weight
      }
    }
  }

  deg <- rowSums(adj)
  L <- diag(deg) - adj
  L
}

compute_true_effect <- function(n_nodes = 10) {
  n_blocks <- 2
  block_size <- n_nodes / n_blocks
  p_block <- P_BLOCK_SBM
  block_assignments <- rep(1:n_blocks, each = block_size)
  mu_eps <- EPS_UNIF_MAX / 2
  base <- 1

  E_A_minus <- matrix(0, n_nodes, n_nodes)
  E_A_plus <- matrix(0, n_nodes, n_nodes)
  for (i in 1:(n_nodes - 1)) {
    for (j in (i + 1):n_nodes) {
      bi <- block_assignments[i]
      bj <- block_assignments[j]
      p <- p_block[bi, bj]
      w_minus <- base + mu_eps
      w_plus <- base + mu_eps + tau_ij(bi, bj)
      E_A_minus[i, j] <- E_A_minus[j, i] <- p * w_minus
      E_A_plus[i, j] <- E_A_plus[j, i] <- p * w_plus
    }
  }

  deg_m <- rowSums(E_A_minus)
  deg_p <- rowSums(E_A_plus)
  L_minus <- -E_A_minus
  diag(L_minus) <- deg_m
  L_plus <- -E_A_plus
  diag(L_plus) <- deg_p

  geodesic_length <- distance(L_minus, L_plus, optns = list(type = "network"))
  list(left = L_minus, right = L_plus, geodesic_length = geodesic_length)
}

run_single_rep <- function(n, true_effect, data_seed) {
  opt_net <- list(type = "network")
  set.seed(data_seed)
  r <- rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 0.1)
  y <- lapply(r, generate_network)

  res <- grdd(y = y, x = r, cutoff = 0, optns = opt_net)

  transported <- gtm(res$tau$left, res$tau$right, true_effect$left, optns = opt_net)
  bias <- distance(true_effect$right, transported, optns = opt_net)

  bias / true_effect$geodesic_length
}

run_simulation <- function(n, n_reps = 100, true_effect, base_seed, size_idx) {
  foreach(
    i = seq_len(n_reps),
    .combine = c,
    .packages = "truncnorm",
    .options.snow = list(progress = progress)
  ) %dopar% {
    data_seed <- base_seed + (size_idx - 1L) * n_reps + i
    run_single_rep(n, true_effect, data_seed)
  }
}

progress <- function(q) {
  if (q %% 10 == 0) {
    cat(sprintf("%d runs are complete\n", q))
  }
}

base_seed <- 1L
n_reps <- 500L
sample_sizes <- c(100L, 200L, 500L, 1000L)

n_cores <- max(1L, min(50L, parallel::detectCores() - 1L))
message("Using ", n_cores, " parallel workers (adjust via machine cores).")
cl <- parallel::makeCluster(n_cores)
registerDoSNOW(cl)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

clusterExport(
  cl,
  c(
    "P_BLOCK_SBM", "EPS_UNIF_MAX",
    "generate_network", "tau_ij", "run_single_rep",
    "grdd", "distance", "gtm", "bw_select", "lfr", "lfr_net",
    "kerFctn", "local_linear_weights"
  )
)

true_effect <- compute_true_effect()
results <- Map(
  function(n_i, size_idx) {
    run_simulation(n_i, n_reps = n_reps, true_effect = true_effect, base_seed = base_seed, size_idx = size_idx)
  },
  sample_sizes,
  seq_along(sample_sizes)
)
names(results) <- paste0("n=", sample_sizes)

out <- file.path("output", "rdata", "estimation.RData")
save(results, file = out)
message("Saved: ", out)
