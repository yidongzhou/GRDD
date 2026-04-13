# Empirical size / power of grdd_inference() for network outcomes (same DGP as simulation_estimation.R).
# Parallel over bootstrap replications; seeds are deterministic from base_seed.
#
# Inputs:  none.
# Outputs: output/rdata/inference.RData containing `reject_array` (see figure1.R).
# Runtime: very long (bootstrap inside each MC rep); use parallel workers as in simulation_estimation.R.
#
# Paper: Figure 1(b) — run scripts/figure1.R after this (or use shipped inference.RData).
#
# Usage (from repository root):  Rscript scripts/simulation_inference.R

source(file.path("R", "paths.R"))
source(grdd_path("R", "grdd.R"))
source(grdd_path("R", "grdd_inference.R"))

library(truncnorm)
library(foreach)
library(doSNOW)
library(parallel)

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

generate_network <- function(r, delta, n_nodes = 10) {
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
        tau <- delta * tau_ij(block_i, block_j)
        weight <- base + ifelse(r < 0, 0, tau) + eps
        adj[i, j] <- adj[j, i] <- weight
      }
    }
  }

  deg <- rowSums(adj)
  diag(deg) - adj
}

run_inference_rep <- function(delta, n, alpha, B, data_seed, boot_seed) {
  set.seed(data_seed)
  r <- rtruncnorm(n, a = -1, b = 1, mean = 0, sd = 0.1)
  y <- lapply(r, function(ri) generate_network(ri, delta))
  fit <- grdd(y = y, x = r, cutoff = 0, optns = list(type = "network"))
  inf <- grdd_inference(fit, alpha = alpha, B = B, seed = boot_seed)
  inf$p_value <= alpha
}

alpha <- 0.05
delta_grid <- seq(0, 1, by = 0.1)
n_reps <- 500L
sample_sizes <- c(100L, 200L, 500L, 1000L)
B_boot <- 1000L
base_seed <- 1L

reject_array <- array(
  NA,
  dim = c(length(sample_sizes), length(delta_grid), n_reps),
  dimnames = list(
    n = as.character(sample_sizes),
    delta = as.character(delta_grid),
    rep = as.character(seq_len(n_reps))
  )
)

progress <- function(q) {
  if (q %% 10L == 0L) {
    cat(sprintf("%d replications are complete\n", q))
  }
}

n_cores <- max(1L, min(50L, parallel::detectCores() - 1L))
message("Using ", n_cores, " parallel workers.")
cl <- parallel::makeCluster(n_cores)
registerDoSNOW(cl)
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

clusterExport(
  cl,
  c(
    "P_BLOCK_SBM", "EPS_UNIF_MAX",
    "tau_ij", "generate_network", "run_inference_rep",
    "alpha", "B_boot", "base_seed", "n_reps",
    "grdd", "distance", "gtm", "bw_select", "lfr", "lfr_net", "kerFctn",
    "local_linear_weights",
    "grdd_inference", "bootstrap_grdd_inference", "embed_object", "hilbert_inner"
  ),
  envir = environment()
)

for (in_idx in seq_along(sample_sizes)) {
  n_current <- sample_sizes[in_idx]
  cat(sprintf("\n=== Sample size n = %d ===\n", n_current))

  for (id in seq_along(delta_grid)) {
    delta <- delta_grid[id]
    rej_vec <- foreach(
      b = seq_len(n_reps),
      .combine = c,
      .packages = "truncnorm",
      .options.snow = list(progress = progress)
    ) %dopar% {
      data_seed <- base_seed + (in_idx - 1L) * length(delta_grid) * n_reps + (id - 1L) * n_reps + b
      boot_seed <- data_seed + 1000000L
      tryCatch(
        run_inference_rep(
          delta = delta,
          n = n_current,
          alpha = alpha,
          B = B_boot,
          data_seed = data_seed,
          boot_seed = boot_seed
        ),
        error = function(e) {
          warning(
            "Replication failed (n = ", n_current, ", delta = ", delta, ", rep = ", b, "): ",
            conditionMessage(e),
            call. = FALSE
          )
          NA
        }
      )
    }

    reject_array[in_idx, id, ] <- rej_vec
    cat(sprintf("n = %d, delta = %.1f: done (%d reps)\n", n_current, delta, n_reps))
  }
}

out <- grdd_path("output", "rdata", "inference.RData")
save(reject_array, file = out)
message("Saved: ", out)
