# =============================================================================
# Bootstrap inference for GRDD (Hilbert embedding, extrinsic inference)
# =============================================================================
# Source R/grdd.R first so kerFctn, lfr, and distance are available.
#
# Pass the object returned by grdd() as \code{fit}. Estimation settings
# (type, bw, kernel, supp, lower, upper, ...) are taken from \code{fit$optns};
# cutoff estimates from \code{fit$tau}. Pass \code{alpha}, \code{B}, and
# optionally \code{seed} for the bootstrap.
#
# Supported \code{fit$optns$type}: composition, euclidean, function, measure, network, spd.
# Bootstrap: H0: d = 0 vs H1: d > 0.
# =============================================================================

#' Embedding \eqn{y \mapsto \Psi(y)}.
#' For compositional data we use the Riemannian logarithm map at a 
#' fixed base point. The base point is chosen automatically as the sample 
#' Fréchet mean on the sphere (via \code{manifold::frechetMean}) when a list 
#' of points is supplied. For a single point the base is taken to be the point 
#' itself (yielding the zero tangent vector).
#' 
#' @param obj A list of outcomes or one outcome (e.g. \code{nu0} from \code{fit$tau}).
#' @param type Outcome type string (same as \code{optns$type} from \code{grdd()}).
#' @param base Optional base point for the log map when \code{type == "composition"}.
#'   If \code{NULL} and \code{obj} is a list, the base is set to the sample Fréchet mean
#'   of the square-root transformed compositions. If \code{NULL} and \code{obj} is a
#'   single composition, the base defaults to the point itself.
#' @return For list input, the list of embedded tangent vectors; 
#'   for a single input, the embedded tangent vector (in \(\mathbb{R}^{d-1}\)).
logmap_sphere <- function(z, p, tol = 1e-10) {
  if (!is.numeric(z) || !is.numeric(p) || length(z) != length(p)) {
    stop("z and p must be numeric vectors of the same length", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be a positive finite scalar", call. = FALSE)
  }

  # Numerical safeguards
  dot <- sum(p * z)
  dot <- pmin(pmax(dot, -1), 1)
  theta <- acos(dot)

  if (theta < tol) {
    return(rep(0, length(p)))  # zero tangent vector
  }

  # Standard formula: (theta / sin(theta)) * (z - <p,z> p)
  v <- z - dot * p
  norm_v <- sqrt(sum(v^2))
  v <- v * (theta / norm_v)

  # v is orthogonal to p and ||v|| = theta = d_g(p, z)
  v
}

embed_object <- function(obj, type, base = NULL) {
  if (type %in% c("euclidean", "function", "measure", "network", "spd")) {
    return(obj)
  }
  if (type == "composition") {
    # Inputs should be valid simplex compositions: finite, nonnegative, sum to 1.
    comp_ok <- function(comp, tol = 1e-8) {
      if (!is.numeric(comp) || any(!is.finite(comp))) return(FALSE)
      if (any(comp < -tol)) return(FALSE)
      if (abs(sum(comp) - 1) > tol) return(FALSE)
      TRUE
    }
    if (is.list(obj)) {
      # obj is a list of simplex compositions; map to the sphere via square-root transform
      if (!all(vapply(obj, comp_ok, logical(1)))) {
        stop("For type = \"composition\", each element of obj must be a valid composition (nonnegative, finite, sum to 1).",
             call. = FALSE)
      }
      obj_sphere <- lapply(obj, sqrt)
      if (is.null(base)) {
        base <- c(manifold::frechetMean(
          mfd  = structure(1, class = 'Sphere'),
          X    = matrix(unlist(obj_sphere), ncol = length(obj_sphere)),
          maxit = 1e04
        ))
      }
      emb <- lapply(obj_sphere, function(y) logmap_sphere(y, base))
      attr(emb, "base") <- base
      return(emb)
    } else {
      # single point: use itself as base (yields the zero tangent vector)
      if (!comp_ok(obj)) {
        stop("For type = \"composition\", obj must be a valid composition (nonnegative, finite, sum to 1).",
             call. = FALSE)
      }
      obj_sphere <- sqrt(obj)
      if (is.null(base)) {
        base <- obj_sphere
      }
      return(logmap_sphere(obj_sphere, base))
    }
  }
  
  stop("Unsupported or unimplemented type in embed_object().", call. = FALSE)
}

#' Frobenius / L2 inner product in the embedding space.
#' @param x,y Numeric vectors (same length) in embedded coordinates.
#' @param optns Options list from \code{grdd()} (uses \code{type}, and \code{supp} for \code{function}).
hilbert_inner <- function(x, y, optns) {
  if (optns$type %in% c("composition", "euclidean", "network", "spd")) {
    return(sum(x * y))
  }
  if (optns$type == "function") {
    if (!requireNamespace("pracma", quietly = TRUE)) {
      stop("Package \"pracma\" is required for functional data inference.", call. = FALSE)
    }
    return(pracma::trapz(optns$supp, x * y))
  }
  if (optns$type == "measure") {
    u <- seq(0, 1, by = 1 / length(x))[-1L]
    if (!requireNamespace("pracma", quietly = TRUE)) {
      stop("Package \"pracma\" is required for distributional inference.", call. = FALSE)
    }
    return(pracma::trapz(u, x * y))
  }
  stop("Unsupported type in hilbert_inner().", call. = FALSE)
}

#' Draw bootstrap replicates of \eqn{T_b} and \eqn{L_b}.
#' @param Psi \eqn{n \times d} embedded data matrix.
#' @param s0,s1 Length-\eqn{n} normalized local-linear side weights:
#'   \eqn{s_{t,i} = n_t w_i / \sum_j w_j} on side \eqn{t}, and 0 off-side.
bootstrap_grdd_inference <- function(Psi, s0, s1, Delta_hat, optns, B, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n <- nrow(Psi)
  n0 <- sum(s0)
  n1 <- sum(s1)
  bw <- optns$bw
  nh <- n * bw
  # Equivalent to weighted.mean(Psi[, j], w = w_side) with normalized s-side notation.
  hat_nu0 <- colSums(s0 * Psi) / n0
  hat_nu1 <- colSums(s1 * Psi) / n1
  T_boot <- numeric(B)
  L_boot <- numeric(B)
  for (b in seq_len(B)) {
    xi <- stats::rnorm(n)
    g0 <- sqrt(nh) / n0 * (as.vector(t(Psi) %*% (xi * s0)) - hat_nu0 * sum(xi * s0))
    g1 <- sqrt(nh) / n1 * (as.vector(t(Psi) %*% (xi * s1)) - hat_nu1 * sum(xi * s1))
    gdiff <- g1 - g0
    T_boot[b] <- hilbert_inner(gdiff, gdiff, optns)
    L_boot[b] <- 2 * abs(hilbert_inner(Delta_hat, gdiff, optns))
  }
  list(T_boot = T_boot, L_boot = L_boot, hat_nu0_Psi = hat_nu0, hat_nu1_Psi = hat_nu1)
}

#' Bootstrap inference for the GRDD treatment effect magnitude at the cutoff
#'
#' Takes the fitted object from \code{grdd()}. Bandwidth, kernel, and
#' \code{type} come from \code{fit$optns}; intrinsic estimators at the cutoff
#' from \code{fit$tau$left} and \code{fit$tau$right}.
#'
#' @param fit Object of class \code{"grdd"} from \code{grdd()} (\code{tau}, \code{y}, \code{x}, \code{cutoff}, \code{optns}).
#' @param alpha Nominal level for the one-sided test and CI.
#' @param B Number of bootstrap replicates.
#' @param seed Optional; passed to \code{set.seed()} before the bootstrap loop.
#'
#' @return A list with effect magnitude, test, p-value, CI, and \code{fit} (input).
grdd_inference <- function(fit, alpha = 0.05, B = 1000L, seed = NULL) {
  if (!inherits(fit, "grdd")) {
    stop('fit must be an object of class "grdd" (output of grdd()).', call. = FALSE)
  }

  optns <- fit$optns
  if (!optns$type %in% c("composition", "euclidean", "function", "measure", "network", "spd")) {
    stop(
      "Inference supports fit$optns$type in {composition, euclidean, function, measure, network, spd}.",
      call. = FALSE
    )
  }

  B <- as.integer(B)

  if (B < 200L) {
    warning("B < 200: bootstrap p-values and CI quantiles may be unstable.", call. = FALSE)
  }

  y <- fit$y
  x <- as.numeric(fit$x)
  cutoff <- fit$cutoff
  n <- length(x)
  if (n != length(y)) {
    stop("length(fit$x) must equal length(fit$y).", call. = FALSE)
  }
  if (cutoff < min(x) || cutoff > max(x)) {
    stop("fit$cutoff must lie in the range of fit$x.", call. = FALSE)
  }
  if (sum(x < cutoff) < 3L || sum(x >= cutoff) < 3L) {
    stop("Need at least 3 observations on each side of the cutoff.", call. = FALSE)
  }

  bw <- optns$bw
  if (length(bw) == 0L || any(!is.finite(bw))) {
    stop("fit$optns$bw must be set (run grdd() so bandwidth is selected or supplied).", call. = FALSE)
  }

  kernel <- optns$kernel
  if (is.null(kernel)) {
    kernel <- "triangular"
  }

  if (optns$type == "function" && is.null(optns$supp)) {
    stop("fit$optns$supp is required when fit$optns$type is \"function\".", call. = FALSE)
  }

  idx0 <- which(x < cutoff)
  idx1 <- which(x >= cutoff)
  n0 <- length(idx0)
  n1 <- length(idx1)

  nu0 <- fit$tau$left
  nu1 <- fit$tau$right

  if (optns$type == "measure") {
    y <- harmonize_measure(y)
  }
  y <- embed_object(y, optns$type)
  base <- NULL
  if (optns$type == "composition") {
    base <- attr(y, "base", exact = TRUE)
  }
  Psi <- do.call(rbind, lapply(y, function(yi) as.numeric(yi)))
  nu0_psi <- embed_object(nu0, optns$type, base = base)
  nu1_psi <- embed_object(nu1, optns$type, base = base)

  Delta_hat <- nu1_psi - nu0_psi
  d_hat_sq <- hilbert_inner(Delta_hat, Delta_hat, optns)
  d_hat <- sqrt(d_hat_sq)

  # Normalize side-specific local linear weights before bootstrap.
  w_left <- local_linear_weights(x[idx0], cutoff, bw, kernel)
  sw_left <- sum(w_left)
  if (!is.finite(sw_left) || sw_left <= 0) {
    stop("Local linear weights sum to a non-positive value on the left side.", call. = FALSE)
  }
  s0 <- rep(0, n)
  s0[idx0] <- n0 * w_left / sw_left

  w_right <- local_linear_weights(x[idx1], cutoff, bw, kernel)
  sw_right <- sum(w_right)
  if (!is.finite(sw_right) || sw_right <= 0) {
    stop("Local linear weights sum to a non-positive value on the right side.", call. = FALSE)
  }
  s1 <- rep(0, n)
  s1[idx1] <- n1 * w_right / sw_right

  boot <- bootstrap_grdd_inference(Psi, s0, s1, Delta_hat, optns, B, seed = seed)
  T_boot <- boot$T_boot
  L_boot <- boot$L_boot

  nh <- n * bw
  T_n <- nh * d_hat_sq
  p_value <- (1 + sum(T_boot >= T_n)) / (B + 1)
  c_alpha <- stats::quantile(L_boot, probs = 1 - alpha, names = FALSE, type = 1)
  ci_lower <- sqrt(max(d_hat_sq - c_alpha / sqrt(nh), 0))
  ci_upper <- sqrt(d_hat_sq + c_alpha / sqrt(nh))

  out <- list(
    fit = fit,
    effect_magnitude = d_hat,
    test_statistic = T_n,
    p_value = p_value,
    ci = c(lower = ci_lower, upper = ci_upper),
    alpha = alpha,
    B = B,
    seed = seed
  )
  class(out) <- c("grdd_inference", "list")
  out
}

print.grdd_inference <- function(x, digits = 4, ...) {
  reject <- x$p_value <= x$alpha
  fit <- x$fit
  type <- fit$optns$type
  n <- length(fit$x)
  n_left <- sum(fit$x < fit$cutoff)
  n_right <- sum(fit$x >= fit$cutoff)
  bw <- fit$optns$bw
  if (length(bw) > 1L) {
    bw <- bw[1L]
  }
  cat("GRDD inference (Hilbert embedding bootstrap)\n")
  cat("  type (from grdd): ", type, "\n", sep = "")
  cat("  n = ", n, ", n_left = ", n_left, ", n_right = ", n_right,
      ", h = ", format(bw, digits = digits), "\n", sep = "")
  cat("  effect magnitude d_hat: ", format(x$effect_magnitude, digits = digits), "\n", sep = "")
  cat("  test statistic T_n:    ", format(x$test_statistic, digits = digits), "\n", sep = "")
  cat("  p-value (one-sided):  ", format(x$p_value, digits = digits), "\n", sep = "")
  cat("  ", 100 * (1 - x$alpha), "% CI for d: [",
      format(x$ci[1], digits = digits), ", ",
      format(x$ci[2], digits = digits), "]\n", sep = "")
  cat("  reject H0 (alpha = ", x$alpha, "): ", reject, "\n", sep = "")
  invisible(x)
}
