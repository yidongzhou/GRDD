# Local linear regression weights (product kernel).
# Source after kerFctn.R (see scripts that source R/*.R in order). Used by lfr_* and grdd_inference.

#' Local linear weights at one evaluation point
#'
#' Predictors use a product kernel \eqn{\prod_j K((x_{ij}-a_j)/h_j)}. Returns
#' unnormalized weights \eqn{w_i} so that a weighted mean of \eqn{y_i} with
#' weights \eqn{w_i} gives the local linear fit at \eqn{a}.
#'
#' @param x \eqn{n \times p} matrix of predictors, or a numeric vector of length
#'   \eqn{n} when \eqn{p = 1} (e.g. running variable in GRDD).
#' @param a Evaluation point, length \eqn{p} (\code{length(a)} must match the
#'   number of predictors after \code{x} is coerced to a matrix).
#' @param bw Bandwidth vector, length \eqn{p}.
#' @param kernel Kernel name for \code{kerFctn()}.
#' @param ridge If \code{TRUE}, add a small ridge to the kernel-weighted second
#'   moment of \eqn{x - a} so the local-linear system is never singular.
#'
#' @return Numeric vector of length \eqn{n} (unnormalized weights).
ll_weights <- function(x, a, bw, kernel, ridge = TRUE) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  } else if (!is.matrix(x)) {
    x <- cbind(as.numeric(x))
  } else {
    x <- as.matrix(x)
  }
  n <- nrow(x)
  p <- ncol(x)
  a <- as.numeric(a)
  bw <- as.numeric(bw)
  if (length(a) != p) {
    stop("length(a) must equal ncol(x).", call. = FALSE)
  }
  if (length(bw) != p) {
    stop("length(bw) must equal ncol(x).", call. = FALSE)
  }
  if (n < 1L) {
    stop("x must have at least one row.", call. = FALSE)
  }

  Kern <- kerFctn(kernel)
  Kprod <- function(d) {
    v <- 1
    for (j in seq_len(p)) {
      v <- v * Kern(d[j] / bw[j])
    }
    v
  }

  zm <- sweep(x, 2L, a, "-")
  ki <- apply(zm, 1L, Kprod)
  mu1 <- colMeans(sweep(zm, 1L, ki, `*`))
  mu2 <- matrix(0, p, p)
  for (i in seq_len(n)) {
    mu2 <- mu2 + ki[i] * tcrossprod(zm[i, , drop = FALSE])
  }
  mu2 <- mu2 / n

  if (ridge) {
    tr <- sum(diag(mu2))
    eps_reg <- if (is.finite(tr) && tr > 0) 1e-12 * tr / p else 1e-12
    mu2 <- mu2 + eps_reg * diag(p)
  }
  mu1v <- matrix(mu1, ncol = 1L)
  wc <- t(solve(mu2, mu1v))
  ki * (1 - c(zm %*% t(wc)))
}
