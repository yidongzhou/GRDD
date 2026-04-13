#' @title Geodesic Regression Discontinuity Design (GRDD)
#'
#' @description
#' Implements regression discontinuity design (RDD) for outcomes that are random objects 
#' taking values in a general metric space.
#'
#' @param y A list of length \code{n}, where \code{y[[i]]} contains the \eqn{i}th random object.
#' @param x A numeric vector of length \code{n} representing the running (forcing) variable.
#' @param cutoff A scalar specifying the cutoff value on the running variable \code{x};
#'   defaults to \code{cutoff = 0}.
#' @param optns A list of control parameters specified via \code{list(name = value)}. 
#'   See the \strong{Details} section below.
#'
#' @details
#' Available options in \code{optns} include:
#' \describe{
#'   \item{\code{type}}{Type of random object. Supported types include: 
#'   \code{"composition"}, \code{"euclidean"}, \code{"function"}, \code{"measure"}, and 
#'   \code{"network"}, \code{"spd"}.}
#'   \item{\code{fuzzy}}{Optional binary vector of length \code{n} indicating treatment 
#'   assignment for a fuzzy RDD. If omitted, the default is a sharp design.}
#'   \item{\code{bw}}{Bandwidth parameter used for local estimation. If not provided, 
#'   it is selected using a mean squared error (MSE)-optimal bandwidth selector.}
#'   \item{\code{bw_optns}}{Optional list of control parameters passed to
#'   \code{bw_select} when selecting \code{bw}. Supported fields include
#'   \code{n_bw}, \code{min_obs}, \code{max_eval}, \code{delta}, \code{delta_min},
#'   \code{delta_max}, \code{bw_min}, and \code{bw_max}. If \code{delta} is
#'   provided, \code{delta_min} and \code{delta_max} are ignored.}
#'   \item{\code{kernel}}{Kernel function for local estimation. Choices include 
#'   \code{"triangular"} (default), \code{"epanechnikov"}, and \code{"uniform"}.}
#'   \item{\code{supp}}{For functional data, a numeric vector indicating the domain (support)
#'   over which the functions are defined.}
#'   \item{\code{lower}}{For distributional data, the lower bound of the support of the measure. Default is \code{NULL}.}
#'   \item{\code{upper}}{For distributional data, the upper bound of the support of the measure. Default is \code{NULL}.}
#' }
#'
#' @return
#' An object of class \code{"grdd"}, which is a list containing local average treatment 
#' effect estimates on the metric space.

# Source required functions (run R from the GRDD repository root, or set GRDD_ROOT).
grdd_repo_root <- Sys.getenv("GRDD_ROOT", "")
if (!nzchar(grdd_repo_root)) grdd_repo_root <- getwd()
source(file.path(grdd_repo_root, "R", "paths.R"))
source(grdd_path("R", "kerFctn.R"))
source(grdd_path("R", "local_linear.R"))
source(grdd_path("R", "lfr_com.R"))
source(grdd_path("R", "lfr_fun.R"))
source(grdd_path("R", "lfr_mea.R"))
source(grdd_path("R", "lfr_net.R"))
source(grdd_path("R", "lfr_spd.R"))
source(grdd_path("R", "lfr_euc.R"))
source(grdd_path("R", "lcm.R"))

grdd <- function(y, x, cutoff = 0, optns = list()) {
  # Input validation
  if (is.null(y)) {
    stop("requires the input of y")
  }
  if (!is.list(y)) {
    stop("y must be a list")
  }
  if (!is.numeric(x)) {
    stop("x must be a numeric vector")
  }
  if (length(x) != length(y)) {
    stop("the length of x must be the same as the length of y")
  }
  
  # Validate cutoff value
  if (cutoff < min(x) || cutoff > max(x)) {
    stop("Cutoff must be within the range of x")
  }
  
  # Check for enough observations
  if (sum(x < cutoff) < 2 || sum(x >= cutoff) < 2) {
    stop("Not enough observations on both sides of the cutoff")
  }
  
  # Set default type if not specified
  if (is.null(optns$type)) {
    optns$type <- 'euclidean'
  }
  
  # Validate function type specific requirements
  if (optns$type == 'function') {
    if (is.null(optns$supp)) {
      stop("for functional data, the support must be specified in optns$supp")
    }
  }
  n <- length(y)
  if (optns$type == "measure") {
    y <- harmonize_measure(y)
  }
  if (optns$type == "composition") {
    # Standardize compositions up front so estimation and inference use the same inputs.
    tol <- 1.5e-8
    sums <- vapply(y, function(yi) sum(as.numeric(yi)), numeric(1))
    if (any(!is.finite(sums)) || any(sums <= 0)) {
      stop("For type = \"composition\", each y[[i]] must have a positive finite sum.", call. = FALSE)
    }
    needs_norm <- any(abs(sums - 1) > tol)
    if (needs_norm) {
      y <- lapply(y, function(yi) {
        s <- sum(as.numeric(yi))
        yi / s
      })
      warning("Each composition in y has been standardized to enforce sum equal to 1", call. = FALSE)
    }
  }
  
  # Set default options
  if (is.null(optns[["fuzzy"]])) {
    optns[["fuzzy"]] <- FALSE
  }
  if (is.null(optns[["kernel"]])) {
    optns[["kernel"]] <- "triangular"
  }
  # Use exact matching for option names (avoid $ partial matching, e.g. bw vs bw_optns)
  if (is.null(optns[["bw"]])) {
    bw_ctrl <- if (is.null(optns[["bw_optns"]])) list() else optns[["bw_optns"]]
    if (!is.list(bw_ctrl)) {
      stop("optns$bw_optns must be a list (e.g., list(max_eval = 100L, bw_max = 2))")
    }
    # Avoid accidental positional matching / typos by only passing known args
    bw_ctrl <- bw_ctrl[names(bw_ctrl) %in% names(formals(bw_select))]
    bw_args <- c(list(y = y, x = x, cutoff = cutoff, optns = optns), bw_ctrl)
    optns[["bw"]] <- do.call(bw_select, bw_args)
  }
  # Perform regression
  res_left <- lfr(y[x < cutoff], x[x < cutoff], xOut = cutoff, optns = optns)
  res_right <- lfr(y[x >= cutoff], x[x >= cutoff], xOut = cutoff, optns = optns)
  
  # Return results
  res <- list(tau = list(left = res_left$predict[[1]], 
                         right = res_right$predict[[1]]), 
              y = y, 
              x = x, 
              cutoff = cutoff,
              optns = optns)
  class(res) <- c("grdd", "list")
  res
}

lfr <- function(y, x, xOut, optns) {
  switch(
    optns$type,
    "composition" = lfr_com,
    "function" = lfr_fun,
    "measure" = lfr_mea,
    "network" = lfr_net,
    "spd" = lfr_spd,
    "euclidean" = lfr_euc,
    stop("Unsupported type in `optns$type`.")
  )(y, x, xOut, optns)
}

#' Bandwidth selection for GRDD using boundary-focused Fréchet cross-validation
#' (in the spirit of Ludwig & Miller (2007) and Imbens & Kalyanaraman (2012, §4.5))
#'
#' Chooses \code{delta} in \eqn{[\code{delta_min},\,\code{delta_max}]} so the
#' evaluation count is as close as possible to \code{max_eval} (default \code{100}):
#' use \code{delta_max} if even at \code{delta_max} there are at most
#' \code{max_eval} points; use \code{delta_min} if at \code{delta_min} there are
#' already at least \code{max_eval} points; otherwise bisect and pick the endpoint
#' nearer to \code{max_eval}. Defaults \code{delta_min = 0.05},
#' \code{delta_max = 0.5}.
bw_select <- function(y, x, cutoff, optns, 
                      n_bw = 50,
                      min_obs = 10,
                      max_eval = 100L,
                      delta = NULL,
                      delta_min = 0.05,
                      delta_max = 0.5,
                      bw_min = NULL,
                      bw_max = NULL) {
  
  if (!is.null(delta)) {
    if (!(is.numeric(delta) && length(delta) == 1 && delta > 0 && delta <= 1)) {
      stop("delta must be a scalar in (0, 1]")
    }
  } else {
    if (!(delta_min < delta_max && delta_min > 0 && delta_max <= 1)) {
      stop("need 0 < delta_min < delta_max <= 1")
    }
  }
  
  # Sort by running variable
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  
  # Identify observations on each side
  left_idx  <- which(x < cutoff)
  right_idx <- which(x >= cutoff)
  N_left    <- length(left_idx)
  N_right   <- length(right_idx)
  
  if (N_left < 10 || N_right < 10) {
    stop("Not enough observations on both sides of the cutoff for CV")
  }
  # k-th smallest |x - cutoff| on each side (same k). Cap k by min_obs, sample
  # size, and ~ floor(N/3) per side so the order statistic stays moderate vs.
  # default bw_max (half the shorter x-span), avoiding bw_min > bw_max.
  half_left  <- max(3L, N_left %/% 3L)
  half_right <- max(3L, N_right %/% 3L)
  min_obs_eff <- min(as.integer(min_obs), half_left, half_right)
  
  eval_count <- function(del) {
    tl <- quantile(x[left_idx],  probs = 1 - del, type = 1)
    tr <- quantile(x[right_idx], probs = del,     type = 1)
    length(which(x >= tl & x <= tr))
  }
  
  K <- as.integer(max_eval)
  if (K < 10L) {
    warning("max_eval < 10; using max_eval = 10 for CV stability")
    K <- 10L
  }
  if (is.null(delta)) {
    n_lo <- eval_count(delta_min)
    n_hi <- eval_count(delta_max)
    if (n_hi <= K) {
      delta <- delta_max
    } else if (n_lo >= K) {
      delta <- delta_min
    } else {
      # Largest delta with count <= K
      lo <- delta_min
      hi <- delta_max
      for (iter in seq_len(50L)) {
        if (hi - lo < 1e-6) break
        mid <- (lo + hi) / 2
        if (eval_count(mid) <= K) lo <- mid else hi <- mid
      }
      d_le <- lo
      n_le <- eval_count(d_le)
      # Smallest delta with count >= K
      lo <- delta_min
      hi <- delta_max
      for (iter in seq_len(50L)) {
        if (hi - lo < 1e-6) break
        mid <- (lo + hi) / 2
        if (eval_count(mid) >= K) hi <- mid else lo <- mid
      }
      d_ge <- hi
      n_ge <- eval_count(d_ge)
      delta <- if (abs(n_le - K) <= abs(n_ge - K)) d_le else d_ge
    }
  }
  theta_left  <- quantile(x[left_idx], probs = 1 - delta, type = 1)
  theta_right <- quantile(x[right_idx], probs = delta, type = 1)
  eval_idx <- which(x >= theta_left & x <= theta_right)
  
  # Step 2: Bandwidth grid
  if (is.null(bw_min)) {
    bw_min <- max(max(diff(x)), 
                  sort(abs(x[x < cutoff] - cutoff))[min_obs_eff] + 1e-8, 
                  sort(abs(x[x >= cutoff] - cutoff))[min_obs_eff] + 1e-8)
  }
  if (is.null(bw_max)) {
    bw_max <- min(max(x) - cutoff, cutoff - min(x)) / 2
  }
  if (!(is.numeric(bw_min) && length(bw_min) == 1 && bw_min > 0)) {
    stop("bw_min must be a positive scalar")
  }
  if (!(is.numeric(bw_max) && length(bw_max) == 1 && bw_max > 0)) {
    stop("bw_max must be a positive scalar")
  }
  if (bw_min >= bw_max) {
    bw_floor <- max(max(diff(x)) + 1e-8, bw_max * 0.01)
    bw_min <- min(bw_floor, bw_max * 0.99)
  }
  if (bw_min >= bw_max) {
    stop("could not enforce bw_min < bw_max; try smaller min_obs or set bw_min/bw_max in bw_optns")
  }
  bw_grid <- seq(bw_min, bw_max, length.out = n_bw)
  
  # optns with each candidate bandwidth
  optns_by_bw <- lapply(bw_grid, function(h) utils::modifyList(optns, list(bw = h)))
  
  # Step 3: CV score = sum of squared errors over the boundary window
  # Outer loop over evaluation points (training data do not depend on h);
  # inner loop over bandwidths.
  cv_vals <- numeric(length(bw_grid))
  
  for (k in seq_along(eval_idx)) {
    idx <- eval_idx[k]
    xj  <- x[idx]
    
    # Strictly one-sided training set (key feature from IK 2012 §4.5)
    if (xj < cutoff) {
      train_idx <- which(x < xj)
    } else {
      train_idx <- which(x > xj)
    }
    
    if (length(train_idx) < 3) {
      next
    }
    
    y_train <- y[train_idx]
    x_train <- x[train_idx]
    
    for (i in seq_along(bw_grid)) {
      pred_obj <- lfr(y = y_train, x = x_train, xOut = xj, optns = optns_by_bw[[i]])
      pred <- pred_obj$predict[[1]]
      cv_vals[i] <- cv_vals[i] + distance(y[[idx]], pred, optns)^2
    }
  }
  
  # Step 4: Select optimal bandwidth
  optimal_bw <- bw_grid[which.min(cv_vals)]
  
  # Return single bandwidth (used on both sides)
  optimal_bw
}

distance <- function(x, y, optns) {
  if (optns$type == 'composition') {
    acos(max(min(sum(sqrt(x) * sqrt(y)), 1), -1))
  } else if (optns$type == 'function') {
    sqrt(pracma::trapz(x = optns$supp, y = (x - y)^2))
  } else if (optns$type == 'measure') {
    # sqrt(mean((sort(x) - sort(y)) ^ 2))
    sqrt(pracma::trapz(x = seq(0, 1, by = 1 / length(x))[-1], y = (x - y)^2))
  } else {# euclidean or network or spd
    sqrt(sum((x - y) ^ 2))
  }
}

# geodesic transport
gtm <- function(alpha, beta, omega, optns) {
  if (optns$type == 'measure') {
    if(alpha[1] > omega[1]) {
      omega[1] <- alpha[1]
    }
    if(alpha[length(alpha)] > omega[length(omega)]) {
      omega[length(omega)] <- alpha[length(omega)]
    }
    zeta <- pracma::spinterp(x = alpha, y = beta, xp = omega)
    zeta[1] <- optns$lower
    zeta[length(zeta)] <- optns$upper
  } else if (optns$type %in% c('composition', 'mvmeasure')) {
    # log--parallel--exp on the unit sphere (both use norm-standardized vectors)
    ab <- sum(alpha * beta)
    ao <- sum(alpha * omega)
    if (1 + ao < 1e-10) stop("The minimizing geodesic from alpha to omega is not unique (omega ≈ -alpha).")
    ab <- max(min(ab, 1), -1)
    theta <- acos(ab)
    s <- sqrt(max(1 - ab^2, 0))
    if (theta < 1e-12 || s < 1e-12) {
      zeta0 <- omega
    } else {
      u <- (beta - ab * alpha) / s
      v <- theta * u
      vo <- sum(v * omega)
      ptv <- v - (vo / (1 + ao)) * (alpha + omega)
      nv <- sqrt(sum(ptv^2))
      if (nv < 1e-12) {
        zeta0 <- omega
      } else {
        zeta0 <- cos(nv) * omega + sin(nv) * (ptv / nv)
      }
    }
    zeta <- pmax(zeta0, 0)
    if (sum(zeta^2) < 1e-15) {
      zeta <- numeric(length(zeta))
      zeta[which.max(zeta0)] <- 1
    } else {
      zeta <- zeta / sqrt(sum(zeta^2))
    }
  } else if (optns$type %in% c('function', 'network', 'spd')) {
    zeta <- omega + beta - alpha
  }
  zeta
}

harmonize_measure <- function(y) {
  n <- length(y)
  N <- sapply(y, length)
  if (length(unique(N)) == 1L) {
    return(y)
  }
  M <- min(plcm(N), n * max(N), 5000)
  lapply(seq_len(n), function(i) {
    residual <- M %% N[i]
    if (residual) {
      sort(c(rep(y[[i]], each = M %/% N[i]), sample(y[[i]], residual)))
    } else {
      sort(rep(y[[i]], each = M %/% N[i]))
    }
  })
}