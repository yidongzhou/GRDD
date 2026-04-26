#' @title Local Fréchet Regression for Compositional Data
#'
#' @description
#' Performs local Fréchet regression for compositional data using the geodesic metric on the sphere.
#'
#' @param y A list of length \code{n}, where \code{y[[i]]} contains the \eqn{i}th compositional observation.
#' @param x An \eqn{n \times p} matrix or data frame of predictor values. If \eqn{p = 1}, 
#'   \code{x} can be provided as a numeric vector of length \code{n}.
#' @param xOut An optional \eqn{n_{\mathrm{out}} \times p} matrix or data frame specifying the predictor 
#'   levels at which to evaluate the fitted regression. Defaults to \code{x}.
#' @param optns A list of control parameters specified as \code{list(name = value)}. 
#'   See \strong{Details} for available options.
#'
#' @details
#' The following options may be specified in \code{optns}:
#' \describe{
#'   \item{\code{bw}}{Bandwidth for local regression (required; selected in \code{grdd()} in our workflow).}
#'   \item{\code{kernel}}{Kernel function to use. Available options include \code{"gaussian"} (default), 
#'     \code{"triangular"}, \code{"uniform"}, \code{"epanechnikov"}, \code{"gausvar"}, and \code{"quartic"}.}
#' }
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{predict}}{A list of predicted compositions at the predictor levels \code{xOut}.}
#'   \item{\code{xOut}}{The predictor levels at which predictions are evaluated.}
#'   \item{\code{y}}{The input compositional data.}
#'   \item{\code{x}}{The input predictor values.}
#'   \item{\code{optns}}{The control options used in the estimation.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' # result <- lfr_com(y, x, xOut, optns = list(bw = 0.3, kernel = "triangular"))
#' }
#'
#' @export

lfr_com <- function(y = NULL, x = NULL, xOut = NULL, optns = list()){
  if (is.null(y) || is.null(x)) {
    stop("requires the input of both y and x")
  }
  if (!is.list(y)) {
    stop("y must be a list")
  }
  if(!all(sapply(y, is.numeric))) {
    stop('all elements of y must be numeric')
  }
  if(length(unique(sapply(y, length))) != 1) {
    stop('all elements of y must have the same length')
  }
  ctg <- names(y[[1]])
  yM <- matrix(unlist(y), nrow = length(y), byrow = TRUE)
  if (!is.matrix(x)) {
    if (is.data.frame(x) | is.vector(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix or a data frame or a vector")
    }
  }
  n <- nrow(x)# number of observations
  p <- ncol(x)# number of covariates
  if (p > 2) {
    stop("local method is designed to work in low dimensional case (p is either 1 or 2)")
  }
  if (nrow(yM) != n) {
    stop("the number of rows in x must be the same as the number of rows in y")
  }
  if (!is.null(xOut)) {
    if (!is.matrix(xOut)) {
      if (is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else if (is.vector(xOut)) {
        if (p == 1) {
          xOut <- as.matrix(xOut)
        } else {
          xOut <- t(xOut)
        }
      } else {
        stop("xOut must be a matrix or a data frame or a vector")
      }
    }
    if (ncol(xOut) != p) {
      stop("x and xOut must have the same number of columns")
    }
  } else {
    xOut <- x
  }
  nOut <- nrow(xOut) # number of predictions
  if (is.null(optns$kernel)) {
    optns$kernel <- "gaussian"
  }
  if (any(abs(rowSums(yM) - 1) > 1.5e-8)){
    yM <- yM / rowSums(yM)
    warning("Each row of y has been standardized to enforce sum equal to 1")
  }
  yM <- sqrt(yM)# from simplex to sphere
  yOut <- lapply(1:nOut, function(i) {
    yOuti <- lfr_com0(x, xOut[i, ], yM, optns$bw, optns$kernel)
    yOuti <- yOuti^2# from sphere to simplex
    names(yOuti) <- ctg
    yOuti
  })
  res <- list(predict = yOut, xOut = xOut, y = y, x = x, optns = optns)
  return(res)
}

# a: output predictor level
lfr_com0 <- function(x, a, y, bw, kernel) {
  n <- nrow(x)
  w <- ll_weights(x, a, bw, kernel, ridge = TRUE)
  y0 <- apply(y, 2, weighted.mean, w)# initial guess
  y0 <- y0 / l2norm(y0)
  
  objFctn <- function(y0) {
    if (abs(l2norm(y0) - 1) > 1.5e-8) {
      return(list(value = Inf))
    }
    f <- mean(w * sapply(1:n, function(i) SpheGeoDist(y[i, ], y0)^2))
    # f = mean_i w_i d_i^2  =>  grad = (2/n) sum_i w_i d_i (nabla d_i); each column i below is w_i d_i nabla d_i
    D <- ncol(y)
    G <- vapply(seq_len(n), function(i) {
      w[i] * SpheGeoDist(y[i, ], y0) * SpheGeoGrad(y[i, ], y0)
    }, numeric(D))
    g <- 2 * rowMeans(G)
    res <- sapply(1:n, function(i) {
      grad_i <- SpheGeoGrad(y[i, ], y0)
      return((grad_i %*% t(grad_i) + SpheGeoDist(y[i, ], y0) * SpheGeoHess(y[i, ], y0)) * w[i])
    }, simplify = "array")
    h <- 2 * apply(res, 1:2, mean)
    return(list(value = f, gradient = g, hessian = h))
  }
  
  res <- tryCatch(
    trust::trust(objFctn, y0, 0.1, 1e5),
    error = function(e) list(argument = y0)
  )
  yOut0 <- res$argument / l2norm(res$argument)
  # project
  yOut <- pmax(yOut0, 0)
  if (sum(yOut^2) < 1e-15) {# all components are negative
    yOut <- numeric(length(yOut))
    yOut[which.max(yOut0)] <- 1
  } else {
    yOut <- yOut / sqrt(sum(yOut^2))
  }
  yOut
}

#'@title Geodesic distance on spheres.
#'@param y1,y2 Two unit vectors, i.e., with \eqn{L^2} norm equal to 1, of the same length.
#'@return A scalar holding the geodesic distance between \code{y1} and \code{y2}.
#'@examples
#'d <- 3
#'y1 <- rnorm(d)
#'y1 <- y1 / sqrt(sum(y1^2))
#'y2 <- rnorm(d)
#'y2 <- y2 / sqrt(sum(y2^2))
#'dist <- SpheGeoDist(y1,y2)
#'@export
SpheGeoDist <- function(y1, y2) {
  if (length(y1) != length(y2)) {
    stop("y1 and y2 should be of the same length")
  }
  y1 <- y1 / l2norm(y1)
  y2 <- y2 / l2norm(y2)
  return(acos(max(min(sum(y1 * y2), 1), -1)))
}

#' Compute gradient w.r.t. y of the geodesic distance \eqn{\arccos(y1^\top y2)} on a unit hypersphere
#' @param y1, y2 Two unit vectors.
#' @return A vector holding gradient w.r.t. \code{y2} of the geodesic distance between \code{y1} and \code{y2}.
#' @export
SpheGeoGrad <- function(y1, y2) {
  dot <- sum(y1 * y2)
  dot <- max(min(dot, 1), -1)
  # tmp = 1 - dot^2 vanishes when y1 || y2 (geodesic gradient singular); floor for trust / FP safety
  tmp <- max(1 - dot^2, 1e-10)
  return(-tmp^(-0.5) * y1)
}

#' Hessian \eqn{\partial^2/\partial y \partial y^\top} of the geodesic distance \eqn{\arccos(x^\top y)} on a unit hypersphere
#' @param x,y Two unit vectors.
#' @return A Hessian matrix.
#' @export
SpheGeoHess <- function(x, y) {
  dot <- sum(x * y)
  dot <- max(min(dot, 1), -1)
  tmp <- max(1 - dot^2, 1e-10)
  return(-dot * tmp^(-1.5) * x %*% t(x))
}

# L2 norm
l2norm <- function(x){
  #sqrt(sum(x^2))
  as.numeric(sqrt(crossprod(x)))
}