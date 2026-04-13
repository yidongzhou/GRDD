#' @title Local Fréchet Regression for Functional Data
#'
#' @description
#' Performs local Fréchet regression for functional data with respect to the \eqn{L^2} metric.
#'
#' @param y A list of length \code{n}, where \code{y[[i]]} contains the measurements 
#'   for the \eqn{i}th functional observation.
#' @param x An \eqn{n \times p} matrix or data frame of predictor values. If \eqn{p = 1}, 
#'   \code{x} can be a numeric vector of length \code{n}.
#' @param xOut An optional \eqn{n_{\mathrm{out}} \times p} matrix or data frame of predictor levels 
#'   at which predictions are evaluated. Defaults to \code{x}.
#' @param optns A list of control parameters specified as \code{list(name = value)}.
#'   See \strong{Details}.
#'
#' @details
#' Available options in \code{optns} include:
#' \describe{
#'   \item{\code{bw}}{Bandwidth for local regression. If omitted, it is selected via cross-validation.}
#'   \item{\code{kernel}}{Kernel function to use. Available options are \code{"gaussian"} (default), 
#'     \code{"triangular"}, \code{"uniform"}, \code{"epanechnikov"}, \code{"gausvar"}, and \code{"quartic"}.}
#' }
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{predict}}{A list of predicted functional values at \code{xOut}.}
#'   \item{\code{xOut}}{The predictor levels at which predictions are made.}
#'   \item{\code{y}}{The input functional responses.}
#'   \item{\code{x}}{The predictor variables used.}
#'   \item{\code{optns}}{The control options used.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example:
#' # result <- lfr_fun(y, x, xOut, optns = list(bw = 0.5, kernel = "gaussian"))
#' }
#'
#' @export

lfr_fun <- function(y = NULL, x = NULL, xOut = NULL, optns = list()){
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
  y <- matrix(unlist(y), nrow = length(y), byrow = TRUE)
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
  if (nrow(y) != n) {
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
  # if (is.null(optns$bw)) {
  #   optns$bw <- bwCV_fun(y = y, x = x, xOut = xOut, optns = optns)
  # } else {
  #   if (length(optns$bw) != p) {
  #     stop("dimension of bandwidth does not agree with x")
  #   }
  #   # Check if any bandwidth is too small
  #   isSmallBw <- logical(p)
  #   for (j in 1:p) {
  #     isSmallBw[j] <- optns$bw[j] < max(diff(sort(x[, j])))
  #   }
  #   # Reset bandwidth if needed using cross validation
  #   if (all(isSmallBw) && optns$kernel %in% c("rect", "quar", "epan")) {
  #     warning("optns$bw was set too small and is reset to be chosen by cross-validation")
  #     optns$bw <- bwCV_fun(y = y, x = x, xOut = xOut, optns = optns)
  #   }
  # }
  
  yOut <- lapply(seq_len(nOut), function(i) {
    w <- local_linear_weights(x, xOut[i, ], optns$bw, optns$kernel)
    apply(y, 2, weighted.mean, w)
  })
  res <- list(predict = yOut, xOut = xOut, y = y, x = x, optns = optns)
  return(res)
}

# # bandwidth selection via cross validation
# bwCV_fun <- function(y, x, xOut, optns) {
#   n <- nrow(x)
#   p <- ncol(x)
#   
#   # k-fold
#   objFctn <- function(h) {
#     numFolds <- ifelse(n > 30, 10, n)# leave-one-out or 10-fold cross-validation
#     folds <- sample(c(rep.int(1:numFolds, n%/%numFolds), seq_len(n%%numFolds)))
#     
#     cv <- 0
#     for (foldidx in seq_len(numFolds)) {
#       # nn by M
#       testidx <- which(folds == foldidx)
#       for (j in testidx) {
#         yOutj <- lfr_fun0(x[-testidx, , drop = FALSE], x[j, ], y[-testidx, , drop = FALSE], h, optns$kernel)
#         cv <- cv + sqrt(sum((y[j, ] - yOutj)^2))
#       }
#     }
#     cv
#   }
#   
#   if (p == 1) {
#     aux <- setBwRange(x = x[, 1], xOut = xOut[, 1], kernel_type = optns$kernel)
#     bwRange <- matrix(c(aux$min, aux$max), nrow = 2, ncol = 1)
#   } else {
#     aux <- setBwRange(x = x[, 1], xOut = xOut[, 1], kernel_type = optns$kernel)
#     aux2 <- setBwRange(x = x[, 2], xOut = xOut[, 2], kernel_type = optns$kernel)
#     bwRange <- as.matrix(cbind(c(aux$min, aux$max), c(aux2$min, aux2$max)))
#   }
#   
#   if (p == 1) {
#     res <- optimize(f = objFctn, interval = bwRange[, 1])$minimum
#   } else {
#     res <- optim(par = colMeans(bwRange), fn = objFctn, lower = bwRange[1, ], 
#                  upper = bwRange[2, ], method = "L-BFGS-B")$par
#   }
#   res
# }
