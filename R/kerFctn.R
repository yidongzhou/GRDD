kerFctn <- function(kernel_type){
  if (kernel_type == 'triangular') {
    ker <- function(x) {
      (1 - abs(x)) * (abs(x) <= 1)
    }
  } else if (kernel_type == 'epanechnikov') {
    ker <- function(x) {
      (3/4) * (1 - x^2) * (abs(x) <= 1)
    }
  } else if (kernel_type == 'uniform') {
    ker <- function(x) {
      (1/2) * (abs(x) <= 1)
    }
  } else if (kernel_type == 'gaussian'){
    ker <- function(x){
      dnorm(x) #exp(-x^2 / 2) / sqrt(2*pi)
    }
  } else if(kernel_type == 'gausvar'){
    ker <- function(x) {
      dnorm(x)*(1.25-0.25*x^2)
    }
  } else if(kernel_type=='quartic'){
    ker <- function(x) {
      (15/16)*(1-x^2)^2 * (abs(x)<=1)
    }
  } else {
    stop('Unavailable kernel')
  }
  return(ker)
}