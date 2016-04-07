#'Kernel functions for kernscr
#'
#'Compute kernel TODO
#'
#'@param x TODO
#'@param y TODO
#'@param sigma variance parameter for the gaussain kernel.
#'@rdname kernels
#'@name kernels
#'@keywords internal
#'@export
GaussianKernel <- function(x, y, sigma = 1){
  return(exp(-t(x-y)%*%(x-y)/(2*sigma^2)))
}

#'@param rho TODO correlation parameter for the linear kernel.
#'Default is 0.
#'@rdname kernels
#'@name kernels
#'@keywords internal
#'@export
LinearKernel <- function(x, y, rho = 0){
  return(t(x)%*%y + rho)
}

#'@param a TODO of the polynomial for the \code{PolyKernel}. Default is \code{0}.
#'@param d degree of the polynomial for the \code{PolyKernel}. Default is \code{2} (Quadratic kernel)
#'@keywords internal
#'@rdname kernels
#'@name kernels
#'@export
PolyKernel <- function(x, y, a=0, d=2){
  return((t(x)%*%y + a)^d)
}