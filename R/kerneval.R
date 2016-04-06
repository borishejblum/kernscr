#'Evaluation of kernels
#'
#'@param K which kernel is evaluated by \code{kerneval}.Possible values include currently implemented kernels
#'designated by a character string \code{"linear"}, \code{"poly"} and \code{"gaussian"}. Otherwise can be a user defined function.
#'
#'
#'@param Z \code{p x n} size Gram matrix.
#'
#'@param ... other arguments to be passed to the kenrel evaluated.
#'
#'@details \code{kernelEval} kernel eval works only for gaussian, polynomial and linear kernels currently.
#'To be implemented for other kernels too...
#'
#'
#'@rdname kernelEval
#'@keywords internal
#'@export
kernelEval <- function(Z, K = c("linear", "poly", "gaussian"), ... ){
  if (K  == "linear"){
    return(linKernelEval(Z))
  }
  else if (K == "poly"){
    return(polyKernelEval(Z, ...))
  }
  else if(K == "gaussian"){
    return(gaussKernelEval(Z, ...))
  }
  else{
    return(genericKernelEval(Z, K, ...))
  }
}


#'@rdname kernelEval
#'@keywords internal
linKernelEval <- function(Z){ Z ## p x n matrix
  Z <- as.matrix(Z) # this is the size Gram matrix
  m <- ncol(Z)
  ## inefficient for loop:
  # G <- matrix(data = 0, nrow = m, ncol = m); for (i in 1:m){G[i,] <- as.vector((K(Z, Z[,i], ...)))}
  G <-  t(Z)%*%Z
  return(G)
}


#'@param sigma sigma variance parameter for the \code{"gaussian"} kernel.
#'@rdname kernelEval
#'@keywords internal
gaussKernelEval <- function(Z, sigma = 1){

  Z <- as.matrix(Z) # this is the size Gram matrix

  m <- ncol(Z)
  # computing the Gram matrix
  Z.rep <- Z[,rep(1:m, rep(m, m))] - Z[,rep(1:m, m)]
  G <- matrix(exp(-colSums(Z.rep * Z.rep)/(2*sigma^2)), m, byrow = T)

  return(G)
}

#'@rdname kernelEval
#'@keywords internal
#'@export
gaussKernelEval_multipleSigmas <- function(Z, sigma = 1){
  Z <- as.matrix(Z) # this is the size Gram matrix
  m <- ncol(Z)

  # computing the Gram matrix
  Z.rep <- Z[,rep(1:m, rep(m, m))] - Z[,rep(1:m, m)]
  # remeber row are the matrix of interest
  # to obtain a particular matrix of interest for sigma[i] for e.g.
  # matrix(G[i,], length(dim))
  G <- exp(-matrix(colSums(Z.rep*Z.rep)[rep(1:ncol(Z.rep), length(sigma))],
                   nrow=length(sigma), byrow=T)/(sigma))
  return(G)
}

#'@param a TODO of the polynomial for the \code{"poly"}. Default is \code{0}.
#'@param d degree of the polynomial for the \code{"poly"}. Default is \code{2} (Quadratic kernel)
#'@rdname kernelEval
#'@keywords internal
polyKernelEval <- function(Z, a = 0, d = 1){
  Z <- as.matrix(Z) # this is the size Gram matrix
  m <- ncol(Z)
  return((matrix(as.vector(t(Z)%*%Z) + rep(a,rep(m*m)), ncol = m*m, byrow = T))^d)
}


#'@rdname kernelEval
#'@keywords internal
#'@export
polyKernelEval_multiple <- function(Z, a = 0, d = 1){
  Z <- as.matrix(Z)# this is the size Gram matrix
  m <- ncol(Z)
  return((matrix(as.vector(t(Z)%*%Z) + rep(a,rep(m*m,length(a))), ncol = m*m, byrow = T))^d)
}


#'@rdname kernelEval
#'@details \code{genericKernelEval}
genericKernelEval <- function(Z, K, ...){
  ## this function is too slow; don't know how to optimize it as of today...
  ## maybe somekind of apply to the vector columns
  if (!is.function(K)){
    cat("K must be a kernel function defined for the columns of X")
    return(0)
  }

  Z <- as.matrix(Z) # this is the size Gram matrix
  m <- ncol(Z)

  # computing the Gram matrix
  # these for loops here are not efficient
  G <- matrix(data = 0, nrow = m, ncol = m)
  for (i in 1:m){
    for (j in 1:m){
      G[i,j] <- K(as.vector(Z[,i]), as.vector(Z[,j]), ...)
    }
  }

  return(G)
}


