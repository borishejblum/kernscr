#'@param X  a matrix
#'@keywords internal
#'@aliases generate_dataset_advanced
#'@rdname generate_dataset_advanced
#'@name generate_dataset_advanced
#'@export
feat.sum.1.5 <- function(X){
  return(0.2*(X[,1] + X[,2] + X[,3] + X[,4] + X[,5]))
}

#'@param X  a matrix
#'@rdname generate_dataset_advanced
#'@keywords internal
#'@export
feat.sum.6.10 <- function(X){
  return(-0.1*(X[,6] + X[,7] + X[,8] + X[,9] + X[,10]))
}

#'@param X  a matrix
#'@rdname generate_dataset_advanced
#'@keywords internal
#'@export
neg.feat.sum.6.10 <- function(X){
  return(0.1*rowSums(X[,6:10]))
}

#'@param X  a matrix
#'@rdname generate_dataset_advanced
#'@keywords internal
#'@export
feat.complex <- function(X){
  return(X[,1]^2 + X[,2]^2 + rowSums(sin(3*X[,3:5])))
}
#'@param X  a matrix
#'@rdname generate_dataset_advanced
#'@keywords internal
#'@export
feat.square <- function(X){
  return(rowMeans(X[,1:5]^2))
}
#'@param X  a matrix
#'@rdname generate_dataset_advanced
#'@keywords internal
#'@export
feat.square.neg <- function(X){
  return(-rowMeans(X[,1:5]^2)/2)
}

#'Functions for generating simulated data
#'
#'@rdname generate_dataset_advanced
#'@keywords internal
#'@export
generate_dataset_advanced <- function(with.feat = F, data.size = 100, ncol.gene.mat = 100, feat.m = NA, feat.d = NA, mu.cen = 1.3,
                                      gamma1 = c(0,0),gamma2 = c(0,0),lam.m = 1/15, lam.d = 1/20, norm.vcov = c(1,.5,.5,1), cov = 0){
  # the models will be generated in the following way
  # log(Tm) = -log(small.lambda0) + gamma_1'X + h(Z) + eps_1
  # log(Td) = -log(small.lambda0) + gamma_2'X + h(Z) + eps_2
  # Td should always be greater than Tm

  # the Zs
  #Zgen <- rnorm(data.size*ncol.gene.mat)#.3*rnorm(10000) + .7*rcauchy(10000)
  #Zgen.mat <- matrix(Zgen, data.size, ncol.gene.mat)
  S <- matrix(cov, ncol.gene.mat, ncol.gene.mat);   diag(S) <- 1
  Zgen.mat <- mvtnorm::rmvnorm(data.size, sigma = S)

  # the X
  BMI <- rnorm(data.size, 21.7, 2.56)
  # distance between eyes in cm (haha!)
  DBE <- rnorm(data.size,5,1)

  eps <- mvtnorm::rmvnorm(n=data.size,mean=c(0,0),sigma=matrix(norm.vcov,2,2))

  eps <- log(-log(pnorm(eps)))

  if (with.feat == F){
    Tm <- exp(-log(lam.m) - gamma1%*%rbind(BMI, DBE) + eps[,1])
    Td <- exp(-log(lam.d) - gamma2%*%rbind(BMI, DBE) + eps[,2])
  }
  else{
    Tm <- exp(-log(lam.m) - gamma1%*%rbind(BMI, DBE) + feat.m(Zgen.mat) + eps[,1])
    Td <- exp(-log(lam.d) - gamma2%*%rbind(BMI, DBE) + feat.d(Zgen.mat) + eps[,2])
  }

  ## not sure if this is not violating the independent censoring assumption
  ## just exponential with fixed mean
  Tc <- rexp(data.size, 1/mu.cen)
  Tmat1 <- cbind(as.vector(Tc), as.vector(Tm), as.vector(Td))
  Tmin1 <- apply(Tmat1, 1, min)
  # sets up the right time, having in mind the coding in the write up
  delta1 <- apply(Tmat1 == Tmin1, 1, which) - 1

  Tmat2 <- cbind(as.vector(Tc), as.vector(Td))
  Tmin2 <- apply(Tmat2, 1, min)
  # sets up the right time, having in mind the coding in the write up
  delta2 <- apply(Tmat2 == Tmin2, 1, which) - 1

  ## adding new stuff to the dataset
  ## I've added the Tmin1, and the (delta1 != 0) as the failures, which practically gives
  ## an event if min(Tm, Td) <= Tc; which is exactly what other people do in this situation

  # this gives a data frame
  data <- data.frame(Tmin1, Tmin2, delta1, delta2, Zgen.mat, BMI, DBE, Tmin1, (delta1 != 0))

  return(data)
}