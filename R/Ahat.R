# THE SAME AS sum.exp.in.risk.set
#'@export
PI_0 <- function(t, all_times, gamma_vec, U){

  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time.mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event.time.mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- matrix(rep(1:length(all_times), length(t)), length(all_times), length(t))*(event.time.mat >= time.mat)

  almost.the.right.res <- as.vector(cumsum(exp(as.matrix(U[indexes, ])%*%gamma_vec))[cumsum(colSums(event.time.mat >= time.mat))])
  res <- almost.the.right.res - c(0, almost.the.right.res[-length(almost.the.right.res)])
  # discards the result for the -1
  return(res[-1])
}

# returns a n x
#'@export
PI_1 <- function(t, all_times, gamma_vec, U, U.cust = U){

  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time.mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event.time.mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- matrix(rep(1:length(all_times), length(t)), length(all_times), length(t))*(event.time.mat >= time.mat)

  almost.the.right.res <- apply(as.matrix(as.vector(exp(as.matrix(U[indexes, ])%*%gamma_vec))*U.cust[indexes, ], ncol = length(gamma_vec)), 2, cumsum)[cumsum(colSums(event.time.mat >= time.mat)), ]
  almost.the.right.res <- as.matrix(almost.the.right.res, nrow = length(t))
  res <- almost.the.right.res - rbind(0, almost.the.right.res[-dim(almost.the.right.res)[1], , drop = F])
  # discards the result for the -1
  return(matrix(res[-1,], nrow = length(t) - 1))
}

#'@export
PI_2 <- function(t, all_times, gamma_vec, U){
  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time.mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event.time.mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- matrix(rep(1:length(all_times), length(t)), length(all_times), length(t))*(event.time.mat >= time.mat)

  # this is a matrix dim(U)[2]^2 x sum(event.time.mat >= time.mat)
  # it contains the matrices of multiplying the transposed rows of U by themselves (cols of t(U))
  # each column of this matrix is the combined rows of a matrices described above
  l <- dim(as.matrix(U))[2]
  simplified.mat <- matrix(U[indexes, rep(1:l, l)]*U[indexes, rep(1:l, rep(l,l))], ncol = l^2)

  almost.the.right.res <- apply(as.vector(exp(as.matrix(U[indexes, ])%*%gamma_vec))*(simplified.mat), 2, cumsum)[cumsum(colSums(event.time.mat >= time.mat)), ]
  almost.the.right.res <- as.matrix(almost.the.right.res, ncol = l^2)
  res <- almost.the.right.res - rbind(0, as.matrix(almost.the.right.res[-dim(almost.the.right.res)[1], , drop = F]))
  # discards the result for the -1
  return(matrix(res[-1,], ncol = l^2))
}

# Ahat computation
# The idea behind this calculation is:
# dEN_i(t) is estimated by the mean number of failure times before i.e. multiply by 1/n (typically dN_i(t) = 1 only for 1 out of n times)
# The second derivative of the log partial likelihood
Ahat <- function(all_times, failures, gamma_vec, U){
  # be careful with the definitions here
  n <- length(all_times)

  fail_times <- all_times[failures == 1]
  #
  pi_0 <- PI_0(fail_times, all_times, gamma_vec, U)/n
  #
  pi_1 <- PI_1(fail_times, all_times, gamma_vec, U)/n

  l <- dim(as.matrix(pi_1))[2]
  pi_1_x2 <- as.matrix(pi_1)[,rep(1:l, l)]*as.matrix(pi_1)[,rep(1:l, rep(l,l))]

  pi_2 <- PI_2(fail_times, all_times, gamma_vec, U)/n

  res <- colSums(pi_0^(-2)*(pi_2*pi_0 - pi_1_x2))

  return(matrix(as.vector(t(res)), ncol = dim(U)[2]))
}

