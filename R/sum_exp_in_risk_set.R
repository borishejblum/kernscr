# SAME AS PI.0
#'@export
sum_exp_in_risk_set <- function(t, all_times, gamma_vec, U){

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



