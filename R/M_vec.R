#'Compute the original gammas
#'
#'@keywords internal
#'@export
gamma_star <- function(perturb.mat, all_times, failures, gamma_vec, U){
  return(matrix(rep(gamma_vec, dim(perturb.mat)[2]), length(gamma_vec)))
}


#'@export
M_vec <- function(t, all_times, failures, gamma_vec, U){
  res <- c()
  # the failure times
  fail_times <- all_times[failures == 1]

  # the trick with the 0 infront is done if we have a failure before events
  lambdas_of_fail_times <- c(0, lambda(fail_times, all_times, failures, gamma_vec, U))

  time_mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event_time_mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- (event_time_mat <= time_mat) & (failures == 1)

  # this returns the max index of the max failure time prior to min(t, all_times[i])
  # here it throw an error if we don't have a single failure!!!
  # should handle that in a more delicate way
  fail_indexes <- apply(indexes,2,cumsum)

  # the + 1 here of the indexes is added if we ave failures before events refer to the logic above
  s <- as.vector(exp(U%*%gamma_vec))*(lambdas_of_fail_times[(fail_indexes + 1)])

  ans <- matrix(as.vector(indexes) - s, length(t), length(all_times), byrow = T)
  return(ans)
}


## the M vector perturbation
#'@export
M_vec_pert <- function(perturb_mat, t, all_times, failures, gamma_vec, U){
  res <- c()
  # the failure times
  fail_times <- all_times[failures == 1]

  ## comment the last two arguments for log lambda pert
  lambdas_of_fail_times <- rbind(rep(0, dim(perturb_mat)[2]), lambda_pert(fail_times, perturb_mat, all_times, failures, gamma_vec, U))

  time_mat <- matrix(rep(t, length(all_times)), length(all_times), length(t), byrow = T)
  event_time_mat <- matrix(rep(all_times, length(t)), length(all_times), length(t))

  indexes <- (event_time_mat <= time_mat) & (failures == 1)

  fail_indexes <- apply(indexes,2,cumsum)

  # the + 1 here of the indexes is added if we ave failures before events refer to the logic above
  lambdas_of_fail_times_prev <- c(0, lambda(fail_times, all_times, failures, gamma_vec, U))
  prev_s <- as.vector(exp(U%*%gamma_vec))*(lambdas_of_fail_times_prev[(fail_indexes + 1)])

  gamma_star <- gamma_star(perturb_mat, all_times, failures, gamma_vec, U)
  exp_matr <- as.matrix(exp(U%*%gamma_star))
  s <- exp_matr[rep(1:dim(exp_matr)[1], length(t)),]*(lambdas_of_fail_times[(fail_indexes + 1),])

  ans <- (as.vector(indexes) - s)*perturb_mat[rep(1:dim(perturb_mat)[1], length(t)),]
  perturbed_part_of_M <- s

  M_v <- matrix(as.vector(indexes) - prev_s, length(t), length(all_times), byrow = T)

  M_v_pert <- as.vector(t(M_v))*(perturb_mat)[rep(1:dim(perturb_mat)[1], dim(M_v)[1]),]

  ans <- M_v_pert - (perturbed_part_of_M - prev_s)
  return(ans)
}

