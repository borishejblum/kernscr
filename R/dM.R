#'@export
dM <- function(all_times, failures, gamma_vec, U){
  fail_times.indexes <- which(failures == 1)
  fail_times <- all_times[fail_times.indexes]

  almost_ans <- M_vec(fail_times, all_times, failures, gamma_vec, U)

  ans <- almost_ans - rbind(rep(0, length(all_times)), almost_ans[-dim(almost_ans)[1],])
  return(ans)
}