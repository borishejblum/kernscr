#'@importFrom survival Surv coxph
#'@keywords internal
#'@export
get_models <- function(dat_fr, set_U){
  # note that I am always expecting the indexes to be in this order!!!
  times1 <- dat_fr[,1]
  times2 <- dat_fr[,2]
  # following the notation of the Surv function
  # in the write up 1 means dead and 0 means alive
  events1 <- dat_fr[,3] == 1
  events2 <- dat_fr[,4]
  model1 <- survival::Surv(times1, events1)
  model2 <- survival::Surv(times2, events2)
  gamma1 <- survival::coxph(model1 ~ as.matrix(dat_fr[,set_U]))$coef
  gamma2 <- survival::coxph(model2 ~ as.matrix(dat_fr[,set_U]))$coef
  return(cbind(gamma1, gamma2))
}