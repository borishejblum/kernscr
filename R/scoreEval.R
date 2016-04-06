
#'@export
scoreEval <- function(data, set_U, Cov_e_M_e_D, K = c("linear", "poly", "gaussian"), ...){
  data_M <- data[order(data[,1]),]
  data_D <- data_M[order(data_M[,2]),]

  M_M <- M_vec(Inf, data_M[,1], data_M[,3], c(0,0), as.matrix(data_M[,set_U]))
  M_D <- M_vec(Inf, data_D[,2], data_D[,4],  c(0,0), as.matrix(data_D[,set_U]))

  M_M <- M_M[,order(data_M[,2]),drop = F]

  K_rho <- kernelEval(t(data_D[,c(c(-1:-4), -set_U, -dim(data_D)[2] + 1, -dim(data_D)[2])]), K, ...)

  return ((M_M)%*%K_rho%*%t(M_M) + (M_D)%*%K_rho%*%t(M_D))
}

#'@export
scoreEval_pert <- function(num_perts, data, set_U, Cov_e_M_e_D, K = c("linear", "poly", "gaussian") , ...){

  Ws <- rnorm(dim(data)[1]*num_perts)
  Ws_mat <- matrix(Ws, dim(data)[1], num_perts)

  data_M <- data[order(data[,1]),]
  data_D <- data_M[order(data_M[,2]),]

  Ws_mat_M <- Ws_mat[order(data[,1]),]
  Ws_mat_D <- Ws_mat[order(data_M[,2]),]

  M_M <- M_vec_pert(Ws_mat_M, Inf, data_M[,1], data_M[,3], c(0,0), as.matrix(data_M[,set_U]))
  M_D <- M_vec_pert(Ws_mat_D, Inf, data_D[,2], data_D[,4], c(0,0), as.matrix(data_D[,set_U]))

  M_M <- M_M[order(data_M[,2]),]


  K_rho <- kernelEval(t(data_D[,c(c(-1:-4), -set_U, -dim(data_D)[2] + 1, -dim(data_D)[2])]), K, ...)

  return(diag(t(M_M)%*%K_rho%*%(M_M) + t(M_D)%*%K_rho%*%(M_D)))
}


#'@export
scoreEval_min_model <- function(data, set_U, Cov_e_M_e_D, K = c("linear", "poly", "gaussian"), ...){
  data_new <- data[order(data[,dim(data)[2]-1]),]

  M_new <- M_vec(Inf, data_new[,dim(data_new)[2]-1], data_new[,dim(data_new)[2]], c(0,0), as.matrix(data_new[,set_U]))
  K_rho <- kernelEval(t(data_new[,c(c(-1:-4), -set_U, -dim(data_new)[2] + 1, -dim(data_new)[2])]), K, ...)

  return((M_new)%*%K_rho%*%t(M_new))
}



#'@export
scoreEval_pert_min_model <- function(num_perts, data, set_U, Cov_e_M_e_D, K = c("linear", "poly", "gaussian"), ...){

  Ws <- rnorm(dim(data)[1]*num_perts)
  Ws_mat <- matrix(Ws, dim(data)[1], num_perts)

  data_new <- data[order(data[,dim(data)[2]-1]),]

  Ws_mat_new <- Ws_mat[order(data[,dim(data)[2]-1]),]

  M_new <- M_vec_pert(Ws_mat_new, Inf, data_new[,dim(data_new)[2]-1], data_new[,dim(data_new)[2]], c(0,0), as.matrix(data_new[,set_U]))

  K_rho <- kernelEval(t(data_new[,c(c(-1:-4), -set_U, -dim(data_new)[2] + 1, -dim(data_new)[2])]), K, ...)

  return(diag(t(M_new)%*%K_rho%*%(M_new)))
}