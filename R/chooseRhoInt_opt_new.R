## choose rho in the range such that the eigenvalue decay rate is between 1.1 and 2.2 ##
## 1.1 to 4
chooseRhoInt_opt_new <- function(Z, rho, rho0, l0 = NA, kernel = c("gaussian", "poly"), d = NA, rate.range=c(1.5,4)){
  n <- nrow(Z)
  if (kernel == "gaussian"){
    G <- gaussKernelEval.multipleSigmas(Z, sigma = rho)
  }
  if (kernel == "poly"){
    G <- polyKernelEval.multipleSigmas(Z, a = rho, d = d)
  }

  slope_eigen <- function(ind, G. = G, Z. = Z, rho0. = rho0){
    eigenvals <- eigen(matrix(G.[ind,], ncol(Z.)), symmetric = T, only.values=T)$values
    log_eig <-  log(eigenvals[1:sum(cumsum(eigenvals)/sum(eigenvals) <= rho0.)])
    log_j <-  log(1:length(log_eig))
    slope <- try(-MASS::rlm(log_eig~log_j)$coef[2], silent=TRUE)
    if(substring(slope[1],1,5)=="Error"){
      slope <- NA
    }else{
      print(c(rho[ind],slope,length(log_j)))
    }
    return(c(slope,length(log_j)))
  }

  slope_vec <-  na.omit(cbind("rho"=rho,"slope"=t(sapply(1:length(rho), slope_eigen))))
  colnames(slope_vec)[2:3] <-  c("slope","m.keep")
  print(slope_vec)

  rho <- slope_vec[,"rho"] ## changed
  res_rho_opt <- c(max(c(min(rho),slope_vec[slope_vec[,"slope"]<=rate.range[1],"rho"])),min(c(max(rho), slope_vec[slope_vec[,"slope"]>=rate.range[2],"rho"]))) ## changed
  return(res_rho_opt)
}