#'sup_stat_both
#'
#'sup_stat_both
#'
#'@param num_perts number of perturbations used
#'@param Ws
#'@param data
#'@param set_U
#'@param Cov_e_M_e_D
#'@parma rho Default is \code{1:40}.
#'@param l0
#'@param kernel a character string indicating which kernel should be used. Currently implemented are
#'\code{"linear"}, \code{"gaussian"} or \code{"poly"}.
#'@param d
#'@param est_gamma logical flag indicating whether \code{gamma} should be estimated. Default is \code{FALSE},
#'in which case \code{0} value is used for gamma.
#'@param pca_thres the threshold to be used for PCA. Default is \code{NULL}, in which case no PCA is performed.
#'
#'
#'@return a \code{data.frame}
#'
#'@export
sup_stat_both <- function(num_perts, Ws=NULL, data, set_U, Cov_e_M_e_D, rho = 1:40,
                          l0 = NA, kernel = c("gaussian", "poly"), d = NA, est_gamma=FALSE, pca_thres=NULL){

  if(is.null(Ws[1])){
    Ws <- rnorm(dim(data)[1]*num_perts)
  }
  n <- nrow(data)
  data_M <- data[order(data[,1]),]
  data_D <- data_M[order(data_M[,2]),]
  p_gene <- ncol(data) - 6

  ind.gene = 1:p_gene + 4 ## data_M sorted by X_M (col #1), data_D sorted by X_D (col #2)

  if(est_gamma){
    gamhat_M = coxph(Surv(data_M[,1],1*(data_M[,3]==1))~as.matrix(data_M[,set_U]))$coef
    gamhat_D = coxph(Surv(data_D[,2],1*(data_D[,4]==1))~as.matrix(data_D[,set_U]))$coef
  }else{
    gamhat_M = gamhat_D = rep(0,length(set_U))
  }
  M_Mc <- M_vec(Inf, data_M[,1], 1*(data_M[,3]==1), gamhat_M, as.matrix(data_M[,set_U])) ## cause specific hazard for M
  M_Dc <- M_vec(Inf, data_M[,1], 1*(data_M[,3]==2), gamhat_M, as.matrix(data_M[,set_U])) ## cause specific hazard for D
  M_Dm <- M_vec(Inf, data_D[,2], 1*(data_D[,4]==1), gamhat_D, as.matrix(data_D[,set_U])) ## marginal hazard for D
  M_Mc <- M_Mc[, order(data_M[,2]), drop = F] ## this ensures that sorting is done according to sorting of data_D (col #2) ##
  M_Dc <- M_Dc[, order(data_M[,2]), drop = F] ## this ensures that sorting is done according to sorting of data_D (col #2) ##

  if(kernel=="linear"){K_rho=kernelEval(t(data_D[,ind.gene]), K = "linear"); K_rho = matrix(K_rho,nrow=1) }
  else if(kernel=="gaussian"){K_rho=gaussKernelEval_multipleSigmas(t(data_D[,ind.gene]), sigma=rho)}
  else if(kernel=="poly"){K_rho=polyKernelEval_multiple( t(data_D[,ind.gene]), a=rho, d=d)}
  sqrt_ncol_K_rho <- sqrt(ncol(K_rho))
  K_rho <- matrix(c(t(K_rho)), ncol = sqrt_ncol_K_rho, byrow = T)
  #dim(K_rho)

  # do PCA
  if(is.null(pca_thres)){
    warning("Not performing PCA: potential loss of power, especially on finite samples")
  }else{
    if(kernel=="linear"){
      K_rho_eig <- eigen(K_rho)
      ncomp <- which(cumsum(K_rho_eig$values/sum(K_rho_eig$values))>=pca_thres)[1]
      K_rho <- tcrossprod(matrix(rep(sqrt(K_rho_eig$values[1:ncomp]), sqrt_ncol_K_rho),
                                 nrow = sqrt_ncol_K_rho, byrow = T)*K_rho_eig$vectors[, 1:ncomp])
    }else{
      K_rho_new <- NULL
      for (i in 1:length(rho)){
        K_rho_eig <- eigen(K_rho[1:sqrt_ncol_K_rho+(i-1)*sqrt_ncol_K_rho, ])
        ncomp <- which(cumsum(K_rho_eig$values/sum(K_rho_eig$values))>=pca_thres)[1]
        K_rho_new <- rbind(K_rho_new,
                           tcrossprod(matrix(rep(sqrt(K_rho_eig$values[1:ncomp]), sqrt_ncol_K_rho),
                                             nrow = sqrt_ncol_K_rho, byrow = T)*K_rho_eig$vectors[, 1:ncomp])
        )
      }
      K_rho <- K_rho_new
    }
  }

  # returns a function to obtaing the K_rhos for a fixed matrix K_rho of the type above
  obtain.K_rhos <- function(K_rho, M_vec){
    temp <- tcrossprod(matrix(tcrossprod(M_vec, K_rho), ncol = sqrt_ncol_K_rho, byrow = T), M_vec)
  }
  ## including cause specific hazard for M, D as well as marginal for D ##
  stats.all <- cbind("Mc"=obtain.K_rhos(K_rho, M_Mc), "Dc"=obtain.K_rhos(K_rho, M_Dc), "Dm"=obtain.K_rhos(K_rho, M_Dm))

  Ws.mat <- matrix(Ws, dim(data)[1], num_perts); Ws.mat.M <- Ws.mat[order(data[,1]),]; Ws.mat.D <- Ws.mat[order(data_M[,2]),]

  ### using the new method based on cause-specific hazard model ###
  M_Mc_pert <- M_vec_pert(Ws.mat.M, Inf, data_M[,1], 1*(data_M[,3]==1), gamhat_M, as.matrix(data_M[,set_U])); ## cause specific for M
  M_Dc_pert <- M_vec_pert(Ws.mat.M, Inf, data_M[,1], 1*(data_M[,3]==2), gamhat_M, as.matrix(data_M[,set_U])); ## cause specific for D
  M_Mc_pert <- M_Mc_pert[order(data_M[,2]),]; ## this ensures that sorting is done according to sorting of data_D
  M_Dc_pert <- M_Dc_pert[order(data_M[,2]),]; ## this ensures that sorting is done according to sorting of data_D
  M_Dm_pert <- M_vec_pert(Ws.mat.D, Inf, data_D[,2], 1*(data_D[,4]==1), gamhat_D, as.matrix(data_D[,set_U])); ## marginal for D
  part_1.Mc <- t(M_Mc_pert)%*%t(K_rho); part_1.Dc <- t(M_Dc_pert)%*%t(K_rho); part_1.Dm <- t(M_Dm_pert)%*%t(K_rho)
  part_2.Mc <- c(t(M_Mc_pert))*c(part_1.Mc); part_2.Dc <- c(t(M_Dc_pert))*c(part_1.Dc); part_2.Dm <- c(t(M_Dm_pert))*c(part_1.Dm)
  part_3.list <- list("Mc"=matrix(part_2.Mc, ncol = num_perts, byrow = T), ## cause specific M, n*K_rho rows: (1:n), n+(1:n), 2n+(1:n)
                      "Dc"=matrix(part_2.Dc, ncol = num_perts, byrow = T), ## cause specific D;
                      "Dm"=matrix(part_2.Dm, ncol = num_perts, byrow = T)) ## marginal D
  stat_pert_std.list = as.list(1:3); names(stat_pert_std.list) = names(part_3.list); sd.all= NULL
  for(l in 1:3){tmpres = NULL; for(i in 0:(length(rho) - 1)){tmpres = cbind(tmpres, colSums(part_3.list[[l]][i*n + 1:n,]))}; ## num_perts x n.rho matrix
  tmp.sd = apply(tmpres,2,sd); stat_pert_std.list[[l]]=as.matrix(tmpres/VTM(tmp.sd,num_perts)); sd.all = cbind(sd.all, tmp.sd)}
  stats.all_std = stats.all/sd.all; colnames(stats.all_std) = names(part_3.list) ## standardized statistic for Mc, Dc and D ##

  ## sum over cause M, cause D and marg D, then take max over rho ##
  stats.sum3_std = max(apply(stats.all_std,1,sum))
  perts_sum3_std = apply(stat_pert_std.list[["Mc"]]+stat_pert_std.list[["Dc"]]+stat_pert_std.list[["Dm"]],1,max)
  pvals_sum3_std = mean(perts_sum3_std > stats.sum3_std)
  ## max over rho, then over cause M, cause D and marg D ##
  stats.max3_std = max(apply(stats.all_std,2,max)) ## maximum statistic across rho and then outcome
  perts_max3_std = apply(matrix(unlist(lapply(stat_pert_std.list,function(xx){apply(xx,1,max)})),ncol=3),1,max)
  pvals_max3_std = mean(perts_max3_std > stats.max3_std)
  ## max over rho, then over cause M, cause D and marg D ##
  stats.max2_std = max(apply(stats.all_std[,1:2,drop=F],2,max)) ## maximum statistic across rho and then outcome
  perts_max2_std = apply(matrix(unlist(lapply(stat_pert_std.list,function(xx){apply(xx,1,max)})),ncol=3)[,1:2],1,max)
  pvals_max2_std = mean(perts_max2_std > stats.max2_std)
  ## sum over cause M and cause D then take max over rho (check how much additional info gained by marg D) ##
  stats.McDc_std = max(apply(stats.all_std[,c("Mc","Dc"),drop=F],1,sum))
  perts_McDc_std = apply(stat_pert_std.list[["Mc"]]+stat_pert_std.list[["Dc"]],1,max)
  pvals_McDc_std = mean(perts_McDc_std > stats.McDc_std)
  ## sum over cause M and cause D then take max over rho (check how much additional info gained by marg D) ##
  stats.McDM_std = max(apply(stats.all_std[,c("Mc","Dm"),drop=F],1,sum))
  perts_McDM_std = apply(stat_pert_std.list[["Mc"]]+stat_pert_std.list[["Dm"]],1,max)
  pvals_McDM_std = mean(perts_McDM_std > stats.McDM_std)
  ## take the maximum statistic across rho for each of the outcome, cause M, cause D and marg D
  stats.max_std = apply(stats.all_std,2,max) ## maximum statistic across rho
  perts_max_std = matrix(unlist(lapply(stat_pert_std.list,function(xx){apply(xx,1,max)})),ncol=3) ## num_perts x 3 matrix
  pvals_max_std = apply(perts_max_std > VTM(stats.max_std,num_perts),2,mean)
  names(pvals_max_std) = c("Mc","Dc","Dm")



  ### using time to min(death,progression) with a single cox model ###
  data_Min <- data[order(data[,dim(data)[2]-1]),]
  if(est_gamma){
    gamhat_min <- coxph(Surv(data_Min[,ncol(data_Min)-1],
                            1*(data_Min[,ncol(data_Min)]))~as.matrix(data_Min[,set_U]))$coef
  }else{
    gamhat_min <- rep(0,length(set_U))
  }
  Ws.mat_min <- Ws.mat[order(data[,dim(data)[2]-1]),]
  M_min <- M_vec(Inf, data_Min[,dim(data_Min)[2]-1],
                 data_Min[,dim(data_Min)[2]], gamhat_min, as.matrix(data_Min[,set_U]))

  if(kernel=="linear"){
    K_rho=kernelEval(t(data_Min[,ind.gene]), K = "linear")
    K_rho = matrix(K_rho,nrow=1)
  }
  if(kernel=="gaussian"){
    K_rho=gaussKernelEval_multipleSigmas(t(data_Min[,ind.gene]), sigma=rho)
  }
  if(kernel=="poly"){
    K_rho=polyKernelEval_multiple(t(data_Min[,ind.gene]), a=rho, d=d)
  }

  sqrt_ncol_K_rho <- sqrt(ncol(K_rho))
  K_rho <- matrix(c(t(K_rho)), ncol = sqrt_ncol_K_rho, byrow = T)

  # do PCA
  if(kernel=="linear"){
    K_rho_eig <- eigen(K_rho)
    ncomp <- which(cumsum(K_rho_eig$values/sum(K_rho_eig$values))>=pca_thres)[1]
    K_rho <- tcrossprod(matrix(rep(sqrt(K_rho_eig$values[1:ncomp]),
                                   sqrt_ncol_K_rho), nrow = sqrt_ncol_K_rho, byrow = T)*K_rho_eig$vectors[, 1:ncomp])
  }else{
    K_rho_new <- NULL
    for (i in 1:length(rho)){
      K_rho_eig <- eigen(K_rho[1:sqrt_ncol_K_rho+(i-1)*sqrt_ncol_K_rho, ])
      ncomp <- which(cumsum(K_rho_eig$values/sum(K_rho_eig$values))>=pca_thres)[1]
      K_rho_new <- rbind(K_rho_new,
                         tcrossprod(matrix(rep(sqrt(K_rho_eig$values[1:ncomp]), sqrt_ncol_K_rho),
                                           nrow = sqrt_ncol_K_rho, byrow = T)*K_rho_eig$vectors[, 1:ncomp])
      )
    }
    K_rho <- K_rho_new
  }

  stats_min <- obtain.K_rhos(K_rho, M_min)
  M_min_pert <- M_vec_pert(Ws.mat_min, Inf, data_Min[,dim(data_Min)[2]-1], data_Min[,dim(data_Min)[2]], gamhat_min, as.matrix(data_Min[,set_U]))
  part_1_min <- t(M_min_pert)%*%t(K_rho)
  part_2_min <- c(t(M_min_pert))*c(part_1_min)
  part_3_min <- matrix(part_2_min, ncol = num_perts, byrow = T)
  res_min <- c(); for(i in 0:(length(rho) - 1)){res_min <- cbind(res_min, colSums(part_3_min[i*n + 1:n,]))}
  st_devs_min <- sqrt(colSums((res_min - matrix(rep(colMeans(res_min), num_perts), nrow = nrow(res_min), byrow=T))^2)/(num_perts - 1))
  pval_min = sum(apply(t(res_min)/st_devs_min,2,max) >= max(stats_min/st_devs_min))/num_perts

  if(est_gamma){
    res_gamma=c("M"=gamhat_M, "D"=gamhat_D, "min"=gamhat_min)
  }else{
    res_gamma=NULL
  }

  res_pvals <- c("sum3" = pvals_sum3_std,
                 "McDc" = pvals_McDc_std,
                 "McDm" = pvals_McDM_std,
                 "max3" = pvals_max3_std,
                 "max.McDc" = pvals_max2_std,
                 pvals_max_std, #TODO
                 "PFS" = pval_min)

  res_pertbs <- cbind.data.frame("sum3" =  1 - rank(perts_sum3_std)/num_perts,
                                 "McDc" = 1 - rank(perts_McDc_std)/num_perts,
                                 "McDm" = 1 - rank(perts_McDM_std)/num_perts,
                                 "max3" =  1 - rank(perts_max3_std)/num_perts,
                                 "max.McDc" = 1 - rank(perts_max2_std)/num_perts,
                                 "Mc" = 1 - rank(perts_max_std[,1])/num_perts,
                                 "Dc" = 1 - rank(perts_max_std[,2])/num_perts,
                                 "Dm" = 1 - rank(perts_max_std[,3])/num_perts,
                                 "PFS" = 1 - rank(apply(t(res_min)/st_devs_min,2,max))/num_perts)

  return(list("raw_pvals"=res_pvals, "null_pvals_pertbs"=res_pertbs, "gammahat"=res_gamma))
}