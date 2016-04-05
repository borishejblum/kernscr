sup_stat_both <- function(num.perts, Ws=NULL, data, set_U, Cov_e_M_e_D, rho = 1:40,
                          l0 = NA, kernel = c("gaussian", "poly"), d = NA, est_gamma=FALSE, pca_thres=NULL){
  if(is.null(Ws[1])){
    Ws <- rnorm(dim(data)[1]*num.perts)
  }
  n <- nrow(data)
  data.M <- data[order(data[,1]),]
  data.D <- data.M[order(data.M[,2]),]
  p.gene <- ncol(data) - 6 # TODO is that correct ?
  ind.gene = 1:p.gene + 4 ## data.M sorted by X_M (col #1), data.D sorted by X_D (col #2)
  if(est_gamma){
    gamhat.M = coxph(Surv(data.M[,1],1*(data.M[,3]==1))~as.matrix(data.M[,set_U]))$coef
    gamhat.D = coxph(Surv(data.D[,2],1*(data.D[,4]==1))~as.matrix(data.D[,set_U]))$coef
  }else{
    gamhat.M = gamhat.D = rep(0,length(set_U))
  }
  M.Mc <- M.vec(Inf, data.M[,1], 1*(data.M[,3]==1), gamhat.M, as.matrix(data.M[,set_U])) ## cause specific hazard for M
  M.Dc <- M.vec(Inf, data.M[,1], 1*(data.M[,3]==2), gamhat.M, as.matrix(data.M[,set_U])) ## cause specific hazard for D
  M.Dm <- M.vec(Inf, data.D[,2], 1*(data.D[,4]==1), gamhat.D, as.matrix(data.D[,set_U])) ## marginal hazard for D
  M.Mc <- M.Mc[, order(data.M[,2]), drop = F] ## this ensures that sorting is done according to sorting of data.D (col #2) ##
  M.Dc <- M.Dc[, order(data.M[,2]), drop = F] ## this ensures that sorting is done according to sorting of data.D (col #2) ##

  if(kernel=="linear"){K.rho=kernelEval(t(data.D[,ind.gene]), K = "linear"); K.rho = matrix(K.rho,nrow=1) }
  else if(kernel=="gaussian"){K.rho=gaussKernelEval.multipleSigmas(t(data.D[,ind.gene]), sigma=rho)}
  else if(kernel=="poly"){K.rho=polyKernelEval.multipleSigmas( t(data.D[,ind.gene]), a=rho, d=d)}
  sqrt.ncol.K.rho <- sqrt(ncol(K.rho))
  K.rho <- matrix(c(t(K.rho)), ncol = sqrt.ncol.K.rho, byrow = T)
  #dim(K.rho)

  # do PCA
  if(is.null(pca_thres)){
    warning("Not performing PCA: potential loss of power, especially on finite samples")
  }else{
    if(kernel=="linear"){
      K.rho_eig <- eigen(K.rho)
      ncomp <- which(cumsum(K.rho_eig$values/sum(K.rho_eig$values))>=pca_thres)[1]
      K.rho <- tcrossprod(matrix(rep(sqrt(K.rho_eig$values[1:ncomp]), sqrt.ncol.K.rho), nrow = sqrt.ncol.K.rho, byrow = T)*K.rho_eig$vectors[, 1:ncomp])
    }else{
      K.rho_new <- NULL
      for (i in 1:length(rho)){
        K.rho_eig <- eigen(K.rho[1:sqrt.ncol.K.rho+(i-1)*sqrt.ncol.K.rho, ])
        ncomp <- which(cumsum(K.rho_eig$values/sum(K.rho_eig$values))>=pca_thres)[1]
        K.rho_new <- rbind(K.rho_new,
                           tcrossprod(matrix(rep(sqrt(K.rho_eig$values[1:ncomp]), sqrt.ncol.K.rho), nrow = sqrt.ncol.K.rho, byrow = T)*K.rho_eig$vectors[, 1:ncomp])
        )
      }
      K.rho <- K.rho_new
    }
  }

  # returns a function to obtaing the K.rhos for a fixed matrix K.rho of the type above
  obtain.K.rhos <- function(K.rho, M.vec){
    temp <- tcrossprod(matrix(tcrossprod(M.vec, K.rho), ncol = sqrt.ncol.K.rho, byrow = T), M.vec)
  }
  ## including cause specific hazard for M, D as well as marginal for D ##
  stats.all <- cbind("Mc"=obtain.K.rhos(K.rho, M.Mc), "Dc"=obtain.K.rhos(K.rho, M.Dc), "Dm"=obtain.K.rhos(K.rho, M.Dm))

  Ws.mat <- matrix(Ws, dim(data)[1], num.perts); Ws.mat.M <- Ws.mat[order(data[,1]),]; Ws.mat.D <- Ws.mat[order(data.M[,2]),]

  ### using the new method based on cause-specific hazard model ###
  M.Mc.pert <- M.vec.pert(Ws.mat.M, Inf, data.M[,1], 1*(data.M[,3]==1), gamhat.M, as.matrix(data.M[,set_U])); ## cause specific for M
  M.Dc.pert <- M.vec.pert(Ws.mat.M, Inf, data.M[,1], 1*(data.M[,3]==2), gamhat.M, as.matrix(data.M[,set_U])); ## cause specific for D
  M.Mc.pert <- M.Mc.pert[order(data.M[,2]),]; ## this ensures that sorting is done according to sorting of data.D
  M.Dc.pert <- M.Dc.pert[order(data.M[,2]),]; ## this ensures that sorting is done according to sorting of data.D
  M.Dm.pert <- M.vec.pert(Ws.mat.D, Inf, data.D[,2], 1*(data.D[,4]==1), gamhat.D, as.matrix(data.D[,set_U])); ## marginal for D
  part.1.Mc <- t(M.Mc.pert)%*%t(K.rho); part.1.Dc <- t(M.Dc.pert)%*%t(K.rho); part.1.Dm <- t(M.Dm.pert)%*%t(K.rho)
  part.2.Mc <- c(t(M.Mc.pert))*c(part.1.Mc); part.2.Dc <- c(t(M.Dc.pert))*c(part.1.Dc); part.2.Dm <- c(t(M.Dm.pert))*c(part.1.Dm)
  part.3.list <- list("Mc"=matrix(part.2.Mc, ncol = num.perts, byrow = T), ## cause specific M, n*K.rho rows: (1:n), n+(1:n), 2n+(1:n)
                      "Dc"=matrix(part.2.Dc, ncol = num.perts, byrow = T), ## cause specific D;
                      "Dm"=matrix(part.2.Dm, ncol = num.perts, byrow = T)) ## marginal D
  stat.pert.std.list = as.list(1:3); names(stat.pert.std.list) = names(part.3.list); sd.all= NULL
  for(l in 1:3){tmpres = NULL; for(i in 0:(length(rho) - 1)){tmpres = cbind(tmpres, colSums(part.3.list[[l]][i*n + 1:n,]))}; ## num.perts x n.rho matrix
  tmp.sd = apply(tmpres,2,sd); stat.pert.std.list[[l]]=as.matrix(tmpres/VTM(tmp.sd,num.perts)); sd.all = cbind(sd.all, tmp.sd)}
  stats.all.std = stats.all/sd.all; colnames(stats.all.std) = names(part.3.list) ## standardized statistic for Mc, Dc and D ##

  ## sum over cause M, cause D and marg D, then take max over rho ##
  stats.sum3.std = max(apply(stats.all.std,1,sum))
  perts.sum3.std = apply(stat.pert.std.list[["Mc"]]+stat.pert.std.list[["Dc"]]+stat.pert.std.list[["Dm"]],1,max)
  pvals.sum3.std = mean(perts.sum3.std > stats.sum3.std)
  ## max over rho, then over cause M, cause D and marg D ##
  stats.max3.std = max(apply(stats.all.std,2,max)) ## maximum statistic across rho and then outcome
  perts.max3.std = apply(matrix(unlist(lapply(stat.pert.std.list,function(xx){apply(xx,1,max)})),ncol=3),1,max)
  pvals.max3.std = mean(perts.max3.std > stats.max3.std)
  ## max over rho, then over cause M, cause D and marg D ##
  stats.max2.std = max(apply(stats.all.std[,1:2,drop=F],2,max)) ## maximum statistic across rho and then outcome
  perts.max2.std = apply(matrix(unlist(lapply(stat.pert.std.list,function(xx){apply(xx,1,max)})),ncol=3)[,1:2],1,max)
  pvals.max2.std = mean(perts.max2.std > stats.max2.std)
  ## sum over cause M and cause D then take max over rho (check how much additional info gained by marg D) ##
  stats.McDc.std = max(apply(stats.all.std[,c("Mc","Dc"),drop=F],1,sum))
  perts.McDc.std = apply(stat.pert.std.list[["Mc"]]+stat.pert.std.list[["Dc"]],1,max)
  pvals.McDc.std = mean(perts.McDc.std > stats.McDc.std)
  ## sum over cause M and cause D then take max over rho (check how much additional info gained by marg D) ##
  stats.McDm.std = max(apply(stats.all.std[,c("Mc","Dm"),drop=F],1,sum))
  perts.McDm.std = apply(stat.pert.std.list[["Mc"]]+stat.pert.std.list[["Dm"]],1,max)
  pvals.McDm.std = mean(perts.McDm.std > stats.McDm.std)
  ## take the maximum statistic across rho for each of the outcome, cause M, cause D and marg D
  stats.max.std = apply(stats.all.std,2,max) ## maximum statistic across rho
  perts.max.std = matrix(unlist(lapply(stat.pert.std.list,function(xx){apply(xx,1,max)})),ncol=3) ## num.perts x 3 matrix
  pvals.max.std = apply(perts.max.std > VTM(stats.max.std,num.perts),2,mean)
  names(pvals.max.std) = c("Mc","Dc","Dm")



  ### using time to min(death,progression) with a single cox model ###
  data.min <- data[order(data[,dim(data)[2]-1]),]
  if(est_gamma){gamhat.min = coxph(Surv(data.min[,ncol(data.min)-1],1*(data.min[,ncol(data.min)]))~as.matrix(data.min[,set_U]))$coef}else{gamhat.min = rep(0,length(set_U))}
  Ws.mat.min <- Ws.mat[order(data[,dim(data)[2]-1]),]
  M.min <- M.vec(Inf, data.min[,dim(data.min)[2]-1], data.min[,dim(data.min)[2]], gamhat.min, as.matrix(data.min[,set_U]))

  if(kernel=="linear"){K.rho=kernelEval(t(data.min[,ind.gene]), K = "linear"); K.rho = matrix(K.rho,nrow=1) }
  if(kernel=="gaussian"){K.rho=gaussKernelEval.multipleSigmas(t(data.min[,ind.gene]),sigma=rho)}
  if(kernel=="poly"){    K.rho=polyKernelEval.multipleSigmas( t(data.min[,ind.gene]),a=rho,d=d)}
  sqrt.ncol.K.rho <- sqrt(ncol(K.rho))
  K.rho <- matrix(c(t(K.rho)), ncol = sqrt.ncol.K.rho, byrow = T)
  #dim(K.rho)

  # do PCA
  if(kernel=="linear"){
    K.rho_eig <- eigen(K.rho)
    ncomp <- which(cumsum(K.rho_eig$values/sum(K.rho_eig$values))>=pca_thres)[1]
    K.rho <- tcrossprod(matrix(rep(sqrt(K.rho_eig$values[1:ncomp]), sqrt.ncol.K.rho), nrow = sqrt.ncol.K.rho, byrow = T)*K.rho_eig$vectors[, 1:ncomp])
  }else{
    K.rho_new <- NULL
    for (i in 1:length(rho)){
      K.rho_eig <- eigen(K.rho[1:sqrt.ncol.K.rho+(i-1)*sqrt.ncol.K.rho, ])
      ncomp <- which(cumsum(K.rho_eig$values/sum(K.rho_eig$values))>=pca_thres)[1]
      K.rho_new <- rbind(K.rho_new,
                         tcrossprod(matrix(rep(sqrt(K.rho_eig$values[1:ncomp]), sqrt.ncol.K.rho), nrow = sqrt.ncol.K.rho, byrow = T)*K.rho_eig$vectors[, 1:ncomp])
      )
    }
    K.rho <- K.rho_new
  }

  stats.min <- obtain.K.rhos(K.rho, M.min)
  M.min.pert <- M.vec.pert(Ws.mat.min, Inf, data.min[,dim(data.min)[2]-1], data.min[,dim(data.min)[2]], gamhat.min, as.matrix(data.min[,set_U]))
  part.1.min <- t(M.min.pert)%*%t(K.rho); part.2.min <- c(t(M.min.pert))*c(part.1.min)
  part.3.min <- matrix(part.2.min, ncol = num.perts, byrow = T)
  res.min <- c(); for(i in 0:(length(rho) - 1)){res.min <- cbind(res.min,colSums(part.3.min[i*n + 1:n,]))}
  st.devs.min <- sqrt(colSums((res.min - matrix(rep(colMeans(res.min),num.perts),nrow = nrow(res.min),byrow=T))^2)/(num.perts - 1))
  pval.min = sum(apply(t(res.min)/st.devs.min,2,max) >= max(stats.min/st.devs.min))/num.perts

  if(est_gamma){
    res_gamma=c("M"=gamhat.M, "D"=gamhat.D, "min"=gamhat.min)
  }else{
    res_gamma=NULL
  }
  res_pvals <- c("sum3"=pvals.sum3.std,
                 "McDc"=pvals.McDc.std,
                 "McDm"=pvals.McDm.std,
                 "max3"=pvals.max3.std,
                 "max.McDc"=pvals.max2.std,
                 pvals.max.std,
                 "PFS"=pval.min)
  res_pertbs <- cbind.data.frame("sum3" =  1 - rank(perts.sum3.std)/num.perts,
                                 "McDc" = 1 - rank(perts.McDc.std)/num.perts,
                                 "McDm" = 1 - rank(perts.McDm.std)/num.perts,
                                 "max3" =  1 - rank(perts.max3.std)/num.perts,
                                 "max.McDc" = 1 - rank(perts.max2.std)/num.perts,
                                 "Mc" = 1 - rank(perts.max.std[,1])/num.perts,
                                 "Dc" = 1 - rank(perts.max.std[,2])/num.perts,
                                 "Dm" = 1 - rank(perts.max.std[,3])/num.perts,
                                 "PFS" = 1 - rank(apply(t(res.min)/st.devs.min,2,max))/num.perts)
  return(list("raw_pvals"=res_pvals, "null_pvals_pertbs"=res_pertbs, "gammahat"=res_gamma))
}