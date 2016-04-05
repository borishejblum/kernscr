
# SAME AS PI.0
sum.exp.in.risk.set <- function(t, all.times, gamma.vec, U){

  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time.mat <- matrix(rep(t, length(all.times)), length(all.times), length(t), byrow = T)
  event.time.mat <- matrix(rep(all.times, length(t)), length(all.times), length(t))

  indexes <- matrix(rep(1:length(all.times), length(t)), length(all.times), length(t))*(event.time.mat >= time.mat)

  almost.the.right.res <- as.vector(cumsum(exp(as.matrix(U[indexes, ])%*%gamma.vec))[cumsum(colSums(event.time.mat >= time.mat))])
  res <- almost.the.right.res - c(0, almost.the.right.res[-length(almost.the.right.res)])
  # discards the result for the -1
  return(res[-1])
}

lambda <- function(t, all.times, failures, gamma.vec, U){

  fail.times <- all.times[failures == 1]

  sums.fail.times <- sum.exp.in.risk.set(fail.times, all.times, gamma.vec, U)
  denoms <- 1/sums.fail.times

  time.mat <- matrix(rep(t, length(fail.times)), length(fail.times), length(t), byrow = T)
  event.time.mat <- matrix(rep(fail.times, length(t)), length(fail.times), length(t))

  return(t(event.time.mat <= time.mat)%*%denoms)
}

M.vec <- function(t, all.times, failures, gamma.vec, U){
  res <- c()
  # the failure times
  fail.times <- all.times[failures == 1]

  # the trick with the 0 infront is done if we have a failure before events
  lambdas.of.fail.times <- c(0, lambda(fail.times, all.times, failures, gamma.vec, U))

  time.mat <- matrix(rep(t, length(all.times)), length(all.times), length(t), byrow = T)
  event.time.mat <- matrix(rep(all.times, length(t)), length(all.times), length(t))

  indexes <- (event.time.mat <= time.mat) & (failures == 1)

  # this returns the max index of the max failure time prior to min(t, all.times[i])
  # here it throw an error if we don't have a single failure!!!
  # should handle that in a more delicate way
  fail.indexes <- apply(indexes,2,cumsum)

  # the + 1 here of the indexes is added if we ave failures before events refer to the logic above
  s <- as.vector(exp(U%*%gamma.vec))*(lambdas.of.fail.times[(fail.indexes + 1)])

  ans <- matrix(as.vector(indexes) - s, length(t), length(all.times), byrow = T)
  return(ans)
}







scoreEval <- function(data, set.U, Cov.e.M.e.D, K = c("linear", "poly", "gaussian"), ...){
  data.M <- data[order(data[,1]),]
  data.D <- data.M[order(data.M[,2]),]

  M.M <- M.vec(Inf, data.M[,1], data.M[,3], c(0,0), as.matrix(data.M[,set.U]))
  M.D <- M.vec(Inf, data.D[,2], data.D[,4],  c(0,0), as.matrix(data.D[,set.U]))

  M.M <- M.M[,order(data.M[,2]),drop = F]

  K.rho <- kernelEval(t(data.D[,c(c(-1:-4), -set.U, -dim(data.D)[2] + 1, -dim(data.D)[2])]), K, ...)

  return ((M.M)%*%K.rho%*%t(M.M) + (M.D)%*%K.rho%*%t(M.D))
}

dM <- function(all.times, failures, gamma.vec, U){
  fail.times.indexes <- which(failures == 1)
  fail.times <- all.times[fail.times.indexes]

  almost.ans <- M.vec(fail.times, all.times, failures, gamma.vec, U)

  ans <- almost.ans - rbind(rep(0, length(all.times)), almost.ans[-dim(almost.ans)[1],])
  return(ans)
}

# THE SAME AS sum.exp.in.risk.set
PI.0 <- function(t, all.times, gamma.vec, U){

  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time.mat <- matrix(rep(t, length(all.times)), length(all.times), length(t), byrow = T)
  event.time.mat <- matrix(rep(all.times, length(t)), length(all.times), length(t))

  indexes <- matrix(rep(1:length(all.times), length(t)), length(all.times), length(t))*(event.time.mat >= time.mat)

  almost.the.right.res <- as.vector(cumsum(exp(as.matrix(U[indexes, ])%*%gamma.vec))[cumsum(colSums(event.time.mat >= time.mat))])
  res <- almost.the.right.res - c(0, almost.the.right.res[-length(almost.the.right.res)])
  # discards the result for the -1
  return(res[-1])
}

# returns a n x
PI.1 <- function(t, all.times, gamma.vec, U, U.cust = U){

  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time.mat <- matrix(rep(t, length(all.times)), length(all.times), length(t), byrow = T)
  event.time.mat <- matrix(rep(all.times, length(t)), length(all.times), length(t))

  indexes <- matrix(rep(1:length(all.times), length(t)), length(all.times), length(t))*(event.time.mat >= time.mat)

  almost.the.right.res <- apply(as.matrix(as.vector(exp(as.matrix(U[indexes, ])%*%gamma.vec))*U.cust[indexes, ], ncol = length(gamma.vec)), 2, cumsum)[cumsum(colSums(event.time.mat >= time.mat)), ]
  almost.the.right.res <- as.matrix(almost.the.right.res, nrow = length(t))
  res <- almost.the.right.res - rbind(0, almost.the.right.res[-dim(almost.the.right.res)[1], , drop = F])
  # discards the result for the -1
  return(matrix(res[-1,], nrow = length(t) - 1))
}

PI.2 <- function(t, all.times, gamma.vec, U){
  # to avoid having to deal with 0 indexes in the initial matrix
  t <- c(-1,t)

  # here we return the sum of the exponentiated gammas mult by the U's of people still in the risk set
  time.mat <- matrix(rep(t, length(all.times)), length(all.times), length(t), byrow = T)
  event.time.mat <- matrix(rep(all.times, length(t)), length(all.times), length(t))

  indexes <- matrix(rep(1:length(all.times), length(t)), length(all.times), length(t))*(event.time.mat >= time.mat)

  # this is a matrix dim(U)[2]^2 x sum(event.time.mat >= time.mat)
  # it contains the matrices of multiplying the transposed rows of U by themselves (cols of t(U))
  # each column of this matrix is the combined rows of a matrices described above
  l <- dim(as.matrix(U))[2]
  simplified.mat <- matrix(U[indexes, rep(1:l, l)]*U[indexes, rep(1:l, rep(l,l))], ncol = l^2)

  almost.the.right.res <- apply(as.vector(exp(as.matrix(U[indexes, ])%*%gamma.vec))*(simplified.mat), 2, cumsum)[cumsum(colSums(event.time.mat >= time.mat)), ]
  almost.the.right.res <- as.matrix(almost.the.right.res, ncol = l^2)
  res <- almost.the.right.res - rbind(0, as.matrix(almost.the.right.res[-dim(almost.the.right.res)[1], , drop = F]))
  # discards the result for the -1
  return(matrix(res[-1,], ncol = l^2))
}

# Ahat computation
# The idea behind this calculation is:
# dEN_i(t) is estimated by the mean number of failure times before i.e. multiply by 1/n (typically dN_i(t) = 1 only for 1 out of n times)
# The second derivative of the log partial likelihood
Ahat <- function(all.times, failures, gamma.vec, U){
  # be careful with the definitions here
  n <- length(all.times)

  fail.times <- all.times[failures == 1]
  #
  pi.0 <- PI.0(fail.times, all.times, gamma.vec, U)/n
  #
  pi.1 <- PI.1(fail.times, all.times, gamma.vec, U)/n

  l <- dim(as.matrix(pi.1))[2]
  pi.1.x2 <- as.matrix(pi.1)[,rep(1:l, l)]*as.matrix(pi.1)[,rep(1:l, rep(l,l))]

  pi.2 <- PI.2(fail.times, all.times, gamma.vec, U)/n

  res <- colSums(pi.0^(-2)*(pi.2*pi.0 - pi.1.x2))

  return(matrix(as.vector(t(res)), ncol = dim(U)[2]))
}

# returns the original gammas!!!
gamma.star <- function(perturb.mat, all.times, failures, gamma.vec, U){
  return(matrix(rep(gamma.vec, dim(perturb.mat)[2]), length(gamma.vec)))
}

## fake lambda pert; W_{gamma,i} = 0
lambda.pert <- function(t, perturb.mat, all.times, failures, gamma.vec, U){
  fail <- which(failures == 1)
  fail.times <- all.times[fail]

  time.mat <- matrix(rep(t, length(fail.times)), length(fail.times), length(t), byrow = T)
  event.time.mat <- matrix(rep(fail.times, length(t)), length(fail.times), length(t))

  fail.times.less.true.false <- (event.time.mat <= time.mat)

  n <- length(all.times)

  dM.mat <- dM(all.times, failures, gamma.vec, U)
  precomp.exp.sum <- PI.0(fail.times, all.times, gamma.vec, U)
  dM.over.precomp.exp.sum <- dM.mat/precomp.exp.sum

  perturbated.mat.across.fail.times <- t(perturb.mat)%*%t(dM.over.precomp.exp.sum)

  fail.times.less.true.false <- (event.time.mat <= time.mat)

  # this is exactly the lambda
  denom.sum.vec <- t(fail.times.less.true.false)%*%(1/precomp.exp.sum)

  # some of the observations are not less than the given time
  # we are also summing all of the other perturbed sums SUM(xi_i dM_i(t)/PI.0(t))
  perturbated.with.only.needed.obs <- t(fail.times.less.true.false)%*%t(perturbated.mat.across.fail.times)

  #res <- matrix(rep(denom.sum.vec, dim(perturb.mat)[2]), length(t)) + perturbated.with.only.needed.obs

  # just repeat the lambda vector
  first.part <- matrix(rep(denom.sum.vec, dim(perturb.mat)[2]), length(t))

  # the sum of SUM(xi_i dM_i(t)/PI.0(t))
  second.part <- perturbated.with.only.needed.obs

  ## THE COMPUTATION OF THE THIRD PART
  res <- exp(log(first.part) + (second.part)/(as.vector(denom.sum.vec)))

  return(res)
}

## the M vector perturbation
M.vec.pert <- function(perturb.mat, t, all.times, failures, gamma.vec, U){
  res <- c()
  # the failure times
  fail.times <- all.times[failures == 1]

  ## comment the last two arguments for log lambda pert
  lambdas.of.fail.times <- rbind(rep(0, dim(perturb.mat)[2]), lambda.pert(fail.times, perturb.mat, all.times, failures, gamma.vec, U))

  time.mat <- matrix(rep(t, length(all.times)), length(all.times), length(t), byrow = T)
  event.time.mat <- matrix(rep(all.times, length(t)), length(all.times), length(t))

  indexes <- (event.time.mat <= time.mat) & (failures == 1)

  fail.indexes <- apply(indexes,2,cumsum)

  # the + 1 here of the indexes is added if we ave failures before events refer to the logic above
  lambdas.of.fail.times.prev <- c(0, lambda(fail.times, all.times, failures, gamma.vec, U))
  prev.s <- as.vector(exp(U%*%gamma.vec))*(lambdas.of.fail.times.prev[(fail.indexes + 1)])

  gamma.star <- gamma.star(perturb.mat, all.times, failures, gamma.vec, U)
  exp.matr <- as.matrix(exp(U%*%gamma.star))
  s <- exp.matr[rep(1:dim(exp.matr)[1], length(t)),]*(lambdas.of.fail.times[(fail.indexes + 1),])

  ans <- (as.vector(indexes) - s)*perturb.mat[rep(1:dim(perturb.mat)[1], length(t)),]
  perturbed.part.of.M <- s

  M.v <- matrix(as.vector(indexes) - prev.s, length(t), length(all.times), byrow = T)

  M.v.pert <- as.vector(t(M.v))*(perturb.mat)[rep(1:dim(perturb.mat)[1], dim(M.v)[1]),]

  ans <- M.v.pert - (perturbed.part.of.M - prev.s)
  return(ans)
}

scoreEval.pert <- function(num.perts, data, set.U, Cov.e.M.e.D, K = c("linear", "poly", "gaussian") , ...){

  Ws <- rnorm(dim(data)[1]*num.perts)
  Ws.mat <- matrix(Ws, dim(data)[1], num.perts)

  data.M <- data[order(data[,1]),]
  data.D <- data.M[order(data.M[,2]),]

  Ws.mat.M <- Ws.mat[order(data[,1]),]
  Ws.mat.D <- Ws.mat[order(data.M[,2]),]

  M.M <- M.vec.pert(Ws.mat.M, Inf, data.M[,1], data.M[,3], c(0,0), as.matrix(data.M[,set.U]))
  M.D <- M.vec.pert(Ws.mat.D, Inf, data.D[,2], data.D[,4], c(0,0), as.matrix(data.D[,set.U]))

  M.M <- M.M[order(data.M[,2]),]


  K.rho <- kernelEval(t(data.D[,c(c(-1:-4), -set.U, -dim(data.D)[2] + 1, -dim(data.D)[2])]), K, ...)

  return(diag(t(M.M)%*%K.rho%*%(M.M) + t(M.D)%*%K.rho%*%(M.D)))
}

scoreEval.min.model <- function(data, set.U, Cov.e.M.e.D, K = c("linear", "poly", "gaussian"), ...){
  data.new <- data[order(data[,dim(data)[2]-1]),]

  M.new <- M.vec(Inf, data.new[,dim(data.new)[2]-1], data.new[,dim(data.new)[2]], c(0,0), as.matrix(data.new[,set.U]))
  K.rho <- kernelEval(t(data.new[,c(c(-1:-4), -set.U, -dim(data.new)[2] + 1, -dim(data.new)[2])]), K, ...)

  return((M.new)%*%K.rho%*%t(M.new))
}

scoreEval.pert.min.model <- function(num.perts, data, set.U, Cov.e.M.e.D, K = c("linear", "poly", "gaussian"), ...){

  Ws <- rnorm(dim(data)[1]*num.perts)
  Ws.mat <- matrix(Ws, dim(data)[1], num.perts)

  data.new <- data[order(data[,dim(data)[2]-1]),]

  Ws.mat.new <- Ws.mat[order(data[,dim(data)[2]-1]),]

  M.new <- M.vec.pert(Ws.mat.new, Inf, data.new[,dim(data.new)[2]-1], data.new[,dim(data.new)[2]], c(0,0), as.matrix(data.new[,set.U]))

  K.rho <- kernelEval(t(data.new[,c(c(-1:-4), -set.U, -dim(data.new)[2] + 1, -dim(data.new)[2])]), K, ...)

  return(diag(t(M.new)%*%K.rho%*%(M.new)))
}

do.simulations <- function(num.sim = 100, with.feat = F, data.size = 100, ncol.gene.mat = 5, feat.m = NA, feat.d  = NA, num.pert = 1000, kernel = c("linear", "poly", "gaussian"), ...){

  res <- c()
  res.new <- c()

  for (i in 1:num.sim){
    data <- generate.dataset.advanced(with.feat = with.feat, data.size = data.size, ncol.gene.mat = ncol.gene.mat, feat.m = feat.m, feat.d  = feat.d)
    num <- scoreEval(data, c(5 + ncol.gene.mat,6 + ncol.gene.mat), matrix(0,data.size, data.size), kernel, ...)
    res <- c(res, sum(scoreEval.pert(1000, data, c(5 + ncol.gene.mat,6 + ncol.gene.mat), matrix(0,data.size, data.size), kernel, ...) >= as.numeric(num))/1000)
    num.new <- scoreEval.min.model(data, c(5 + ncol.gene.mat,6 + ncol.gene.mat), matrix(0,data.size, data.size), kernel, ...)
    res.new <- c(res.new, sum(scoreEval.pert.min.model(1000, data, c(5 + ncol.gene.mat,6 + ncol.gene.mat), matrix(0,data.size, data.size), kernel, ...) >= as.numeric(num.new))/1000)
  }

  print(sum(res <= 0.05)/num.sim)
  print(sum(res.new <= 0.05)/num.sim)

  return(list(res,res.new))
}



## choose rho in the range such that the eigenvalue decay rate is between 1.1 and 2.2 ##
## 1.1 to 4
chooseRhoInt.opt.new <- function(Z, rho, rho0, l0 = NA, kernel = c("gaussian", "poly"), d = NA, rate.range=c(1.5,4)){
  n <- nrow(Z)
  if (kernel == "gaussian"){
    G <- gaussKernelEval.multipleSigmas(Z, sigma = rho)
  }
  if (kernel == "poly"){
    G <- polyKernelEval.multipleSigmas(Z, a = rho, d = d)
  }

  slope.eigen <- function(ind, G. = G, Z. = Z, rho0. = rho0){
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

  slope.vec <-  na.omit(cbind("rho"=rho,"slope"=t(sapply(1:length(rho),slope.eigen))))
  colnames(slope.vec)[2:3] <-  c("slope","m.keep")
  print(slope.vec)

  rho <- slope.vec[,"rho"] ## changed
  res_rho.opt <- c(max(c(min(rho),slope.vec[slope.vec[,"slope"]<=rate.range[1],"rho"])),min(c(max(rho), slope.vec[slope.vec[,"slope"]>=rate.range[2],"rho"]))) ## changed
  return(res_rho.opt)
}


VTM <- function(vc, dm){matrix(vc, ncol=length(vc), nrow=dm, byrow=T)}

sup.stat.both_boris <- function(num.perts, Ws=NULL, data, set.U, Cov.e.M.e.D, rho = 1:40, l0 = NA, kernel = c("gaussian", "poly"), d = NA, est.gamma=F, pca.thres=NULL){
  if(is.null(Ws[1])){
    Ws <- rnorm(dim(data)[1]*num.perts)
  }
  n <- nrow(data)
  data.M <- data[order(data[,1]),]
  data.D <- data.M[order(data.M[,2]),]
  p.gene <- ncol(data) - 6 # TODO is that correct ?
  ind.gene = 1:p.gene + 4 ## data.M sorted by X_M (col #1), data.D sorted by X_D (col #2)
  if(est.gamma){
    gamhat.M = coxph(Surv(data.M[,1],1*(data.M[,3]==1))~as.matrix(data.M[,set.U]))$coef
    gamhat.D = coxph(Surv(data.D[,2],1*(data.D[,4]==1))~as.matrix(data.D[,set.U]))$coef
  }else{
    gamhat.M = gamhat.D = rep(0,length(set.U))
  }
  M.Mc <- M.vec(Inf, data.M[,1], 1*(data.M[,3]==1), gamhat.M, as.matrix(data.M[,set.U])) ## cause specific hazard for M
  M.Dc <- M.vec(Inf, data.M[,1], 1*(data.M[,3]==2), gamhat.M, as.matrix(data.M[,set.U])) ## cause specific hazard for D
  M.Dm <- M.vec(Inf, data.D[,2], 1*(data.D[,4]==1), gamhat.D, as.matrix(data.D[,set.U])) ## marginal hazard for D
  M.Mc <- M.Mc[, order(data.M[,2]), drop = F] ## this ensures that sorting is done according to sorting of data.D (col #2) ##
  M.Dc <- M.Dc[, order(data.M[,2]), drop = F] ## this ensures that sorting is done according to sorting of data.D (col #2) ##

  if(kernel=="linear"){K.rho=kernelEval(t(data.D[,ind.gene]), K = "linear"); K.rho = matrix(K.rho,nrow=1) }
  else if(kernel=="gaussian"){K.rho=gaussKernelEval.multipleSigmas(t(data.D[,ind.gene]), sigma=rho)}
  else if(kernel=="poly"){K.rho=polyKernelEval.multipleSigmas( t(data.D[,ind.gene]), a=rho, d=d)}
  sqrt.ncol.K.rho <- sqrt(ncol(K.rho))
  K.rho <- matrix(c(t(K.rho)), ncol = sqrt.ncol.K.rho, byrow = T)
  #dim(K.rho)

  # do PCA
  if(is.null(pca.thres)){
    warning("Not performing PCA: potential loss of power, especially on finite samples")
  }else{
    if(kernel=="linear"){
      K.rho_eig <- eigen(K.rho)
      ncomp <- which(cumsum(K.rho_eig$values/sum(K.rho_eig$values))>=pca.thres)[1]
      K.rho <- tcrossprod(matrix(rep(sqrt(K.rho_eig$values[1:ncomp]), sqrt.ncol.K.rho), nrow = sqrt.ncol.K.rho, byrow = T)*K.rho_eig$vectors[, 1:ncomp])
    }else{
      K.rho_new <- NULL
      for (i in 1:length(rho)){
        K.rho_eig <- eigen(K.rho[1:sqrt.ncol.K.rho+(i-1)*sqrt.ncol.K.rho, ])
        ncomp <- which(cumsum(K.rho_eig$values/sum(K.rho_eig$values))>=pca.thres)[1]
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
  M.Mc.pert <- M.vec.pert(Ws.mat.M, Inf, data.M[,1], 1*(data.M[,3]==1), gamhat.M, as.matrix(data.M[,set.U])); ## cause specific for M
  M.Dc.pert <- M.vec.pert(Ws.mat.M, Inf, data.M[,1], 1*(data.M[,3]==2), gamhat.M, as.matrix(data.M[,set.U])); ## cause specific for D
  M.Mc.pert <- M.Mc.pert[order(data.M[,2]),]; ## this ensures that sorting is done according to sorting of data.D
  M.Dc.pert <- M.Dc.pert[order(data.M[,2]),]; ## this ensures that sorting is done according to sorting of data.D
  M.Dm.pert <- M.vec.pert(Ws.mat.D, Inf, data.D[,2], 1*(data.D[,4]==1), gamhat.D, as.matrix(data.D[,set.U])); ## marginal for D
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
  if(est.gamma){gamhat.min = coxph(Surv(data.min[,ncol(data.min)-1],1*(data.min[,ncol(data.min)]))~as.matrix(data.min[,set.U]))$coef}else{gamhat.min = rep(0,length(set.U))}
  Ws.mat.min <- Ws.mat[order(data[,dim(data)[2]-1]),]
  M.min <- M.vec(Inf, data.min[,dim(data.min)[2]-1], data.min[,dim(data.min)[2]], gamhat.min, as.matrix(data.min[,set.U]))

  if(kernel=="linear"){K.rho=kernelEval(t(data.min[,ind.gene]), K = "linear"); K.rho = matrix(K.rho,nrow=1) }
  if(kernel=="gaussian"){K.rho=gaussKernelEval.multipleSigmas(t(data.min[,ind.gene]),sigma=rho)}
  if(kernel=="poly"){    K.rho=polyKernelEval.multipleSigmas( t(data.min[,ind.gene]),a=rho,d=d)}
  sqrt.ncol.K.rho <- sqrt(ncol(K.rho))
  K.rho <- matrix(c(t(K.rho)), ncol = sqrt.ncol.K.rho, byrow = T)
  #dim(K.rho)

  # do PCA
  if(kernel=="linear"){
    K.rho_eig <- eigen(K.rho)
    ncomp <- which(cumsum(K.rho_eig$values/sum(K.rho_eig$values))>=pca.thres)[1]
    K.rho <- tcrossprod(matrix(rep(sqrt(K.rho_eig$values[1:ncomp]), sqrt.ncol.K.rho), nrow = sqrt.ncol.K.rho, byrow = T)*K.rho_eig$vectors[, 1:ncomp])
  }else{
    K.rho_new <- NULL
    for (i in 1:length(rho)){
      K.rho_eig <- eigen(K.rho[1:sqrt.ncol.K.rho+(i-1)*sqrt.ncol.K.rho, ])
      ncomp <- which(cumsum(K.rho_eig$values/sum(K.rho_eig$values))>=pca.thres)[1]
      K.rho_new <- rbind(K.rho_new,
                         tcrossprod(matrix(rep(sqrt(K.rho_eig$values[1:ncomp]), sqrt.ncol.K.rho), nrow = sqrt.ncol.K.rho, byrow = T)*K.rho_eig$vectors[, 1:ncomp])
      )
    }
    K.rho <- K.rho_new
  }

  stats.min <- obtain.K.rhos(K.rho, M.min)
  M.min.pert <- M.vec.pert(Ws.mat.min, Inf, data.min[,dim(data.min)[2]-1], data.min[,dim(data.min)[2]], gamhat.min, as.matrix(data.min[,set.U]))
  part.1.min <- t(M.min.pert)%*%t(K.rho); part.2.min <- c(t(M.min.pert))*c(part.1.min)
  part.3.min <- matrix(part.2.min, ncol = num.perts, byrow = T)
  res.min <- c(); for(i in 0:(length(rho) - 1)){res.min <- cbind(res.min,colSums(part.3.min[i*n + 1:n,]))}
  st.devs.min <- sqrt(colSums((res.min - matrix(rep(colMeans(res.min),num.perts),nrow = nrow(res.min),byrow=T))^2)/(num.perts - 1))
  pval.min = sum(apply(t(res.min)/st.devs.min,2,max) >= max(stats.min/st.devs.min))/num.perts

  if(est.gamma){
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