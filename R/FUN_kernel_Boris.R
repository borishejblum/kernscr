
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








