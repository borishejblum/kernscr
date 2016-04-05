# #R setup
# rm(list=ls())
# library(mvtnorm)
# library(MASS)
# library(reshape2)
# source("~/Dropbox/HSPH/Kernel Score Test-Matey/Code/BorisCode/FUN_kernel_Boris.R")
# setwd("~/Dropbox/HSPH/Kernel Score Test-Matey/Code/BorisCode/BreastCancer/")
#
# #loading data
# load("c2.cp.v2.5.GS.RData")
# load("BreastCancer-Vijver.RData")
# load("BreastCancer-Vijver-Clinical.RData")
# nm.pathway_BC <- as.matrix(read.csv("BC-pathway-list.csv", header=F))
# nm.pathway_CaiToniniLin <- as.matrix(read.table("Cai_Tonini_Lin_List.txt", header=FALSE))
#
# #processing/preparing data
# BC.dat.clin <- get("clin.data")
# BC.dat.exp <- get("llid_norm.mat")
# PCA.percent <- 0.9
# pvalue.all.Lin <- NULL
# pvalue.all.Gau <- NULL
# sup.hat.Gau <- NULL
# sup.hat.Lin <-  NULL
# sup.ptb.Gau <- list(NULL,NULL,NULL)
# sup.ptb.Lin <-  list(NULL,NULL)
#
#
#
# #Setting up the simulation settings
# nm.pathway <- nm.pathway_CaiToniniLin
# kern2test <- c("linear", "gaussian")#("linear", "poly", "gaussian")
#
#
# #Description of the data
# library(survival)
# fit <- coxph(Surv(tsurv, 1-evdeath) ~ 1, data=BC.dat.clin)
# dim(BC.dat.clin)
# survfit(fit)
# apply(BC.dat.clin[, -c(1, ncol(BC.dat.clin))], 2, mean, na.rm=TRUE)
# length(nm.pathway)
# summary(sapply(GS.llids[nm.pathway], function(v){length(intersect(rownames(BC.dat.exp), v))}))
# tt <- table(BC.dat.clin[,c("evdeath", "evmeta")])/nrow(BC.dat.clin)
# round(tt, digits=2)
# sum(round(tt, digits=2))
#
#
# #
# set.seed(1202)
# B0 <- 5000
# num.perts <- B0
# Gi.mat0 <- matrix(rnorm(B0*nrow(BC.dat.clin)), nrow=nrow(BC.dat.clin))
# obs2keep <- rep(TRUE, 260)
#
#
# # run the actual test ----
# for (mykern in kern2test){
#   for(ip in 1:length(nm.pathway)){
#     tmpind <- match(nm.pathway[ip],names(GS.llids))
#     gene.llids <- rownames(BC.dat.exp)
#     pathway.llids <- GS.llids[[tmpind]]
#     pathway.llids <- intersect(gene.llids, pathway.llids)
#     pathway.expression <- t(BC.dat.exp[pathway.llids, ])
#     BC.dat <- data.frame("Tmin1"=BC.dat.clin$tmeta,
#                          "Tmin2"=BC.dat.clin$tsurv,
#                          "delta1"=(BC.dat.clin$evmeta == 0)*(BC.dat.clin$evdeath == 1)*2 + BC.dat.clin$evmeta,
#                          "delta2"=BC.dat.clin$evdeath)
#     BC.dat <- data.frame(BC.dat,
#                          pathway.expression[match(BC.dat.clin$ID,row.names(pathway.expression)),])
#     BC.dat <- data.frame(BC.dat,
#                          BC.dat$Tmin1, (BC.dat$delta1 != 0))
#     mydata <- BC.dat
#     na.ind <- apply(is.na(mydata),1,sum)>0
#     obs2keep <- (obs2keep & !na.ind)
#     mydata <- mydata[!na.ind,]
#     #Gi.mat <- Gi.mat0[!na.ind,]
#     #nn <- nrow(mydata)
#
#     p.gene <- ncol(mydata) - 6
#
#     if(mykern != "linear"){
#       #print(apply(mydata[,c("delta1","delta2")]==1, 2, mean))
#
#       #Jen changes:
#       #rho.max <- round(p.gene*2*sum(apply(mydata[,1:p.gene+4], 2, sd))) #old version before Jen change
#       rho.max <- 200 #new version after Jen change
#       #rho.set <- chooseRhoInt.opt.new(t(mydata[,1:p.gene+4]), rho=seq(1, rho.max, by=2), rho0=.95, l0=NA, kernel=mykern, d=2) #old version before Jen change
#       rho.set_temp <- chooseRhoInt.opt.new(t(mydata[,1:p.gene+4]), rho=seq(1, rho.max, by=2), rho0=PCA.percent, l0=NA, kernel=mykern, d=2) #new version after Jen change
#       #e <- try(rho.set <- seq(rho.set[1], 1, length=20)) #old version before Jen change
#       e <- try(rho.set <-  exp(seq(log(rho.set_temp[1]),log(rho.set_temp[2]), length=30))) #new version after Jen change
#       if(class(e)=="try-error"){
#         rho.set = 1
#       }
#       #print(rho.set)
#     }else{
#       rho.set = 1
#     }
#
#     tmpout <- sup.stat.both_boris(num.perts=B0, data=mydata, set.U=(5:(ncol(mydata)-2)),
#                                   Cov.e.M.e.D=0, rho=rho.set, l0=NA, kernel=mykern, d=2,
#                                   est.gamma=FALSE, pca.thres=PCA.percent)
#
#     pvals_out <- c(tmpout[["raw_pvals"]], mykern, ip, nm.pathway[ip])
#     write(pvals_out, file="raw_pvals_CaiToniniLin.txt", append=TRUE, ncol=500)
#
#     pertb_nullpvals_out <- cbind(round(tmpout[["null_pvals_pertbs"]],16), rep(mykern, B0), rep(ip, B0), rep(nm.pathway[ip], B0))
#     write.table(pertb_nullpvals_out, file="pertb_nullpvals_CaiToniniLin_newrho.txt", append=TRUE,
#                 row.names = FALSE, col.names=FALSE, quote=TRUE)
#
#     cat(mykern, " ", ip, "/", length(nm.pathway), ":\n", paste(pvals_out, collapse = " "), "\n\n", sep="")
#   }
# }
#
#
#
# res_raw <- read.table("raw_pvals_CaiToniniLin.txt", quote="\"", comment.char="")
# colnames(res_raw) <- c("sum3", "McDc", "McDm", "max3", "max.McDc", "Mc", "Dc", "Dm", "PFS", "kern", "pathway_num", "pathway_name")
#
# nullpvals_ptb <- read.table("pertb_nullpvals_CaiToniniLin.txt")
# colnames(nullpvals_ptb) <- c("sum3", "McDc", "McDm", "max3", "max.McDc", "Mc", "Dc", "Dm", "PFS", "kern", "pathway_num", "pathway_name")
#
# res_adj <- res_raw
# for(i in 1:(ncol(res_raw)-3)){
#   for (mykern in kern2test){
#     null_pv_temp <- nullpvals_ptb[(nullpvals_ptb$kern==mykern),c(i, ncol(res_raw)-2:0)]
#     colnames(null_pv_temp)[1] <- "pval"
#     null_pv_temp <- cbind.data.frame(null_pv_temp, "ptb_num"=paste0("ptb",rep(1:5000, length(nm.pathway))))
#     null_pv_temp <- dcast(null_pv_temp, pathway_num + pathway_name + kern ~ptb_num, value.var="pval")
#     null_pv_temp.min <- apply(1-null_pv_temp[,-c(1:3)], 2, FUN=min)
#     res_adj[res_adj$kern==mykern,i] <- sapply(res_raw[res_raw$kern==mykern,i], function(p){mean(p > null_pv_temp.min)})
#   }
#   cat(i, "/", ncol(res_raw)-3, "\n", sep="")
# }
#
# res_adj_melted <- reshape2::melt(res_adj, id.vars = c("kern", "pathway_num", "pathway_name"))
# res_raw_melted <- reshape2::melt(res_raw, id.vars = c("kern", "pathway_num", "pathway_name"))
#
#
#
#
#
#
#
#
#
#
# pdf(file=paste("Dm_Gauss_kernPCA",".pdf",sep=""),width=5,height=10)
# plot_res_2011(raw_melted= res_raw_melted, adj_melted = res_adj_melted, kernel="gaussian", method="Dm")
# dev.off()
#
# pdf(file=paste("McDm_PFS_GausskernPCA",".pdf",sep=""),width=9,height=10)
# plot_res_2011(raw_melted= res_raw_melted, adj_melted = res_adj_melted, kernel="gaussian", method=c("PFS", "McDm"), title="")
# dev.off()
#
# pdf(file=paste("SCR_PFS_GausskernPCA_2015",".pdf",sep=""),width=9,height=10)
# plot_res_2015(raw_melted= res_raw_melted, adj_melted = res_adj_melted, kernel="gaussian", method=c("McDm", "PFS"), title="")
# dev.off()
#
# adj_temp <- cbind("Dm"=res_adj_melted[res_adj_melted$kern=="gaussian" & res_adj_melted$variable=="Dm", "value"],
#                   "PFS"=res_adj_melted[res_adj_melted$kern=="gaussian" & res_adj_melted$variable=="PFS", "value"],
#                   "McDm"=res_adj_melted[res_adj_melted$kern=="gaussian" & res_adj_melted$variable=="McDm", "value"])
# colSums(adj_temp<=0.05)
#
# raw_temp <- cbind("Dm"=res_raw_melted[res_raw_melted$kern=="gaussian" & res_raw_melted$variable=="Dm", "value"],
#                   "PFS"=res_raw_melted[res_raw_melted$kern=="gaussian" & res_raw_melted$variable=="PFS", "value"],
#                   "McDm"=res_raw_melted[res_raw_melted$kern=="gaussian" & res_raw_melted$variable=="McDm", "value"])
# colSums(raw_temp<0.05)
# sum(raw_temp[, "McDm"] > raw_temp[, "PFS"])
# sum(adj_temp[, "McDm"] > adj_temp[, "PFS"])
#
# res_raw_melted[res_raw_melted$kern=="gaussian" & res_raw_melted$variable=="Dm" & res_raw_melted$pathway_name=="AKAP13PATHWAY", "value"]
