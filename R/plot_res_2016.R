plot_res_2016 <- function(raw_melted, adj_melted, kernel, method, pathway_names=TRUE, title=NULL){
  range <- seq(70,1)

  if(pathway_names){
    nplot <- length(method) + 1
  }else{
    nplot <- length(method)
  }

  tmpy <- log10(raw_melted[raw_melted$kern==kernel & raw_melted$variable==method[1], "value"])
  tmpy = pmax(tmpy, -4)


  layout(matrix(1:nplot, nrow=1), c(2.5, rep(2, nplot-1)), rep(1,nplot))

  par(mar=c(5,0,1,0))
  if(pathway_names){
    plot(tmpy,range,type="n",pch=15,ylab="",xlab="",
         yaxt="n",xaxt="n",xlim=c(-3,0),bty="n",cex=0.7)
    txt.pathway <- adj_melted[adj_melted$kern==kernel & adj_melted$variable==method[1], "pathway_name"]
    for(jj in range){
      text(0,jj, txt.pathway[70 - jj + 1], srt = 0, xpd = TRUE, pos = 2, cex=0.8)
    }
    method_legend <- method
    if(sum(method_legend=="McDm")>0){
      method_legend[method_legend=="McDm"] <- "SCR"
    }
    legend("bottom", pch=c(4,15), rev(method_legend))
  }

  par(mar=c(5,0,1,4))
  plot(tmpy,range,pch=15,ylab="",xlab="log10 raw p-value",
       yaxt="n",xaxt="n",xlim=c(-4,0),bty="n",cex=0.7, cex.lab=1.1)
  axis(side=1, at = c(-4, -3, -2, -1, 0),labels=c("< -4","-3", "-2","-1","0"))

  tmpy = pmax(log10(raw_melted[raw_melted$kern==kernel & raw_melted$variable==method[2], "value"]),-4)
  points(tmpy,range,pch=4,cex=0.8)
  abline(v=log10(0.05),col="gray");
  #mtext(side=3,text="Raw p-values", cex=0.95,line=0.7, adj=0.6)

  tmpy <- log10(adj_melted[adj_melted$kern==kernel & adj_melted$variable==method[1], "value"])
  tmpy = pmax(tmpy, -4)
  plot(tmpy,range,pch=15,ylab="",xlab="log10 adjusted p-value",
       yaxt="n",xaxt="n",xlim=c(-4,0),bty="n",cex=0.7, cex.lab=1.1)
  axis(side=1, at = c(-4, -3, -2, -1, 0),labels=c("< -4","-3", "-2","-1","0"))

  tmpy = pmax(log10(adj_melted[adj_melted$kern==kernel & adj_melted$variable==method[2], "value"]),-4)
  points(tmpy,range,pch=4,cex=0.8)
  abline(v=log10(0.05),col="gray")
  #mtext(side=3,text="Ajusted p-values", cex=0.95,line=0.7, adj=0.6)

}