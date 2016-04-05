plot_res_2011 <- function(raw_melted, adj_melted, kernel, method,
                          pathway_names=TRUE, title=NULL, lower_threshold=-4){
  range <- seq(70,1)
  if(pathway_names){
    nplot <- length(kernel)+1
  }else{
    nplot= length(kernel)
  }

  if(is.null(title)){
    title <- method
  }

  layout(matrix(1:nplot, nrow=1), c(3, rep(2, nplot-1)), rep(1, nplot))
  if(title==""){
    par(mar=c(5,0,1,0))
  }else{
    par(mar=c(5,0,5,0))
  }

  dolegend <- TRUE

  for(k in kernel){

    tmpy <- log10(adj_melted[adj_melted$kern==k & adj_melted$variable==method, "value"])
    tmpy <- pmax(tmpy, lower_threshold)

    if(pathway_names){
      plot(tmpy,range,type="n",pch=15,ylab="",xlab="",
           yaxt="n",xaxt="n",xlim=c(-3,0),bty="n",cex=0.7)
      txt.pathway <- as.character(adj_melted[adj_melted$kern==k & adj_melted$variable==method, "pathway_name"])
      #Cai 2011 paper order:
      #ordralpha <- order(txt.pathway)
      #txt.pathway <- txt.pathway[ordralpha]
      #ordr <- rev(order(sapply(txt.pathway, nchar)))
      #txt.pathway <- txt.pathway[ordr]
      #txt.pathway <- txt.pathway[c(1:13, 15:16, 14, 17:19, 21:22, 20, 23:24, 29:31, 25:28, 34:35, 32:33, 36:70)]
      for(jj in range){
        text(0,jj, txt.pathway[70 - jj + 1], srt = 0, xpd = TRUE, pos = 2, cex=0.57)
      }
      pathway_names=FALSE
    }

    if(dolegend){
      legend("bottom", pch=c(4,15), c("raw p-val.", "adj. p-val."), cex=0.85)
      dolegend <- FALSE
    }

    if(title==""){
      par(mar=c(5,0,1,4))
    }else{
      par(mar=c(5,0,5,4))
    }
    plot(tmpy,range,pch=15,ylab="",xlab="log10 Pvalue",
         yaxt="n",xaxt="n",xlim=c(lower_threshold,0),bty="n",cex=0.)
    axis(side=1, at = c(lower_threshold, -3, -2, -1, 0),
         labels=c(paste0("< ", lower_threshold),"-3", "-2","-1","0"))
    tmpy = pmax(log10(raw_melted[raw_melted$kern==k & raw_melted$variable==method, "value"]),
                lower_threshold)
    points(tmpy,range,pch=4,cex=0.7)
    tmpy = pmax(log10(adj_melted[adj_melted$kern==k & raw_melted$variable==method, "value"]),
                lower_threshold)
    points(tmpy,range,pch=15,cex=0.62)
    abline(v=log10(0.05),col="gray")
    if(title!=""){
      mtext(side=3, text=paste(paste0(toupper(substring(k,1,1)), substring(k,2)), "kernel"),
            cex=0.95,line=0.7, adj=0.6)
    }
  }
  if(title!=""){
    mtext(side=3,text=title, cex=1,line=3, at=-5)
  }
  #mtext(side=1,text="log10 Pvalue",cex=0.9,line=3.5, adj=0.5, at=-5)
}