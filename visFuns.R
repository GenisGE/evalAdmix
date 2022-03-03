plotCorRes <- function(cor_mat, pop=NULL, ord=NULL, superpop=NULL,
                       title="Correlation of residuals", min_z=NA,max_z=NA, 
                       cex.main=1.5, cex.lab=1.5, cex.legend=1.5, color_palette=c("#001260", "#EAEDE9", "#601200"),
                       pop_labels = c(T,T), plot_legend = T, adjlab = 0.1, rotatelabpop=0, rotatelabsuperpop=0,lineswidth=1, lineswidthsuperpop=2,
                       adjlabsuperpop=0.16,cex.lab.2 = 1.5){
  
  op <- par(mfrow=c(1,1) ,mar=c(5,4,4,2) +0.1,xpd=F, oma=c(0,0,0,0))
  on.exit(par(op))
  
  N <- dim(cor_mat)[1]
  

  if(is.null(ord)&!is.null(pop)) ord <- order(pop)
  if(is.null(ord)&is.null(pop)) ord <- 1:nrow(cor_mat)

  if(is.null(pop)){
      pop <- rep(" ", nrow(cor_mat))
      lineswidth <- 0
  }
    
  pop<-pop[ord]
  
  N_pop <- vapply(unique(pop[ord]), function(x) sum(pop==x),1)
  
  cor_mat <- cor_mat[ord,ord]
  
  ## Set lower part of matrix as population mean correlation
  mean_cors <- matrix(ncol=length(unique(pop)), nrow=length(unique(pop)))
  colnames(mean_cors) <- unique(pop)
  rownames(mean_cors) <- unique(pop)
  
  for(i1 in 1:(length(unique(pop)))){
    for(i2 in 1:(length(unique(pop)))){
      p1 <- unique(pop)[i1]
      p2 <- unique(pop)[i2]
      mean_cors[i1,i2]<- mean(cor_mat[which(pop==p1),
                                      which(pop==p2)][!is.na(cor_mat[which(pop==p1),
                                                                     which(pop==p2)])])
      
    }
  }
  
  for(i1 in 1:(N-1)){
    for(i2 in (i1+1):N){
      cor_mat[i1, i2] <- mean_cors[pop[i2], pop[i1]]
      
    }
  }
  
  z_lims <- c(min_z, max_z)
  
    if(all(is.na(z_lims))) z_lims <- c(-max(abs(cor_mat[!is.na(cor_mat)])),
                                         max(abs(cor_mat[!is.na(cor_mat)])))
  #if(all(is.null(z_lims))) max_z <- max(abs(cor_mat[!is.na(cor_mat)]))

    
  if(any(is.na(z_lims))) z_lims <- c(-z_lims[!is.na(z_lims)], z_lims[!is.na(z_lims)])
  #if(any(is.null(z_lims))) max_z <- z_lims[!is.null(z_lims)]

    min_z <- z_lims[1]
    max_z <- z_lims[2]
    
  diag(cor_mat) <- 10
  nHalf <- 10
  
  # make sure col palette is centered on 0
  Min <- min_z
  Max <- max_z
  Thresh <- 0
  
  ## Make vector of colors for values below threshold
  rc1 <- colorRampPalette(colors = color_palette[1:2], space="Lab")(nHalf)    
  ## Make vector of colors for values above threshold
  rc2 <- colorRampPalette(colors = color_palette[2:3], space="Lab")(nHalf)
  rampcols <- c(rc1, rc2)
  
  rampcols[c(nHalf, nHalf+1)] <- rgb(t(col2rgb(color_palette[2])), maxColorValue=256) 
  
  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)
  
  rlegend <- as.raster(matrix(rampcols, ncol=1)[length(rampcols):1,])
  if(plot_legend){
    layout(matrix(1:2,ncol=2), width = c(4,1),height = c(1,1))
    par(mar=c(5,4,4,0),oma=c(1,4.5,2,0))
  }else
    par(mar=c(5,4,4,5),oma=c(1,4.5,2,0))
  image(t(cor_mat), col=rampcols, breaks=rampbreaks,
        yaxt="n",xaxt="n", zlim=c(min_z,max_z),useRaster=T,
        main=title, 
        oldstyle=T,cex.main=cex.main,xpd=NA)
  image(ifelse(t(cor_mat>max_z),1,NA),col="darkred",add=T)
  if(min(cor_mat)<min_z) image(ifelse(t(cor_mat<min_z),1,NA),col="darkslateblue",add=T)
  image(ifelse(t(cor_mat==10),1,NA),col="black",add=T)
  
  # put pop info
  if(pop_labels[2])
    text(sort(tapply(1:length(pop),pop,mean)/length(pop)),-adjlab,unique(pop),xpd=NA,cex=cex.lab, srt=rotatelabpop)
  if(pop_labels[1])
    text(-adjlab,sort(tapply(1:length(pop),pop,mean)/length(pop)),unique(pop),xpd=NA, cex=cex.lab,srt=90-rotatelabpop)
  abline(v=grconvertX(cumsum(sapply(unique(pop),function(x){sum(pop==x)}))/N,"npc","user"),
         col=1,lwd=lineswidth,xpd=F)
  abline(h=grconvertY(cumsum(sapply(unique(pop),function(x){sum(pop==x)}))/N, "npc", "user"),
         col=1,lwd=lineswidth,xpd=F)
  
  # put superpop if not null
    if(!is.null(superpop)){
        superpop <- superpop[ord]
    if(pop_labels[2])
      text(sort(tapply(1:length(superpop),superpop,mean)/length(superpop)),-adjlabsuperpop,unique(superpop),xpd=NA,cex=cex.lab.2, srt=rotatelabsuperpop, font=2)
    if(pop_labels[1])
      text(-adjlabsuperpop,sort(tapply(1:length(superpop),superpop,mean)/length(superpop)),unique(superpop),xpd=NA, cex=cex.lab.2,srt=90-rotatelabsuperpop,font=2)
    abline(v=grconvertX(cumsum(sapply(unique(superpop),function(x){sum(superpop==x)}))/N,"npc","user"),
           col=1,lwd=lineswidthsuperpop,xpd=F)
    abline(h=grconvertY(cumsum(sapply(unique(superpop),function(x){sum(superpop==x)}))/N, "npc", "user"),
           col=1,lwd=lineswidthsuperpop,xpd=F)
  }
  
  if(plot_legend){
    par(mar=c(5,0.5,4,2))
    plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')    
    
    rasterImage(rlegend, 0, 0.25, 0.4,0.75)
    text(x=0.8, y = c(0.25,0.5, 0.75),
         labels = c(-max(abs(min_z),abs(max_z)), 0, max(abs(min_z),abs(max_z))),
         cex=cex.legend,xpd=NA)
  }
}


orderInds <- function(q=NULL, pop=NULL, popord=NULL){
  # Function to order individuals for admixture and evalAdmix plots. 
  # recommended is to use pop, then if q is given it will order within pop by admixture proporiton. poporder allows to pre-specify order of populations
  # if only q is given will group individuals by main cluster they are assigned
  
  ordpop <- function(x, pop, q){
    idx <- which(pop==x)
    main_k <- which.max(apply(as.matrix(q[idx,]),2,mean))
    ord <- order(q[idx,main_k])
    idx[ord]
  } 
  
  if(!is.null(pop)){
    
    if(is.null(popord)) popord <- unique(pop)
    
    if(!is.null(q)){ 
      
      ord <- unlist(sapply(popord, ordpop, pop=pop, q=q))
      
    } else if (is.null(q)) {
      
      ord <- unlist(sapply(popord, function(x) which(pop==x)))
      
    }
  } else if (is.null(pop)&!is.null(q)) {
    
    # get index of k with max value per individual
    main_k <- apply(q,1, which.max)
    
    # get max q per indivdiual
    main_q <- q[cbind(1:nrow(q),main_k)]
    
    ord <- order(main_k, main_q)
    
  } else {stop("Need at least an argument to order.")}

  return(ord)
  
}


orderK <- function(q, refinds= NULL,refpops = NULL, pop=NULL){
  # Function to order ancestral populations, useful to keep cluster colors in admix plot the same when comparing results across different k values
  # if you give refinds will use maximum Q value of each individual to define clusters
  # if you give refpops (must also give pops) will use maximum mean admixture proportions within inds from pop to define clusters
  # if any refpops or refinds have same cluster as maximum, the admixture plot will look really bad (you will lose a cluster and another will be twice)
  
  k <- ncol(q)
  kord <- integer(0)
  
  if(is.null(refinds)){
  refpops <- refpops[1:k]
  
  for(p in refpops){
    
    kord <- c(kord, which.max(apply(q[pop==p,],2,mean)))
    
  }
  } else {
    
    refinds <- refinds[1:k]
    
    for(i in refinds){
      
      kord <- c(kord, which.max(q[i,]))
    }
  }
  
    # if(any(rowSums(q[,kord]!=1))) warning("reordered admixture proportions don't sum to 1, make sure every refind or refpop defines a unique cluster.")

    return(kord)
}


plotAdmix <- function(q, pop=NULL, ord=NULL, inds=NULL,
                      colorpal= c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"),
                      main=paste("Admixture proportions assuming K =",k),
                      cex.main=1.5, cex.lab=1, rotatelab=0,padj=0, cex.inds=1,
                      drawindslines=TRUE){
  # simple function to plot admixture proprotions, just to make sure the ordering of individuals is handled as in plotCorRes.
  
  k <- ncol(q)
  
  if(k>length(colorpal))
    warning("not enought colors for all Ks in palette.")
  
  # if(!is.null(ord)) if(!ord) ord <- 1:nrow(q)
  
  if(is.null(ord)&!is.null(pop)) ord <- order(pop)
  if(is.null(ord)&is.null(pop)) ord <- 1:nrow(q)
  
  barplot(t(q)[,ord], col=colorpal, space=0, border=NA, cex.axis=1.2,cex.lab=1.8,
          ylab="Admixture proportions", xlab="", main=main, cex.main=cex.main,xpd=NA)
  
  if(!is.null(inds)){
    text(x = 1:nrow(q) - 0.5,-0.1, inds[ord],xpd=NA,srt=90, cex=cex.inds)
  }
  
  if(!is.null(pop)){
    
    text(sort(tapply(1:length(pop),pop[ord],mean)),-0.05-padj,unique(pop[ord]),xpd=NA, srt=rotatelab, cex=cex.lab)
    if(drawindslines) abline(v=1:nrow(q), col="white", lwd=0.2)
    abline(v=cumsum(sapply(unique(pop[ord]),function(x){sum(pop[ord]==x)})),col=1,lwd=1.2)
    
  }
  
}
