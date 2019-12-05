plotCorRes <- function(cor_mat, pop, 
                       title="Correlation of residuals", min_z=NULL,max_z=NULL, 
                       is.ord=F, cex.main=1.5, cex.lab=1.5, cex.legend=1.5, color_palette=c("blue", "green", "red"),
                       pop_labels = c(T,T), plot_legend = T){
    
    
    N <- dim(cor_mat)[1]
    
    if(!is.ord){
        ord <- order(pop)
    }else{
        ord <- 1:N
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
            cor_mat[i1, i2] <- mean_cors[pop[i1], pop[i2]]
            
        }
    }
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
    
    
    if(is.null(min_z)) min_z <- min(cor_mat[!is.na(cor_mat)])
    if(is.null(max_z)) max_z <- max(cor_mat[!is.na(cor_mat)])
    
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

    if(pop_labels[2])
        text(sort(tapply(1:length(pop),pop,mean)/length(pop)),-0.1,unique(pop),xpd=NA,cex=cex.lab)
    if(pop_labels[1])
        text(-0.1,sort(tapply(1:length(pop),pop,mean)/length(pop)),unique(pop),xpd=NA, cex=cex.lab,srt=90)
    abline(v=grconvertX(cumsum(sapply(unique(pop),function(x){sum(pop==x)}))/N,"npc","user"),
           col=1,lwd=1.2,xpd=F)
    abline(h=grconvertY(cumsum(sapply(unique(pop),function(x){sum(pop==x)}))/N, "npc", "user"),
           col=1,lwd=1.2,xpd=F)
     if(plot_legend){
         par(mar=c(5,0.5,4,2))
         plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')    
         
         rasterImage(rlegend, 0, 0.25, 0.4,0.75)
         text(x=0.8, y = c(0.25,0.5, 0.75),
              labels = c(-max(abs(min_z),abs(max_z)), 0, max(abs(min_z),abs(max_z))),
              cex=cex.legend,xpd=NA)
    }
}
