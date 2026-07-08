## =============================================================================
## evalAdmix visualization functions (base R graphics)
##
## Public functions:
##   plotCorRes      - single-panel correlation of residuals (rewritten, same API)
##   plotCorResMulti - multi-panel comparison of several correlation matrices
##                     with shared color scale and outer-margin labels
##   orderInds       - order individuals for plots (unchanged)
##   orderK          - order ancestral populations (unchanged)
##   plotAdmix       - barplot of admixture proportions (unchanged)
##
## Internal helpers (prefix "."):
##   .corres_prepare          - ordering + triangle replacement
##   .corres_resolve_z        - resolve NA min_z/max_z
##   .corres_palette           - build ramp palette + breaks
##   .corres_draw_grid         - draw one heatmap panel + diagonal mask
##   .corres_draw_labels       - draw pop/superpop labels and separator lines
##   .corres_draw_legend       - draw the color legend panel
## =============================================================================



## ---------------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------------

.corres_prepare <- function(cor_mat, pop, ord, summary_triangle="population_mean",
                             summary_group=NULL) {
  # Ordering + (optional) replacement of one triangle with population means.
  # summary_group: optional vector (same length as pop) specifying the grouping
  #   for computing triangle means. If NULL, uses pop itself.
  # Returns list(cor_mat=ordered matrix with diag untouched, pop=ordered pop).
  N <- nrow(cor_mat)

  if (is.null(ord) & !is.null(pop)) ord <- order(pop)
  if (is.null(ord) & is.null(pop)) ord <- 1:N

  if (is.null(pop)) {
    pop <- rep(" ", N)
  }

  pop <- pop[ord]
  cor_mat <- cor_mat[ord, ord]

  # Grouping for the triangle means: defaults to pop, but can be a
  # finer grouping (e.g. site within transect) passed via summary_group.
  grp <- if (is.null(summary_group)) pop else summary_group[ord]

  if (summary_triangle == "population_mean" && length(unique(grp)) > 0) {
    ups <- unique(grp)
    mean_cors <- matrix(NA, nrow=length(ups), ncol=length(ups))
    rownames(mean_cors) <- ups
    colnames(mean_cors) <- ups

    for (i1 in seq_along(ups)) {
      for (i2 in seq_along(ups)) {
        p1 <- ups[i1]; p2 <- ups[i2]
        sub <- cor_mat[which(grp==p1), which(grp==p2)]
        mean_cors[i1, i2] <- mean(sub[!is.na(sub)])
      }
    }

    for (i1 in 1:(N-1)) {
      for (i2 in (i1+1):N) {
        cor_mat[i1, i2] <- mean_cors[grp[i2], grp[i1]]
      }
    }
  }

  list(cor_mat=cor_mat, pop=pop)
}


.corres_resolve_z <- function(cor_mat, min_z, max_z, symmetric=FALSE) {
  # Resolve NA min_z/max_z preserving original evalAdmix semantics.
  # symmetric=TRUE: mirror whichever bound is given about 0.
  z_lims <- c(min_z, max_z)

  if (all(is.na(z_lims))) {
    m <- max(abs(cor_mat[!is.na(cor_mat)]))
    z_lims <- c(-m, m)
  }
  if (any(is.na(z_lims))) {
    known <- z_lims[!is.na(z_lims)]
    z_lims <- c(-known, known)
  }

  z_lims
}


.corres_palette <- function(min_z, max_z, color_palette, nHalf=10) {
  Min <- min_z; Max <- max_z; Thresh <- 0

  rc1 <- colorRampPalette(colors=color_palette[1:2], space="Lab")(nHalf)
  rc2 <- colorRampPalette(colors=color_palette[2:3], space="Lab")(nHalf)
  rampcols <- c(rc1, rc2)

  rampcols[c(nHalf, nHalf+1)] <-
    rgb(t(col2rgb(color_palette[2])), maxColorValue=256)

  rb1 <- seq(Min, Thresh, length.out=nHalf+1)
  rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks <- c(rb1, rb2)

  rlegend <- as.raster(matrix(rampcols, ncol=1)[length(rampcols):1,])

  list(cols=rampcols, breaks=rampbreaks, rlegend=rlegend)
}


.corres_draw_grid <- function(cm, pal, min_z, max_z, title, cex.main) {
  # Draws one heatmap panel using palette `pal` (list from .corres_palette).
  # cm should have diag set to 10 for the black diagonal mask.
  image(t(cm), col=pal$cols, breaks=pal$breaks,
        yaxt="n", xaxt="n", zlim=c(min_z, max_z), useRaster=TRUE,
        main=title, oldstyle=TRUE, cex.main=cex.main, xpd=NA)
  image(ifelse(t(cm > max_z), 1, NA), col="darkred", add=TRUE)
  if (min(cm, na.rm=TRUE) < min_z)
    image(ifelse(t(cm < min_z), 1, NA), col="darkslateblue", add=TRUE)
  image(ifelse(t(cm == 10), 1, NA), col="black", add=TRUE)
}


.corres_draw_labels <- function(pop, superpop, N,
                                label_y=TRUE, label_x=TRUE,
                                draw_lines=TRUE,
                                adjlab=0.1, adjlabsuperpop=0.16,
                                cex.lab=1.5, cex.lab.2=1.5,
                                rotatelabpop=0, rotatelabsuperpop=0,
                                lineswidth=1, lineswidthsuperpop=2,
                                summary_group=NULL, summarylineswidth=0.5) {
  # Draws pop labels on the left (y) and/or bottom (x) margins,
  # the corresponding separator lines, and optionally superpop labels.
  # If summary_group is provided, draws thin separator lines at its
  # boundaries (inside the plot area) before the pop/superpop lines.
  if (!is.null(summary_group) && draw_lines && length(unique(summary_group)) > 1) {
    abline(v=grconvertX(cumsum(sapply(unique(summary_group),
                                     function(x) sum(summary_group==x)))/N,
                        "npc", "user"),
           col="grey50", lwd=summarylineswidth, xpd=FALSE)
    abline(h=grconvertY(cumsum(sapply(unique(summary_group),
                                     function(x) sum(summary_group==x)))/N,
                        "npc", "user"),
           col="grey50", lwd=summarylineswidth, xpd=FALSE)
  }

  if (!is.null(pop)) {
    if (label_x)
      text(sort(tapply(1:length(pop), pop, mean)/length(pop)),
           -adjlab, unique(pop), xpd=NA, cex=cex.lab, srt=rotatelabpop)
    if (label_y)
      text(-adjlab, sort(tapply(1:length(pop), pop, mean)/length(pop)),
           unique(pop), xpd=NA, cex=cex.lab, srt=90-rotatelabpop)
    if (draw_lines) {
      abline(v=grconvertX(cumsum(sapply(unique(pop),
                                       function(x) sum(pop==x)))/N,
                          "npc", "user"),
             col=1, lwd=lineswidth, xpd=FALSE)
      abline(h=grconvertY(cumsum(sapply(unique(pop),
                                       function(x) sum(pop==x)))/N,
                          "npc", "user"),
             col=1, lwd=lineswidth, xpd=FALSE)
    }
  }

  if (!is.null(superpop)) {
    if (label_x)
      text(sort(tapply(1:length(superpop), superpop,
                      mean)/length(superpop)),
           -adjlabsuperpop, unique(superpop), xpd=NA,
           cex=cex.lab.2, srt=rotatelabsuperpop, font=2)
    if (label_y)
      text(-adjlabsuperpop, sort(tapply(1:length(superpop), superpop,
                                       mean)/length(superpop)),
           unique(superpop), xpd=NA, cex=cex.lab.2,
           srt=90-rotatelabsuperpop, font=2)
    if (draw_lines) {
      abline(v=grconvertX(cumsum(sapply(unique(superpop),
                                       function(x) sum(superpop==x)))/N,
                          "npc", "user"),
             col=1, lwd=lineswidthsuperpop, xpd=FALSE)
      abline(h=grconvertY(cumsum(sapply(unique(superpop),
                                       function(x) sum(superpop==x)))/N,
                          "npc", "user"),
             col=1, lwd=lineswidthsuperpop, xpd=FALSE)
    }
  }
}


.corres_draw_legend <- function(rlegend, min_z, max_z, cex.legend,
                                mar=c(5, 0.5, 4, 2)) {
  # Draws the legend panel. Caller controls par(mar).
  par(mar=mar)
  plot(c(0,1), c(0,1), type='n', axes=FALSE, xlab='', ylab='', main='')
  rasterImage(rlegend, 0, 0.25, 0.4, 0.75)
  text(x=0.8, y=c(0.25, 0.5, 0.75),
       labels=c(-max(abs(min_z), abs(max_z)), 0,
                max(abs(min_z), abs(max_z))),
       cex=cex.legend, xpd=NA)
}


## ---------------------------------------------------------------------------
## plotCorRes (single-panel, rewritten to use helpers; same public API)
## ---------------------------------------------------------------------------

plotCorRes <- function(cor_mat, pop=NULL, ord=NULL, superpop=NULL,
                       title="Correlation of residuals", min_z=NA, max_z=NA,
                       cex.main=1.5, cex.lab=1.5, cex.legend=1.5,
                       color_palette=c("#001260", "#EAEDE9", "#601200"),
                       pop_labels = c(T,T), plot_legend = T, adjlab = 0.1,
                       rotatelabpop=0, rotatelabsuperpop=0,
                       lineswidth=1, lineswidthsuperpop=2,
                       adjlabsuperpop=0.16, cex.lab.2 = 1.5, multipanel=FALSE,
                       summary_triangle = "population_mean",
                       summary_group = NULL) {

  N <- nrow(cor_mat)
  lineswidth_local <- lineswidth

  if (is.null(pop)) {
    pop <- rep(" ", N)
    lineswidth_local <- 0
  }

  prep <- .corres_prepare(cor_mat, pop, ord, summary_triangle, summary_group)
  cor_mat <- prep$cor_mat
  pop <- prep$pop
  superpop <- if (!is.null(superpop)) superpop[ord] else NULL

  z_lims <- .corres_resolve_z(cor_mat, min_z, max_z)
  min_z <- z_lims[1]; max_z <- z_lims[2]

  pal <- .corres_palette(min_z, max_z, color_palette)

  cm <- cor_mat
  diag(cm) <- 10

  if (!multipanel) {
    if (plot_legend) {
      layout(matrix(1:2, ncol=2), width=c(4,1), height=c(1,1))
      par(mar=c(5,4,4,0), oma=c(1,4.5,2,0))
    } else {
      par(mar=c(5,4,4,5), oma=c(1,4.5,2,0))
    }
  }

  .corres_draw_grid(cm, pal, min_z, max_z, title, cex.main)

  .corres_draw_labels(pop, superpop, N,
                      label_y=pop_labels[1], label_x=pop_labels[2],
                      draw_lines=TRUE,
                      adjlab=adjlab, adjlabsuperpop=adjlabsuperpop,
                      cex.lab=cex.lab, cex.lab.2=cex.lab.2,
                      rotatelabpop=rotatelabpop,
                      rotatelabsuperpop=rotatelabsuperpop,
                      lineswidth=lineswidth_local,
                      lineswidthsuperpop=lineswidthsuperpop)

  if (plot_legend && !multipanel) {
    .corres_draw_legend(pal$rlegend, min_z, max_z, cex.legend)
  }
}


## ---------------------------------------------------------------------------
## plotCorResMulti (multi-panel comparison with shared scale and outer labels)
## ---------------------------------------------------------------------------

plotCorResMulti <- function(
    cor_mats,
    titles = if (!is.null(names(cor_mats))) names(cor_mats)
             else paste("Panel", seq_along(cor_mats)),
    pop = NULL, ord = NULL, superpop = NULL,
    nrow = NULL, ncol = NULL,
    shared_scale = TRUE,
    label_outer = TRUE,
    summary_triangle = "population_mean",
    summary_group = NULL,
    min_z = NA, max_z = NA,
    cex.main = 1.5, cex.lab = 1.5, cex.legend = 1.5, cex.lab.2 = 1.5,
    color_palette = c("#001260", "#EAEDE9", "#601200"),
    pop_labels = c(TRUE, TRUE),
    plot_legend = TRUE,
    adjlab = 0.1, adjlabsuperpop = 0.16,
    rotatelabpop = 0, rotatelabsuperpop = 0,
    lineswidth = 1, lineswidthsuperpop = 2,
    summarylineswidth = 0.5,
    oma = c(4.5, 4.5, 2, 2),
    mar = c(1.5, 1.5, 2, 1.5)
) {

  if (!is.list(cor_mats))
    stop("cor_mats must be a list of matrices.")
  nPanels <- length(cor_mats)
  if (nPanels < 1)
    stop("cor_mats must contain at least one matrix.")
  if (is.null(titles))
    titles <- paste("Panel", seq_len(nPanels))
  if (length(titles) != nPanels)
    titles <- rep_len(titles, nPanels)

  N <- nrow(cor_mats[[1]])

  if (is.null(ord) & !is.null(pop)) ord <- order(pop)
  if (is.null(ord) & is.null(pop)) ord <- 1:N

  if (is.null(pop))
    pop <- rep(" ", N)
  pop_local <- pop[ord]
  superpop_local <- if (!is.null(superpop)) superpop[ord] else NULL
  summary_group_local <- if (!is.null(summary_group)) summary_group[ord] else NULL

  # Prepare every panel (ordering + triangle replacement) once.
  prep_list <- lapply(cor_mats, function(m) {
    p <- .corres_prepare(m, pop, ord, summary_triangle, summary_group)
    p$cor_mat
  })

  # Resolve shared scale.
  if (shared_scale) {
    if (all(is.na(c(min_z, max_z)))) {
      mm <- max(abs(unlist(lapply(prep_list,
                                  function(m) m[!is.na(m)]))))
      z_lims <- c(-mm, mm)
    } else {
      z_lims <- .corres_resolve_z(NULL, min_z, max_z)
    }
  } else {
    z_lims <- NULL  # resolved per-panel inside the loop below
  }
  if (!is.null(z_lims)) { min_z <- z_lims[1]; max_z <- z_lims[2] }
  shared_pal <- if (!is.null(z_lims))
    .corres_palette(min_z, max_z, color_palette) else NULL

  # Layout grid.
  if (is.null(nrow) && is.null(ncol)) {
    g <- grDevices::n2mfrow(nPanels)
    nr <- g[1]; nc <- g[2]
  } else if (is.null(nrow)) {
    nc <- ncol
    nr <- ceiling(nPanels / nc)
  } else if (is.null(ncol)) {
    nr <- nrow
    nc <- ceiling(nPanels / nr)
  } else {
    nr <- nrow; nc <- ncol
  }

  # Build layout matrix: panels 1..nPanels padded with 0; legend column at right.
  panel_mat <- matrix(seq_len(nr * nc), nrow=nr, ncol=nc, byrow=TRUE)
  panel_mat[panel_mat > nPanels] <- 0

  has_legend <- isTRUE(plot_legend)
  if (has_legend) {
    legend_id <- nPanels + 1
    layout_mat <- cbind(panel_mat, rep(legend_id, nr))
    widths  <- c(rep(4, nc), 1)
  } else {
    layout_mat <- panel_mat
    widths  <- rep(4, nc)
  }
  heights <- rep(1, nr)

  op <- par(no.readonly=TRUE)
  on.exit(par(op))

  # Outer margins host left (y) and bottom (x) labels and titles.
  par(oma=oma, mar=mar)
  layout(layout_mat, widths=widths, heights=heights)

  for (k in seq_len(nPanels)) {
    cm <- prep_list[[k]]
    if (is.null(shared_pal)) {
      cur_z <- .corres_resolve_z(cm, min_z, max_z)
      cur_pal <- .corres_palette(cur_z[1], cur_z[2], color_palette)
      cur_min <- cur_z[1]; cur_max <- cur_z[2]
    } else {
      cur_pal  <- shared_pal
      cur_min  <- min_z; cur_max <- max_z
    }

    diag(cm) <- 10

    # Panel position within the grid.
    r <- (k - 1) %/% nc + 1  # 1-based row index
    c <- (k - 1) %%  nc + 1  # 1-based column index

    # Determine which labels to draw on this panel.
    if (label_outer) {
      draw_y <- c == 1
      draw_x <- r == nr
    } else {
      draw_y <- TRUE
      draw_x <- TRUE
    }

    .corres_draw_grid(cm, cur_pal, cur_min, cur_max,
                      titles[k], cex.main)

    .corres_draw_labels(
      pop=pop_local, superpop=superpop_local, N=N,
      label_y=pop_labels[1] && draw_y,
      label_x=pop_labels[2] && draw_x,
      draw_lines=TRUE,
      adjlab=adjlab, adjlabsuperpop=adjlabsuperpop,
      cex.lab=cex.lab, cex.lab.2=cex.lab.2,
      rotatelabpop=rotatelabpop,
      rotatelabsuperpop=rotatelabsuperpop,
      lineswidth=lineswidth,
      lineswidthsuperpop=lineswidthsuperpop,
      summary_group=summary_group_local,
      summarylineswidth=summarylineswidth
    )
  }

  if (has_legend) {
    # First draw using the shared scale; legend panel.
    if (!is.null(shared_pal)) {
      .corres_draw_legend(shared_pal$rlegend, min_z, max_z, cex.legend,
                          mar=c(0.2, 0.2, 1.5, 2))
    } else {
      # No shared scale: take the last panel's scale as a reasonable fallback.
      last <- prep_list[[nPanels]]
      last_z <- .corres_resolve_z(last, NA, NA)
      lp <- .corres_palette(last_z[1], last_z[2], color_palette)
      .corres_draw_legend(lp$rlegend, last_z[1], last_z[2], cex.legend,
                          mar=c(0.2, 0.2, 1.5, 2))
    }
  }

  invisible(NULL)
}


## ---------------------------------------------------------------------------
## orderInds (unchanged from original)
## ---------------------------------------------------------------------------

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


## ---------------------------------------------------------------------------
## orderK (unchanged from original)
## ---------------------------------------------------------------------------

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


## ---------------------------------------------------------------------------
## plotAdmix (unchanged from original)
## ---------------------------------------------------------------------------

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