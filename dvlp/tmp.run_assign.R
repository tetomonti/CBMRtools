assign_heatmap<-function(x,
                         gene_scores,
                         sample_scores,
                         sample_probs,
                         posteriors,
                         inclusion_cutoff,
                         gene_cols=rep('green',nrow(x)),
                         drop_outliers=TRUE,
                         fancy_order=TRUE,
                         colSideCols=NULL,
                         main='',
                         xlab="Samples",
                         ylab="Signature")
{
    ## support functions
    ##
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    ## ordering the heatmap
    gene_order<-order(gene_scores,decreasing=T)
    sample_order<-order(sample_scores,decreasing=T)   
    x<-x[gene_order,sample_order]
    
    if (fancy_order){
        fancy_score<-colSums(sweep(x,1,gene_scores,'*'))
        fancy_order<-order(fancy_score,decreasing=T)
        x<-x[,fancy_order]
        sample_scores<-sample_scores[fancy_order]
        sample_probs<-sample_probs[fancy_order]   
    }
    rm <- rowMeans(x, na.rm = T)
    x <- sweep(x, 1, rm)
    sx <- apply(x, 1, sd, na.rm = T)
    x <- sweep(x, 1, sx, "/")
    
    if (drop_outliers){
        x[x>4]   <- 4
        x[x<(-4)]<- -4
    }
    
    ## defining the canvas
    lhei <- c(1.5, 4)
    lwid <- c(1.5, 4)
    lmat <- rbind(4:3, 2:1)
    
    if (!is.null(colSideCols)) {
        lmat <- rbind(lmat[1,] + 1, c(0, 1), lmat[2, ] +  1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!is.null(colSideCols)) {
        par(mar = c(0.1, 1.1, 0.5, 2))
        image(cbind(1:ncol(x)), col = colSideCols[sample_order], axes = FALSE)
    }
    
    ## heatmap
    col<-brewer.pal(11,"RdBu")[11:1]
    breaks<-seq(min(x, na.rm = T), max(x, na.rm = T), length = length(col)+1)
    par(mar=c(2,1.1,0.5,2))
    image(t(x[nrow(x):1,]), 
          axes = FALSE, 
          xlab = '', 
          ylab = ylab, 
          col = col,
          breaks=breaks)
    mtext(xlab,1,cex=1.5,line=0.5)
    mtext(ylab,4,cex=1.5,line=0.5)
    
    gene_cols[posteriors<inclusion_cutoff | gene_scores <0]<-'grey'   
    bar_ord<-order(gene_scores)
    par(mar=c(0.9,4,0,0))
    barplot(-1*abs(gene_scores[bar_ord]),
            horiz=T,
            border=NA,
            axes=F,
            col=gene_cols[bar_ord],
            space=0)
    mtext('Gene scores',2,line=-1.5)
  
    ## top histogram
    par(mar=c(0,0,7,1.2))
    barplot(sort(sample_probs,decreasing=T),
            border=NA,
            col='grey',
            axes=F,
            main='',
            space=0)
    barplot(sort(sample_scores,decreasing=T),
            border=NA,
            col='purple',
            space=0,
          axes=F,
          add=T)
    axis(2,at=seq(0,1,0.2),tick=T,line=-0.8)
    axis(1,labels=F)
    legend('topright',
           legend=c('Score', 'Activation Prob.'),
           col=c('purple','grey'),
           pch=15,
           cex=1,
           pt.cex=1)
    
    title(main, cex.main = 2.5)
  
    ## key
    par(mar = c(5, 2, 2, 3), cex = 0.75)
    tmpbreaks <- breaks
    min.raw <- min(x, na.rm = TRUE)
    max.raw <- max(x, na.rm = TRUE)
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    mtext(side = 1, "Row Z-Score", line = 2)
    
    h <- hist(x, plot = FALSE, breaks = breaks)
    hx <- scale01(breaks, min.raw, max.raw)
    hy <- c(h$counts, h$counts[length(h$counts)])
    lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
          col = "cyan")
    axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
    title("Color Key\nand Histogram")
    par(cex = 0.5)
}
