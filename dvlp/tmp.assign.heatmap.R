assign_heatmap <- function
(
 expression_matrix,
 kappa_pos,
 S_pos,
 delta_pos,
 gamma_pos=NULL,
 posterior_cutoff=0.95,
 xlab='Samples',
 ylab='Signature',
 assign_score_lab='Assign Score',
 S_pos_lab='Gene scores',
 gene_cols=rep('green',nrow(expression_matrix)),
 drop_outliers=T,
 colSideCols=NULL
 )
{
    ## reduce expression dataset to the gene set space
    x<-expression_matrix
    x<-x[rownames(x)%in%names(S_pos),]
    
    ## match sample and gene names
    x<-x[match(names(S_pos),rownames(x)),]
    x<-x[,match(names(kappa_pos),colnames(x))]
    
    ## order the matrix according to the scores
    gene_order<-order(S_pos,decreasing=T)
    sample_order<-order(kappa_pos,decreasing=T)
    
    x<-x[gene_order,sample_order]
    S_pos<-S_pos[gene_order]
    delta_pos<-delta_pos[gene_order]
    kappa_pos<-kappa_pos[sample_order]
    
    if(!is.null(gamma_pos)){
        gamma_pos<-gamma_pos[sample_order]
    }
    ## row scaling of the expression data
    rm <- rowMeans(x, na.rm = T)
    x <- sweep(x, 1, rm)
    sx <- apply(x, 1, sd, na.rm = T)
    x <- sweep(x, 1, sx, "/")
    
    ## removing the outliers by reducing the row scaled log2 data to a range of +/-4
    if (drop_outliers){
        x[x>4]   <- 4
        x[x<(-4)]<- -4
    }
    ## heatmap
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
    col<-c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
    breaks<-seq(min(x, na.rm = T), max(x, na.rm = T), length = length(col)+1)
    orig.mar <- par(mar=c(2,1.1,0.5,2))
    image(t(x[nrow(x):1,]),
          axes = FALSE,
          xlab = '',
          ylab = "",
          col = col,
          breaks=breaks)
    mtext(xlab,1,cex=1,line=0.5)
    mtext(ylab,4,cex=1,line=0.5)
    
    ## bar plot on the left
    gene_cols[delta_pos<posterior_cutoff]<-'grey'
    bar_ord<-order(S_pos)
    par(mar=c(0.9,7,0,0))
    barplot(-1*abs(S_pos[bar_ord]),
            axisnames=F,
            horiz=T,
            border=NA,
            axes=F,
            col=gene_cols[bar_ord],
            space=0)
    mtext('Gene scores',2,line=-1.5)
    
    ## bar plot on the top
    
    par(mar=c(0,0,7,1.2))
    names(kappa_pos)<-''
    add=F
    if(!is.null(gamma_pos)){
        names(gamma_pos)<-''
        barplot(gamma_pos,
                border=NA,
                col='grey',
                axes=F,
                main='',
                space=0)
        add=T
    }
    barplot(kappa_pos,
            border=NA,
            col='purple',
            space=0,
            axes=F,
            add=add)
    
    axis(2,at=seq(0,1,0.2),tick=T)
    axis(1,labels=F)
    if(!is.null(gamma_pos)){
        legend('topright',
               legend=c('kappa_pos', 'gamma_pos'),
               col=c('purple','grey'),
               pch=15,
               cex=1,
               pt.cex=1)
    }
    ## key
    par(mar = c(5, 2, 2, 3), cex = 0.75)
    tmpbreaks <- breaks
    min.raw <- min(x, na.rm = TRUE)
    max.raw <- max(x, na.rm = TRUE)
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- (as.numeric(lv) - min.raw)/(max.raw - min.raw)
    axis(1, at = xv, labels = lv)
    mtext(side = 1, "Row Z-Score", line = 2)
    
    
    h <- hist(x, plot = FALSE, breaks = breaks)
    hx <- (breaks - min.raw)/(max.raw - min.raw)
    
    hy <- c(h$counts, h$counts[length(h$counts)])
    lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
          col = "cyan")
    axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
    title("Color Key\nand Histogram")
    par(cex = 0.5)
}
