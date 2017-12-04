#require(ASSIGN)
#require(Biobase)
#require(RColorBrewer)

## function: RUN ASSIGN
##
## This function runs ASSIGN given an expression set object
     
run_assign<-function(ES,         # ExpressionSet object 
                     eSig,       # Expression signature used by ASSIGN (vector of valid gene symbols)
                     oDir,       # Output directory (folder is created automatically)
                     iter=3000,  # Number of MCMC iterations
                     beta=TRUE,  # Whether or not to use mixture beta when running assign 
                     do.log=TRUE,# log2-transform the data
                     geneID="gene_symbol"
)
{
  cat("Number of genes in signature: ",length(eSig),"\n")
  cat("Number of genes from signature that overlap with features in ESet: ",sum(fData(ES)$gene_symbol %in% eSig),"\n")
  cat("Running ASSIGN with beta mixture modelling set to ",beta," ..\n")
  cat("")
  cat("Saving output to directory: ",eval(oDir),"\n")

  ## checks
  if ( !(geneID %in% colnames(fData(ES))) ) stop( "geneID not found: ",geneID)
  
  ## Create output directory
  dir.create(oDir,showWarnings=F)
  
  mat<-exprs(ES)[fData(ES)[,geneID] %in% eSig,]
  rownames(mat)<-fData(ES)[,geneID][fData(ES)[,geneID] %in% eSig]
  gene_list<-rownames(mat)
  if(length(gene_list) > length(unique(gene_list)))
    stop(paste("Duplicate entries present in matrix when filtering in for gene list..\n Please either remove these genes or consolidate before projecting onto dataset:\n",
               names(which(table(gene_list)>1)),sep="\n"))
  
  ## Take log transform (if only normalized RNASeq data), adding pseudocount to avoid NA's
  if (do.log) {
    if ( max(mat)<100 ) warning("data migth be log-transformed twice")
    cat("Taking log2-transorm of expression values ..\n\n")
    mat<-log2(mat+1)
  }
  ## Start assign
  processed.data <- assign.preprocess(trainingData=NULL,
                                      testData=mat,
                                      trainingLabel=NULL,
                                      geneList=list(gene_list),
                                      n_sigGene=NA,
                                      theta0=0.05,
                                      theta1=0.9)
  
  #Run mcmc
  mcmc.chain <- assign.mcmc(Y=processed.data$testData_sub, 
                            Bg = processed.data$B_vector, 
                            X=processed.data$S_matrix, 
                            Delta_prior_p = processed.data$Pi_matrix, 
                            iter = iter, 
                            adaptive_B=TRUE, 
                            adaptive_S=TRUE, 
                            mixture_beta=beta)
  
  mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=1000, 
                                  iter=iter, adaptive_B=TRUE, 
                                  adaptive_S=TRUE,mixture_beta=beta)
  
  assign.output(processed.data=processed.data, 
                mcmc.pos.mean.testData=mcmc.pos.mean, 
                trainingData=NULL, 
                testData=mat, 
                trainingLabel=NULL, 
                testLabel=NULL, 
                geneList=list(gene_list), 
                adaptive_B=TRUE, 
                adaptive_S=TRUE, 
                mixture_beta=beta, 
                outputDir=oDir)
  
  output.data <- list( processed.data=processed.data, mcmc.pos.mean.testData=mcmc.pos.mean )
  output.data <- annotateAssignOutput( output.data )
  
  save( output.data, file=paste(oDir,"output.rda",sep='/') )
  cat("Done!\n")
  return( output.data)
}
## function: ANNOTATE ASSIGN OUTPUT
##
annotateAssignOutput <- function
(
  assignOutput,
  pathways=NULL
)
{
  ## Gxk matrix (or vector) of gene weights
  ##
  Pi_matrix <- assignOutput$processed.data$Pi_matrix

  if ( is.null(dim(Pi_matrix)) ) {
    names(assignOutput$processed.data$B_vector) <-
      names(assignOutput$processed.data$S_matrix) <- 
        names(Pi_matrix)
  }
  ## gene annotation
  rownames(assignOutput$mcmc.pos.mean.testData$S_pos) <-
    rownames(assignOutput$mcmc.pos.mean.testData$Delta_pos) <-
      if ( is.null(dim(Pi_matrix)) ) names(Pi_matrix) else rownames(Pi_matrix)

  ## pathway annotation
  if ( !is.null(pathways) && !is.null(names(pathways)) )
    colnames(assignOutput$mcmc.pos.mean.testData$S_pos) <-
        colnames(assignOutput$mcmc.pos.mean.testData$Delta_pos) <- names(pathways)

  ## sample annotation
  rownames(assignOutput$mcmc.pos.mean.testData$kappa_pos) <-
    rownames(assignOutput$mcmc.pos.mean.testData$gamma_pos) <- 
        colnames(assignOutput$processed.data$testData_sub)

  assignOutput
}
## function: ASSIGN HEATMAP
##
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
## function: PLOT ALL
##
plotAll<-function
(
    eSet,
    output.data,
    title,
    gene_list,
    inclusion_cutoff=0.75,
    colSideCols=NULL,
    do.log=TRUE,
    xlab="Samples",
    ylab="Signature",
    geneID="gene_symbol"
)
{
  ## checks
  if ( !(geneID %in% colnames(fData(eSet))) ) stop( "geneID not found: ",geneID)

  ## Take log transform, adding pseudocount to avoid NA's
  if (do.log) {
    if ( max(exprs(eSet))<100 ) warning("data migth be log-transformed twice")
    cat("Taking log2-transorm of expression values ..\n\n")
    exprs(eSet) <- log2( exprs(eSet)+1 )
  }
  assay.data<-log2(exprs(eSet)+1)
  assay.data<-assay.data[fData(eSet)[,geneID] %in% gene_list,]
  rownames(assay.data)<-fData(eSet)[,geneID][fData(eSet)[,geneID] %in% gene_list]
  gene_list<-rownames(assay.data)
  
  ## all
  hc01.col <- hclust(dist(t(assay.data)),method="ward.D")
  hc01.row <- hclust(as.dist(1-cor(t(assay.data))),method="ward.D")
  
  if (!is.null(colSideCols)){
    heatmap.2(assay.data,
              scale='row',
              trace='none',
              main=title,
              Colv=as.dendrogram(hc01.col),
              Rowv=as.dendrogram(hc01.row),
              col=brewer.pal(11,"RdBu")[11:1],
              ColSideColors=colSideCols)
  }else{
    heatmap.2(assay.data,
              scale='row',
              trace='none',
              main=title,
              Colv=as.dendrogram(hc01.col),
              Rowv=as.dendrogram(hc01.row),
              col=brewer.pal(11,"RdBu")[11:1])   
  }
  
  gene_scores <- output.data$mcmc.pos.mean.testData$S_pos[match.nona(rownames(assay.data),names(output.data$processed.data$Pi_matrix))]
  posteriors <- output.data$mcmc.pos.mean.testData$Delta_pos[match.nona(rownames(assay.data),names(output.data$processed.data$Pi_matrix))]
  sample_scores<-output.data$mcmc.pos.mean.testData$kappa_pos
  sample_probs<-output.data$mcmc.pos.mean.testData$gamma_pos
  
  assign_heatmap(assay.data,
                 gene_scores,
                 sample_scores,
                 sample_probs,
                 posteriors,
                 inclusion_cutoff,
                 colSideCols=colSideCols,
                 main=title,xlab=xlab,ylab=ylab)
}

#Read arguments from command line
#args<-commandArgs(TRUE)
#ES<-args[1]
#sig<-args[2]
#out_dir<-args[3]

#Run ASSIGN
#run_assign(ES=ES,eSig=sig,oDir=out_dir,iter=3000)
