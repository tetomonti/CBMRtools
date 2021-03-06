#' geneSetProjection
#' 
#' \code{geneSetProjection} maps a gene-level dataset to a
#' geneset-level dataset based on control-treatment pairing (i.e., the
#' ks score is done with respect to the 'control vs. treatment'
#' phenotype).  It can be used to carry out 'single sample projection'
#' (the average of all samples is used as controls), as well as with
#' 'time-series' data of response to treatment, where for each time
#' point there are several replicates for both the control and the
#' treated.
#'
#' @param dat ExpressionSet data object ( require(biobase) )
#' @param pairing list with pairing between treatment and control (see format below)
#' @param gsets GeneSet data object (see GeneSet.R)
#' @param method one of {"ks","logistic","median","mean"} (only ks implemented so far)
#' @param collapse collapse multiple replicates into a single output vector (TRUE)
#' @param weighted gsea-like weighting of KS score (TRUE)
#' @param absolute use absolute values when calculating the KS score ignoring up and down (FALSE)
#' @param verbose verbosity on/off (TRUE)
#'
#' @return an expression set 
#'
#' @details
#'
#' The format of 'pairing' is a list of lists (of sample names), of the form:
#' 
#'   list(group1=list(control=c(s_111,..,s_11n),treatment=c(s_121,..,s_12n)),
#'        group2=list(control=c(s_211,..,s_21n),treatment=c(s_221,..,s_22n)),
#'        ...)
#' 
#' And the output ExpressionSet will have as many columns as either
#' the number of groups (if collapse=TRUE), or the number of samples
#' across all treatment groups (if collapse=FALSE)
#' 
#' @examples
#' # the very basic steps are:
#' 
#' # 1) load the data (expression and genesets)
#' # 2) rename the dataset by replacing probesetIDs by gene symbols ==> DB
#' # 3) create the list of lists pairing ==> PAIRS
#' # 4) run geneSetProjection(dat=DB,pairing=PAIRS,gset.db=GSET, ...)
#'
#' ## 1) LOAD the data
#' ##
#' data(gspData)
#' if (is.null(gsp.eSet)) stop("is.null(gsp.eSet)")
#' if (is.null(gsp.GeneSet)) stop("is.null(gsp.GeneSet)")
#' 
#' ## 2) RENAME the dataset rows (w/ gene symbols)
#'
#' DAT1 <- gsp.eSet[fData(gsp.eSet)[,"symbol"]!="",]
#' exprs(DAT1) <- log2(exprs(DAT1)+1)
#' featureNames(DAT1) <- toupper(fData(DAT1)[,"symbol"])
#' 
#' ## 3) CREATE the list of lists pairing (one group only)
#' 
#' ANidx <- pData(DAT1)[,"tissue_type"]=="AN"
#' PAIRS <- list(OSCC=list(control=sampleNames(DAT1)[ANidx],
#'                         treatment=sampleNames(DAT1)[!ANidx]))
#' 
#' ## 4) RUN geneSetProjection ..
#' ##
#' 
#' GSPdir <- geneSetProjection(dat=DAT1,
#'                             pairing=PAIRS,
#'                             GS=gsp.GeneSet,
#'                             collapse=FALSE, # single-sample projection
#'                             weighted=FALSE, # standard KS score, no weighting
#'                             absolute=FALSE, # keep sign of enrichment scores
#'                             min.gset=5,
#'                             verbose=TRUE)
#' 
#' ## shortening pathway names for better visualization
#' featureNames(GSPdir) <- c("WNT_BCATENIN","BCATENIN.CSNK1E","BCATENIN.CTBP1","TGF_BETA")
#'
#' gradeID <- 'my_grade'
#' stageID <- 'my_stage'
#' p2 <- heatmap.ggplot2(eSet=GSPdir,col.clust=TRUE,row.clust=TRUE,
#'                       col.lab=c(gradeID,stageID),row.lab="",
#'                       heatmap.y.text=TRUE, heatmap.x.text=FALSE,
#'                       heatmap.colorlegend.name="RNASeq_expression",
#'                       title.text="TCGA BRCA log2 gene set projection",
#'                       col.legend.name=c(gradeID,stageID), row.legend.name="", 
#'                       row.scaling="none",z.norm=FALSE, 
#'                       cuttree.col=0, cuttree.row=0,
#'                       verbose=FALSE, show=TRUE)
#' grid.arrange(p2)
#'
#' @export

## END documentation support
#####################################################################################
## GENESET PROJECTION
##
## function to map a gene-level dataset to a geneset-level dataset
## based on control-treatment pairing (i.e., the ks score is done
## with respect to the 'control vs. treatment' phenotype).
## This function is ideal for 'time-series' data of response to
## treatment, where for each time point there are several replicates
## for both the control and the treated.
## 
geneSetProjection <- function
(
 dat,           # expression data
 pairing,       # list with pairing between treatment and control (see format below)
 GS,            # gene set object
 method=c("ks","logistic","median","mean"), # only ks implemented at the moment
 collapse=TRUE, # collapse multiple replicates into a single output vector
 weighted=TRUE, # gsea-like weighting of KS score
 absolute=FALSE,# use absolute values when calculating the KS score ignoring up and down
 min.gset=5,    # min genes in a geneset allowed
 verbose=TRUE
 )
{
  ## FORMAT: pairing is a list of lists (of sample names), of the form:
  ##   list(group1=list(control=c(s_111,..,s_11n),treatment=c(s_121,..,s_12n)),
  ##        group2=list(control=c(s_211,..,s_21n),treatment=c(s_221,..,s_22n)),
  ##        ...)
  
  ## checks on inputs
  ##
  method <- match.arg(method)
  if ( class(dat)!='ExpressionSet' )
    stop( 'dat expected to be an ExpressionSet object: ', class(dat) )
  if ( class(GS)!='GeneSet' )
    stop( 'GS expected to be a GeneSet object: ', class(GS) )
  if ( class(pairing)!='list' )
    stop( 'pairing expected to be a list: ', class(pairing) )      
  if ( method!="ks" )
    stop( "method not implemented yet: ", method )
  if ( any(unlist(sapply(pairing,function(z) is.na(match(unlist(z),sampleNames(dat)))))) )
    stop( "mismatch btw pairing and data" )
  if ( any(sapply(pairing,function(z) any(names(z)!=c("control","treatment")))) )
    stop( "pairing groups must be named 'control' and 'treatment'" )

  ## drop the genes that are in the gene sets and not in the dataset
  setGeneSet(GS) <- lapply(getGeneSet(GS),intersect,toupper(featureNames(dat)))
  if ( length(GS)<1 ) stop( "no genesets overlap w/ data's genes" )

  VERBOSE(verbose,"Projecting",length(pairing),"groups ..\n")
  PRJ <- NULL
  for ( i in 1:length(pairing) )
  {
    VERBOSE(verbose,"  > group",names(pairing)[i],".. ")
    grpI <- pairing[[i]]
    ctlI <- match(grpI$control,sampleNames(dat))
    trtI <- match(grpI$treatment,sampleNames(dat))
    VERBOSE(verbose,"[",length(ctlI),"/",length(trtI),"]")
    if ( length(ctlI)<2 ) {
      VERBOSE(verbose,"insufficient control samples: ",length(ctlI),", skipping.\n",sep="")
      next
    }
    if ( length(trtI)<1 ) {
      VERBOSE(verbose,"insufficient treatment samples: ",length(trtI),", skipping\n",sep="")
      next 
    }
    datI <- dat[,c(ctlI,trtI)]
    clsI <- rep(0:1,times=c(length(ctlI),length(trtI))); levels(clsI) <- paste(c("ctl","trt"),0:1,sep=".")
    control <- !collapse && length(intersect(grpI$control,grpI$treatment))==0
    prjI <- ksProjection(exprs(datI),cls=clsI,gsets=getGeneSet(GS),collapse=collapse,weighted=weighted,absolute=absolute,control=control)
    PRJ <- cbind(PRJ,prjI)

    VERBOSE(verbose," done.\n")
  }
  outdat <- NULL
  TMP <- pData(dat)[as.vector(unlist(sapply(pairing,function(z) z$treatment))),,drop=FALSE]
  
  if ( collapse ) {
    ## *** WARNING: THIS HAS NOT BEEN ADEQUATELY TESTED ***
    outdat <- ExpressionSet(assayData=PRJ,
                            phenoData=AnnotatedDataFrame(pData(dat)[as.vector(unlist(sapply(pairing,function(z) z$treatment[1]))),]))
    #outdat <- dat[,as.vector(unlist(sapply(pairing,function(z) z$treatment[1])))]
    sampleNames(outdat) <- names(pairing)
  }
  else {
    outdat <- ExpressionSet(assayData=PRJ,
                            phenoData=AnnotatedDataFrame(pData(dat)[match(colnames(PRJ),sampleNames(dat)),]))
  }
  outdat
}
## function: KS PROJECTION
##
## project one or more samples from gene space onto geneset ks-score space
## the function can either:
## - turn many-cases vs. many-controls into a single ks-score vector (collapse=T, cls!=NULL)
## - turn many-cases vs. many-controls into many ks-score vectors (collapse=FALSE, cls!=NULL)
## - turn many-cases into many ks-score vectors (cls==NULL)
##
ksProjection <- function
(
    dat,              # [genes x samples] data matrix
    cls=NULL,         # class template (expected coding: 0=control; 1=treatment)
    gsets,            # list of (named) genesets
    collapse=FALSE,   # collapse multiple replicates
    control=!collapse,# include normals (only when collapse==FALSE)
    weighted=TRUE,
    do.check=TRUE,
    min.gset=10,
    verbose=TRUE,
    absolute=FALSE
)
{
  ## CHECKS ON INPUTS
  ##
  if ( !is.matrix(dat) ) {
    stop( 'dat must be a matrix: ', class(dat) )
  }
  if ( collapse && control )
     warning( 'inclusion of controls ignored when collapse==TRUE' )
  if ( !is.null(cls) && length(cls)!=ncol(dat) ) {
    stop( "length(cls)!=ncol(dat)" )
  }
  if ( !is.null(cls) && length(unique(cls))!=2 ) {
    stop( "cls must be binary" )
  }
  if ( do.check && length(cmn <- intersect(unique(unlist(gsets)),toupper(rownames(dat))))<min.gset ) {
    stop( "genes in common are less then ", min.gset, ": ",length(cmn) )
  }
  ## restrict to the geneset genes represented in the data (and to genesets w/ at least min.gset genes)
  ##
  if ( do.check ) {
      VERBOSE( verbose, "\n\tChecking genesets .." )
      gsets <- lapply(lapply(gsets,unique),intersect,cmn)
      gsets <- gsets[nz <- sapply(gsets,length)>=min.gset]
      if ( length(gsets)<1 ) stop( "no genesets w/ min # of genes required" )
      VERBOSE( verbose && sum(!nz)>0, "\t", sum(!nz), "genesets removed because of too few genes.\n")
      VERBOSE( verbose, "\tDone checking.\n")
  }
  ## map from gene names to gene positions within dataset
  ## (for max efficiency, this could ideally be brought out of the function)
  ##
  gidx <- lapply( gsets, function(Z) sort(match.nona(Z,toupper(rownames(dat)))) )
  
  ksProject <- function( dat, cls, var.equal=FALSE,absolute=FALSE,verbose=TRUE)
  {
    ## WARNING: gidx is a 'global' variable (not passed as argument)
    ##
    RNK <- rank(SCR <- tscore.4ks(dat,cls,var.equal=var.equal),ties.method="first")
    KS <- sapply(gidx, function(z)
        ksGenescore(n.x=nrow(dat),y=RNK[z],do.pval=FALSE,weight=if(weighted) sort(SCR),absolute=absolute))
    VERBOSE(verbose,".")
    return(KS)
  }
  ## return single vector from replicates
  ##
  if ( collapse ) {
    KS <- ksProject(dat,cls,var.equal=FALSE,absolute=absolute)
    return(KS)
  }
  ## return as many vectors as 'replicates'
  ##
  else if ( !is.null(cls) ) {
    ctl.idx <- which(cls==0)
    trt.idx <- which(cls==1)
    if (control) {
        KSctl <- sapply(ctl.idx,function(z) ksProject(dat=dat[,c(ctl.idx,z)],cls=c(cls[ctl.idx],1),
                                                      var.equal=TRUE,absolute=absolute))
        colnames(KSctl) <- colnames(dat)[ctl.idx]
    }
    KStrt <- sapply(trt.idx,function(z) ksProject(dat=dat[,c(ctl.idx,z)],cls=cls[c(ctl.idx,z)],
                                                  var.equal=TRUE,absolute=absolute))
    colnames(KStrt) <- colnames(dat)[trt.idx]
    return ( if (control) cbind(KSctl,KStrt) else KStrt )
  }
  ## no control group provided
  ##
  else {
    stop( "not implemented yet" )
  }
}
## function: TSCORE 4 KS
##
## fast computation of t-score based on matrix multiplications
## it also takes care of the case when one of the two groups is a singleton
## return value: positive => up in cls==0; negative => up in cls==1
##
tscore.4ks <- function( x, cls, var.equal=FALSE, robust=FALSE )
{
  lev <- sort(unique(cls))
  if (length(lev)!=2) stop( "cls must be binary" )
  CLS <- cbind( as.numeric(cls==lev[1]), as.numeric(cls==lev[2]) )
  nc <- apply(CLS,2,sum)

  ## both groups have multiple observations
  ##
  if ( all(nc>1) ) # use mean and stdev
  {
    s <- s2 <- NULL
    if ( !robust || all(nc<3) )
    {
      s  <- x %*% CLS
      s2 <- x^2 %*% CLS
      
      s2[,1] <- (s2[,1] - (s[,1]^2)/nc[1]) / (nc[1]-1)
      s2[,2] <- (s2[,2] - (s[,2]^2)/nc[2]) / (nc[2]-1)
      s[,1] <- s[,1]/nc[1]
      s[,2] <- s[,2]/nc[2]
    }     
    else           # use median and mad instead
    { 
      s <- sapply(lev,function(z) apply(x[,cls==z],1,median))
      s2 <- sapply(lev,function(z) apply(x[,cls==z],1,mad))^2
    }
    stderr <- if (var.equal)
      sqrt( (((nc[1]-1)*s2[,1] + (nc[2]-1)*s2[,2])/(nc[1]+nc[2]-2)) * (1/nc[1]+1/nc[2]) )
    else
      sqrt( s2[,1]/nc[1] + s2[,2]/nc[2] )
    
    (s[,1]-s[,2]) / stderr
  }
  ## one of the two groups has a single observation
  ##
  else if ( any(nc>1) )
  {
    ## determine the group w/ multiple observations
    iN <- which(nc>1)
    ncN <- nc[iN]
    clsN <- CLS[,iN,drop=FALSE]
    sgn <- if (iN==1) 1 else -1 # needed to keep consistent sign

    ## compute mean and standard deviation for that group
    s <- s2 <- NULL
    if ( !robust || ncN<3 )
    {
      s  <- drop(x %*% clsN)
      s2 <- drop(x^2 %*% clsN)
      s2 <- (s2 - (s^2)/ncN) / (ncN-1)
      s <- s/ncN
    }
    else {
      s <- apply(x[,cls==lev[iN]],1,median)
      s2 <- apply(x[,cls==lev[iN]],1,mad)^2     
    }
    ## I found a case where s2 was negative for a single value it was
    ## a very small value e-14 so I assume it was a rounding error,
    ## but sqrt of a negative number doesn't work so this is a rather
    ## crude way to fix this:
    s2<-abs(s2)
    
    ## compute score
    sgn * ( s - x[,cls!=lev[iN]] ) / sqrt(s2)
  }
  ## both groups have a single observation (just compute the difference)
  ##
  else {
    x[,cls==lev[1]] - x[,cls==lev[2]]
  }
}
## CONTROL STANDARDIZE
##
control.based.standardize <- function
(
 dat,            # ExpressionSet object
 pairing,        # list with pairing between treatment and control (see format below)
 robust=FALSE,   # use median/mad instead of mean/sd
 collapse=FALSE, # not implemented yet (but would collapse multiple treatment replicates into one)
 min.sd=.01,     # can be a scalar or a vector of length nrow(dat)
 verbose=TRUE
 )
{
  ## FORMAT: pairing is a list of lists, of the form:
  ##   list(group1=list(control=c(s_111,..,s_11n),treatment=c(s_121,..,s_12n)),
  ##        group2=list(control=c(s_211,..,s_21n),treatment=c(s_221,..,s_22n)),
  ##        ...)

  ## CHECKS on inputs
  ##
  if ( length(min.sd)==1 )
    min.sd <- rep(min.sd,nrow(dat))
  if ( length(min.sd)!=nrow(dat) )
    stop( "length(min.sd)!=nrow(dat)" )
  if ( collapse )
    stop( "collapse not implemented yet, ignored" )
  
  VERBOSE(verbose,"Control-standardizing",length(pairing),"groups ..\n")
  NRM <- NULL
  nrm.names <- NULL
  for ( i in 1:length(pairing) )
  {
    VERBOSE(verbose,"  > group",names(pairing)[i],".. ")
    grpI <- pairing[[i]]
    if ( any(is.na(ctlI <- match(grpI$control,sampleNames(dat)))) ) stop( "control samples missing from dat" )
    if ( any(is.na(trtI <- match(grpI$treatment,sampleNames(dat)))) ) stop( "treatment samples missing from dat" )
    VERBOSE(verbose,"[",length(trtI),"/",length(ctlI),"]")
    if ( length(ctlI)<2 ) {
      VERBOSE(verbose,"insufficient control samples: ",length(ctlI),", skipping.\n",sep="")
      next
    }
    if ( length(trtI)<1 ) {
      VERBOSE(verbose,"insufficient treatment samples: ",length(trtI),", skipping\n",sep="")
      next
    }
    mnI <- if (robust) apply(exprs(dat)[,ctlI],1,median) else drop(fast.mean(exprs(dat)[,ctlI]))
    sdI <- if (robust) apply(exprs(dat)[,ctlI],1,mad) else drop(fast.sd(exprs(dat)[,ctlI]))
    VERBOSE(verbose, " (thresholding ",round(100*sum(sdI<min.sd)/length(sdI),2),"% SDs)",sep="")
    sdI[sdI<min.sd] <- min.sd[sdI<min.sd]

    nrmI <- (dat@signal[,trtI] - mnI) / sdI
    NRM <- cbind(NRM,nrmI)
    nrm.names <- c(nrm.names,grpI$treatment)

    VERBOSE(verbose," done.\n")
  }
  if ( ncol(NRM)!=length(nrm.names) ) stop( "ncol(NRM)!=length(nrm.names):",ncol(NRM),length(nrm.names) )
  colnames(NRM) <- if ( collapse ) names(pairing) else nrm.names
  outdat <- dat[,match.nona(colnames(NRM),sampleNames(dat))]
  exprs(outdat) <- NRM
  return( outdat )
}
if ( FALSE )
{
  CBMGIT <- Sys.getenv('CBMGIT')
  if (CBMGIT=="") stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )
  CBMMLAB <- Sys.getenv('CBMMLAB')
  if (CBMMLAB=="") stop( "Use 'setenv CBMMLAB ..' to set CBMrepository's base directory" )
  source( paste(CBMGIT, "scripts/R/CBMRtools/R/broad.file.formats.R", sep="/") )
  source( paste(CBMGIT, "scripts/R/CBMRtools/R/GeneSet.R", sep="/") )
  require(Biobase)
  
  ## CREATION OF 'TOY' DATA
  ##
  DAT <- readRDS('~/Research/Projects/oralcancer/tcga/firehose_2014_12_06/TCGA_OSCC_mRNA.annotated.rds')
  DAT1 <- DAT[pData(featureData(DAT))[,'symbol']!='',]
  DAT1 <- DAT1[match(unique(pData(featureData(DAT1))[,'symbol']),pData(featureData(DAT1))[,'symbol']),]
  featureNames(DAT1) <- pData(featureData(DAT1))[,'symbol']

  GS <- GeneSet(paste(CBMMLAB,'/annot/c2.cp.v4.0.symbols.gmt',sep=''))
  GS1 <- GS
  setGeneSet(GS1) <- getGeneSet(GS1)[unlist(sapply(c('catenin','emt'),grep,names(getGeneSet(GS1)),ignore.case=TRUE))]
  sapply( lapply(lapply(getGeneSet(GS1),toupper),intersect,toupper(featureNames(DAT1))), length)
  geneSetName(GS1) <- '4pathways'

  gsp.eSet <- DAT1[intersect(unique(unlist(getGeneSet(GS1))),featureNames(DAT1)),]
  gsp.GeneSet <- GS1

  save(gsp.eSet,gsp.GeneSet,file=paste(CBMGIT, "scripts/R/CBMRtools/data/gspData.rda",sep="/"))

  INP <- paste(CBMMLAB,"/annot/c2.cp.v4.0.symbols.gmt",sep="")
  OUT <- paste(CBMGIT, "/scripts/R/CBMRtools/data/gsp.GeneSet.gmt",sep="")
  system( paste("grep CATENIN", INP, ">", OUT) )
  system( paste("grep _EMT_", INP, ">>", OUT) )
  
  ## EXAMPLE OF USE
  ##
  ## the very basic steps are:
  ##
  ## 1) load the data (expression and genesets)
  ## 2) rename the dataset by replacing probesetIDs by gene symbols ==> DB
  ## 3) create the list of lists pairing ==> PAIRS
  ## 4) run geneSetProjection(dat=DB,pairing=PAIRS,gset.db=GSET, ...)
  ##
  CBMGIT <- Sys.getenv('CBMGIT')
  if (CBMGIT=="") stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )
  source( paste(CBMGIT, "scripts/R/CBMRtools/R/misc.R", sep="/") )
  source( paste(CBMGIT, "scripts/R/CBMRtools/R/misc.math.R", sep="/") )
  source( paste(CBMGIT, "scripts/R/CBMRtools/R/broad.file.formats.R", sep="/") )
  source( paste(CBMGIT, "scripts/R/CBMRtools/R/GeneSet.R", sep="/") )
  source( paste(CBMGIT, "scripts/R/CBMRtools/R/ks.score.R", sep="/") )
  source( paste(CBMGIT, "scripts/R/CBMRtools/R/geneSetProjection.R", sep="/") )
  #source( paste(CBMGIT, "scripts/R/CBMRtools/R/heatmap.ggplot.R", sep="/") )
  require(Biobase)
  require(RColorBrewer)
  require(ggdendro)
  require(ggplot2)
  require(grid)
  require(gridExtra)
  require(gtable)
  require(reshape2)
  require(scales)
  require(stats)
  require(cba)
  
  CBMMLAB <- Sys.getenv('CBMMLAB')
  if (CBMMLAB=="") stop( "Use 'setenv CBMMLAB ..' to set CBMrepository's base directory" )
  source( paste(CBMMLAB, "R/heatmap.R", sep="/") )
  
  setwd( paste(CBMGIT, "/scripts/R/CBMRtools/", sep="") )

  ## 1) LOAD the data
  ##
  load("data/gspData.rda")
  if (is.null(gsp.eSet)) stop("is.null(gsp.eSet)")
  if (is.null(gsp.GeneSet)) stop("is.null(gsp.GeneSet)")

  ## 2) RENAME the dataset rows (w/ gene symbols)
  ##
  DAT1 <- gsp.eSet[pData(featureData(gsp.eSet))[,"symbol"]!="",]
  featureNames(DAT1) <- toupper(pData(featureData(DAT1))[,"symbol"])

  ## 3) CREATE the list of lists pairing
  ##
  PAIRS <- list(OSCC=list(control=sampleNames(DAT1)[pData(DAT1)[,"tissue_type"]=="AN"],
                          treatment=sampleNames(DAT1)[pData(DAT1)[,"tissue_type"]!="AN"]))

  ## 4) RUN geneSetProjection ..
  ##
  GSPdir <- geneSetProjection(dat=DAT1,
                              pairing=PAIRS,
                              GS=gsp.GeneSet,
                              collapse=FALSE,
                              weighted=FALSE,
                              absolute=FALSE, # keep the directionality of the enrichment scores
                              min.gset=5,
                              verbose=TRUE)

  gradeID <- 'my_grade'
  stageID <- 'my_stage'
  short.names <- sapply(gsub("REACTOME_","",featureNames(GSPdir)),function(z)
      paste(unlist(strsplit(z,"_"))[1:3],collapse="_"))
  
 # p2 <- heatmap.ggplot2(eSet=GSPdir,col.clust=TRUE,row.clust=TRUE,col.lab=c(gradeID,stageID),row.lab="",
  #                      heatmap.y.text=FALSE, heatmap.x.text=FALSE,heatmap.colorlegend.name="RNASeq_expression",
   #                     title.text="TCGA BRCA log2 gene set projection",
    #                    col.legend.name=c(gradeID,stageID), 
     #                   row.legend.name="", 
      #                  row.scaling="none",
       #                 z.norm=FALSE, 
        #                cuttree.col=0, cuttree.row=0,
         #               verbose=FALSE, show=TRUE)
}
