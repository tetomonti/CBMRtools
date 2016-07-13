#####################################################################################
#' cbmGSEA
#' 
#' \code{cbmGSEA} re-implementation of Broad's GSEA, with
#' gene-specific permutation-based p-value calculation and FDR
#' correction by the BH method
#'
#' @param eSet                 an [m-genes x n-samples] matrix
#' @param rnkScore             provide a ranking score instead
#' @param pheno                class ID
#' @param tag                  id of gene to use as template (nearest neighbor mode)
#' @param cond                 treatment label (e.g., TCDD)
#' @param control              control label (e.g., DMSO)
#' @param geneId               fData column containing the gene symbols
#' @param gSet                 a named list of genesets (i.e. vectors of gene IDs)
#' @param minGset              minimum length of accepted genesets
#' @param nperm                number of permutation iterations
#' @param score                one of {"t.score","delta"}: differential score used when cls is specified
#' @param do.abs               rank genes by absolute value of score
#' @param robust               compute robust score
#' @param method               one of {"pearson","spearman","euclidean"}: score used when tag is specified
#' @param alternative          one of {"two.sided","greater","less"}
#' @param weighted             use weighted KS score as in GSEA
#' @param weight.p             weight's exponent
#' @param smoother             smooth p-values
#' @param clsLev               phenotype labels (e.g., c('normal','tumor'))
#' @param confound             confounder variable (must be categorical)
#' @param seed                 random seed (for reproducible results)
#' @param do.pval              compute asymptotic p-values
#' @param exhaustive           generate all possible permutations
#' @param do.plot              generate ks plots 
#' @param plot.name            file stub where to save ks plots
#' @param plot.dev             display device
#' @param topG                 number of genesets to display per direction
#' @param verbose              verbosity TRUE/FALSE
#'
#' @return a data.frame with one geneset per row and columns reporting KS score, p-value, q-value, etc.
#'
#' @examples
#'
#' data(cbmGSEAdata)
#' if ( is.null(cbmGSEA.eSet) || is.null(cbmGSEA.gSet) ) stop( "missing objects" )
#'
#' table(pData(cbmGSEA.eSet)$tissue_type) # expecting 8 AN's + 8 Tumor's
#'
#' # carry out a simple t-score calculation (notice the '-', to achieve a decreasing=TRUE sort)
#' rnkS <- -t.score(x=exprs(cbmGSEA.eSet),cls=factor(pData(cbmGSEA.eSet)$tissue_type))
#'
#' # 1) phenotype shuffling
#'
#' set.seed(123)
#' OUT1 <- cbmGSEA(eSet=cbmGSEA.eSet,gSet=cbmGSEA.gSet,pheno="tissue_type",cond="Tumor",control="AN",
#'                 nperm=100,verbose=TRUE,do.plot=TRUE,topG=1)
#' print(OUT1[,1:5])
#' 
#' # 2) rank shuffling (based on pre-computed ranked list)
#'
#' set.seed(123)
#' OUT2 <- cbmGSEA(rnkScore=rnkS,gSet=cbmGSEA.gSet,nperm=100,verbose=TRUE,do.plot=TRUE,topG=1)
#' print(OUT2[,1:5])
#' all.equal(OUT1[,"score"],OUT2[,"score"]) # same results expected
#'
#' # 3) gene tag shuffling (ranking by nearest neighbor)
#'
#' set.seed(123)
#' OUT3 <- cbmGSEA(eSet=cbmGSEA.eSet,tag="MAML2",gSet=cbmGSEA.gSet,nperm=100,verbose=TRUE,do.plot=TRUE,topG=1)
#' print(OUT3[,1:5])
#' 
#' @export
#' 
#####################################################################################
cbmGSEA <- function
(
 eSet=NULL,                  # an (m-genes x n-samples) matrix
 rnkScore=NULL,              # provide a ranking score instead
 pheno=NULL,                 # class ID
 tag=NULL,                   # id of gene to use as template (nearest neighbor mode)
 cond=NULL,                  # 'condition' label
 control=NULL,               # 'control' label
 geneId=NULL,                # fData column containing the gene symbols
 gSet,                       # a named list of genesets (i.e., vectors of gene IDs)
 minGset=5,                  # minimum length of accepted genesets
 nperm=100,                  # number of permutation iterations
 score=c("t.score","delta"), # differential score used when cls is specified
 do.abs=FALSE,               # rank genes by absolute value of score
 robust=FALSE,               # compute robust score
 method=c("pearson","spearman","euclidean"),
                             # score used when tag is specified
 alternative=c("two.sided","greater","less"),
 weighted=FALSE,             # use weighted KS score as in GSEA
 weight.p=1,                 # weight's exponent
 smoother=1,                 # smooth p-values
 clsLev=c("neg","pos"),      # phenotype labels (e.g., c('normal','tumor'))
 confound=NULL,              # confounder variable (must be categorical)
 seed=NULL,                  # random seed (for reproducible results)
 do.pval=!weighted,          # compute asymptotic p-values
 exhaustive=FALSE,           # generate all possible permutations
 do.plot=!is.null(plot.name),# generate ks plots 
 plot.name=NULL,             # file stub where to save ks plots
 plot.dev=pdf,
 topG=20,                    # max number of genesets to display per direction
 verbose=FALSE,
 ...
)
{
    ## BEGIN checks on inputs
    ##
    if ( !xor(is.null(eSet),is.null(rnkScore)) )
        stop( "only one of {eSet,rnkScore} can be specified" )
    if ( class(gSet)!="list" ) stop( "gSet must be a list: ", class(gSet) )
    if ( is.null(names(gSet)) ) stop( "gSet must ba a *named* list" )
    if ( is.null(rnkScore) ) {
        if ( !xor(is.null(pheno),is.null(tag)) )
            stop( "only one of {pheno, tag} must be specified" )
        if ( !is.null(geneId) && !(geneId %in% colnames(fData(eSet))) )
            stop( "geneId must be a valid column name in fData(eSet)" )
        if ( !is.null(pheno) && (is.null(cond) || is.null(control)) )
            stop( "must specify cond and control when using pheno" )
        if ( nrow(eSet)<100 )
            warning( "less than 100 genes in the eSet:", nrow(eSet) )
    }
    else {
        if ( is.null(names(rnkScore)) ) stop( "rnkScore must be annotated with gene symbols" )
    }
    score <- match.arg(score)
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    ##
    ## END checks

    alt.dir <- if (alternative=="greater") 1 else -1
    Ndat <- NULL
    
    ## ExpressionSet supplied in input
    ##
    if ( !is.null(eSet) )
    {    
        Ndat <- nrow(eSet)
        SCORE <- {
            if ( is.null(tag) )
                match.fun(score)
            else
                switch(method,
                       euclidean = function(x,y){sqrt(sum((x-y)^2))},
                       pearson =   function(x,y){1 - cor(x,y)},
                       spearman =  function(x,y){1 - cor(rank(x),rank(y))})
        }
        if ( !is.null(pheno) ) {
            condIdx <- which(pData(eSet)[,pheno]==cond);
            if (length(condIdx)<1) stop( sprintf("no '%s' samples",condition) )
            cntlIdx <- which(pData(eSet)[,pheno]==control);
            if (length(cntlIdx)<1) stop( sprintf("no '%s' samples",control) )
            #eDat <- eDat[,c(cntlIdx,condIdx)]
            eSet <- eSet[,c(cntlIdx,condIdx)]
            cls <- factor(pData(eSet)[c(cntlIdx,condIdx),pheno],levels=c(control,cond))
            clsLev <- levels(cls)
        }
        else if ( !is.null(tag) ) {
            tagI <- matchIndex(tag,featureNames(eSet))
            tag <- exprs(eSet)[tagI,]
            eSet <- eSet[-tagI,]
        }
        eDat <- exprs(eSet)

        ## compute observed ranking scores
        ##
        VERBOSE( verbose, "  computing observed scores .. " )
        rnkScore <- {
            if ( is.null(pheno) ) {
                VERBOSE( verbose, "(based on NN) .. ")
                apply(eDat,1,SCORE,y=tag)
            }
            else {
                VERBOSE( verbose, "(based on class template) .. ")
                # this is to have t.score and similar agree with run_limma
                -SCORE(eDat, cls=cls, robust=robust)
            }
        }
        VERBOSE( verbose, "done.\n" )
        geneSymbols <- if (is.null(geneId)) rownames(eDat) else fData(eSet)[,geneId]
    }
    ## ELSE: ranking score supplied in input
    ##
    else {
        Ndat <- length(rnkScore)
        geneSymbols <- names(rnkScore)
    }
    if (do.abs) {
        rnkScore <- abs(rnkScore)
    }
    ## Map gene symbols to indices
    gsetIdx <- glist2idx( gSet, gnames=geneSymbols, minGset=minGset )
    ghitLen <- sapply( gsetIdx, length )
    gsetLen <- sapply( gSet[match(names(gsetIdx),names(gSet))], length )
    ngset <- if (is.list(gSet)) length(gSet) else 1

    ## rank genes by scores (decreasing==TRUE is achieved by taking -score)
    geneRanks <- rank(-rnkScore,ties.method="first")

    ## compute observed KS scores
    VERBOSE( verbose, "  computing observed enrichment score(s) .." )
    ksObs <- t(sapply(names(gsetIdx), function(z)
        ksGenescore(Ndat, y=geneRanks[gsetIdx[[z]]], do.pval=do.pval, bare=TRUE, 
                    weight=if (weighted) sort(rnkScore), weight.p=weight.p,
                    do.plot=FALSE, cls.lev=clsLev, main=z,...)))
    if (!do.pval) {
        ksObs <- t(ksObs)
        colnames(ksObs) <- "score"
    }
    if ( do.plot ) {
        if ( !is.null(plot.name) ) {
            plot.dev( plot.name )
        }
        tmp <- ksObs[order(ksObs[,"score"],decreasing=TRUE),,drop=FALSE]
        topG1 <- min(topG,sum(ksObs[,"score"]>0))
        topG2 <- min(topG,sum(ksObs[,"score"]<0))
        gsetIdx1 <- c(rownames(tmp)[1:topG1],
                      rev(rownames(tmp))[1:topG2])
        if ( any(!(gsetIdx1 %in% names(gsetIdx))) ) stop( "gsetIdx1 %in% names(gsetIdx)" )
        tmp <- t(sapply(gsetIdx1, function(z)
            ksGenescore(Ndat, y=geneRanks[gsetIdx[[z]]], do.pval=FALSE, bare=TRUE, 
                        weight=if (weighted) sort(rnkScore), weight.p=weight.p,
                        do.plot=TRUE, cls.lev=clsLev, main=z,...)))
        if ( !is.null(plot.name) ) {
            dev.off()
            VERBOSE( verbose, "\n  (plots saved to file '", plot.name, "')\n", sep="" )
        }
    }
    VERBOSE( verbose, " done.\n" )

    ## return scores if no permutation requested
    if (nperm<1) return(ksObs)

    ## PERMUTATION-BASED SIGNIFICANCE TESTING ..

    VERBOSE( verbose, "  permutation testing ")
    ## 1) by shuffling the phenotype
    ##
    if ( !is.null(eSet) && !is.null(pheno) ) {
        clsPerm <- permute.binarray(cls,nperm=nperm,control=confound,exhaustive=exhaustive,
                                    verbose=verbose)
        if ( nperm>nrow(clsPerm) ) smoother=0 # case when exhaustive is less than required
        nperm <- nrow(clsPerm)
        permScores <- apply(clsPerm, 1, function(Z) SCORE(eDat, cls=Z, robust=robust))
        VERBOSE( verbose, "(by phenotype shuffling):\n" )
    }
    ## 2) by shuffling the reference gene tag
    ##
    else if ( !is.null(eSet) && !is.null(tag) ) {
        clsPerm <- sapply(1:nperm, function(z) sample(tag))
        permScores <- apply(clsPerm, 2, function(y) t(apply(eDat,1,SCORE,y=y)))
        VERBOSE( verbose, "(by tag shuffling):\n" )
    }
    ## 3) by shuffling the preranked input list
    ##
    else if ( !is.null(rnkScore) ) {
        permScores <- sapply(1:nperm, function(z) sample(rnkScore))
        VERBOSE( verbose, "(by rank shuffling):\n" )
    }
    else stop( "unexpected permutation option" )

    ## run ksGenescore on the permuted scores ..
    ##
    VERBOSE( verbose, "\tcomputing permuted ks scores .. " )
    percent <- pctstep <- max( .1, round(1/nperm,2) )

    xPerm <- matrix( 0, length(gsetIdx), length(perm.2side.names),
                     dimnames=list(names(gsetIdx),perm.2side.names) )
 
    for ( i in 1:nperm )
    {
        permScore <- permScores[,i]
        if (do.abs) permScore <- abs(permScore)
        scoreRank <- rank(permScore,ties.method="first")
        weight <- if (weighted) sort(permScore)
        
        ksNul <- sapply(gsetIdx, function(z)
            ksGenescore(Ndat,y=scoreRank[z],weight=weight,weight.p=weight.p,do.pval=FALSE) )

        tmp <- perm.2side.online(ksObs[,1], ksNul, dir=alt.dir)
        if ( !all(dim(tmp)==dim(xPerm)) ) stop("incompatible dimensions for summary")
        xPerm <- xPerm + tmp
        
        if ( verbose & i>=nperm*percent ) {
            VERBOSE( verbose, percent*100,"% ", sep="" )
            percent <- round(percent + pctstep,1)
        }
    }
    VERBOSE( verbose, "\n" )

    ## complete computation of some p-values (see perm.2side.summary 
    ## ..for the "offline" computation of these statistics

    OUT <- p2ss.add( x.prm=xPerm, x.obs=ksObs[,1], nperm=nperm, smoother=smoother )

    bind.fun <- if (ngset>1) cbind else c
    names(ghitLen) <- names(gsetLen) <- NULL
    if ( do.pval ) OUT <- bind.fun(OUT,
                                   asymptotic.p=signif(ksObs[,2]),
                                   asymptotic.fdr=signif(pval2fdr(ksObs[,2])))
    OUT <- bind.fun( OUT, hits=ghitLen, size=gsetLen )
    
    return( OUT )
}
## FUNCTION: GLIST 2 IDX
##
glist2idx <- function( glist, gnames, minGset, verbose=TRUE )
{
  if ( !is.list(glist) ) stop( "glist must be a list: ",class(glist) )
      
  ## replace list of geneset names with list of geneset indices
  ##
  glistIdx <- lapply( glist, match.nona, gnames, na.rm=TRUE, suppressWarnings=TRUE )
  glistLen <- sapply( glistIdx, length )
  glistIdx <- glistIdx[glistLen>=minGset]

  if ( (Nrm <- sum(glistLen<minGset))>0 ) {
      VERBOSE( verbose, "  Removed ", Nrm, " genesets because too short (<", minGset,")\n", sep="" )
  }
  if ( length(glistIdx)==0 )
      stop( "no geneset with the required minimum of genetags present in the dataset" )

  return( glistIdx )
}
##################################################################
## code to create and test the data objects used in the examples #
##################################################################
if ( FALSE )
{
    require(Biobase)
    require(CBMRtools)
    CBMRtools <- Sys.getenv("CBMRtools")
    source(paste(CBMRtools,"CBMRtools/R/diffanal.scores.R",sep="/"))
    source(paste(CBMRtools,"CBMRtools/R/ks.score.R",sep="/"))
    source(paste(CBMRtools,"CBMRtools/R/misc.R",sep="/"))
    source(paste(CBMRtools,"CBMRtools/R/perm.2side.R",sep="/"))
    source(paste(CBMRtools,"CBMRtools/R/permute.array.R",sep="/"))
    source(paste(CBMRtools,"CBMRtools/R/cbmGSEA.R",sep="/"))
    
    gSet <- getGeneSet(new("GeneSet","~/Research/CancerGenomeAnalysisData/annot/h.all.v5.1.symbols.gmt"))
    brca <- readRDS("~/Research/Projects/tcga/brca/BRCA_2015_02_04_ES.rds")
    exprs(brca) <- log2(exprs(brca)+1)
    
    ## prepare GEP data (extract small subset of samples and subset of genes)
    
    brca <- brca[!is.na(fData(brca)$gene_symbol),]
    featureNames(brca) <- fData(brca)$gene_symbol
    gSet <- gSet[order(sapply(gSet,length))[1:3]]
    set.seed(123)
    brca <- brca[,c(sample(which(pData(brca)$my_stage=="AN"),size=8),
                    sample(which(pData(brca)$my_stage=="stage i"),size=8))]
    brca <- brca[!apply(exprs(brca)==0,1,all),]
    brca1 <- brca[c(unique(unlist(gSet)),
                    sample(setdiff(featureNames(brca),unique(unlist(gSet))),size=100)),]
    pData(brca1) <- pData(brca1)[,"tissue_type",drop=FALSE]
    
    ## t.score returns negative scores for the cond class and positive for the control class (opposite run_limma)
    rnkS <- -t.score(x=exprs(brca1),cls=factor(pData(brca1)$tissue_type))
    
    ## add two custom "true positive" genesets 
    gSet <- c(list(tumor50=names(rnkS)[order(rnkS,decreasing=TRUE)[1:50]],
                   normal50=names(rnkS)[order(rnkS,decreasing=FALSE)[1:50]]),
              gSet)
    
    set.seed(123)
    OUT1 <- cbmGSEA(eSet=brca1,gSet=gSet,pheno="tissue_type",cond="Tumor",control="AN",nperm=100,verbose=TRUE,plot.name="cbmGSEA.plots1.pdf")
    
    set.seed(123)
    OUT2 <- cbmGSEA(rnkScore=rnkS,gSet=gSet,nperm=100,verbose=TRUE,plot.name="cbmGSEA.plots2.pdf")
    
    all.equal(OUT1[,"score"],OUT2[,"score"])
    
    set.seed(123)
    OUT3 <- cbmGSEA(eSet=brca1,tag="MAML2",gSet=gSet,nperm=100,verbose=TRUE,plot.name="cbmGSEA.plots3.pdf")
    
    cbmGSEA.eSet <- brca1
    cbmGSEA.gSet <- gSet
    save(cbmGSEA.eSet,cbmGSEA.gSet,file=file.path(CBMRtools, "CBMRtools/data/cbmGSEAdata.rda"))
}
