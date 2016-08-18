#######################################################################
## function: CBM GSEA 2 QMATRIX
##
## Take a list of outputs from multiple calls to cbmGSEA and generate
## a summary geneset-by-signature matrix (and optionally a heatmap) of
## the genesets' q-values
##
#######################################################################
cbmGSEA2qmatrix <- function
(
    cgsea,           # list of data.frames, the up and down genesets, as output by cbmGSEA
    fdr=c(.05,.01),  # FDR thresholds (must be in decreasing order)
    method=c("union","intersect"),
                     # union or intersection of genesets across signatures
    do.sort=TRUE,    # sort matrices by HC
    do.heat=FALSE,   # display heatmap
    rm.zero=TRUE,    # remove genesets/rows w/ no hits
    pvalID="fdr",    # which significance measure to use (either "p2" or "fdr")
    globalMHT=FALSE, # correct for MHT *across* signatures (default is within)
    na.col="gray",   # color for missing values in the heatmap
    verbose=TRUE,    # verbose output
    ...              # extra arguments to my.heatmap
)
{
    ## each element of the variable cgsea corresponds to a distinct cbmGSEA
    ## run, and it contains a data.frame object

    ## input checks
    if ( length(fdr)>1 && any(diff(fdr)>0) ) stop( "fdr must be in decreasing order" )
    if ( globalMHT ) pvalID <- "p2"
    if ( !(pvalID %in% colnames(cgsea[[1]])) ) stop( "unrecognized pvalID: ", pvalID)
    method <- match.arg(method)
    combineFun <- match.fun(method)

    VERBOSE(verbose,"Generating qmatrix based on",method,"..")
    
    ## extract all the geneset IDs in common among all cbmGSEA runs
    gID <- Reduce(combineFun,lapply(cgsea,rownames))

    ## extract the FDRs of all genesets for each cbmGSEA run 
    mx <- matrix(NA,length(gID),length(cgsea),dimnames=list(gID,names(cgsea)))
    for ( i in 1:ncol(mx) ) {
        z <- cgsea[[i]]
        z[z[,"score"]<0,pvalID] <- -z[z[,"score"]<0,pvalID]
        if ( method=="intersect" ) 
            mx[,i] <- z[match.nona(rownames(mx),rownames(z)),pvalID]
        else if ( method=="union" )
            mx[match.nona(rownames(z),rownames(mx)),i] <- z[,pvalID]
        else
            stop("unrecognized method:",method)
    }    
    VERBOSE(verbose, " done, [", paste(dim(mx),collapse=","),"] matrix generated.\n",sep="")

    ## if global multiple hypothesis correction (MHT), take the uncorrected
    ## p-values across signatures and carry out a global FDR correction
    ##
    if ( globalMHT ) {
        VERBOSE(verbose, "carrying out global MHT ..")
        absMX <- abs(mx)
        absMX[,] <- p.adjust(as.vector(absMX),method="BH")
        mx[is.na(mx)] <- 999 # remove NA's (necessary for 'mx<0' test)
        absMX[mx<0] <- -absMX[mx<0]
        mx <- absMX
        VERBOSE(verbose, " done, min(fdr):", min(abs(mx),na.rm=TRUE), "max(fdr):", max(mx,na.rm=TRUE),"\n")
    }
    return( qmatrix2heatmap(mx=mx,fdr=fdr,do.sort=do.sort,do.heat=do.heat,rm.zero=rm.zero,na.col=na.col,...) )
}
#######################################################################
## function: GSEA 2 QMATRIX
##
## Take a list of outputs of multiple calls to runGSEA or
## runGSEApreranked and generate a summary geneset-by-signature matrix
## (and optionally a heatmap) of the genesets' q-values
##
#######################################################################
gsea2qmatrix <- function
(
    gsea,               # list of lists (each with 2 data.frames, the up and down genesets, as output by gsea)
    fdr=c(.05,.01),     # FDR thresholds (must be in decreasing order)
    method=c("union","intersect"),
                        # union or intersection of genesets across signatures
    do.sort=TRUE,       # sort matrices by HC
    do.heat=FALSE,      # display heatmap
    rm.zero=TRUE,       # remove genesets/rows w/ no hits
    pvalID="FDR q-val", # which significance measure to use (see GSEA output for choices)
    na.col="gray",      # color for missing values in the heatmap
    verbose=TRUE,       # verbose output
    ...                 # extra arguments to qmatrix2heatmap
)
{
    ## each element of the list gsea corresponds to a distinct GSEA
    ## run, and it contains a two-element list corresponding to the
    ## files 'gsea_report_for_na_{pos,neg}_*.xls' output by GSEA
    ##
    if ( length(fdr)>1 && any(diff(fdr)>0) ) stop( "fdr must be in decreasing order" )
    if ( !(pvalID %in% colnames(gsea[[1]][[1]])) ) stop( "unrecognized pvalID: ", pvalID)
    if ( !all(sapply(gsea,function(Z1) sapply(Z1,function(Z2) is.character(Z2[,'NAME'])))) )
        stop( "column 'NAME' expected to be of type character")
        
    method <- match.arg(method)
    combineFun <- match.fun(method)

    VERBOSE(verbose,"Generating qmatrix based on",method,"..")

    ## extract all the geneset IDs in common among all GSEA runs
    gID <- Reduce(combineFun,lapply(gsea,function(Z) c(Z[[1]][,'NAME'],Z[[2]][,'NAME'])))

    ## extract the FDRs of all genesets for each GSEA run (i.e., each list item)
    mx <- matrix(NA,length(gID),length(gsea),dimnames=list(gID,names(gsea)))
    for ( i in 1:ncol(mx) ) {
        z <- gsea[[i]]
        tmp <- rbind(z[[1]][,c('NAME',pvalID),drop=FALSE],
                     data.frame(z[[2]][,'NAME',FALSE],-z[[2]][,pvalID,FALSE],
                                check.names=FALSE,stringsAsFactors=FALSE))
        if ( method=="intersect" ) 
            mx[,i] <- tmp[match.nona(rownames(mx),tmp[,'NAME']),pvalID]
        else if ( method=="union" )
            mx[match.nona(tmp[,'NAME'],rownames(mx)),i] <- tmp[,pvalID]
        else
            stop("unrecognized method:",method)
    }
    return( qmatrix2heatmap(mx=mx, fdr=fdr, do.sort=do.sort, do.heat=do.heat, rm.zero=rm.zero, na.col=na.col,...) )
}
#######################################################################
## function: QMATRIX 2 HEATMAP
##
## (Function called by both gsea2qmatrix and cbmGSEA2qmatrix)
## Take a geneset-by-signature matrix of q-values, and generate:
## 1) a 'discretized' matrix based on the input threshold(s)
## 2) (optionally) a color-coded heatmap
#######################################################################
qmatrix2heatmap <- function
(
    mx,             # geneset-by-signature matrix of q-values
    fdr=c(.05,.01), # FDR thresholds (must be in decreasing order)
    do.sort=TRUE,   # sort matrices by HC
    do.heat=FALSE,  # display heatmap
    rm.zero=TRUE,   # remove genesets/rows w/ no hits
    verbose=TRUE,   # extra arguments to my.heatmap
                    # hclust methods
    method=c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"),
    na.col="gray",  # color for missing values in the heatmap
    ...
)
{
    method <- match.arg(method)
    
    levs <- 0:length(fdr)
    zero <- -0.000001
    mx01 <- suppressWarnings(matrix(cut(as.vector(mx),
                                        breaks=c(-1,-fdr,zero,rev(fdr),1),
                                        labels=as.numeric(c(-levs,rev(levs))),include.lowest=TRUE),
                                    nrow=nrow(mx),ncol=ncol(mx)))
    mx01 <- apply(mx01,2,as.numeric)
    dimnames(mx01) <- dimnames(mx)

    if ( rm.zero ) {
        rm.idx <- apply(mx01!=0,1,any,na.rm=TRUE)
        if ( sum(rm.idx)==0 ) {
            VERBOSE(verbose,"no significant genesets found, not removing any")
        }
        else {
            mx01 <- mx01[rm.idx,,drop=FALSE]
            mx <- mx[rm.idx,,drop=FALSE]
            VERBOSE(verbose, "Removed",sum(!rm.idx),"non-significant genesets.\n")
        }
    }
    ## sort rows and columns by HC
    if ( do.sort || do.heat ) {
        hc.col <- hcopt(dist(t(mx01),method="euclidean"),method=method)
        hc.row <- hcopt(dist(mx01,method="euclidean"),method=method)

        if ( do.heat ) {
            ncolors <- length(fdr)*2+1
            COL <- col.gradient(c("blue","white","red"),length=ncolors)
            
            ## yucky handling of color-coding (but necessary 'cause the default is not smart enough)
            ## (if anybody has a better solution, feel free to amend)
            COL <- COL[sort(unique(as.vector(mx01)))+levs[length(levs)]+1]
                     
            ## add color used for the NA's
            if ( any(is.na(mx)) ) {
                COL <- c(COL,col.gradient(na.col,1))
                mx01[is.na(mx01)] <- max(mx01,na.rm=TRUE)+1
            }            
            my.heatmap(mx01,Rowv=as.dendrogram(hc.row),Colv=as.dendrogram(hc.col),
                       scale="none",col=COL,revC=TRUE,...)
        }
        if ( do.sort ) {
            mx <- mx[hc.row$order,hc.col$order]
            mx01 <- mx01[hc.row$order,hc.col$order]
        }
    }
    return( list(mx=mx,mx01=mx01) )
}
#######################################################################
## function: HYPER 2 QMATRIX
##
## Take the output of hyperEnrichment (called on multiple signatures)
## and generate a summary geneset-by-signature matrix (and optionally
## a heatmap) of the genesets' q-values
#######################################################################
hyper2qmatrix <- function
(
    hyper,          # output of hyperEnrichment
    fdr=c(.05,.01), # FDR thresholds (must be in decreasing order)
    do.sort=TRUE,   # sort matrices by HC
    do.heat=FALSE,  # display heatmap
    rm.zero=TRUE,   # remove genesets/rows w/ no hits
                    # hclust methods
    method=c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"),
    ...             # extra arguments to my.heatmap
)
{
    if ( length(fdr)>1 && any(diff(fdr)>0) )
        stop( "fdr must be in decreasing order" )
    if ( length(unique(hyper[,"set"]))==1 )
        stop( "hyper must contain results for at least two signatures")
    method <- match.arg(method)
    
    SIG <- unique(hyper[,"set"])
    gsets <- hyper[hyper[,"set"]==SIG[1],"category"]
    mx <- sapply(SIG,function(z) {
        tmp <- hyper[hyper[,"set"]==z,c("fdr","category")]
        as.numeric(tmp[match(gsets,tmp[,"category"]),"fdr"])
    });
    rownames(mx) <- gsets
    mx01 <- mx
    for ( i in 1:length(fdr) ) {
        mx01[mx<=fdr[i]] <- i
    }
    mx01[mx>fdr[1]] <- 0

    if ( rm.zero ) {
        keep.idx <- apply(mx01!=0,1,any)
        if ( sum(keep.idx)==0 ) {
            warning("no significant genesets found")
        }
        else {
            mx01 <- mx01[keep.idx,,drop=FALSE]
            mx <- mx[keep.idx,,drop=FALSE]
        }
    }
    ## sort rows and columns by HC
    if ( do.sort || do.heat ) {
        hc.col <- hcopt(dist(t(mx01),method="euclidean"),method="ward.D")
        hc.row <- hcopt(dist(mx01,method="euclidean"),method="ward.D")

        if ( do.heat ) {
            ncolors <- length(unique(as.vector(mx01)))
            my.heatmap(mx01,Rowv=as.dendrogram(hc.row),Colv=as.dendrogram(hc.col),scale="none",
                       col=col.gradient(c("white","red"),length=ncolors),revC=TRUE,...)
        }
        if ( do.sort ) {
            mx <- mx[hc.row$order,hc.col$order]
            mx01 <- mx01[hc.row$order,hc.col$order]
        }
    }
    return( list(mx01=mx01,mx=mx) )
}

##################################################################
## code to create and test the data objects used in the examples #
##################################################################
if ( FALSE )
{
    require(CBMRtools)
    CBMRtools <- Sys.getenv("CBMRtools")
    source(file.path(CBMRtools,"dvlp/heatmap.R"))
    source(file.path(CBMRtools,"dvlp/enrich2qmatrix.R"))
    source(file.path(CBMRtools,"dvlp/diffanalOverlap.R"))
    ## path in my desktop
    PATH <- "~/Research/Projects/AhR/AHR_CYP1B1_CB799_microarrays2016"
    GSPATH <- "~/Research/CancerGenomeAnalysisData/annot"
    
    ## paths on SCC:
    #PATH <- "/restricted/projectnb/montilab-p/projects/AhR/AHR_CYP1B1_CB799_microarrays2016"
    #GSPATH <- "/restricted/projectnb/montilab-p/CBMrepositoryData/annot/"

    ## upload MSigDB's hallmark geneset compendium (stored locally)
    gSet <- new("GeneSet",file.path(GSPATH,"h.all.v5.1.symbols.gmt"))

    ## loading pre-computed list of diffanal results
    DIF <- readRDS(file.path(PATH,"results/rds/DIF2.RDS"))
    names(DIF)
    ##  [1] "MDA.AhR" "MDA.CYP" "SUM.AhR" "SUM.CYP" "AhR"     "CYP"     "CB113"  
    ##  [8] "BaP"     "PYO"     "CH"     

    ## loading pre-computed list of cbmGSEA results (this takes a long time to compute)
    cGSEA <- readRDS(file.path(PATH,"results/rds/cbmGSEA.RDS"))
    names(cGSEA)
    ## [1] "h.all"  "c2.cp"  "c3.all" "c2.cgp"

    ## extract the hallmarks results (hGSEA is a list of data.frames, one per signature)
    hGSEA <- cGSEA[['h.all']]
    names(hGSEA)
    ## [1] "MDA.AhR" "MDA.CYP" "SUM.AhR" "SUM.CYP" "AhR"     "CYP"     "CB113"  
    ## [8] "BaP"     "PYO"     "CH"     

    png("ahr.hallmarks.png")
    OUT <- cbmGSEA2qmatrix(hGSEA,fdr=c(0.10,0.05,0.01),do.heat=TRUE,globalMHT=TRUE,margins=c(6,15))
    dev.off()

    ## let's test with different genesets
    hGSEA1 <- hGSEA
    for ( i in 1:length(hGSEA1) ) # remove a (different) geneset from each data.frame
        hGSEA1[[i]] <- hGSEA[[i]][-i,]

    ## take the intersection
    OUT1 <- cbmGSEA2qmatrix(hGSEA1,fdr=c(0.10,0.05,0.01),do.heat=TRUE,globalMHT=TRUE,margins=c(6,15),method="intersect")

    ## take the union (there should be rows with missing values, color-coded 'gray' by default)
    OUT2 <- cbmGSEA2qmatrix(hGSEA1,fdr=c(0.10,0.05,0.01),do.heat=TRUE,globalMHT=TRUE,margins=c(6,15),method="union")

    ## TEST HYPER-ENRICHMENT VISUALIZATION

    SIG <- diffanal2signatures(DIF, maxFDR=0.05, minFC=1.5)
    # notice it will return a warning re 'multiple matches'. The latest version of matchIndex eliminates it.
    
    gSet <- new("GeneSet","~/Research/CancerGenomeAnalysisData/annot/h.all.v5.1.symbols.gmt")

    hOut <- hyperEnrichment(SIG,categories=getGeneSet(gSet),ntotal=22981)

    tmp <- hyper2qmatrix(hOut,fdr=c(.05,.01),do.heat=TRUE,margins=c(9,20))

    my.write.table(tmp$mx,names="GeneSet",file="HYPERtable.xls")
}
