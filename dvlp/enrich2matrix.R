#######################################################################
## function: HYPER 2 MATRIX
##
## Take the output of hyperEnrichment (called on multiple signatures)
## and generate a summary geneset-by-signature matrix (and optionally
## a heatmap) of the genesets' q-values
#######################################################################
hyper2matrix <- function
(
    hyper,          # output of hyperEnrichment
    fdr=c(.05,.01), # FDR thresholds (must be in decreasing order)
    do.sort=TRUE,   # sort matrices by HC
    do.heat=FALSE,  # display heatmap
    rm.zero=TRUE,   # remove genesets/rows w/ no hits
    ...             # extra arguments to my.heatmap
)
{
    if ( length(fdr)>1 && any(diff(fdr)>0) ) stop( "fdr must be in decreasing order" )
    if ( length(unique(hyper[,"set"]))==1 ) stop( "hyper must contain results for at least two signatures")
    
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
        rm.idx <- apply(mx01!=0,1,any)
        if ( sum(rm.idx)==0 ) {
            warning("no significant genesets found")
        }
        else {
            mx01 <- mx01[rm.idx,,drop=FALSE]
            mx <- mx[rm.idx,,drop=FALSE]
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
#######################################################################
## function: GSEA 2 MATRIX
##
## Take a list of outputs of multiple calls to runGSEA or
## runGSEApreranked and generate a summary geneset-by-signature matrix
## (and optionally a heatmap) of the genesets' q-values
##
#######################################################################
gsea2matrix <- function
(
    gsea,               # list of lists (each with 2 data.frames, the up and down genesets, as output by gsea)
    fdr=c(.05,.01),     # FDR thresholds (must be in decreasing order)
    do.sort=TRUE,       # sort matrices by HC
    do.heat=FALSE,      # display heatmap
    rm.zero=TRUE,       # remove genesets/rows w/ no hits
    pvalID="FDR q-val", # which significance measure to use (see GSEA output for choices)
    ...                 # extra arguments to fdrMatrix2heatmap
)
{
    ## each element of the variable gsea correspond to a distinct GSEA
    ## run, and it contains a two-element list corresponding to the
    ## files 'gsea_report_for_na_{pos,neg}_*.xls' output by GSEA
    ##
    if ( length(fdr)>1 && any(diff(fdr)>0) ) stop( "fdr must be in decreasing order" )
    if ( !(pvalID %in% colnames(gsea[[1]][[1]])) ) stop( "unrecognized pvalID: ", pvalID)

    ## extract all the geneset IDs
    gID <- c(gsea[[1]][[1]][,'NAME'],gsea[[1]][[2]][,'NAME'])

    ## extract the FDRs of all genesets for each GSEA run (i.e., each list item)
    mx <- sapply( gsea, function(z) {
        tmp <- rbind(z[[1]][,c('NAME',pvalID)],
                     data.frame(z[[2]][,'NAME',FALSE],-z[[2]][,pvalID,FALSE],
                                check.names=FALSE,stringsAsFactors=FALSE))
        tmp[ match.nona(gID,tmp[,'NAME']), pvalID ]
    })
    rownames(mx) <- gID

    return( fdrMatrix2heatmap(mx=mx, fdr=fdr, do.sort=do.sort, do.heat=do.heat, rm.zero=rm.zero, ...) )
}
#######################################################################
## function: FDR MATRIX 2 HEATMAP
##
## (Function called by both gsea2matrix and cbmGSEA2matrix)
## Take a geneset-by-signature matrix of q-values, and generate:
## 1) a 'discretized' matrix based on the input threshold(s)
## 2) (optionally) a color-coded heatmap
#######################################################################
fdrMatrix2heatmap <- function
(
    mx,             # geneset-by-signature matrix of q-values
    fdr=c(.05,.01), # FDR thresholds (must be in decreasing order)
    do.sort=TRUE,   # sort matrices by HC
    do.heat=FALSE,  # display heatmap
    rm.zero=TRUE,   # remove genesets/rows w/ no hits
    verbose=TRUE,   # extra arguments to my.heatmap
    ...
)
{
    levs <- 0:length(fdr)
    zero <- -0.000001
    mx01 <- suppressWarnings(matrix(cut(as.vector(mx),
                                        breaks=c(-1,-fdr,zero,rev(fdr),1),
                                        labels=as.numeric(c(-levs,rev(levs))),include.lowest=TRUE),
                                    nrow=nrow(mx),ncol=ncol(mx)))
    mx01 <- apply(mx01,2,as.numeric)
    dimnames(mx01) <- dimnames(mx)

    if ( rm.zero ) {
        rm.idx <- apply(mx01!=0,1,any)
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
        hc.col <- hcopt(dist(t(mx01),method="euclidean"),method="ward.D")
        hc.row <- hcopt(dist(mx01,method="euclidean"),method="ward.D")

        if ( do.heat ) {
            ncolors <- length(fdr)*2+1
            COL <- col.gradient(c("blue","white","red"),length=ncolors)
            COL <- COL[sort(unique(as.vector(mx01)))+levs[length(levs)]+1]
            my.heatmap(mx01,Rowv=as.dendrogram(hc.row),Colv=as.dendrogram(hc.col),scale="none",col=COL,revC=TRUE,...)
        }
        if ( do.sort ) {
            mx <- mx[hc.row$order,hc.col$order]
            mx01 <- mx01[hc.row$order,hc.col$order]
        }
    }
    return( list(mx=mx,mx01=mx01) )
}
#######################################################################
## function: CBM GSEA 2 MATRIX
##
## Take a list of outputs from multiple calls to cbmGSEA and generate
## a summary geneset-by-signature matrix (and optionally a heatmap) of
## the genesets' q-values
##
#######################################################################
cbmGSEA2matrix <- function
(
    cgsea,           # list of data.frames, the up and down genesets, as output by cbmGSEA
    fdr=c(.05,.01),  # FDR thresholds (must be in decreasing order)
    do.sort=TRUE,    # sort matrices by HC
    do.heat=FALSE,   # display heatmap
    rm.zero=TRUE,    # remove genesets/rows w/ no hits
    pvalID="fdr",    # which significance measure to use (either "p2" or "fdr")
    globalMHT=FALSE, # correct for MHT *across* signatures (default is within)
    verbose=TRUE,    # verbose output
    ...              # extra arguments to my.heatmap
)
{
    ## each element of the variable gsea correspond to a distinct cbmGSEA
    ## run, and it contains a data.frame object
    ##
    if ( length(fdr)>1 && any(diff(fdr)>0) ) stop( "fdr must be in decreasing order" )
    if ( globalMHT ) pvalID <- "p2"
    
    ## extract all the geneset IDs from the first data.frame
    gID <- rownames(cgsea[[1]])
    if ( !(pvalID %in% colnames(cgsea[[1]])) ) stop( "unrecognized pvalID: ", pvalID)

    ## extract the FDRs of all genesets for each cbmGSEA run 
    VERBOSE(verbose, "extracting into matrix ..")
    mx <- sapply( cgsea, function(z) {
        tmp <- z[match.nona(gID,rownames(z)),,drop=FALSE]
        tmp[tmp[,"score"]<0,pvalID] <- -tmp[tmp[,"score"]<0,pvalID]
        tmp[,pvalID]})
    rownames(mx) <- gID
    VERBOSE(verbose, " done, [", paste(dim(mx),collapse=","),"] matrix generated.\n",sep="")

    ## if global multiple hypothesis correction (MHT), take the uncorrected
    ## p-values across signatures and carry out a global FDR correction
    ##
    if ( globalMHT ) {
        VERBOSE(verbose, "carrying out global MHT ..")
        absMX <- abs(mx)
        absMX[,] <- p.adjust(as.vector(absMX),method="BH")
        absMX[mx<0] <- -absMX[mx<0]
        mx <- absMX
        VERBOSE(verbose, " done, min(fdr):", min(abs(mx)), "max(fdr):", max(mx),"\n")
    }
    return( fdrMatrix2heatmap(mx=mx,fdr=fdr,do.sort=do.sort,do.heat=do.heat,rm.zero=rm.zero, ...) )
}
