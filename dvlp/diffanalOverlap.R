#######################################################################
## function: DIFFANAL OVERLAP
##
## Take a list of diffanal results (each list item being the output of
## limma_wrapper, or edgeR_wrapper, or similar), and return a list of
## correspoding gene lists (selected based on significance criteria),
## and a table with pairwise signature comparison (overlap
## significance and overlap counts).
#######################################################################

diffanalOverlap <- function
(
    diffanalList,      # list of diffanal results (each a data.frame, see e.g. limma_wrapper)
    key="hgnc_symbol", # column id of gene identifiers in the table
    maxFDR=1.0,        # max FDR for genes to include in signature
    minFCG=1.0,        # min fold-change for genes to include in signature
    ntotal=nrow(diffanalList[[1]]),
    fdrKey="adj.P.Val",
    fcgKey="limma.fold.change",
    scrKey="t",
    adjust=TRUE,
    outfile=NULL
)
{
    ## First, extract signatures (i.e., gene lists) from diffanal results
    SIGlist <- diffanal2signatures(diffanalList,key=key,maxFDR=maxFDR,minFCG=minFCG,fdrKey=fdrKey,fcgKey=fcgKey,scrKey=scrKey)
    
    ## Then, compare the signatures (by 'Venn-Diagramming')
    signatureOverlap( SIGlist, ntotal=ntotal, adjust=adjust, outfile=outfile )
}
#######################################################################
## function: DIFFANAL 2 SIGNATURES
##
## Take a list of diffanal results (each list item being the output of
## limma_wrapper, or edgeR_wrapper, or ..), and return a list of
## 'signatures' (each being a list of gene identifiers selected based
## on various criteria).
#######################################################################

diffanal2signatures <- function
(
    diffanalList,               # list of diffanal results (possibly, data frames)
    maxFDR=1.0,                 # max FDR to be included in the signatures
    minFCG=1.0,                 # min fold-change to be included
    key="hgnc_symbol",          # column header where to find gene symbols
    fdrKey="adj.P.Val",         # column IDs (defaults are run_limma columns)
    fcgKey="limma.fold.change", # fold-change key
    scrKey="t"                  # score key
)
{
    SIGlist <- NULL

    ## check inputs
    if ( is.null(names(diffanalList)) ) stop( "is.null(names(diffanalList))" )
    if ( maxFDR>1 | maxFDR<=0 ) stop( "maxFDR must be in (0,1]: ", maxFDR )
    if ( minFCG<1 ) stop( "minFCG must be >= 1: ", minFCG )
    
    ## establish key mapping (and throw error if any mapping fails)
    cnames <- colnames(diffanalList[[1]])
    keyI <- matchIndex(key,cnames)
    fdrI <- matchIndex(fdrKey,cnames)
    fcgI <- matchIndex(fcgKey,cnames)
    scrI <- matchIndex(scrKey,cnames)
    
    for ( i in 1:length(diffanalList) )
    {
        DIFi <- diffanalList[[i]]
        SIGi <- list(up=DIFi[as.numeric(DIFi[,fdrI])<=maxFDR &
                             as.numeric(DIFi[,fcgI])>=minFCG &
                             as.numeric(DIFi[,scrI])>0,keyI],
                     dn=DIFi[as.numeric(DIFi[,fdrI])<=maxFDR &
                             1/as.numeric(DIFi[,fcgI])>=minFCG &
                             as.numeric(DIFi[,scrI])<0,keyI])
        names(SIGi) <- paste( names(diffanalList)[i], names(SIGi), sep="." )
        SIGlist <- c( SIGlist, SIGi )
                  
    }
    return( SIGlist )
}
#######################################################################
## function: SIGNATURE OVERLAP
##
## Take a list of signatures (each being a list of gene identifiers),
## and generate a table with pairwise signature overlap (significance and
## counts).
#######################################################################

signatureOverlap <- function
(
    SIGlist,     # a list of lists of gene identifiers
    ntotal,      # background population size for the hyper test
    adjust=TRUE, # adjust p-values across all signature comparisons
    outfile=NULL # write summary table to output file
)
{
    OVLP <- NOVLP <- JOVLP <- matrix(NA,length(SIGlist),length(SIGlist),dimnames=list(names(SIGlist),names(SIGlist)))
    for ( i in 1:(length(SIGlist)-1) ) {
        for ( j in (i+1):length(SIGlist) ) {
            n1 <- length(SIGlist[[i]])
            n2 <- length(SIGlist[[j]])
            ov <- length(intersect(SIGlist[[i]],SIGlist[[j]]))
            OVLP[i,j] <- test.overlap(novlp=ov,nset1=n1,nset2=n2,ntotal=ntotal)
            NOVLP[i,j] <- paste(ov,n1,n2,sep="|")
            JOVLP[i,j] <- ov/length(union(SIGlist[[i]],SIGlist[[j]]))
        }
    }
    if ( adjust ) {
        OVLP[upper.tri(OVLP)] <- p.adjust(OVLP[upper.tri(OVLP)],method="BH")
    }
    if ( !is.null(outfile) ) {
        my.write.table(OVLP,names="overlap",file=outfile)
        cat("\n",file=outfile,append=TRUE)
        my.write.table(NOVLP,names="Noverlap",file=outfile,append=TRUE)
        cat("\n",file=outfile,append=TRUE)
        my.write.table(JOVLP,names="Jaccard",file=outfile,append=TRUE)
    }
    return( list(signatures=SIGlist,overlap=OVLP,noverlap=NOVLP,jaccard=JOVLP) )
}
#######################################################################
## support functions
#######################################################################

## function: TEST OVERLAP
##
## hyper-geometric test of two signatures overlap

test.overlap <- function(novlp,nset1,nset2,ntotal)
{
  if ( novlp<1 ) return (1.0)
  phyper(q=novlp-1,k=nset1,m=nset2,n=ntotal-nset2,lower.tail=FALSE)
}
## function: LIST 2 TABLE
##
list2table <- function( obj, fill=NA )
{
  if ( !is.list(obj) )
    stop( "input object must be a list" )

  mx <- matrix( fill, nrow=max(sapply(obj,length)), ncol=length(obj), dimnames=list(NULL,names(obj)) )
  for ( i in 1:length(obj) ) {
    if ( length(obj[[i]])<1 ) next
    mx[1:length(obj[[i]]),i] <- obj[[i]]
  }
  mx
}
