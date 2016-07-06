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
    diffanalList,      # list of diffanal results (each a data.frame, see limma_wrapper)
    key="hgnc_symbol", # column id of gene identifiers in the table
    maxFDR=1.0,        # max FDR for genes to include in signature
    minFC=1.0,         # min fold-change for genes to include in signature
    ntotal=nrow(diffanalList[[1]]),
    fdrKey="adj.P.Val",
    fcKey="limma.fold.change",
    tKey="t",
    adjust=TRUE,
    outfile=NULL
)
{
    ## First, extract signature (i.e., gene lists) from diffanal results
    SIGlist <- diffanal2signatures(diffanalList,key=key,maxFDR=maxFDR,minFC=minFC,fdrKey=fdrKey,fcKey=fcKey,tKey=tKey)

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
    diffanalList,             # list of diffanal results (possibly, data frames)
    key="hgnc_symbol",        # column header where to find gene symbols
    maxFDR=1.0,               # max FDR to be included in the signatures
    minFC=1.0,                # min fold-change to be included
    fdrKey="adj.P.Val",       # column IDs (defaults are run_limma columns)
    fcKey="limma.fold.change",# ..
    tKey="t"                  # ..
)
{
    SIGlist <- NULL
    if ( !(key %in% colnames(diffanalList[[1]])) ) {
        stop( "column key",key,"not found")
    }
    for ( i in 1:length(diffanalList) )
    {
        DIFi <- diffanalList[[i]]
        SIGi <- list(up=DIFi[as.numeric(DIFi[,fdrKey])<=maxFDR &
                             as.numeric(DIFi[,fcKey])>=minFC &
                             as.numeric(DIFi[,tKey])>0,key],
                     dn=DIFi[as.numeric(DIFi[,fdrKey])<=maxFDR &
                             1/as.numeric(DIFi[,fcKey])>=minFC &
                             as.numeric(DIFi[,tKey])<0,key])
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
    OVLP <- matrix("",length(SIGlist),length(SIGlist),dimnames=list(names(SIGlist),names(SIGlist)))
    NOVLP <- matrix("",length(SIGlist),length(SIGlist),dimnames=list(names(SIGlist),names(SIGlist)))
    for ( i in 1:(length(SIGlist)-1) ) {
        for ( j in (i+1):length(SIGlist) ) {
            n1 <- length(SIGlist[[i]])
            n2 <- length(SIGlist[[j]])
            ov <- length(intersect(SIGlist[[i]],SIGlist[[j]]))
            OVLP[i,j] <- test.overlap(novlp=ov,nset1=n1,nset2=n2,ntotal=ntotal)
            NOVLP[i,j] <- paste(ov,n1,n2,sep="|")
        }
    }
    if ( adjust ) {
        OVLP[upper.tri(OVLP)] <- p.adjust(OVLP[upper.tri(OVLP)],method="BH")
    }
    if ( !is.null(outfile) ) {
        my.write.table(OVLP,names="overlap",file=outfile)
        cat("\n",file=outfile,append=TRUE)
        my.write.table(NOVLP,names="Noverlap",file=outfile,append=TRUE)
    }
    return( list(signatures=SIGlist,overlap=OVLP,noverlap=NOVLP) )
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
