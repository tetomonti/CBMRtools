#####################################################################################
## BEGIN documentation support (what follows are keyworded entries
## from which documentation pages will be extracted automatically)

#' Carry out set enrichment test based on hyper-geometric distribution
#'
#' @param drawn One or more sets of 'drawn' items (e.g., genes). Basically, a list of signatures.
#' @param categories list of gene sets (e.g., MSigDB c2)
#' @param ntotal background population, i.e., the total no. of items from which items are supposed to have been drawn
#' @param min.drawsize min no. of drawn items that must be among categories' items
#' @param mht correct for multiple hypothesis testing across multiple 'draws'
#'
#' @return a data.frame with rows indexed by the signature(s) tested
#'
#' @examples
#'
#' # load objects hyperSig (a list of signatures) and hyperGsets (a GeneSet object)
#' # and run hyper-enrichment test
#'
#' data(hyper) # contains objects hyperSig and hyperGsets
#' hyperE <- hyperEnrichment(drawn=hyperSig,categories=getGeneSet(hyperGsets),ntotal=10000)
#' head(hyperE)
#' 
#' @export 

## END documentation support
#####################################################################################

## function: HYPER ENRICHMENT
##
hyperEnrichment <- function
(
   drawn,          # one or more sets of 'drawn' items (e.g., genes). Basically, a list of signatures.
   categories,     # gene sets (list of gene sets)
   ntotal=length(unique(unlist(categories))),
                   # background population, i.e., the total no. of items from which
                   # ..items are supposed to have been drawn
   min.drawsize=4, # min no. of drawn items that must be among categories' items
   mht=TRUE,       # correct for MHT across multiple 'draws'
   verbose=TRUE
)
{
   ## checks on inputs
   ##
   if (!is(categories, "list") ) {
     stop( "categories expected to be a list of gene sets" )
   }
   gene.names<-unique(unlist(categories))
   if ( is.list(drawn) && is.null(names(drawn)) ) {
     stop( "drawn must have non-null names when a list" )
   }
   if ( ntotal<length(unique(unlist(categories)))) {
     warning( "background population's size less than unique categories' items: ", ntotal,"<",length(gene.names))
   }
   ##
   ## end checks
   
   cnames <-
      c("pval","fdr","set annotated","set size","category annotated","total annotated","category","hits")
   
   ## handling of multiple 'draws'
   ##
   if ( is.list(drawn) ) 
   {
     ncat <- length(categories)
     ENRICH <- NULL
     
     VERBOSE(verbose,"Testing",length(drawn),"drawsets on",ncat,"categories and",
             length(gene.names),"total items ..\n")
     
     percent <- 0.1
     base <- 0
     ntst <- 0
     for ( i in 1:length(drawn) )
     {
       VERBOSE(verbose,"*** Testing", names(drawn)[i], ".. " )
       dset <- drawn[[i]]
       tmp <- hyperEnrichment(dset,categories,ntotal=ntotal,verbose=verbose)
       if (is.null(tmp)) {
         VERBOSE(verbose,"not enough items drawn\n")
         next
       }
       ntst <- ntst+1
       rng <- (base+1):(base+ncat)
       if (any(!is.na(enrich[rng,]))) stop( "something wrong")
       
       enrich[rng,] <- cbind(set=rep(names(drawn)[i],ncat),tmp)

       ENRICH <- rbind(ENRICH,cbind(set=names(drawn)[i],tmp,stringsAsFactors=FALSE))
                       
       base <- base+ncat
       if (FALSE && i>=round(length(drawn)*percent)) {
         VERBOSE(verbose, round(100*percent),"% ",sep="")
         percent <- percent+0.1
       }
       VERBOSE(verbose," (min fdr: ", signif(min(as.numeric(tmp[,"fdr"])),2),")\n",sep="")
     }
     VERBOSE(verbose,"done.\n")
     colnames(enrich) <- c("set",cnames)
      
     enrich <- enrich[1:base,,drop=FALSE]
     if (mht) {
       VERBOSE(verbose,"MHT-correction across multiple draws ..")
       enrich[,"fdr"] <- pval2fdr(as.numeric(enrich[,"pval"]))
       VERBOSE(verbose,"done.\n")
     }    
     VERBOSE(verbose,
             "Categories tested: ",rjust(length(categories),4),"\n",
             "Candidate sets:    ",rjust(length(drawn),4),"\n",
             "Sets tested:       ",rjust(ntst,4),"\n",
             "Items tested:      ",rjust(sum(sapply(drawn,length)),4)," (min,med,max: ",
             paste(quantile(sapply(drawn,length),probs=c(0,.5,1)),collapse=","),")\n",
             "N(FDR<=0.25):      ",rjust(sum(enrich[,"fdr"]<=.25),4),"\n",
             "N(FDR<=0.05):      ",rjust(sum(enrich[,"fdr"]<=.05),4),"\n",
             "N(FDR<=0.01):      ",rjust(sum(enrich[,"fdr"]<=.01),4),"\n",
             sep="")
     return(enrich)
   }
   ## handling of a single draw
   ##
   m.idx <- drawn[drawn %in% gene.names] 
   
   if ( length(m.idx)<min.drawsize ) {
     VERBOSE(verbose,"insufficient annotated genes in the drawn set: ",
             paste(gene.names[m.idx],collapse=","),"\n")
     return(NULL)
   }
   VERBOSE(verbose,length(m.idx),"/",length(drawn), " annotated genes found",sep="")
   
   nhits <-sapply(categories, function(x,y) length(intersect(x,y)), m.idx)
   ndrawn <- length(drawn) # length(m.idx)
   ncats <- sapply(categories,length)
   nleft <- ntotal-ncats
   
   ## compute P[X>=nhits]
   enrich <- phyper(q=nhits-1,m=ncats,n=nleft,k=ndrawn,lower.tail=F)
   ##enrich <- cbind(pval=enrich,
   enrich <- data.frame(pval=enrich,
                        fdr=pval2fdr(enrich),
                        nhits=nhits,
                        ndrawn=ndrawn,
                        ncats=ncats,
                        ntot=ntotal,
                        category=names(categories),
                        hits=sapply(categories,function(x,y)paste(intersect(x,y),collapse=','),m.idx),
                        stringsAsFactors=FALSE)  

   ord <- order(as.numeric(enrich[,"pval"]))
   enrich <- enrich[ord,,drop=FALSE]
   enrich[,"pval"] <- signif(as.numeric(enrich[,"pval"]),2)
   enrich[,"fdr"] <- signif(as.numeric(enrich[,"fdr"]),2)
   
   colnames(enrich) <- cnames
   rownames(enrich) <- names(categories)[ord]
   
   return(enrich)
}
## GENERATION OF THE DATA FOR THE EXAMPLE
##
if ( FALSE )
{
  ## data generation
  ##
  CBMMLAB <- Sys.getenv('CBMMLAB')
  CBMGIT <- Sys.getenv('CBMGIT')
  CBMRtools <- Sys.getenv('CBMRtools')
  source(paste(CBMRtools,'R/misc.R',sep='/'))
  source(paste(CBMRtools,'R/GeneSet.R',sep='/'))
  source(paste(CBMRtools,'R/hyperEnrichment.R',sep='/'))
  source(paste(CBMRtools,'R/broad.file.formats.R',sep='/'))

  GS <- GeneSet(paste(CBMMLAB,'/annot/c2.cp.v4.0.symbols.gmt',sep=''))

  SIG <- read.tab.delim('~/Research/Projects/oralcancer/taz_yap_dpagt1/results/SIGtab.xls')
  hyperSig <- table2list(SIG,fill="")

  HYP <- hyperEnrichment(drawn=hyperSig,categories=getGeneSet(GS),ntotal=10000)

  hyperGsets <- GS
  hyperGsets@source.file <- GS@name
  hyperGsets@name <- basename(GS@name)
  save(hyperSig,hyperGsets,file=paste(CBMRtools,'data/hyper.rda',sep='/'))

  ## testing the script outside the package
  ##
  load(file=paste(CBMRtools,'data/hyper.rda',sep='/'))
  HYP <- hyperEnrichment(drawn=hyperSig,categories=getGeneSet(hyperGsets),ntotal=10000)
}
