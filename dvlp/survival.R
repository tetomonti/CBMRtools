##################################################################
#                      Coxph-based gene ranking
##################################################################

coxph.pvalue <- function( x, digits=3, logtest=FALSE )
{
  # this is copied from the procedure print.coxph. It simply implements 
  # the computation of the wald statistic z and corresponding p-value.
  #
  # INPUT:  x is the object returned from a call to coxph()
  # OUTPUT: the wald statistic and the corresponding p-value

  coef <- x$coef
  se <- sqrt(diag(x$var))
  if(is.null(coef) | is.null(se))
    stop("Input is not valid")
  z <- coef/se
  
  # p-value of logtest
  #
  if ( logtest )
  {
    lt <- -2 * (x$loglik[1] - x$loglik[2])
    df <- if (is.null(x$df)) sum(!is.na(coef)) else round(sum(x$df),2)
    p <- signif(1 - pchisq(lt, df),digits)
  }
  # alternatively, p-value of coefficient
  #
  else {
    p <- signif(1 - pchisq((coef/se)^2, 1), digits )
  }
  return ( c(z,p) )
}
gene.coxph <- function( x, y, bare=FALSE, verbose=FALSE, singular.ok=FALSE )
{
  # INPUT:
  #   - x      (m x n) matrix for m genes and n experiments
  #   - y      Surv object
  #   - bare   if true, return p-value only
  #
  # OUPUT:
  #   (m x 2) matrix reporting z-statistic and p-value as
  #   computed by coxph for each gene
  #
  if ( !is.Surv(y) )
    stop( "y must be of type Surv" )

  z <- t(apply( x,1,function(w){ coxph.pvalue( coxph(y~w,singular.ok=singular.ok) ) } ))

  colnames(z) <- c("z","p")

  if ( bare )
    return( z[,2] )
  else
    return(z)
}
perm.coxph <- function( x, y, nperm=100, seed=NULL, verbose=FALSE, online=TRUE, two.sided=FALSE )
{
  if ( two.sided )
  {
    VERBOSE( verbose, "Two-sided test.\n" )
    score <- function(x,cls,verbose) { gene.coxph( x, y=cls, bare=FALSE, verbose=verbose )[,1] }
    x.perm <-
        perm.2side( x, y, nperm=nperm, score=score, seed=seed, online=online, rnd=3, verbose=verbose )$pval
  }
  else
  {
    VERBOSE( verbose, "One-sided test.\n" )
    score <- function(x,cls,verbose) { gene.coxph( x, y=cls, bare=FALSE, verbose=verbose )[,2] }
    x.perm <-
        perm.1side( x, y, nperm=nperm, score=score, seed=seed, online=online, rnd=3, verbose=verbose )$pval
  }
  x.perm <- x.perm[,c("score","p2","fdr","maxT")]
  pval <- gene.coxph(x,y,bare=TRUE)
  x.perm <- cbind(x.perm,asymp.p=pval,asymp.fdr=pval2fdr(pval))
  x.perm
}
##################################################################
#                      KM-based gene ranking
##################################################################

survdiff.pvalue <- function( sobj )
{
  if ( !is.list(sobj) | names(sobj)[5]!="chisq" ) stop( "Survdiff object expected" )
  return( 1-pchisq(sobj$chisq,length(sobj$n)-1) )
}
survdiff.matrix <- function(sobj)
{
  out <- cbind(N=sobj$n,
               Observed=round(sobj$obs,2),
               Expected=round(sobj$exp,2),
               '(O-E)^2/E'=round((sobj$obs-sobj$exp)^2/sobj$exp,3),
               '(O-E)^2/B'=round(diag((sobj$obs-sobj$exp)^2/sobj$var),3))
  out <- cbind(out,chisq=round(sobj$chis,2))
  out <- cbind(out,pval=signif(survdiff.pvalue(sobj),2))
  out
}
km.survfit <- function( x, S, nmin=1, probs=NULL, digits=3, debug=FALSE )
{
  # INPUT:
  #   - x        n-sized vector of predictive values
  #   - S        n-sized Surv object
  #   - nmin     min acceptable number of observations in a split
  #   - probs    prequantize x by calling quantile(x,probs=probs)
  #
  # OUPUT:
  #   3-tuple with:
  #   - thresh   the selected threshold in the range of x
  #   - p.value  the p.value returned by survfit regressed on 'x>thresh'
  #   - geq      the number of x's for which x>thresh is true
  #
  if ( !is.Surv(S) )
    stop( "S must be a Surv object" ) 
  if ( nmin>(length(x)-1)/2 )
    stop( paste("nmin too large", nmin) )
  if ( nmin<1 )
    stop( paste("nmin too small:", nmin) )
  if ( !is.null(probs) & nmin>1 )
    warning( "nmin>1 with quantization deprecated" )
  n <- length(unique(x))
  if ( n<nmin*2 )
    stop( paste("not enough x values for given nmin:",n,"(", nmin*2,"needed )") )
  if ( !is.null(probs) & n<length(probs) )
    stop( paste("not enough x values for given quantiles:",n,"(", length(probs),"needed )") )
  
  idx <- !is.na(x)
  x <- x[idx]
  S <- S[idx,]

  # define the set of candidate thresholds
  #
  x.srt <- sort(unique(x))
  x.thresh <- if ( is.null(probs) ) 
    ( x.srt[-1] + x.srt[-length(x.srt)] ) / 2
  else
    quantile(x.srt,probs)
  if (is.null(probs))
    x.thresh <- x.thresh[nmin:(length(x.thresh)-nmin+1)]

  # x.split is a matrix with as many columns as thresholds, with each
  # column indicating which x-values are less than or equal to (FALSE) and
  # greater than (TRUE) the corresponding threshold
  #
  x.split  <- t(apply(matrix(x,length(x),length(x.thresh)),1,'>',x.thresh))

  x.p <- apply( x.split,2,function(z){(1-pchisq(survdiff(S~z)$chisq,1))} )
  best <- which.min( x.p )
  if ( !is.null(digits) ) x.p <- signif(x.p,digits)
  
  if ( !debug )
    return( c( x.thresh[best],      # best trheshold
               x.p[best],           # p-value associated to best threshold
               sum(x.split[,best]), # no. of obs greater than threshold
               sum(S[x.split[,best],1])>sum(S[!x.split[,best],1]))
                                    # does greater than threshold predict good prognosis
           )
  else
    return( cbind(x.thresh, x.p, apply(x.split,2,sum)) )
}
gene.km.survfit <- function( x, surv, nmin=1, srt=FALSE, probs=NULL, bare=FALSE, debug=FALSE, ... )
{
  x.km <- t(apply( x,1,km.survfit,surv,nmin=nmin,probs=probs,debug=debug))

  if ( debug ) return(x.km)

  # ELSE ..
  #
  colnames(x.km) <- c( "thresh", "p.value", "gt", "prognosis" )
  
  if ( srt ) {
    x.ord <- order(x.km[,2])
    x.km <- cbind( x.km[x.ord,],x.ord )
    colnames(x.km)[ncol(x.km)] <- "order"
  }
  if (bare)
    return (x.km[,2])
  else
    return (x.km)
}
plot.km.survfit <- function( gct, surv, km, r, col=c("blue","red"), lty=c(1,2) )
{
  i.pv <- match( "p.value", colnames(km) ); if (is.na(i.pv)) stop("p.value not found" )
  i.th <- match( "thresh", colnames(km) );  if (is.na(i.th)) stop("thresh not found" )
  i.gt <- match( "gt", colnames(km) );      if (is.na(i.gt)) stop("gt not found" )
  i.od <- match( "order", colnames(km) );
  if ( is.na(i.od) )
  {
    warning( "order not found, genes sorted automatically\n" )
    ord  <- order(km[,i.pv])
    km   <- cbind( km[ord,], ord )
    colnames(km)[ncol(km)] <- "order"
    i.od <- ncol(km)
  }
  gname <- rownames(km)[r]
  n   <- nrow(surv)
  nx  <- max(surv[,1])
  off <- .05*nx
  
  x.t <- as.vector(sign(gct[km[r,i.od],]>km[r,i.th]))
  plot( survfit(surv~x.t), col=col, lty=lty )  
  legend( off, .01, col=col, lty=lty, yjust=0,
          legend=paste( gname, c("<=",">"), round(km[r,i.th],2), sep="") )
  legend( nx-off,.01,col=col,lty=lty,xjust=1, yjust=0,
          legend=paste("n=",c(nrow(surv)-km[r,i.gt],km[r,i.gt])))
  text( x=nx-off, y=1, labels=paste("p-value =",km[r,i.pv]), adj=c(1,0) )
  #cat("offset:", off, "\n" )
  #cat("p-value =", km[r,i.pv], "\n" )
}
perm.survfit <- function( x, y, nperm, nmin=2, probs=NULL, seed=NULL, online=TRUE, verbose=FALSE )
{
  # define the score function
  #
  score <- function(x,y,nmin,probs,verbose)
  {
    gene.km.survfit( x, y, nmin=nmin, probs=probs, bare=TRUE, verbose=FALSE )
  }
  # call the generic permutation routine
  #
  x.perm <- perm.1side( x, y, nperm=nperm, score=score, seed=seed,
                        online=online, verbose=verbose, nmin=nmin, probs=probs )
  (x.perm)
}
debug.perm.survfit <- function()
{
  source("read.res.R")
  library(survival)
  gct.x <- read.gct("~/wi/scripts/tests/t_Survfit.gct" )
  gct.y <- read.table("~/wi/scripts/tests/t_Survfit.txt",header=TRUE)[,c(4,3)]
  gct.s <- Surv(gct.y[,1],gct.y[,2])
  x.km <- gene.km.survfit( gct.x, gct.s, probs=seq(.1,.9,.1),nmin=1 )
  x.prm <- perm.survfit( gct.x, gct.s, nperm=100, probs=seq(.2,.8,.1), nmin=1, verbose=TRUE )
}
# source("~/dvlp/R/survival.R")
