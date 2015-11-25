########################################################################
## CV CLASSIFY
########################################################################

cvClassify <- function
(
 dat,              # M x N matrix (M rows/cases, N columns/predictors)
 cls,              # M-vector of class labels
 model,            # classifier type (must be a list w/ fields 'estimate' and 'predict')
 nfeats,           # number of features (can be a vector, to return learning curve)
 nfolds=nrow(dat), # CV folds (default is loocv)
 fs.score=c("t.score","snr","mean.diff"),
 fs.balanced=T,    # pick 1/2 features for one class and 1/2 features for the other (only for binary phenotype)
 verbose=T,
 do.debug=F,
 ...               # model's specs
)
{
  ## BEGIN checks on inputs
  ##
  if ( any((match(cls,my.levels(cls))-1)!=cls) )
    stop("class labels must be contiguous integers starting at 0")
  if ( nrow(dat)!=length(cls) )
    stop( "number of samples must be same as length of class vector" )
  if ( nfolds==0 )
    nfolds=nrow(dat)
  if ( nfolds>nrow(dat) )
    stop( "nfolds cannot be greater than number of sample: ", nfolds )
  fs.score <- match.arg(fs.score)
  ##
  ## END checks on inputs
  
  VERBOSE(verbose,
          "Cross validation w/ following specs:\n",
          "dat:     [", paste(dim(dat),collapse="x"), "] matrix.\n",
          "cls:     {", paste(tabulate(as.factor(cls)),collapse=","), "}\n",
          "nfolds: ", nfolds, "\n",
          "nfeats:  {", paste(nfeats,collapse=","), "}\n",
          "fs.score: ", fs.score, "\n")

  levs <- levels(cls); nlevs <- length(levs)
  snames <- c( "ERR", "RERR", "ROC", "LS", "NERR", paste("NERR",levs,sep=".") )
  
  # define CV folds
  #
  folds <- if (nfolds!=1) xval.select( nrow(dat), nfolds, cls=cls ) else rep(1,nrow(dat))
  nfolds <- sort(unique(folds))
  
  # where to save prediction results
  #
  out.predict <- matrix(NA,nrow(dat),length(levs),dimnames=list(rownames(dat),paste("P(C=",levs,")",sep="")))
  out.predict <- lapply( 1:length(nfeats), function(z) out.predict )
  names(out.predict) <- paste( nfeats,"feats",sep="." )
  FEATS <- matrix( 0, ncol(dat), length(nfeats), dimnames=list(colnames(dat),names(out.predict)))
  NFEATS <- matrix( 0, length(nfolds), length(nfeats) )
  colnames(NFEATS) <- names(out.predict)
  rownames(NFEATS) <- paste("fold",1:nrow(NFEATS),sep=".")
  rownames(FEATS) <- colnames(dat)
  
  MODELS <- NULL
  # iterate over folds
  #
  percent <- .1
  for ( fold in nfolds )
  {
    #VERBOSE( verbose, "Fold", fold, ": " )
    #f.idx <- folds==fold
    tst.idx <- folds==fold
    trn.idx <- if (length(nfolds)==1) tst.idx else !tst.idx

    feats <- feature.select(t(dat[trn.idx,,drop=F]),
                            cls=subset.cls(cls,trn.idx),
                            nfeat=nfeats,
                            score=fs.score,
                            balanced=fs.balanced)
    if (length(nfeats)==1) feats <- list(feats)
    models <- NULL
    
    # iterate over number of features
    #
    for ( i in 1:length(nfeats) )
    {
      #VERBOSE( verbose, nfeats[i], ".", sep="" )
      model.fit <- model$estimate(dat=dat[trn.idx,feats[[i]],drop=F],
                                  cls=subset.cls(cls,trn.idx), ...)
      #VERBOSE( verbose, "." )
      model.prd <- model$predict(dat=dat[tst.idx,feats[[i]],drop=F],
                                 model=model.fit )
      out.predict[[i]][tst.idx,] <- model.prd
      
      fidx <- feats[[i]]
      if (!is.null(model$predictors)) {
        fidx <- match(model$predictors(model.fit),colnames(dat))
        if (any(is.na(fidx))) stop( "something wrong: missing predictors in dataset" )
      }
      if (do.debug) {
        models <- c( models, list(model.fit) )
      }
      FEATS[fidx,i] <- FEATS[fidx,i]+1
      NFEATS[fold,i] <- length(fidx)
    }
    if (do.debug) {
      MODELS <- c( MODELS, list(models) )
    }
    #VERBOSE( verbose, "\n" )

    if ( fold>=round(length(nfolds)*percent) ) {
      VERBOSE(verbose," ",round(100*percent),"%",sep="")
      percent <- percent+.1
    }
  }
  VERBOSE(verbose,"\n")
  OUT.predict <- matrix( NA, length(nfeats), length(snames) )
  for ( i in 1:length(nfeats) )
  {
    #if (i==5)
    #  debug(prediction.summary)
    OUT.predict[i,] <- prediction.summary( out.predict[[i]], cls )
  }
  nfeats[nfeats==0] <- ncol(dat)
  OUT.predict <- cbind( nfeats, OUT.predict )
  rownames(OUT.predict) <- names(out.predict)
  colnames(OUT.predict) <- c("Nfeats",snames)

  FEATS <- cbind( FEATS, tot=apply(FEATS,1,sum) )
  FEATS <- FEATS[FEATS[,"tot"]>0,]
  return( list(summary=OUT.predict,details=out.predict,features=FEATS,nfeatures=NFEATS,models=MODELS) )
}
cls.summary <- function( x, cls, eps=.Machine$double.neg.eps )
{
  nlevs <- length(levels(cls))
  x[x<=eps] <- eps
  x <- t(apply( x, 1, function(z) {z[z>=1] <- 1-sum(z[z==eps]); z}))
  prd <- apply( x, 1, which.max )-1
  cnf <- my.ftable( cls, prd, x.lev=0:(nlevs-1), y.lev=0:(nlevs-1) )
  ROC <- NA
  return(c(ERR=1-sum(diag(cnf))/sum(cnf),
           RERR=1-mean(diag(cnf)/apply(cnf,1,sum)),
           ROC=if(nlevs==2) AUC(rocdemo.sca(cls,x[,2],dxrule.sca)) else NA,
           LS=sum(log(x[cbind(1:nrow(x),cls+1)])),
           NERR=sum(cnf)-sum(diag(cnf)),
           apply(cnf,1,sum)-diag(cnf)))
  
}
prediction.summary <- function( x, cls, eps=.Machine$double.neg.eps, thresh=0.5, detailed=F )
{
  if ( is.null(nlevs <- nlevels(cls)) )
    stop( "cls must have a non-null levels attribute" )
  
  x[x<=eps] <- eps
  x <- t(apply( x, 1, function(z) {z[z>=1] <- 1-sum(z[z==eps]); z}))
  prd <- {
    if ( nlevels(cls)==2 )
      as.numeric(x[,2]>thresh)
    else
      apply( x, 1, which.max )-1
  }
  cnf <- my.ftable( cls, prd, x.lev=0:(nlevs-1), y.lev=0:(nlevs-1) )
  ROC <- NA
  if (nlevs==2) {
    ROC <- if ( all(x[,2]>(1-1.0e-12)) ) 1.0 else AUC(rocdemo.sca(cls,x[,2],dxrule.sca))
  }  
  out <- c(ERR=1-sum(diag(cnf))/sum(cnf),
           RERR=1-mean(diag(cnf)/apply(cnf,1,sum)),
           ROC=ROC,
           LS=sum(log(x[cbind(1:nrow(x),cls+1)])),
           NERR=sum(cnf)-sum(diag(cnf)),
           apply(cnf,1,sum)-diag(cnf))

  return( if (detailed) list(out=out,confusion=cnf) else out)
}
# function: BEST NF
#
bestNF <- function(dat,
                   score1=c("ERR","RERR","ROC","LS"),
                   score2=c("LS","ERR","RERR","ROC","none"),
                   occam=T,
                   verbose=T)
{
  score1 <- match.arg( score1 )
  score2 <- match.arg( score2 )

  VERBOSE( verbose, "  1st selection based on ", score1, "\n", sep="" )
  VERBOSE( verbose, "  2nd selection based on ", score2, "\n", sep="" )
  
  dat[,c("LS","ROC",if(!occam) "Nfeats")] <- -dat[,c("LS","ROC",if(!occam) "Nfeats")]
  
  ord <- order(dat[,score1],
               if(score2!="none") dat[,score2] else dat[,score1],
               dat[,"Nfeats"],
               decreasing=F)
  
  dat[,c("LS","ROC",if(!occam) "Nfeats")] <- -dat[,c("LS","ROC",if(!occam) "Nfeats")]
  
  return( list(nf=dat[ord[1],"Nfeats"], score=dat[ord[1],score1]) )
}
cvselect.predict <- function
(
 trn.dat,
 trn.cls,
 tst.dat,
 tst.cls=NULL,
 nfolds=nrow(dat),
 nfeats,
 fs.score=c("t.score","snr","mean.diff"),
 #feats=NULL,
 model=c("naive","knn","lps","lrgP","hclust","lda","qda","svm"),
 fs1=c("ERR","RERR","ROC","LS"), # level-1 criterion for choosing number of features
 fs2=c("LS","ERR","RERR","ROC"), # level-2 criterion for choosing number of features
 fs.balanced=T,
 verbose=T,
 ...
 )
{
  fs.score <- match.arg(fs.score)
  model <- match.arg(model)
  fs1 <- match.arg(fs1)
  fs2 <- match.arg(fs2)

  model.fun <- switch(model,
                      naive=list(estimate=naive.estimate,
                                 predict=naive.predict),
                      knn=list(estimate=knn.estimate,
                               predict=knn.predict),
                      lps=list(estimate=if (length(nlevels(trn.cls))>2) lpsM.estimate else lpsT.estimate,
                               predict= if (length(nlevels(trn.cls))>2) lpsM.predict else lpsT.predict),
                      lrgP=list(estimate=largeP.estimate,
                                predict=largeP.predict),
                      hclust=list(estimate=hc.estimate,
                                  predict=hc.predict),
                      lda=list(estimate=xda.pca.estimate,
                               predict=xda.pca.predict),
                      svm=list(estimate=function(dat,cls,...) {
                                 svm(x=dat,y=as.factor(cls),probability=T,...)
                               },
                               predict=function(dat,model) {
                                 tmp <- attr(predict(object=model,newdata=dat,probability=T),"probabilities");
                                 tmp[,order(colnames(tmp)),drop=F]
                               }))

  VERBOSE(verbose, "*** Model selection ***\n")
  CV <- cv.classify(dat=trn.dat,cls=trn.cls,model=model.fun,nfeats=nfeats,nfolds=nfolds,
                    fs.score=fs.score,...)

  bestN <- bestNF( CV$summary, score1=fs1, score2=fs2 )
  VERBOSE(verbose, "  Best num. of features: ", bestN$nf, "\n",sep="")

  VERBOSE(verbose, "*** Prediction ***\n")
  FSI <- feature.select(t(trn.dat),trn.cls,score=fs.score,nfeat=bestN$nf,balanced=fs.balanced,verbose=verbose)

  CLSfit <- model.fun$estimate(dat=trn.dat[,FSI],cls=trn.cls,...)
  CLStst <- model.fun$predict(dat=tst.dat[,FSI],model=CLSfit)
  rownames(CLStst) <- rownames(tst.dat)

  tst.err <- NULL
  #bst.thr <- if ( nlevels(trn.cls)==2 )
  #  calibrate.classifier(model.fun$predict(dat=trn.dat[,FSI],cls=trn.cls,...)[,2],cls=trn.cls)
  
  if ( !is.null(tst.cls) ) {
    VERBOSE(verbose, "Test results:\n")
    tst.err <- prediction.summary(CLStst,cls=tst.cls,detailed=T)
    if (verbose) print(tst.err)
  }
  list(prediction=CLStst,model=CLSfit,summary=tst.err,selection=CV,bestF=colnames(trn.dat)[FSI])
}

calibrate.classifier <- function
(
 probs,  # probabilties output by the model
 cls     # true classification labels
 )
{
  if ( nlevels(cls)!=2 )
    stop( "calibration implemented for binary classification only" )
  if ( length(probs)!=length(cls) )
    stop( "probs and cls must be of same length" )

  ord <- order(probs)
  probs <- probs[ord]
  cls[1:length(cls)] <- cls[ord]
  
  THRESH <- probs[-length(probs)]+diff(probs)/2

  min.err <- Inf
  bst.thr <- NULL
  for ( thresh in THRESH )
  {
    smry <- prediction.summary(cbind(1-probs,probs),cls=cls,thresh=thresh)
    VERBOSE(verbose, round(thresh,4), ": ", round(smry[1],4), "\n", sep="")

    if (smry[1]<min.err) {
      bst.thr <- thresh
      min.err <- smry[1]
    }
  }
  bst.thr
}
