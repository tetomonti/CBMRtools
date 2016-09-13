library(limma)
library(plyr)

#' wrapper for running limma with additional columns attached
#' @import Biobase limma plyr
#' @param eset expression set object with gene expression
#' @param treatment column name in pData(eset) specifying the conditions for differential expression
#' @param cond value in treatment specifying the positive condition
#' @param control value in treatment specifying the control condition
#' @param verbose whether to display additional columns in report table
#' @param sort.by column to sort by
#' @param decreasing direction for sorting
#' @param cutoff threshold value for keeping results
#' @examples
#'	data(eSet.brca.100)
#'	eset<-eSet.brca.100
#'	eset<-eset[, eset$ER_status %in% c("Positive", "Negative")]
#'	res<-run_limma(eset = eset, treatment = "ER_status", cond = "Positive", 
#'		control = "Negative", verbose = TRUE, sort.by = "adj.P.Val", decreasing = FALSE, cutoff = 0.05)
#' @export

#assumes values in eset are log2 transformed
run_limma <- function(eset, treatment, cond, control, verbose = TRUE,
                      sort.by = "adj.P.Val", decreasing = FALSE, cutoff = NA)
{
    if(!(verbose %in% c(TRUE, FALSE)))
        stop("verbose must in be c(TRUE, FALSE)")
    if(ncol(fData(eset))<1)
        stop("empty fData")
    
    cat("running limma\n")
    test_name<-paste(control, "_vs_", cond, sep = "")
    treatmentvec<-pData(eset)[, treatment]
    eset<-eset[, treatmentvec %in% c(cond, control)]
    treatmentvec<-pData(eset)[, treatment]
    design <- model.matrix(~ 0 + factor(treatmentvec))
    colnames(design) <- levels( factor(treatmentvec))
    fit <- lmFit(eset, design)
    command_str<-paste("makeContrasts(",
                       "(", cond , "-", control, ")", 
                       ",levels = design)", sep = "")
    contrast.matrix<-eval(parse(text =command_str)) 
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
                                        #return full table without sorting
    fit2.table<-topTable(fit2, coef=1, adjust="BH", number =length(fit2) ,
                         sort.by = "none")
    
    if (verbose){
        ## appending additional columns
        cat("running accessory function\n")
        fit2.accessary<-run_limma_accessory(eset, treatment, cond, control)
        
        cat("joining table\n")
        fdat.colnames<-colnames(fData(eset))
        
        fit2.table<-join(x=fit2.table, y=fit2.accessary, by = fdat.colnames, type = "left",
                         match = "first")
        
        fit2.table$limma.fold.change<-2^fit2.table$logFC
        colnames(fit2.table)[which(colnames(fit2.table) == "logFC")]<-"limma.logFC"
#	 ordered.names<-c(fdat.colnames,"limma.logFC", "limma.fold.change","t", "P.Value",
#	                  "adj.P.Val","actual.mean0","actual.mean1","high.class",
#			  "actual.fold.change","actual.sd0","actual.sd1","n0","n1","actual.logFC")
        ordered.names<-c(fdat.colnames,
                         "high.class","t", "P.Value","adj.P.Val","actual.fold.change",
                         "limma.fold.change","actual.mean0","actual.mean1",
                         "actual.sd0","actual.sd1","n0","n1",
                         "limma.logFC", "actual.logFC")
        fit2.table<-fit2.table[, ordered.names]	
    }
    if(!is.na(sort.by)){
        if(!(sort.by %in% colnames(fit2.table))){
            stop(paste(sort.by, " is not in colnames of output table\n", sep =""))
        } else {
            fit2.table<-fit2.table[order(fit2.table[, sort.by], decreasing = decreasing),]
        }
    }
    if(!is.na(cutoff)){
        if(decreasing){
            fit2.table<-fit2.table[fit2.table[,sort.by] > cutoff,]
        }
        else {
            fit2.table<-fit2.table[fit2.table[,sort.by] < cutoff,]
        }
    }
    return(fit2.table)
}

##assumes values in eset are log2 transformed
#' @export
run_limma_accessory<-function(eset, treatment, cond, control){

	test_name<-paste(control, "_vs_", cond, sep = "")
	treatmentvec<-pData(eset)[, treatment]
	eset<-eset[, treatmentvec %in% c(cond, control)]
	treatmentvec<-pData(eset)[, treatment]

	exprs(eset)<-2^exprs(eset)
	eset.sub.cond<-eset[, treatmentvec %in% c(cond)]
	eset.sub.control<-eset[, treatmentvec %in% c(control)]

	exprs.cond<-exprs(eset.sub.cond)
	exprs.control<-exprs(eset.sub.control)

	mean.control<-apply(exprs.control,1, mean)
	sd.control<-apply(exprs.control,1, sd)
	n.control<-rep(ncol(exprs.control), nrow(exprs.cond))
	
	mean.cond<-apply(exprs.cond,1, mean)
	sd.cond<-apply(exprs.cond,1, sd)
	n.cond<-rep(ncol(exprs.cond), nrow(exprs.cond))

	fc<-mean.cond/mean.control
	logfc<-log2(fc)

	high.class<-mean.cond >= mean.control
	high.class[high.class == TRUE]<- cond
	high.class[high.class == FALSE]<- control

	res<-data.frame(actual.mean0 = mean.control, actual.mean1 = mean.cond,
		 actual.sd0 = sd.control,  actual.sd1 = sd.cond,
		 n0 = n.control, n1 = n.cond, 
		 actual.fold.change = fc, 
		 actual.logFC = logfc, 
		 high.class = high.class)

	res<-data.frame(res, fData(eset))
	return(res)
}

#test function, do not export
test_limma_wrapper<-function(){
	#source("limma_wrapper.R")
	#library(CBMRtools)

	data(eSet.brca.100)
	eset<-eSet.brca.100
	eset<-eset[, eset$ER_status %in% c("Positive", "Negative")]
	res<-run_limma(eset = eset, treatment = "ER_status", cond = "Positive", 
		control = "Negative", verbose = TRUE, sort.by = "adj.P.Val", decreasing = FALSE, cutoff = 0.05)

}	
