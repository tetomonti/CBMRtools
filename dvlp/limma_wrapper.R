#######################################################################
## function: LIMMA ACROSS
##
## Generalization of run_limma to deal with confounders
#######################################################################

limma.across <- function
(
    eset,
    clines,
    cond,
    control,
    sort.by="adj.P.Val",
    decreasing=FALSE
)
{
    eset <- eset[,(pData(eset)$cellline %in% clines) &
                  (pData(eset)$treatment %in% c(control,cond))]
    cline <- factor(pData(eset)$cellline)
    treat <- factor(pData(eset)$treatment,levels=c(control,cond))
    design <- model.matrix(~ cline + treat)
    fit <- lmFit(eset, design)
    fit2 <- eBayes(fit)
    fit2.table <- topTable(fit2, coef=paste("treat",cond,sep=""), adjust="BH",
                           number=length(fit2), sort.by="none")
    fit2.table$limma.fold.change <- 2^fit2.table$logFC
    fit2.accessory <- run_limma_accessory(eset, "treatment", cond, control)

    fdat.colnames <- colnames(fData(eset))
    fit2.table <- join(x=fit2.table, y=fit2.accessory, by=fdat.colnames, type="left", match="first")
    colnames(fit2.table)[which(colnames(fit2.table)=="logFC")]<-"limma.logFC"

    ordered.names <- c(fdat.colnames,
                      "high.class","t", "P.Value","adj.P.Val","actual.fold.change",
                      "limma.fold.change","actual.mean0","actual.mean1",
                      "actual.sd0","actual.sd1","n0","n1",
                      "limma.logFC", "actual.logFC")
    fit2.table <- fit2.table[, ordered.names]
    if ( !is.null(sort.by) ) {
        if( !(sort.by %in% colnames(fit2.table)) ) {
            stop(paste(sort.by, " is not in colnames of output table\n", sep =""))
        } else {
            fit2.table <- fit2.table[order(fit2.table[,sort.by], decreasing=decreasing),]
        }
    }
    return(fit2.table)
}
