#' Run edgeR on eset object
#' @import edgeR
#' @param eset expression set object with raw count data
#' @param class_id column name in pData(eset) specifying the conditions for differential expression
#' @param treatment value in class_id specifying the treatment condition
#' @param control value in class_id specifying the control condition
#' @param verbose add additional information and perform significance filtering
#' @param sortByPValue whether to sort results by nominal p-value, if verbose = TRUE
#' @param FDRcutoff FDR threshold for filtering results, if verbose = TRUE
#' @export

#Assumes values in eset are raw count values
run_edgeR<-function(eset, class_id, control, treatment, verbose = TRUE, sortByPValue = TRUE, FDRcutoff = NULL){
  require(edgeR)
  
  if(!(verbose %in% c(TRUE, FALSE)))
    stop("verbose must in be c(TRUE, FALSE)")
  
  library(edgeR)
  
  control_inds<-which(pData(eset)[, class_id] == control)
  treatment_inds<-which(pData(eset)[, class_id] == treatment)
  
  eset.control<-eset[, control_inds]
  eset.treatment<-eset[, treatment_inds]
  eset.compare<-eset[, c(control_inds, treatment_inds)]
  condition<-factor(pData(eset.compare)[, class_id], levels = c(control, treatment))
  
  cat("running edgeR\n")
  es <- DGEList(counts=exprs(eset.compare), group = condition)
  es <- calcNormFactors(es)
  es <- estimateDisp(es)
  es <- exactTest(es)
  res<-topTags(es, n = nrow(eset.compare),  sort.by = "none")
  if(verbose == TRUE){
    cat("Adding additional information and filtering results\n")
    res <- format_results_edgeR(edgeR_out = res, eset = eset, class_id = class_id, control = control, treatment = treatment, sortByPValue = sortByPValue,  FDRcutoff = FDRcutoff)
  }
  
  return(res)
}


## Function for formatting output of edgeR
#' @export
format_results_edgeR<-function(edgeR_out, eset, class_id, control, treatment, sortByPValue, FDRcutoff){
  require(edgeR)
  
  edgerOut <- edgeR_out@.Data[[1]]
  
  eset <- eset[rownames(edgerOut),]
  control_inds<-which(pData(eset)[, class_id] == control)
  treatment_inds<-which(pData(eset)[, class_id] == treatment)
  eset.compare<-eset[, c(control_inds, treatment_inds)]
  condition<-factor(pData(eset.compare)[, class_id], levels = c(control, treatment))
  
  #make edgeR compliant dataset
  es <- DGEList(counts=exprs(eset.compare), group = condition)
  
  #Normalize the expression data
  es <- calcNormFactors(es)
  NormExp <- cpm(es, normalized.lib.sizes = TRUE, log=FALSE)
  
  #calculate base mean
  edgerOut$baseMean<-rowMeans(NormExp)
  
  #Calculate the within group means and standard deviations
  NormCont <- NormExp [, pData(eset.compare)[, class_id] == control]
  NormTreat <- NormExp [, pData(eset.compare)[, class_id] == treatment]
  
  edgerOut$contMean <- rowMeans(NormCont)
  edgerOut$contSD <- apply(NormCont, 1, sd)
  edgerOut$treatMean <- rowMeans(NormTreat)
  edgerOut$treatSD <- apply(NormTreat, 1, sd)
  
  # Transform log2 fold change to regular fold change
  edgerOut$FoldChange <- 2^edgerOut$logFC
  
  # Get high class
  edgerOut$HighClass <- NA
  edgerOut$HighClass[edgerOut$logFC >= 0 ] <- treatment
  edgerOut$HighClass[edgerOut$logFC < 0 ] <- control
  
  # Change column order and rename columns to group names
  edgerOut<-edgerOut[, c("baseMean", 
                         "logCPM", 
                         "FoldChange", 
                         "logFC", 
                         "treatMean", 
                         "contMean", 
                         "treatSD", 
                         "contSD",
                         "HighClass",
                         "PValue", 
                         "FDR")]
  
  colnames(edgerOut)[2] <- "baseLogCPM"
  colnames(edgerOut)[4] <- "log2FoldChange"
  
  colnames(edgerOut)<-sub("cont", paste(control, "_", sep =""), colnames(edgerOut))
  colnames(edgerOut)<-sub("treat", paste(treatment, "_", sep =""), colnames(edgerOut))
  
  # Sort by nominal p-value
  if(sortByPValue == TRUE) edgerOut <- edgerOut[order(edgerOut$PValue),]
  
  # Filter out non-significant results (IF FDRcutoff != NULL)
  if(!is.null(FDRcutoff)) {edgerOut <- edgerOut[edgerOut$FDR <= FDRcutoff,]
  edgerOut <- edgerOut[!is.na(edgerOut$FDR),]
  }
  
  return(edgerOut)}