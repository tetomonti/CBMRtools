#' Run DESeq2 on eset object
#' @import DESeq2
#' @param eset expression set object with raw count data
#' @param class_id column name in pData(eset) specifying the conditions for differential expression
#' @param treatment value in class_id specifying the treatment condition
#' @param control value in class_id specifying the control condition
#' @param verbose add additional information and perform filtering
#' @param sortByPValue whether to sort results by nominal p-value, if verbose = TRUE
#' @param FDRcutoff FDR threshold for filtering results, if verbose = TRUE
#' @export

#Assumes values in eset are raw count values
run_deseq<-function(eset, class_id, control, treatment, verbose = TRUE, sortByPValue = TRUE, FDRcutoff = NULL){
  require(DESeq2)
  
  if(!(verbose %in% c(TRUE, FALSE)))
    stop("verbose must in be c(TRUE, FALSE)")
  
  control_inds<-which(pData(eset)[, class_id] == control)
  treatment_inds<-which(pData(eset)[, class_id] == treatment)
  eset.control<-eset[, control_inds]
  eset.treatment<-eset[, treatment_inds]
  eset.compare<-eset[, c(control_inds, treatment_inds)]
  
  #make deseq2 compliant dataset
  colData<-data.frame(condition=as.character(pData(eset.compare)[, class_id]))
  dds<-DESeqDataSetFromMatrix(exprs(eset.compare), colData, formula( ~ condition))
  
  #set reference to control, otherwise default is alphabetical order
  dds$condition <- factor(dds$condition, levels=c(control,treatment))
  
  cat("running DESeq2\n")
  #run deseq2
  #3 steps:
  #1.) estimate size factors
  #2.) estimate dispersion
  #3.) negative binomial GLM fitting and wald test
  dds_res<-DESeq(dds)
  res <- results(dds_res)
  
  if(verbose == TRUE){
    cat("Adding additional information and filtering results\n")
    res <- format_results_deseq(deseq_out = res, eset = eset, class_id = class_id, control = control, treatment = treatment, sortByPValue = sortByPValue,  FDRcutoff = FDRcutoff)
  }
  return(res)
}

## Function for formatting output of DESeq2
#' @export
format_results_deseq<-function(deseq_out, eset, class_id, control, treatment, sortByPValue, FDRcutoff){
  require(DESeq2)
  
  desOut <- as.data.frame(deseq_out)
  
  eset <- eset[rownames(desOut),]
  control_inds<-which(pData(eset)[, class_id] == control)
  treatment_inds<-which(pData(eset)[, class_id] == treatment)
  eset.compare<-eset[, c(control_inds, treatment_inds)]
  
  #make deseq2 compliant dataset
  colData<-data.frame(condition=as.character(pData(eset.compare)[, class_id]))
  dds<-DESeqDataSetFromMatrix(exprs(eset.compare), colData, formula( ~ condition))
  
  #Normalize the expression data
  dds <- estimateSizeFactors(dds)
  NormExp <- counts(dds, normalized = TRUE)
  
  #Calculate the within group means and standard deviations
  NormCont <- NormExp [, pData(eset.compare)[, class_id] == control]
  NormTreat <- NormExp [, pData(eset.compare)[, class_id] == treatment]
  
  desOut$contMean <- rowMeans(NormCont)
  desOut$contSD <- apply(NormCont, 1, sd)
  desOut$treatMean <- rowMeans(NormTreat)
  desOut$treatSD <- apply(NormTreat, 1, sd)
  
  # Transform log2 fold change to regular fold change
  desOut$FoldChange <- 2^desOut$log2FoldChange
  
  # Get high class
  desOut$HighClass <- NA
  desOut$HighClass[desOut$log2FoldChange >= 0 ] <- treatment
  desOut$HighClass[desOut$log2FoldChange < 0 ] <- control
  
  # Change column order and rename columns to group names
  desOut<-desOut[, c("baseMean", 
                     "FoldChange", 
                     "log2FoldChange", 
                     "lfcSE", 
                     "treatMean", 
                     "contMean", 
                     "treatSD", 
                     "contSD", 
                     "stat",
                     "HighClass",
                     "pvalue", 
                     "padj")]
  
  colnames(desOut)[11] <- "PValue"
  colnames(desOut)[12] <- "FDR"
  
  colnames(desOut)<-sub("cont", paste(control, "_", sep =""), colnames(desOut))
  colnames(desOut)<-sub("treat", paste(treatment, "_", sep =""), colnames(desOut))
  
  # Sort by nominal p-value
  if(sortByPValue == TRUE) desOut <- desOut[order(desOut$PValue),]
  
  # Filter out non-significant results (IF FDRcutoff != NULL)
  if(!is.null(FDRcutoff)) {desOut <- desOut[desOut$FDR <= FDRcutoff,]
    desOut <- desOut[!is.na(desOut$FDR),]
  }
  return(desOut)}