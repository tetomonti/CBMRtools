runGSEA_biDir <- function(GSigs, # List of "GeneSig" objects
                          pcut = 0.2, # FDR cutoff
                          minGenes = 10, # Min Genes in Gene Set
                          maxGenes = 500, # Max Gene in Gene Set
                          nperm = 1000, # Number of permutations
                          outDir # Output directory for writing gene sets
                          ){
  require(CBMRtools)
  
  GSEAlist <- list()
  for(i in 1:length(GSigs)){
    
    GSigI <- GSigs[[i]]
    SigI <- GSigI@Signature
    MDI <- GSigI@MetaData
    NameI <- as.character( MDI["Name",] )
    
    rnkScoreI <- qnorm(abs(SigI$P.Value)/2, 0, 1,lower.tail = FALSE) * sign(SigI$Score)
    names(rnkScoreI) <- SigI$ID
    
    J <- (1:length(GSigs))
    
    listJ <-list()
    indexJ <- 0
    for(j in J){
      GSigJ <- GSigs[[j]]
      SigJ <- GSigJ@Signature
      MDJ <- GSigJ@MetaData
      NameJ <- as.character( MDJ["Name",] )
      
      # Remove observations with missings IDs
      SigJ <- SigJ[!is.na(SigJ$ID),]
      SigJ <- SigJ[!is.na(SigJ$P.Value),]
      SigJ$P.adj <- p.adjust(SigJ$P.Value, method = "BH")
      SigJ <- SigJ[order(SigJ$P.Value),]
      
      # Separate positive and negatively enriched genes
      upJ <- SigJ[SigJ$Score > 0,]
      downJ <- SigJ[SigJ$Score < 0,]
      
      # Find number of positively enriched genes
      upCount <- sum(upJ$P.adj < pcut)
      if(upCount >= minGenes & upCount <= maxGenes)  upJ <- upJ[upJ$P.adj < pcut,]
      if(upCount < minGenes) upJ <- upJ[1:minGenes,]
      if(upCount > maxGenes) upJ <- upJ[1:maxGenes,]
      
      downCount <- sum(downJ$P.adj < pcut)
      if(downCount >= minGenes & downCount <= maxGenes)  downJ <- downJ[downJ$P.adj < pcut,]
      if(downCount < minGenes) downJ <- downJ[1:minGenes,]
      if(downCount > maxGenes) downJ <- downJ[1:maxGenes,]
      
      #Write gene sets to file
      upOut <- GeneSig(
        SigType = "gene_list",
        MetaData = MDJ,
        Signature = data.frame(upJ$ID, row.names = NULL)
      )
      write.GeneSig(upOut, paste(NameJ, "up", sep="_"), outDir)
      
      downOut <- GeneSig(
        SigType = "gene_list",
        MetaData = MDJ,
        Signature = data.frame(downJ$ID, row.names = NULL)
      )
      write.GeneSig(downOut, paste(NameJ, "down", sep="_"), outDir)
      
      
      
      upJ <- upJ$ID
      upJ <- unique( upJ[!is.na(upJ)] )
      
      downJ <- downJ$ID
      downJ <- unique( downJ[!is.na(downJ)] )
      
      indexJ <- indexJ + 1
      listJ[[indexJ]] <- upJ
      names(listJ)[indexJ] <- paste(NameJ, "up", sep="_")
      
      indexJ <- indexJ + 1
      listJ[[indexJ]] <- downJ
      names(listJ)[indexJ] <- paste(NameJ, "down", sep="_")
    }
    
    GSEAlistI <- cbmGSEA( rnkScore = rnkScoreI,
                          gSet = listJ,
                          nperm = nperm,
                          minGset=1,
                          seed = 1
    )
    GSEAlistI$dir <- 1
    GSEAlistI$dir[ grepl("_down", rownames(GSEAlistI) )] <- -1
    matched <- sign(GSEAlistI$score) == GSEAlistI$dir
    
    #Make p2 a one-sided
    GSEAlistI$p2 <- GSEAlistI$p2/2
    GSEAlistI$p2 <- (!matched)*1 - GSEAlistI$p2
    GSEAlistI$p2 <- GSEAlistI$p2 * sign(GSEAlistI$p2)
    GSEAlistI$fdr <- p.adjust(GSEAlistI$p2, method = "BH")
    
    #Make asymptotic.p one-directional
    GSEAlistI$asymptotic.p <- GSEAlistI$asymptotic.p/2
    GSEAlistI$asymptotic.p <- (!matched)*1 - GSEAlistI$asymptotic.p
    GSEAlistI$asymptotic.p <- GSEAlistI$asymptotic.p * sign(GSEAlistI$asymptotic.p)
    GSEAlistI$asymptotic.fdr <- p.adjust(GSEAlistI$asymptotic.p, method = "BH")
    
    GSEAlist[[i]] <- GSEAlistI
    names(GSEAlist)[i] <- NameI
  }
  GSEAfram <- do.call(rbind, GSEAlist)
  
  SigName <- unlist( strsplit( rownames(GSEAfram), "[.]" ) )[seq(1, (nrow(GSEAfram)*2), by = 2)]
  
  SetName <- unlist( strsplit( rownames(GSEAfram), "[.]" ) )[seq(2, (nrow(GSEAfram)*2), by = 2)]
  SetName <- factor(SetName, levels = paste( rep( unique(SigName), each = 2 ), c("up", "down"), sep = "_") )
  
  SigName <-factor(SigName, levels = rev( unique(SigName) ) )
  
  rownames(GSEAfram) <- NULL
  
  GSEAout <- cbind(Signature = SigName, Set = SetName, GSEAfram)
  
  
  return(GSEAout)
}