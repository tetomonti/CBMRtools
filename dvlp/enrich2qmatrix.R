#######################################################################
## function: CBM GSEA 2 QMATRIX
##
## Take a list of outputs from multiple calls to cbmGSEA and generate
## a summary geneset-by-signature matrix (and optionally a heatmap) of
## the genesets' q-values
##
#######################################################################
cbmGSEA2qmatrix <- function
(
    cgsea,           # list of data.frames, the up and down genesets, as output by cbmGSEA
    fdr=c(.05,.01),  # FDR thresholds (must be in decreasing order)
    method=c("union","intersect"),
                     # union or intersection of genesets across signatures
    do.sort=TRUE,    # sort matrices by HC
    do.heat=FALSE,   # display heatmap
    rm.zero=TRUE,    # remove genesets/rows w/ no hits
    pvalID="fdr",    # which significance measure to use (either "p2" or "fdr")
    globalMHT=FALSE, # correct for MHT *across* signatures (default is within)
    na.col="gray",   # color for missing values in the heatmap
    verbose=TRUE,    # verbose output
    outfile=NULL,    # save qmatrix output to workbook file
    sheetName="sheet1", # add a new worksheet to workbook
    annotation_col=NULL, # data frame of GSEA result annotations
    annotation_row=NULL, # data frame of pathway annotations
    annotation_colors=NULL, # List of names colors for annotations
    aggMethod=c("ward.D2","ward.D","single","complete","average","mcquitty","median","centroid"),
    ...              
)
{
    ## each element of the variable cgsea corresponds to a distinct cbmGSEA
    ## run, and it contains a data.frame object

    ## input checks
    if ( length(fdr)>1 && any(diff(fdr)>0) ) stop( "fdr must be in decreasing order" )
    if ( globalMHT ) pvalID <- "p2"
    if ( !(pvalID %in% colnames(cgsea[[1]])) ) stop( "unrecognized pvalID: ", pvalID)
    method <- match.arg(method)
    combineFun <- match.fun(method)

    VERBOSE(verbose,"Generating qmatrix based on",method,"..")
    
    ## extract all the geneset IDs in common among all cbmGSEA runs
    gID <- Reduce(combineFun,lapply(cgsea,rownames))

    ## extract the FDRs of all genesets for each cbmGSEA run 
    mx <- matrix(NA,length(gID),length(cgsea),dimnames=list(gID,names(cgsea)))
    for ( i in 1:ncol(mx) ) {
        z <- cgsea[[i]]
        z[z[,"score"]<0,pvalID] <- -z[z[,"score"]<0,pvalID]
        if ( method=="intersect" ) 
            mx[,i] <- z[match.nona(rownames(mx),rownames(z)),pvalID]
        else if ( method=="union" )
            mx[match.nona(rownames(z),rownames(mx)),i] <- z[,pvalID]
        else
            stop("unrecognized method:",method)
    }    
    VERBOSE(verbose, " done, [", paste(dim(mx),collapse=","),"] matrix generated.\n",sep="")

    ## if global multiple hypothesis correction (MHT), take the uncorrected
    ## p-values across signatures and carry out a global FDR correction
    ##
    if ( globalMHT ) {
        VERBOSE(verbose, "carrying out global MHT ..")
        absMX <- abs(mx)
        absMX[,] <- p.adjust(as.vector(absMX),method="BH")
        mx[is.na(mx)] <- 999 # remove NA's (necessary for 'mx<0' test)
        absMX[mx<0] <- -absMX[mx<0]
        mx <- absMX
        VERBOSE(verbose, " done, min(fdr):", min(abs(mx),na.rm=TRUE), "max(fdr):", max(mx,na.rm=TRUE),"\n")
    }
    q2h <- qmatrix2heatmap(mx=mx,fdr=fdr,do.sort=do.sort,do.heat=do.heat,rm.zero=rm.zero,na.col=na.col,annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, aggMethod = aggMethod, ...)
    
    if ( !is.null(outfile) )
    {
        out <- qmatrix2workbook(q2h$mx,q2h$mx01,col=c("blue","white","red"), sheetName = sheetName, outfile = outfile, annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
        saveWorkbook(out,file=outfile,overwrite=TRUE)
        VERBOSE(verbose,"Workbook saved to '",outfile,"'\n",sep="")
    }
    return( q2h )
}
#######################################################################
## function: GSEA 2 QMATRIX
##
## Take a list of outputs of multiple calls to runGSEA or
## runGSEApreranked and generate a summary geneset-by-signature matrix
## (and optionally a heatmap) of the genesets' q-values
##
#######################################################################
gsea2qmatrix <- function
(
    gsea,               # list of lists (each with 2 data.frames, the up and down genesets, as output by gsea)
    fdr=c(.05,.01),     # FDR thresholds (must be in decreasing order)
    method=c("union","intersect"),
                        # union or intersection of genesets across signatures
    do.sort=TRUE,       # sort matrices by HC
    do.heat=FALSE,      # display heatmap
    rm.zero=TRUE,       # remove genesets/rows w/ no hits
    pvalID="FDR q-val", # which significance measure to use (see GSEA output for choices)
    na.col="gray",      # color for missing values in the heatmap
    verbose=TRUE,       # verbose output
    outfile=NULL,       # save qmatrix output to workbook file
    zero=1.0e-10,       # min p-value 
    ...                 # extra arguments to qmatrix2heatmap
)
{
    ## each element of the list gsea corresponds to a distinct GSEA
    ## run, and it contains a two-element list corresponding to the
    ## files 'gsea_report_for_na_{pos,neg}_*.xls' output by GSEA
    ##
    if ( length(fdr)>1 && any(diff(fdr)>0) ) stop( "fdr must be in decreasing order" )
    if ( !(pvalID %in% colnames(gsea[[1]][[1]])) ) stop( "unrecognized pvalID: ", pvalID)
    if ( !all(sapply(gsea,function(Z1) sapply(Z1,function(Z2) is.character(Z2[,'NAME'])))) )
        stop( "column 'NAME' expected to be of type character")
        
    method <- match.arg(method)
    combineFun <- match.fun(method)

    VERBOSE(verbose,"Generating qmatrix based on",method,"..")

    ## extract all the geneset IDs in common among all GSEA runs
    gID <- Reduce(combineFun,lapply(gsea,function(Z) c(Z[[1]][,'NAME'],Z[[2]][,'NAME'])))

    ## extract the FDRs of all genesets for each GSEA run (i.e., each list item)
    mx <- matrix(NA,length(gID),length(gsea),dimnames=list(gID,names(gsea)))
    for ( i in 1:ncol(mx) ) {
        z <- gsea[[i]]
        up <- z[[1]][,c('NAME',pvalID),drop=FALSE] # positive sign p-values
        up[,pvalID] <- pmax(zero,up[,pvalID])      # ..
        dn <- z[[2]][,c('NAME',pvalID),drop=FALSE] # negative sign p-values
        dn[,pvalID] <- -pmax(zero,dn[,pvalID])     # ..
        tmp <- rbind(up,dn)
        
        if ( method=="intersect" ) 
            mx[,i] <- tmp[match.nona(rownames(mx),tmp[,'NAME']),pvalID]
        else if ( method=="union" )
            mx[match.nona(tmp[,'NAME'],rownames(mx)),i] <- tmp[,pvalID]
        else
            stop("unrecognized method:",method)
    }
    q2h <- qmatrix2heatmap(mx=mx, fdr=fdr, do.sort=do.sort, do.heat=do.heat, rm.zero=rm.zero, na.col=na.col, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, ...)
    
    if ( !is.null(outfile) )
    {
        out <- qmatrix2workbook(q2h$mx,fdr=fdr,col=c("blue","white","red"), sheetName = sheetName, outfile = outfile)
        saveWorkbook(out,file=outfile,overwrite=TRUE)
        VERBOSE(verbose,"Workbook saved to '",outfile,"'\n",sep="")
    }
    return( q2h )
}
#######################################################################
## function: QMATRIX 2 HEATMAP
##
## (Function called by both gsea2qmatrix and cbmGSEA2qmatrix)
## Take a geneset-by-signature matrix of q-values, and generate:
## 1) a 'discretized' matrix based on the input threshold(s)
## 2) (optionally) a color-coded heatmap
#######################################################################
qmatrix2heatmap <- function
(
    mx,             # geneset-by-signature matrix of q-values
    fdr=c(.05,.01), # FDR thresholds (must be in decreasing order)
    do.sort=TRUE,   # sort matrices by HC
    do.heat=FALSE,  # display heatmap
    rm.zero=TRUE,   # remove genesets/rows w/ no hits
    verbose=TRUE,   # extra arguments to my.heatmap
                    # hclust methods
    aggMethod=c("ward.D2","ward.D","single","complete","average","mcquitty","median","centroid"),
    na.col="gray",  # color for missing values in the heatmap
    annotation_col=NULL, # data frame of GSEA result annotations
    annotation_row=NULL, # data frame of pathway annotations
    annotation_colors=NULL, # List of names colors for annotations
    ...
)
{
    require(pheatmap)
    aggMethod <- match.arg(aggMethod)
    
    levs <- 0:length(fdr)
    zero <- -0.000001
    mx01 <- suppressWarnings(matrix(cut(as.vector(mx),
                                        breaks=c(-1,-fdr,zero,rev(fdr),1),
                                        labels=as.numeric(c(-levs,rev(levs))),include.lowest=TRUE),
                                    nrow=nrow(mx),ncol=ncol(mx)))
    mx01 <- apply(mx01,2,as.numeric)
    dimnames(mx01) <- dimnames(mx)

    if ( rm.zero ) {
        rm.idx <- apply(mx01!=0,1,any,na.rm=TRUE)
        if ( sum(rm.idx)==0 ) {
            VERBOSE(verbose,"no significant genesets found, not removing any")
        }
        else {
            mx01 <- mx01[rm.idx,,drop=FALSE]
            mx <- mx[rm.idx,,drop=FALSE]
            VERBOSE(verbose, "Removed",sum(!rm.idx),"non-significant genesets.\n")
        }
    }
    heat <- NULL
    ## sort rows and columns by HC
    if ( do.sort || do.heat )
    {
        ## Remove gene sets that cause crash - no overlap with another gene set
        dist.row <-dist(mx01,method="euclidean")
        rmC <- 0
        while(sum(is.na(dist.row)) > 0){
            rmC <- rmC + 1
            remGS <- names( sort( colSums( is.na( as.matrix(dist.row) ) ), decreasing = TRUE ) )[1]
            mx01 <- mx01[rownames(mx01)!=remGS,]
            mx <- mx[rownames(mx)!=remGS,]
            dist.row  <- dist(mx01,method="euclidean")
        }
        if( rmC > 0 ){
            VERBOSE(verbose, "Removed", rmC ,"genesets due to no overlap with another gene set(s).\n")
        }
        hc.row <- hcopt(dist.row,method=aggMethod)
        hc.col <- hcopt(dist(t(mx01),method="euclidean"),method=aggMethod)
        
        if ( do.heat ) {
            
            mxMin <- min(mx01, na.rm = T)
            mxMax <- max(mx01, na.rm = T)
            mxU <- seq(mxMin, mxMax, by = 1)
 
            mxAbsMax <- max(abs(mxU))
            ncolors <- mxAbsMax*2+1
            
            COL <- col.gradient(c("blue","white","red"),length=ncolors)
            names(COL) <- as.character(seq(-mxAbsMax, mxAbsMax, by = 1))
            COL <- COL[names(COL) %in% mxU]
  
            
            # Create heatmap
          heat <- pheatmap(mx01,
                     color = COL,
                     cluster_rows = hc.row,
                     cluster_col = hc.col,
                     annotation_col = annotation_col,
                     annotation_row = annotation_row,
                     annotation_colors = annotation_colors,
                     legend = FALSE
                     )

        }
        if ( do.sort ) {
            mx <- mx[hc.row$order,hc.col$order]
            mx01 <- mx01[hc.row$order,hc.col$order]
        }
    }
    return( list(mx=mx,mx01=mx01,heat=heat) )
}
#######################################################################
## function: HYPER 2 QMATRIX
##
## Take the output of hyperEnrichment (called on multiple signatures)
## and generate a summary geneset-by-signature matrix (and optionally
## a heatmap) of the genesets' q-values
#######################################################################
hyper2qmatrix <- function
(
    hyper,          # output of hyperEnrichment
    fdr=c(.05,.01), # FDR thresholds (must be in decreasing order)
    do.sort=TRUE,   # sort matrices by HC
    do.heat=FALSE,  # display heatmap
    rm.zero=TRUE,   # remove genesets/rows w/ no hits
    globalMHT = FALSE, # correct for MHT *across* signatures (default is within)
    aggMethod=c("ward.D2","ward.D","single","complete","average","mcquitty","median","centroid"), # hclust methods
    outfile=NULL,   # save qmatrix output to workbook file
    sheetName = "sheet1",
    verbose = TRUE,
    annotation_col=NULL, # data frame of GSEA result annotations
    annotation_row=NULL, # data frame of pathway annotations
    annotation_colors=NULL, # List of names colors for annotations
    ...             
)
{
    if ( length(fdr)>1 && any(diff(fdr)>0) )
        stop( "fdr must be in decreasing order" )
    if ( length(unique(hyper[,"set"]))==1 )
        stop( "hyper must contain results for at least two signatures")
    aggMethod <- match.arg(aggMethod)
    
    SIG <- unique(hyper[,"set"])
    gsets <- hyper[hyper[,"set"]==SIG[1],"category"]
    mx <- sapply(SIG,function(z) {
        tmp <- hyper[hyper[,"set"]==z,c("fdr","category")]
        as.numeric(tmp[match.nona(gsets,tmp[,"category"]),"fdr"])
    });
    rownames(mx) <- gsets
    
    ## if global multiple hypothesis correction (MHT), take the uncorrected
    ## p-values across signatures and carry out a global FDR correction
    ##
    if ( globalMHT ) {
      VERBOSE(verbose, "carrying out global MHT ..")
      absMX <- abs(mx)
      absMX[,] <- p.adjust(as.vector(absMX),method="BH")
      mx[is.na(mx)] <- 999 # remove NA's (necessary for 'mx<0' test)
      absMX[mx<0] <- -absMX[mx<0]
      mx <- absMX
      VERBOSE(verbose, " done, min(fdr):", min(abs(mx),na.rm=TRUE), "max(fdr):", max(mx,na.rm=TRUE),"\n")
    }
    q2h <- qmatrix2heatmap(mx=mx,fdr=fdr,do.sort=do.sort,do.heat=do.heat,rm.zero=rm.zero,na.col=na.col,annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, aggMethod = aggMethod, ...)
    
    if ( !is.null(outfile) )
    {
      out <- qmatrix2workbook(q2h$mx,q2h$mx01,col=c("blue","white","red"), sheetName = sheetName, outfile = outfile, annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
      saveWorkbook(out,file=outfile,overwrite=TRUE)
      VERBOSE(verbose,"Workbook saved to '",outfile,"'\n",sep="")
    }
    return( q2h )
}
#######################################################################
## function: QMATRIX 2 WORKBOOK
##
## Save the qmatrix into a Excel workboox conditionally-formatted to
## match the heatmap.
#######################################################################
qmatrix2workbook <- function
(
    qmatrix,        # matrix of signed q-values
    imatrix,        # matrix of signed integers
    col=c("blue","white","red"),
    fcol="black",    # font color
    sheetName,
    outfile,
    annotation_row = NULL,
    annotation_col = NULL, 
    annotation_colors = NULL,
    ...
)
{
  
    require(openxlsx)
    # Get start row
    startRow=1
    if(!is.null(annotation_col)){
      startRow = 2+ncol(annotation_col)
      cNames <- colnames(annotation_col)
      annotation_col<- as.data.frame(annotation_col[colnames(qmatrix),])
      rownames(annotation_col) <- colnames(qmatrix); colnames(annotation_col) <- cNames
    }
    
    # Get start column
    startCol=1
    if(!is.null(annotation_row)){
      startCol = 2+ncol(annotation_row)
      cNames <- colnames(annotation_row)
      annotation_row<- as.data.frame(annotation_row[rownames(qmatrix),])
      rownames(annotation_row) <- rownames(qmatrix); colnames(annotation_row) <- cNames
    }
    
    
    # Write data
    if(file.exists(outfile)){
      wb <- loadWorkbook(outfile)
      shNames <- names(wb)
      if(sheetName %in% shNames){
        removeWorksheet(wb, sheetName)
      }
      
    } else {
      wb <- createWorkbook()
    }
    addWorksheet(wb,sheetName)
    
    # Write row names and column names
    writeData(wb,sheetName,rownames(qmatrix),startCol=1,startRow=startRow+1, colNames = F, rowNames = F)
    writeData(wb,sheetName,t(colnames(qmatrix)),startCol=startCol+1,startRow=1, colNames = F, rowNames = F)
    
    # Write data
    writeData(wb,sheetName,qmatrix,startCol=startCol + 1,startRow=startRow + 1, colNames = F, rowNames = F)
    
    # Add colors to fill in heatmap
    mxMin <- min(imatrix, na.rm = T)
    mxMax <- max(imatrix, na.rm = T)
    mxU <- seq(mxMin, mxMax, by = 1)
    
    mxAbsMax <- max(abs(mxU))
    ncolors <- mxAbsMax*2+1
    
    COL <- col.gradient(col,length=ncolors)
    names(COL) <- as.character(seq(-mxAbsMax, mxAbsMax, by = 1))
    COL <- COL[names(COL) %in% mxU]

    ## Fill in colors based on significance
    for ( i in 1:length(COL) )
    {
        colInd <- which(imatrix == as.numeric(names(COL)[i]), arr.ind=TRUE) 
        colInd[,1] <- colInd[,1] + startRow
        colInd[,2] <- colInd[,2] + startCol
        if ( length(colInd)>0 ){
            addStyle(wb, sheetName, style=createStyle(fgFill=COL[i]), rows=colInd[,1], cols=colInd[,2])
        }
    }
    
    # Add values to row annotation
    if(!is.null(annotation_row)){
      writeData(wb,sheetName,annotation_row,startCol=startCol,startRow=startRow,rowNames=FALSE)
      if(!is.null(annotation_colors)){
        for(i in 1:ncol(annotation_row)){
          anCol <- annotation_colors[[colnames(annotation_row)[i]]]
          for(j in 1:length(anCol)){
            rowInd <- which(annotation_row[,i] == names(anCol)[j]) + startRow
            if ( length(rowInd)>0 ){
              addStyle(wb, sheetName, style=createStyle(fgFill=anCol[j]), rows=rowInd, cols=(startCol + (i-1)))
            }
          }
        }
      }
    }
    
    # Add values to column annotation
    if(!is.null(annotation_col)){
      writeData(wb,sheetName, t(annotation_col),startCol=startCol+1, startRow=i+1, colNames = F, rowNames = F)
      writeData(wb,sheetName, colnames(annotation_col), startCol=startCol+nrow(annotation_col)+1, startRow=i+1, colNames = F, rowNames = F)
      if(!is.null(annotation_colors)){
        for(i in 1:ncol(annotation_col)){
          anCol <- annotation_colors[[colnames(annotation_col)[i]]]
          for(j in 1:length(anCol)){
            colInd <- which(annotation_col[,i] == names(anCol)[j]) + startCol
            if ( length(colInd)>0 ){
              addStyle(wb, sheetName, style=createStyle(fgFill=anCol[j]), rows=i+1, cols=colInd)
            }
          }
        }
      }
    }
    
    # Freeze header cells
    freezePane(wb, sheetName, firstActiveRow = startRow+1, firstActiveCol = startCol+1)
    
    # Save workbook
    saveWorkbook(wb, outfile, overwrite = T)
    return( wb )
}
## this was an alternative version
##
qmatrix2workbook2 <- function
(
    qmatrix,        # matrix of signed q-values
    fdr=c(.05,.01), # FDR thresholds (must be in decreasing order)
    col=c("blue","white","red"),
    fcol="black"    # font color
)
{
    bi <- any(qmatrix<0)        # bi- or uni-directional q-values?
    if ( !bi && length(col)>2 ) # two colors expected when qmatrix is uni-directional
        warning( "more than two colors input with unidirectional qmatrix?" )
    
    levs <- 0:length(fdr); nlevs <- length(levs)
    ncolors <- if ( bi ) length(fdr)*2+1 else length(fdr)+1
    COL <- col.gradient(col,length=ncolors)
    
    posCol <- if (bi) COL[length(COL):(nlevs+1)] else rev(COL[-1])
    negCol <- if (bi) COL[1:(nlevs-1)]
    
    shName <- "sheet1"
    wb <- createWorkbook()
    addWorksheet(wb,shName)
    writeData(wb,shName,qmatrix,startCol=1,startRow=1,rowNames=TRUE)

    ## handling of the positive q-values
    ##
    fdr <- c(0.0,rev(fdr),1.0)
    if ( length(fdr)!=length(posCol)+2 ) stop( "length(fdr)!=length(posCol)+2" )
    for ( i in 1:length(posCol) )
    {
        colInd <- which(qmatrix>fdr[i] & qmatrix<=fdr[i+1], arr.ind=TRUE)+1  
        if ( length(colInd)>0 ){
            addStyle(wb, shName, style=createStyle(fgFill=posCol[i]), rows=colInd[,1], cols=colInd[,2])
        }
    }
    ## if bi-directional, handling of the negative q-values
    ##
    if ( bi )
    {
        if ( length(fdr)!=length(negCol)+2 ) stop( "length(fdr)!=length(negCol)+2" )
        for ( i in 1:length(fdr) )
        {
            colInd <- which(qmatrix < -fdr[i], arr.ind=TRUE)+1  
            if ( length(colInd)>0 ){
                addStyle(wb, shName, style=createStyle(fgFill=negCol[i]), rows=colInd[,1], cols=colInd[,2])
            }
        }
    }
    wb
}

##################################################################
## code to create and test the data objects used in the examples #
##################################################################
if ( FALSE )
{
    require(CBMRtools)
    CBMRtools <- Sys.getenv("CBMRtools")
    source(file.path(CBMRtools,"dvlp/heatmap.R"))
    source(file.path(CBMRtools,"dvlp/enrich2qmatrix.R"))
    source(file.path(CBMRtools,"dvlp/diffanalOverlap.R"))
    ## path in my desktop
    PATH <- "~/Research/Projects/AhR/AHR_CYP1B1_CB799_microarrays2016"
    GSPATH <- "~/Research/CancerGenomeAnalysisData/annot"
    
    ## paths on SCC:
    #PATH <- "/restricted/projectnb/montilab-p/projects/AhR/AHR_CYP1B1_CB799_microarrays2016"
    #GSPATH <- "/restricted/projectnb/montilab-p/CBMrepositoryData/annot/"

    ## upload MSigDB's hallmark geneset compendium (stored locally)
    gSet <- new("GeneSet",file.path(GSPATH,"h.all.v5.1.symbols.gmt"))
    
    ## loading pre-computed list of diffanal results
    DIF <- readRDS(file.path(PATH,"results/rds/DIF2.RDS"))
    names(DIF)
    ##  [1] "MDA.AhR" "MDA.CYP" "SUM.AhR" "SUM.CYP" "AhR"     "CYP"     "CB113"  
    ##  [8] "BaP"     "PYO"     "CH"     
    
    ## loading pre-computed list of cbmGSEA results (this takes a long time to compute)
    cGSEA <- readRDS(file.path(PATH,"results/rds/cbmGSEA.RDS"))
    names(cGSEA)
    ## [1] "h.all"  "c2.cp"  "c3.all" "c2.cgp"
    
    ## extract the hallmarks results (hGSEA is a list of data.frames, one per signature)
    hGSEA <- cGSEA[['h.all']]
    names(hGSEA)
    ## [1] "MDA.AhR" "MDA.CYP" "SUM.AhR" "SUM.CYP" "AhR"     "CYP"     "CB113"  
    ## [8] "BaP"     "PYO"     "CH"     
    
    # Run without annotations
    png("ahr.hallmarks.png")
    OUT <- cbmGSEA2qmatrix(hGSEA,fdr=c(0.10,0.05,0.01),do.heat=TRUE,globalMHT=TRUE)
    dev.off()
    
    # Run with row and column annotations
    
    # Row annotation
    annot_row <- data.frame(source = rep("Hallmarks", length(names(gSet@geneset))), row.names = names(gSet@geneset))
    
    # Column annotation (This is a meaningless example)
    annot_col <- data.frame(pert = sub("MDA[.]|SUM[.]", "",names(hGSEA)), row.names = names(hGSEA))
    annot_col$first_letter <- substr(names(hGSEA), 1, 1)
    
    # Add colors
    annot_colors <- list(source = "green",
                         pert = brewer.pal(length(unique(annot_col$pert)), "Set1"),
                         first_letter = brewer.pal(length(unique(annot_col$first_letter)), "Set2"))
    names(annot_colors$source) <- "Hallmarks"
    names(annot_colors$pert) <- unique(annot_col$pert)
    names(annot_colors$first_letter) <- unique(annot_col$first_letter)
    
    png("ahr.hallmarks_annotated.png")
    OUTa <- cbmGSEA2qmatrix(hGSEA,fdr=c(0.10,0.05,0.01),do.heat=TRUE,globalMHT=TRUE, annotation_col = annot_col, annotation_row = annot_row, annotation_colors = annot_colors)
    dev.off()
    
    # Add sheet w/out annotations
    wOUT1 <- qmatrix2workbook(OUT$mx,OUT$mx01,outfile = "wOUT1.xlsx", sheet = "TESTSHEET1")
    
    wOUT2 <- qmatrix2workbook(OUT$mx,OUT$mx01,outfile = "wOUT1.xlsx", sheet = "TESTSHEET2", annotation_row = annot_row, annotation_col = annot_col, annotation_colors = annot_colors)
    
    
    ## let's test with different genesets
    hGSEA1 <- hGSEA
    for ( i in 1:length(hGSEA1) ) # remove a (different) geneset from each data.frame
      hGSEA1[[i]] <- hGSEA[[i]][-i,]
    
    ## take the intersection
    OUT1 <- cbmGSEA2qmatrix(hGSEA1,fdr=c(0.10,0.05,0.01),do.heat=TRUE,globalMHT=TRUE,method="intersect")
    
    ## take the union (there should be rows with missing values, color-coded 'gray' by default)
    OUT2 <- cbmGSEA2qmatrix(hGSEA1,fdr=c(0.10,0.05,0.01),do.heat=TRUE,globalMHT=TRUE,method="union")
    
    ## TEST HYPER-ENRICHMENT VISUALIZATION
    
    SIG <- diffanal2signatures(DIF, maxFDR=0.05, minFC=1.5)
    hOut <- hyperEnrichment(SIG,categories=getGeneSet(gSet),ntotal=22981)
    tmp <- hyper2qmatrix(hOut,fdr=c(.05,.01),do.heat=TRUE)
    my.write.table(tmp$mx,names="GeneSet",file="HYPERtable.xls")
    
    ## TESTING qmatrix2workbook
    
    require(openxlsx)
    wb1 <- qmatrix2workbook(tmp$mx,tmp$mx01, outfile = "wOUT1.xlsx", sheetName = "TESTSHEET3")
    
    ## Test again but with annotations
    # Column annotation (This is a meaningless example)
    annot_col <- data.frame(dir = gsub("[[:alnum:]]*[.]", "",unique(hOut[,1])), row.names = unique(hOut[,1]))
    
    # Add colors
    annot_colors <- list(source = "green",
                         dir = c("red", "blue"))
    
    names(annot_colors$source) <- "Hallmarks"
    names(annot_colors$dir) <- c("up", "dn")
    
    # Create qmatrix, heatmap, and worksheet
    tmp1 <- hyper2qmatrix(hOut,fdr=c(.05,.01),do.heat=TRUE,outfile = "wOUT1.xlsx", sheetName = "TESTSHEET4", annotation_col = annot_col, annotation_row = annot_row, annotation_colors = annot_colors)
}
