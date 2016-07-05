## Function to annotate probesets with gene symbols, locations, etc.
## probesets (filters) can be ENSEMBLE IDs, Affy IDs, ...

probesetAnnotation <- function
(
    eset,
    use.gs=FALSE,
    filters="ensembl_gene_id",
    verbose=TRUE,
    na.rm=TRUE
)
{
  ## FEATURES' annotation
  ##
  if ( grepl("_at",featureNames(eset)[1]) && filters=="ensembl_gene_id" ) {
      featureNames(eset) <- gsub('_at','',featureNames(eset))
      VERBOSE(verbose,"postfix '_at' removed from featureNames\n")
  }
  mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

  ## retrieve probes' annotation
  VERBOSE(verbose,"retrieving feature ID mappings ..")
  GS <- getBM(filters=filters,
              attributes= c("ensembl_gene_id","entrezgene","hgnc_symbol","chromosome_name",
                  "band","strand","start_position","end_position","description"),
              values=featureNames(eset),
              mart= mart)
  VERBOSE(verbose, " done.\n")
  
  ## make sure that all ensemble replicates correspond to same gene symbol
  RP <- table(GS[,"ensembl_gene_id"]); RP <- RP[RP>1]
  if ( !all(CHK <- sapply(names(RP),function(Z) length(unique(GS[grep(Z,GS[,1]),3])))) )
    warning( "not all ensemble replicates correspond to same gene symbol" )
  
  ## pick 1st instance of each ensemble occurrence (duplicates are due
  ## to multiple entrez IDs matching the same ensemble ID)
  ##
  GS <- GS[match(unique(GS[,"ensembl_gene_id"]),GS[,"ensembl_gene_id"]),]
  rownames(GS) <- GS[,'ensembl_gene_id']

  ## remove NA's
  if ( na.rm ) {
      GS <- GS[!is.na(GS[,"hgnc_symbol"]),]
  }
  ## reformat eset data and add probes' annotation
  match.idx <- match.nona(rownames(GS),featureNames(eset))
  eset1 <- eset[match.idx,]

  ## if option selected, use gene symbols as featureNames
  if ( use.gs ) {
    if ( length(GS[,"hgnc_symbol"]!=length(unique(GS[,"hgnc_symbol"]))) )
      stop( "replicate gene symbols" )
    featureNames(eset1) <- rownames(GS) <- GS[,"hgnc_symbol"]
  }
  else {
    fData(eset1) <- GS
  }
  eset1
}
