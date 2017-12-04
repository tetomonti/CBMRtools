#######################################################################
## class: GeneSig
##
## Define a format for gene signatures
##
## Takes an input for three slots
##    1. SigType slot => either gene_rank or gene_list
##      a. gene_rank if gene signature is a ranked list
##      b. gene_list if gene signature is a gene set
##    2. MetaData slot => A data.frame with one column
##          and row names given the category of each row
##    3. Signature slot => A data.frame that must contain 4 columns
##      a. ID => Gene identifiers
##      b. P.Value => P.Values from DGE testing
##      c. Score => A bidirectional score to assign direction
##          to ranking
##      d. high.class => An identifier for the group with the
##          larger mean expression
##      e. Any additional column allowed
##
#######################################################################

GeneSig <- setClass(
    "GeneSig",
    slots = c(SigType = "character",
              MetaData = "data.frame",
              Signature = "data.frame"),
    validity = function(object)
    {
        errMessage = c()
        if( sum( c("gene_rank", "gene_list") %in% object@SigType) != 1) {
            errMessage <- c(errMessage, "SigType must be on of: 'gene_rank' or 'gene_list'.")
        }
        if( ncol(object@MetaData) != 1 ) {
            errMessage <- c( errMessage, "MetaData must be 1 column." )
        }
        geneCheck <- c("ID", "P.Value", "Score", "high.class") %in% colnames(object@Signature)
        if( !all(geneCheck) & object@SigType == "gene_rank" ) {
            errMessage <- c(errMessage,
                            paste( "Column(s): '",
                                  paste( c("ID", "P.Value", "Score", "high.class")[!geneCheck], collapse = "', '"),
                                  "' missing. Required when SigType = 'gene_rank'." ,
                                  sep = ""))
        }
        if( ncol(object@Signature) !=1 & object@SigType == "gene_list") {
            errMessage <- c(errMessage, "If SigType = 'gene_list', SigType must be a single column.")
        }
        if( length(errMessage) == 0 ) {
            return(TRUE)
        }
        else {
            errMessage <- paste(errMessage, collapse = "\n")
            errMessage <- paste("\nCould not create GeneSig object:", errMessage, sep ="\n")
            return(errMessage)
        }
    }
)

#######################################################################
## function: write.GeneSig
##
## Write a GeneSig file to .gsg file and optional RDS file
##
## .gsg file includes a commented (#) header with tab delimited
##      @Signature slot
##
#######################################################################

write.GeneSig <- function(SigObj, fileName, dir = getwd(), writeRDS FALSE){

  if( class(SigObj) != "GeneSig" ) stop("SigObj must be of class 'GeneSig'.")
  outFile <- file.path(dir, fileName)

  MD <- SigObj@MetaData
  Sig <- SigObj@Signature

  colnames(MD) <- paste("#", "Metadata", sep = "")
  rownames(MD) <- paste("#", rownames(MD), sep = "")

  write.table( MD, paste(outFile, "gsg", sep="."), row.names=TRUE, quote=FALSE, sep = "\t")
  suppressWarnings( write.table( Sig, paste(outFile, "gsg", sep="."), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t", append=TRUE) )

  if (writeRDS) saveRDS(SigObj, file=paste(outFile, "rds", sep = "."))
}


#######################################################################
## function: read.GeneSig
##
## Read .gsg file and format to object of class GeneSig
##
#######################################################################
read.GeneSig <- function(file)
{
  # Get meta data
  metaRows <- readLines(file)[-1]
  metaRows <- metaRows[grepl("^#", metaRows)]
  metaRows <- gsub("#", "", metaRows)
  metaRows <- strsplit(metaRows, "\t")
  MetaData <- data.frame( Info = unlist( lapply(metaRows, function(x) x[2]) ) )
  rownames(MetaData) <- unlist( lapply(metaRows, function(x) x[1]) )

  # Get Signature information
  Signature <- read.table(file, sep="\t", header = T, stringsAsFactors = FALSE)

  # Decide between "gene_rank" and"gene_list"
  if(ncol(Signature)==1) SigType <- "gene_list" else SigType <- "gene_rank"

  SigObj <- GeneSig(MetaData = MetaData,
                    Signature = Signature,
                    SigType = SigType)

  return(SigObj)
}



