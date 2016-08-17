
# Create class "GeneSig"
GeneSig <- setClass("GeneSig", 
                    slots = c(SigType = "character",
                              MetaData = "data.frame",
                              Signature = "data.frame"),
                    validity = function(object){
                      errMessage = c()
                      if( sum( c("gene_rank", "gene_list") %in% object@SigType) != 1) errMessage <- c(errMessage,
                                                                                     "SigType must be on of: 'gene_rank' or 'gene_list'.")
                      if( ncol(object@MetaData) != 1 ) errMessage <- c( errMessage, 
                                                                      "MetaData must be 1 column." )
                      geneCheck <- c("ID", "P.Value", "Score") %in% colnames(object@Signature)
                      if( !all(geneCheck) & object@SigType == "gene_rank" ) errMessage <- c( errMessage ,
                                                            paste( "Column(s): '", 
                                                                  paste( c("ID", "P.Value", "Score")[!geneCheck],
                                                                        collapse = "', '"), 
                                                                  "' missing. Required when SigType = 'gene_rank'." ,
                                                                  sep = "") )
                      if( ncol(object@Signature) !=1 & object@SigType == "gene_list") errMessage <- c(errMessage, 
                                                                                                               "If SigType = 'gene_list', SigType must be a single column.")
                      if( length(errMessage) == 0 ){
                        return(TRUE)
                      } else {
                        errMessage <- paste(errMessage, collapse = "\n")
                        errMessage <- paste("\nCould not create GeneSig object:", errMessage, sep ="\n")
                        return(errMessage)
                      }
                      
                    }
                    )


# Create file that writes GeneSig to text file
write.GeneSig <- function(SigObj, fileName, dir = getwd()){
  
  if( class(SigObj) != "GeneSig" ) stop("SigObj must be of class 'GeneSig'.")
  outFile <- file.path(dir, fileName)
  
  MD <- SigObj@MetaData
  Sig <- SigObj@Signature
  
  colnames(MD) <- paste("#", "Metadata", sep = "")
  rownames(MD) <- paste("#", rownames(MD), sep = "")
  
  write.table( MD, outFile, row.names = TRUE, quote = FALSE, sep = "\t")
  suppressWarnings( write.table( Sig, outFile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", append = TRUE) )

  }

# Read in txt file
read.GeneSig <- function(file){
  
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
  
  

