#' wrapper script for microarray gene ST chips processing into expression set object 
#'
#' @import Biobase oligo oligoClasses
#' @param indir input directory
#' @param outdir output directory
#' @param header of outfile
#' @return the normalized expressionSet object
#' 
#' @examples
#' #running example
#' #indir<-"/Users/amyli/Desktop/temp_dir/example_script/raw"
#' #outdir<-"/Users/amyli/Desktop/temp_dir/example_script/processed"
#' #dat<-oligo_wrapper(indir = indir, outdir = outdir, header = "my_eset")
#'
#' 
#' @export 

#require(oligo)

oligo_wrapper<-function(indir, outdir, header){
	celFiles <- list.celfiles(indir, full.names = TRUE)
	affyRaw <- read.celfiles(celFiles)
	eset <- rma(affyRaw, target = "core")
	fileout <- paste(outdir, "/", header, ".RDS", sep = "")
	cat(paste("Saving to ", fileout, "\n", sep = ""))
	saveRDS(eset, file = fileout)
	cat("Done!\n")
	return(eset)
}

#indir<-"/Users/amyli/Desktop/temp_dir/example_script/raw"
#outdir<-"/Users/amyli/Desktop/temp_dir/example_script/processed"
#dat<-oligo_wrapper(indir = indir, outdir = outdir, header = "my_eset")



