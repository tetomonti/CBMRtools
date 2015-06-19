
#' \code{save.xlsx} save to .xlsx  
#' @import XLConnect
#' @param x R object to save, must be a named list of data.frames, row limit =65535 for each data frame
#' @param f.dir directory to save to, default = NA, current working directory
#' @param f.header file header
#' @param f.name file name (appended to file header)
#' @param java.param max java memory, default =  "-Xmx3g"
#' @export

save.xlsx<-function(x, f.dir = NA, f.header = "", f.name, java.param = "-Xmx3g"){
	options(java.parameters = java.param )
	#require(XLConnect)
	# Load workbook (create if not existing)
	if(is.na(f.dir)) 
		wb.name<-paste(f.name, ".xlsx", sep = "")
	else
		wb.name<-paste(f.dir, "/", f.header, f.name, ".xlsx", sep = "")
	
	if (file.exists(wb.name)){
		file.remove(wb.name)
	}
	
	wb <- loadWorkbook(wb.name, create = TRUE)
	
	for (i in names(x)){ 
		createSheet(wb, name = i)
		writeWorksheet(wb, x[[i]], sheet = i, startRow = 1, startCol = 1)
	}
	saveWorkbook(wb)
	xlcFreeMemory()
}

