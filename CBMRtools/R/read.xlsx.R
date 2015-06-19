
#' \code{read.xlsx} read from .xlsx  
#' @import XLConnect
#' @param file the xlsx file to read from, may contain multiple sheets
#' @param java.param max java memory, default =  "-Xmx3g"
#' @export
read.xlsx <- function( file, java.param = "-Xmx3g")
{
	options(java.parameters = java.param )
	wb <- loadWorkbook(file)
	wb.names <- getSheets(wb)
	res<-list()
	for (i in wb.names){
		res[[i]]<-readWorksheet(wb, sheet = i)
	}
	xlcFreeMemory()
	return(res)
}

