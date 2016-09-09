require(openxlsx)

qmatrix <- readRDS("qmatrix.RDS")
nCol <- ncol(qmatrix)+1
nRow <- nrow(qmatrix)+1
wb <- createWorkbook()
addWorksheet(wb,"sheet1")
writeData(wb,"sheet1",qmatrix,startCol=1,startRow=1,rowNames=TRUE)

## creating cell styles 
style05 <- createStyle(fontColour = "black", bgFill = "pink")
style01 <- createStyle(fontColour = "black", bgFill = "red")
style10 <- createStyle(fontColour = "black", bgFill = "white")

## conditional formatting based on q-values

conditionalFormatting(wb, "sheet1",cols=2:nCol,rows=2:nRow,rule=">=0",style=style01)
saveWorkbook(wb,"test.xlsx",overwrite=TRUE) # check the file: correct (all red, as expected)

conditionalFormatting(wb, "sheet1",cols=2:nCol,rows=2:nRow,rule=">0.01",style=style05)
saveWorkbook(wb,"test.xlsx",overwrite=TRUE) # check the file: correct (red and pink, as expected)

conditionalFormatting(wb, "sheet1",cols=2:nCol,rows=2:nRow,rule=">0.05",style=style10)
saveWorkbook(wb,"test.xlsx",overwrite=TRUE) # check the file: incorrect (pink is gone and only red and white left)
## there's a bug in 'conditionalFormatting'
