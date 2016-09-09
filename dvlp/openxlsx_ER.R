require(openxlsx)

qmatrix <- readRDS("qmatrix.RDS")
nCol <- ncol(qmatrix)+1
nRow <- nrow(qmatrix)+1
wb <- createWorkbook()
addWorksheet(wb,"sheet1")
writeData(wb,"sheet1",qmatrix,startCol=1,startRow=1,rowNames=TRUE)


## cell colors and p-value cutoffs
bgFill = c("red", "pink", "white")
bgVal = c(0, 0.01, 0.05)

## for each color and cutoff, fill in excell sheet accordingly
for(i in 1:length(bgFill)){
  colInd <- which(qmatrix >= bgVal[i], arr.ind = TRUE)+1
  
  if(length(colInd) > 0){
      addStyle(wb, "sheet1", style = createStyle(fgFill = bgFill[i]), rows = colInd[,1], cols = colInd[,2])
  }
}

saveWorkbook(wb,"test.xlsx",overwrite=TRUE)
