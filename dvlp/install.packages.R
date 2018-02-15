pkgs <- c("mclust",
          "getopt",
          "optparse",
          "randomForest",
          "randomForestSRC",
          "pamr",
          "e1071",
          "combinat",
          "heatmap.plus",
#          "missForest",
          "rgl",
          "ggplot2",
          "ggdendro",
          "gridExtra",
          "NanoStringNorm",
          "gplots",
          "rjags",
#          "h5r",
#          "RCurl",
#          "WGCNA",
          "PCIT",
          "rjson",
          "dynamicTreeCut",
          "cba",
          "devtools",
          "roxygen2",
          "rmarkdown",
#          "xlsx",
          "XLConnect",
          "openxlsx",
          "data.table",
          "dendextend",
          "pheatmap")
install.packages(pkgs,repos="http://cran.r-project.org")
require(devtools)
install_github('hadley/staticdocs')
install_github('hadley/pkgdown')

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("biomaRt","ROC","pathifier","ConsensusClusterPlus","ASSIGN","Biobase","oligo",
           "oligoClasses","limma","frma","GSEAlm"))
biocLite("DESeq2","edgeR")
biocLite("SCAN.UPC")
biocLite("GSVA")
