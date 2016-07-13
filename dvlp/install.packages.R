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
          "xlsx",
          "dynamicTreeCut",
          "cba",
          "devtools",
          "roxygen2",
          "rmarkdown",
          "XLConnect",
          "data.table",
          "dendextend")
install.packages(pkgs,repos="http://cran.r-project.org")
require(devtools)
install_github('hadley/staticdocs')

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("biomaRt","ROC","pathifier","ConsensusClusterPlus","ASSIGN","Biobase","oligo","oligoClasses",
           "limma","frma","GSEAlm"))
biocLite("DESeq2","edgeR")
