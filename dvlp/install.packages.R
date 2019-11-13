pkgs <- c("cba",
          "caret",
          "combinat",
          "data.table",
          "dendextend",
          "devtools",
          "doMC",
          "dplyr",
          "dynamicTreeCut",
          "e1071",
          "ff",
          "flexdashboard",
          "GEOquery",
          "getopt",
          "ggdendro",
          "ggplot2",
          "ggpubr",
          "glmnet",
          "gplots",
          "gridExtra",
          "kableExtra",
          "heatmaply",
          "heatmap.plus",
          "manhattanly",
          "mclust",
          "msigdbr",
#          "NanoStringNorm",
          "openxlsx",
          "optparse",
          "pamr",
          "pheatmap",
          "PCIT",
          "plotly",
          "Rmisc",
          "randomForest",
          "randomForestSRC",
          "rgl",
          "rjags",
          "rjson",
          "rmarkdown",
          "roxygen2",
          "rvest",
          "statmod",
          "tsne",
          "VennDiagram",
          "webshot"
          )

install.packages(pkgs,repos="http://cran.r-project.org")
install.packages(pkgs,repos="http://cran.r-project.org",type="source")

## Few packages not available through CRAN
require(devtools)
install_github('hadley/staticdocs')
install_github('hadley/pkgdown')
## if using R version < 3.5
devtools::install_github("montilab/hypeR")

## bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos="http://cran.r-project.org")
BiocManager::install()
BiocManager::install(c("GEOquery","Biobase","biomaRt","ROC","ConsensusClusterPlus","ASSIGN",
                       "limma","frma","DESeq2","edgeR","GSVA","RNASeqPower","pdInfoBuilder"))


## old version
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("biomaRt","ROC","pathifier","ConsensusClusterPlus","ASSIGN","Biobase","oligo",
           "oligoClasses","limma","frma"))
biocLite("limma")
biocLite("DESeq2")
biocLite("edgeR")
biocLite("SCAN.UPC")
biocLite("GSVA")
biocLite(c("car","pvclust","infotheo","RTN"))
biocLite("igraph")
biocLite("RTN")
biocLite('pdInfoBuilder')
biocLite('RNASeqPower')

