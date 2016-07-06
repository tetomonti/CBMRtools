require(Biobase)
require(CBMRtools)
CBMRtools <- Sys.getenv("CBMRtools")
source(paste(CBMRtools,"CBMRtools/R/diffanal.scores.R",sep="/"))
source(paste(CBMRtools,"CBMRtools/R/ks.score.R",sep="/"))
source(paste(CBMRtools,"CBMRtools/R/misc.R",sep="/"))
source(paste(CBMRtools,"CBMRtools/R/perm.2side.R",sep="/"))
source(paste(CBMRtools,"CBMRtools/R/permute.array.R",sep="/"))
source(paste(CBMRtools,"dvlp/cbmGSEA.R",sep="/"))

gSet <- getGeneSet(new("GeneSet","~/Research/CancerGenomeAnalysisData/annot/h.all.v5.1.symbols.gmt"))
brca <- readRDS("~/Research/Projects/tcga/brca/BRCA_2015_02_04_ES.rds")
exprs(brca) <- log2(exprs(brca)+1)

## prepare GEP data (extract small subset of samples and subset of genes)

brca <- brca[!is.na(fData(brca)$gene_symbol),]
featureNames(brca) <- fData(brca)$gene_symbol
gSet <- gSet[order(sapply(gSet,length))[1:3]]

set.seed(123)
brca1 <- brca[c(unique(unlist(gSet)),
                sample(setdiff(featureNames(brca1),unique(unlist(gSet))),size=100)),
              c(sample(which(pData(brca)$my_stage=="AN"),size=8),
                sample(which(pData(brca)$my_stage=="stage i"),size=8))]
pData(brca1) <- pData(brca1)[,"tissue_type",drop=FALSE]

## t.score returns negative scores for the cond class and positive for the control class (opposite run_limma)
rnkS <- -t.score(x=exprs(brca1),cls=factor(pData(brca1)$tissue_type))

## add two custom "true positive" genesets 
gSet <- c(list(tumor50=rownames(diff)[order(diff[,"t",],decreasing=TRUE)[1:50]],
               normal50=rownames(diff)[order(diff[,"t"],decreasing=FALSE)[1:50]]),
          gSet)

set.seed(123)
OUT1 <- cbmGSEA(eSet=brca1,gSet=gSet,pheno="tissue_type",cond="Tumor",control="AN",nperm=100,verbose=TRUE,plot.name="cbmGSEA.plots1.pdf")

set.seed(123)
OUT2 <- cbmGSEA(rnkScore=-rnkS,gSet=gSet,nperm=100,verbose=TRUE,plot.name="cbmGSEA.plots2.pdf")

all.equal(OUT1[,"score"],OUT2[,"score"])

set.seed(123)
OUT3 <- cbmGSEA(eSet=brca1,tag="MAML2",gSet=gSet,nperm=100,verbose=TRUE,plot.name="cbmGSEA.plots3.pdf")

cbmGSEA.eSet <- brca1
cbmGSEA.gSet <- gSet
save(cbmGSEA.eSet,cbmGSEA.gSet,file=file.path(CBMRtools, "CBMRtools/data/cbmGSEAdata.rda"))
