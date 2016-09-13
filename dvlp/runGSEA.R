#####################################################################################
#' runGSEA
#' 
#' \code{runGSEA} wrapper to call the jar version of Broad's GSEA
#'
#' @param eSet      
#' @param eSet
#' @param covariate
#' @param test
#' @param control
#' @param nperm=1000
#' @param baseName
#' @param gmt.file                geneset compendium file (e.g., c2.cp)
#' @param topX=20                 generate gsea plots for tox X genesets in each direction
#' @param rnkFile=".RANKfile.rnk" temp file where to write the input ranking needed by gsea
#' @param maxG=500                max geneset size
#' @param minG=15                 min geneset size
#' @param jarFile="~/Research/Tools/java/gsea2-2.2.2.jar"
#' @param chip="gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip"
#' @param verbose=TRUE
#'
#' @examples
#'
#' # see end of file
#' # <CODE SNIPPET>
#' 
#' @export
#####################################################################################
runGSEA <- function
(
    eSet,
    covariate,
    test,
    control,
    nperm=1000,
    outdir,                  # where to save the output
    rptLabel,                # prefix label for the output directory
    gmtFile,                 # geneset compendium file (e.g., c2.cp)
    topX=20,                 # generate gsea plots for tox X genesets in each direction
    maxG=500,                # max geneset size
    minG=15,                 # min geneset size
    jarFile="~/Research/Tools/java/gsea2-2.2.2.jar",
    chip="gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip",
    gepFile="expr.txt",
    clsFile="pheno.cls",
    do.rm=TRUE,
    verbose=TRUE
)
{
   if ( !dir.exists(outdir) ) stop( "undefined outdir: ",outdir )
   
   testname <- paste(rptLabel,test,'vs',control,sep='_')
   testname <- gsub(' ','_',testname)

   ## use only the samples the have either test or control classes
   rdata <- eSet[,pData(eSet)[[covariate]] %in% c(test,control)]
   
   ## ordering them in a way that we have first test and then control
   rdata <- rdata[,c(grep(paste('^',test,'$',sep=''),pData(rdata)[[covariate]]),
                   grep(paste('^',control,'$',sep=''),pData(rdata)[[covariate]]))]
   
   ## dump expression matrix into a tab-delimited txt file
   VERBOSE(verbose,"Writing expression data to file:",gepFile)
   write.table(rbind(c('symbols',colnames(rdata)),
                     cbind(featureNames(rdata),exprs(rdata))),
               file=gepFile,
               quote=FALSE,
               sep='\t',
               col.names=FALSE,
               row.names=FALSE)
   VERBOSE(verbose,", done.\n")
   
   ## dump phenotypes into a cls file
   VERBOSE(verbose,"Writing class template to file:", clsFile)
   pheno <- as.character(pData(rdata)[[covariate]])
   con <- file(clsFile,open='w')
   write(paste(length(pheno),'2 1'),con)
   write(paste('# ',test,' ',control,sep=''),con)
   classes <- ''
   for (i in 1:length(pheno)){
      classes <- paste(classes,pheno[i])
   }
   write(classes,con)
   close(con)         
   VERBOSE(verbose,", done.\n")
   
   ## keep track of the directory content (will see why below)
   lsBefore <- if ( dir.exists(outdir) ) list.files(outdir)

   CMD <- paste("java -Xmx3072m -cp ", jarFile,
                " xtools.gsea.Gsea -res ", gepFile, " -cls ", clsFile, "#", test, "_versus_", control," -gmx ", gmtFile,
                " -collapse false -nperm ", nperm,
                " -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label ", testname, 
                " -metric Signal2Noise -sort real -order descending -include_only_symbols false",
                " -make_sets true -median false -num 100 -plot_top_x ", topX,
                " -rnd_seed timestamp -save_rnd_lists false -set_max ", maxG, " -set_min ", minG, 
                " -zip_report false -out ", outdir, " -gui false", sep="")
   
   ## call java gsea version
   VERBOSE(verbose,CMD,"\n")
   system(CMD)
   if (do.rm) unlink(c(clsFile,gepFile))

   ## determine the name of the generated directory
   lsAfter <- list.files(outdir)
   OUT <- setdiff(lsAfter,lsBefore)
   if ( length(OUT)>1 ) stop( "more than one output directories: ", paste(OUT,sep=",") )
   VERBOSE(verbose,"output saved to:",OUT,"\n")
   
   ## determine the names and read the tabular results files' content, to be returned
   posFile <- system( paste("ls ",outdir,"/",OUT,"/gsea_report_for_",test,"_*.xls",sep=""), intern=TRUE )
   negFile <- system( paste("ls ",outdir,"/",OUT,"/gsea_report_for_",control,"_*.xls",sep=""), intern=TRUE )
   VERBOSE(verbose,"posFile:",posFile,"\n")
   VERBOSE(verbose,"negFile:",negFile,"\n")
   POS <- read.delim(posFile,check.names=FALSE,stringsAsFactors=FALSE)
   NEG <- read.delim(negFile,check.names=FALSE,stringsAsFactors=FALSE)
   return( list(pos=POS,neg=NEG) )
}
#####################################################################################
#' runGSEApreranked
#' 
#' \code{runGSEApreranked} wrapper to call the jar version of Broad's GSEApreranked
#'
#' @param ranking      2-column (gene IDs and scores) matrix
#' @param gmtFile      geneset compendium file (e.g., c2.cp)
#' @param outdir       where to save the output
#' @param rptLabel     prefix label for the output directory
#' @param nperm        number of permutation iterations
#' @param topX         generate gsea plots for tox X genesets in each direction
#' @param maxG         max geneset size
#' @param minG         min geneset size
#' @param rnkFile      temp file where to write the input ranking needed by gsea
#' @param jarFile      file with the jar executable
#' @param chip         '.chip' file of probeset-to-genesymbol mapping
#' @param verbose
#'
#' @examples
#'
#' # see end of file
#' # <CODE SNIPPET>
#' 
#' @export
#' 
#####################################################################################

runGSEApreranked <- function
(
    ranking,                 # 2-column (gene IDs and scores) matrix
    gmtFile,                 # geneset compendium file (e.g., c2.cp)
    outdir,                  # where to save the output
    rptLabel,                # prefix label for the output directory
    nperm=1000,              # number of permutation iterations
    topX=20,                 # generate gsea plots for tox X genesets in each direction
    maxG=500,                # max geneset size
    minG=15,                 # min geneset size
    rnkFile=".RANKfile.rnk", # temp file where to write the input ranking needed by gsea
    jarFile="~/Research/Tools/java/gsea2-2.2.2.jar",
    chip="gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip",
    verbose=TRUE
)
{
    if ( !dir.exists(outdir) ) stop( "undefined outdir: ",outdir )
    
    if ( ncol(ranking)!=2 ) stop( "ranking must be a 2-column matrix:", ncol(ranking) )
    if ( !is.numeric(ranking[,2]) ) stop ( "2nd column of ranking must be numeric:", class(ranking[,2]) )

    ## write the ranking data to a temp file
    my.write.table(ranking,rnkFile,row.names=FALSE)

    ## keep track of the directory content (will see why below)
    lsBefore <- if ( dir.exists(outdir) ) list.files(outdir)

    ## run gsea
    CMD <- paste("java -Xmx3072m  -cp", jarFile,
                 "xtools.gsea.GseaPreranked -gmx", gmtFile,
                 "-collapse false -mode Max_probe -norm meandiv -nperm", nperm, "-rnk", rnkFile,
                 "-scoring_scheme weighted -rpt_label", rptLabel,
                 "-chip", chip,
                 "-include_only_symbols true -make_sets true -plot_top_x", topX,
                 "-rnd_seed timestamp -zip_report false -set_max", maxG,"-set_min", minG,
                 "-out", outdir, "-gui false", sep=" ")
    VERBOSE(verbose,CMD,"\n")
    system( CMD )

    ## remove the temp file
    file.remove(rnkFile)

    ## determine the name of the generated directory
    lsAfter <- list.files(outdir)
    OUT <- setdiff(lsAfter,lsBefore)
    if ( length(OUT)>1 ) stop( "more than one output directories: ", paste(OUT,sep=",") )
    
    VERBOSE(verbose,"output saved to:",OUT,"\n")

    ## determine the names and read the tabular results files' content, to be returned
    posFile <- system( paste("ls ",outdir,"/",OUT,"/gsea_report_for_na_pos_*.xls",sep=""), intern=TRUE )
    negFile <- system( paste("ls ",outdir,"/",OUT,"/gsea_report_for_na_neg_*.xls",sep=""), intern=TRUE )
    VERBOSE(verbose,"posFile:",posFile,"\n")
    VERBOSE(verbose,"negFile:",negFile,"\n")
    POS <- read.delim(posFile,check.names=FALSE,stringsAsFactors=FALSE)
    NEG <- read.delim(negFile,check.names=FALSE,stringsAsFactors=FALSE)
    return( list(pos=POS,neg=NEG) )
}

##################################################################
## code to create and test the data objects used in the examples #
##################################################################
if ( FALSE )
{
    CBMRtools <- Sys.getenv("CBMRtools")
    source(file.path(CBMRtools,"CBMRtools/R/diffanal.scores.R"))
    source(file.path(CBMRtools,"dvlp/runGSEA.R"))

    ## use the toy data created for cbmGSEA
    ##
    load(file.path(CBMRtools,"CBMRtools/data/cbmGSEAdata.rda"))
    if ( is.null(cbmGSEA.eSet) || is.null(cbmGSEA.gSet) ) stop( "missing data objects" )

    ## write toy geneset compendium to file
    ##
    write.gmt(cbmGSEA.gSet,gmtfile="gsets.gmt")

    ## testing full-fledged GSEA
    ##
    OUT1 <- runGSEA(eSet=cbmGSEA.eSet, covariate="tissue_type",test="Tumor",control="AN",
                    outdir="test",rptLabel="testGSEA", gmtFile="gsets.gmt")

    ## testing GSEApreranked
    ##
    rnkS <- -t.score(x=exprs(cbmGSEA.eSet),cls=factor(pData(cbmGSEA.eSet)$tissue_type))
    rnkS <- data.frame(GeneID=names(rnkS),score=rnkS)
    OUT2 <- runGSEApreranked( ranking=rnkS, gmtFile="gsets.gmt", outdir="test", rptLabel="testPreranked" )
}
