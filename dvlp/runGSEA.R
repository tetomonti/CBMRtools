runGSEA <- function
(
    eSet,
    covariate,
    test,
    control,
    nperm=1000,
    base.name,
    gmt.file,
    topX=20,                 # generate gsea plots for tox X genesets in each direction
    rnkFile=".RANKfile.rnk", # temp file where to write the input ranking needed by gsea
    maxG=500,                # max geneset size
    minG=15,                 # min geneset size
    jarFile="~/Research/Tools/java/gsea2-2.2.2.jar",
    chip="gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip",
    verbose=TRUE
)
{
   testname<-paste(base.name,test,'vs',control,sep='_')
   testname<-gsub(' ','_',testname)

   ## use only the samples the have either test or control classes
   rdata<-eSet[,pData(eSet)[[covariate]]%in%c(test,control)]
   
   ## ordering them in a way that we have first test and then control
   rdata<-rdata[,c(grep(paste('^',test,'$',sep=''),pData(rdata)[[covariate]]),
                   grep(paste('^',control,'$',sep=''),pData(rdata)[[covariate]]))]
   
   ## dump expression matrix into a tab-delimited txt file
   write.table(rbind(c('symbols',colnames(rdata)),
                     cbind(featureNames(rdata),exprs(rdata))),
               file='expr.txt',
               quote=FALSE,
               sep='\t',
               col.names=FALSE,
               row.names=FALSE)
   
   ## dump phenotypes into a cls file
   pheno<-as.character(pData(rdata)[[covariate]])
   con<-file('pheno.cls',open='w')
   write(paste(length(pheno),'2 1'),con)
   write(paste('# ',test,' ',control,sep=''),con)
   classes<-''
   for (i in 1:length(pheno)){
      classes<-paste(classes,pheno[i])
   }
   write(classes,con)
   close(con)         
   
   ## call java gsea version
   system(paste("java -Xmx3072m -cp ", jarFile,
                " xtools.gsea.Gsea -res expr.txt -cls pheno.cls#", test, "_versus_",c ontrol," -gmx ", gmtFile,
                " -collapse false -nperm ", nperm,
                " -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label ", testname, 
                " -metric Signal2Noise -sort real -order descending -include_only_symbols false",
                " -make_sets true -median false -num 100 -plot_top_x ", topX,
                " -rnd_seed timestamp -save_rnd_lists false -set_max ", maxG,
                "-set_min 5 -zip_report false -out htmls -gui false",sep=""))
   unlink(c("pheno.cls","expr.txt"))
   
}
runGSEApreranked <- function
(
    ranking,                 # 2-column (IDs and scores) matrix
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
    if (ncol(ranking)!=2) stop( "ranking must be a 2-column matrix:", ncol(ranking) )

    ## write the ranking data to a temp file
    my.write.table(ranking,rnkFile,row.names=FALSE)

    ## keep track of the directory content (will see why below)
    lsBefore <- if ( dir.exists(outdir) ) system(paste("ls",outdir),intern=TRUE)

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
    lsAfter <- system(paste("ls",outdir),intern=TRUE)
    OUT <- setdiff(lsAfter,lsBefore)
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
