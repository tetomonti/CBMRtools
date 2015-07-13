###this script reduces expression set objects to rows with unique transcript IDS

##single always fetches the probe with maximum of summary function (e.g. probe with maximum according to MAD statistic)
#' \code{collapse.probes} collapse eSet object by non-unique row feature IDs
#' @import Biobase data.table
#' @param eSet expression set object
#' @param fdat.id feature id to collapse on
#' @param method = "single" or "summarize" followed by "." followed by "mean", "median", "max" or "min", e.g. "single.mean"
#'  "single" chooses one probe for ambiguous probes with maximum of summary function across row, 
#'  "summarize" aggregates summary metric across multiple rows for each column index
#' @export

collapse.probes<-function(eset, 
	fdat.id = "SYMBOL",
	method = c("single.mad",
	"summarize.mean", 
	"summarize.median", 
	"summarize.max",
	"summarize.min")){

	library(data.table)
#	library(Biobase)
	#aggregate information across multiple "probes" using summary function
	if(grepl("summarize", method)){
		x<-data.table(exprs(eset))
		x$y<-fData(eset)[, fdat.id]

		fxn<-strsplit(method, split = "\\.")[[1]][2]
		summarize.valid<-c("mean", "median", "max", "min")
		
		if(!(fxn %in% summarize.valid)) 
		stop (paste("invalid criteria, summarize function must be one of ",
			paste(summarize.valid, collapse = ","), sep = ""))

		cmdstr<-paste("x[, lapply(.SD,", fxn, 
			"), by= y]")

		df.summary<-eval(parse(text =cmdstr))
		df.summary<-data.frame(df.summary)

		df.id<-df.summary$y
		df.exprs<-df.summary[, colnames(df.summary) != "y"]

		fdat<-fData(eset)
		fdat.unique.id<-make.unique(as.character(fdat[, fdat.id]))
		new.eset<-eset[match(df.id, fdat.unique.id), ]
		exprs(new.eset)<-as.matrix(df.exprs)
		return(new.eset)
	}

	#using single probe per unique probe id with criteria function
	else {
		x<-exprs(eset)
		fxn<-strsplit(method, split = "\\.")[[1]][2]

		single.valid<-c("median", "max", "min", "mad")
		
		if(!(fxn %in% single.valid)) 
		stop (paste("invalid criteria, single function must be one of ",
			paste(single.valid, collapse = ","), sep = ""))

		cmdstr<-paste("apply(x, 1,", fxn, ", na.rm = TRUE",
			")")
		dat<-data.table(data.frame(id = fData(eset)[, fdat.id],
			summary = eval(parse(text =cmdstr))))
		dat$ind<-1:nrow(dat)
	
		dat.reduced<-dat[dat[, .I[summary == max(summary)][1] , by=id][,V1]]
		summary.idx<-dat.reduced$ind

		return(eset[summary.idx,])
	}
}