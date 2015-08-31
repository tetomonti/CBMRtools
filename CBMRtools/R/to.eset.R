#' construct an eSet object 
#' @import Biobase
#' @param mat expression matrix
#' @param pdat phenoData data frame
#' @param fdat featureData data frame
#' @export

to.eSet<-function(mat, pdat, fdat){
	#require(Biobase)
	mat<-as.matrix(mat)

	#checking data type and dimensions
	if (!is.data.frame(pdat))
		stop("pdat must be a data frame")
	
	if (!is.data.frame(fdat))
		stop("fdat must be a data frame")
	
	if ( nrow(fdat) != nrow(mat))
		stop("nrow(fdat) must equal nrow(mat)")
	
	if( nrow(pdat) != ncol(mat))
		stop("nrow(pdat) must equal ncol(mat)")

	if (!all(rownames(fdat) == rownames(mat))){
		warning("fdat rownames and mat rownames do not match, setting fdat rownames to mat rownames")
		rownames(fdat) <- rownames(mat)
	}

	if (!all(rownames(pdat) == colnames(mat))){
		warning("pdat rownames and mat colnames do not match, setting fdat rownames to mat rownames")
		rownames(pdat) <- colnames(mat)
	}

	fMetaData<-data.frame(labelDescription = colnames(fdat), row.names = colnames(fdat))
	featureData<-new("AnnotatedDataFrame", data= fdat, varMetadata=fMetaData) 

	pMetaData<-data.frame(labelDescription = colnames(pdat), row.names = colnames(pdat))
	phenoData<-new("AnnotatedDataFrame", data= pdat, varMetadata=pMetaData) 
	
	eSet<-ExpressionSet(assayData=mat, featureData = featureData, phenoData = phenoData,  annotation = "")
	return(eSet)
}