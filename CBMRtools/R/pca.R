require(Biobase)

pca_plot<-function(eset, label, colors, madF = NA){

	if (!is.na(madF)){

		eset<--variationFilter(eset, 
			score= "mad", dir = "top", ngenes = madF)
	}

	emat<-exprs(eset)
	emat.pca<-prcomp(t(emat))
	res<-data.frame(emat.pca$x)
	res<-data.frame(t(t(res)/apply(res,2,max)))
	res$class<-pData(eset)[, label]
	res$class <- factor(res$class)

	#colors.n<-length(unique(as.character(res$class)))
	#Naïve (white), Vh (green), Rosi (Red), TBT (Blue) and TPP (purple)


	p<-ggplot(res, aes(x=PC1, y = PC2, color=class, size = 3))+
	geom_point(size = 3) + theme_bw() + ggtitle("PCA Analysis by Sample Day and Treatment")
	
	# +
	#scale_colour_manual(values = colors)

	return(p)
}