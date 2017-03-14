#example for usage of CBMRtools::ggheat.continuous.single
library(CBMRtools)

#load expression data
data(tcga.subset.400g.200s)
dat

#hclust for rows and columns
hc<-clust_eset(dat)

pData(dat)$hclust.groups<-as.factor(cutree(hc$hc, k = 4))

#scale expression by row
dat.scaled<-scale_row(dat)

subtypelevels<-levels(dat$subtype)

#color legends for column labels
col_legend<-list(subtype = list(col_breaks = subtypelevels,
								col_values = brewer.pal(length(subtypelevels),"Set1"),
								col_labels = subtypelevels), 
		hclust.groups = list(col_breaks = levels(pData(dat)$hclust.groups),
								col_values = sapply(c("pink", "orange", "yellow", "cyan"), to.hex),
								col_labels = levels(pData(dat)$hclust.groups)))

#heatmap fill gradient
hmcolors<-function(... ) scale_fill_gradient2(low = "blue", mid = "white",
       high = "red", midpoint = 0, limits=c(-3,3), oob=squish, ...)

p<-ggheat.continuous.single(eset = dat.scaled, 
	hc = hc$hc, 
	hr = hc$hr, 
	hmcolors = hmcolors,
	hmtitle = "row-zscore GE",
	col_lab = c("subtype", "hclust.groups"), 
	col_legend = col_legend,
	ylabstr = "",
	fout = NA, 
	p.heights = c(1.5, 0.5, 5),
	xsize = 0,
	ysize = 0, 
	ysizelab = 7,
	xright = 0.18)