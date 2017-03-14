data(tcga.subset.400g.200s)
dat

subtypelevels<-levels(dat$subtype)

grps<-ggheat.make.groups(eset = dat, 
	labelcol = "subtype",  #column name for grouping in pData(eset)
	labelvals = subtypelevels, #values to group on (e.g. factor levels of pData(eset)[, labelcol])
	clustFUN = clust_eset #clustering function for eset
	)

col_legend<-list(subtype = list(col_breaks = subtypelevels,
								col_values = brewer.pal(length(subtypelevels),"Set1"),
								col_labels = subtypelevels))

esetlist<-lapply(grps$esetlist, function(i) scale_row(i))

hmcolors<-function(... ) scale_fill_gradient2(low = "blue", mid = "white",
       high = "red", midpoint = 0, limits=c(-3,3), oob=squish, ...)

p2<-ggheat.continuous.group(esetlist, 
	grps$hclist, 
	grps$hrlist, 
	hmcolors,
	hmtitle = "row-zscore GE",
	col_lab = "subtype", 
	col_legend = col_legend,
	ylabstr = "",
	fout  = NA, 
	p.heights = c(1.5, 0.5, 5),
	xsize = 0,
	ysize = 0, 
	ysizelab = 7,
	xleft = 0.10, 
	xright = 0.24)

