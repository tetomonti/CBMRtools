
library(Rtsne)
library(reshape2)
library(ggplot2)

run_Rtsne<-function(mat, col, col.color = NA, col.label = "type", k = 5, 
	strmain = "", seed = NA, fout){
	
	if(!is.na(seed))
		set.seed(seed)

	res<-Rtsne(X = mat, dims = k)$Y
	colnames(res)<-paste("tSNE", 1:ncol(res), sep = "")
	res.m<-melt(res)
	res.m<-merge(x = res.m, y = res.m, by = "Var1")
	col<-data.frame(col)
	colnames(col)<-col.label
	obj<-cbind(res.m, col[res.m$Var1,col.label])
	colnames(obj)<-c(colnames(res.m), col.label)

	p<-ggplot(obj, aes_string(x="value.x", 
					y="value.y", color=col.label))+
				geom_point(size = 1, alpha = 0.5) + theme_bw() + 
				facet_grid(Var2.x~Var2.y)+
				ggtitle(strmain)+ 
				theme(plot.title = element_text(size = 15, face = "bold"))+
				xlab("")+ylab("")
	
	isna<-function(i){	
		if(length(i)>1) return(FALSE)
		else return(is.na(i))
	}

	if(!isna(col.color))
		p<-p+scale_colour_manual(values = col.color, breaks = names(col.color))

	ggsave(file = fout, p, width = 20, height = 20, unit = "in")
	return(obj)
}