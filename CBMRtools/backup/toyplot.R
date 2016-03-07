# source("~/Desktop/git_projects/montilab/CBMRtools/CBMRtools/backup/toyplot.R")


source("/Users/amyli/Desktop/git_projects/montilab/CBMRtools/CBMRtools/backup/heatmap.ggplot.discrete.R")

library(CBMRtools)

f.dir.in<-"/Users/amyli/Desktop/monti_lab/year3/network_analysis/data/data_dump_2.2.16/res_ok"
f.dir.out<-"/Users/amyli/Desktop/monti_lab/year3/network_analysis/results/manuscript_network_2.2.16_results"


load(paste(f.dir.in, "/", "data_RI.RData", sep = ""))

load(paste(f.dir.in, "/", "groups_colors_RI.RData", sep = ""))


##make data into expression set

col_carc[col_carc == "purple"]<-"Carcinogen"
col_carc[col_carc == "yellow"]<-"Noncarcinogen"
col_carc[col_carc == "white"]<-NA


col_geno[col_geno == "green"]<-"Nongenotoxic"
col_geno[col_geno == "pink"]<-"Genotoxic"
col_geno[col_geno == "white"]<-NA


pheno<-data.frame(carcinogenicity = col_carc, genotoxicity = col_geno, chem = ch_ok, group = sel_groups)
dat<-to.eSet(mat = rand_matrix_adj, pdat =pheno, fdat = pheno)

mat.orig<-exprs(dat)


mat<-exprs(dat)
mat<-log10(mat)
datnum<-as.numeric(mat)
mat<-replace(mat, mat == 0, max(datnum[datnum != 0]))
exprs(dat)<-mat

meta.c.color.string<-c("yellow", "green", "pink", "purple")
meta.c.color<-as.character(sapply(meta.c.color.string, to.hex))
names(meta.c.color)<-c("Noncarcinogen", "Nongenotoxic", "Genotoxic", "Carcinogen")


m2.colors<-c(
"black" ,     "blue",        "brown",       "green" ,      "greenyellow",
 "magenta" ,    "pink",        "purple",      "red" ,        "salmon"    , 
"tan" ,        "turquoise" ,  "yellow"    )

m2.groups<-c(5, 9, 1, 13, 4,
	7, 3, 6, 2, 11,
	10, 8, 12)
m2.groups<-paste("G", m2.groups, sep = "" )

dat$groups<-factor(sapply(1:length(dat$group), function(i){
	m2.groups[which(m2.colors == dat$group[i])]
	}))

dat$groups <- factor(dat$groups, 
	levels = paste("G", 1:13, sep = ""))

m2<-sapply(m2.colors, to.hex)
names(m2)<-m2.groups


HCM3<-HCM2
HCM3$height <-HCM3$height - min(HCM3$height) +0.01


jpeg(paste(f.dir.out, "RI_heatmap_reds_discrete.jpg", sep = "/"),
	, width = 4, height = 6, units = 'in', res = 300)



discrete.colors.breaks<-sapply(seq(0, 1, by = 1/5),
	function(i) quantile(as.numeric(exprs(dat)), i))
discrete.colors.breaks[1]<-(-Inf)
discrete.colors.breaks[length(discrete.colors.breaks)]<-(Inf)


labels<-sapply(seq(0, 1, by = 1/5),
	function(i) round(quantile(as.numeric(exprs(dat)), i),2))

sndlast<-labels[length(labels)-1]
last<-labels[length(labels)]
maxchar<-12

newlast<-paste(labels[length(labels)-1], 
	paste(rep(" ", maxchar - nchar(sndlast) - nchar(last)), collapse = ""), 
	labels[length(labels)])


labels2<-c(sapply(labels[1:(length(labels)-2)], function(i){
	paste(i, paste(rep(" ", maxchar - nchar(i)), collapse = ""), sep = "")
	}) , newlast)

#discrete.colors.breaks = c(-Inf, -3, -2.5, -2, -1.5, -1 ,Inf)

discrete.colors.values =  colorRampPalette(brewer.pal(9,"Reds"))(5)
discrete.colors.labels = as.character(labels2)

curr.heights = c(1, 1.5, 0.8, 6.5, 1.6, 1.5, 1)

grid.newpage()
p1<-heatmap.ggplot2.discrete(
	eSet = dat[order.dendrogram(as.dendrogram(HCM3)),],
	#eSet= dat[rev(order.dendrogram(as.dendrogram(HCM3))),], 
 	col.clust = TRUE, row.clust = FALSE,
    col.clust.hc = HCM3, row.clust.hc = NA,
    col.legend.brewer = c(meta.c.color),
    row.legend.brewer = "",

    discrete.colors = TRUE,
	discrete.colors.breaks = discrete.colors.breaks,
	discrete.colors.values =  discrete.colors.values,
	discrete.colors.labels = discrete.colors.labels,
    col.lab = c("genotoxicity", "carcinogenicity"), 
    row.lab = "",
    heatmap.y.text = TRUE, heatmap.x.text = TRUE,
    heatmap.y.text.size = 3,
    heatmap.x.text.size = 3,
    heatmap.colorlegend.name = "log10_RI",
    title.text = "",
    col.legend.name =c("genotoxicity", "carcinogenicity"),
    row.legend.name = "",
    cuttree.col =0, cuttree.row =0,
    verbose = FALSE, show = FALSE,
    grid.heights = curr.heights/sum(curr.heights))

 	
grid.arrange(p1)
dev.off()


jpeg(paste(f.dir.out, "RI_heatmap_reds_discrete_with_groups.jpg", sep = "/"),
	, width = 4, height = 6, units = 'in', res = 300)


grid.newpage()
p2<-heatmap.ggplot2.discrete(
	eSet = dat[order.dendrogram(as.dendrogram(HCM3)),],
	#eSet= dat[rev(order.dendrogram(as.dendrogram(HCM3))),], 
 	col.clust = TRUE, row.clust = FALSE,
    col.clust.hc = HCM3, row.clust.hc = NA,
    col.legend.brewer = c(meta.c.color, m2),
    row.legend.brewer = "",

    discrete.colors = TRUE,
	discrete.colors.breaks = discrete.colors.breaks,
	discrete.colors.values =  discrete.colors.values,
	discrete.colors.labels = discrete.colors.labels,
    col.lab = c("genotoxicity", "carcinogenicity", "groups"), 
    row.lab = "",
    heatmap.y.text = TRUE, heatmap.x.text = TRUE,
    heatmap.y.text.size = 3,
    heatmap.x.text.size = 3,
    heatmap.colorlegend.name = "log10_RI",
    title.text = "",
    col.legend.name =c("genotoxicity", "carcinogenicity", "groups"),
    row.legend.name = "",
    cuttree.col =0, cuttree.row =0,
    verbose = FALSE, show = FALSE,
    grid.heights = curr.heights/sum(curr.heights))

grid.arrange(p2)
dev.off()

