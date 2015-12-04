require(gplots)
require(Biobase)
require(RColorBrewer)
source('heatmap.2g.R')


plotDiagnostic<-function(bwh,unmc,signature,mode){
   #reduce to signature genes
   signature<-sort(signature,decreasing=F)
   bwh<-bwh[names(signature),]
   unmc<-unmc[names(signature),]
   
   bwh$CLASS<-as.factor(as.character(bwh$CLASS))
   
   #extract all relevant scores and classifications
   bwh<-bwh[,order(bwh[[paste0('Diagnostic_score_M',mode)]],decreasing = F)]
   bwh_scores<-bwh[[paste0('Diagnostic_score_M',mode)]]
   bwh_labels<-bwh[[paste0('Diagnostic_M',mode)]]
   
   #get the unmc scores in the same order
   unmc_scores<-unmc[[paste0('Diagnostic_score_M',mode)]]
   unmc_labels<-unmc[[paste0('Diagnostic_M',mode)]]
   names(unmc_scores)<-names(unmc_labels)<-unmc$id
   unmc_scores<-unmc_scores[bwh$id]
   unmc_labels<-unmc_labels[bwh$id]
 
   #use the sample names instead of the file names
   colnames(bwh)<-bwh$id   

   #scatterplot
   plot(bwh_scores,unmc_scores,pch=19,main='Correlation of the scores between BWH/UNMC')

   #graphics
   COLS=c('purple','red','black','blue','green','#FF33FF','yellow','grey')
   COLS<-COLS[1:length(levels(bwh$CLASS))]
   x<-as.matrix(exprs(bwh))
   nr<-nrow(x)
   nc<-ncol(x)
   rowInd <- nr:1
   colInd <- 1:nc
   
   #remove outliers
   outliers<-t(apply(x,1,quantile,c(0.01,0.99),na.rm=T))
   for (i in 1:nrow(x)){
      x[i,x[i,]>outliers[i,2]]<-outliers[i,2]
      x[i,x[i,]<outliers[i,1]]<-outliers[i,1]      
   }
   
   #row scaling
   rm <- rowMeans(x)
   x <- sweep(x, 1, rm)
   sx <- apply(x, 1, sd)
   x <- sweep(x, 1, sx, "/")
   
   #layout
   op <- par(no.readonly = TRUE)
   rightM<-5
   bottomM<-6
   lmat <- rbind(c(8,5), c(7,4), c(6,3), 2:1)
   layout(lmat, widths = c(0.1,0.9), heights = c(0.2,0.2,0.05,0.85), respect = FALSE)
   
   par(mar = c(bottomM, 1, 1, rightM))
   breaks=seq(min(x,na.rm=T), max(x,na.rm=T),length = 12)
   image(t(x[nr:1,nc:1]), 
         axes = F, 
         col = brewer.pal(11,"RdBu")[11:1],
         breaks=breaks)
   
   axis(1, 
        ((1:nc)-1)/(nc-1), 
        labels = colnames(x)[ncol(x):1], 
        las = 2,  
        tick = 0,
        cex.axis=0.7,
        line=-0.8)
   axis(4, 
        ((1:nr)-1)/(nr-1), 
        labels = rownames(x)[nrow(x):1], 
        las = 1,  
        tick = 0, 
        cex.axis=1,
        line=-0.8)
   
   par(mar = c(bottomM, 2, 1, 0))
   coef<-signature
   names(coef)<-NULL
   barplot.g(-abs(coef[nr:1]),horiz=T,axes=F,border=NA)
   mtext('Signature coefficients',side=2,cex=1.2)
   
   #truth
   par(mar = c(0, 1, 1, rightM))
   truth<-as.numeric(as.factor(bwh$CLASS))
   image(as.matrix(truth[nc:1]),
         axes=F,
         col=COLS)
   
   #BWH 
   cols<-as.numeric(bwh_labels)
   cols[is.na(cols)]<-3
   cols<-c('black','lightgrey','white')[cols]
   
   par(mar = c(0, 1, 2, rightM))
   barplot.g(bwh_scores[nc:1],
            axes=F,
            col=cols[nc:1],
            axisnames=F)
   mtext('Diagnostic score - BWH',side=3,cex=1)
   
   legend('right', 
          pch = 19, 
          levels(bwh$CLASS),
          col=COLS, 
          cex=.75)
   
   #UNMC
   cols<-as.numeric(unmc_labels)
   cols[is.na(cols)]<-3
   cols<-c('black','lightgrey','white')[cols]
   
   par(mar = c(0, 1, 2, rightM))
   barplot.g(unmc_scores[nc:1],
            axes=F,
            col=cols[nc:1],
            axisnames=F)
   mtext('Diagnostic score - UNMC',side=3,cex=1)
   
   par(op)  
}



##############################################################################
#using MYCKEY 1
##############################################################################
bwh<-readRDS('BWH_elements_M13_reference.RDS')
unmc<-readRDS('UNMC_elements_M13_reference.RDS')

#check which sames do not overlap
overlap<-intersect(bwh$id,unmc$id)
sort(setdiff(bwh$id,overlap))
sort(setdiff(unmc$id,overlap))

#reduce to overlapping samples
bwh<-bwh[,bwh$id%in%overlap]
unmc<-unmc[,unmc$id%in%overlap]

signature<-readRDS('diagnostic_coefficients_MYCKEY13.RDS')[-1]

#all samples
pdf('Diagnostic_classifier_all___M13_ref.pdf',width=14)
plotDiagnostic(bwh,unmc,signature,mode='13')
dev.off()


