rm(list = ls())

source("DESeq_Model.R")

##########################################################################################
gene_symbol<-mapIds(org.Mm.eg.db,keys=rownames(tab),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")

dir.create(path = "GeneExpressionPlots/", showWarnings = FALSE)

tiff("GeneExpressionPlots/Trem2.tiff", width = 1000, height = 1500,res = 150)

Trem2<-plotCounts(dm,gene=rownames(tab)[pmatch("Trem2",gene_symbol)], intgroup =c("Grouping","Genotype"),returnData=TRUE)
ggplot(Trem2, aes(x=Grouping, y=log10(count), color=Genotype)) +
  scale_y_continuous(name = "Counts", labels = scales::math_format(10^.x)) +
  #geom_text(aes(colour=Genotype),label=rownames(design))+
  annotation_logticks(sides = "l") +
  geom_point(position=position_jitter(width=.1,height=0), size=3)+
  ggtitle("Trem2 Expression")+
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

tiff("GeneExpressionPlots/Adgre1.tiff", width = 1000, height = 1500,res = 150)
Adgre1<-plotCounts(dm,gene=rownames(tab)[pmatch("Adgre1",gene_symbol)], intgroup =c("Grouping","Genotype"),returnData=TRUE)
ggplot(Adgre1, aes(x=Grouping, y=log10(count), color=Genotype)) +
  scale_y_continuous(name = "Counts", labels = scales::math_format(10^.x)) +
  #geom_text(aes(colour=Genotype),label=rownames(design))+
  annotation_logticks(sides = "l") +
  geom_point(position=position_jitter(width=.1,height=0), size=3)+
  ggtitle("Adgre1 Expression")+
  theme(plot.title = element_text(hjust = 0.5))

dev.off()


save.image(file = paste0("environments/GeneExpPlots_",Sys.Date(),".RData",sep=""))