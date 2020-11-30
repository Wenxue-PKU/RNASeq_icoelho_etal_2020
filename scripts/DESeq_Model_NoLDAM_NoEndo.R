rm(list=ls())

library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("org.Mm.eg.db")
library("ggplot2")
library("genefilter")
library("xlsx")
library("stringr")
library("vsn")

theme_set(theme_bw())

######################################################################################
#setwd("//filessmb/folders/III/EXTRA/ONGOING/ASBARROS/RNASeq_icoelho/R_Analysis/")

source("scripts/DESeq_Model.R",echo = T)
dm_sub<-dm[,!(dm$Grouping %in% c("LDAM","Endo_CD26+"))]
dm_sub$Grouping<-droplevels(dm_sub$Grouping)

dvst_sub<-vst(dm_sub,blind=F)

ddist_sub<-dist(t(assay(dvst_sub)))

tiff("GeneralPlots/Heatmap_NoLDAM_NoEndo.tiff", width = 1500, height = 1000,res = 150)
heatmap1<-pheatmap(ddist_sub,cluster_rows=T,cluster_cols=T,annotation_col = design[,c("Batch","Grouping","Genotype")])
heatmap1
dev.off()

pca1<-plotPCA(dvst_sub, intgroup = c("Grouping","Genotype"),returnData=T)
percentVar<-round(100*attr(pca1,"percentVar"))

data1<-ggplot(pca1,aes(PC1,PC2,COlour=Grouping,shape=Genotype))+
  #geom_point(aes(colour=Grouping),size=3)+
  geom_text(aes(colour=Grouping),label=rownames(colData(dm_sub)))+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom",text=element_text(size=18))

tiff("GeneralPlots/PCA_NoLDAM_NoEndo_Names.tiff", width = 1500, height = 1000,res = 150)
data1
dev.off()

data2<-ggplot(pca1,aes(PC1,PC2,COlour=Grouping,shape=Genotype))+
  geom_point(aes(colour=Grouping),size=5)+
  #geom_text(aes(colour=Grouping),label=rownames(design))+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5,size=18),legend.position="bottom",text=element_text(size=16))

tiff("GeneralPlots/PCA_NoLDAM_NoEndo_Grouping&Genotype.tiff", width = 2300, height = 1500,res = 200)
data2
dev.off()

save.image(file = paste0("environments/DESeq_Model_NoLDAM_NoEndo",Sys.Date(),".RData",sep=""))
