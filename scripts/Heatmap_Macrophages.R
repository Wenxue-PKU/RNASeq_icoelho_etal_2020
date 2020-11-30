rm(list = ls())

require(VennDiagram)
library(xlsx)
library(gplots)
setwd("C:/Users/asbarros/Google Drive/RNASeq_icoelho/R_Analysis/")

##########################
#Running the DESeq model #
##########################

source("scripts/DESeq_Model.R",echo = T)

#######################
#Creating the folders #
#######################
dir.create(path = "Heatmap_Macrophages/", showWarnings = FALSE)
my_palette <- colorRampPalette(c("blue","white","red"))(n = 50000)

genes<-read.csv(file = "Heatmap_Macrophages/KC associated genes.csv",header=T,sep=";")
dm_sub<-dm[,!(dm$Grouping %in% c("Endo_CD26+","LDAM"))] #Excluding Endo & LDAM because they are endothelial, Excluding KC_APAP due to not enough sampling
dm_sub$Grouping<-droplevels(dm_sub$Grouping) #Dropping levels, in order for the model to work

tab<-colData(dm_sub)[order(colData(dm_sub)$Grouping),]
sample_order<-rownames(tab)

dvst_sub<-vst(dm_sub,blind=F)
sep<-assay(dvst_sub)

sep_filt<-sep[rownames(sep) %in% genes$Gene.stable.ID,]
genes<-genes[genes$Gene.stable.ID %in% rownames(sep_filt),]

sep_filt<-sep_filt[,sample_order]
all(colnames(sep_filt)==sample_order)

sep_filt<-sep_filt[order(match(rownames(sep_filt),genes$Gene.stable.ID)),]
all(rownames(sep_filt)==genes$Gene.stable.ID)

tiff("Heatmap_Macrophages/heatmap_KC_AssociatedGenes_col&rowcluster.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6))
dev.off()

tiff("Heatmap_Macrophages/heatmap_KC_AssociatedGenes_rowcluster.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='row',Colv = NA)
dev.off()

out<-heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='row',Colv = NA)

tiff("Heatmap_Macrophages/heatmap_KC_AssociatedGenes_nocluster.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram = 'none', Rowv=NA,
          Colv = NA)
dev.off()

genes_df<-data.frame(order=c(1:115),genes_list=genes$Gene.stable.ID)

write.xlsx(genes_df,file = "Heatmap_Macrophages/GenesOrder_Dendogram_KC_AssociatedGenes.xlsx",sheetName = "NoRowCluster",append = F)

genes_df<-genes_df[order(match(rownames(genes_df),out$rowInd)),]

write.xlsx(genes_df,file = "Heatmap_Macrophages/GenesOrder_Dendogram_KC_AssociatedGenes.xlsx",sheetName = "WithRowCluster",append = T)

########## Changed Scales ############
tiff("Heatmap_Macrophages/heatmap_KC_AssociatedGenes_col&rowcluster_changedscale.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",breaks=seq(0,20,length=50001),density.info='none',col=my_palette,trace="none",margins=c(20,6))
dev.off()

tiff("Heatmap_Macrophages/heatmap_KC_AssociatedGenes_rowcluster_changedscale.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",breaks=seq(0,20,length=50001),density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='row',Colv = NA)
dev.off()

out<-heatmap.2(sep_filt,labRow = F,notecol="black",breaks=seq(0,17,length=50001),density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='row',Colv = NA)
dev.off()

tiff("Heatmap_Macrophages/heatmap_KC_AssociatedGenes_nocluster_changedscale.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",breaks=seq(0,20,length=50001),density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram = 'none', Rowv=NA,
          Colv = NA)
dev.off()
