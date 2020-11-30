rm(list = ls())

setwd("C:/Users/asbarros/Google Drive/RNASeq_icoelho/R_Analysis/")

source("scripts/DESeq_Model.R",echo=T)
require(VennDiagram)
require(gplots)

dev.off()

##########################################################################################
#Creating the subdirectories where to save all the results
dir.create(path = "KC_LDAM_Endo/", showWarnings = FALSE)
dir.create(path = "KC_LDAM_Endo/DE_Results_Filt/", showWarnings = FALSE)
dir.create(path = "KC_LDAM_Endo/Plots/", showWarnings = FALSE)

#subsetting the dm object to include only KC, Endo & LDAM
dm_sub<-dm[,dm$Grouping %in% c("KC","Endo_CD26+","LDAM")]
dm_sub$Grouping<-droplevels(dm_sub$Grouping)
dm_sub$Grouping<-relevel(dm_sub$Grouping,ref = "KC")

dvst_sub<-vst(dm_sub,blind=F)


pca1<-plotPCA(dvst_sub, intgroup = c("Grouping","Genotype"),returnData=T)
percentVar<-round(100*attr(pca1,"percentVar"))


data1<-ggplot(pca1,aes(PC1,PC2,colour=Grouping,shape=Genotype))+
  geom_point(aes(colour=Grouping),size=3)+
  #geom_text(aes(colour=Grouping),label=pca1$name)+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

tiff("KC_LDAM_Endo/Plots/PCA_KC_LDAM_ENDO.tiff", width = 1500, height = 1000,res = 150)
data1
dev.off()

#dm_sub<-dm
design(dm_sub)<- ~ Grouping
dm_sub<-DESeq(dm_sub)

resultsNames(dm_sub)

res1<-results(dm_sub,contrast = c("Grouping", "LDAM", "KC"),alpha = 0.05)
res2<-results(dm_sub,contrast = c("Grouping", "Endo_CD26.", "KC"),alpha = 0.05)
res3<-results(dm_sub,contrast = c("Grouping", "Endo_CD26.", "LDAM"),alpha = 0.05)

print(summary(res1))
print(summary(res2))
print(summary(res3))

res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de1<-lfcShrink(dm_sub,contrast = c("Grouping", "LDAM", "KC"),type="ashr",res=res1)
ashr_de2<-lfcShrink(dm_sub,contrast = c("Grouping", "Endo_CD26.", "KC"),type="ashr",res=res2)
ashr_de3<-lfcShrink(dm_sub,contrast = c("Grouping", "Endo_CD26.", "LDAM"),type="ashr",res=res3)

res1<-res1[,-4]
res2<-res2[,-4]
res3<-res3[,-4]

check1<-vector()
check2<-vector()
check3<-vector()

for (i in 1:ncol(ashr_de1)){
  check1[i]<-all(as.data.frame(res1)[,i]==as.data.frame(ashr_de1)[,i],na.rm=T)
  check2[i]<-all(as.data.frame(res2)[,i]==as.data.frame(ashr_de2)[,i],na.rm=T)
  check3[i]<-all(as.data.frame(res3)[,i]==as.data.frame(ashr_de3)[,i],na.rm=T)
}

print(check1)
print(check2)
print(check3)

de1<-data.frame(GeneSymbol=res1$symbol, GeneName=res1$geneName, 
                P_adj=res1$padj,log_ashr=ashr_de1$log2FoldChange)

de2<-data.frame(GeneSymbol=res1$symbol, GeneName=res1$geneName, 
                P_adj=res2$padj,log_ashr=ashr_de2$log2FoldChange)

de3<-data.frame(GeneSymbol=res1$symbol, GeneName=res1$geneName, 
                P_adj=res3$padj,log_ashr=ashr_de3$log2FoldChange)

de1<-de1[!is.na(de1$P_adj),]
de1<-de1[de1$P_adj < 0.05, ]

de2<-de2[!is.na(de2$P_adj),]
de2<-de2[de2$P_adj < 0.05, ]

de3<-de3[!is.na(de3$P_adj),]
de3<-de3[de3$P_adj < 0.05, ]

write.xlsx(de1,file = "KC_LDAM_Endo/DE_Results_Filt/LDAM_vs_KC_DiffExp_filt.xls",sheetName = "Summary")
write.xlsx(de2,file = "KC_LDAM_Endo/DE_Results_Filt/Endo_vs_KC_DiffExp_filt.xls",sheetName = "Summary")
write.xlsx(de3,file = "KC_LDAM_Endo/DE_Results_Filt/Endo_vs_LDAM_DiffExp_filt.xls",sheetName = "Summary")


de1_filt<-de1[de1$log_ashr >= 2, ]
de2_filt<-de2[de2$log_ashr >= 2, ]

n_join<-length(rownames(de1_filt)[rownames(de1_filt) %in% rownames(de2_filt)])
genes_join<-rownames(de1_filt)[rownames(de1_filt) %in% rownames(de2_filt)]

n_ldam<-length(rownames(de1_filt)[!(rownames(de1_filt) %in% rownames(de2_filt))])
genes_ldam<-rownames(de1_filt)[!(rownames(de1_filt) %in% rownames(de2_filt))]
                               
n_endo<-length(rownames(de2_filt)[!(rownames(de2_filt) %in% rownames(de1_filt))])
genes_endo<-rownames(de2_filt)[!(rownames(de2_filt) %in% rownames(de1_filt))]

tab_endo<-data.frame(Ensembl=genes_endo,
                          GeneSymbol=mapIds(org.Mm.eg.db,keys=genes_endo,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"))
tab_ldam<-data.frame(Ensembl=genes_ldam,
                     GeneSymbol=mapIds(org.Mm.eg.db,keys=genes_ldam,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"))

tab_join<-data.frame(Ensembl=genes_join,
                     GeneSymbol=mapIds(org.Mm.eg.db,keys=genes_join,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"))


write.xlsx(tab_join,file = "KC_LDAM_Endo/DE_Results_Filt/LDAM & Endo vs KC_list.xls",sheetName = "Join",append = F)
write.xlsx(tab_endo,file = "KC_LDAM_Endo/DE_Results_Filt/LDAM & Endo vs KC_list.xls",sheetName = "Endo",append=T)
write.xlsx(tab_ldam,file = "KC_LDAM_Endo/DE_Results_Filt/LDAM & Endo vs KC_list.xls",sheetName = "LDAM",append=T)


dev.off()

tiff("KC_LDAM_Endo/Plots/Venn_KC_LDAM_ENDO.tiff", width = 1500, height = 1000,res = 150)
draw.pairwise.venn(area1 = 1250, area2 = 546, cross.area = 429,
                   category = c("LDAM", "Endo/CD26+"), lty = "blank",
                 fill = c("#E69F00", "#56B4E9"))
dev.off()


sum_de1<-de1_filt[rownames(de1_filt) %in% genes_join,]
sum_de2<-de2_filt[rownames(de2_filt) %in% genes_join,]

all(rownames(sum_de1)==rownames(sum_de2))

sum<-data.frame(GeneSymbol=sum_de1$GeneSymbol,GeneName=sum_de1$GeneName,
                Padj_LDAM_vs_KC=sum_de1$P_adj, Log2FC_LDAM_vs_KC=sum_de1$log_ashr,
                Padj_Endo_vs_KC=sum_de2$P_adj, Log2FC_Endo_vs_KC=sum_de2$log_ashr)

rownames(sum)<-rownames(sum_de1)

write.xlsx(sum,file = "KC_LDAM_Endo/DE_Results_Filt/LDAM & Endo vs KC_join_filt.xls",sheetName = "Summary")

sep<-assay(dvst_sub)
sep_filt<-sep[rownames(sep) %in% genes_join,]

my_palette <- colorRampPalette(c("blue","white","red"))(n = 5000)


sep_filt<-sep_filt[,c("9-WT-KC-control","1-TREM2KO-KC-control",
                      "4-WT-LDAMs-APAP","WTAPAP1LDAMsB", "WTAPAP2LDAMsB",
                      "8-TREM2KO-LDAMs-APAP","Trem2KOAPAP1LDAMsA","Trem2KOAPAP2LDAMsA",
                      "WTcontrol1Endothelia-CD26posA","WTcontrol2Endothelia-CD26posA","Trem2KOcontrol1Endothelia-CD26posA",
                      "Trem2KOcontrol2Endothelia-CD26posA")]

tiff("KC_LDAM_Endo/Plots/heatmap_vst_nosamplecluster.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='row',Colv=NA)
dev.off()

tiff("KC_LDAM_Endo/Plots/heatmap_vst.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6))
dev.off()

out<-heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(5,6),dendrogram='row',Colv=NA)

genes_df<-data.frame(rownames(sep_filt),origin=c(rep("Overlap",429)))
write.xlsx(genes_df,file = "KC_LDAM_Endo/Plots/Genes_Order_Dendogram.xlsx",sheetName = "NoRowCluster",append = F)

sep_filt<-sep_filt[out$rowInd,]

genes_df<-genes_df[order(match(rownames(genes_df),out$rowInd)),]

write.xlsx(genes_df,file = "KC_LDAM_Endo/Plots/Genes_Order_Dendogram.xlsx",sheetName = "WithRowCluster",append = T)

sep<-assay(dm_sub)
sep_filt<-sep[rownames(sep) %in% genes_join,]

tiff("KC_LDAM_Endo/Plots/heatmap_raw_nosamplecluster.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='row',Colv=NA)
dev.off()

tiff("KC_LDAM_Endo/Plots/heatmap_raw.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='none',Colv=NA,Rowv = NA)
dev.off()


#######################################################
#Heatmap with the full dataset from KC vs LDAM vs Endo#
#######################################################
genes_list<-c(genes_ldam, genes_join,genes_endo)

sep<-assay(dvst_sub)
sep_filt<-sep[rownames(sep) %in% genes_list,]
sep_filt<-sep_filt[order(match(rownames(sep_filt),genes_list)),]
sep_filt<-sep_filt[,c("9-WT-KC-control","1-TREM2KO-KC-control",
                      "4-WT-LDAMs-APAP","WTAPAP1LDAMsB", "WTAPAP2LDAMsB",
                      "8-TREM2KO-LDAMs-APAP","Trem2KOAPAP1LDAMsA","Trem2KOAPAP2LDAMsA",
                      "WTcontrol1Endothelia-CD26posA","WTcontrol2Endothelia-CD26posA","Trem2KOcontrol1Endothelia-CD26posA",
                      "Trem2KOcontrol2Endothelia-CD26posA")]

all(rownames(sep_filt) == genes_list)

tiff("KC_LDAM_Endo/Plots/heatmap_fulldataset_col&rowcluster.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6))
dev.off()

tiff("KC_LDAM_Endo/Plots/heatmap_fulldataset_rowcluster.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='row',Colv=NA)
dev.off()

out<-heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='row',Colv = NA)

tiff("KC_LDAM_Endo/Plots/heatmap_fulldataset_nocluster.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram = 'none', Rowv=NA, Colv=NA )
dev.off()

genes_df<-data.frame(genes_list,origin=c(rep("LDAM Only", 821),rep("Overlap",429),rep("Endo Only",117)))

write.xlsx(genes_df,file = "KC_LDAM_Endo/Plots/GenesOrder_Dendogram_FullDataset.xlsx",sheetName = "NoRowCluster",append = F)

genes_df<-genes_df[order(match(rownames(genes_df),out$rowInd)),]

write.xlsx(genes_df,file = "KC_LDAM_Endo/Plots/GenesOrder_Dendogram_FullDataset.xlsx",sheetName = "WithRowCluster",append = T)

############################ Heatmaps with changed scales ##############################################
tiff("KC_LDAM_Endo/Plots/heatmap_fulldataset_col&rowcluster_changedscale.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",breaks=seq(0,17,length=50001),density.info='none',col=my_palette,trace="none",margins=c(20,6))
dev.off()

tiff("KC_LDAM_Endo/Plots/heatmap_fulldataset_rowcluster_changedscale.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",breaks=seq(0,17,length=50001),density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram='row',Colv=NA)
dev.off()


tiff("KC_LDAM_Endo/Plots/heatmap_fulldataset_nocluster_changedscale.tiff", width = 1200, height = 1200,res = 150)
heatmap.2(sep_filt,labRow = F,notecol="black",breaks=seq(0,17,length=50001),density.info='none',col=my_palette,trace="none",margins=c(20,6),dendrogram = 'none', Rowv=NA, Colv=NA )
dev.off()