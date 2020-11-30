rm(list = ls())

require(VennDiagram)
library(xlsx)
setwd("C:/Users/asbarros/Google Drive/RNASeq_icoelho/R_Analysis/")

##########################
#Running the DESeq model #
##########################

source("scripts/DESeq_Model.R",echo = T)

#######################
#Creating the folders #
#######################
dir.create(path = "DE_LDAM&Endo&KC_APAP/", showWarnings = FALSE)
dir.create(path = "DE_LDAM&Endo&KC_APAP/Plots/", showWarnings = FALSE)
dir.create(path = "DE_LDAM&Endo&KC_APAP/Diff_Exp/", showWarnings = FALSE)
dir.create(path = "DE_LDAM&Endo&KC_APAP/Diff_Exp/DE_Results", showWarnings = FALSE)
dir.create(path = "DE_LDAM&Endo&KC_APAP/Diff_Exp/VD", showWarnings = FALSE)

dir.create(path = "DE_LDAM&Endo&KC_APAP/Up_Regulated/", showWarnings = FALSE)
dir.create(path = "DE_LDAM&Endo&KC_APAP/Up_Regulated/DE_Results", showWarnings = FALSE)
dir.create(path = "DE_LDAM&Endo&KC_APAP/Up_Regulated/VD", showWarnings = FALSE)

dir.create(path = "DE_LDAM&Endo&KC_APAP/Down_Regulated/", showWarnings = FALSE)
dir.create(path = "DE_LDAM&Endo&KC_APAP/Down_Regulated/DE_Results", showWarnings = FALSE)
dir.create(path = "DE_LDAM&Endo&KC_APAP/Down_Regulated/VD", showWarnings = FALSE)

#############################################
#Subsetting the DESeq object for Macrophages#
#############################################

dm_sub<-dm[,dm$Grouping %in% c("Endo_CD26+","LDAM","KC_APAP")] #Excluding Endo & LDAM because they are endothelial, Excluding KC_APAP due to not enough sampling
dm_sub$Grouping<-droplevels(dm_sub$Grouping) #Dropping levels, in order for the model to work
dm_sub$Grouping<-relevel(dm_sub$Grouping,ref = "KC_APAP") #define KC_APAP as the baseline.

###########################################
#Subsetting the DESeq object for Wild Type#
###########################################
dm_wt<-dm_sub[,dm_sub$Genotype=="WT" | dm_sub$Grouping=="KC_APAP"] #Subset our subset, to inclue ONLY WT's AND the KC_APAP set
dm_wt$Genotype<-droplevels(dm_wt$Genotype) #Dropping levels in order for the model to work

dvst_sub<-vst(dm_wt,blind=F) #running the VST normalization in our subset
data1<-ggplot(pca1,aes(PC1,PC2,colour=Grouping,shape=Genotype))+
  geom_point(aes(colour=Grouping),size=3)+
  #geom_text(aes(colour=Grouping),label=pca1$name)+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

tiff("DE_LDAM&Endo&KC_APAP/Plots/PCA_WT.tiff", width = 1500, height = 1000,res = 150)
data1
dev.off()

design(dm_wt)<- ~ Grouping
dm_wt<-DESeq(dm_wt)

#################### Diff. Expressed #######################

comp<-levels(dm_wt$Grouping)
check1<-vector()
check2<-vector()
check3<-vector()

res1<-results(dm_wt,contrast = c("Grouping", comp[2],comp[1]),alpha = 0.05) #Endo vs KC_APAP
res2<-results(dm_wt,contrast = c("Grouping", comp[3],comp[1]),alpha = 0.05) #LDAM vs KC_APAP
res3<-results(dm_wt,contrast = c("Grouping", comp[3],comp[2]),alpha = 0.05) #LDAM vs Endo

res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de1<-lfcShrink(dm_wt,contrast = c("Grouping", comp[2],comp[1]),type="ashr",res=res1)
ashr_de2<-lfcShrink(dm_wt,contrast = c("Grouping", comp[3],comp[1]),type="ashr",res=res2)
ashr_de3<-lfcShrink(dm_wt,contrast = c("Grouping", comp[3],comp[2]),type="ashr",res=res3)

res1_sub<-res1[,-4]
res2_sub<-res2[,-4]
res3_sub<-res3[,-4]

for (j in 1:ncol(ashr_de1)) {
  check1[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)
  check2[j]<-all(as.data.frame(res2_sub)[,j] == as.data.frame(ashr_de2)[,j],na.rm=T)
  check3[j]<-all(as.data.frame(res3_sub)[,j] == as.data.frame(ashr_de3)[,j],na.rm=T)
}

print(check1)
print(check2)
print(check3)

de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                PAdj=res1$padj,Log2FC_ashr=ashr_de1$log2FoldChange)

de2<-data.frame(Gene_Symbol = res2$symbol, Gene_Name=res2$geneName,
                PAdj=res2$padj,Log2FC_ashr=ashr_de2$log2FoldChange)

de3<-data.frame(Gene_Symbol = res3$symbol, Gene_Name=res3$geneName,
                PAdj=res3$padj,Log2FC_ashr=ashr_de3$log2FoldChange)

de1_filt<-de1[!is.na(de1$PAdj),]
de1_filt<-de1_filt[de1_filt$PAdj < 0.05,]

de2_filt<-de2[!is.na(de2$PAdj),]
de2_filt<-de2_filt[de2_filt$PAdj < 0.05,]

de3_filt<-de3[!is.na(de3$PAdj),]
de3_filt<-de3_filt[de3_filt$PAdj < 0.05,]

write.xlsx(de1_filt,file = "DE_LDAM&Endo&KC_APAP/Diff_Exp/DE_Results/WT_DE_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[2],"_vs_",comp[1],"_filt"),append = F)
write.xlsx(de2_filt,file = "DE_LDAM&Endo&KC_APAP/Diff_Exp/DE_Results/WT_DE_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[1],"_filt"),append = T)
write.xlsx(de3_filt,file = "DE_LDAM&Endo&KC_APAP/Diff_Exp/DE_Results/WT_DE_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[2],"_filt"),append = T)

wt_de_list<-list(de1_filt,de2_filt,de3_filt)
names(wt_de_list)<-c(paste0(comp[2],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[2]))

#################### Up-Regulated #######################

up1_filt<-de1_filt[de1_filt$Log2FC_ashr > 0, ]
up2_filt<-de2_filt[de2_filt$Log2FC_ashr > 0, ]
up3_filt<-de3_filt[de3_filt$Log2FC_ashr > 0, ]

write.xlsx(up1_filt,file = "DE_LDAM&Endo&KC_APAP/Up_Regulated/DE_Results/WT_Up_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[2],"_vs_",comp[1],"_up_filt"),append = F)
write.xlsx(up2_filt,file = "DE_LDAM&Endo&KC_APAP/Up_Regulated/DE_Results/WT_Up_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[1],"_up_filt"),append = T)
write.xlsx(up3_filt,file = "DE_LDAM&Endo&KC_APAP/Up_Regulated/DE_Results/WT_Up_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[2],"_up_filt"),append = T)

wt_up_list<-list(up1_filt,up2_filt,up3_filt)
names(wt_up_list)<-c(paste0(comp[2],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[2]))

#################### Down-Regulated #######################

down1_filt<-de1_filt[de1_filt$Log2FC_ashr < 0, ]
down2_filt<-de2_filt[de2_filt$Log2FC_ashr < 0, ]
down3_filt<-de3_filt[de3_filt$Log2FC_ashr < 0, ]

write.xlsx(down1_filt,file = "DE_LDAM&Endo&KC_APAP/Down_Regulated/DE_Results/WT_Down_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[2],"_vs_",comp[1],"_down_filt"),append = F)
write.xlsx(down2_filt,file = "DE_LDAM&Endo&KC_APAP/Down_Regulated/DE_Results/WT_Down_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[1],"_down_filt"),append = T)
write.xlsx(down3_filt,file = "DE_LDAM&Endo&KC_APAP/Down_Regulated/DE_Results/WT_Down_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[2],"_down_filt"),append = T)

wt_down_list<-list(down1_filt,down2_filt,down3_filt)
names(wt_down_list)<-c(paste0(comp[2],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[2]))

####################################
#Subsetting the DESeq object for KO#
####################################
dm_ko<-dm_sub[,!(dm_sub$Genotype=="WT") | dm_sub$Grouping=="KC_APAP"]
dm_ko$Grouping<-droplevels(dm_ko$Grouping)
dm_ko$Genotype<-droplevels(dm_ko$Genotype)

dvst_sub<-vst(dm_ko,blind=F)

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

tiff("DE_LDAM&Endo&KC_APAP/Plots/PCA_Trem2KO.tiff", width = 1500, height = 1000,res = 150)
data1
dev.off()

design(dm_ko)<- ~ Grouping
dm_ko<-DESeq(dm_ko)


#################### Diff. Expressed #######################

comp<-levels(dm_ko$Grouping) #Remove the levels KC_APAP (baseline) & KC (one sample per genotype)
check1<-vector()
check2<-vector()
check3<-vector()

res1<-results(dm_ko,contrast = c("Grouping", comp[2],comp[1]),alpha = 0.05)
res2<-results(dm_ko,contrast = c("Grouping", comp[3],comp[1]),alpha = 0.05)
res3<-results(dm_ko,contrast = c("Grouping", comp[3],comp[2]),alpha = 0.05)

res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")

res3$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res3),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res3$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res3), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de1<-lfcShrink(dm_ko,contrast = c("Grouping", comp[2],comp[1]),type="ashr",res=res1)
ashr_de2<-lfcShrink(dm_ko,contrast = c("Grouping", comp[3],comp[1]),type="ashr",res=res2)
ashr_de3<-lfcShrink(dm_ko,contrast = c("Grouping", comp[3],comp[2]),type="ashr",res=res3)

res1_sub<-res1[,-4]
res2_sub<-res2[,-4]
res3_sub<-res3[,-4]

for (j in 1:ncol(ashr_de1)) {
  check1[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)
  check2[j]<-all(as.data.frame(res2_sub)[,j] == as.data.frame(ashr_de2)[,j],na.rm=T)
  check3[j]<-all(as.data.frame(res3_sub)[,j] == as.data.frame(ashr_de3)[,j],na.rm=T)
}

print(check1)
print(check2)
print(check3)

de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                PAdj=res1$padj,Log2FC_ashr=ashr_de1$log2FoldChange)

de2<-data.frame(Gene_Symbol = res2$symbol, Gene_Name=res2$geneName,
                PAdj=res2$padj,Log2FC_ashr=ashr_de2$log2FoldChange)

de3<-data.frame(Gene_Symbol = res3$symbol, Gene_Name=res3$geneName,
                PAdj=res3$padj,Log2FC_ashr=ashr_de3$log2FoldChange)

de1_filt<-de1[!is.na(de1$PAdj),]
de1_filt<-de1_filt[de1_filt$PAdj < 0.05,]

de2_filt<-de2[!is.na(de2$PAdj),]
de2_filt<-de2_filt[de2_filt$PAdj < 0.05,]

de3_filt<-de3[!is.na(de3$PAdj),]
de3_filt<-de3_filt[de3_filt$PAdj < 0.05,]

write.xlsx(de1_filt,file = "DE_LDAM&Endo&KC_APAP/Diff_Exp/DE_Results/KO_DE_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[2],"_vs_",comp[1],"_filt"),append = F)
write.xlsx(de2_filt,file = "DE_LDAM&Endo&KC_APAP/Diff_Exp/DE_Results/KO_DE_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[1],"_filt"),append = T)
write.xlsx(de3_filt,file = "DE_LDAM&Endo&KC_APAP/Diff_Exp/DE_Results/KO_DE_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[2],"_filt"),append = T)

ko_de_list<-list(de1_filt,de2_filt,de3_filt)
names(ko_de_list)<-c(paste0(comp[2],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[2]))

#################### Up-Regulated #######################

up1_filt<-de1_filt[de1_filt$Log2FC_ashr > 0, ]
up2_filt<-de2_filt[de2_filt$Log2FC_ashr > 0, ]
up3_filt<-de3_filt[de3_filt$Log2FC_ashr > 0, ]

write.xlsx(up1_filt,file = "DE_LDAM&Endo&KC_APAP/Up_Regulated/DE_Results/KO_Up_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[2],"_vs_",comp[1],"_up_filt"),append = F)
write.xlsx(up2_filt,file = "DE_LDAM&Endo&KC_APAP/Up_Regulated/DE_Results/KO_Up_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[1],"_up_filt"),append = T)
write.xlsx(up3_filt,file = "DE_LDAM&Endo&KC_APAP/Up_Regulated/DE_Results/KO_Up_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[2],"_up_filt"),append = T)

ko_up_list<-list(up1_filt,up2_filt,up3_filt)
names(ko_up_list)<-c(paste0(comp[2],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[2]))

#################### Down-Regulated #######################

down1_filt<-de1_filt[de1_filt$Log2FC_ashr < 0, ]
down2_filt<-de2_filt[de2_filt$Log2FC_ashr < 0, ]
down3_filt<-de3_filt[de3_filt$Log2FC_ashr < 0, ]

write.xlsx(down1_filt,file = "DE_LDAM&Endo&KC_APAP/Down_Regulated/DE_Results/KO_Down_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[2],"_vs_",comp[1],"_down_filt"),append = F)
write.xlsx(down2_filt,file = "DE_LDAM&Endo&KC_APAP/Down_Regulated/DE_Results/KO_Down_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[1],"_down_filt"),append = T)
write.xlsx(down3_filt,file = "DE_LDAM&Endo&KC_APAP/Down_Regulated/DE_Results/KO_Down_LDAM&Endo&KC_APAP_Results.xlsx",
           sheetName = paste0(comp[3],"_vs_",comp[2],"_down_filt"),append = T)

ko_down_list<-list(down1_filt,down2_filt,down3_filt)
names(ko_down_list)<-c(paste0(comp[2],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[1]),paste0(comp[3],"_vs_",comp[2]))

##################################################
#Assessing the overlap & Performing Venn Diagrams#
##################################################

#### Diff Expressed ####

for (i in 1:length(wt_de_list)) {
  join<-rownames(wt_de_list[[i]])[rownames(wt_de_list[[i]]) %in% rownames(ko_de_list[[i]])]
  wt<-rownames(wt_de_list[[i]])[!(rownames(wt_de_list[[i]]) %in% rownames(ko_de_list[[i]]))]
  ko<-rownames(ko_de_list[[i]])[!(rownames(ko_de_list[[i]]) %in% rownames(wt_de_list[[i]]))]
  
  join_df<-data.frame(#Genes=join,
    GeneSymbol=mapIds(org.Mm.eg.db,keys=join,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
    GeneName = mapIds(org.Mm.eg.db,keys=join,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  wt_df<-data.frame(#Genes=wt,
    GeneSymbol=mapIds(org.Mm.eg.db,keys=wt,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
    GeneName = mapIds(org.Mm.eg.db,keys=wt,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  
  ko_df<-data.frame(#Genes=ko,
    GeneSymbol=mapIds(org.Mm.eg.db,keys=ko,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
    GeneName = mapIds(org.Mm.eg.db,keys=ko,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  
  write.xlsx(join_df,file = paste0("DE_LDAM&Endo&KC_APAP/Diff_Exp/VD/",names(wt_de_list)[i],
                                   "_Genes_VD.xlsx"),sheetName = "Genes_Overlap",append = F)
  write.xlsx(wt_df,file = paste0("DE_LDAM&Endo&KC_APAP/Diff_Exp/VD/",names(wt_de_list)[i],
                                 "_Genes_VD.xlsx"),sheetName = "Genes_WT_Only",append = T)
  write.xlsx(ko_df,file = paste0("DE_LDAM&Endo&KC_APAP/Diff_Exp/VD/",names(wt_de_list)[i],
                                 "_Genes_VD.xlsx"),sheetName = "Genes_KO_Only",append = T)
  
  
  tiff(paste0("DE_LDAM&Endo&KC_APAP/Diff_Exp/VD/", names(wt_de_list)[i],"_DE.tiff",sep=""), width = 1500, height = 1000,res = 150)
  
  draw.pairwise.venn(area1 = length(wt) + length(join), 
                     area2 = length(ko) + length(join),
                     cross.area = length(join),
                     category = c("WT", "KO"), lty = "blank",
                     fill = c("#E69F00", "#56B4E9"))
  dev.off()
  
}

for (i in 1:length(wt_up_list)) {
  join<-rownames(wt_up_list[[i]])[rownames(wt_up_list[[i]]) %in% rownames(ko_up_list[[i]])]
  wt<-rownames(wt_up_list[[i]])[!(rownames(wt_up_list[[i]]) %in% rownames(ko_up_list[[i]]))]
  ko<-rownames(ko_up_list[[i]])[!(rownames(ko_up_list[[i]]) %in% rownames(wt_up_list[[i]]))]
  
  join_df<-data.frame(#Genes=join,
    GeneSymbol=mapIds(org.Mm.eg.db,keys=join,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
    GeneName = mapIds(org.Mm.eg.db,keys=join,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  wt_df<-data.frame(#Genes=wt,
    GeneSymbol=mapIds(org.Mm.eg.db,keys=wt,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
    GeneName = mapIds(org.Mm.eg.db,keys=wt,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  
  ko_df<-data.frame(#Genes=ko,
    GeneSymbol=mapIds(org.Mm.eg.db,keys=ko,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
    GeneName = mapIds(org.Mm.eg.db,keys=ko,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  
  write.xlsx(join_df,file = paste0("DE_LDAM&Endo&KC_APAP/Up_Regulated/VD/",names(wt_de_list)[i],
                                   "_Genes_VD.xlsx"),sheetName = "Genes_Overlap",append = F)
  write.xlsx(wt_df,file = paste0("DE_LDAM&Endo&KC_APAP/Up_Regulated/VD/",names(wt_de_list)[i],
                                 "_Genes_VD.xlsx"),sheetName = "Genes_WT_Only",append = T)
  write.xlsx(ko_df,file = paste0("DE_LDAM&Endo&KC_APAP/Up_Regulated/VD/",names(wt_de_list)[i],
                                 "_Genes_VD.xlsx"),sheetName = "Genes_KO_Only",append = T)
  
  
  tiff(paste0("DE_LDAM&Endo&KC_APAP/Up_Regulated/VD/", names(wt_de_list)[i],"_Up.tiff",sep=""), width = 1500, height = 1000,res = 150)
  
  draw.pairwise.venn(area1 = length(wt) + length(join), 
                     area2 = length(ko) + length(join),
                     cross.area = length(join),
                     category = c("WT", "KO"), lty = "blank",
                     fill = c("#E69F00", "#56B4E9"))
  dev.off()
  
}

for (i in 1:length(wt_down_list)) {
  join<-rownames(wt_down_list[[i]])[rownames(wt_down_list[[i]]) %in% rownames(ko_down_list[[i]])]
  wt<-rownames(wt_down_list[[i]])[!(rownames(wt_down_list[[i]]) %in% rownames(ko_down_list[[i]]))]
  ko<-rownames(ko_down_list[[i]])[!(rownames(ko_down_list[[i]]) %in% rownames(wt_down_list[[i]]))]
  
  join_df<-data.frame(#Genes=join,
    GeneSymbol=mapIds(org.Mm.eg.db,keys=join,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
    GeneName = mapIds(org.Mm.eg.db,keys=join,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  wt_df<-data.frame(#Genes=wt,
    GeneSymbol=mapIds(org.Mm.eg.db,keys=wt,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
    GeneName = mapIds(org.Mm.eg.db,keys=wt,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  
  ko_df<-data.frame(#Genes=ko,
    GeneSymbol=mapIds(org.Mm.eg.db,keys=ko,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
    GeneName = mapIds(org.Mm.eg.db,keys=ko,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  
  write.xlsx(join_df,file = paste0("DE_LDAM&Endo&KC_APAP/Down_Regulated/VD/",names(wt_de_list)[i],
                                   "_Genes_VD.xlsx"),sheetName = "Genes_Overlap",append = F)
  write.xlsx(wt_df,file = paste0("DE_LDAM&Endo&KC_APAP/Down_Regulated/VD/",names(wt_de_list)[i],
                                 "_Genes_VD.xlsx"),sheetName = "Genes_WT_Only",append = T)
  write.xlsx(ko_df,file = paste0("DE_LDAM&Endo&KC_APAP/Down_Regulated/VD/",names(wt_de_list)[i],
                                 "_Genes_VD.xlsx"),sheetName = "Genes_KO_Only",append = T)
  
  
  tiff(paste0("DE_LDAM&Endo&KC_APAP/Down_Regulated/VD/", names(wt_de_list)[i],"_Down.tiff",sep=""), width = 1500, height = 1000,res = 150)
  
  draw.pairwise.venn(area1 = length(wt) + length(join), 
                     area2 = length(ko) + length(join),
                     cross.area = length(join),
                     category = c("WT", "KO"), lty = "blank",
                     fill = c("#E69F00", "#56B4E9"))
  dev.off()
  
}
