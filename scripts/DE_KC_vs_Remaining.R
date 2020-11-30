rm(list = ls())

require(VennDiagram)
setwd("C:/Users/asbarros/Google Drive/RNASeq_icoelho/R_Analysis/")

##########################
#Running the DESeq model #
##########################

source("scripts/DESeq_Model.R",echo = T)

#######################
#Creating the folders #
#######################

dir.create(path = "KC_APAP_vs_Macrophages/", showWarnings = FALSE)
dir.create(path = "KC_APAP_vs_Macrophages/DE_Results_Filt/", showWarnings = FALSE)
dir.create(path = "KC_APAP_vs_Macrophages/Plots/", showWarnings = FALSE)

#############################################
#Subsetting the DESeq object for Macrophages#
#############################################

dm_sub<-dm[,!(dm$Grouping %in% c("Endo_CD26+","LDAM"))] #Excluding Endo & LDAM because they are endothelial, Excluding KC_APAP due to not enough sampling
dm_sub$Grouping<-droplevels(dm_sub$Grouping) #Dropping levels, in order for the model to work
dm_sub$Grouping<-relevel(dm_sub$Grouping,ref = "KC_APAP") #define KC_APAP as the baseline.

###########################################
#Subsetting the DESeq object for Wild Type#
###########################################
dm_wt<-dm_sub[,dm_sub$Genotype=="WT" | dm_sub$Grouping=="KC_APAP"] #Subset our subset, to inclue ONLY WT's AND the KC_APAP set
dm_wt$Genotype<-droplevels(dm_wt$Genotype) #Dropping levels in order for the model to work

dvst_sub<-vst(dm_wt,blind=F) #running the VST normalization in our subset

pca1<-plotPCA(dvst_sub, intgroup = c("Grouping","Genotype"),returnData=T) #getting the plot
percentVar<-round(100*attr(pca1,"percentVar")) #Extracting the percentage of variance explained 


data1<-ggplot(pca1,aes(PC1,PC2,colour=Grouping,shape=Genotype))+
  geom_point(aes(colour=Grouping),size=3)+
  #geom_text(aes(colour=Grouping),label=pca1$name)+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

tiff("KC_APAP_vs_Macrophages/Plots/PCA_WT.tiff", width = 1500, height = 1000,res = 150)
data1
dev.off()

design(dm_wt)<- ~ Grouping
dm_wt<-DESeq(dm_wt)

comp<-levels(dm_wt$Grouping)[-c(1,2)] #Remove the levels KC_APAP (baseline) & KC (one sample per genotype)

check<-vector()

#wb<-createWorkbook(type="xlsx")
sep_wt<-list()

for (i in comp){
  print(paste0("Performing DE for KC_APAP vs ", i, sep=""))
  
  res1<-results(dm_wt,contrast = c("Grouping", i, "KC_APAP"),alpha = 0.05)
  res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
  res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")
  print(summary(res1))
  
  ashr_de1<-lfcShrink(dm_wt,contrast = c("Grouping", i, "KC_APAP"),type="ashr",res=res1)
  res1_sub<-res1[,-4]
  for (j in 1:ncol(ashr_de1)) {check[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)}
  
  print(check)
  
  de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                  PAdj=res1$padj,Log2FC_ashr=ashr_de1$log2FoldChange)
  
  de1_filt<-de1[!is.na(de1$PAdj),]
  de1_filt<-de1_filt[de1_filt$PAdj < 0.05,]
  
  j<-which(comp==i)
  sep_wt[[j]]<-de1_filt
  names(sep_wt)[j]<-i
  
  #sheet<-createSheet(wb, sheetName = paste0(i,"_vs_KC_APAP_DE",sep=""))
  #addDataFrame(de1_filt, sheet, startRow=1, startColumn=1)
  
}

#saveWorkbook(wb, "KC_APAP_vs_Macrophages/DE_Results_Filt/WT_KC_APAPvsAll_DE_filt.xlsx")

###########################################
#Subsetting the DESeq object for KO#
###########################################
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

tiff("KC_APAP_vs_Macrophages/Plots/PCA_Trem2KO.tiff", width = 1500, height = 1000,res = 150)
data1
dev.off()

design(dm_ko)<- ~ Grouping
dm_ko<-DESeq(dm_ko)

comp<-levels(dm_wt$Grouping)[-c(1,2)] #Remove the levels KC_APAP (baseline) & KC (one sample per genotype)

check<-vector()

#wb<-createWorkbook(type="xlsx")
sep_ko<-list()

for (i in comp){
  print(paste0("Performing DE for KC_APAP vs ", i, sep=""))
  
  res1<-results(dm_ko,contrast = c("Grouping", i, "KC_APAP"),alpha = 0.05)
  res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
  res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")
  print(summary(res1))
  
  ashr_de1<-lfcShrink(dm_ko,contrast = c("Grouping", i, "KC_APAP"),type="ashr",res=res1)
  res1_sub<-res1[,-4]
  for (j in 1:ncol(ashr_de1)) {check[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)}
  
  print(check)
  
  de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                  PAdj=res1$padj,Log2FC_ashr=ashr_de1$log2FoldChange)
  
  de1_filt<-de1[!is.na(de1$PAdj),]
  de1_filt<-de1_filt[de1_filt$PAdj < 0.05,]
  
  j<-which(comp==i)
  sep_ko[[j]]<-de1_filt
  names(sep_ko)[j]<-i
  
  #sheet<-createSheet(wb, sheetName = paste0(i,"_vs_KC_APAP_DE",sep=""))
  #addDataFrame(de1_filt, sheet, startRow=1, startColumn=1)
  
}

#saveWorkbook(wb, "KC_APAP_vs_Macrophages/DE_Results_Filt/Trem2KO_KC_APAPvsAll_DE_filt.xlsx")



##################################################
#Assessing the overlap & Performing Venn Diagrams#
##################################################
dir.create(path = "KC_APAP_vs_Macrophages/Overlapping_Results", showWarnings = FALSE)
dir.create(path = "KC_APAP_vs_Macrophages/Venn", showWarnings = FALSE)

for (i in 1:length(sep_wt)) {
  print(paste0("Performing Overlapp Analysis for ", names(sep_wt)[i], sep=""))
  
  wb<-createWorkbook(type="xlsx")
  
  join<-rownames(sep_wt[[i]])[rownames(sep_wt[[i]]) %in% rownames(sep_ko[[i]])]
  wt<-rownames(sep_wt[[i]])[!(rownames(sep_wt[[i]]) %in% rownames(sep_ko[[i]]))]
  ko<-rownames(sep_ko[[i]])[!(rownames(sep_ko[[i]]) %in% rownames(sep_wt[[i]]))]
  
  print(length(join))
  print(length(wt))
  print(length(ko))
  
  join_df<-data.frame(#Genes=join,
                      GeneSymbol=mapIds(org.Mm.eg.db,keys=join,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
                      GeneName = mapIds(org.Mm.eg.db,keys=join,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  wt_df<-data.frame(#Genes=wt,
                      GeneSymbol=mapIds(org.Mm.eg.db,keys=wt,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
                      GeneName = mapIds(org.Mm.eg.db,keys=wt,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  
  ko_df<-data.frame(#Genes=ko,
                      GeneSymbol=mapIds(org.Mm.eg.db,keys=ko,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"),
                      GeneName = mapIds(org.Mm.eg.db,keys=ko,column="GENENAME", keytype="ENSEMBL", multiVals="first"))
  
  sheet<-createSheet(wb, sheetName = "Overlapping")
  addDataFrame(join_df, sheet, startRow=1, startColumn=1)
  sheet<-createSheet(wb, sheetName = "WT_only")
  addDataFrame(wt_df, sheet, startRow=1, startColumn=1)
  sheet<-createSheet(wb, sheetName = "KO_only")
  addDataFrame(ko_df, sheet, startRow=1, startColumn=1)
  
  saveWorkbook(wb, paste0("KC_APAP_vs_Macrophages/Overlapping_Results/",names(sep_wt)[i],"_Overlapp.xlsx",sep=""))
  
  tiff(paste0("KC_APAP_vs_Macrophages/Venn/VD_", names(sep_wt)[i],".tiff",sep=""), width = 1500, height = 1000,res = 150)
  
  draw.pairwise.venn(area1 = length(wt) + length(join), 
                     area2 = length(ko) + length(join),
                     cross.area = length(join),
                     category = c("WT", "KO"), lty = "blank",
                     fill = c("#E69F00", "#56B4E9"))
  dev.off()
}









