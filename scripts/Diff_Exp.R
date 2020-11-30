rm(list = ls())

setwd("C:/Users/asbarros/Google Drive/RNASeq_icoelho/R_Analysis/")

source("scripts/DESeq_Model.R")

##########################################################################################
dir.create(path = "DiffExp_Results/", showWarnings = FALSE)

types<-levels(dm$Grouping)
types<-types[-c(2,3)] #Remove KC_APAP, since it only has WT's and KC, because only have one sample per group
check<-vector()

for (i in types){
  print(paste0("Performing DE for ", i, sep=""))
  dm_sub<-dm[,dm$Grouping %in% c(i)]
  dm_sub$Grouping<-droplevels(dm_sub$Grouping)
  
  design(dm_sub)<- ~ Genotype
  dm_sub<-DESeq(dm_sub)
  
  res1<-results(dm_sub,name = "Genotype_Trem2_KO_vs_WT",alpha = 0.05)
  res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
  res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")
  
  print(summary(res1))
  ashr_de1<-lfcShrink(dm_sub,coef =  "Genotype_Trem2_KO_vs_WT",type="ashr",res=res1)
  res1_sub<-res1[,-4]
  for (j in 1:ncol(ashr_de1)) {check[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)}
  
  check
  
  apeglm_de1<-lfcShrink(dm_sub,coef = "Genotype_Trem2_KO_vs_WT",type="apeglm",res = res1)
  for (j in 1:ncol(apeglm_de1)) {check[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(apeglm_de1)[,j],na.rm=T)}
  
  de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                  baseMean=res1$baseMean,Pval_Adj=res1$padj,
                  NoCorr_logFC=res1$log2FoldChange,ashr_logFC=ashr_de1$log2FoldChange,
                  apeglm_logFC=apeglm_de1$log2FoldChange)
  
  sum_de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,Pval_Adj=res1$padj,logFC=apeglm_de1$log2FoldChange)
  
  all(de1$apeglm_de1$log2FoldChange == sum_de1$logFC)
  
  #write.xlsx(x = sum_de1,file = paste0("DiffExp_Results/",i,"_DiffExp.xls",sep=""),sheetName = "Summary", append = F)
  #write.xlsx(x = de1,file = paste0("DiffExp_Results/",i,"_DiffExp.xls",sep=""),sheetName = "Extended", append = T)
  }

############ Filtering ###########################
dir.create(path = "DiffExp_Results_filt/", showWarnings = FALSE)

for (i in types){
  print(paste0("Performing DE for ", i, sep=""))
  dm_sub<-dm[,dm$Grouping %in% c(i)]
  dm_sub$Grouping<-droplevels(dm_sub$Grouping)
  
  design(dm_sub)<- ~ Genotype
  dm_sub<-DESeq(dm_sub)
  
  res1<-results(dm_sub,name = "Genotype_Trem2_KO_vs_WT",alpha = 0.05)
  res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
  res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")
  
  print(summary(res1))
  ashr_de1<-lfcShrink(dm_sub,coef =  "Genotype_Trem2_KO_vs_WT",type="ashr",res=res1)
  res1_sub<-res1[,-4]
  for (j in 1:ncol(ashr_de1)) {check[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)}
  
  check
  
  apeglm_de1<-lfcShrink(dm_sub,coef = "Genotype_Trem2_KO_vs_WT",type="apeglm",res = res1)
  for (j in 1:ncol(apeglm_de1)) {check[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(apeglm_de1)[,j],na.rm=T)}
  
  de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                  baseMean=res1$baseMean,Pval_Adj=res1$padj,
                  NoCorr_logFC=res1$log2FoldChange,ashr_logFC=ashr_de1$log2FoldChange,
                  apeglm_logFC=apeglm_de1$log2FoldChange)
  
  de1<-de1[!is.na(de1$Pval_Adj),]
  de1<-de1[de1$Pval_Adj < 0.05,]
  
  sum_de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,Pval_Adj=res1$padj,logFC=apeglm_de1$log2FoldChange)
  sum_de1<-sum_de1[!is.na(sum_de1$Pval_Adj), ]
  sum_de1<-sum_de1[sum_de1$Pval_Adj< 0.05,]
  
  all(de1$apeglm_de1$log2FoldChange == sum_de1$logFC)
  
  #write.xlsx(x = sum_de1,file = paste0("DiffExp_Results_filt/",i,"_DiffExp_filt.xls",sep=""),sheetName = "Summary", append = F)
  #write.xlsx(x = de1,file = paste0("DiffExp_Results_filt/",i,"_DiffExp_filt.xls",sep=""),sheetName = "Extended", append = T)
}
