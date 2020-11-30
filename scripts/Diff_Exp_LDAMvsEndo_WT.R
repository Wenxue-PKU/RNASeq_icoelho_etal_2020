rm(list = ls())

source("scripts/DESeq_Model.R")

##########################################################################################
dm_sub<-dm[,dm$Genotype == "WT"]
dm_sub$Genotype<-droplevels(dm_sub$Genotype)

dm_sub$Grouping<-relevel(x = dm_sub$Grouping,ref = "LDAM")

design(dm_sub)<- ~ Grouping
dm_sub<-DESeq(dm_sub)

resultsNames(dm_sub)
check<-vector()

res1<-results(dm_sub,name = "Grouping_Endo_CD26._vs_LDAM",alpha = 0.05)
res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

print(summary(res1))
ashr_de1<-lfcShrink(dm_sub,coef =  "Grouping_Endo_CD26._vs_LDAM",type="ashr",res=res1)
res1_sub<-res1[,-4]
for (j in 1:ncol(ashr_de1)) {check[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)}

print(check)

apeglm_de1<-lfcShrink(dm_sub,coef = "Grouping_Endo_CD26._vs_LDAM",type="apeglm",res = res1)
for (j in 1:ncol(apeglm_de1)) {check[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(apeglm_de1)[,j],na.rm=T)}

print(check)


de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                baseMean=res1$baseMean,Pval_Adj=res1$padj,
                NoCorr_logFC=res1$log2FoldChange,ashr_logFC=ashr_de1$log2FoldChange,
                apeglm_logFC=apeglm_de1$log2FoldChange)

de1<-de1[!is.na(de1$Pval_Adj),]
de1<-de1[de1$Pval_Adj < 0.05,]

sum_de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,Pval_Adj=res1$padj,logFC=apeglm_de1$log2FoldChange)
sum_de1<-sum_de1[!is.na(sum_de1$Pval_Adj), ]
sum_de1<-sum_de1[sum_de1$Pval_Adj< 0.05,]

#write.xlsx(x = sum_de1,file = paste0("DiffExp_Results_filt/","Endo_CD26+_vs_LDAM","_DiffExp.xls",sep=""),sheetName = "Summary", append = F)
#write.xlsx(x = de1,file = paste0("DiffExp_Results_filt/","Endo_CD26+_vs_LDAM","_DiffExp.xls",sep=""),sheetName = "Extended", append = T)