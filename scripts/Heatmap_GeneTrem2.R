rm(list = ls())

require(VennDiagram)
library(xlsx)
library(gplots)
library(org.Mm.eg.db)
setwd("C:/Users/asbarros/Google Drive/RNASeq_icoelho/R_Analysis/")

##########################
#Running the DESeq model #
##########################

source("scripts/DESeq_Model.R",echo = T)

##### Importing the genes #####
dir.create(path = "Heatmap_GenesTrem2/", showWarnings = FALSE)

genes<-read.xlsx("Genes_Trem2.xlsx",sheetIndex=1,header=T)
my_palette <- colorRampPalette(c("blue","white","red"))(n = 50000)

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

design(dm_wt)<- ~ Grouping
dm_wt<-DESeq(dm_wt)

comp<-levels(dm_wt$Grouping)[-c(1,2)]
check1<-vector()

res1<-results(dm_wt,contrast = c("Grouping", comp[3],comp[2]),alpha = 0.05)
res1$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res1),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res1$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res1), column="GENENAME", keytype="ENSEMBL", multiVals="first")

ashr_de1<-lfcShrink(dm_wt,contrast = c("Grouping", comp[3],comp[2]),type="ashr",res=res1)

res1_sub<-res1[,-4]

for (j in 1:ncol(ashr_de1)) {
  check1[j]<-all(as.data.frame(res1_sub)[,j] == as.data.frame(ashr_de1)[,j],na.rm=T)
}

print(check1)

de1<-data.frame(Gene_Symbol = res1$symbol, Gene_Name=res1$geneName,
                PAdj=res1$padj,Log2FC_ashr=ashr_de1$log2FoldChange)


####################################
#Subsetting the DESeq object for KO#
####################################
dm_ko<-dm_sub[,!(dm_sub$Genotype=="WT") | dm_sub$Grouping=="KC_APAP"]
dm_ko$Grouping<-droplevels(dm_ko$Grouping)
dm_ko$Genotype<-droplevels(dm_ko$Genotype)

design(dm_ko)<- ~ Grouping
dm_ko<-DESeq(dm_ko)

comp<-levels(dm_ko$Grouping)[-c(1,2)] #Remove the levels KC_APAP (baseline) & KC (one sample per genotype)
check2<-vector()
res2<-results(dm_ko,contrast = c("Grouping", comp[3],comp[2]),alpha = 0.05)

res2$symbol<-mapIds(org.Mm.eg.db,keys=rownames(res2),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
res2$geneName <- mapIds(org.Mm.eg.db, keys=rownames(res2), column="GENENAME", keytype="ENSEMBL", multiVals="first")


ashr_de2<-lfcShrink(dm_ko,contrast = c("Grouping", comp[3],comp[2]),type="ashr",res=res2)

res2_sub<-res2[,-4]

for (j in 1:ncol(ashr_de2)) {
  check2[j]<-all(as.data.frame(res2_sub)[,j] == as.data.frame(ashr_de2)[,j],na.rm=T)
}

print(check2)

de2<-data.frame(Gene_Symbol = res2$symbol, Gene_Name=res2$geneName,
                PAdj=res2$padj,Log2FC_ashr=ashr_de2$log2FoldChange)


############# make full table ######

res<-data.frame(ENSEMBL=rownames(res2),GeneSymbol=de2$Gene_Symbol,
                WT_Log2FC=de1$Log2FC_ashr, KO_Log2FC=de2$Log2FC_ashr)
res<-res[res$ENSEMBL %in% genes$Ensembl,]

all(res$ENSEMBL==as.character(genes$Ensembl))

res<-res[order(match(res$ENSEMBL,genes$Ensembl)),]
all(res$ENSEMBL==as.character(genes$Ensembl))

res<-as.matrix(data.frame(res[,3:4],row.names = res$GeneSymbol))

tiff("Heatmap_GenesTrem2/heatmap_GenesTrem2.tiff", width = 1500, height = 1700,res = 150)
heatmap.2(res,labRow = rownames(res),notecol="black",density.info='none',col=my_palette,trace="none",margins=c(15,20), dendrogram = 'none', Rowv=NA,Colv = NA,cexCol = 2,cexRow = 2)
dev.off()

res<-res[!(rownames(res) %in% c("Cst3","Cst7")),]

my_palette2 <- colorRampPalette(c("white","red"))(n = 50000)


tiff("Heatmap_GenesTrem2/heatmap_GenesTrem2_nocontrols.tiff", width = 1500, height = 1700,res = 150)
heatmap.2(res,labRow = rownames(res),notecol="black",density.info='none',col=my_palette2,trace="none",margins=c(15,20), dendrogram = 'none', Rowv=NA,Colv = NA,cexCol = 2,cexRow = 2)
dev.off()

