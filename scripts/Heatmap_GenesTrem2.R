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

genes<-read.xlsx("Genes_Trem2.xlsx",sheetIndex=1,header=F)
temp<-data.frame(ENSEMBL=rownames(tab),
                 gene_symbol=mapIds(org.Mm.eg.db,keys=rownames(tab),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first"))


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
