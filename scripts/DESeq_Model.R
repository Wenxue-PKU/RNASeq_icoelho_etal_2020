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
setwd("C:/Users/asbarros/Google Drive/RNASeq_icoelho/R_Analysis/")
dir.create(path = "SessionInfo/", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()),paste0("SessionInfo/SessionInfo_",Sys.Date(),".txt",sep=""))

#Import the files
gc_tab_nov18<-"gene_count_Nov18/"
files_nov18<-paste(gc_tab_nov18,dir(gc_tab_nov18,pattern="tab"),sep="")

gc_tab_march19<-"gene_count_March19/"
files_march19<-paste(gc_tab_march19,dir(gc_tab_march19,pattern="tab"),sep="")

file.names<-c(files_nov18,files_march19)

outfile<-list()
for (i in file.names) {
  outfile[[i]]<-read.table(i,row.names=1)
}

#Remove the first four lines of the files, that contain unmappeds, dubious, and other artifacts
for (i in 1:length(outfile)) {
  outfile[[i]]<-outfile[[i]][-c(1:4),]
}

#Select one of the columns to use as the expression data (normally, select the one w/ most reads)
tab<-data.frame(outfile[[1]][,1])
rownames(tab)<-rownames(outfile[[1]])

for (i in 2:length(outfile)) {
  tab[,i]<-outfile[[i]][,1]
  rownames(tab)<-rownames(outfile[[1]])
}

#define the column names
colnames(tab)<-gsub(".*/","",names(outfile))
colnames(tab)<-gsub("_.*","",colnames(tab))

#import the design matrix
design<-read.xlsx("metadata.xlsx",sheetIndex = 1,header=T,row.names = 1)

all(rownames(design)==colnames(tab)) #check if the rownames of design match exactly w/ the column names

design<-design[order(match(rownames(design),colnames(tab))),]

all(rownames(design)==colnames(tab)) #check if the rownames of design match exactly w/ the column names

design$Genotype<-relevel(x = design$Genotype,ref = "WT")
design$LibrarySize<-colSums(tab)

design<-design[!(rownames(design) %in% c("3-WT-KC-control","2-TREM2KO-LDAMs-APAP","7-TREM2KO-KC-control","10-WT-LDAMs-APAP")),]
tab<-tab[,!(colnames(tab) %in% c("3-WT-KC-control","2-TREM2KO-LDAMs-APAP","7-TREM2KO-KC-control","10-WT-LDAMs-APAP")),]

all(rownames(design)==colnames(tab)) #check if the rownames of design match exactly w/ the column names

#DESeqResults
#########################################################################################asbarros
dir.create(path = "GeneralPlots/", showWarnings = FALSE)

tiff("GeneralPlots/LibrarySize_Plot.tiff", width = 1500, height = 1000,res = 150)

ggplot(data=design, aes(x=rownames(design),y=(design$LibrarySize/1E6)))+
  geom_bar(stat="identity",fill="dodgerblue4")+
  ylab("Total Reads aligned to Genes (Million Reads)")+
  xlab("Samples")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

dm <- DESeqDataSetFromMatrix(countData = tab, colData = design, design = ~ Batch)
dm<-dm[rowSums(counts(dm)) > 0 , ] #remove the genes w/ 0 reads
dm<-estimateSizeFactors(dm) #each sample has a Size Factor, which will help further on in the normalization step


dvst<-vst(dm,blind=F)
dm<-estimateDispersions(dm)
plotDispEsts(dm,ylim=c(1e-6, 1e2))

ddist<-dist(t(assay(dvst)))
heatmap1<-pheatmap(ddist,cluster_rows=T,cluster_cols=T,annotation_col = design[,c("Batch","Grouping","Genotype")])

pca1<-plotPCA(dvst, intgroup = c("Grouping","Batch"),returnData=T)
percentVar<-round(100*attr(pca1,"percentVar"))

data1<-ggplot(pca1,aes(PC1,PC2,COlour=Grouping,shape=Batch))+
  #geom_point(aes(colour=Grouping),size=3)+
  geom_text(aes(colour=Grouping),label=rownames(design))+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

data1

##################### We Need to remove a specific sample #####################################################
design<-design[!(rownames(design) %in% c("WTAPAP3LDAMa")),]
tab<-tab[,!(colnames(tab) %in% c("WTAPAP3LDAMa")),]

all(rownames(design)==colnames(tab)) #check if the rownames of design match exactly w/ the column names

dm <- DESeqDataSetFromMatrix(countData = tab, colData = design, design = ~ Grouping)
dm<-dm[rowSums(counts(dm)) > 0 , ] #remove the genes w/ 0 reads
dm<-estimateSizeFactors(dm) #each sample has a Size Factor, which will help further on in the normalization step


dvst<-vst(dm,blind=F)
dm<-estimateDispersions(dm)
plotDispEsts(dm,ylim=c(1e-6, 1e2))
dev.off()

ddist<-dist(t(assay(dvst)))

tiff("GeneralPlots/Heatmap_AllSamples.tiff", width = 1500, height = 1000,res = 150)
heatmap1<-pheatmap(ddist,cluster_rows=T,cluster_cols=T,annotation_col = design[,c("Batch","Grouping","Genotype")])
print(heatmap1)

dev.off()

pca1<-plotPCA(dvst, intgroup = c("Grouping","Batch"),returnData=T)
percentVar<-round(100*attr(pca1,"percentVar"))

data1<-ggplot(pca1,aes(PC1,PC2,COlour=Grouping,shape=Batch))+
  #geom_point(aes(colour=Grouping),size=3)+
  geom_text(aes(colour=Grouping),label=rownames(design))+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

tiff("GeneralPlots/PCA_AllSamples_Names.tiff", width = 1500, height = 1000,res = 150)
data1
dev.off()

data2<-ggplot(pca1,aes(PC1,PC2,COlour=Grouping,shape=Batch))+
  geom_point(aes(colour=Grouping),size=3)+
  #geom_text(aes(colour=Grouping),label=rownames(design))+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

tiff("GeneralPlots/PCA_AllSamples_Grouping&Batch.tiff", width = 1500, height = 1000,res = 150)
data2
dev.off()

pca1<-plotPCA(dvst, intgroup = c("Grouping","Genotype"),returnData=T)
data3<-ggplot(pca1,aes(PC1,PC2,Colour=Grouping,shape=Genotype))+
  geom_point(aes(colour=Grouping),size=3)+
  #geom_text(aes(colour=Grouping),label=rownames(design))+
  xlab(paste0("PC1: ",percentVar[1],"% Variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% Variance"))+
  coord_fixed()+
  ggtitle("Principal Component Analysis based on Overall Genetic Expression across all samples")+
  theme(plot.title = element_text(hjust = 0.5),legend.position="bottom")

tiff("GeneralPlots/PCA_AllSamples_Grouping&Genotype.tiff", width = 1500, height = 1000,res = 150)
data3
dev.off()

design(dm)<- ~ Grouping
dm<-DESeq(dm)

dir.create(path = "variables/", showWarnings = FALSE)
dir.create(path = "environments/", showWarnings = FALSE)
saveRDS(dm, paste0("variables/dm_",Sys.Date(),".rds",sep="")) #remember that you need to use readRDS function to read this file
save.image(file = paste0("environments/DESeq_Model_",Sys.Date(),".RData",sep=""))
