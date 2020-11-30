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
dir.create(path = "DE_Macrophages/VulcanoPlots/", showWarnings = FALSE)
dir.create(path = "DE_Macrophages/VulcanoPlots/WT", showWarnings = FALSE)
dir.create(path = "DE_Macrophages/VulcanoPlots/KO", showWarnings = FALSE)

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

#################### Diff. Expressed #######################

comp<-levels(dm_wt$Grouping)[-c(1,2)]
check1<-vector()
check2<-vector()
check3<-vector()

res1<-results(dm_wt,contrast = c("Grouping", comp[2],comp[1]),alpha = 0.05)
res2<-results(dm_wt,contrast = c("Grouping", comp[3],comp[1]),alpha = 0.05)
res3<-results(dm_wt,contrast = c("Grouping", comp[3],comp[2]),alpha = 0.05)

print(summary(res1))
print(summary(res2))
print(summary(res3))

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

####### Vulcano Plot for de1 #####
de1$PAdj[is.na(de1$PAdj)]<-1
de1$threshold<-"NoVariance"
de1$threshold[de1$PAdj < 0.05 & de1$Log2FC_ashr < 0] <- c("Down-Regulated")
de1$threshold[de1$PAdj < 0.05 & de1$Log2FC_ashr > 0]<- c("Up-Regulated")

de1$threshold<-as.factor(de1$threshold)
de1$threshold<-factor(de1$threshold, levels(de1$threshold)
                      [c(2,1,3)])
tiff(paste0("DE_Macrophages/VulcanoPlots/WT/", paste0("Vulcano_Plot_WT_", comp[2]," _vs_",comp[1],".tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de1, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-40, 40)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[2]," vs ",comp[1], " (WT)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =18, label = length(de1$threshold[de1$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -12, y =18, label = length(de1$threshold[de1$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

de1$PAdj[-log(de1$PAdj,base=10) > 20]<-1E-20
tiff(paste0("DE_Macrophages/VulcanoPlots/WT/", paste0("Vulcano_Plot_WT_", comp[2]," _vs_",comp[1],"_Threshold.tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de1, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-40, 40)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[2]," vs ",comp[1]," (WT)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =18, label = length(de1$threshold[de1$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -12, y = 18, label = length(de1$threshold[de1$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

####### Vulcano Plot for de2 #####
de2<-data.frame(Gene_Symbol = res2$symbol, Gene_Name=res2$geneName,
                PAdj=res2$padj,Log2FC_ashr=ashr_de2$log2FoldChange)

de2$PAdj[is.na(de2$PAdj)]<-1
de2$threshold<-"NoVariance"
de2$threshold[de2$PAdj < 0.05 & de2$Log2FC_ashr < 0] <- c("Down-Regulated")
de2$threshold[de2$PAdj < 0.05 & de2$Log2FC_ashr > 0]<- c("Up-Regulated")

de2$threshold<-as.factor(de2$threshold)
de2$threshold<-factor(de2$threshold, levels(de2$threshold)
                      [c(2,1,3)])
tiff(paste0("DE_Macrophages/VulcanoPlots/WT/", paste0("Vulcano_Plot_WT_", comp[3]," _vs_",comp[1],".tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de2, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-45, 45)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[3]," vs ",comp[1], " (WT)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =18, label = length(de2$threshold[de2$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -12, y = 18, label = length(de2$threshold[de2$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

de2$PAdj[-log(de2$PAdj,base=10) > 20]<-1E-20
tiff(paste0("DE_Macrophages/VulcanoPlots/WT/", paste0("Vulcano_Plot_WT_", comp[3]," _vs_",comp[1],"_Threshold.tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de2, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-45, 45)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[3]," vs ",comp[1], " (WT)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =18, label = length(de2$threshold[de2$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -12, y = 18, label = length(de2$threshold[de2$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

####### Vulcano Plot for de3 #####
de3<-data.frame(Gene_Symbol = res3$symbol, Gene_Name=res3$geneName,
                PAdj=res3$padj,Log2FC_ashr=ashr_de3$log2FoldChange)

de3$PAdj[is.na(de3$PAdj)]<-1
de3$threshold<-"NoVariance"
de3$threshold[de3$PAdj < 0.05 & de3$Log2FC_ashr < 0] <- c("Down-Regulated")
de3$threshold[de3$PAdj < 0.05 & de3$Log2FC_ashr > 0]<- c("Up-Regulated")

de3$threshold<-as.factor(de3$threshold)
de3$threshold<-factor(de3$threshold, levels(de3$threshold)
                      [c(2,1,3)])
tiff(paste0("DE_Macrophages/VulcanoPlots/WT/", paste0("Vulcano_Plot_WT_", comp[3]," _vs_",comp[2],".tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de3, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-35, 35)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[3]," vs ",comp[2], " (WT)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =22, label = length(de3$threshold[de3$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -12, y = 22, label = length(de3$threshold[de3$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

de3$PAdj[-log(de3$PAdj,base=10) > 20]<-1E-20
tiff(paste0("DE_Macrophages/VulcanoPlots/WT/", paste0("Vulcano_Plot_WT_", comp[3]," _vs_",comp[2],"_Threshold.tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de3, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-35, 35)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[3]," vs ",comp[2], " (WT)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =18, label = length(de3$threshold[de3$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -12, y = 18, label = length(de3$threshold[de3$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

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

tiff("KC_APAP_vs_Macrophages/Plots/PCA_Trem2KO.tiff", width = 1500, height = 1000,res = 150)
data1
dev.off()

design(dm_ko)<- ~ Grouping
dm_ko<-DESeq(dm_ko)


#################### Diff. Expressed #######################

comp<-levels(dm_ko$Grouping)[-c(1,2)] #Remove the levels KC_APAP (baseline) & KC (one sample per genotype)
check1<-vector()
check2<-vector()
check3<-vector()

res1<-results(dm_ko,contrast = c("Grouping", comp[2],comp[1]),alpha = 0.05)
res2<-results(dm_ko,contrast = c("Grouping", comp[3],comp[1]),alpha = 0.05)
res3<-results(dm_ko,contrast = c("Grouping", comp[3],comp[2]),alpha = 0.05)

print(summary(res1))
print(summary(res2))
print(summary(res3))


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


####### Vulcano Plot for de1 #####
de1$PAdj[is.na(de1$PAdj)]<-1
de1$threshold<-"NoVariance"
de1$threshold[de1$PAdj < 0.05 & de1$Log2FC_ashr < 0] <- c("Down-Regulated")
de1$threshold[de1$PAdj < 0.05 & de1$Log2FC_ashr > 0]<- c("Up-Regulated")

de1$threshold<-as.factor(de1$threshold)
de1$threshold<-factor(de1$threshold, levels(de1$threshold)
                      [c(2,1,3)])
tiff(paste0("DE_Macrophages/VulcanoPlots/KO/", paste0("Vulcano_Plot_KO_", comp[2]," _vs_",comp[1],".tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de1, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-25, 25)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[2]," vs ",comp[1], " (KO)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 15, y =10, label = length(de1$threshold[de1$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -15, y = 10, label = length(de1$threshold[de1$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

de1$PAdj[-log(de1$PAdj,base=10) > 20]<-1E-20
tiff(paste0("DE_Macrophages/VulcanoPlots/KO/", paste0("Vulcano_Plot_KO_", comp[2]," _vs_",comp[1],"_Threshold.tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de1, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-25, 25)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[2]," vs ",comp[1]," (KO)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 15, y =10, label = length(de1$threshold[de1$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -15, y = 10, label = length(de1$threshold[de1$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

####### Vulcano Plot for de2 #####
de2<-data.frame(Gene_Symbol = res2$symbol, Gene_Name=res2$geneName,
                PAdj=res2$padj,Log2FC_ashr=ashr_de2$log2FoldChange)

de2$PAdj[is.na(de2$PAdj)]<-1
de2$threshold<-"NoVariance"
de2$threshold[de2$PAdj < 0.05 & de2$Log2FC_ashr < 0] <- c("Down-Regulated")
de2$threshold[de2$PAdj < 0.05 & de2$Log2FC_ashr > 0]<- c("Up-Regulated")

de2$threshold<-as.factor(de2$threshold)
de2$threshold<-factor(de2$threshold, levels(de2$threshold)
                      [c(2,1,3)])
tiff(paste0("DE_Macrophages/VulcanoPlots/KO/", paste0("Vulcano_Plot_KO_", comp[3]," _vs_",comp[1],".tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de2, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[3]," vs ",comp[1], " (KO)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =18, label = length(de2$threshold[de2$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -12, y = 18, label = length(de2$threshold[de2$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

de2$PAdj[-log(de2$PAdj,base=10) > 20]<-1E-20
tiff(paste0("DE_Macrophages/VulcanoPlots/KO/", paste0("Vulcano_Plot_KO_", comp[3]," _vs_",comp[1],"_Threshold.tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de2, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[3]," vs ",comp[1], " (KO)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 12, y =18, label = length(de2$threshold[de2$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -12, y = 18, label = length(de2$threshold[de2$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

####### Vulcano Plot for de3 #####
de3<-data.frame(Gene_Symbol = res3$symbol, Gene_Name=res3$geneName,
                PAdj=res3$padj,Log2FC_ashr=ashr_de3$log2FoldChange)

de3$PAdj[is.na(de3$PAdj)]<-1
de3$threshold<-"NoVariance"
de3$threshold[de3$PAdj < 0.05 & de3$Log2FC_ashr < 0] <- c("Down-Regulated")
de3$threshold[de3$PAdj < 0.05 & de3$Log2FC_ashr > 0]<- c("Up-Regulated")

de3$threshold<-as.factor(de3$threshold)
de3$threshold<-factor(de3$threshold, levels(de3$threshold)
                      [c(2,1,3)])
tiff(paste0("DE_Macrophages/VulcanoPlots/KO/", paste0("Vulcano_Plot_KO_", comp[3]," _vs_",comp[2],".tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de3, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-30, 30)) +
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[3]," vs ",comp[2], " (KO)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 15, y =12, label = length(de3$threshold[de3$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -15, y = 12, label = length(de3$threshold[de3$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()

de3$PAdj[-log(de3$PAdj,base=10) > 20]<-1E-20
tiff(paste0("DE_Macrophages/VulcanoPlots/KO/", paste0("Vulcano_Plot_KO_", comp[3]," _vs_",comp[2],"_Threshold.tiff",sep="")), width = 1500, height = 1000,res = 150)

ggplot(data=de3, aes(x=Log2FC_ashr, y=-log(PAdj,base=10), colour=threshold)) +
  geom_point(alpha=0.4, size=3) +
  xlim(c(-30, 30))+
  xlab("log2 fold change") + ylab("-log10(p-value)") +
  scale_color_manual(values=c("lightgrey","blue", "red")) +
  geom_hline(yintercept = -log(0.05,base=10),colour="blue", linetype = "longdash") +
  labs(title = paste0("Vulcano plot for ", comp[3]," vs ",comp[2], " (KO)",sep="")) + 
  labs(colour="Expression Level") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate("text", x = 15, y =12, label = length(de3$threshold[de3$threshold=="Up-Regulated"]),size=6)+
  annotate("text", x = -15, y = 12, label = length(de3$threshold[de3$threshold=="Down-Regulated"]),size=6)+ 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 1),
        plot.title=element_text(hjust=0.5,size=18), 
        legend.position = "bottom", legend.key = element_blank(),text=element_text(size=16))
dev.off()