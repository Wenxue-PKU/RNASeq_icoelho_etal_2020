library(biomaRt)
library(org.Mm.eg.db)
require(VennDiagram)
library(xlsx)
library(gplots)

########################
my_palette <- colorRampPalette(c("blue","white","red"))(n = 50000)
######################

Trst<-dm[,dm$Grouping=="Transition"]

Trst<-counts(Trst,normalized = FALSE)

ensembl_list <- rownames(Trst)
ensembl = useEnsembl("ensembl",dataset="mmusculus_gene_ensembl"
                     #,mirror = "uswest"
                     ,version=95)

gene_coords=getBM(attributes=c("ensembl_gene_id","start_position","end_position","mgi_symbol"), 
                  filters="ensembl_gene_id", values=ensembl_list, mart=ensembl,)

gene_coords$size=gene_coords$end_position - gene_coords$start_position
gene_coords$kilo=gene_coords$size / 1e3

gene_coords<-gene_coords[order(match(gene_coords$ensembl_gene_id,rownames(Trst))),]

all(rownames(Trst) == gene_coords$ensembl_gene_id)

Trst_rpk<-matrix(nrow = nrow(Trst),ncol = ncol(Trst))

for (i in 1:nrow(Trst)) {
  Trst_rpk[i,]<-Trst[i,] / gene_coords$kilo[i]
}

rownames(Trst_rpk)<-rownames(Trst)
colnames(Trst_rpk)<-colnames(Trst)

CSums<-colSums(Trst_rpk)

Trst_tpm<-matrix(nrow = nrow(Trst),ncol = ncol(Trst))

#Transform the data into TPM
for (i in 1:ncol(Trst)) {
  Trst_tpm[,i]<-1e6*(Trst_rpk[,i]/CSums[i])
}

#Re-name the rownames & colnames of the TPM Trst 
rownames(Trst_tpm)<-rownames(Trst)
colnames(Trst_tpm)<-colnames(Trst)

htmap<-as.data.frame(Trst_tpm)
htmap$symbol<-mapIds(org.Mm.eg.db,keys=rownames(htmap),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
htmap<-htmap[,c("11-WT-Transition-APAP","5-WT-Transition-APAP","WTAPAP1TransitionA","WTAPAP2TransitionB", "WTAPAP3TransitionA",
                "Trem2KOAPAP1TransitionB","Trem2KOAPAP2TransitionA")]
write.xlsx(htmap,file = "NicheNet/GeneExpression_Transition_TPM.xlsx",
           sheetName = "TPM_Expression",append = T)

htmap<-htmap[htmap$symbol %in% best_upstream_ligands, ]
htmap<-htmap[order(match(htmap$symbol,best_upstream_ligands)),]
rownames(htmap)<-htmap$symbol
htmap<-htmap[,-ncol(htmap)]
htmap<-htmap[,c("11-WT-Transition-APAP","5-WT-Transition-APAP","WTAPAP1TransitionA","WTAPAP2TransitionB", "WTAPAP3TransitionA",
                    "Trem2KOAPAP1TransitionB","Trem2KOAPAP2TransitionA")]

htmap<-as.matrix(htmap)

tiff("GeneralPlots/heatmap.tiff", width = 2700, height = 2400,res = 300)
heatmap.2(htmap,labCol=colnames(htmap),labRow = rownames(htmap),notecol="black",density.info='none',col=my_palette,trace="none",
          dendrogram = "none",Rowv = F,margins=c(20,6),Colv = "Rowv")
dev.off()
