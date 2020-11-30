library(biomaRt)
library(org.Mm.eg.db)

ensembl_list <- rownames(df)
ensembl = useEnsembl("ensembl",dataset="mmusculus_gene_ensembl"
                     #,mirror = "uswest"
                     ,version=95)

gene_coords=getBM(attributes=c("ensembl_gene_id","start_position","end_position","mgi_symbol"), 
                  filters="ensembl_gene_id", values=ensembl_list, mart=ensembl,)

gene_coords$size=gene_coords$end_position - gene_coords$start_position
gene_coords$kilo=gene_coords$size / 1e3

gene_coords<-gene_coords[order(match(gene_coords$ensembl_gene_id,rownames(df))),]

all(rownames(df) == gene_coords$ensembl_gene_id)

df_rpk<-matrix(nrow = nrow(df),ncol = ncol(df))

for (i in 1:nrow(df)) {
  df_rpk[i,]<-df[i,] / gene_coords$kilo[i]
}

rownames(df_rpk)<-rownames(df)
colnames(df_rpk)<-colnames(df)

CSums<-colSums(df_rpk)

df_tpm<-matrix(nrow = nrow(df),ncol = ncol(df))

#Transform the data into TPM
for (i in 1:ncol(df)) {
  df_tpm[,i]<-1e6*(df_rpk[,i]/CSums[i])
}

#Re-name the rownames & colnames of the TPM df 
rownames(df_tpm)<-rownames(df)
colnames(df_tpm)<-colnames(df)