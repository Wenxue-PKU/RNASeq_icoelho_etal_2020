#Import the relevant libraries
library("org.Mm.eg.db")
library("xlsx")

#Import the file with the gene symbols concerning the genes of interest
genes<-read.xlsx("Gene_NicheNet_070120.xlsx",sheetIndex = 1,header=T)

#Create a new list
ensembl<-list()

#Save the ensembl ids from the genes of interest file
for(i in 1:3){
  a<-genes[!is.na(genes[,i]),i]
  ensembl[[i]]<-mapIds(org.Mm.eg.db,keys=as.character(a),column="ENSEMBL",keytype = "SYMBOL",multiVals = "first")
}

#Rename the list
names(ensembl)<-colnames(genes)
