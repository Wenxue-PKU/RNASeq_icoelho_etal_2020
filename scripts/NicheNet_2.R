
Trst_wt<-dm[,dm$Grouping_Genotype=="Transition_WT"]

df<-counts(Trst_wt)
colSums(df)

df_tpm<-as.data.frame(matrix(nrow = nrow(df),ncol = ncol(df)))
CSums<-colSums(df)/1E6

for (i in 1:ncol(df)) {
  df_tpm[,i]<-df[,i]/CSums[i]
}
rownames(df_tpm)<-rownames(df)
colnames(df_tpm)<-colnames(df)

crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=4] %>% names()

sep<-data.frame(V1=seq(1:20))
for (i in 1:20){
  crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=i] %>% names()
  sep[i,2]<-length(crit)
}

p1<-ggplot(data = sep,aes(x=V1,y=V2))+
  geom_point()

crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=4] %>% names()

sender<-df[rownames(df) %in% crit,]
genes_sender<-mapIds(org.Mm.eg.db,keys=rownames(sender),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
genes_sender<-genes_sender[!is.na(genes_sender)]

############################################ 
LDAM<-dm[,dm$Grouping=="LDAM"]
df<-counts(LDAM)
colSums(df)

df_tpm<-as.data.frame(matrix(nrow = nrow(df),ncol = ncol(df)))
CSums<-colSums(df)/1E6

for (i in 1:ncol(df)) {
  df_tpm[,i]<-df[,i]/CSums[i]
}
rownames(df_tpm)<-rownames(df)
colnames(df_tpm)<-colnames(df)

crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=4] %>% names()

sep<-data.frame(V1=seq(1:20))
for (i in 1:20){
  crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=i] %>% names()
  sep[i,2]<-length(crit)
}

p1<-ggplot(data = sep,aes(x=V1,y=V2))+
  geom_point()

crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=4] %>% names()

receivers<-df[rownames(df) %in% crit,]
genes_receiver<-mapIds(org.Mm.eg.db,keys=rownames(receivers),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
genes_receiver<-genes_receiver[!is.na(genes_receiver)]

################################
source("scripts/Import_Genes_NicheNet.R",echo=T)
geneset_oi<-as.character(genes$LDECs[genes$LDECs %in% rownames(ligand_target_matrix)])
background_expressed_genes = as.character(genes_receiver[genes_receiver %in% rownames(ligand_target_matrix)])

source("scripts/NicheNet_MainLigands.R")

write.xlsx(ligand_activities,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS_MainLigands.xlsx",
           sheetName = "Full Program",append = F)


tiff("NicheNet/Transition_WT/LDECs/Plots/PCC_TrsntWT_LDECS_Hist.tiff", width = 1500, height = 1000,res = 150)
p_hist_lig_activity
dev.off()

####### Best 20 ######
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
source("scripts/NicheNet_Ligands_to_Genes.R")

tiff("NicheNet/Transition_WT/LDECs/Plots/GenesHeatmap_TrsntWT_LDECS_Hist_Best20_DB500.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_500

write.xlsx(active_ligand_target_links_500,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS_LigandstoGenes.xlsx",
           sheetName = "Best20_DB500",append = F)
dev.off()

tiff("NicheNet/Transition_WT/LDECs/Plots/GenesHeatmap_TrsntWT_LDECS_Hist_Best20_DB250.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_250
dev.off()

write.xlsx(active_ligand_target_links_250,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS_LigandstoGenes.xlsx",
           sheetName = "Best20_DB250",append = T)

tiff("NicheNet/Transition_WT/LDECs/Plots/GenesHeatmap_TrsntWT_LDECS_Hist_Best20_DB100.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_100

write.xlsx(active_ligand_target_links_100,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS_LigandstoGenes.xlsx",
           sheetName = "Best20_DB100",append = T)
dev.off()

########### 
source("scripts/NicheNet_Ligands_to_Receptors.R")
tiff("NicheNet/Transition_WT/LDECs/Plots/ReceptHeatmap_TrsntWT_LDECS_Hist_Best20.tiff", width = 2000, height = 1000,res = 150)
p_ligand_receptor_network
write.xlsx(vis_ligand_receptor_network,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS_LigandstoReceptors.xlsx",
           sheetName = "Best20",append = F)
dev.off()

####### Best 10 ######
best_upstream_ligands = ligand_activities %>% top_n(10, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
source("scripts/NicheNet_Ligands_to_Genes.R")

tiff("NicheNet/Transition_WT/LDECs/Plots/GenesHeatmap_TrsntWT_LDECS_Hist_Best10_DB500.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_500

write.xlsx(active_ligand_target_links_500,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS_LigandstoGenes.xlsx",
           sheetName = "Best10_DB500",append = T)
dev.off()

tiff("NicheNet/Transition_WT/LDECs/Plots/GenesHeatmap_TrsntWT_LDECS_Hist_Best10_DB250.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_250
dev.off()

write.xlsx(active_ligand_target_links_250,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS_LigandstoGenes.xlsx",
           sheetName = "Best10_DB250",append = T)

tiff("NicheNet/Transition_WT/LDECs/Plots/GenesHeatmap_TrsntWT_LDECS_Hist_Best10_DB100.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_100

write.xlsx(active_ligand_target_links_100,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS_LigandstoGenes.xlsx",
           sheetName = "Best10_DB100",append = T)
dev.off()

########### 
source("scripts/NicheNet_Ligands_to_Receptors.R")
tiff("NicheNet/Transition_WT/LDECs/Plots/ReceptHeatmap_TrsntWT_LDECS_Hist_Best10.tiff", width = 2000, height = 1000,res = 150)
p_ligand_receptor_network
write.xlsx(vis_ligand_receptor_network,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS_LigandstoReceptors.xlsx",
           sheetName = "Best10",append = F)
dev.off()

############################################ 
dir.create(path = "NicheNet/Transition_WT", showWarnings = FALSE)
dir.create(path = "NicheNet/Transition_WT/Endo", showWarnings = FALSE)
dir.create(path = "NicheNet/Transition_WT/Endo/Plots", showWarnings = FALSE)
dir.create(path = "NicheNet/Transition_WT/Endo/Results", showWarnings = FALSE)

Endo<-dm[,dm$Grouping=="Endo_CD26+"]

df<-counts(Endo)
colSums(df)

df_tpm<-as.data.frame(matrix(nrow = nrow(df),ncol = ncol(df)))
CSums<-colSums(df)/1E6

for (i in 1:ncol(df)) {
  df_tpm[,i]<-df[,i]/CSums[i]
}
rownames(df_tpm)<-rownames(df)
colnames(df_tpm)<-colnames(df)

crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=4] %>% names()

sep<-data.frame(V1=seq(1:20))
for (i in 1:20){
  crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=i] %>% names()
  sep[i,2]<-length(crit)
}

p1<-ggplot(data = sep,aes(x=V1,y=V2))+
  geom_point()

crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=4] %>% names()

receivers<-df[rownames(df) %in% crit,]
genes_receiver<-mapIds(org.Mm.eg.db,keys=rownames(receivers),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
genes_receiver<-genes_receiver[!is.na(genes_receiver)]

################################
source("scripts/Import_Genes_NicheNet.R",echo=T)
geneset_oi<-as.character(genes$Endothelial[genes$Endothelial %in% rownames(ligand_target_matrix)])
background_expressed_genes = as.character(genes_receiver[genes_receiver %in% rownames(ligand_target_matrix)])

source("scripts/NicheNet_MainLigands.R")

write.xlsx(ligand_activities,file = "NicheNet/Transition_WT/Endo/Results/TrsntWT_Endo_MainLigands.xlsx",
           sheetName = "Full Program",append = F)


tiff("NicheNet/Transition_WT/Endo/Plots/PCC_TrsntWT_Endo_Hist.tiff", width = 1500, height = 1000,res = 150)
p_hist_lig_activity
dev.off()

####### Best 20 ######
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
source("scripts/NicheNet_Ligands_to_Genes.R")

tiff("NicheNet/Transition_WT/Endo/Plots/GenesHeatmap_TrsntWT_Endo_Hist_Best20_DB500.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_500

write.xlsx(active_ligand_target_links_500,file = "NicheNet/Transition_WT/Endo/Results/TrsntWT_Endo_LigandstoGenes.xlsx",
           sheetName = "Best20_DB500",append = F)
dev.off()

tiff("NicheNet/Transition_WT/Endo/Plots/GenesHeatmap_TrsntWT_Endo_Hist_Best20_DB250.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_250
dev.off()

write.xlsx(active_ligand_target_links_250,file = "NicheNet/Transition_WT/Endo/Results/TrsntWT_Endo_LigandstoGenes.xlsx",
           sheetName = "Best20_DB250",append = T)

tiff("NicheNet/Transition_WT/Endo/Plots/GenesHeatmap_TrsntWT_Endo_Hist_Best20_DB100.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_100

write.xlsx(active_ligand_target_links_100,file = "NicheNet/Transition_WT/Endo/Results/TrsntWT_Endo_LigandstoGenes.xlsx",
           sheetName = "Best20_DB100",append = T)
dev.off()

#### 
source("scripts/NicheNet_Ligands_to_Receptors.R")
tiff("NicheNet/Transition_WT/Endo/Plots/ReceptHeatmap_TrsntWT_Endo_Hist_Best20.tiff", width = 2000, height = 1000,res = 150)
p_ligand_receptor_network
write.xlsx(vis_ligand_receptor_network,file = "NicheNet/Transition_WT/Endo/Results/TrsntWT_Endo_LigandstoReceptors.xlsx",
           sheetName = "Best20",append = F)
dev.off()

####### Best 10 ######
best_upstream_ligands = ligand_activities %>% top_n(10, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
source("scripts/NicheNet_Ligands_to_Genes.R")

tiff("NicheNet/Transition_WT/Endo/Plots/GenesHeatmap_TrsntWT_Endo_Hist_Best10_DB500.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_500

write.xlsx(active_ligand_target_links_500,file = "NicheNet/Transition_WT/Endo/Results/TrsntWT_Endo_LigandstoGenes.xlsx",
           sheetName = "Best10_DB500",append = T)
dev.off()

tiff("NicheNet/Transition_WT/Endo/Plots/GenesHeatmap_TrsntWT_Endo_Hist_Best10_DB250.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_250
dev.off()

write.xlsx(active_ligand_target_links_250,file = "NicheNet/Transition_WT/Endo/Results/TrsntWT_Endo_LigandstoGenes.xlsx",
           sheetName = "Best10_DB250",append = T)

tiff("NicheNet/Transition_WT/Endo/Plots/GenesHeatmap_TrsntWT_Endo_Hist_Best10_DB100.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_100

write.xlsx(active_ligand_target_links_100,file = "NicheNet/Transition_WT/Endo/Results/TrsntWT_Endo_LigandstoGenes.xlsx",
           sheetName = "Best10_DB100",append = T)
dev.off()

########### 
source("scripts/NicheNet_Ligands_to_Receptors.R")
tiff("NicheNet/Transition_WT/Endo/Plots/ReceptHeatmap_TrsntWT_Endo_Hist_Best10.tiff", width = 2000, height = 1000,res = 150)
p_ligand_receptor_network
write.xlsx(vis_ligand_receptor_network,file = "NicheNet/Transition_WT/Endo/Results/TrsntWT_Endo_LigandstoReceptors.xlsx",
           sheetName = "Best10",append = F)
dev.off()


############################################ 
dir.create(path = "NicheNet/Transition_WT", showWarnings = FALSE)
dir.create(path = "NicheNet/Transition_WT/Common", showWarnings = FALSE)
dir.create(path = "NicheNet/Transition_WT/Common/Plots", showWarnings = FALSE)
dir.create(path = "NicheNet/Transition_WT/Common/Results", showWarnings = FALSE)

Com<-dm[,dm$Grouping=="Endo_CD26+" | dm$Grouping=="LDAM"]

df<-counts(Com)
colSums(df)

df_tpm<-as.data.frame(matrix(nrow = nrow(df),ncol = ncol(df)))
CSums<-colSums(df)/1E6

for (i in 1:ncol(df)) {
  df_tpm[,i]<-df[,i]/CSums[i]
}
rownames(df_tpm)<-rownames(df)
colnames(df_tpm)<-colnames(df)

crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=4] %>% names()

sep<-data.frame(V1=seq(1:20))
for (i in 1:20){
  crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=i] %>% names()
  sep[i,2]<-length(crit)
}

p1<-ggplot(data = sep,aes(x=V1,y=V2))+
  geom_point()

crit<- df_tpm %>% apply(1,function(x){log2(mean(x) + 1)}) %>% .[. >=4] %>% names()

receivers<-df[rownames(df) %in% crit,]
genes_receiver<-mapIds(org.Mm.eg.db,keys=rownames(receivers),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
genes_receiver<-genes_receiver[!is.na(genes_receiver)]

################################
source("scripts/Import_Genes_NicheNet.R",echo=T)
geneset_oi<-as.character(genes$Common[genes$Common %in% rownames(ligand_target_matrix)])
background_expressed_genes = as.character(genes_receiver[genes_receiver %in% rownames(ligand_target_matrix)])

source("scripts/NicheNet_MainLigands.R")

write.xlsx(ligand_activities,file = "NicheNet/Transition_WT/Common/Results/TrsntWT_Common_MainLigands.xlsx",
           sheetName = "Full Program",append = F)


tiff("NicheNet/Transition_WT/Common/Plots/PCC_TrsntWT_Common_Hist.tiff", width = 1500, height = 1000,res = 150)
p_hist_lig_activity
dev.off()

####### Best 20 ######
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
source("scripts/NicheNet_Ligands_to_Genes.R")

tiff("NicheNet/Transition_WT/Common/Plots/GenesHeatmap_TrsntWT_Common_Hist_Best20_DB500.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_500

write.xlsx(active_ligand_target_links_500,file = "NicheNet/Transition_WT/Common/Results/TrsntWT_Common_LigandstoGenes.xlsx",
           sheetName = "Best20_DB500",append = F)
dev.off()

tiff("NicheNet/Transition_WT/Common/Plots/GenesHeatmap_TrsntWT_Common_Hist_Best20_DB250.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_250
dev.off()

write.xlsx(active_ligand_target_links_250,file = "NicheNet/Transition_WT/Common/Results/TrsntWT_Common_LigandstoGenes.xlsx",
           sheetName = "Best20_DB250",append = T)

tiff("NicheNet/Transition_WT/Common/Plots/GenesHeatmap_TrsntWT_Common_Hist_Best20_DB100.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_100

write.xlsx(active_ligand_target_links_100,file = "NicheNet/Transition_WT/Common/Results/TrsntWT_Common_LigandstoGenes.xlsx",
           sheetName = "Best20_DB100",append = T)
dev.off()

#### 
source("scripts/NicheNet_Ligands_to_Receptors.R")
tiff("NicheNet/Transition_WT/Common/Plots/ReceptHeatmap_TrsntWT_Common_Hist_Best20.tiff", width = 2000, height = 1000,res = 150)
p_ligand_receptor_network
write.xlsx(vis_ligand_receptor_network,file = "NicheNet/Transition_WT/Common/Results/TrsntWT_Common_LigandstoReceptors.xlsx",
           sheetName = "Best20",append = F)
dev.off()

####### Best 10 ######
best_upstream_ligands = ligand_activities %>% top_n(10, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
source("scripts/NicheNet_Ligands_to_Genes.R")

tiff("NicheNet/Transition_WT/Common/Plots/GenesHeatmap_TrsntWT_Common_Hist_Best10_DB500.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_500

write.xlsx(active_ligand_target_links_500,file = "NicheNet/Transition_WT/Common/Results/TrsntWT_Common_LigandstoGenes.xlsx",
           sheetName = "Best10_DB500",append = T)
dev.off()

tiff("NicheNet/Transition_WT/Common/Plots/GenesHeatmap_TrsntWT_Common_Hist_Best10_DB250.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_250
dev.off()

write.xlsx(active_ligand_target_links_250,file = "NicheNet/Transition_WT/Common/Results/TrsntWT_Common_LigandstoGenes.xlsx",
           sheetName = "Best10_DB250",append = T)

tiff("NicheNet/Transition_WT/Common/Plots/GenesHeatmap_TrsntWT_Common_Hist_Best10_DB100.tiff", width = 1700, height = 1000,res = 150)
p_ligand_target_network_100

write.xlsx(active_ligand_target_links_100,file = "NicheNet/Transition_WT/Common/Results/TrsntWT_Common_LigandstoGenes.xlsx",
           sheetName = "Best10_DB100",append = T)
dev.off()

########### 
source("scripts/NicheNet_Ligands_to_Receptors.R")
tiff("NicheNet/Transition_WT/Common/Plots/ReceptHeatmap_TrsntWT_Common_Hist_Best10.tiff", width = 2000, height = 1000,res = 150)
p_ligand_receptor_network
write.xlsx(vis_ligand_receptor_network,file = "NicheNet/Transition_WT/Common/Results/TrsntWT_Common_LigandstoReceptors.xlsx",
           sheetName = "Best10",append = F)
dev.off()
























