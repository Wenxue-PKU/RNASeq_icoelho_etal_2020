####### Libraries & Setwd #####
library(nichenetr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

library("org.Mm.eg.db")

#### Perform the DESeq analysis to get the dm file
setwd("C:/Users/asbarros/Google Drive/RNASeq_icoelho/R_Analysis/")
source("scripts/DESeq_Model.R",echo=T)
dir.create(path = "NicheNet/", showWarnings = FALSE)

###### 1st) Import the NicheNet network & convert to mouse symbols #######
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
head(ligand_target_matrix)

colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

dim(ligand_target_matrix)

###########
Trst_wt<-dm[,dm$Grouping_Genotype=="Transition_WT"]

dat<-counts(Trst_wt)

sep<-data.frame(V1=seq(10,500,10))
for (i in 1:50){
  j<-i*10
  means<-rowMeans(dat)[rowMeans(dat)>j]
  min<-rowMin(dat)[rowMin(dat)>j]
  sep[i,2]<-length(means)
  sep[i,3]<-length(min)
}

p1<-ggplot(data = sep,aes(x=V1,y=V2))+
  geom_point()+
  geom_point(aes(x=V1,y=V3),color="blue")

p1
sender<-dat[rowMin(dat)>10,]

genes_sender<-mapIds(org.Mm.eg.db,keys=rownames(sender),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
genes_sender<-genes_sender[!is.na(genes_sender)]

#### 4th) Define set of possible ligands #####
#Get the NicheNet network information
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

lr_network$from = lr_network$from %>% convert_human_to_mouse_symbols() 
lr_network$to = lr_network$to %>% convert_human_to_mouse_symbols() 

#Get the unique ligands for the expressed genes for CAF cells
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,genes_sender)

################################
source("scripts/Import_Genes_NicheNet.R",echo=T)
dir.create(path = "NicheNet/Transition_WT", showWarnings = FALSE)
dir.create(path = "NicheNet/Transition_WT/LDECs", showWarnings = FALSE)
dir.create(path = "NicheNet/Transition_WT/LDECs/Plots", showWarnings = FALSE)
dir.create(path = "NicheNet/Transition_WT/LDECs/Results", showWarnings = FALSE)

LDAM<-dm[,dm$Grouping=="LDAM"]

dat<-counts(LDAM)
crit<-rowMeans(dat)[rowMin(dat)>10]
length(crit)
receiver<-dat[rowMin(dat)>10,]

genes_receiver<-mapIds(org.Mm.eg.db,keys=rownames(receiver),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
genes_receiver<-genes_receiver[!is.na(genes_receiver)]

geneset_oi<-as.character(genes$LDECs[genes$LDECs %in% rownames(ligand_target_matrix)])
background_expressed_genes = as.character(genes_receiver[genes_receiver %in% rownames(ligand_target_matrix)])

#Get the unique ligads whose receptors are expressed in the malignant cells
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,background_expressed_genes)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

#All of the expressed ligands from CAF are going to be considered potential ligands
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

#################################
#Now, filter our network data to include only the entries related with the elements from the previous functions
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
lr_network_expressed = lr_network %>% filter(!is.na(from) & !is.na(to)) 

head(lr_network_expressed)

#### 5th) Perform NicheNet's ligand activity prediction on the gene set of interest ####
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)


ligand_activities %>% arrange(-pearson) 

write.xlsx(ligand_activities,file = "NicheNet/Transition_WT/LDECs/Results/TrsntWT_LDECS.xlsx",
           sheetName = "Full Program",append = F)

# In addition, we can also assess the threshold to use the prioritization of top-ligands, using a histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="seagreen4")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(10, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="blue", linetype="dashed", size=1) +
  labs(x="ligand activity (PCC)", y = "# ligands") +
  ggtitle("Histogram for Ligand Distribution according with Pearson Correlation Coefficient (PCC)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

tiff("NicheNet/Transition_WT/LDECs/Plots/PCC_TrsntWT_LDECS_Hist.tiff", width = 1500, height = 1000,res = 150)
p_hist_lig_activity
dev.off()

#### 6th) Infer target genes of top-ranked ligands & visualize them in a heatmap ####

#In this step, we are going to visualize the interaction of the top-20 ranked ligands & both the target genes & 250 top genes (according with the prior model)
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

#For visualization purposes, we are going to atribute a score of zero for each potential score below the 1st quartile (0.25)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

#In this case, we are going to order the ligands as well as the genes for the regulatory potential value & nr. of targets
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized CAF-ligands","p-EMT genes in malignant cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network


####### Best 20 ######
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

#In this step, we are going to visualize the interaction of the top-20 ranked ligands & both the target genes & 250 top genes (according with the prior model)
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

#For visualization purposes, we are going to atribute a score of zero for each potential score below the 1st quartile (0.25)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

#In this case, we are going to order the ligands as well as the genes for the regulatory potential value & nr. of targets
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Transition ligands (WT)","LDEC's Transcript Program", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

tiff("NicheNet/Transition_WT/LDECs/Plots/GenesHeatmap_TrsntWT_LDECS_Best20.tiff", width = 1500, height = 1000,res = 150)
p_ligand_target_network
dev.off()

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

weighted_networks$lr_sig$from = weighted_networks$lr_sig$from %>% convert_human_to_mouse_symbols() 
weighted_networks$lr_sig$to = weighted_networks$lr_sig$to %>% convert_human_to_mouse_symbols() 

lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized CAF-ligands","Receptors expressed by malignant cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network























order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()










ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")

p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized Transition Ligands (WT)","", color = "darkorange",legend_position = "top", x_axis_position = "top")+
  scale_fill_gradient(low = "white",high="blue",breaks=seq(from = 0,to=.08,by = 0.02))+
  ggtitle("Pearson Correlation Coefficient\non predicting target gene expression")+
  theme(axis.text.x = element_blank(),legend.title=element_blank(),plot.title = element_text(hjust = 0.5),axis.ticks = element_blank())

tiff("NicheNet/Transition_WT/LDECs/Plots/PCC_TrsntWT_LDECS_Best20.tiff", width = 400, height = 1000,res = 150)
p_ligand_pearson
dev.off()

######### Best 10 ##########
best_upstream_ligands = ligand_activities %>% top_n(10, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()

ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")

p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized Transition Ligands (WT)","", color = "darkorange",legend_position = "top", x_axis_position = "top")+
  scale_fill_gradient(low = "white",high="blue",breaks=seq(from = 0,to=.08,by = 0.02))+
  ggtitle("Pearson Correlation Coefficient\non predicting target gene expression")+
  theme(axis.text.x = element_blank(),legend.title=element_blank(),plot.title = element_text(hjust = 0.5),axis.ticks = element_blank())

tiff("NicheNet/Transition_WT/LDECs/Plots/PCC_TrsntWT_LDECS_Best10.tiff", width = 400, height = 1000,res = 150)
p_ligand_pearson
dev.off()





#### 6th) Infer target genes of top-ranked ligands & visualize them in a heatmap ####

#In this step, we are going to visualize the interaction of the top-20 ranked ligands & both the target genes & 250 top genes (according with the prior model)
active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links,geneset = geneset_oi, 
         ligand_target_matrix = ligand_target_matrix,n=200) %>% bind_rows()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

#For visualization purposes, we are going to atribute a score of zero for each potential score below the 1st quartile (0.25)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                 ligand_target_matrix = ligand_target_matrix, cutoff = 0.1)

#In this case, we are going to order the ligands as well as the genes for the regulatory potential value & nr. of targets
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% .[. %in% rownames(ligand_target_matrix)] %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized Transition Ligands","Transcriptional Program LDECs", 
                      color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.01,0.02)) +
  ggtitle("Transition on LDECs (WT Genotype)")+
  theme(axis.text.x = element_text(face = "italic"),plot.title = element_text(hjust = 0.5))

p_ligand_target_network







#################################################################################################
Endo<-dm[,dm$Grouping=="Endo_CD26+"]
dat<-counts(Endo)
crit<-rowMin(dat)[rowMeans(dat)>10]
length(crit)
receiver<-dat[rowMin(dat)>10,]

genes_receiver<-mapIds(org.Mm.eg.db,keys=rownames(receiver),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
genes_receiver<-genes_receiver[!is.na(genes_receiver)]

geneset_oi<-as.character(genes$LDECs[genes$LDECs %in% rownames(ligand_target_matrix)])
background_expressed_genes = as.character(genes_receiver[genes_receiver %in% rownames(ligand_target_matrix)])

#Get the unique ligads whose receptors are expressed in the malignant cells
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,background_expressed_genes)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

#All of the expressed ligands from CAF are going to be considered potential ligands
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

#################################
#Now, filter our network data to include only the entries related with the elements from the previous functions
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
lr_network_expressed = lr_network %>% filter(!is.na(from) & !is.na(to)) 

head(lr_network_expressed)

#### 5th) Perform NicheNet's ligand activity prediction on the gene set of interest ####
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)


ligand_activities %>% arrange(-pearson) 

# In addition, we can also assess the threshold to use the prioritization of top-ligands, using a histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity


best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

#### 6th) Infer target genes of top-ranked ligands & visualize them in a heatmap ####

#In this step, we are going to visualize the interaction of the top-20 ranked ligands & both the target genes & 250 top genes (according with the prior model)
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

#For visualization purposes, we are going to atribute a score of zero for each potential score below the 1st quartile (0.25)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

#In this case, we are going to order the ligands as well as the genes for the regulatory potential value & nr. of targets
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized Transition Ligands","Transcriptional Program LDECs", 
                      color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.01,0.02)) +
  ggtitle("Transition on Endothelial (WT Genotype)")+
  theme(axis.text.x = element_text(face = "italic"),plot.title = element_text(hjust = 0.5))

p_ligand_target_network

###############################################################
Common<-dm[,dm$Grouping=="LDAM" | dm$Grouping=="Endo_CD26+"]
dat<-counts(Common)
crit<-rowMin(dat)[rowMin(dat)>10]
length(crit)
receiver<-dat[rowMin(dat)>10,]

genes_receiver<-mapIds(org.Mm.eg.db,keys=rownames(receiver),column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
genes_receiver<-genes_receiver[!is.na(genes_receiver)]

geneset_oi<-as.character(genes$LDECs[genes$LDECs %in% rownames(ligand_target_matrix)])
background_expressed_genes = as.character(genes_receiver[genes_receiver %in% rownames(ligand_target_matrix)])

#Get the unique ligads whose receptors are expressed in the malignant cells
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,background_expressed_genes)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

#All of the expressed ligands from CAF are going to be considered potential ligands
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

#################################
#Now, filter our network data to include only the entries related with the elements from the previous functions
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
lr_network_expressed = lr_network %>% filter(!is.na(from) & !is.na(to)) 

head(lr_network_expressed)

#### 5th) Perform NicheNet's ligand activity prediction on the gene set of interest ####
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)


ligand_activities %>% arrange(-pearson) 

# In addition, we can also assess the threshold to use the prioritization of top-ligands, using a histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity


best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)

#### 6th) Infer target genes of top-ranked ligands & visualize them in a heatmap ####

#In this step, we are going to visualize the interaction of the top-20 ranked ligands & both the target genes & 250 top genes (according with the prior model)
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

#For visualization purposes, we are going to atribute a score of zero for each potential score below the 1st quartile (0.25)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

#In this case, we are going to order the ligands as well as the genes for the regulatory potential value & nr. of targets
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized Transition Ligands","Transcriptional Program LDECs", 
                      color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.01,0.02)) +
  ggtitle("Transition on LDECs & Endothelial (WT Genotype)")+
  theme(axis.text.x = element_text(face = "italic"),plot.title = element_text(hjust = 0.5))

p_ligand_target_network

