rm(list=ls())
####### Libraries & Setwd #####
library(nichenetr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("org.Mm.eg.db")
library("ggplot2")
library("genefilter")
library("xlsx")
library("stringr")
library("vsn")

####### Import the necessary variables ########
setwd("C:/Users/asbarros/Google Drive/RNASeq_icoelho/R_Analysis/")
dm<-readRDS("variables/dm_2020-01-21.rds")
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

# Now, it's important to transform the gene information from human to mouse
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

lr_network$from = lr_network$from %>% convert_human_to_mouse_symbols() 
lr_network$to = lr_network$to %>% convert_human_to_mouse_symbols() 

weighted_networks$lr_sig$from = weighted_networks$lr_sig$from %>% convert_human_to_mouse_symbols() 
weighted_networks$lr_sig$to = weighted_networks$lr_sig$to %>% convert_human_to_mouse_symbols() 

############ Run the respective central scripts ####################
source("scripts/NicheNet_scripts/NicheNet_TrstWT.R",echo=T)
source("scripts/NicheNet_scripts/NicheNet_TrstKO.R",echo=T)
