#In order to plot the ligands in the proper order, we need to reverse the order of the best_ligands objects
order_ligands = best_upstream_ligands %>% rev()

#create a new matrix with the PCC for ligands
ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% set_rownames(ligand_activities$test_ligand)

#Structure the matrix to get the ligands in order
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")

#Now, produce the heatmap
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized Ligands ","", color = "darkorange",legend_position = "top", x_axis_position = "top")+
  scale_fill_gradient(low = "white",high="blue")+
  ggtitle("Pearson Correlation Coefficient\non predicting target gene expression")+
  theme(axis.text.x = element_blank(),legend.title=element_blank(),plot.title = element_text(hjust = 0.5),axis.ticks = element_blank(),legend.key.width = unit(1,"cm"))

