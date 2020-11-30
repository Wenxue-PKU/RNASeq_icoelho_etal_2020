########### DB @ 500 ##############
#In this step, we are going to visualize the interaction of the top-20 ranked ligands & both the target genes & 250 top genes (according with the prior model)
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 500) %>% bind_rows()
active_ligand_target_links_df<-active_ligand_target_links_df[!is.na(active_ligand_target_links_df$target),]

#For visualization purposes, we are going to atribute a score of zero for each potential score below the 1st quartile (0.25)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.0)

active_ligand_target_links_500<-active_ligand_target_links

#In this case, we are going to order the ligands as well as the genes for the regulatory potential value & nr. of targets
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network_500 = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Ligands","Transcriptional Program", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple") + theme(axis.text.x = element_text(face = "italic"),legend.key.width = unit(1,"cm"))

########### DB @ 250 ##############
#In this step, we are going to visualize the interaction of the top-20 ranked ligands & both the target genes & 250 top genes (according with the prior model)
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df<-active_ligand_target_links_df[!is.na(active_ligand_target_links_df$target),]

#For visualization purposes, we are going to atribute a score of zero for each potential score below the 1st quartile (0.25)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.0)

active_ligand_target_links_250<-active_ligand_target_links

#In this case, we are going to order the ligands as well as the genes for the regulatory potential value & nr. of targets
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network_250 = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Ligands","Transcriptional Program", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple") + theme(axis.text.x = element_text(face = "italic"),legend.key.width = unit(1,"cm"))

####### DB @ 100 ##################
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 100) %>% bind_rows()
active_ligand_target_links_df<-active_ligand_target_links_df[!is.na(active_ligand_target_links_df$target),]

#For visualization purposes, we are going to atribute a score of zero for each potential score below the 1st quartile (0.25)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.0)

active_ligand_target_links_100<-active_ligand_target_links

#In this case, we are going to order the ligands as well as the genes for the regulatory potential value & nr. of targets
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network_100 = vis_ligand_target %>% make_heatmap_ggplot("Prioritized Ligands","Transcriptional Program", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple") + theme(axis.text.x = element_text(face = "italic"),legend.key.width = unit(1,"cm"))
