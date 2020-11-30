#Intersect the list of the senders' genes & the ligands of NicheNet db
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,genes_sender)

#Intersect the expressed genes from receivers cells & the receptors/endpoints of NicheNet db
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,background_expressed_genes)

#To making it easier, we are going to work on the subset of the network that contains only the expressed ligands & expressed receptors
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)

#The db may contain a couple of NA's; this is better to be removed
lr_network_expressed = lr_network_expressed %>% filter(!is.na(from) & !is.na(to)) 

#All of the expressed ligands from senders are going to be considered potential ligands
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

#Now, we perform the prediction of the ligand activities concerning the gene set of interest
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

#Printing the df of the prediction of ligand activies, sorted in decreasing order using the Pearson Correlation Coef. (PCC)
ligand_activities %>% arrange(-pearson) 

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
