#Distance Matrix of full Menu based on Hill numbers----

# metric
# A vector containing "C", "U", "V" or "S". C: Sørensen-type overlap or complement. U: Jaccard-type overlap or complement. V: Sørensen-type turnover or complement. S: Jaccard-type turnover or complement. See hilldiv wiki for further information.

# "C" The Sørensen-type overlap quantifies the effective average proportion of a sub‐systems OTUs (or lineages in the case of phylodiversities) that is shared across all subsystems. This is thus a metric that quantifies overlap from the subsystems perspective. Its corresponding dissimilarity measure (1 - CqN) quantifies the effective average proportion of nonshared OTUs or lineages in a system. CqN is integrated in the functions beta.dis() and pair.dis().

# "V" The Sørensen-type turnover-complement is the complement of the Sørensen‐type turnover, which quantifies the normalized OTU turnover rate with respect to the average subsystem (i.e., alpha), thus provides the proportion of a typical subsystem that changes across subsystems.

#"U" The Jaccard-type overlap quantifies the effective proportion of OTUs or lineages in a system that are shared across all subsystems. Hence, this metric quantifies overlap from the perspective of the overall system. Its corresponding dissimilarity (1 - UqN) quantifies the effective proportion of nonshared OTUs or lineages in the overall system.

# "S" The Jaccard-type turnover-complement is the complement of the Jaccard‐type turnover, which quantifies the normalized OTU turnover rate with respect to the whole system (i.e. gamma).


## prep taxonomic (taxa-merged) species-named eDNA index for hilldiv----
Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps<-readRDS("../10_Phyloseq/Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps.RDS")

#make otu table and ultrametric tree for Menu_taxmerge_eDNAindex ASVs(taxa)
Menu_taxmerge_eDNAindex_otudf<-cbind(data.frame(otu_table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))
#Make ultrametric tree for otu object
Menu_taxmerge_eDNAindex_tree<-Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@phy_tree
Menu_taxmerge_eDNAindex_ultra_tree<-chronos(Menu_taxmerge_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Menu_taxmerge_eDNAindex_ultra_tree, 'Menu_taxmerge_eDNAindex_ultra.tre')
Menu_taxmerge_eDNAindex_ultra_tree <- ape::read.tree('Menu_taxmerge_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
namechange<-rownames(Menu_taxmerge_eDNAindex_otudf)
namechange<-gsub(pattern = " ",replacement = "_",namechange)
rownames(Menu_taxmerge_eDNAindex_otudf)<-namechange

#make envvars for later
Menu_taxmerge_eDNAindex_envvars<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))

#Calculate distances for various Hill numbers----

##Phylogenetic Hill q=0 (Faith's PD) distance matrix with taxa merged object ----

Menu_taxmerged_Faiths_Hilldistmat<-pair_dis(countable = Menu_taxmerge_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1,tree = Menu_taxmerge_eDNAindex_ultra_tree)

##Phylogenetic Hill q=1 (Allen's H) distance matrix with taxa merged object ----

Menu_taxmerged_Allens_Hilldistmat<-pair_dis(countable = Menu_taxmerge_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1,tree = Menu_taxmerge_eDNAindex_ultra_tree)

##Non-Phylogenetic Hill q=0 (Richness) distance matrix with taxa merged object ----

Menu_taxmerged_Richness_Hilldistmat<-pair_dis(countable = Menu_taxmerge_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)


##Non-Phylogenetic Hill q=1 (exponent of Shannon's) distance matrix with taxa merged object ----

Menu_taxmerged_Shannons_Hilldistmat<-pair_dis(countable = Menu_taxmerge_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)

#Get Jaccard-type turnover complement (with respect to gamma) as a symmetrical distance matrix for each----
Menu_taxmerged_Faiths_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

Menu_taxmerged_Allens_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_Richness_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

Menu_taxmerged_Shannons_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#Get Sørensen type turnover complement (with respect to alpha) as a symmetrical distance matrix for each----
Menu_taxmerged_Faiths_Hill_Sor_dissim<-as.dist(Menu_taxmerged_Faiths_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_Allens_Hill_Sor_dissim<-as.dist(Menu_taxmerged_Allens_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_Richness_Hill_Sor_dissim<-as.dist(Menu_taxmerged_Richness_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_Shannons_Hill_Sor_dissim<-as.dist(Menu_taxmerged_Shannons_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)


##Do the same for Menu with hammerheads removed from the phyloseq object

#First remove hammerheads so they don't confound the calculation.
Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps <- Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps %>% subset_taxa(species != "Sphyrna lewini")

## prep taxonomic (taxa-merged) species-named eDNA index rm_hh for hilldiv----


Menu_taxmerge_rmhh_eDNAindex_otudf<-cbind(data.frame(otu_table(Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps)))
#Make ultrametric tree for otu object
Menu_taxmerge_rmhh_eDNAindex_tree<-Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps@phy_tree
Menu_taxmerge_rmhh_eDNAindex_ultra_tree<-chronos(Menu_taxmerge_rmhh_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Menu_taxmerge_rmhh_eDNAindex_ultra_tree, 'Menu_taxmerge_rmhh_eDNAindex_ultra.tre')
Menu_taxmerge_rmhh_eDNAindex_ultra_tree <- ape::read.tree('Menu_taxmerge_rmhh_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
namechange_rmhh<-rownames(Menu_taxmerge_rmhh_eDNAindex_otudf)
namechange_rmhh<-gsub(pattern = " ",replacement = "_",namechange_rmhh)
rownames(Menu_taxmerge_rmhh_eDNAindex_otudf)<-namechange_rmhh

#make sample table for later
HH_rm_eDNAindex_envvars<-as.data.frame(cbind(sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps)))

##Hill distance matrix with taxa merged rmmhh object ----
Menu_taxmerged_rmhh_Faiths_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rmhh_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1,tree = Menu_taxmerge_rmhh_eDNAindex_ultra_tree)

Menu_taxmerged_rmhh_Allens_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rmhh_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1,tree = Menu_taxmerge_rmhh_eDNAindex_ultra_tree)

Menu_taxmerged_rmhh_Richness_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rmhh_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)

Menu_taxmerged_rmhh_Shannons_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rmhh_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)


## Jaccard-type turnover-complement of hammerheads-removed menu object----
#Get Jaccard-type turnover complement (with respect to gamma) as a symmetrical distance matrix for each
Menu_taxmerged_rmhh_Faiths_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rmhh_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

Menu_taxmerged_rmhh_Allens_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rmhh_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rmhh_Richness_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rmhh_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

Menu_taxmerged_rmhh_Shannons_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rmhh_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#Get Sørensen type turnover complement (with respect to alpha) as a symmetrical distance matrix for each----
Menu_taxmerged_rmhh_Faiths_Hill_Sor_dissim<-as.dist(Menu_taxmerged_rmhh_Faiths_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rmhh_Allens_Hill_Sor_dissim<-as.dist(Menu_taxmerged_rmhh_Allens_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rmhh_Richness_Hill_Sor_dissim<-as.dist(Menu_taxmerged_rmhh_Richness_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rmhh_Shannons_Hill_Sor_dissim<-as.dist(Menu_taxmerged_rmhh_Shannons_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)




##Do the same for carcharhinus----


#First remove carcharhinus reads so they don't confound the calculation.
Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps <- Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps %>% subset_taxa(species != "Carcharhinus genus")

tax_table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)
## prep taxonomic (taxa-merged) species-named eDNA index rm_carch for hilldiv----


Menu_taxmerge_rm_carch_eDNAindex_otudf<-cbind(data.frame(otu_table(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps)))
#Make ultrametric tree for otu object
Menu_taxmerge_rm_carch_eDNAindex_tree<-Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps@phy_tree
Menu_taxmerge_rm_carch_eDNAindex_ultra_tree<-chronos(Menu_taxmerge_rm_carch_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Menu_taxmerge_rm_carch_eDNAindex_ultra_tree, 'Menu_taxmerge_rm_carch_eDNAindex_ultra.tre')
Menu_taxmerge_rm_carch_eDNAindex_ultra_tree <- ape::read.tree('Menu_taxmerge_rm_carch_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
namechange_rm_carch<-rownames(Menu_taxmerge_rm_carch_eDNAindex_otudf)
namechange_rm_carch<-gsub(pattern = " ",replacement = "_",namechange_rm_carch)
rownames(Menu_taxmerge_rm_carch_eDNAindex_otudf)<-namechange_rm_carch

#make sample table for later
rm_carch_eDNAindex_envvars<-as.data.frame(cbind(sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps)))

##Hill distance matrix with taxa merged rm_carch object for all three measures ----

Menu_taxmerged_rm_carch_Faiths_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rm_carch_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1,tree = Menu_taxmerge_rm_carch_eDNAindex_ultra_tree)

Menu_taxmerged_rm_carch_Allens_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rm_carch_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1,tree = Menu_taxmerge_rm_carch_eDNAindex_ultra_tree)

Menu_taxmerged_rm_carch_Richness_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rm_carch_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)

Menu_taxmerged_rm_carch_Shannons_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rm_carch_eDNAindex_otudf,qvalue =10,hierarchy=menu_sample_site_hierarchy,level = 1)

#get Jaccard dissimilarity matrix (turnover-complement with respect to gamma)

Menu_taxmerged_rm_carch_Faiths_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rm_carch_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rm_carch_Allens_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rm_carch_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rm_carch_Richness_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rm_carch_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rm_carch_Shannons_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rm_carch_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#get Sørensen dissimilarity matrix (turnover-complement with respect to alpha)

Menu_taxmerged_rm_carch_Faiths_Hill_Sor_dissim<-as.dist(Menu_taxmerged_rm_carch_Faiths_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rm_carch_Allens_Hill_Sor_dissim<-as.dist(Menu_taxmerged_rm_carch_Allens_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rm_carch_Richness_Hill_Sor_dissim<-as.dist(Menu_taxmerged_rm_carch_Richness_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_rm_carch_Shannons_Hill_Sor_dissim<-as.dist(Menu_taxmerged_rm_carch_Shannons_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
