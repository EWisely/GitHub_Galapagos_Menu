#following Stier et al----
library(tidyverse)
library(vegan)
library(hilldiv)
library(phyloseq)
library(metagMisc)
library(ape)
library(metafor)
library(BiodiversityR)
library(ggplot2)
library(ggsci)

#load untransformed dataset----
Shark_Menu_taxa_merge_species_named_traits.ps<-readRDS("../10_Phyloseq/Menu_combined2/Menu_combined2_untransformed_traits_added_env_added_taxa_named_Phyloseq.rds")

tax_table(Shark_Menu_taxa_merge_species_named_traits.ps)

#remove each shark separately 
Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps<-subset_taxa(Shark_Menu_taxa_merge_species_named_traits.ps, species!="Sphyrna lewini")
Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps = prune_taxa(taxa_sums(Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps) > 0, Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps) 
Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps<-prune_samples(sample_sums(Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps)>0, Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps)

Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps<-subset_taxa(Shark_Menu_taxa_merge_species_named_traits.ps, species!="Carcharhinus genus")
Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps = prune_taxa(taxa_sums(Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps) > 0, Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps) 
Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps<-prune_samples(sample_sums(Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps)>0, Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps)

#Make Fish and Crustacean datasets for each shark
#Hammerhead Fish
Fish_taxa_merge_species_named_traits_rm_hh.ps<-subset_taxa(Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps,f_or_c=="Fish")
Fish_taxa_merge_species_named_traits_rm_hh.ps<-prune_taxa(taxa_sums(Fish_taxa_merge_species_named_traits_rm_hh.ps)>0,Fish_taxa_merge_species_named_traits_rm_hh.ps)
Fish_taxa_merge_species_named_traits_rm_hh.ps<-prune_samples(sample_sums(Fish_taxa_merge_species_named_traits_rm_hh.ps)>0, Fish_taxa_merge_species_named_traits_rm_hh.ps)

#Hammerhead Crustaceans
Crustacean_taxa_merge_species_named_traits_rm_hh.ps<-subset_taxa(Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps,f_or_c=="Crustacean")
Crustacean_taxa_merge_species_named_traits_rm_hh.ps<-prune_taxa(taxa_sums(Crustacean_taxa_merge_species_named_traits_rm_hh.ps)>0,Crustacean_taxa_merge_species_named_traits_rm_hh.ps)
Crustacean_taxa_merge_species_named_traits_rm_hh.ps<-prune_samples(sample_sums(Crustacean_taxa_merge_species_named_traits_rm_hh.ps)>0, Crustacean_taxa_merge_species_named_traits_rm_hh.ps)

#Carcharhinus Fish
Fish_taxa_merge_species_named_traits_rm_carch.ps<-subset_taxa(Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps,f_or_c=="Fish")
Fish_taxa_merge_species_named_traits_rm_carch.ps<-prune_taxa(taxa_sums(Fish_taxa_merge_species_named_traits_rm_carch.ps)>0,Fish_taxa_merge_species_named_traits_rm_carch.ps)
Fish_taxa_merge_species_named_traits_rm_carch.ps<-prune_samples(sample_sums(Fish_taxa_merge_species_named_traits_rm_carch.ps)>0, Fish_taxa_merge_species_named_traits_rm_carch.ps)


#Carcharhinus Crustaceans
Crustacean_taxa_merge_species_named_traits_rm_carch.ps<-subset_taxa(Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps,f_or_c=="Crustacean")
Crustacean_taxa_merge_species_named_traits_rm_carch.ps<-prune_taxa(taxa_sums(Crustacean_taxa_merge_species_named_traits_rm_carch.ps)>0,Crustacean_taxa_merge_species_named_traits_rm_carch.ps)
Crustacean_taxa_merge_species_named_traits_rm_carch.ps<-prune_samples(sample_sums(Crustacean_taxa_merge_species_named_traits_rm_carch.ps)>0, Crustacean_taxa_merge_species_named_traits_rm_carch.ps)

taxa_sums(Crustacean_taxa_merge_species_named_traits_rm_carch.ps)


##1000 rarefaction replicates to get every sample to the same coverage----
#Rarefy Hammerhead Menu, Fish, and Crustaceans, then Carcharhinus Menu Fish and Crustaceans

#This creates a new phyloseq object after averaging the 1000 rarefaction iterations
multi_rare_avg_rm_hh.ps<-phyloseq_mult_raref_avg(physeq = Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps,iter = 1000,verbose = TRUE, parallel = TRUE,SampSize = min(sample_sums(Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps)))

otu_table(Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps)
otu_table(multi_rare_avg_rm_hh.ps)

##Rarefy Hammerhead Fish
fish_multi_rare_avg_rm_hh.ps<-phyloseq_mult_raref_avg(physeq = Fish_taxa_merge_species_named_traits_rm_hh.ps,iter = 1000,verbose = TRUE, parallel = TRUE,SampSize = min(sample_sums(Fish_taxa_merge_species_named_traits_rm_hh.ps)))

otu_table(Fish_taxa_merge_species_named_traits_rm_hh.ps)
otu_table(fish_multi_rare_avg_rm_hh.ps)


##Rarefy Hammerhead Crustaceans
crustacean_multi_rare_avg_rm_hh.ps<-phyloseq_mult_raref_avg(physeq = Crustacean_taxa_merge_species_named_traits_rm_hh.ps,iter = 1000,verbose = TRUE, parallel = TRUE,SampSize = min(sample_sums(Crustacean_taxa_merge_species_named_traits_rm_hh.ps)))

otu_table(Crustacean_taxa_merge_species_named_traits_rm_hh.ps)
otu_table(crustacean_multi_rare_avg_rm_hh.ps)


##Rarefy Carcharhinus Fish
fish_multi_rare_avg_rm_carch.ps<-phyloseq_mult_raref_avg(physeq = Fish_taxa_merge_species_named_traits_rm_carch.ps,iter = 1000,verbose = TRUE, parallel = TRUE,SampSize = min(sample_sums(Fish_taxa_merge_species_named_traits_rm_carch.ps)))

otu_table(Fish_taxa_merge_species_named_traits_rm_carch.ps)
otu_table(fish_multi_rare_avg_rm_carch.ps)

##Rarefy Carcharhinus Crustaceans
crustacean_multi_rare_avg_rm_carch.ps<-phyloseq_mult_raref_avg(physeq = Crustacean_taxa_merge_species_named_traits_rm_carch.ps,iter = 1000,verbose = TRUE, parallel = TRUE,SampSize = min(sample_sums(Crustacean_taxa_merge_species_named_traits_rm_carch.ps)))

otu_table(Crustacean_taxa_merge_species_named_traits_rm_carch.ps)
otu_table(crustacean_multi_rare_avg_rm_carch.ps)

###eDNA index transform rarefied and unrarefied phyloseq objects----
#Hammerhead Menu unrarefied
Shark_Menu_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps<-phyloseq_standardize_otu_abundance(Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps, method = "wisconsin")
otu_table(Shark_Menu_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps)

#Hammerhead Menu rarefied
multi_rare_avg_rm_hh.eDNAindex.ps<-phyloseq_standardize_otu_abundance(multi_rare_avg_rm_hh.ps, method = "wisconsin")
otu_table(multi_rare_avg_rm_hh.eDNAindex.ps)

#Hammerhead Fish unrarefied
Fish_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps<-phyloseq_standardize_otu_abundance(Fish_taxa_merge_species_named_traits_rm_hh.ps, method = "wisconsin")
otu_table(Fish_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps)

#Hammerhead Fish rarefied
fish_multi_rare_avg_rm_hh.eDNAindex.ps<-phyloseq_standardize_otu_abundance(fish_multi_rare_avg_rm_hh.ps, method = "wisconsin")
otu_table(fish_multi_rare_avg_rm_hh.eDNAindex.ps)

#Hammerhead Crustaceans unrarefied
Crustacean_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps<-phyloseq_standardize_otu_abundance(Crustacean_taxa_merge_species_named_traits_rm_hh.ps, method = "wisconsin")
otu_table(Crustacean_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps)

#Hammerhead Crustaceans rarefied
crustacean_multi_rare_avg_rm_hh.eDNAindex.ps<-phyloseq_standardize_otu_abundance(crustacean_multi_rare_avg_rm_hh.ps, method = "wisconsin")
otu_table(crustacean_multi_rare_avg_rm_hh.eDNAindex.ps)


#Carcharhinus Menu unrarefied
Shark_Menu_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps<-phyloseq_standardize_otu_abundance(Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps, method = "wisconsin")
otu_table(Shark_Menu_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps)

#Carcharhinus Menu rarefied
multi_rare_avg_rm_carch.eDNAindex.ps<-phyloseq_standardize_otu_abundance(multi_rare_avg_rm_carch.ps, method = "wisconsin")
otu_table(multi_rare_avg_rm_carch.eDNAindex.ps)

#Carcharhinus Fish unrarefied
Fish_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps<-phyloseq_standardize_otu_abundance(Fish_taxa_merge_species_named_traits_rm_carch.ps, method = "wisconsin")
otu_table(Fish_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps)

#Carcharhinus Fish rarefied
fish_multi_rare_avg_rm_carch.eDNAindex.ps<-phyloseq_standardize_otu_abundance(fish_multi_rare_avg_rm_carch.ps, method = "wisconsin")
otu_table(fish_multi_rare_avg_rm_carch.eDNAindex.ps)

#Carcharhinus Crustaceans unrarefied
Crustacean_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps<-phyloseq_standardize_otu_abundance(Crustacean_taxa_merge_species_named_traits_rm_carch.ps, method = "wisconsin")
otu_table(Crustacean_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps)

#Carcharhinus Crustaceans rarefied
crustacean_multi_rare_avg_rm_carch.eDNAindex.ps<-phyloseq_standardize_otu_abundance(crustacean_multi_rare_avg_rm_carch.ps, method = "wisconsin")
otu_table(crustacean_multi_rare_avg_rm_carch.eDNAindex.ps)


##create input objects for hilldiv from unrarefied and rarefied eDNA indexed taxa-merged phyloseq----

#12 Versions: 2 predators (hammerheads and Carcharhinus), 3 groupings (Full Menu, Fish, Crustaceans), rarefied and unrarefied (2 rarefaction states)----

###1 Unrarefied Hammerhead Menu----

Menu_taxmerge_rm_hh_unraref_eDNAindex_otudf<-cbind(data.frame(otu_table(Shark_Menu_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps)))
#Make ultrametric tree for otu object 
Menu_taxmerge_rm_hh_unraref_eDNAindex_tree<-Shark_Menu_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps@phy_tree
Menu_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree<-chronos(Menu_taxmerge_rm_hh_unraref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Menu_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree, 'Menu_taxmerge_rm_hh_unraref_eDNAindex_ultra.tre')
Menu_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree <- ape::read.tree('Menu_taxmerge_rm_hh_unraref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
unrarefnamechange<-rownames(Menu_taxmerge_rm_hh_unraref_eDNAindex_otudf)
unrarefnamechange<-gsub(pattern = " ",replacement = "_",unrarefnamechange)
rownames(Menu_taxmerge_rm_hh_unraref_eDNAindex_otudf)<-unrarefnamechange

match_data(Menu_taxmerge_rm_hh_unraref_eDNAindex_otudf,Menu_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

#make envvars for later
Menu_taxmerge_rm_hh_unraref_eDNAindex_envvars<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps)))



###2 Rarefied Hammerhead Menu----
Menu_taxmerge_rm_hh_raref_eDNAindex_otudf<-cbind(data.frame(otu_table(multi_rare_avg_rm_hh.ps)))
#Make ultrametric tree for otu object 
Menu_taxmerge_rm_hh_raref_eDNAindex_tree<-multi_rare_avg_rm_hh.ps@phy_tree
Menu_taxmerge_rm_hh_raref_eDNAindex_ultra_tree<-chronos(Menu_taxmerge_rm_hh_raref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Menu_taxmerge_rm_hh_raref_eDNAindex_ultra_tree, 'Menu_taxmerge_rm_hh_raref_eDNAindex_ultra.tre')
Menu_taxmerge_rm_hh_raref_eDNAindex_ultra_tree <- ape::read.tree('Menu_taxmerge_rm_hh_raref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
rarefnamechange<-rownames(Menu_taxmerge_rm_hh_raref_eDNAindex_otudf)
rarefnamechange<-gsub(pattern = " ",replacement = "_",rarefnamechange)
rownames(Menu_taxmerge_rm_hh_raref_eDNAindex_otudf)<-rarefnamechange

match_data(Menu_taxmerge_rm_hh_raref_eDNAindex_otudf,Menu_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

#make envvars for later
Menu_taxmerge_rm_hh_raref_eDNAindex_envvars<-cbind(data.frame(sample_data(multi_rare_avg_rm_hh.ps)))

###3 Unrarefied Hammerhead Fish----

Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf<-cbind(data.frame(otu_table(Fish_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps)))
#Make ultrametric tree for otu object 
Fish_taxmerge_rm_hh_unraref_eDNAindex_tree<-Fish_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps@phy_tree
Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree<-chronos(Fish_taxmerge_rm_hh_unraref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree, 'Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra.tre')
Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree <- ape::read.tree('Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
unrarefnamechange<-rownames(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf)
unrarefnamechange<-gsub(pattern = " ",replacement = "_",unrarefnamechange)
rownames(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf)<-unrarefnamechange

match_data(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

#make envvars for later
Fish_taxmerge_rm_hh_unraref_eDNAindex_envvars<-cbind(data.frame(sample_data(Fish_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps)))

###4 Rarefied Hammerhead Fish----
Fish_taxmerge_rm_hh_raref_eDNAindex_otudf<-cbind(data.frame(otu_table(fish_multi_rare_avg_rm_hh.ps)))
#Make ultrametric tree for otu object 
Fish_taxmerge_rm_hh_raref_eDNAindex_tree<-fish_multi_rare_avg_rm_hh.ps@phy_tree
Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree<-chronos(Fish_taxmerge_rm_hh_raref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree, 'Fish_taxmerge_rm_hh_raref_eDNAindex_ultra.tre')
Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree <- ape::read.tree('Fish_taxmerge_rm_hh_raref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
rarefnamechange<-rownames(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf)
rarefnamechange<-gsub(pattern = " ",replacement = "_",rarefnamechange)
rownames(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf)<-rarefnamechange

match_data(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

#make envvars for later
Fish_taxmerge_rm_hh_raref_eDNAindex_envvars<-cbind(data.frame(sample_data(fish_multi_rare_avg_rm_hh.ps)))





###5 Unrarefied Hammerhead Crustaceans----

Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf<-cbind(data.frame(otu_table(Crustacean_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps)))
#Make ultrametric tree for otu object 
Crustacean_taxmerge_rm_hh_unraref_eDNAindex_tree<-Crustacean_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps@phy_tree
Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree<-chronos(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree, 'Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra.tre')
Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree <- ape::read.tree('Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
unrarefnamechange<-rownames(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf)
unrarefnamechange<-gsub(pattern = " ",replacement = "_",unrarefnamechange)
rownames(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf)<-unrarefnamechange

match_data(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

#make envvars for later
Crustacean_taxmerge_rm_hh_unraref_eDNAindex_envvars<-cbind(data.frame(sample_data(Crustacean_taxa_merge_species_named_traits_rm_hh_eDNAindex.ps)))

###6 Rarefied Hammerhead Crustaceans----
Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf<-cbind(data.frame(otu_table(crustacean_multi_rare_avg_rm_hh.ps)))
#Make ultrametric tree for otu object 
Crustacean_taxmerge_rm_hh_raref_eDNAindex_tree<-crustacean_multi_rare_avg_rm_hh.ps@phy_tree
Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree<-chronos(Crustacean_taxmerge_rm_hh_raref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree, 'Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra.tre')
Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree <- ape::read.tree('Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
rarefnamechange<-rownames(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf)
rarefnamechange<-gsub(pattern = " ",replacement = "_",rarefnamechange)
rownames(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf)<-rarefnamechange

match_data(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

#make envvars for later
Crustacean_taxmerge_rm_hh_raref_eDNAindex_envvars<-cbind(data.frame(sample_data(crustacean_multi_rare_avg_rm_hh.ps)))





###7 Unrarefied Carcharhinus Menu----

Menu_taxmerge_rm_carch_unraref_eDNAindex_otudf<-cbind(data.frame(otu_table(Shark_Menu_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps)))
#Make ultrametric tree for otu object 
Menu_taxmerge_rm_carch_unraref_eDNAindex_tree<-Shark_Menu_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps@phy_tree
Menu_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree<-chronos(Menu_taxmerge_rm_carch_unraref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Menu_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree, 'Menu_taxmerge_rm_carch_unraref_eDNAindex_ultra.tre')
Menu_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree <- ape::read.tree('Menu_taxmerge_rm_carch_unraref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
unrarefnamechange<-rownames(Menu_taxmerge_rm_carch_unraref_eDNAindex_otudf)
unrarefnamechange<-gsub(pattern = " ",replacement = "_",unrarefnamechange)
rownames(Menu_taxmerge_rm_carch_unraref_eDNAindex_otudf)<-unrarefnamechange

match_data(Menu_taxmerge_rm_carch_unraref_eDNAindex_otudf,Menu_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

#make envvars for later
Menu_taxmerge_rm_carch_unraref_eDNAindex_envvars<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps)))



###8 Rarefied Carcharhinus Menu----
Menu_taxmerge_rm_carch_raref_eDNAindex_otudf<-cbind(data.frame(otu_table(multi_rare_avg_rm_carch.ps)))
#Make ultrametric tree for otu object 
Menu_taxmerge_rm_carch_raref_eDNAindex_tree<-multi_rare_avg_rm_carch.ps@phy_tree
Menu_taxmerge_rm_carch_raref_eDNAindex_ultra_tree<-chronos(Menu_taxmerge_rm_carch_raref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Menu_taxmerge_rm_carch_raref_eDNAindex_ultra_tree, 'Menu_taxmerge_rm_carch_raref_eDNAindex_ultra.tre')
Menu_taxmerge_rm_carch_raref_eDNAindex_ultra_tree <- ape::read.tree('Menu_taxmerge_rm_carch_raref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
rarefnamechange<-rownames(Menu_taxmerge_rm_carch_raref_eDNAindex_otudf)
rarefnamechange<-gsub(pattern = " ",replacement = "_",rarefnamechange)
rownames(Menu_taxmerge_rm_carch_raref_eDNAindex_otudf)<-rarefnamechange

match_data(Menu_taxmerge_rm_carch_raref_eDNAindex_otudf,Menu_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

#make envvars for later
Menu_taxmerge_rm_carch_raref_eDNAindex_envvars<-cbind(data.frame(sample_data(multi_rare_avg_rm_carch.ps)))

###9 Unrarefied Carcharhinus Fish----

Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf<-cbind(data.frame(otu_table(Fish_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps)))
#Make ultrametric tree for otu object 
Fish_taxmerge_rm_carch_unraref_eDNAindex_tree<-Fish_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps@phy_tree
Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree<-chronos(Fish_taxmerge_rm_carch_unraref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree, 'Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra.tre')
Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree <- ape::read.tree('Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
unrarefnamechange<-rownames(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf)
unrarefnamechange<-gsub(pattern = " ",replacement = "_",unrarefnamechange)
rownames(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf)<-unrarefnamechange

match_data(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

#make envvars for later
Fish_taxmerge_rm_carch_unraref_eDNAindex_envvars<-cbind(data.frame(sample_data(Fish_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps)))

###10 Rarefied Carcharhinus Fish----
Fish_taxmerge_rm_carch_raref_eDNAindex_otudf<-cbind(data.frame(otu_table(fish_multi_rare_avg_rm_carch.ps)))
#Make ultrametric tree for otu object 
Fish_taxmerge_rm_carch_raref_eDNAindex_tree<-fish_multi_rare_avg_rm_carch.ps@phy_tree
Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree<-chronos(Fish_taxmerge_rm_carch_raref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree, 'Fish_taxmerge_rm_carch_raref_eDNAindex_ultra.tre')
Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree <- ape::read.tree('Fish_taxmerge_rm_carch_raref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
rarefnamechange<-rownames(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf)
rarefnamechange<-gsub(pattern = " ",replacement = "_",rarefnamechange)
rownames(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf)<-rarefnamechange

match_data(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

#make envvars for later
Fish_taxmerge_rm_carch_raref_eDNAindex_envvars<-cbind(data.frame(sample_data(fish_multi_rare_avg_rm_carch.ps)))




###11 Unrarefied Carcharhinus Crustaceans----

Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf<-cbind(data.frame(otu_table(Crustacean_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps)))
#Make ultrametric tree for otu object 
Crustacean_taxmerge_rm_carch_unraref_eDNAindex_tree<-Crustacean_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps@phy_tree
Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree<-chronos(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree, 'Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra.tre')
Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree <- ape::read.tree('Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
unrarefnamechange<-rownames(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf)
unrarefnamechange<-gsub(pattern = " ",replacement = "_",unrarefnamechange)
rownames(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf)<-unrarefnamechange

match_data(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

#make envvars for later
Crustacean_taxmerge_rm_carch_unraref_eDNAindex_envvars<-cbind(data.frame(sample_data(Crustacean_taxa_merge_species_named_traits_rm_carch_eDNAindex.ps)))

###12 Rarefied Carcharhinus Crustaceans----
Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf<-cbind(data.frame(otu_table(crustacean_multi_rare_avg_rm_carch.ps)))
#Make ultrametric tree for otu object 
Crustacean_taxmerge_rm_carch_raref_eDNAindex_tree<-crustacean_multi_rare_avg_rm_carch.ps@phy_tree
Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree<-chronos(Crustacean_taxmerge_rm_carch_raref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree, 'Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra.tre')
Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree <- ape::read.tree('Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
rarefnamechange<-rownames(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf)
rarefnamechange<-gsub(pattern = " ",replacement = "_",rarefnamechange)
rownames(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf)<-rarefnamechange

match_data(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

#make envvars for later
Crustacean_taxmerge_rm_carch_raref_eDNAindex_envvars<-cbind(data.frame(sample_data(crustacean_multi_rare_avg_rm_carch.ps)))





#Calculate Hill alpha diversity metrics ----
#16 Versions: Just the 8 Fish and Crustacean versions above done with Faith's and Allen's Hill Alpha diversity ----

###Menu: Get Hill diversity q=0 considering phylogenetic diversity for each taxa-merged sample----

#Hammerhead Menu Faith's PD----

##unrarefied----
match_data(Menu_taxmerge_rm_hh_unraref_eDNAindex_otudf,Menu_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

rm_hh_unraref_FaithsPD_taxa_sample<-hill_div(Menu_taxmerge_rm_hh_unraref_eDNAindex_otudf,0,Menu_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)
rm_hh_unraref_FaithsPD_taxa_sample

rm_hh_unraref_Faiths.df<-rownames_to_column(as.data.frame(rm_hh_unraref_FaithsPD_taxa_sample),var="Sample")

##rarefied----
match_data(Menu_taxmerge_rm_hh_raref_eDNAindex_otudf,Menu_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

rm_hh_raref_FaithsPD_taxa_sample<-hill_div(Menu_taxmerge_rm_hh_raref_eDNAindex_otudf,0,Menu_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
rm_hh_raref_FaithsPD_taxa_sample

rm_hh_raref_Faiths.df<-rownames_to_column(as.data.frame(rm_hh_raref_FaithsPD_taxa_sample),var="Sample")

#Carcharhinus Menu Faith's PD----
##unrarefied----
match_data(Menu_taxmerge_rm_carch_unraref_eDNAindex_otudf,Menu_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

rm_carch_unraref_FaithsPD_taxa_sample<-hill_div(Menu_taxmerge_rm_carch_unraref_eDNAindex_otudf,0,Menu_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)
rm_carch_unraref_FaithsPD_taxa_sample

rm_carch_unraref_Faiths.df<-rownames_to_column(as.data.frame(rm_carch_unraref_FaithsPD_taxa_sample),var="Sample")

##rarefied ----
match_data(Menu_taxmerge_rm_carch_raref_eDNAindex_otudf,Menu_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

rm_carch_raref_FaithsPD_taxa_sample<-hill_div(Menu_taxmerge_rm_carch_raref_eDNAindex_otudf,0,Menu_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)
rm_carch_raref_FaithsPD_taxa_sample

rm_carch_raref_Faiths.df<-rownames_to_column(as.data.frame(rm_carch_raref_FaithsPD_taxa_sample),var="Sample")

# #Hammerhead Menu Allen’s H using eDNAindex transformed taxa-merged read abundances----
# rm_hh_raref_AllensH_taxa_sample_eDNAindex<-hill_div(Menu_taxmerge_rm_hh_raref_eDNAindex_otudf,1,Menu_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
# rm_hh_raref_AllensH_taxa_sample_eDNAindex
# 
# rm_hh_raref_Allens.df<-rownames_to_column(as.data.frame(rm_hh_raref_AllensH_taxa_sample_eDNAindex),var="Sample")
# 

#1 Hammerhead Fish Unrarefied Faith's PD ----

match_data(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

rm_hh_unraref_fish_FaithsPD_taxa_sample<-hill_div(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,0,Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)
rm_hh_unraref_fish_FaithsPD_taxa_sample

rm_hh_unraref_fish_Faiths.df<-rownames_to_column(as.data.frame(rm_hh_unraref_fish_FaithsPD_taxa_sample),var="Sample")



#2 Hammerhead Fish Rarefied Faith's PD ----

match_data(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

rm_hh_raref_fish_FaithsPD_taxa_sample<-hill_div(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,0,Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
rm_hh_raref_fish_FaithsPD_taxa_sample

rm_hh_raref_fish_Faiths.df<-rownames_to_column(as.data.frame(rm_hh_raref_fish_FaithsPD_taxa_sample),var="Sample")


#1.5 Hammerhead Fish Unrarefied Shannon's exponent ----

match_data(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

rm_hh_unraref_fish_ShannonsExp_taxa_sample<-hill_div(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,1)#,Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)
rm_hh_unraref_fish_ShannonsExp_taxa_sample

rm_hh_unraref_fish_Shannons.df<-rownames_to_column(as.data.frame(rm_hh_unraref_fish_ShannonsExp_taxa_sample),var="Sample")



#2.5 Hammerhead Fish Rarefied Shannon's exponent ----

match_data(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

rm_hh_raref_fish_ShannonsExp_taxa_sample<-hill_div(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,1)#,Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
rm_hh_raref_fish_ShannonsExp_taxa_sample

rm_hh_raref_fish_Shannons.df<-rownames_to_column(as.data.frame(rm_hh_raref_fish_ShannonsExp_taxa_sample),var="Sample")


#1.8 Hammerhead Fish Unrarefied Richness ----

match_data(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

rm_hh_unraref_fish_Richness_taxa_sample<-hill_div(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,0)#,Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)
rm_hh_unraref_fish_Richness_taxa_sample

rm_hh_unraref_fish_Richness.df<-rownames_to_column(as.data.frame(rm_hh_unraref_fish_Richness_taxa_sample),var="Sample")



#2.8 Hammerhead Fish Rarefied Richness ----

match_data(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

rm_hh_raref_fish_Richness_taxa_sample<-hill_div(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,1)#,Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
rm_hh_raref_fish_Richness_taxa_sample

rm_hh_raref_fish_Richness.df<-rownames_to_column(as.data.frame(rm_hh_raref_fish_Richness_taxa_sample),var="Sample")





#3 Hammerhead Fish Unrarefied Allen's H ----

match_data(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

rm_hh_unraref_fish_AllensH_taxa_sample<-hill_div(Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,1,Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)
rm_hh_unraref_fish_AllensH_taxa_sample

rm_hh_unraref_fish_Allens.df<-rownames_to_column(as.data.frame(rm_hh_unraref_fish_AllensH_taxa_sample),var="Sample")



#4 Hammerhead Fish Rarefied Allen's H ----

match_data(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

rm_hh_raref_fish_AllensH_taxa_sample<-hill_div(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,1,Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
rm_hh_raref_fish_AllensH_taxa_sample

rm_hh_raref_fish_Allens.df<-rownames_to_column(as.data.frame(rm_hh_raref_fish_AllensH_taxa_sample),var="Sample")


#5 Hammerhead Crustacean Unrarefied Faith's PD ----

match_data(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

rm_hh_unraref_crustacean_FaithsPD_taxa_sample<-hill_div(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,0,Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)
rm_hh_unraref_crustacean_FaithsPD_taxa_sample

rm_hh_unraref_crustacean_Faiths.df<-rownames_to_column(as.data.frame(rm_hh_unraref_crustacean_FaithsPD_taxa_sample),var="Sample")



#6 Hammerhead Crustacean Rarefied Faith's PD ----

match_data(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

rm_hh_raref_crustacean_FaithsPD_taxa_sample<-hill_div(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,0,Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
rm_hh_raref_crustacean_FaithsPD_taxa_sample

rm_hh_raref_crustacean_Faiths.df<-rownames_to_column(as.data.frame(rm_hh_raref_crustacean_FaithsPD_taxa_sample),var="Sample")

#5.5 Hammerhead Crustacean Unrarefied Shannon's exponent ----

match_data(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

rm_hh_unraref_crustacean_ShannonsExp_taxa_sample<-hill_div(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,1)#,Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)
rm_hh_unraref_crustacean_ShannonsExp_taxa_sample

rm_hh_unraref_crustacean_Shannons.df<-rownames_to_column(as.data.frame(rm_hh_unraref_crustacean_ShannonsExp_taxa_sample),var="Sample")



#6.5 Hammerhead Crustacean Rarefied Shannon's exponent ----

match_data(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

rm_hh_raref_crustacean_ShannonsExp_taxa_sample<-hill_div(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,1)#,Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
rm_hh_raref_crustacean_ShannonsExp_taxa_sample

rm_hh_raref_crustacean_Shannons.df<-rownames_to_column(as.data.frame(rm_hh_raref_crustacean_ShannonsExp_taxa_sample),var="Sample")

#5.8 Hammerhead Crustacean Unrarefied Richness ----

match_data(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

rm_hh_unraref_crustacean_Richness_taxa_sample<-hill_div(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,0)#,Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)
rm_hh_unraref_crustacean_Richness_taxa_sample

rm_hh_unraref_crustacean_Richness.df<-rownames_to_column(as.data.frame(rm_hh_unraref_crustacean_Richness_taxa_sample),var="Sample")



#6.8 Hammerhead Crustacean Rarefied Richness ----

match_data(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

rm_hh_raref_crustacean_Richness_taxa_sample<-hill_div(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,0)#,Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
rm_hh_raref_crustacean_Richness_taxa_sample

rm_hh_raref_crustacean_Richness.df<-rownames_to_column(as.data.frame(rm_hh_raref_crustacean_Richness_taxa_sample),var="Sample")




#7 Hammerhead Crustacean Unrarefied Allen's H ----

match_data(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

rm_hh_unraref_crustacean_AllensH_taxa_sample<-hill_div(Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,1,Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)
rm_hh_unraref_crustacean_AllensH_taxa_sample

rm_hh_unraref_crustacean_Allens.df<-rownames_to_column(as.data.frame(rm_hh_unraref_crustacean_AllensH_taxa_sample),var="Sample")



#8 Hammerhead Crustacean Rarefied Allen's H ----

match_data(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

rm_hh_raref_crustacean_AllensH_taxa_sample<-hill_div(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,1,Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)
rm_hh_raref_crustacean_AllensH_taxa_sample

rm_hh_raref_crustacean_Allens.df<-rownames_to_column(as.data.frame(rm_hh_raref_crustacean_AllensH_taxa_sample),var="Sample")


##9 Carcharhinus Fish Unrarefied Faith's PD ----

match_data(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

rm_carch_unraref_fish_FaithsPD_taxa_sample<-hill_div(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,0,Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)
rm_carch_unraref_fish_FaithsPD_taxa_sample

rm_carch_unraref_fish_Faiths.df<-rownames_to_column(as.data.frame(rm_carch_unraref_fish_FaithsPD_taxa_sample),var="Sample")



#10 Carcharhinus Fish Rarefied Faith's PD ----

match_data(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

rm_carch_raref_fish_FaithsPD_taxa_sample<-hill_div(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,0,Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)
rm_carch_raref_fish_FaithsPD_taxa_sample

rm_carch_raref_fish_Faiths.df<-rownames_to_column(as.data.frame(rm_carch_raref_fish_FaithsPD_taxa_sample),var="Sample")

##9.5 Carcharhinus Fish Unrarefied Shannon's exponent ----

match_data(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

rm_carch_unraref_fish_ShannonsExp_taxa_sample<-hill_div(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,1)#,Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)
rm_carch_unraref_fish_ShannonsExp_taxa_sample

rm_carch_unraref_fish_Shannons.df<-rownames_to_column(as.data.frame(rm_carch_unraref_fish_ShannonsExp_taxa_sample),var="Sample")



#10.5 Carcharhinus Fish Rarefied Shannon's exponent ----

match_data(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

rm_carch_raref_fish_ShannonsExp_taxa_sample<-hill_div(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,1)#,Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)
rm_carch_raref_fish_ShannonsExp_taxa_sample

rm_carch_raref_fish_Shannons.df<-rownames_to_column(as.data.frame(rm_carch_raref_fish_ShannonsExp_taxa_sample),var="Sample")


##9.8 Carcharhinus Fish Unrarefied Richness ----

match_data(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

rm_carch_unraref_fish_Richness_taxa_sample<-hill_div(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,1)#,Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)
rm_carch_unraref_fish_Richness_taxa_sample

rm_carch_unraref_fish_Richness.df<-rownames_to_column(as.data.frame(rm_carch_unraref_fish_Richness_taxa_sample),var="Sample")



#10.8 Carcharhinus Fish Rarefied Richness ----

match_data(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

rm_carch_raref_fish_Richness_taxa_sample<-hill_div(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,1)#,Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)
rm_carch_raref_fish_Richness_taxa_sample

rm_carch_raref_fish_Richness.df<-rownames_to_column(as.data.frame(rm_carch_raref_fish_Richness_taxa_sample),var="Sample")





#11 Carcharhinus Fish Unrarefied Allen's H ----

match_data(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

rm_carch_unraref_fish_AllensH_taxa_sample<-hill_div(Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,1,Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)
rm_carch_unraref_fish_AllensH_taxa_sample

rm_carch_unraref_fish_Allens.df<-rownames_to_column(as.data.frame(rm_carch_unraref_fish_AllensH_taxa_sample),var="Sample")



#12 Carcharhinus Fish Rarefied Allen's H ----

match_data(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

rm_carch_raref_fish_AllensH_taxa_sample<-hill_div(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,1,Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)
rm_carch_raref_fish_AllensH_taxa_sample

rm_carch_raref_fish_Allens.df<-rownames_to_column(as.data.frame(rm_carch_raref_fish_AllensH_taxa_sample),var="Sample")


#13 Carcharhinus Crustacean Unrarefied Faith's PD ----

match_data(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

rm_carch_unraref_crustacean_FaithsPD_taxa_sample<-hill_div(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,0,Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)
rm_carch_unraref_crustacean_FaithsPD_taxa_sample

rm_carch_unraref_crustacean_Faiths.df<-rownames_to_column(as.data.frame(rm_carch_unraref_crustacean_FaithsPD_taxa_sample),var="Sample")



#14 Carcharhinus Crustacean Rarefied Faith's PD ----

match_data(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

rm_carch_raref_crustacean_FaithsPD_taxa_sample<-hill_div(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,0,Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)
rm_carch_raref_crustacean_FaithsPD_taxa_sample

rm_carch_raref_crustacean_Faiths.df<-rownames_to_column(as.data.frame(rm_carch_raref_crustacean_FaithsPD_taxa_sample),var="Sample")


#13.5 Carcharhinus Crustacean Unrarefied Shannon's exponent ----

match_data(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

rm_carch_unraref_crustacean_ShannonsExp_taxa_sample<-hill_div(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,1)#,Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)
rm_carch_unraref_crustacean_ShannonsExp_taxa_sample

rm_carch_unraref_crustacean_Shannons.df<-rownames_to_column(as.data.frame(rm_carch_unraref_crustacean_ShannonsExp_taxa_sample),var="Sample")



#14.5 Carcharhinus Crustacean Rarefied Shannon's exponent ----

match_data(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

rm_carch_raref_crustacean_ShannonsExp_taxa_sample<-hill_div(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,1)#,Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)
rm_carch_raref_crustacean_ShannonsExp_taxa_sample

rm_carch_raref_crustacean_Shannons.df<-rownames_to_column(as.data.frame(rm_carch_raref_crustacean_ShannonsExp_taxa_sample),var="Sample")


#13.8 Carcharhinus Crustacean Unrarefied Richness ----

match_data(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

rm_carch_unraref_crustacean_Richness_taxa_sample<-hill_div(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,0)#,Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)
rm_carch_unraref_crustacean_Richness_taxa_sample

rm_carch_unraref_crustacean_Richness.df<-rownames_to_column(as.data.frame(rm_carch_unraref_crustacean_Richness_taxa_sample),var="Sample")



#14.8 Carcharhinus Crustacean Rarefied Richness ----

match_data(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

rm_carch_raref_crustacean_Richness_taxa_sample<-hill_div(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,0)#,Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)
rm_carch_raref_crustacean_Richness_taxa_sample

rm_carch_raref_crustacean_Richness.df<-rownames_to_column(as.data.frame(rm_carch_raref_crustacean_Richness_taxa_sample),var="Sample")


#15 Carcharhinus Crustacean Unrarefied Allen's H ----

match_data(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

rm_carch_unraref_crustacean_AllensH_taxa_sample<-hill_div(Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,1,Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)
rm_carch_unraref_crustacean_AllensH_taxa_sample

rm_carch_unraref_crustacean_Allens.df<-rownames_to_column(as.data.frame(rm_carch_unraref_crustacean_AllensH_taxa_sample),var="Sample")



#16 Carcharhinus Crustacean Rarefied Allen's H ----

match_data(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

rm_carch_raref_crustacean_AllensH_taxa_sample<-hill_div(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,1,Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)
rm_carch_raref_crustacean_AllensH_taxa_sample

rm_carch_raref_crustacean_Allens.df<-rownames_to_column(as.data.frame(rm_carch_raref_crustacean_AllensH_taxa_sample),var="Sample")


#Calculate Beta Diversity (Jaccard-style following Stier et al.)----
#Same 16 as the Alpha Diversity above Faith's and Allens

##Hammerhead Unrarefied Menu Faith's PD----

Menu_taxmerged_rm_hh_unraref_Faiths_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rm_hh_unraref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Menu_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

Menu_taxmerged_rm_hh_unraref_Faiths_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rm_hh_unraref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Menu_rm_hh_unraref_betadisp<-betadisper(d = Menu_taxmerged_rm_hh_unraref_Faiths_Hill_Jaccard_dissim,group = Menu_taxmerge_rm_hh_unraref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Menu_rm_hh_unraref_betadisp)
vegan::permutest(Faiths_Menu_rm_hh_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Faiths_Menu_rm_hh_unraref_betadisp)

Faiths_Menu_Beta_rm_hh_unraref<-rownames_to_column(as.data.frame(Faiths_Menu_rm_hh_unraref_betadisp$distances),var="Sample")

#Hammerhead Rarefied Menu Faiths PD----

Menu_taxmerged_rm_hh_raref_Faiths_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rm_hh_raref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Menu_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

Menu_taxmerged_rm_hh_raref_Faiths_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rm_hh_raref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Menu_rm_hh_raref_betadisp<-betadisper(d = Menu_taxmerged_rm_hh_raref_Faiths_Hill_Jaccard_dissim,group = Menu_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Menu_rm_hh_raref_betadisp)
vegan::permutest(Faiths_Menu_rm_hh_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Faiths_Menu_rm_hh_raref_betadisp)

Faiths_Menu_Beta_rm_hh_raref<-rownames_to_column(as.data.frame(Faiths_Menu_rm_hh_raref_betadisp$distances),var="Sample")






##1 Hammerhead Unrarefied Fish Faith's PD----

Fish_taxmerged_rm_hh_unraref_Faiths_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_hh_unraref_Faiths_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_hh_unraref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Fish_rm_hh_unraref_betadisp<-betadisper(d = Fish_taxmerged_rm_hh_unraref_Faiths_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_hh_unraref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Fish_rm_hh_unraref_betadisp)
vegan::permutest(Faiths_Fish_rm_hh_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Faiths_Fish_rm_hh_unraref_betadisp)

Faiths_Fish_Beta_rm_hh_unraref<-rownames_to_column(as.data.frame(Faiths_Fish_rm_hh_unraref_betadisp$distances),var="Sample")

##2 Hammerhead Rarefied Fish Faith's PD----

Fish_taxmerged_rm_hh_raref_Faiths_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_hh_raref_Faiths_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_hh_raref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Fish_rm_hh_raref_betadisp<-betadisper(d = Fish_taxmerged_rm_hh_raref_Faiths_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Fish_rm_hh_raref_betadisp)
vegan::permutest(Faiths_Fish_rm_hh_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Faiths_Fish_rm_hh_raref_betadisp)

Faiths_Fish_Beta_rm_hh_raref<-rownames_to_column(as.data.frame(Faiths_Fish_rm_hh_raref_betadisp$distances),var="Sample")

##1.5 Hammerhead Unrarefied Fish Shannon's exponent----

Fish_taxmerged_rm_hh_unraref_Shannons_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_hh_unraref_Shannons_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_hh_unraref_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Shannons_Fish_rm_hh_unraref_betadisp<-betadisper(d = Fish_taxmerged_rm_hh_unraref_Shannons_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_hh_unraref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Shannons_Fish_rm_hh_unraref_betadisp)
vegan::permutest(Shannons_Fish_rm_hh_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Shannons_Fish_rm_hh_unraref_betadisp)

Shannons_Fish_Beta_rm_hh_unraref<-rownames_to_column(as.data.frame(Shannons_Fish_rm_hh_unraref_betadisp$distances),var="Sample")

##2.5 Hammerhead Rarefied Fish Shannon's exponent----

Fish_taxmerged_rm_hh_raref_Shannons_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_hh_raref_Shannons_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_hh_raref_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Shannons_Fish_rm_hh_raref_betadisp<-betadisper(d = Fish_taxmerged_rm_hh_raref_Shannons_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Shannons_Fish_rm_hh_raref_betadisp)
vegan::permutest(Shannons_Fish_rm_hh_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Shannons_Fish_rm_hh_raref_betadisp)

Shannons_Fish_Beta_rm_hh_raref<-rownames_to_column(as.data.frame(Shannons_Fish_rm_hh_raref_betadisp$distances),var="Sample")




##1.8 Hammerhead Unrarefied Fish Richness----

Fish_taxmerged_rm_hh_unraref_Richness_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_hh_unraref_Richness_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_hh_unraref_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Richness_Fish_rm_hh_unraref_betadisp<-betadisper(d = Fish_taxmerged_rm_hh_unraref_Richness_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_hh_unraref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Richness_Fish_rm_hh_unraref_betadisp)
vegan::permutest(Richness_Fish_rm_hh_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Richness_Fish_rm_hh_unraref_betadisp)

Richness_Fish_Beta_rm_hh_unraref<-rownames_to_column(as.data.frame(Richness_Fish_rm_hh_unraref_betadisp$distances),var="Sample")

##2.8 Hammerhead Rarefied Fish Richness----

Fish_taxmerged_rm_hh_raref_Richness_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_hh_raref_Richness_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_hh_raref_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Richness_Fish_rm_hh_raref_betadisp<-betadisper(d = Fish_taxmerged_rm_hh_raref_Richness_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Richness_Fish_rm_hh_raref_betadisp)
vegan::permutest(Richness_Fish_rm_hh_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Richness_Fish_rm_hh_raref_betadisp)

Richness_Fish_Beta_rm_hh_raref<-rownames_to_column(as.data.frame(Richness_Fish_rm_hh_raref_betadisp$distances),var="Sample")


##3 Allen's Unrarefied Fish Hammerhead distance matrix)----

Fish_taxmerged_rm_hh_unraref_Allens_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_hh_unraref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Fish_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_hh_unraref_Allens_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_hh_unraref_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Allens_Fish_rm_hh_unraref_betadisp<-betadisper(d = Fish_taxmerged_rm_hh_unraref_Allens_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_hh_unraref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Allens_Fish_rm_hh_unraref_betadisp)
vegan::permutest(Allens_Fish_rm_hh_unraref_betadisp)

plot(Allens_Fish_rm_hh_unraref_betadisp)

Allens_Fish_Beta_rm_hh_unraref<-rownames_to_column(as.data.frame(Allens_Fish_rm_hh_unraref_betadisp$distances),var="Sample")

##4 Allen's Rarefied Fish Hammerhead distance matrix)----

Fish_taxmerged_rm_hh_raref_Allens_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_hh_raref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Fish_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_hh_raref_Allens_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_hh_raref_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Allens_Fish_rm_hh_raref_betadisp<-betadisper(d = Fish_taxmerged_rm_hh_raref_Allens_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Allens_Fish_rm_hh_raref_betadisp)
vegan::permutest(Allens_Fish_rm_hh_raref_betadisp)
plot(Allens_Fish_rm_hh_raref_betadisp)

Allens_Fish_Beta_rm_hh_raref<-rownames_to_column(as.data.frame(Allens_Fish_rm_hh_raref_betadisp$distances),var="Sample")

##5 Hammerhead Unrarefied Crustacean Faith's PD----

Crustacean_taxmerged_rm_hh_unraref_Faiths_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_hh_unraref_Faiths_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_hh_unraref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Crustacean_rm_hh_unraref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_hh_unraref_Faiths_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Crustacean_rm_hh_unraref_betadisp)
vegan::permutest(Faiths_Crustacean_rm_hh_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Faiths_Crustacean_rm_hh_unraref_betadisp)

Faiths_Crustacean_Beta_rm_hh_unraref<-rownames_to_column(as.data.frame(Faiths_Crustacean_rm_hh_unraref_betadisp$distances),var="Sample")

##6 Hammerhead Rarefied Crustacean Faith's PD----

Crustacean_taxmerged_rm_hh_raref_Faiths_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_hh_raref_Faiths_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_hh_raref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Crustacean_rm_hh_raref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_hh_raref_Faiths_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Crustacean_rm_hh_raref_betadisp)
vegan::permutest(Faiths_Crustacean_rm_hh_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Faiths_Crustacean_rm_hh_raref_betadisp)

Faiths_Crustacean_Beta_rm_hh_raref<-rownames_to_column(as.data.frame(Faiths_Crustacean_rm_hh_raref_betadisp$distances),var="Sample")

##5.2 Hammerhead Unrarefied Crustacean Richness----

Crustacean_taxmerged_rm_hh_unraref_Richness_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_hh_unraref_Richness_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_hh_unraref_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Richness_Crustacean_rm_hh_unraref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_hh_unraref_Richness_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Richness_Crustacean_rm_hh_unraref_betadisp)
vegan::permutest(Richness_Crustacean_rm_hh_unraref_betadisp)

plot(Richness_Crustacean_rm_hh_unraref_betadisp)

Richness_Crustacean_Beta_rm_hh_unraref<-rownames_to_column(as.data.frame(Richness_Crustacean_rm_hh_unraref_betadisp$distances),var="Sample")

##6.2 Hammerhead Rarefied Crustacean Richness----

Crustacean_taxmerged_rm_hh_raref_Richness_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_hh_raref_Richness_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_hh_raref_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Richness_Crustacean_rm_hh_raref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_hh_raref_Richness_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Richness_Crustacean_rm_hh_raref_betadisp)
vegan::permutest(Richness_Crustacean_rm_hh_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Richness_Crustacean_rm_hh_raref_betadisp)

Richness_Crustacean_Beta_rm_hh_raref<-rownames_to_column(as.data.frame(Richness_Crustacean_rm_hh_raref_betadisp$distances),var="Sample")

##5.5 Hammerhead Unrarefied Crustacean Shannon's exponent----

Crustacean_taxmerged_rm_hh_unraref_Shannons_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_hh_unraref_Shannons_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_hh_unraref_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Shannons_Crustacean_rm_hh_unraref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_hh_unraref_Shannons_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Shannons_Crustacean_rm_hh_unraref_betadisp)
vegan::permutest(Shannons_Crustacean_rm_hh_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Shannons_Crustacean_rm_hh_unraref_betadisp)

Shannons_Crustacean_Beta_rm_hh_unraref<-rownames_to_column(as.data.frame(Shannons_Crustacean_rm_hh_unraref_betadisp$distances),var="Sample")

##6.5 Hammerhead Rarefied Crustacean Shannon's exponent----

Crustacean_taxmerged_rm_hh_raref_Shannons_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_hh_raref_Shannons_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_hh_raref_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Shannons_Crustacean_rm_hh_raref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_hh_raref_Shannons_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Shannons_Crustacean_rm_hh_raref_betadisp)
vegan::permutest(Shannons_Crustacean_rm_hh_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Shannons_Crustacean_rm_hh_raref_betadisp)

Shannons_Crustacean_Beta_rm_hh_raref<-rownames_to_column(as.data.frame(Shannons_Crustacean_rm_hh_raref_betadisp$distances),var="Sample")






##7 Allen's Unrarefied Crustacean Hammerhead distance matrix)----

Crustacean_taxmerged_rm_hh_unraref_Allens_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_hh_unraref_Allens_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_hh_unraref_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Allens_Crustacean_rm_hh_unraref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_hh_unraref_Allens_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_hh_unraref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Allens_Crustacean_rm_hh_unraref_betadisp)
vegan::permutest(Allens_Crustacean_rm_hh_unraref_betadisp)
plot(Allens_Crustacean_rm_hh_unraref_betadisp)

Allens_Crustacean_Beta_rm_hh_unraref<-rownames_to_column(as.data.frame(Allens_Crustacean_rm_hh_unraref_betadisp$distances),var="Sample")

##8 Allen's Rarefied Crustacean Hammerhead distance matrix)----

Crustacean_taxmerged_rm_hh_raref_Allens_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Crustacean_taxmerge_rm_hh_raref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_hh_raref_Allens_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_hh_raref_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Allens_Crustacean_rm_hh_raref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_hh_raref_Allens_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Allens_Crustacean_rm_hh_raref_betadisp)
vegan::permutest(Allens_Crustacean_rm_hh_raref_betadisp)
plot(Allens_Crustacean_rm_hh_raref_betadisp)

Allens_Crustacean_Beta_rm_hh_raref<-rownames_to_column(as.data.frame(Allens_Crustacean_rm_hh_raref_betadisp$distances),var="Sample")

##9 Carcharhinus Unrarefied Fish Faith's PD----

Fish_taxmerged_rm_carch_unraref_Faiths_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_carch_unraref_Faiths_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_carch_unraref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Fish_rm_carch_unraref_betadisp<-betadisper(d = Fish_taxmerged_rm_carch_unraref_Faiths_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_carch_unraref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Fish_rm_carch_unraref_betadisp)
vegan::permutest(Faiths_Fish_rm_carch_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Faiths_Fish_rm_carch_unraref_betadisp)

Faiths_Fish_Beta_rm_carch_unraref<-rownames_to_column(as.data.frame(Faiths_Fish_rm_carch_unraref_betadisp$distances),var="Sample")

##10 Carcharhinus Rarefied Fish Faith's PD----

Fish_taxmerged_rm_carch_raref_Faiths_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_carch_raref_Faiths_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_carch_raref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Fish_rm_carch_raref_betadisp<-betadisper(d = Fish_taxmerged_rm_carch_raref_Faiths_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Fish_rm_carch_raref_betadisp)
vegan::permutest(Faiths_Fish_rm_carch_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Faiths_Fish_rm_carch_raref_betadisp)

Faiths_Fish_Beta_rm_carch_raref<-rownames_to_column(as.data.frame(Faiths_Fish_rm_carch_raref_betadisp$distances),var="Sample")

##11 Allen's Unrarefied Fish Carcharhinus distance matrix)----

Fish_taxmerged_rm_carch_unraref_Allens_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_carch_unraref_Allens_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_carch_unraref_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)


#calculate dispersion around the centroid

Allens_Fish_rm_carch_unraref_betadisp<-betadisper(d = Fish_taxmerged_rm_carch_unraref_Allens_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_carch_unraref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Allens_Fish_rm_carch_unraref_betadisp)
vegan::permutest(Allens_Fish_rm_carch_unraref_betadisp)
plot(Allens_Fish_rm_carch_unraref_betadisp)

Allens_Fish_Beta_rm_carch_unraref<-rownames_to_column(as.data.frame(Allens_Fish_rm_carch_unraref_betadisp$distances),var="Sample")

##12 Allen's Rarefied Fish Carcharhinus distance matrix)----

Fish_taxmerged_rm_carch_raref_Allens_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_carch_raref_Allens_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_carch_raref_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Allens_Fish_rm_carch_raref_betadisp<-betadisper(d = Fish_taxmerged_rm_carch_raref_Allens_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Allens_Fish_rm_carch_raref_betadisp)
vegan::permutest(Allens_Fish_rm_carch_raref_betadisp)
plot(Allens_Fish_rm_carch_raref_betadisp)

Allens_Fish_Beta_rm_carch_raref<-rownames_to_column(as.data.frame(Allens_Fish_rm_carch_raref_betadisp$distances),var="Sample")

##13 Carcharhinus Unrarefied Crustacean Faith's PD----

Crustacean_taxmerged_rm_carch_unraref_Faiths_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_carch_unraref_Faiths_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_carch_unraref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Crustacean_rm_carch_unraref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_carch_unraref_Faiths_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Crustacean_rm_carch_unraref_betadisp)
vegan::permutest(Faiths_Crustacean_rm_carch_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Faiths_Crustacean_rm_carch_unraref_betadisp)

Faiths_Crustacean_Beta_rm_carch_unraref<-rownames_to_column(as.data.frame(Faiths_Crustacean_rm_carch_unraref_betadisp$distances),var="Sample")

##14 Carcharhinus Rarefied Crustacean Faith's PD----

Crustacean_taxmerged_rm_carch_raref_Faiths_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_carch_raref_Faiths_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_carch_raref_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Faiths_Crustacean_rm_carch_raref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_carch_raref_Faiths_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faiths_Crustacean_rm_carch_raref_betadisp)
vegan::permutest(Faiths_Crustacean_rm_carch_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Faiths_Crustacean_rm_carch_raref_betadisp)

Faiths_Crustacean_Beta_rm_carch_raref<-rownames_to_column(as.data.frame(Faiths_Crustacean_rm_carch_raref_betadisp$distances),var="Sample")

##15 Allen's Unrarefied Crustacean Carcharhinus distance matrix)----

Crustacean_taxmerged_rm_carch_unraref_Allens_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_carch_unraref_Allens_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_carch_unraref_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Allens_Crustacean_rm_carch_unraref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_carch_unraref_Allens_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Allens_Crustacean_rm_carch_unraref_betadisp)
vegan::permutest(Allens_Crustacean_rm_carch_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object
plot(Allens_Crustacean_rm_carch_unraref_betadisp)

Allens_Crustacean_Beta_rm_carch_unraref<-rownames_to_column(as.data.frame(Allens_Crustacean_rm_carch_unraref_betadisp$distances),var="Sample")

##16 Allen's Rarefied Crustacean Carcharhinus distance matrix)----

Crustacean_taxmerged_rm_carch_raref_Allens_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1, tree = Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_carch_raref_Allens_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_carch_raref_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Allens_Crustacean_rm_carch_raref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_carch_raref_Allens_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Allens_Crustacean_rm_carch_raref_betadisp)
vegan::permutest(Allens_Crustacean_rm_carch_raref_betadisp)
plot(Allens_Crustacean_rm_carch_raref_betadisp)

Allens_Crustacean_Beta_rm_carch_raref<-rownames_to_column(as.data.frame(Allens_Crustacean_rm_carch_raref_betadisp$distances),var="Sample")







##17 Carcharhinus Unrarefied Crustacean Richness----

Crustacean_taxmerged_rm_carch_unraref_Richness_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_carch_unraref_Richness_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_carch_unraref_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Richness_Crustacean_rm_carch_unraref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_carch_unraref_Richness_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Richness_Crustacean_rm_carch_unraref_betadisp)
vegan::permutest(Richness_Crustacean_rm_carch_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Richness_Crustacean_rm_carch_unraref_betadisp)

Richness_Crustacean_Beta_rm_carch_unraref<-rownames_to_column(as.data.frame(Richness_Crustacean_rm_carch_unraref_betadisp$distances),var="Sample")

##18 Carcharhinus Rarefied Crustacean Richness----

Crustacean_taxmerged_rm_carch_raref_Richness_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_carch_raref_Richness_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_carch_raref_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Richness_Crustacean_rm_carch_raref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_carch_raref_Richness_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Richness_Crustacean_rm_carch_raref_betadisp)
vegan::permutest(Richness_Crustacean_rm_carch_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Richness_Crustacean_rm_carch_raref_betadisp)

Richness_Crustacean_Beta_rm_carch_raref<-rownames_to_column(as.data.frame(Richness_Crustacean_rm_carch_raref_betadisp$distances),var="Sample")

##19 Shannon's Unrarefied Crustacean Carcharhinus distance matrix)----

Crustacean_taxmerged_rm_carch_unraref_Shannons_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_carch_unraref_Shannons_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_carch_unraref_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Shannons_Crustacean_rm_carch_unraref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_carch_unraref_Shannons_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_carch_unraref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Shannons_Crustacean_rm_carch_unraref_betadisp)
vegan::permutest(Shannons_Crustacean_rm_carch_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object
plot(Shannons_Crustacean_rm_carch_unraref_betadisp)

Shannons_Crustacean_Beta_rm_carch_unraref<-rownames_to_column(as.data.frame(Shannons_Crustacean_rm_carch_unraref_betadisp$distances),var="Sample")

##20 Shannon's Rarefied Crustacean Carcharhinus distance matrix)----

Crustacean_taxmerged_rm_carch_raref_Shannons_Hilldistmat<-pair_dis(countable = Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Crustacean_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

Crustacean_taxmerged_rm_carch_raref_Shannons_Hill_Jaccard_dissim<-as.dist(Crustacean_taxmerged_rm_carch_raref_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Shannons_Crustacean_rm_carch_raref_betadisp<-betadisper(d = Crustacean_taxmerged_rm_carch_raref_Shannons_Hill_Jaccard_dissim,group = Crustacean_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Shannons_Crustacean_rm_carch_raref_betadisp)
vegan::permutest(Shannons_Crustacean_rm_carch_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object
plot(Shannons_Crustacean_rm_carch_raref_betadisp)

Shannons_Crustacean_Beta_rm_carch_raref<-rownames_to_column(as.data.frame(Shannons_Crustacean_rm_carch_raref_betadisp$distances),var="Sample")

##21 Carcharhinus Unrarefied Fish Richness----

Fish_taxmerged_rm_carch_unraref_Richness_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_carch_unraref_Richness_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_carch_unraref_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Richness_Fish_rm_carch_unraref_betadisp<-betadisper(d = Fish_taxmerged_rm_carch_unraref_Richness_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_carch_unraref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Richness_Fish_rm_carch_unraref_betadisp)
vegan::permutest(Richness_Fish_rm_carch_unraref_betadisp)
#not significantly different dispersions for hh and non-hh with the unrarefied object.
plot(Richness_Fish_rm_carch_unraref_betadisp)

Richness_Fish_Beta_rm_carch_unraref<-rownames_to_column(as.data.frame(Richness_Fish_rm_carch_unraref_betadisp$distances),var="Sample")

##22 Carcharhinus Rarefied Fish Richness----

Fish_taxmerged_rm_carch_raref_Richness_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_carch_raref_Richness_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_carch_raref_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)
#calculate dispersion around the centroid

Richness_Fish_rm_carch_raref_betadisp<-betadisper(d = Fish_taxmerged_rm_carch_raref_Richness_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Richness_Fish_rm_carch_raref_betadisp)
vegan::permutest(Richness_Fish_rm_carch_raref_betadisp)
#not significantly different dispersions for hh and non-hh with the rarefied object.
plot(Richness_Fish_rm_carch_raref_betadisp)

Richness_Fish_Beta_rm_carch_raref<-rownames_to_column(as.data.frame(Richness_Fish_rm_carch_raref_betadisp$distances),var="Sample")

##23 Shannon's Unrarefied Fish Carcharhinus distance matrix)----

Fish_taxmerged_rm_carch_unraref_Shannons_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_carch_unraref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Fish_taxmerge_rm_carch_unraref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_carch_unraref_Shannons_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_carch_unraref_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)


#calculate dispersion around the centroid

Shannons_Fish_rm_carch_unraref_betadisp<-betadisper(d = Fish_taxmerged_rm_carch_unraref_Shannons_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_carch_unraref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Shannons_Fish_rm_carch_unraref_betadisp)
vegan::permutest(Shannons_Fish_rm_carch_unraref_betadisp)
plot(Shannons_Fish_rm_carch_unraref_betadisp)

Shannons_Fish_Beta_rm_carch_unraref<-rownames_to_column(as.data.frame(Shannons_Fish_rm_carch_unraref_betadisp$distances),var="Sample")

##24 Shannon's Rarefied Fish Carcharhinus distance matrix)----

Fish_taxmerged_rm_carch_raref_Shannons_Hilldistmat<-pair_dis(countable = Fish_taxmerge_rm_carch_raref_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)#, tree = Fish_taxmerge_rm_carch_raref_eDNAindex_ultra_tree)

Fish_taxmerged_rm_carch_raref_Shannons_Hill_Jaccard_dissim<-as.dist(Fish_taxmerged_rm_carch_raref_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid

Shannons_Fish_rm_carch_raref_betadisp<-betadisper(d = Fish_taxmerged_rm_carch_raref_Shannons_Hill_Jaccard_dissim,group = Fish_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Shannons_Fish_rm_carch_raref_betadisp)
vegan::permutest(Shannons_Fish_rm_carch_raref_betadisp)
plot(Shannons_Fish_rm_carch_raref_betadisp)

Shannons_Fish_Beta_rm_carch_raref<-rownames_to_column(as.data.frame(Shannons_Fish_rm_carch_raref_betadisp$distances),var="Sample")





#Make Hammerhead Dataframes for Fish and Crustaceans----

#Faith's Fish unrarefied alpha
#Faith's Fish rarefied alpha
Hammerheads_Faith_Fish_alpha.df<-full_join(rm_hh_unraref_fish_Faiths.df,rm_hh_raref_fish_Faiths.df)
#Faith's Fish unrarefied beta
#Faith's Fish rarefied beta
Hammerheads_Faith_Fish_beta.df<-full_join(Faiths_Fish_Beta_rm_hh_unraref,Faiths_Fish_Beta_rm_hh_raref)

#Richness Fish unrarefied alpha
#Richness Fish rarefied alpha
Hammerheads_Richness_Fish_alpha.df<-full_join(rm_hh_unraref_fish_Richness.df,rm_hh_raref_fish_Richness.df)
#Richness Fish unrarefied beta
#Richness Fish rarefied beta
Hammerheads_Richness_Fish_beta.df<-full_join(Richness_Fish_Beta_rm_hh_unraref,Richness_Fish_Beta_rm_hh_raref)

#Shannons's Fish unrarefied alpha
#Shannons's Fish rarefied alpha
Hammerheads_Shannons_Fish_alpha.df<-full_join(rm_hh_unraref_fish_Shannons.df,rm_hh_raref_fish_Shannons.df)
#Shannons's Fish unrarefied beta
#Shannons's Fish rarefied beta
Hammerheads_Shannons_Fish_beta.df<-full_join(Shannons_Fish_Beta_rm_hh_unraref,Shannons_Fish_Beta_rm_hh_raref)

#Allen's Fish unrarefied alpha
#Allen's Fish rarefied alpha
Hammerheads_Allen_Fish_alpha.df<-full_join(rm_hh_unraref_fish_Allens.df,rm_hh_raref_fish_Allens.df)
#Allen's Fish unrarefied beta
#Allen's Fish rarefied beta
Hammerheads_Allen_Fish_beta.df<-full_join(Allens_Fish_Beta_rm_hh_unraref,Allens_Fish_Beta_rm_hh_raref)

#Faith's Crustacean unrarefied alpha
#Faith's Crustacean rarefied alpha
Hammerheads_Faith_Crustacean_alpha.df<-full_join(rm_hh_unraref_crustacean_Faiths.df,rm_hh_raref_crustacean_Faiths.df)
#Allen's Crustacean unrarefied alpha
#Allen's Crustacean rarefied alpha
Hammerheads_Allen_Crustacean_alpha.df<-full_join(rm_hh_unraref_crustacean_Allens.df,rm_hh_raref_crustacean_Allens.df)

#Faith's Crustacean unrarefied beta
#Faith's Crustacean rarefied beta
Hammerheads_Faith_Crustacean_beta.df<-full_join(Faiths_Crustacean_Beta_rm_hh_unraref,Faiths_Crustacean_Beta_rm_hh_raref)
#Allen's Crustacean unrarefied beta
#Allen's Crustacean rarefied beta
Hammerheads_Allen_Crustacean_beta.df<-full_join(Allens_Crustacean_Beta_rm_hh_unraref,Allens_Crustacean_Beta_rm_hh_raref)

#Shannon's Crustacean unrarefied alpha
#Shannon's Crustacean rarefied alpha
Hammerheads_Shannon_Crustacean_alpha.df<-full_join(rm_hh_unraref_crustacean_Shannons.df,rm_hh_raref_crustacean_Shannons.df)
#Shannon's Crustacean unrarefied beta
#Shannon's Crustacean rarefied beta
Hammerheads_Shannon_Crustacean_beta.df<-full_join(Shannons_Crustacean_Beta_rm_hh_unraref,Shannons_Crustacean_Beta_rm_hh_raref)


#Richness Crustacean unrarefied alpha
#Richness Crustacean rarefied alpha
Hammerheads_Richness_Crustacean_alpha.df<-full_join(rm_hh_unraref_crustacean_Richness.df,rm_hh_raref_crustacean_Richness.df)
#Richness Crustacean unrarefied beta
#Richness Crustacean rarefied beta
Hammerheads_Richness_Crustacean_beta.df<-full_join(Richness_Crustacean_Beta_rm_hh_unraref,Richness_Crustacean_Beta_rm_hh_raref)

#Hammerheads Fish Faith and Allen alpha
Hammerheads_Fish_alpha.df<-full_join(Hammerheads_Faith_Fish_alpha.df,Hammerheads_Allen_Fish_alpha.df)
#Hammerheads Fish Faith and Allen beta
Hammerheads_Fish_beta.df<-full_join(Hammerheads_Faith_Fish_beta.df,Hammerheads_Allen_Fish_beta.df)


#Hammerheads Fish Shannon and Allen alpha
Hammerheads_Fish_alpha.df<-full_join(Hammerheads_Shannons_Fish_alpha.df,Hammerheads_Allen_Fish_alpha.df)
#Hammerheads Fish Shannon and Allen beta
Hammerheads_Fish_beta.df<-full_join(Hammerheads_Shannons_Fish_beta.df,Hammerheads_Allen_Fish_beta.df)

#Hammerheads Fish Richness and Shannon alpha
Hammerheads_Fish_alpha.df<-full_join(Hammerheads_Richness_Fish_alpha.df,Hammerheads_Shannons_Fish_alpha.df)
#Hammerheads Fish Richness and Shannon beta
Hammerheads_Fish_beta.df<-full_join(Hammerheads_Richness_Fish_beta.df,Hammerheads_Shannons_Fish_beta.df)

#Hammerheads Fish Richness and Faith alpha
Hammerheads_Fish_alpha.df<-full_join(Hammerheads_Richness_Fish_alpha.df,Hammerheads_Faith_Fish_alpha.df)
#Hammerheads Fish Richness and Faith beta
Hammerheads_Fish_beta.df<-full_join(Hammerheads_Richness_Fish_beta.df,Hammerheads_Faith_Fish_beta.df)

#Hammerheads Crustacean Faith and Allen alpha
Hammerheads_Crustacean_alpha.df<-full_join(Hammerheads_Faith_Crustacean_alpha.df,Hammerheads_Allen_Crustacean_alpha.df)
#Hammerheads Crustacean Faith and Allen beta
Hammerheads_Crustacean_beta.df<-full_join(Hammerheads_Faith_Crustacean_beta.df,Hammerheads_Allen_Crustacean_beta.df)


#Hammerheads Crustacean Shannon and Allen alpha
Hammerheads_Crustacean_alpha.df<-full_join(Hammerheads_Shannon_Crustacean_alpha.df,Hammerheads_Allen_Crustacean_alpha.df)
#Hammerheads Crustacean Shannon and Allen beta
Hammerheads_Crustacean_beta.df<-full_join(Hammerheads_Shannon_Crustacean_beta.df,Hammerheads_Allen_Crustacean_beta.df)

#Hammerheads Crustacean Richness and Shannon alpha
Hammerheads_Crustacean_alpha.df<-full_join(Hammerheads_Richness_Crustacean_alpha.df,Hammerheads_Shannon_Crustacean_alpha.df)
#Hammerheads Crustacean Richness and Shannon beta
Hammerheads_Crustacean_beta.df<-full_join(Hammerheads_Richness_Crustacean_beta.df,Hammerheads_Shannon_Crustacean_beta.df)

#Hammerheads Crustacean Richness and Faith alpha
Hammerheads_Crustacean_alpha.df<-full_join(Hammerheads_Richness_Crustacean_alpha.df,Hammerheads_Faith_Crustacean_alpha.df)

#Hammerheads Crustacean Richness and Faith beta
Hammerheads_Crustacean_alpha.df<-full_join(Hammerheads_Richness_Crustacean_beta.df,Hammerheads_Faith_Crustacean_beta.df)

#Hammerheads Environmental Variables

Hammerheads_Envvars.df<-rownames_to_column(Menu_taxmerge_rm_hh_unraref_eDNAindex_envvars, var = "Sample")%>%
  dplyr::select(Sample, DD_lat, DD_long, Site_Name, Microhabitat,hammerhead_site,blacktip_site)


#Hammerheads Fish
Hammerheads_Fish.df<-full_join(Hammerheads_Fish_alpha.df,Hammerheads_Fish_beta.df)

Hammerheads_Fish.df1<-inner_join(Hammerheads_Envvars.df,Hammerheads_Fish.df)


#Hammerheads Crustaceans
Hammerheads_Crustacean.df<-full_join(Hammerheads_Crustacean_alpha.df,Hammerheads_Crustacean_beta.df)

Hammerheads_Crustacean.df1<-inner_join(Hammerheads_Envvars.df,Hammerheads_Crustacean.df)

# #Calculate means and sds and n for Fish ----
# Hammerheads_Fish_Means.df<-Hammerheads_Fish.df1%>%
#   group_by(hammerhead_site)%>%
#   summarize(mean_Fish.FaithsPD = mean(rm_hh_unraref_fish_FaithsPD_taxa_sample),
#             sd_Fish.FaithsPD=sd(rm_hh_unraref_fish_FaithsPD_taxa_sample),
#             mean_rarefied.Fish.FaithsPD = mean(rm_hh_raref_fish_FaithsPD_taxa_sample),
#             sd_rarefied.Fish.FaithsPD=sd(rm_hh_raref_fish_FaithsPD_taxa_sample),
#             mean_Fish.AllensH=mean(rm_hh_unraref_fish_AllensH_taxa_sample),
#             sd_Fish.AllensH=sd(rm_hh_unraref_fish_AllensH_taxa_sample),
#             mean_rarefied.Fish.AllensH=mean(rm_hh_raref_fish_AllensH_taxa_sample),
#             sd_rarefied.Fish.AllensH=sd(rm_hh_raref_fish_AllensH_taxa_sample),
#             
#             mean_beta.Fish.FaithsPD = mean(Faiths_Fish_rm_hh_unraref_betadisp$distances),
#             sd_beta.Fish.FaithsPD=sd(`Faiths_Fish_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Fish.FaithsPD = mean(`Faiths_Fish_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Fish.FaithsPD=sd(`Faiths_Fish_rm_hh_raref_betadisp$distances`),
#             mean_beta.Fish.AllensH=mean(`Allens_Fish_rm_hh_unraref_betadisp$distances`),
#             sd_beta.Fish.AllensH=sd(`Allens_Fish_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Fish.AllensH=mean(`Allens_Fish_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Fish.AllensH=sd(`Allens_Fish_rm_hh_raref_betadisp$distances`),
#             n=n())

# #Calculate means and sds and n for Fish ----
# Hammerheads_Fish_Means.df<-Hammerheads_Fish.df1%>%
#   group_by(hammerhead_site)%>%
#   summarize(mean_Fish.ShannonsExp = mean(rm_hh_unraref_fish_ShannonsExp_taxa_sample),
#             sd_Fish.ShannonsExp=sd(rm_hh_unraref_fish_ShannonsExp_taxa_sample),
#             mean_rarefied.Fish.ShannonsExp = mean(rm_hh_raref_fish_ShannonsExp_taxa_sample),
#             sd_rarefied.Fish.ShannonsExp=sd(rm_hh_raref_fish_ShannonsExp_taxa_sample),
#             mean_Fish.AllensH=mean(rm_hh_unraref_fish_AllensH_taxa_sample),
#             sd_Fish.AllensH=sd(rm_hh_unraref_fish_AllensH_taxa_sample),
#             mean_rarefied.Fish.AllensH=mean(rm_hh_raref_fish_AllensH_taxa_sample),
#             sd_rarefied.Fish.AllensH=sd(rm_hh_raref_fish_AllensH_taxa_sample),
#             
#             mean_beta.Fish.ShannonsExp = mean(Shannons_Fish_rm_hh_unraref_betadisp$distances),
#             sd_beta.Fish.ShannonsExp=sd(`Shannons_Fish_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Fish.ShannonsExp = mean(`Shannons_Fish_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Fish.ShannonsExp=sd(`Shannons_Fish_rm_hh_raref_betadisp$distances`),
#             mean_beta.Fish.AllensH=mean(`Allens_Fish_rm_hh_unraref_betadisp$distances`),
#             sd_beta.Fish.AllensH=sd(`Allens_Fish_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Fish.AllensH=mean(`Allens_Fish_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Fish.AllensH=sd(`Allens_Fish_rm_hh_raref_betadisp$distances`),
#             n=n())

# #Calculate means and sds and n for Fish ----
# Hammerheads_Fish_Means.df<-Hammerheads_Fish.df1%>%
#   group_by(hammerhead_site)%>%
#   summarize(mean_Fish.Richness = mean(rm_hh_unraref_fish_Richness_taxa_sample),
#             sd_Fish.Richness=sd(rm_hh_unraref_fish_Richness_taxa_sample),
#             mean_rarefied.Fish.Richness = mean(rm_hh_raref_fish_Richness_taxa_sample),
#             sd_rarefied.Fish.Richness=sd(rm_hh_raref_fish_Richness_taxa_sample),
#             mean_Fish.ShannonsExp=mean(rm_hh_unraref_fish_ShannonsExp_taxa_sample),
#             sd_Fish.ShannonsExp=sd(rm_hh_unraref_fish_ShannonsExp_taxa_sample),
#             mean_rarefied.Fish.ShannonsExp=mean(rm_hh_raref_fish_ShannonsExp_taxa_sample),
#             sd_rarefied.Fish.ShannonsExp=sd(rm_hh_raref_fish_ShannonsExp_taxa_sample),
#             
#             mean_beta.Fish.Richness = mean(Richness_Fish_rm_hh_unraref_betadisp$distances),
#             sd_beta.Fish.Richness=sd(`Richness_Fish_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Fish.Richness = mean(`Richness_Fish_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Fish.Richness=sd(`Richness_Fish_rm_hh_raref_betadisp$distances`),
#             mean_beta.Fish.ShannonsExp=mean(`Shannons_Fish_rm_hh_unraref_betadisp$distances`),
#             sd_beta.Fish.ShannonsExp=sd(`Shannons_Fish_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Fish.ShannonsExp=mean(`Shannons_Fish_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Fish.ShannonsExp=sd(`Shannons_Fish_rm_hh_raref_betadisp$distances`),
#             n=n())

#Calculate means and sds and n for Fish ----
Hammerheads_Fish_Means.df<-Hammerheads_Fish.df1%>%
  group_by(hammerhead_site)%>%
  summarize(mean_Fish.Richness = mean(rm_hh_unraref_fish_Richness_taxa_sample),
            sd_Fish.Richness=sd(rm_hh_unraref_fish_Richness_taxa_sample),
            mean_rarefied.Fish.Richness = mean(rm_hh_raref_fish_Richness_taxa_sample),
            sd_rarefied.Fish.Richness=sd(rm_hh_raref_fish_Richness_taxa_sample),
            mean_Fish.FaithsPD = mean(rm_hh_unraref_fish_FaithsPD_taxa_sample),
            sd_Fish.FaithsPD=sd(rm_hh_unraref_fish_FaithsPD_taxa_sample),
            mean_rarefied.Fish.FaithsPD = mean(rm_hh_raref_fish_FaithsPD_taxa_sample),
            sd_rarefied.Fish.FaithsPD=sd(rm_hh_raref_fish_FaithsPD_taxa_sample),
            
            mean_beta.Fish.Richness = mean(Richness_Fish_rm_hh_unraref_betadisp$distances),
            sd_beta.Fish.Richness=sd(`Richness_Fish_rm_hh_unraref_betadisp$distances`),
            mean_beta.rarefied.Fish.Richness = mean(`Richness_Fish_rm_hh_raref_betadisp$distances`),
            sd_beta.rarefied.Fish.Richness=sd(`Richness_Fish_rm_hh_raref_betadisp$distances`),
            mean_beta.Fish.FaithsPD = mean(Faiths_Fish_rm_hh_unraref_betadisp$distances),
            sd_beta.Fish.FaithsPD=sd(`Faiths_Fish_rm_hh_unraref_betadisp$distances`),
            mean_beta.rarefied.Fish.FaithsPD = mean(`Faiths_Fish_rm_hh_raref_betadisp$distances`),
            sd_beta.rarefied.Fish.FaithsPD=sd(`Faiths_Fish_rm_hh_raref_betadisp$distances`),
            n=n())

# #Calculate means and sds and n for Crustaceans ----
# Hammerheads_Crustacean_Means.df<-Hammerheads_Crustacean.df1%>%
#   group_by(hammerhead_site)%>%
#   summarize(mean_Crustacean.FaithsPD = mean(rm_hh_unraref_crustacean_FaithsPD_taxa_sample),
#             sd_Crustacean.FaithsPD=sd(rm_hh_unraref_crustacean_FaithsPD_taxa_sample),
#             mean_rarefied.Crustacean.FaithsPD = mean(rm_hh_raref_crustacean_FaithsPD_taxa_sample),
#             sd_rarefied.Crustacean.FaithsPD=sd(rm_hh_raref_crustacean_FaithsPD_taxa_sample),
#             mean_Crustacean.AllensH=mean(rm_hh_unraref_crustacean_AllensH_taxa_sample),
#             sd_Crustacean.AllensH=sd(rm_hh_unraref_crustacean_AllensH_taxa_sample),
#             mean_rarefied.Crustacean.AllensH=mean(rm_hh_raref_crustacean_AllensH_taxa_sample),
#             sd_rarefied.Crustacean.AllensH=sd(rm_hh_raref_crustacean_AllensH_taxa_sample),
#             
#             mean_beta.Crustacean.FaithsPD = mean(Faiths_Crustacean_rm_hh_unraref_betadisp$distances),
#             sd_beta.Crustacean.FaithsPD=sd(`Faiths_Crustacean_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Crustacean.FaithsPD = mean(`Faiths_Crustacean_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Crustacean.FaithsPD=sd(`Faiths_Crustacean_rm_hh_raref_betadisp$distances`),
#             mean_beta.Crustacean.AllensH=mean(`Allens_Crustacean_rm_hh_unraref_betadisp$distances`),
#             sd_beta.Crustacean.AllensH=sd(`Allens_Crustacean_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Crustacean.AllensH=mean(`Allens_Crustacean_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Crustacean.AllensH=sd(`Allens_Crustacean_rm_hh_raref_betadisp$distances`),
#             n=n())



# #Calculate means and sds and n for Crustaceans ----
# Hammerheads_Crustacean_Means.df<-Hammerheads_Crustacean.df1%>%
#   group_by(hammerhead_site)%>%
#   summarize(mean_Crustacean.ShannonsExp = mean(rm_hh_unraref_crustacean_ShannonsExp_taxa_sample),
#             sd_Crustacean.ShannonsExp=sd(rm_hh_unraref_crustacean_ShannonsExp_taxa_sample),
#             mean_rarefied.Crustacean.ShannonsExp = mean(rm_hh_raref_crustacean_ShannonsExp_taxa_sample),
#             sd_rarefied.Crustacean.ShannonsExp=sd(rm_hh_raref_crustacean_ShannonsExp_taxa_sample),
#             mean_Crustacean.AllensH=mean(rm_hh_unraref_crustacean_AllensH_taxa_sample),
#             sd_Crustacean.AllensH=sd(rm_hh_unraref_crustacean_AllensH_taxa_sample),
#             mean_rarefied.Crustacean.AllensH=mean(rm_hh_raref_crustacean_AllensH_taxa_sample),
#             sd_rarefied.Crustacean.AllensH=sd(rm_hh_raref_crustacean_AllensH_taxa_sample),
#             
#             mean_beta.Crustacean.ShannonsExp = mean(Shannons_Crustacean_rm_hh_unraref_betadisp$distances),
#             sd_beta.Crustacean.ShannonsExp=sd(`Shannons_Crustacean_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Crustacean.ShannonsExp = mean(`Shannons_Crustacean_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Crustacean.ShannonsExp=sd(`Shannons_Crustacean_rm_hh_raref_betadisp$distances`),
#             mean_beta.Crustacean.AllensH=mean(`Allens_Crustacean_rm_hh_unraref_betadisp$distances`),
#             sd_beta.Crustacean.AllensH=sd(`Allens_Crustacean_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Crustacean.AllensH=mean(`Allens_Crustacean_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Crustacean.AllensH=sd(`Allens_Crustacean_rm_hh_raref_betadisp$distances`),
#             n=n())


# #Calculate means and sds and n for Crustaceans ----
# Hammerheads_Crustacean_Means.df<-Hammerheads_Crustacean.df1%>%
#   group_by(hammerhead_site)%>%
#   summarize(mean_Crustacean.Richness = mean(rm_hh_unraref_crustacean_Richness_taxa_sample),
#             sd_Crustacean.Richness=sd(rm_hh_unraref_crustacean_Richness_taxa_sample),
#             mean_rarefied.Crustacean.Richness = mean(rm_hh_raref_crustacean_Richness_taxa_sample),
#             sd_rarefied.Crustacean.Richness=sd(rm_hh_raref_crustacean_Richness_taxa_sample),
#             mean_Crustacean.ShannonsExp=mean(rm_hh_unraref_crustacean_ShannonsExp_taxa_sample),
#             sd_Crustacean.ShannonsExp=sd(rm_hh_unraref_crustacean_ShannonsExp_taxa_sample),
#             mean_rarefied.Crustacean.ShannonsExp=mean(rm_hh_raref_crustacean_ShannonsExp_taxa_sample),
#             sd_rarefied.Crustacean.ShannonsExp=sd(rm_hh_raref_crustacean_ShannonsExp_taxa_sample),
#             
#             mean_beta.Crustacean.Richness = mean(Richness_Crustacean_rm_hh_unraref_betadisp$distances),
#             sd_beta.Crustacean.Richness=sd(`Richness_Crustacean_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Crustacean.Richness = mean(`Richness_Crustacean_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Crustacean.Richness=sd(`Richness_Crustacean_rm_hh_raref_betadisp$distances`),
#             mean_beta.Crustacean.ShannonsExp=mean(`Shannons_Crustacean_rm_hh_unraref_betadisp$distances`),
#             sd_beta.Crustacean.ShannonsExp=sd(`Shannons_Crustacean_rm_hh_unraref_betadisp$distances`),
#             mean_beta.rarefied.Crustacean.ShannonsExp=mean(`Shannons_Crustacean_rm_hh_raref_betadisp$distances`),
#             sd_beta.rarefied.Crustacean.ShannonsExp=sd(`Shannons_Crustacean_rm_hh_raref_betadisp$distances`),
#             n=n())


#Calculate means and sds and n for Crustaceans ----
Hammerheads_Crustacean_Means.df<-Hammerheads_Crustacean.df1%>%
  group_by(hammerhead_site)%>%
  summarize(mean_Crustacean.Richness = mean(rm_hh_unraref_crustacean_Richness_taxa_sample),
            sd_Crustacean.Richness=sd(rm_hh_unraref_crustacean_Richness_taxa_sample),
            mean_rarefied.Crustacean.Richness = mean(rm_hh_raref_crustacean_Richness_taxa_sample),
            sd_rarefied.Crustacean.Richness=sd(rm_hh_raref_crustacean_Richness_taxa_sample),
            mean_Crustacean.FaithsPD = mean(rm_hh_unraref_crustacean_FaithsPD_taxa_sample),
            sd_Crustacean.FaithsPD=sd(rm_hh_unraref_crustacean_FaithsPD_taxa_sample),
            mean_rarefied.Crustacean.FaithsPD = mean(rm_hh_raref_crustacean_FaithsPD_taxa_sample),
            sd_rarefied.Crustacean.FaithsPD=sd(rm_hh_raref_crustacean_FaithsPD_taxa_sample),
            
            mean_beta.Crustacean.Richness = mean(Richness_Crustacean_rm_hh_unraref_betadisp$distances),
            sd_beta.Crustacean.Richness=sd(`Richness_Crustacean_rm_hh_unraref_betadisp$distances`),
            mean_beta.rarefied.Crustacean.Richness = mean(`Richness_Crustacean_rm_hh_raref_betadisp$distances`),
            sd_beta.rarefied.Crustacean.Richness=sd(`Richness_Crustacean_rm_hh_raref_betadisp$distances`),
            mean_beta.Crustacean.FaithsPD = mean(Faiths_Crustacean_rm_hh_unraref_betadisp$distances),
            sd_beta.Crustacean.FaithsPD=sd(`Faiths_Crustacean_rm_hh_unraref_betadisp$distances`),
            mean_beta.rarefied.Crustacean.FaithsPD = mean(`Faiths_Crustacean_rm_hh_raref_betadisp$distances`),
            sd_beta.rarefied.Crustacean.FaithsPD=sd(`Faiths_Crustacean_rm_hh_raref_betadisp$distances`),
            n=n())


#Re-shape the Fish dataframe----
hammerheads_fish_diversities_long<-Hammerheads_Fish_Means.df%>%
  pivot_longer(
    cols=-c(hammerhead_site,n),
    names_to = c(".value","measure"),
    names_sep = "_"
  )

#make it wider by hh or non-hh site
hammerheads_fish_metameansd<-hammerheads_fish_diversities_long %>% 
  pivot_longer(c(-measure, -hammerhead_site)) %>%
  pivot_wider(names_from = c(hammerhead_site, name))

hammerheads_fish_metameansd_num<- hammerheads_fish_metameansd%>%
  mutate(across(hh_site_n:`non-hh_site_sd`,as.numeric))

is.numeric(hammerheads_fish_metameansd_num$hh_site_mean)

#Re-shape the Crustacean dataframe----
hammerheads_crustacean_diversities_long<-Hammerheads_Crustacean_Means.df%>%
  pivot_longer(
    cols=-c(hammerhead_site,n),
    names_to = c(".value","measure"),
    names_sep = "_"
  )

#make it wider by hh or non-hh site
hammerheads_crustacean_metameansd<-hammerheads_crustacean_diversities_long %>% 
  pivot_longer(c(-measure, -hammerhead_site)) %>%
  pivot_wider(names_from = c(hammerhead_site, name))

hammerheads_crustacean_metameansd_num<- hammerheads_crustacean_metameansd%>%
  mutate(across(hh_site_n:`non-hh_site_sd`,as.numeric))

is.numeric(hammerheads_crustacean_metameansd_num$hh_site_mean)

#Make Hammerhead Menu Dataframe----
Hammerheads_Menu_alpha.df<-full_join(rm_hh_unraref_Faiths.df,rm_hh_raref_Faiths.df)

Hammerheads_Menu_beta.df<-full_join(Faiths_Menu_Beta_rm_hh_unraref,Faiths_Menu_Beta_rm_hh_raref)

Hammerheads_Menu_diversities<-full_join(Hammerheads_Menu_alpha.df,Hammerheads_Menu_beta.df)
#Hammerheads Environmental Variables

Hammerheads_Envvars.df<-rownames_to_column(Menu_taxmerge_rm_hh_unraref_eDNAindex_envvars, var = "Sample")%>%
  dplyr::select(Sample, DD_lat, DD_long, Site_Name, Microhabitat,hammerhead_site,blacktip_site)

Hammerheads_Menu_Diversity.df<-inner_join(Hammerheads_Envvars.df,Hammerheads_Menu_diversities)

#Calculate means and sds and n for Menu ----
Hammerheads_Menu_Means.df<-Hammerheads_Menu_Diversity.df%>%
  group_by(hammerhead_site)%>%
  summarize(mean_Menu.FaithsPD = mean(rm_hh_unraref_FaithsPD_taxa_sample),
            sd_Menu.FaithsPD=sd(rm_hh_unraref_FaithsPD_taxa_sample),
            mean_rarefied.Menu.FaithsPD = mean(rm_hh_raref_FaithsPD_taxa_sample),
            sd_rarefied.Menu.FaithsPD=sd(rm_hh_raref_FaithsPD_taxa_sample),
            mean_beta.Menu.FaithsPD = mean(Faiths_Menu_rm_hh_unraref_betadisp$distances),
            sd_beta.Menu.FaithsPD=sd(`Faiths_Menu_rm_hh_unraref_betadisp$distances`),
            mean_beta.rarefied.Menu.FaithsPD = mean(`Faiths_Menu_rm_hh_raref_betadisp$distances`),
            sd_beta.rarefied.Menu.FaithsPD=sd(`Faiths_Menu_rm_hh_raref_betadisp$distances`),
            n=n())

#Re-shape the Menu dataframe----
hammerheads_menu_diversities_long<-Hammerheads_Menu_Means.df%>%
  pivot_longer(
    cols=-c(hammerhead_site,n),
    names_to = c(".value","measure"),
    names_sep = "_"
  )

#make it wider by hh or non-hh site
hammerheads_menu_metameansd<-hammerheads_menu_diversities_long %>% 
  pivot_longer(c(-measure, -hammerhead_site)) %>%
  pivot_wider(names_from = c(hammerhead_site, name))

hammerheads_menu_metameansd_num<- hammerheads_menu_metameansd%>%
  mutate(across(hh_site_n:`non-hh_site_sd`,as.numeric))

is.numeric(hammerheads_menu_metameansd_num$hh_site_mean)

#metafor Menu Hammerheads random effects model and graph----

#install.packages("metafor")
library(metafor) #masks permutest from vegan

#calculate yi and vi for hammerheads menu
hh_menu_logresponse_meta<-escalc(measure = "ROM",m1i = hh_site_mean,sd1i = hh_site_sd,m2i = `non-hh_site_mean`,sd2i = `non-hh_site_sd`,n1i = hh_site_n,n2i = `non-hh_site_n`,data = hammerheads_menu_metameansd_num)

# fit random-effects model

hh_menu_res<-rma(yi, vi, data=hh_menu_logresponse_meta, test = "t")
hh_menu_res$data

# predicted pooled risk ratio (with 95% confidence/prediction intervals)
predict(hh_menu_res,transf = exp,digits = 2)

plot.new()
forest(hh_menu_res, atransf = exp,slab = measure,header = "Hammerhead Predator Effects on Fish and Crustaceans",shade="zebra")




#metafor random effects model and graph Hammerheads----
#meta analysis of Fish and Crustaceans
#install.packages("metafor")
library(metafor) #masks permutest from vegan

#calculate yi and vi for hammerheads fish
hh_fish_logresponse_meta<-escalc(measure = "ROM",m1i = hh_site_mean,sd1i = hh_site_sd,m2i = `non-hh_site_mean`,sd2i = `non-hh_site_sd`,n1i = hh_site_n,n2i = `non-hh_site_n`,data = hammerheads_fish_metameansd_num)

# fit random-effects model

hh_fish_res<-rma(yi, vi, data=hh_fish_logresponse_meta, test = "t")
hh_fish_res$data

# predicted pooled risk ratio (with 95% confidence/prediction intervals)
predict(hh_fish_res,transf = exp,digits = 2)

plot.new()
forest(hh_fish_res, atransf = exp,slab = measure,header = "Hammerhead Predator Effects on Fish",shade="zebra")






#calculate yi and vi for hammerheads crustaceans
hh_crustacean_logresponse_meta<-escalc(measure = "ROM",m1i = hh_site_mean,sd1i = hh_site_sd,m2i = `non-hh_site_mean`,sd2i = `non-hh_site_sd`,n1i = hh_site_n,n2i = `non-hh_site_n`,data = hammerheads_crustacean_metameansd_num)

# fit random-effects model

hh_crustacean_res<-rma(yi, vi, data=hh_crustacean_logresponse_meta, test = "t")
hh_crustacean_res$data

# predicted pooled risk ratio (with 95% confidence/prediction intervals)
predict(hh_crustacean_res,transf = exp,digits = 2)

plot.new()
forest(hh_crustacean_res, atransf = exp,slab = measure,header = "Hammerhead Predator Effects on Crustaceans",shade="zebra")



#Make Carcharhinus Dataframes for Fish and Crustaceans----



#Richness Fish unrarefied alpha
#Richness Fish rarefied alpha
Carcharhinus_Richness_Fish_alpha.df<-full_join(rm_carch_unraref_fish_Richness.df,rm_carch_raref_fish_Richness.df)
#Richness Fish unrarefied beta
#Richness Fish rarefied beta
Carcharhinus_Richness_Fish_beta.df<-full_join(Richness_Fish_Beta_rm_carch_unraref,Richness_Fish_Beta_rm_carch_raref)

#Shannons's Fish unrarefied alpha
#Shannons's Fish rarefied alpha
Carcharhinus_Shannons_Fish_alpha.df<-full_join(rm_carch_unraref_fish_Shannons.df,rm_carch_raref_fish_Shannons.df)
#Shannons's Fish unrarefied beta
#Shannons's Fish rarefied beta
Carcharhinus_Shannons_Fish_beta.df<-full_join(Shannons_Fish_Beta_rm_carch_unraref,Shannons_Fish_Beta_rm_carch_raref)

#Richness Crustacean unrarefied alpha
#Richness Crustacean rarefied alpha
Carcharhinus_Richness_Crustacean_alpha.df<-full_join(rm_carch_unraref_crustacean_Richness.df,rm_carch_raref_crustacean_Richness.df)
#Richness Crustacean unrarefied beta
#Richness Crustacean rarefied beta
Carcharhinus_Richness_Crustacean_beta.df<-full_join(Richness_Crustacean_Beta_rm_carch_unraref,Richness_Crustacean_Beta_rm_carch_raref)

#Shannons's Crustacean unrarefied alpha
#Shannons's Crustacean rarefied alpha
Carcharhinus_Shannons_Crustacean_alpha.df<-full_join(rm_carch_unraref_crustacean_Shannons.df,rm_carch_raref_crustacean_Shannons.df)
#Shannons's Crustacean unrarefied beta
#Shannons's Crustacean rarefied beta
Carcharhinus_Shannons_Crustacean_beta.df<-full_join(Shannons_Crustacean_Beta_rm_carch_unraref,Shannons_Crustacean_Beta_rm_carch_raref)


#Carcharhinus Fish Richness and Shannon alpha
Carcharhinus_Fish_alpha.df<-full_join(Carcharhinus_Richness_Fish_alpha.df,Carcharhinus_Shannons_Fish_alpha.df)
#Carcharhinus Fish Richness and Shannon beta
Carcharhinus_Fish_beta.df<-full_join(Carcharhinus_Richness_Fish_beta.df,Carcharhinus_Shannons_Fish_beta.df)



#Carcharhinus Crustacean Richness and Shannon alpha
Carcharhinus_Crustacean_alpha.df<-full_join(Carcharhinus_Richness_Crustacean_alpha.df,Carcharhinus_Shannons_Crustacean_alpha.df)
#Carcharhinus Crustacean Richness and Shannon beta
Carcharhinus_Crustacean_beta.df<-full_join(Carcharhinus_Richness_Crustacean_beta.df,Carcharhinus_Shannons_Crustacean_beta.df)



#Carcharhinus Fish Richness and Faith alpha
Carcharhinus_Fish_alpha.df<-full_join(Carcharhinus_Richness_Fish_alpha.df,Carcharhinus_Faith_Fish_alpha.df)
#Carcharhinus Fish Richness and Faith beta
Carcharhinus_Fish_beta.df<-full_join(Carcharhinus_Richness_Fish_beta.df,Carcharhinus_Faith_Fish_beta.df)



#Carcharhinus Crustacean Richness and Faith alpha
Carcharhinus_Crustacean_alpha.df<-full_join(Carcharhinus_Richness_Crustacean_alpha.df,Carcharhinus_Faith_Crustacean_alpha.df)
#Carcharhinus Crustacean Richness and Faith beta
Carcharhinus_Crustacean_beta.df<-full_join(Carcharhinus_Richness_Crustacean_beta.df,Carcharhinus_Faith_Crustacean_beta.df)



#Carcharhinus Environmental Variables

Carcharhinus_Envvars.df<-rownames_to_column(Menu_taxmerge_rm_carch_unraref_eDNAindex_envvars, var = "Sample")%>%
  dplyr::select(Sample, DD_lat, DD_long, Site_Name, Microhabitat,hammerhead_site,blacktip_site)


#Carcharhinus Fish
Carcharhinus_Fish.df<-full_join(Carcharhinus_Fish_alpha.df,Carcharhinus_Fish_beta.df)

Carcharhinus_Fish.df1<-inner_join(Carcharhinus_Envvars.df,Carcharhinus_Fish.df)


#Carcharhinus Crustaceans
Carcharhinus_Crustacean.df<-full_join(Carcharhinus_Crustacean_alpha.df,Carcharhinus_Crustacean_beta.df)

Carcharhinus_Crustacean.df1<-inner_join(Carcharhinus_Envvars.df,Carcharhinus_Crustacean.df)

#Calculate means and sds and n for Fish ----
Carcharhinus_Fish_Means.df<-Carcharhinus_Fish.df1%>%
  group_by(blacktip_site)%>%
  summarize(mean_Fish.Richness = mean(rm_carch_unraref_fish_Richness_taxa_sample),
            sd_Fish.Richness=sd(rm_carch_unraref_fish_Richness_taxa_sample),
            mean_rarefied.Fish.Richness = mean(rm_carch_raref_fish_Richness_taxa_sample),
            sd_rarefied.Fish.Richness=sd(rm_carch_raref_fish_Richness_taxa_sample),
            mean_Fish.ShannonsExp=mean(rm_carch_unraref_fish_ShannonsExp_taxa_sample),
            sd_Fish.ShannonsExp=sd(rm_carch_unraref_fish_ShannonsExp_taxa_sample),
            mean_rarefied.Fish.ShannonsExp=mean(rm_carch_raref_fish_ShannonsExp_taxa_sample),
            sd_rarefied.Fish.ShannonsExp=sd(rm_carch_raref_fish_ShannonsExp_taxa_sample),
            
            mean_beta.Fish.Richness = mean(Richness_Fish_rm_carch_unraref_betadisp$distances),
            sd_beta.Fish.Richness=sd(`Richness_Fish_rm_carch_unraref_betadisp$distances`),
            mean_beta.rarefied.Fish.Richness = mean(`Richness_Fish_rm_carch_raref_betadisp$distances`),
            sd_beta.rarefied.Fish.Richness=sd(`Richness_Fish_rm_carch_raref_betadisp$distances`),
            mean_beta.Fish.ShannonsExp=mean(`Shannons_Fish_rm_carch_unraref_betadisp$distances`),
            sd_beta.Fish.ShannonsExp=sd(`Shannons_Fish_rm_carch_unraref_betadisp$distances`),
            mean_beta.rarefied.Fish.ShannonsExp=mean(`Shannons_Fish_rm_carch_raref_betadisp$distances`),
            sd_beta.rarefied.Fish.ShannonsExp=sd(`Shannons_Fish_rm_carch_raref_betadisp$distances`),
            n=n())


#Calculate means and sds and n for Crustaceans ----
Carcharhinus_Crustacean_Means.df<-Carcharhinus_Crustacean.df1%>%
  group_by(blacktip_site)%>%
  summarize(mean_Crustacean.Richness = mean(rm_carch_unraref_crustacean_Richness_taxa_sample),
            sd_Crustacean.Richness=sd(rm_carch_unraref_crustacean_Richness_taxa_sample),
            mean_rarefied.Crustacean.Richness = mean(rm_carch_raref_crustacean_Richness_taxa_sample),
            sd_rarefied.Crustacean.Richness=sd(rm_carch_raref_crustacean_Richness_taxa_sample),
            mean_Crustacean.ShannonsExp=mean(rm_carch_unraref_crustacean_ShannonsExp_taxa_sample),
            sd_Crustacean.ShannonsExp=sd(rm_carch_unraref_crustacean_ShannonsExp_taxa_sample),
            mean_rarefied.Crustacean.ShannonsExp=mean(rm_carch_raref_crustacean_ShannonsExp_taxa_sample),
            sd_rarefied.Crustacean.ShannonsExp=sd(rm_carch_raref_crustacean_ShannonsExp_taxa_sample),
            
            mean_beta.Crustacean.Richness = mean(Richness_Crustacean_rm_carch_unraref_betadisp$distances),
            sd_beta.Crustacean.Richness=sd(`Richness_Crustacean_rm_carch_unraref_betadisp$distances`),
            mean_beta.rarefied.Crustacean.Richness = mean(`Richness_Crustacean_rm_carch_raref_betadisp$distances`),
            sd_beta.rarefied.Crustacean.Richness=sd(`Richness_Crustacean_rm_carch_raref_betadisp$distances`),
            mean_beta.Crustacean.ShannonsExp=mean(`Shannons_Crustacean_rm_carch_unraref_betadisp$distances`),
            sd_beta.Crustacean.ShannonsExp=sd(`Shannons_Crustacean_rm_carch_unraref_betadisp$distances`),
            mean_beta.rarefied.Crustacean.ShannonsExp=mean(`Shannons_Crustacean_rm_carch_raref_betadisp$distances`),
            sd_beta.rarefied.Crustacean.ShannonsExp=sd(`Shannons_Crustacean_rm_carch_raref_betadisp$distances`),
            n=n())


#Re-shape the Fish dataframe----
carcharhinus_fish_diversities_long<-Carcharhinus_Fish_Means.df%>%
  pivot_longer(
    cols=-c(blacktip_site,n),
    names_to = c(".value","measure"),
    names_sep = "_"
  )

#make it wider by carch or non-carch site
carcharhinus_fish_metameansd<-carcharhinus_fish_diversities_long %>% 
  pivot_longer(c(-measure, -blacktip_site)) %>%
  pivot_wider(names_from = c(blacktip_site, name))

carcharhinus_fish_metameansd_num<- carcharhinus_fish_metameansd%>%
  mutate(across(carch_site_n:`non-carch_site_sd`,as.numeric))

is.numeric(carcharhinus_fish_metameansd_num$carch_site_mean)

#Re-shape the Crustacean dataframe----
carcharhinus_crustacean_diversities_long<-Carcharhinus_Crustacean_Means.df%>%
  pivot_longer(
    cols=-c(blacktip_site,n),
    names_to = c(".value","measure"),
    names_sep = "_"
  )

#make it wider by carch or non-carch site
carcharhinus_crustacean_metameansd<-carcharhinus_crustacean_diversities_long %>% 
  pivot_longer(c(-measure, -blacktip_site)) %>%
  pivot_wider(names_from = c(blacktip_site, name))

carcharhinus_crustacean_metameansd_num<- carcharhinus_crustacean_metameansd%>%
  mutate(across(carch_site_n:`non-carch_site_sd`,as.numeric))

is.numeric(carcharhinus_crustacean_metameansd_num$carch_site_mean)


#metafor random effects model and graph Carcharhinus----
#meta analysis of Fish and Crustaceans
#install.packages("metafor")
library(metafor) #masks permutest from vegan

#calculate yi and vi for carcharhinus fish
carch_fish_logresponse_meta<-escalc(measure = "ROM",m1i = carch_site_mean,sd1i = carch_site_sd,m2i = `non-carch_site_mean`,sd2i = `non-carch_site_sd`,n1i = carch_site_n,n2i = `non-carch_site_n`,data = carcharhinus_fish_metameansd_num)

# fit random-effects model

carch_fish_res<-rma(yi, vi, data=carch_fish_logresponse_meta, test = "t")
carch_fish_res$data

# predicted pooled risk ratio (with 95% confidence/prediction intervals)
predict(carch_fish_res,transf = exp,digits = 2)

plot.new()
forest(carch_fish_res, atransf = exp,slab = measure,header = "Carcharhinus Predator Effects on Fish",shade="zebra")






#calculate yi and vi for carcharhinus crustaceans
carch_crustacean_logresponse_meta<-escalc(measure = "ROM",m1i = carch_site_mean,sd1i = carch_site_sd,m2i = `non-carch_site_mean`,sd2i = `non-carch_site_sd`,n1i = carch_site_n,n2i = `non-carch_site_n`,data = carcharhinus_crustacean_metameansd_num)

# fit random-effects model

carch_crustacean_res<-rma(yi, vi, data=carch_crustacean_logresponse_meta, test = "t")
carch_crustacean_res$data

# predicted pooled risk ratio (with 95% confidence/prediction intervals)
predict(carch_crustacean_res,transf = exp,digits = 2)

plot.new()
forest(carch_crustacean_res, atransf = exp,slab = measure,header = "Carcharhinus Predator Effects on Crustaceans",shade="zebra")



##### END of log response ratios ----------
#log response ratio https://search.r-project.org/CRAN/refmans/ARPobservation/html/logRespRatio.html

#install.packages("ARPobservation")
library(ARPobservation)
hh_raref_beta_lrr<-logRespRatio(
  Faiths_rm_hh_raref_betadisp$distances,
  phase = Menu_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site,
  base_level = "non-hh_site",
  conf_level = 0.95,
  bias_correct = TRUE,
  exponentiate = FALSE
)
hh_raref_beta_lrr




#try renyi accumulation curves (like species accumulation curves but with hill numbers)----
# data(BCI)
# i <- sample(nrow(BCI), 12)
# mod <- renyi(BCI[i,])
# plot(mod)
# mod <- renyiaccum(BCI[i,])
# plot(mod, as.table=TRUE, col = c(1, 2, 2))
# persp(mod)
# 
# raref_hh_renyi_otudf<-t(Menu_taxmerge_rm_hh_raref_eDNAindex_otudf)
# hhsitelist<-rownames(Menu_taxmerge_rm_hh_raref_eDNAindex_envvars%>%
#   dplyr::filter(hammerhead_site=="hh_site"))
# 
# rarefhhsitemod<-renyiaccum(raref_hh_renyi_otudf,scales = c(0,0.5,1,2),subset = rownames(raref_hh_renyi_otudf)%in% hhsitelist)
# rarefhhsiteplot<-plot(rarefhhsitemod, as.table=TRUE)
# rarefhhsiteplot
# persp(rarefhhsitemod)
# 
# rarefnonhhsitemod<-renyiaccum(raref_hh_renyi_otudf, scales = c(0,0.5,1,2),subset = rownames(raref_hh_renyi_otudf)%ni% hhsitelist)
# 
# 
# rarefnonhhsiteplot<-plot(rarefnonhhsitemod, as.table=TRUE)
# rarefnonhhsiteplot
# 
# ##plot on top of each other 

# # doesn't work
#plot(hhsitemod)
#plot(nonhhsitemod, add = TRUE)
# 

#try this method: it's plain richness like specaccum and overplotting with ggplot2----
#https://rpubs.com/Roeland-KINDT/694066

# data(dune.env)
# data(dune)
# dune.env$site.totals <- apply(dune,1,sum)
# Accum.1 <- accumresult(dune, y=dune.env, scale='site.totals', method='exact', conditioned=TRUE)
# Accum.1
# accumplot(Accum.1)
# 
# Accum.2 <- accumcomp(dune, y=dune.env, factor='Management', method='exact', 
#                      legend=FALSE, conditioned=TRUE, scale='site.totals')
## CLICK IN THE GRAPH TO INDICATE WHERE THE LEGEND NEEDS TO BE PLACED FOR
## OPTION WHERE LEGEND=TRUE (DEFAULT).

## Not run: 
# ggplot2 plotting method
# 
# data(warcom)
# data(warenv)
# 
# Accum.3 <- accumcomp(warcom, y=warenv, factor='population', 
#                      method='exact', conditioned=F, plotit=F)

raref_hh_otudf<-t(Menu_taxmerge_rm_hh_raref_eDNAindex_otudf)
raref_hh_envvars<-Menu_taxmerge_rm_hh_raref_eDNAindex_envvars
raref_hh_envvars$hammerhead_site<-as.factor(raref_hh_envvars$hammerhead_site)

#Hammerheads Menu accum plot----
Menu_hh_accum <- accumcomp(x=raref_hh_otudf, y=raref_hh_envvars, factor='hammerhead_site', 
                     method='exact', conditioned=F, plotit=F)


#Hammerheads Fish accum plot----
Fish_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site<-as.factor(Fish_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site)
trans_raref_hh_fish<-t(Fish_taxmerge_rm_hh_raref_eDNAindex_otudf)

hh_Fish_accum<-accumcomp(x=trans_raref_hh_fish,y=Fish_taxmerge_rm_hh_raref_eDNAindex_envvars, factor ='hammerhead_site',method = 'exact',conditioned=F,plotit = F)

#Hammerheads Crustacean accum plot----
trans_raref_hh_crustaceans<-t(Crustacean_taxmerge_rm_hh_raref_eDNAindex_otudf)
Crustacean_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site<-as.factor(Crustacean_taxmerge_rm_hh_raref_eDNAindex_envvars$hammerhead_site)

hh_Crustacean_accum<-accumcomp(x=trans_raref_hh_crustaceans,y = Crustacean_taxmerge_rm_hh_raref_eDNAindex_envvars, factor = 'hammerhead_site',method = 'exact',conditioned = F, plotit = F)


#Carcharhinus Fish accum plot----
Fish_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site<-as.factor(Fish_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site)
trans_raref_carch_fish<-t(Fish_taxmerge_rm_carch_raref_eDNAindex_otudf)

carch_Fish_accum<-accumcomp(x=trans_raref_carch_fish,y=Fish_taxmerge_rm_carch_raref_eDNAindex_envvars, factor ='blacktip_site',method = 'exact',conditioned=F,plotit = F)

#Carcharhinus Crustacean accum plot----
trans_raref_carch_crustaceans<-t(Crustacean_taxmerge_rm_carch_raref_eDNAindex_otudf)
Crustacean_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site<-as.factor(Crustacean_taxmerge_rm_carch_raref_eDNAindex_envvars$blacktip_site)

carch_Crustacean_accum<-accumcomp(x=trans_raref_carch_crustaceans,y = Crustacean_taxmerge_rm_carch_raref_eDNAindex_envvars, factor = 'blacktip_site',method = 'exact',conditioned = F, plotit = F)







library(ggplot2)

# possibly need for extrafont::loadfonts(device="win") to have Arial
# as alternative, use library(ggThemeAssist)
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

accum.long3 <- accumcomp.long(Menu_hh_accum, ci=NA, label.freq=5)
hh_Fish_accum.long<-accumcomp.long(hh_Fish_accum,ci=NA,label.freq = 5)
hh_Crustacean_accum.long<-accumcomp.long(hh_Crustacean_accum,ci=NA,label.freq = 5)

hammerheads_fish_plot <- ggplot(data=hh_Fish_accum.long, aes(x = Sites, y = Richness, ymax =  UPR, ymin= LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=2) +
  geom_point(data=subset(hh_Fish_accum.long, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping), alpha=0.2, show.legend=FALSE) + 
  BioR.theme +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Rarefied Samples", y = "Fish Richness", colour = "Hammerhead Site Status", shape = "Hammerhead Site Status")

hammerheads_fish_plot

hammerheads_crustacean_plot <- ggplot(data=hh_Crustacean_accum.long, aes(x = Sites, y = Richness, ymax =  UPR, ymin= LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=2) +
  geom_point(data=subset(hh_Crustacean_accum.long, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping), alpha=0.2, show.legend=FALSE) + 
  BioR.theme +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Rarefied Samples", y = "Crustacean Richness", colour = "Hammerhead Site Status", shape = "Hammerhead Site Status")

hammerheads_crustacean_plot


#same for Carcharhinus

carch_Fish_accum.long<-accumcomp.long(carch_Fish_accum,ci=NA,label.freq = 5)
carch_Crustacean_accum.long<-accumcomp.long(carch_Crustacean_accum,ci=NA,label.freq = 5)

carcharhinus_fish_plot <- ggplot(data=carch_Fish_accum.long, aes(x = Sites, y = Richness, ymax =  UPR, ymin= LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=2) +
  geom_point(data=subset(carch_Fish_accum.long, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping), alpha=0.2, show.legend=FALSE) + 
  BioR.theme +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Rarefied Samples", y = "Fish Richness", colour = "Carcharhinus Site Status", shape = "Carcharhinus Site Status")

carcharhinus_fish_plot

carcharhinus_crustacean_plot <- ggplot(data=carch_Crustacean_accum.long, aes(x = Sites, y = Richness, ymax =  UPR, ymin= LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=2) +
  geom_point(data=subset(carch_Crustacean_accum.long, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=5) +
  geom_ribbon(aes(colour=Grouping), alpha=0.2, show.legend=FALSE) + 
  BioR.theme +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Rarefied Samples", y = "Crustacean Richness", colour = "Carcharhinus Site Status", shape = "Carcharhinus Site Status")

carcharhinus_crustacean_plot

pdf(carcharhinus_crustacean_plot, "../11_Vegan/vis/carcharhinus_site_crustacean_rarefaction.pdf")


