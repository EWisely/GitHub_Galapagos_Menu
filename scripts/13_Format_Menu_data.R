#Eldridge Wisely
#6-11-2024
#Read in and format Menu phyloseq objects for gllvm workshop exercises

library(phyloseq)
library(tidyverse)

#starting from reading in Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps (this was what I called that phyloseq object I attached to the email just now)----

## Wisconsin transformation =eDNA index, now with traits and environmental variables added!----


Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps<-readRDS("~/Desktop/Galapagos projects/Menu/10_Phyloseq/Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps.RDS")

#Untransformed data

Shark_Menu_taxa_merge_species_named_traits.ps<-readRDS("~/Desktop/Galapagos projects/Menu/10_Phyloseq/Menu_combined2/Menu_combined2_untransformed_traits_added_env_added_taxa_named_Phyloseq.rds")

plot(Shark_Menu_taxa_merge_species_named_traits.ps@phy_tree)


Menu_envvars<-cbind(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps))
write.csv(Menu_envvars, "../results_csvs/Menu_envvars.csv")

#filter out incomplete field measurement columns and coordinates
Menu_envvars<-Menu_envvars%>%
  select(-c(Water_temp,Salinity,Associated_Blood_Samples,Associated_Fecal_Samples,Date.collected,DD_lat,DD_long))
Menu_envvars1<-Menu_envvars[1:15]

#center and scale and set as factors
Menu_envvars <- data.frame(lapply(Menu_envvars1, function(x)if(is.numeric(x)){scale(x)}else{as.factor(x)}))

## Make abundance table with wisconsin-transformed abundances ----
Menu_abund<-as.data.frame(t(cbind(otu_table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps))))
write.csv(Menu_abund, "../results_csvs/Menu_transformed_abundance.csv")
## And lastly, you’ll find the tree is already in the format that ape uses. ----
Menu_tree<-phy_tree(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)

#append Menu_abund to Menu_envvars table

Menu_envvars<-read.csv("../results_csvs/Menu_envvars.csv")
Menu_envvars<-Menu_envvars%>%
  rename(sample=X)
Menu_abund<-rownames_to_column(Menu_abund, var = "sample")
Menu_for_map<-full_join(Menu_envvars, Menu_abund)
write.csv(Menu_for_map, "../results_csvs/Menu_for_map.csv")

### for the traits analysis because there’s so much missing data and the imputation I tried worked VERY poorly, I’ve just filtered the phyloseq object (many species, and some sites were lost) and created new abundance and envvars tables too. ###

#Make a phyloseq object of taxa that has no missing important traits----
Shark_Menu_taxa_merge_species_named_no_missing_traits.eDNAindex.ps<- subset_taxa(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, FoodTroph != "NA")

Shark_Menu_taxa_merge_species_named_no_missing_traits.eDNAindex.ps <- prune_samples(sample_sums(Shark_Menu_taxa_merge_species_named_no_missing_traits.eDNAindex.ps)>0, Shark_Menu_taxa_merge_species_named_no_missing_traits.eDNAindex.ps)

## Make the traits table----
Menu_traits<-as.data.frame(cbind(tax_table(Shark_Menu_taxa_merge_species_named_no_missing_traits.eDNAindex.ps)))

Menu_traits<-Menu_traits%>%
  select(-DepthRangeDeep,-DepthRangeShallow,-AddRems)

##Environmental variables of only the samples that had species with trait data----
Menu_envvars_tr<-cbind(sample_data(Shark_Menu_taxa_merge_species_named_no_missing_traits.eDNAindex.ps))

#filter out incomplete field measurement columns and coordinates
Menu_envvars_tr<-Menu_envvars_tr%>%
  select(-c(Water_temp,Salinity,Associated_Blood_Samples,Associated_Fecal_Samples,Date.collected,DD_lat,DD_long))
Menu_envvars_tr1<-Menu_envvars_tr[1:15]

#center and scale and set as factors
Menu_traits_envvars <- data.frame(lapply(Menu_envvars_tr1, function(x)if(is.numeric(x)){scale(x)}else{as.factor(x)}))



## Make abundance table from "no missing traits" object----
Menu_traits_abund<-as.data.frame(t(cbind(otu_table(Shark_Menu_taxa_merge_species_named_no_missing_traits.eDNAindex.ps))))

## Menu "no missing traits" tree----
Menu_traits_tree<-phy_tree(Shark_Menu_taxa_merge_species_named_no_missing_traits.eDNAindex.ps)

