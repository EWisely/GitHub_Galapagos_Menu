#Eldridge Wisely
#Date: 6-5-24
#Do analysis and visualization of combined BerryCrust and MiFish phyloseq objects 

library("ggplot2")      # graphics
library("tidyverse")
library("dplyr")        # filter and reformat data frames
library(phyloseq)
library(ggsci)
library(metagMisc)
library(ggiraph)
library(microViz)
library(rfishbase)
library(CopernicusMarine)
library(vegan)
library(ggvegan)
library(pairwiseAdonis)
library(Rphylopars)
library(phyloseqGraphTest)

'%ni%' <- Negate("%in%")

#set variables
Primer<-"Menu_combined2" #options Menu_combined has empty PCRs, Menu_combined1 has had the empty PCRs removed, and Menu_combined2 is aggregated by sample (PCRs from the same sample combined), and empty samples removed, phylogenetic tree added.


#Import Phyloseq object----
Menu.ps<-readRDS(file= paste("../09_Metabar_to_Phyloseq/",Primer,"_Phyloseq.rds", sep = ""))

Menu.ps

##filter out passive mid-bay samples because they didn't seem to work well----
Menu1.ps<- ps_filter(Menu.ps, Microhabitat != "passive mid bay")

Menu1.ps

Passive.ps<-ps_filter(Menu.ps, Microhabitat=="passive mid bay")

plot(Menu1.ps@phy_tree)

#Add information about which samples had Carcharhinus or Sphyrna genus detected.----

#Find which samples had hammerheads and which had blacktips----

Blacktips.ps <- subset_taxa(Menu1.ps, genus =="Carcharhinus")
Blacktips.ps <- prune_samples(sample_sums(Blacktips.ps)>0, Blacktips.ps)
Blacktip_sites<-Blacktips.ps@sam_data$Site_Name
Blacktip_samples<-Blacktips.ps@sam_data$sample_id
Blacktip_samples.ps<-prune_samples(x = Menu1.ps,samples = Blacktip_samples)
Blacktip_samples.ps<-prune_taxa(taxa_sums(Blacktip_samples.ps)>0,Blacktip_samples.ps)


#saveRDS(Blacktip_samples.ps, file=paste0("../11_Vegan/Menu_combined2/Blacktip_samples_Phyloseq.rds"))


Hammerheads.ps<-subset_taxa(Menu1.ps,genus=="Sphyrna")
Hammerheads.ps <- prune_samples(sample_sums(Hammerheads.ps)>0, Hammerheads.ps)
Hammerhead_sites<-Hammerheads.ps@sam_data$Site_Name
Hammerhead_samples<-Hammerheads.ps@sam_data$sample_id
Hammerhead_samples.ps<-prune_samples(x = Menu1.ps,samples = Hammerhead_samples)
Hammerhead_samples.ps<-prune_taxa(taxa_sums(Hammerhead_samples.ps)>0,Hammerhead_samples.ps)


#saveRDS(Hammerhead_samples.ps, file=paste0("../11_Vegan/Menu_combined2/Hammerhead_samples_Phyloseq.rds"))

##Make new columns in the sample data for blacktips detected and hammerheads detected with presence absence.----  

Shark_Menu.ps<-Menu1.ps

#Count things before aggregating further
Menu_tax_table<-as.data.frame(cbind(tax_table(Shark_Menu.ps)))

write.csv(Menu_tax_table, file="../10_Phyloseq/Final_menu_taxonomy_table.csv")
nrow(Menu_tax_table)
#246 ASVs ID'ed to at least family level


Menu_tax_table<-rownames_to_column(Menu_tax_table, var = "ASV")

counts<-Menu_tax_table%>%
  dplyr::count(species)
counts


#67 families

#121 taxa
#count the species-level assignments
Menu_species_table<-Menu_tax_table%>%
  filter(family!="NA")%>%
  filter(phylum!="Arthropoda")


distinct<-n_distinct(Menu_species_table$species)
distinct

nrow(Menu_species_table)

Menu_tax_table1<-Menu_tax_table%>%
  group_by(class)%>%
  count()%>%
  ungroup()

Menu_tax_table1
#class
# Groups:   class [4]
# class              n
# <chr>          <int>
#   1 Actinopteri      118
# 2 Chondrichthyes     3
# 3 Malacostraca     122
# 4 Thecostraca        3


#246 ASVs
#71 ID'ed to species level 71 species
#118 ID'ed to at least genus level (101 genuses)
#121 taxa ID'ed to at least family level (67 families)







##Tax_fix to remove empty cells from taxonomy----
Shark_Menu.taxfix.ps<-Shark_Menu.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )


###rename to each ASV separately so that intraspecific variation remains.----
Shark_Menu_ASVs_taxon_named.ps<-tax_rename(Shark_Menu.taxfix.ps,rank = "species")


#count number of ASVs assigned to each taxon
Shark_Menu_ASVs_taxon_named.df<-psmelt(Shark_Menu_ASVs_taxon_named.ps)



#Merge_species but keep taxa that aren't ID'ed to species----
Shark_Menu_taxa_merge.ps <- phyloseq::tax_glom(Shark_Menu.ps, taxrank = "species", NArm = FALSE)


##tax_fix species-merged object to keep higher classifications in the rename----
Shark_Menu_taxa_merge_tax_fix.ps<-Shark_Menu_taxa_merge.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )
##rename with name of species or higher taxonomic level----
Shark_Menu_taxa_merge_species_named.ps<-tax_rename(Shark_Menu_taxa_merge_tax_fix.ps,rank = "species")

Shark_Menu_taxa_merge_species_named.ps

Menu_tax_table_taxamerge<-as.data.frame(cbind(tax_table(Shark_Menu_taxa_merge_species_named.ps)))

Menu_species_table_taxmerge<-Menu_tax_table_taxamerge%>%
  filter(genus!="NA")
distinct<-n_distinct(Menu_species_table_taxmerge$genus)
distinct


#This is the new full dataset!----

saveRDS(Shark_Menu_taxa_merge_species_named.ps, file = "../10_Phyloseq/Shark_Menu_taxa_merge_species_named_ps.RDS")

#Add environmental variables for each sample----

new_sample_data<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named.ps)))
envvars<-read.csv("../000_environmental_data/Menu_sampling_envvars.csv")

envvars<-envvars%>%
  dplyr::rename(sample_id=X)

new_sample_data<-new_sample_data%>%
  select(sample_id,Site_Name,Microhabitat,Tide,Water_temp,Date.collected,DD_lat,DD_long,Salinity,Associated_Blood_Samples,Associated_Fecal_Samples)

new_sample_data<-new_sample_data%>%
  mutate(hammerheads_detected=
           if_else(sample_id %in% Hammerhead_samples, 
                   "hh_detected",
                   "no_hh_detected"))%>%
  mutate(carcharhinus_detected=
           if_else(sample_id %in% Blacktip_samples, 
                   "carch_detected",
                   "no_carch_detected"))%>%
  mutate(blacktip_site=
           if_else(Site_Name %in% Blacktip_sites, 
                   "carch_site",
                   "non-carch_site"))%>%
  mutate(hammerhead_site=
           if_else(Site_Name %in% Hammerhead_sites,
                   "hh_site",
                   "non-hh_site"))


new_sample_data1<-left_join(new_sample_data,envvars)

##center and scale continuous environmental variables----

new_sample_data1<-new_sample_data1%>%
  mutate(scaled_oxygen = scale(oxygen)[,1],
         scaled_Water_temp= scale(Water_temp)[,1],
         scaled_sst= scale(sst)[,1],
         scaled_Salinity= scale(Salinity)[,1],
         scaled_depth= scale(depth)[,1],
         scaled_silicate= scale(silicate)[,1],
         scaled_phosphate= scale(phosphate)[,1],
         scaled_primaryprod= scale(primary_productivity)[,1],
         scaled_nitrate= scale(nitrate)[,1],
         scaled_chloride= scale(chloride)[,1])

new_sample_data2<-column_to_rownames(new_sample_data1, var="sample_id")


sample_data(Shark_Menu_taxa_merge_species_named.ps)<-new_sample_data2
sample_data(Menu1.ps)<-new_sample_data2

sample_names(Shark_Menu_taxa_merge_species_named.ps)

Shark_Menu_taxa_merge_species_named.ps@sam_data$carcharhinus_detected
Shark_Menu_taxa_merge_species_named.ps@sam_data$Tide
Shark_Menu_taxa_merge_species_named.ps@sam_data$scaled_silicate




#Hellinger transform the read counts (square-root of percent abundances) ----
# see: Proper environmental DNA metabarcoding data transformation reveals temporal stability of fish communities in a dendritic river system
# Martin Laporte, Emilie Reny-Nolin, Victoria Chouinard, Cécilia Hernandez, Eric Normandeau, Bérénice Bougas, Caroline Côté, Sonja Behmel, Louis Bernatchez
# First published: 23 June 2021
# https://doi-org.ezproxy4.library.arizona.edu/10.1002/edn3.224
# 

Shark_Menu_ASVs_taxon_named.hell.ps<-phyloseq_standardize_otu_abundance(Shark_Menu_ASVs_taxon_named.ps, method = "hellinger")

Shark_Menu_taxa_merge_species_named.hell.ps<-phyloseq_standardize_otu_abundance(Shark_Menu_taxa_merge_species_named.ps, method = "hellinger")

#Make sure there are no empty samples etc. in the phyloseq object

Shark_Menu_taxa_merge_species_named.hell.ps<-phyloseq_validate(
  Shark_Menu_taxa_merge_species_named.hell.ps,
  remove_undetected = TRUE,
  min_tax_length = 4,
  verbose = TRUE
)

#Get trait fields from fishbase/sealifebase----

##make species_list of only species names to give to fishbase----
Shark_Menu_just_species.ps <- phyloseq::tax_glom(Shark_Menu.ps, taxrank = "species", NArm = TRUE)


species_names<-as.data.frame(tax_table(Shark_Menu_just_species.ps))%>%
  dplyr::select("species")
#50 ASVs were not identifiable to species level. (because there are 121 species plus other ASVs)
#71 unique species found!

sealife_table<-fb_tbl("species","sealifebase") %>% 
  dplyr::mutate(sci_name = paste(Genus, Species)) %>%
  dplyr::filter(sci_name %in% species_names$species) %>% 
  dplyr::select(sci_name, SpecCode, FBname,DemersPelag,DepthRangeShallow, DepthRangeDeep, Vulnerability, Length)

fish_table<-fb_tbl("species","fishbase")%>% 
  dplyr::mutate(sci_name = paste(Genus, Species)) %>%
  dplyr::filter(sci_name %in% species_names$species) %>% 
  dplyr::select(sci_name, SpecCode, FBname,DemersPelag,DepthRangeShallow, DepthRangeDeep, Vulnerability, Length)

Menu_species_table<-rbind(fish_table, sealife_table)

crustacean_ecology_traits<-ecology(species_list = species_names$species, server = "sealifebase")
crustacean_ecology_traits

crustacean_trophic_traits<-ecology(species_list = species_names$species, server = "sealifebase") %>% 
  dplyr::select(dplyr::matches("Troph"),"Species")

crustacean_ecology_traits<-full_join(crustacean_trophic_traits, crustacean_ecology_traits)

fish_ecology_traits<-ecology(species_list = species_names$species, server = "fishbase")
fish_ecology_traits

#get the FoodTroph, FoodSeTroph, Neritic, Intertidal, Oceanic, Epipelagic, Neritic, Mesopelagic, Estuaries, Mangroves, FeedingType, Solitary, and AddRems columns (and add them to the Menu_species_table)
fish_ecology_traits1<-fish_ecology_traits%>%
  dplyr::select(Species,FoodTroph, FoodSeTroph, Neritic, Intertidal, Oceanic, Epipelagic, Neritic, Mesopelagic, Estuaries, Mangroves, FeedingType, Solitary, AddRems)

crustacean_ecology_traits<-crustacean_ecology_traits%>%
  dplyr::rename(Mesopelagic = mesopelagic)

crustacean_ecology_traits1<-crustacean_ecology_traits%>%
  dplyr::select(Species, FoodTroph, FoodSeTroph, Neritic, Intertidal, Oceanic, Epipelagic, Neritic, Mesopelagic, Estuaries, Mangroves, FeedingType, Solitary, AddRems)

Menu_species_table<-Menu_species_table%>%
  dplyr::rename(Species = sci_name)


##Make a traits table for everything that could be ID'ed to species.----
Menu_traits_table<-left_join(Menu_species_table, fish_ecology_traits1, by = "Species")

Menu_traits_table<-left_join(Menu_traits_table, crustacean_ecology_traits1)


#Add the traits to the taxonomy table----

taxonomy_table<-cbind(data.frame(tax_table(Shark_Menu_taxa_merge_species_named.hell.ps)))

#taxonomy_table<-rownames_to_column(taxonomy_table, var="ID")

Menu_traits_table<-Menu_traits_table%>%
  dplyr::rename(species = Species)

taxo_traits_table<-full_join(taxonomy_table,Menu_traits_table)

taxo_traits_table$ID<-paste0(taxo_traits_table$species," 1")

taxo_traits_table<-column_to_rownames(taxo_traits_table, var = "ID")

taxo_traits_matrix<-as.matrix(taxo_traits_table)

##put traits back into the phyloseq object----

Shark_Menu_taxa_merge_species_named_traits.ps<-Shark_Menu_taxa_merge_species_named.ps

tax_table(Shark_Menu_taxa_merge_species_named_traits.ps)<-taxo_traits_matrix

Shark_Menu_taxa_merge_species_named_traits.ps@tax_table
Shark_Menu_taxa_merge_species_named_traits.ps@sam_data

Shark_Menu_taxa_merge_species_named_traits_taxtable<-as.data.frame(cbind(tax_table(Shark_Menu_taxa_merge_species_named_traits.ps)))

write.csv(Shark_Menu_taxa_merge_species_named_traits_taxtable, file="../10_Phyloseq/Shark_Menu_taxa_merge_species_named_traits_taxonomy_table.csv")

#Add trophic information for the fish and crustaceans that rfishbase couldn't find.  
trophic_added<-read.csv("../Traits/Menu final taxonomy table - Sheet1.csv")

trophic_added$ID<-paste0(trophic_added$Taxon.name," 1")

trophic_added<-column_to_rownames(trophic_added, var = "ID")%>%
  dplyr::select(-Taxon.name)

trophic_added_matrix<-as.matrix(trophic_added)

##put traits back into the phyloseq object----
tax_table(Shark_Menu_taxa_merge_species_named_traits.ps)<-trophic_added_matrix

Shark_Menu_taxa_merge_species_named_traits.ps@tax_table
Shark_Menu_taxa_merge_species_named_traits.ps@sam_data


#This is the very new full dataset with traits and environmental variables!----

Shark_Menu_taxa_merge_species_named_traits.ps

saveRDS(Shark_Menu_taxa_merge_species_named_traits.ps, file = paste0("../10_Phyloseq/",Primer,"/",Primer,"_untransformed_traits_added_env_added_taxa_named_Phyloseq.rds"))




#Hellinger transform the read counts (square-root of percent abundances) ----
# see: Proper environmental DNA metabarcoding data transformation reveals temporal stability of fish communities in a dendritic river system
# Martin Laporte, Emilie Reny-Nolin, Victoria Chouinard, Cécilia Hernandez, Eric Normandeau, Bérénice Bougas, Caroline Côté, Sonja Behmel, Louis Bernatchez
# First published: 23 June 2021
# https://doi-org.ezproxy4.library.arizona.edu/10.1002/edn3.224
# 

Shark_Menu_taxa_merge_species_named_traits.hell.ps<-phyloseq_standardize_otu_abundance(Shark_Menu_taxa_merge_species_named_traits.ps, method = "hellinger")

Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps<-phyloseq_standardize_otu_abundance(Shark_Menu_taxa_merge_species_named_traits.ps, method = "wisconsin")

# #merge by trophic information
# Shark_Menu_Troph_merge.eDNAindex.ps<- phyloseq::tax_glom(Shark_Menu_taxa_merge_species_named_traits.ps, taxrank = "FeedingType", NArm = TRUE)
  
  
#Make sure there are no empty samples etc. in the phyloseq object

Shark_Menu_taxa_merge_species_named_traits.hell.ps<-phyloseq_validate(
  Shark_Menu_taxa_merge_species_named_traits.hell.ps,
  remove_undetected = TRUE,
  min_tax_length = 4,
  verbose = TRUE
)

#this now has NAs in the taxonomy table because of the missing traits data.

saveRDS(Shark_Menu_taxa_merge_species_named_traits.hell.ps, file = paste0("../10_Phyloseq/",Primer,"/",Primer,"_hellinger_transformed_traits_added_taxa_named_Phyloseq.rds"))

Shark_Menu_taxa_merge_species_named_traits.pa.ps<-phyloseq_standardize_otu_abundance(Shark_Menu_taxa_merge_species_named_traits.hell.ps, method = "pa")

saveRDS(Shark_Menu_taxa_merge_species_named_traits.pa.ps, file = paste0("../10_Phyloseq/",Primer,"/",Primer,"_pa_transformed_traits_added_taxa_named_Phyloseq.rds"))

###Also for the ASV phyloseq object----
ASV_taxonomy_table<-cbind(data.frame(tax_table(Shark_Menu_ASVs_taxon_named.hell.ps)))

ASV_taxonomy_table<-rownames_to_column(ASV_taxonomy_table, var="ID")

#Menu_traits_table<-Menu_traits_table%>%
#  dplyr::rename(species = Species)

ASV_traits_table<-left_join(ASV_taxonomy_table,Menu_traits_table, by ="species")

ASV_traits_table<-column_to_rownames(ASV_traits_table, var = "ID")

ASV_traits_matrix<-as.matrix(ASV_traits_table)

Shark_Menu_ASVs_taxon_named_traits.hell.ps<-Shark_Menu_ASVs_taxon_named.hell.ps

tax_table(Shark_Menu_ASVs_taxon_named_traits.hell.ps)<-ASV_traits_matrix

Shark_Menu_ASVs_taxon_named_traits.hell.ps@tax_table


Shark_Menu_ASVs_taxon_named_traits.hell.ps@sam_data$Microhabitat

Shark_Menu_ASVs_taxon_named_traits.hell.ps@sam_data$Tide

#Make a phyloseq object for gllvm----
##Make a phyloseq object of taxa that has DemersPelag info----
Shark_Menu_taxa_merge_species_named_DemersPelag.hell.ps<- subset_taxa(Shark_Menu_taxa_merge_species_named_traits.hell.ps, DemersPelag != "NA")

Shark_Menu_taxa_merge_species_named_DemersPelag.hell.ps@tax_table

##Make a phyloseq object of taxa that has no missing traits----
Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps<- subset_taxa(Shark_Menu_taxa_merge_species_named_DemersPelag.hell.ps, FoodTroph != "NA")

Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps@tax_table

#has all traits complete except DepthRangeShallow and DepthRangeDeep, I'll just delete those columns for the gllvm

##Make an eDNAindex phyloseq object of taxa that has DemersPelag info----
Shark_Menu_taxa_merge_species_named_DemersPelag.eDNAindex.ps<- subset_taxa(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, DemersPelag != "NA")

Shark_Menu_taxa_merge_species_named_DemersPelag.eDNAindex.ps@tax_table

##Make an eDNAindex phyloseq object of taxa that has no missing traits----
Shark_Menu_taxa_merge_species_named_ALL_traits.eDNAindex.ps<- subset_taxa(Shark_Menu_taxa_merge_species_named_DemersPelag.eDNAindex.ps, FoodTroph != "NA")

Shark_Menu_taxa_merge_species_named_ALL_traits.eDNAindex.ps@tax_table



##Remove samples that are missing temperature data----

Shark_Menu_taxa_merge_species_named_ALL_traits_temp.eDNAindex.ps <- Shark_Menu_taxa_merge_species_named_ALL_traits.eDNAindex.ps %>% ps_filter(Water_temp !="NA")

Shark_Menu_taxa_merge_species_named_ALL_traits_temp.eDNAindex.ps



#gllvm using eDNAindex transformed phyloseq----
#https://jenniniku.github.io/gllvm/articles/vignette3.html
library(gllvm)
#Package **mvabund** is loaded with **gllvm** so just load with a function `data()`.
data("spider")
# more info: 
#?spider
fitp <- gllvm(y = spider$abund, family = poisson(), num.lv = 2)
fitnb <- gllvm(y = spider$abund, family = "negative.binomial", num.lv = 2)
AIC(fitp)
## [1] 1761.655
AIC(fitnb)
## [1] 1561.361

#model it off of the demo data
menu_gllvm_data<-spider

#the traits were only available for a few species. The list was built off of this phyloseq object: Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps 

menu_gllvm_data$trait
menu_gllvm_trait<-cbind(data.frame(tax_table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))

menu_gllvm_trait<-rownames_to_column(menu_gllvm_trait, var="ID")
colnames(menu_gllvm_trait)
menu_gllvm_trait<-menu_gllvm_trait%>%
  dplyr::select(ID, FoodTroph, trophic_role)
menu_gllvm_trait

menu_gllvm_data$trait<-menu_gllvm_trait

menu_gllvm_data$abund
traits_otu_matrix<-as(otu_table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps), "matrix")
traits_otu_table<-as.data.frame(traits_otu_matrix)
traits_otu_table1<-rownames_to_column(traits_otu_table, var="taxname")
taxnames<-traits_otu_table1$taxname
taxnames<-gsub(pattern = " 1", replacement = "",taxnames)
traits_otu_table1$taxname<-taxnames
traits_otu_table<-column_to_rownames(traits_otu_table1, var = "taxname")
menu_gllvm_data$abund<-t(traits_otu_table)
menu_gllvm_data$abund


menu_gllvm_data$x
menu_gllvm_env_vars<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))
menu_gllvm_env_vars<-menu_gllvm_env_vars%>%
  dplyr::select(Site_Name,Microhabitat,Tide, scaled_sst,scaled_depth,scaled_silicate,scaled_phosphate,scaled_primaryprod,scaled_nitrate,scaled_chloride,scaled_oxygen)
menu_gllvm_env_vars

menu_gllvm_data$x<-menu_gllvm_env_vars

saveRDS(menu_gllvm_data, "../11_Vegan/menu_gllvm_data_eDNAindex.rds")



fitp <- gllvm(y = menu_gllvm_data$abund, family = poisson(), num.lv = 2)
fitnb <- gllvm(y = menu_gllvm_data$abund, family = "negative.binomial", num.lv = 2)
AIC(fitp)
## [1] 704.4989 eDNAindex 395.7438eDNAindex with traits 1301.177 eDNAindex all species
AIC(fitnb)
## [1] 772.4951 463.8158eDNAindex with traits 1839.015 eDNAindex all species

par(mfrow = c(1,2))
plot(fitnb, which = 1:2)

fitLAp <- gllvm(y=menu_gllvm_data$abund, family = poisson(), method = "LA", num.lv = 2)
fitLAnb <- gllvm(y=menu_gllvm_data$abund, family = "negative.binomial", method = "LA", num.lv = 2)
fitLAzip <- gllvm(y=menu_gllvm_data$abund, family = "ZIP", method = "LA", num.lv = 2)
AIC(fitLAp)
AIC(fitLAnb)
AIC(fitLAzip)

#fitLAp and fitp are tied for the lowest AIC so I'll go with fitp

plot(fitp)
plot.new()
dev.off()
ordiplot(fitp, biplot = TRUE)
abline(h = 0, v = 0, lty=2)
dev.off()


# Parameters:
coef(fitp)

getLV(fitp)

# Standard errors for parameters:
fitp$sd

#add environmental variables
fitpx <- gllvm(y = menu_gllvm_data$abund, X = menu_gllvm_data$x, family = poisson(), seed = 123, num.lv = 2)
fitpx

# Call: #for hellinger-transformed data
#   gllvm(y = menu_gllvm_data$abund, X = menu_gllvm_data$x, family = poisson(), 
#         num.lv = 2, seed = 123)
# family: 
#   [1] "poisson"
# method: 
#   [1] "VA"
# 
# log-likelihood:  -458.2023 
# Residual degrees of freedom:  4720 
# AIC:  5754.405 
# AICc:  8235.43 
# BIC:  22380.98 


# Call: for the eDNA index transformed data
#   gllvm(y = menu_gllvm_data$abund, X = menu_gllvm_data$x, family = poisson(), 
#         num.lv = 2, seed = 123)
# family: 
#   [1] "poisson"
# method: 
#   [1] "VA"
# 
# log-likelihood:  -197.6456 
# Residual degrees of freedom:  4720 
# AIC:  5233.291 
# AICc:  7714.317 
# BIC:  21859.87 

plot.new()
ordiplot(fitpx, biplot = TRUE)
abline(h = 0, v = 0, lty=2)

#coefplot(fitpx, mfrow = c(3,3), cex.ylab = 0.8,cex.xlab = 1.5, mar = c(4,10,2,1))#for hellinger version
coefplot(fitpx, mfrow = c(3,3), cex.ylab = 0.8,cex.xlab = 1.5, mar = c(4,10,2,1))#for eDNAindex version
#coefplot(fitpx, mfrow = c(1,3), cex.ylab = 0.8,cex.xlab = 1.5)


cv<-getResidualCov(fitpx)

cr <- getResidualCor(fitpx)
#install.packages("corrplot")
library(corrplot)
plot.new()
dev.off()
corrplot(cr, diag = FALSE, type = "lower", method = "square", tl.srt = 25, tl.cex = 0.4)#for eDNAindex transformation
#corrplot(cr, diag = FALSE, type = "lower", method = "square", tl.srt = 25, tl.cex = 0.4) #for hellinger transformed object

#for all species residual correlations after all environmental variables are accounted for
corrplot(cr, diag = FALSE, type = "lower", method = "circle", tl.srt = 25, tl.cex = 0.4)


#custom function to calculate matrix of p=values from the data: 
#copied from: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram#google_vignette

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(cr)
head(p.mat[, 1:5])

# mark out the ones with p-values greater than 0.05
corrplot(cr, type="lower", 
         p.mat = p.mat, sig.level = 0.05,tl.cex = 0.4)

# Leave blank on no significant coefficient p<0.05
corrplot(cr, type="lower", 
         p.mat = p.mat, sig.level = 0.05, insig = "blank",tl.cex = 0.4)

#Best correlation plot----  p<0.05
corrplot(cr, type="lower", order = "FPC",
         #addCoef.col = "black",number.cex = .3, # Add coefficient of correlation
         tl.col="black", tl.srt=45,tl.cex = .45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=TRUE,
         #title = "Correlation of Fish and Crustacean eDNA index values after removing environmental effects" 
)

#make a list of taxa in species complex1 and species complex2 with greater than 0.6 correlation coefficients.----  

species_species_correllation_coeffs<-rownames_to_column(as.data.frame(cr), var="taxaname")
Carcharhinus_species_complex_gllvm<-species_species_correllation_coeffs%>%
  filter(`Carcharhinus genus`>0.6)

Sphyrna_lewini_species_complex_gllvm<-species_species_correllation_coeffs%>%
  filter(`Sphyrna lewini`>0.6)

Carch_associated_speciescomplex<-Carcharhinus_species_complex_gllvm$taxaname
Slewini_associated_speciescomplex<-Sphyrna_lewini_species_complex_gllvm$taxaname

length(Carch_associated_speciescomplex)
# 61 species (eDNA index transformation)
length(Slewini_associated_speciescomplex)
# 58 species (eDNA index transformation)


#Add membership to each species complex to the taxonomy table----

taxonomy_table<-cbind(data.frame(tax_table(Shark_Menu_taxa_merge_species_named_traits.ps)))


taxo_sc_table<-taxonomy_table%>%
  mutate(Slew_complex_gllvm_eDNAindex=ifelse(species %in% Slewini_associated_speciescomplex,
                             1,0),
         Carch_complex_gllvm_eDNAindex=ifelse(species %in% Carch_associated_speciescomplex,
                              1,0))


taxo_sc_matrix<-as.matrix(taxo_sc_table)

##put species complex membership info back into the phyloseq object----

#Shark_Menu_taxa_merge_species_named_traits.ps<-Shark_Menu_taxa_merge_species_named.ps

#tax_table(Shark_Menu_taxa_merge_species_named_traits.ps)<-taxo_sc_matrix

Shark_Menu_taxa_merge_species_named_traits.ps@tax_table
Shark_Menu_taxa_merge_species_named_traits.ps@sam_data


#gllvm using hellinger transformed phyloseq----
#https://jenniniku.github.io/gllvm/articles/vignette3.html
library(gllvm)
#Package **mvabund** is loaded with **gllvm** so just load with a function `data()`.
data("spider")
# more info: 
#?spider
fitp <- gllvm(y = spider$abund, family = poisson(), num.lv = 2)
fitnb <- gllvm(y = spider$abund, family = "negative.binomial", num.lv = 2)
AIC(fitp)
## [1] 1761.655
AIC(fitnb)
## [1] 1561.361

#model it off of the demo data
menu_gllvm_data<-spider

#the traits were only available for a few species. The list was built off of this phyloseq object: Shark_Menu_taxa_merge_species_named_traits.hell.ps (because doing comparisons between samples, I don't want to have transformed across samples, just within samples.)

menu_gllvm_data$trait
menu_gllvm_trait<-cbind(data.frame(tax_table(Shark_Menu_taxa_merge_species_named_traits.hell.ps)))

menu_gllvm_trait<-rownames_to_column(menu_gllvm_trait, var="ID")
colnames(menu_gllvm_trait)
menu_gllvm_trait<-menu_gllvm_trait%>%
  select(ID, FoodTroph, Epipelagic, Neritic, Mesopelagic, Solitary, DemersPelag, Intertidal, Estuaries, Oceanic, Mangroves)
menu_gllvm_trait

menu_gllvm_data$trait<-menu_gllvm_trait

menu_gllvm_data$abund
traits_otu_matrix<-as(otu_table(Shark_Menu_taxa_merge_species_named_traits.hell.ps), "matrix")
traits_otu_table<-as.data.frame(traits_otu_matrix)
traits_otu_table1<-rownames_to_column(traits_otu_table, var="taxname")
taxnames<-traits_otu_table1$taxname
taxnames<-gsub(pattern = " 1", replacement = "",taxnames)
traits_otu_table1$taxname<-taxnames
traits_otu_table<-column_to_rownames(traits_otu_table1, var = "taxname")
menu_gllvm_data$abund<-t(traits_otu_table)
menu_gllvm_data$abund


menu_gllvm_data$x
menu_gllvm_env_vars<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)))
menu_gllvm_env_vars<-menu_gllvm_env_vars%>%
  select(Site_Name,Microhabitat,Tide, scaled_sst,scaled_depth,scaled_silicate,scaled_phosphate,scaled_primaryprod,scaled_nitrate,scaled_chloride,scaled_oxygen)
menu_gllvm_env_vars

menu_gllvm_data$x<-menu_gllvm_env_vars

saveRDS(menu_gllvm_data, "../11_Vegan/menu_gllvm_data_hell.rds")



fitp <- gllvm(y = menu_gllvm_data$abund, family = poisson(), num.lv = 2)
fitnb <- gllvm(y = menu_gllvm_data$abund, family = "negative.binomial", num.lv = 2)
AIC(fitp)
## [1] 704.4989 hell 395.7438eDNAindex with traits 1301.177 eDNAindex all species
AIC(fitnb)
## [1] 772.4951 463.8158eDNAindex with traits 1839.015 eDNAindex all species

par(mfrow = c(1,2))
plot(fitnb, which = 1:2)

fitLAp <- gllvm(y=menu_gllvm_data$abund, family = poisson(), method = "LA", num.lv = 2)
fitLAnb <- gllvm(y=menu_gllvm_data$abund, family = "negative.binomial", method = "LA", num.lv = 2)
fitLAzip <- gllvm(y=menu_gllvm_data$abund, family = "ZIP", method = "LA", num.lv = 2)
AIC(fitLAp)
AIC(fitLAnb)
AIC(fitLAzip)

#fitLAp and fitp are tied for the lowest AIC so I'll go with fitp

plot(fitp)
plot.new()
dev.off()
ordiplot(fitp, biplot = TRUE)
abline(h = 0, v = 0, lty=2)
dev.off()


# Parameters:
coef(fitp)

getLV(fitp)

# Standard errors for parameters:
fitp$sd

#add environmental variables
fitpx <- gllvm(y = menu_gllvm_data$abund, X = menu_gllvm_data$x, family = poisson(), seed = 123, num.lv = 2)
fitpx

# Call: #for hellinger-transformed data
#   gllvm(y = menu_gllvm_data$abund, X = menu_gllvm_data$x, family = poisson(), 
#         num.lv = 2, seed = 123)
# family: 
#   [1] "poisson"
# method: 
#   [1] "VA"
# 
# log-likelihood:  -458.2023 
# Residual degrees of freedom:  4720 
# AIC:  5754.405 
# AICc:  8235.43 
# BIC:  22380.98 

plot.new()
ordiplot(fitpx, biplot = TRUE)
abline(h = 0, v = 0, lty=2)

coefplot(fitpx, mfrow = c(3,3), cex.ylab = 0.8,cex.xlab = 1.5, mar = c(4,10,2,1))#for hellinger version
#coefplot(fitpx, mfrow = c(2,3), cex.ylab = 0.8,cex.xlab = 1.5, mar = c(4,9,2,1))#for eDNAindex version
#coefplot(fitpx, mfrow = c(1,3), cex.ylab = 0.8,cex.xlab = 1.5)


cv<-getResidualCov(fitpx)

cr <- getResidualCor(fitpx)
#install.packages("corrplot")
library(corrplot)
plot.new()
dev.off()
#corrplot(cr, diag = FALSE, type = "lower", method = "square", tl.srt = 25, tl.cex = 0.4)#for eDNAindex transformation
corrplot(cr, diag = FALSE, type = "lower", method = "square", tl.srt = 25, tl.cex = 0.4) #for hellinger transformed object

#for all species residual correlations after all environmental variables are accounted for
corrplot(cr, diag = FALSE, type = "lower", method = "circle", tl.srt = 25, tl.cex = 0.4)


#custom function to calculate matrix of p=values from the data: 
#copied from: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram#google_vignette

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(cr)
head(p.mat[, 1:5])

# mark out the ones with p-values greater than 0.05
corrplot(cr, type="lower", 
         p.mat = p.mat, sig.level = 0.05,tl.cex = 0.4)

# Leave blank on no significant coefficient p<0.05
corrplot(cr, type="lower", 
         p.mat = p.mat, sig.level = 0.05, insig = "blank",tl.cex = 0.4)

#Best correlation plot----  p<0.05
corrplot(cr, type="lower", order = "FPC",
         #addCoef.col = "black",number.cex = .3, # Add coefficient of correlation
         tl.col="black", tl.srt=45,tl.cex = .45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank",
         # hide correlation coefficient on the principal diagonal
         diag=TRUE,
         #title = "Correlation of Fish and Crustacean Hellinger transformed values after removing environmental effects" 
)

#make a list of taxa in species complex1 and species complex2 with greater than 0.6 correlation coefficients.----  

species_species_correllation_coeffs<-rownames_to_column(as.data.frame(cr), var="taxaname")
Carcharhinus_species_complex_gllvm<-species_species_correllation_coeffs%>%
  filter(`Carcharhinus genus`>0.6)

Sphyrna_lewini_species_complex_gllvm<-species_species_correllation_coeffs%>%
  filter(`Sphyrna lewini`>0.6)

Carch_associated_speciescomplex<-Carcharhinus_species_complex_gllvm$taxaname
Slewini_associated_speciescomplex<-Sphyrna_lewini_species_complex_gllvm$taxaname


length(Carch_associated_speciescomplex)
#60 species (hellinger transformation)
length(Slewini_associated_speciescomplex)
#51 species (hellinger transformation)

#Add membership to each species complex to the taxonomy table----

taxonomy_table<-cbind(data.frame(tax_table(Shark_Menu_taxa_merge_species_named_traits.ps)))

taxo_sc_table<-taxonomy_table%>%
  mutate(Slew_complex_hell=ifelse(species %in% Slewini_associated_speciescomplex,
                                1,0),
         Carch_complex_hell=ifelse(species %in% Carch_associated_speciescomplex,
                              1,0))


taxo_sc_matrix<-as.matrix(taxo_sc_table)

##put species complex membership info back into the phyloseq object----

#Shark_Menu_taxa_merge_species_named_traits.ps<-Shark_Menu_taxa_merge_species_named.ps

#tax_table(Shark_Menu_taxa_merge_species_named_traits.ps)<-taxo_sc_matrix

Shark_Menu_taxa_merge_species_named_traits.ps@tax_table
Shark_Menu_taxa_merge_species_named_traits.ps@sam_data

#This is the very new full dataset with species complex from gllvm, traits, and environmental variables!----

Shark_Menu_taxa_merge_species_named_traits.ps

#saveRDS(Shark_Menu_taxa_merge_species_named_traits.ps, file = paste0("../10_Phyloseq/",Primer,"/",Primer,"_untransformed_traits_added_env_added_taxa_named_Phyloseq.rds"))




# #Extra gllvm stuff that I don't think I need
# 
# X <- menu_gllvm_data$x%>%
#   select(scaled_sst,scaled_primaryprod,scaled_phosphate,scaled_oxygen)
# fitpx1 <- gllvm(menu_gllvm_data$abund, X, family = poisson(), num.lv = 1)
# coefplot(fitpx1, mfrow = c(1,2), cex.ylab = 0.8)
# #ggsave(filename = "../11_Vegan/Menu_combined2/Tide_gllvm_coeff.jpg")
# coefplot(fitpx1, mfrow = c(2,2), cex.ylab = 0.8,mar=c(2,12,1,1),cex.xlab = 2)
# coefplot(fitpx1, mfrow = c(2,2), cex.ylab = 0.8,cex.xlab = 1.5,mar=c(4,10,2,1))
# #ggsave(filename = "../11_Vegan/Menu_combined2/Microhabitat_and_Tide_gllvm_coeff.jpg")
# 
# #the coefplots all crossed 0 so are not significant for just the ones with traits data hellinger transformed
# 
# ##Include trait information in the gllvm ----
# 
# #https://jenniniku.github.io/gllvm/articles/vignette1.html#incorporating-functional-traits-into-fourth-corner-models
# criteria <- NULL
# for(i in 1:5){
#   fiti <- gllvm(menu_gllvm_data$abund, menu_gllvm_data$x, family = poisson(), num.lv = i, sd.errors = FALSE,
#                 formula = ~ Site_Name + Microhabitat + Water_temp, seed = 1234)
#   criteria[i] <- summary(fiti)$AICc
#   names(criteria)[i] = i
# }
# 
# criteria
# #1 latent variable had the lowest AIC criteria number, so we'll use that.
# 
# fit_env <- gllvm(menu_gllvm_data$abund, menu_gllvm_data$x, family = poisson(), num.lv = 1,
#                  formula = ~Site_Name+Tide+Microhabitat, seed = 1234)
# dev.off()
# coefplot(fit_env, cex.ylab = 0.7, mfrow=c(1,3),mar = c(4,9,2,4))
# #ggsave(filename = "../11_Vegan/GLLVM_species_microhabitat_coeffs.jpg")
# 
# # Fit GLLVM without environmental variables and 1 latent variable:
# fit1lv <- gllvm(menu_gllvm_data$abund, family = poisson(), num.lv = 1, seed = 1234)
# 
# # Correlation matrix
# library(gclus)
# cr0 <- getResidualCor(fit1lv)
# dev.off()
# corrplot(cr0[order.single(cr0), order.single(cr0)], diag = FALSE, type = "lower", 
#          method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")
# 
# ##correlation plots from gllvm----
# # Residual correlation matrix:
# cr <- getResidualCor(fit_env)
# library(corrplot); library(gclus)
# dev.off()
# corrplot(cr[order.single(cr), order.single(cr)], diag = FALSE, type = "lower", 
#          method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")
# 
# 

###add traits finally----

# # #Error in colMeans(x, na.rm = TRUE) : 'x' must be numeric
# # In addition: Warning message:
# #   In gllvm(menu_gllvm_data$abund, menu_gllvm_data$x, menu_gllvm_data$trait,  :
# #              There are rows full of zeros in y.
# 
# fit_4th <- gllvm(menu_gllvm_data$abund, menu_gllvm_data$x, menu_gllvm_data$trait, family = poisson(), num.lv = 1,
#                  formula = menu_gllvm_data$abund ~ (Site_Name + Microhabitat + Tide) +
#                    (Site_Name + Microhabitat + Tide) : (DemersPelag +
#                                                                         Oceanic + Mangroves+ Neritic+Epipelagic+FeedingType), seed = 123,
#                  row.eff = "random", control.start =list(n.init = 3, jitter.var = 0.01),
#                  randomX = ~ Site_Name + Microhabitat + Tide)
# 
# library(lattice)
# coefplot(fit_4th, mar = c(4, 11, 1, 1), cex.ylab = 0.8)
# fourth <- fit_4th$fourth.corner
# a <- 1.5
# colort <- colorRampPalette(c("blue", "white", "red"))
# plot.4th <- levelplot((as.matrix(fourth)), xlab = "Environmental Variables",
#                       ylab = "Species traits", col.regions = colort(100), cex.lab = 1.3,
#                       at = seq(-a, a, length = 100), scales = list(x = list(rot = 45)))
# plot.4th


#Permanovas-----
##calculate Landes distances----

pairwise_landes_distances<-as.matrix(print(phyloseq::distance(Shark_Menu_taxa_merge_species_named_traits.hell.ps, method="l", binary=F, type= "samples", upper=TRUE, diag=TRUE)))

## Site_Name*Microhabitat----
set.seed(200)
landes_permanova<- adonis2(as.dist(pairwise_landes_distances)~Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Site_Name*Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Microhabitat)

landes_permanova
write_delim(landes_permanova, file = "../11_Vegan/Menu_combined2/landes_permanova_results.txt")

#Site is significant p=0.001, Microhabitat p=0.022 less so.
Shark_Menu_taxa_merge_species_named.hell.ps@sam_data$Microhabitat



Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps<-readRDS("~/Desktop/Galapagos projects/Menu/10_Phyloseq/Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps.RDS")

#following https://github.com/pmartinezarbizu/pairwiseAdonis

landes_dist_matrix <- phyloseq::distance(Shark_Menu_taxa_merge_species_named_traits.hell.ps, method ="l")

set.seed(200)
vegan::adonis2(landes_dist_matrix ~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Site_Name+Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Microhabitat+Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Tide
+Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$scaled_sst)
# vegan::adonis2(formula = landes_dist_matrix ~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Site_Name + Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Microhabitat + Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Tide + Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$scaled_sst)
# Df SumOfSqs      R2      F Pr(>F)    
# phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Site_Name  6   1076.4 0.20745 2.4430  0.001 ***
#   Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Microhabitat             2    278.2 0.05362 1.8942  0.027 *  
#   Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Tide                     1     83.0 0.01600 1.1308  0.319    
# Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$scaled_sst               1    226.2 0.04360 3.0806  0.007 ** 
#   Residual                                                                            48   3525.0 0.67933                  
# Total                                                                               58   5188.9 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

site_micro_tide_sst_permanova<-vegan::adonis2(landes_dist_matrix ~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Site_Name+Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Microhabitat+Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$Tide
                                              +Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data$scaled_sst)

site_micro_tide_sst_permanova

###  Faith's site microhabitat tide water temp
site_micro_tide_sst_permanova_hillsfaithjacc<-vegan::adonis2(Menu_taxmerged_Faiths_Hill_Jaccard_dissim~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)$Site_Name+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Microhabitat+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Tide+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$scaled_Water_temp, na.action = na.omit)

site_micro_tide_sst_permanova_hillsfaithjacc

site_lat_long_micro_tide_sst_permanova_hillsfaithjacc<-vegan::adonis2(Menu_taxmerged_Faiths_Hill_Jaccard_dissim~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)$Site_Name+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$DD_lat+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$DD_long+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Microhabitat+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Tide+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$scaled_Water_temp, na.action = na.omit)

site_lat_micro_tide_sst_permanova_hillsfaithjacc


finalotutable<-as.data.frame(t(cbind(otu_table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps))))

specieswithmostcontribution<-simper(finalotutable,ordered=TRUE, permutations = (999))

specieswithmostcontribution

sitenamesimper<-simper(finalotutable, group = Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Site_Name, ordered=TRUE, permutations = (999),dist=Menu_taxmerged_Faiths_Hill_Jaccard_dissim) #no difference with "l" bray, or jaccard

sitenamesimper


hhsitesimper<-simper(finalotutable, group = Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$hammerhead_site, ordered=TRUE, permutations = (999),dist=Menu_taxmerged_Faiths_Hill_Jaccard_dissim)#no diff with Landes or jaccard

hhsitesimper


btsitesimper<-simper(finalotutable, group = Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$blacktip_site, ordered=TRUE, permutations = (999),dist=Menu_taxmerged_Faiths_Hill_Jaccard_dissim)#no diff with Landes or jaccard

btsitesimper

###Just site_name----
set.seed(2)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Site_Name)

#                                   pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1               Puerto Grande vs La Seca  1 121.70783 2.0612459 0.07616966   0.049      1.000    
# 2         Puerto Grande vs Rosa Blanca 1  1 325.58091 4.4087088 0.12108940   0.001      0.021   .
# 3            Puerto Grande vs La Tortuga  1 112.98899 1.9052590 0.07354719   0.076      1.000    
# 4           Puerto Grande vs Cerro Brujo  1 206.64575 3.3079001 0.09641803   0.006      0.126    
# 5         Puerto Grande vs Rosa Blanca 2  1 345.30598 5.4321101 0.17282041   0.002      0.042   .
# 6  Puerto Grande vs Cartago Bay, Isabela  1 182.37361 3.4388761 0.12532860   0.001      0.021   .
# 7               La Seca vs Rosa Blanca 1  1 154.28523 1.3982761 0.09711414   0.213      1.000    
# 8                  La Seca vs La Tortuga  1 100.12202 1.0123135 0.16837337   0.383      1.000    
# 9                 La Seca vs Cerro Brujo  1 106.25893 1.2652226 0.09537892   0.235      1.000    
# 10              La Seca vs Rosa Blanca 2  1 149.67917 1.4471992 0.17132297   0.136      1.000    
# 11       La Seca vs Cartago Bay, Isabela  1 211.83631 3.0788297 0.38109848   0.066      1.000    
# 12           Rosa Blanca 1 vs La Tortuga  1  59.60011 0.5176764 0.04135563   0.799      1.000    
# 13          Rosa Blanca 1 vs Cerro Brujo  1 128.90963 1.2925984 0.06369802   0.253      1.000    
# 14        Rosa Blanca 1 vs Rosa Blanca 2  1 116.14460 1.0093091 0.06724554   0.395      1.000    
# 15 Rosa Blanca 1 vs Cartago Bay, Isabela  1 244.88582 2.3870716 0.16591783   0.037      0.777    
# 16             La Tortuga vs Cerro Brujo  1  47.38782 0.5458521 0.04727690   0.858      1.000    
# 17           La Tortuga vs Rosa Blanca 2  1  46.61667 0.4167577 0.06494833   0.799      1.000    
# 18    La Tortuga vs Cartago Bay, Isabela  1 124.83333 1.7149399 0.30008013   0.200      1.000    
# 19          Cerro Brujo vs Rosa Blanca 2  1 154.70833 1.6980462 0.11552870   0.068      1.000    
# 20   Cerro Brujo vs Cartago Bay, Isabela  1 223.69551 3.0587670 0.21757008   0.002      0.042   .
# 21 Rosa Blanca 2 vs Cartago Bay, Isabela  1 232.36667 2.6778923 0.30858787   0.070      1.000  

set.seed(2)
site_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Site_Name)
write_delim(site_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/site_landes_pairwise_adonis_results_hell.txt")



set.seed(200)
site_taxmerged_Hill_Shannonexp_Jaccard_pairwise_adonis<-pairwise.adonis(Menu_taxmerged_Hill_Shannonexp_Jaccard_dissim, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)$Site_Name)
write_delim(site_taxmerged_Hill_Shannonexp_Jaccard_pairwise_adonis, "../11_Vegan/Menu_combined2/site_taxmerged_Hill_Shannonexp_Jaccard_pairwise_adonis_eDNAindex.txt")




###just microhabitat----
set.seed(200)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Microhabitat)

# > ###just microhabitat----
# > set.seed(200)
# > pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Microhabitat)
# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 mangroves vs middle of bay  1  55.40317 0.5803425 0.01245896   0.800      1.000    
# 2         mangroves vs beach  1 203.61407 2.9285037 0.09468624   0.011      0.033   .
# 3     middle of bay vs beach  1 180.98928 1.9667519 0.04920970   0.035      0.105 

set.seed(200)
microhabitat_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Microhabitat)
write_delim(microhabitat_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/microhabitat_landes_pairwise_adonis_results_hell.txt")

#mangroves vs. beach were significant at and middle of bay vs beach was significant

set.seed(200)
microhabitat_Shannonexpjaccard_pairwise_adonis<-pairwise.adonis(Menu_taxmerged_Hill_Shannonexp_Jaccard_dissim, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)$Microhabitat)
write_delim(microhabitat_Shannonexpjaccard_pairwise_adonis, "../11_Vegan/Menu_combined2/microhabitat_shannonjaccard_pairwise_adonis_results_hell.txt")


##Check for differences between hammerhead sites and non-hammerhead sites.----
#First remove hammerheads so they don't confound the calculation.
Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps <- Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps %>% subset_taxa(species != "Sphyrna lewini")

landes_dist_matrix_rm_hh<-phyloseq::distance(Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps, method = "l")
set.seed(2)
pairwise.adonis(landes_dist_matrix_rm_hh, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps)$hammerhead_site)

# hellinger
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 hh_site vs non-hh_site  1  443.6926 5.385417 0.08632493   0.001      0.001  **

#eDNAindex
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 hh_site vs non-hh_site  1  452.8712 5.450488 0.08727695   0.001      0.001  **

# pairwise.adonis(Menu_taxmerged_rmhh_Hill_Jaccard_dissim, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.hell.ps)$hammerhead_site)
# pairs Df  SumsOfSqs    F.Model          R2 p.value p.adjusted sig
# 1 hh_site vs non-hh_site  1 -0.0227335 -0.6081797 -0.01078489   0.988      0.988    


set.seed(2)
hhsite_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix_rm_hh, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.hell.ps)$hammerhead_site)
write_delim(hhsite_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/hhsite_landes_pairwise_adonis_results_hell.txt")
#significant 0.001

set.seed(200)
pairwise.adonis(Menu_taxmerged_rmhh_Shannons_Hill_Jaccard_dissim, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps)$hammerhead_site)

# > pairwise.adonis(Menu_taxmerged_rmhh_Shannons_Hill_Jaccard_dissim, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps)$hammerhead_site)
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 hh_site vs non-hh_site  1 0.9402211 2.512913 0.04222466   0.001      0.001  **

set.seed(200)
hhsite_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix_rm_hh, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.hell.ps)$hammerheads_detected)
write_delim(hhsite_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/hhdetected_landes_pairwise_adonis_results_hell.txt")
#not significant

rda_ps_rm_hh <- ordinate( physeq = Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps,  method = "RDA",distance = "unifrac",formula = ~ Microhabitat+Tide+hammerhead_site+hammerheads_detected+scaled_sst+scaled_depth+scaled_primaryprod+scaled_oxygen+scaled_phosphate + Condition(Site_Name), scale=TRUE)

summary(rda_ps_rm_hh)
anova(rda_ps_rm_hh, by = "margin")
plot(rda_ps_rm_hh)

cca_ps_rm_hh <- ordinate( physeq = Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps,  method = "CCA", formula = ~ Microhabitat+Tide+hammerhead_site+hammerheads_detected+scaled_sst+scaled_depth+scaled_primaryprod+scaled_oxygen+scaled_phosphate + Condition(Site_Name), scale=FALSE, distance = Menu_taxmerged_rmhh_Hill_Jaccard_dissim)

summary(cca_ps_rm_hh)
anova(cca_ps_rm_hh, by = "margin")
plot(cca_ps_rm_hh)

#dbRDA----
#use hilldiv distance matrix on eDNA transformed otu table to do a dbRDA----
#these objects came from hilldiv in script 11!
#Then proceed with dbrda function from vegan
dbRDA_Faith <- dbrda(Menu_taxmerged_Faiths_Hill_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars,na.action = na.omit,sqrt.dist=FALSE) #dbrda from vegan
summary(dbRDA_Faith)
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total           6.130     1.0000
# Conditioned     1.235     0.2015
# Constrained     1.529     0.2494
# Unconstrained   3.366     0.5491

anova(dbRDA_Faith, by ="margin")
# Model: dbrda(formula = Menu_taxmerged_Faiths_Hill_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)   
# Microhabitat           2   0.2514 1.6060  0.145   
# scaled_Water_temp      1   0.1243 1.5874  0.176   
# Tide                   1   0.0841 1.0745  0.344   
# scaled_nitrate         1   0.4198 5.3629  0.002 **
#   hammerheads_detected   1   0.2566 3.2778  0.025 * 
#   carcharhinus_detected  1   0.0640 0.8174  0.488   
# scaled_primaryprod     1   0.3600 4.5983  0.007 **
#   Residual              43   3.3661                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


plot(dbRDA_Faith)

Shark_Menu_taxa_merge_species_named_traits_temp.eDNAindex.ps<- Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps %>% ps_filter(Water_temp !="NA")


speDNAindex<-t(as.data.frame(cbind(otu_table(Shark_Menu_taxa_merge_species_named_traits_temp.eDNAindex.ps))))
sppscores(dbRDA_Faith)<-decostand(speDNAindex,"range")


plot(dbRDA_Faith)


dbRDA_Richness <- dbrda(Menu_taxmerged_Richness_Hill_Sor_dissim ~Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan
summary(dbRDA_Richness) #Richness
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total          13.897     1.0000
# Conditioned     2.140     0.1540
# Constrained     3.318     0.2387
# Unconstrained   8.439     0.6072
anova(dbRDA_Richness, by ="margin")
# Model: dbrda(formula = Menu_taxmerged_Richness_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)   
# Microhabitat           2   0.7810 1.9899  0.002 **
#   scaled_Water_temp      1   0.4045 2.0611  0.026 * 
#   Tide                   1   0.2413 1.2294  0.239   
# scaled_nitrate         1   0.4968 2.5313  0.004 **
#   hammerheads_detected   1   0.2900 1.4778  0.138   
# carcharhinus_detected  1   0.3847 1.9603  0.026 * 
#   scaled_primaryprod     1   0.3436 1.7507  0.061 . 
# Residual              43   8.4389                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
plot(dbRDA_Richness)

dbRDA_Shannon<-dbrda(Menu_taxmerged_Hill_Shannonexp_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(dbRDA_Shannon)
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total          21.769     1.0000
# Conditioned     3.606     0.1656
# Constrained     4.564     0.2097
# Unconstrained  13.599     0.6247
anova(dbRDA_Shannon, by ="margin")  
# Model: dbrda(formula = Menu_taxmerged_Hill_Shannonexp_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)    
# Microhabitat           2   1.2141 1.9196  0.001 ***
#   scaled_Water_temp      1   0.5609 1.7737  0.016 *  
#   Tide                   1   0.5199 1.6441  0.023 *  
#   scaled_nitrate         1   0.6246 1.9750  0.004 ** 
#   hammerheads_detected   1   0.3810 1.2049  0.196    
# carcharhinus_detected  1   0.6539 2.0678  0.004 ** 
#   scaled_primaryprod     1   0.4792 1.5151  0.032 *  
#   Residual              43  13.5987                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
plot(dbRDA_Shannon)

dbRDA_Allen<-dbrda(Allensdistmat_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(dbRDA_Allen)
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total         0.63439    1.00000
# Conditioned   0.04909    0.07738
# Constrained   0.24341    0.38369
# Unconstrained 0.34189    0.53893
anova(dbRDA_Allen, by ="margin")  
# Model: dbrda(formula = Allensdistmat_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)   
# Microhabitat           2  0.07890 4.9615  0.008 **
#   scaled_Water_temp      1  0.03659 4.6023  0.038 * 
#   Tide                   1  0.00998 1.2548  0.293   
# scaled_nitrate         1  0.06937 8.7248  0.003 **
#   hammerheads_detected   1  0.01000 1.2577  0.286   
# carcharhinus_detected  1  0.01613 2.0282  0.188   
# scaled_primaryprod     1  0.03789 4.7657  0.024 * 
#   Residual              43  0.34189                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Model: dbrda(formula = Allensdistmat_Jaccard_dissim ~ Site_Name + Microhabitat + scaled_Water_temp + hammerheads_detected + carcharhinus_detected + Tide + scaled_primaryprod + scaled_oxygen, data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)   
# Site_Name              6  0.04772 1.0006  0.471   
# Microhabitat           2  0.07865 4.9475  0.005 **
#   scaled_Water_temp      1  0.03665 4.6110  0.043 * 
#   hammerheads_detected   1  0.01002 1.2605  0.302   
# carcharhinus_detected  1  0.01614 2.0302  0.156   
# Tide                   1  0.00996 1.2526  0.295   
# scaled_primaryprod     1  0.04825 6.0706  0.012 * 
#   scaled_oxygen          1  0.06946 8.7389  0.004 **
#   Residual              43  0.34180                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
plot(dbRDA_Allen)


dbRDA1<-dbrda(Menu_taxmerged_Hill_Jaccard_dissim ~ Microhabitat + scaled_Water_temp +Tide+scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(dbRDA1)
anova(dbRDA1, by ="margin")       
plot(dbRDA1)

dbRDA2<-dbrda(Menu_taxmerged_Hill_Shannonexp_Jaccard_dissim ~ Microhabitat + scaled_Water_temp +Tide+scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(dbRDA2)
anova(dbRDA2, by ="margin")       
plot(dbRDA2)


dbRDA3<-dbrda(Allensdistmat_Jaccard_dissim ~ Microhabitat + scaled_Water_temp +Tide+scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(dbRDA3)
anova(dbRDA3, by ="margin")       
plot(dbRDA3)

landes_dist_matrix_eDNAindex<-phyloseq::distance(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, method = "l")
dbRDA4<-dbrda(landes_dist_matrix_eDNAindex ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(dbRDA4)
# Partitioning of mean squared Beta.l distance:
#   Inertia Proportion
# Total           90.15     1.0000
# Conditioned     18.77     0.2082
# Constrained     18.83     0.2089
# Unconstrained   52.55     0.5830
anova(dbRDA4, by ="margin")
# Model: dbrda(formula = landes_dist_matrix_eDNAindex ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df Variance      F Pr(>F)  
# Microhabitat           2    4.134 1.6912  0.049 *
#   scaled_Water_temp      1    1.686 1.3799  0.193  
# Tide                   1    1.421 1.1628  0.308  
# scaled_nitrate         1    2.966 2.4265  0.022 *
#   hammerheads_detected   1    2.248 1.8390  0.081 .
# carcharhinus_detected  1    1.988 1.6267  0.099 .
# scaled_primaryprod     1    2.218 1.8152  0.078 .
# Residual              43   52.553                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
plot(dbRDA4)

rda_ps_rm_hh <- ordinate( physeq = Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps,  method = "RDA", formula = ~ Microhabitat+Tide+hammerhead_site+hammerheads_detected+scaled_sst+scaled_depth+scaled_primaryprod+scaled_oxygen+scaled_phosphate + Condition(Site_Name), scale=FALSE, distance = Menu_taxmerged_rmhh_Hill_Jaccard_dissim)
rda_ps_rm_hh
anova(rda_ps_rm_hh, by ="margin")       


rda_ps_rm_hh1 <- ordinate( physeq = Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps,  method = "RDA", formula = ~ Microhabitat+Tide+hammerhead_site+hammerheads_detected+scaled_sst+scaled_depth+scaled_primaryprod+scaled_oxygen+scaled_phosphate + Condition(Site_Name), scale=FALSE,distance = Menu_taxmerged_rmhh_Hill_Jaccard_dissim)
rda_ps_rm_hh1
anova(rda_ps_rm_hh1, by ="margin")       
plot(rda_ps_rm_hh1)

              
HH_rm_dbRDA <- dbrda(Menu_taxmerged_rmhh_Hill_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), HH_rm_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)#,na.action = na.omit) #dbrda from vegan


summary(HH_rm_dbRDA)
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total           6.407     1.0000
# Conditioned     1.198     0.1870
# Constrained     1.618     0.2526
# Unconstrained   3.591     0.5604
anova(HH_rm_dbRDA,by = 
            "margin")




plot(dbRDA, main = "dbRDA")
dbRDA
set.seed(2)
anova.cca(dbRDA, by="margin")

goodness(object = dbRDA,display = "sites")
goodness(object = cca_ps_rm_hh)

upr <- dbRDA#dbrda(Menu_taxmerged_rmhh_Hill_Jaccard_dissim ~ ., data = HH_rm_eDNAindex_envvars,na.action=na.omit)
lwr <- dbrda(Menu_taxmerged_rmhh_Hill_Jaccard_dissim ~ 1+Condition(Site_Name), data = HH_rm_eDNAindex_envvars,na.action=na.omit,sqrt.dist=FALSE)
set.seed(1)
mods <- ordiR2step(lwr, scope = formula(upr), trace = 0)
mods

mod1<-dbrda(formula = Menu_taxmerged_rmhh_Hill_Jaccard_dissim ~ Condition(Site_Name) + Microhabitat+scaled_sst, data = HH_rm_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)

anova(mod1, by ="margin")

#eDNAindex

# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          14.638     1.0000
# Conditioned     1.700     0.1161
# Constrained     2.694     0.1841
# Unconstrained  10.244     0.6998

# Model: cca(formula = OTU ~ Microhabitat + Tide + hammerhead_site + hammerheads_detected + scaled_sst + scaled_depth + scaled_primaryprod + scaled_oxygen + scaled_phosphate + Condition(Site_Name), data = data)
# Df ChiSquare      F Pr(>F)   
# Microhabitat          2    0.5860 1.2300  0.016 * 
#   Tide                  1    0.2791 1.1717  0.132   
# hammerhead_site       0    0.0000    Inf  1.000   
# hammerheads_detected  1    0.2832 1.1889  0.210   
# scaled_sst            1    0.4086 1.7151  0.006 **
#   scaled_depth          1    0.2566 1.0771  0.392   
# scaled_primaryprod    1    0.2630 1.1040  0.344   
# scaled_oxygen         1    0.2623 1.1012  0.321   
# scaled_phosphate      1    0.2625 1.1017  0.306   
# Residual             43   10.2439                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#hellinger

# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          4.9249     1.0000
# Conditioned    0.6316     0.1282
# Constrained    0.9762     0.1982
# Unconstrained  3.3172     0.6736

# Model: cca(formula = OTU ~ Microhabitat + Tide + hammerhead_site + hammerheads_detected + scaled_sst + scaled_depth + scaled_primaryprod + scaled_oxygen + scaled_phosphate + Condition(Site_Name), data = data)
# Df ChiSquare      F Pr(>F)   
# Microhabitat          2    0.1943 1.2595  0.024 * 
#   Tide                  1    0.0802 1.0395  0.260   
# hammerhead_site       0    0.0000   -Inf  1.000   
# hammerheads_detected  1    0.0844 1.0937  0.265   
# scaled_sst            1    0.1672 2.1668  0.006 **
#   scaled_depth          1    0.0594 0.7699  0.625   
# scaled_primaryprod    1    0.0694 0.8999  0.450   
# scaled_oxygen         1    0.0694 0.8997  0.512   
# scaled_phosphate      1    0.0694 0.9001  0.558   
# Residual             43    3.3172                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#When controlling for Site_Name, Microhabitat was significant, hammerheads_detected was significant, and hammerhead_site couldn't be tested, it seems.  

##Check for differences between Carcharhinus sites and non-carcharhinus sites.----
#First remove Carcharhinus reads so they don't confound the calculation.
Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps <- Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps %>% subset_taxa(genus != "Carcharhinus")

landes_dist_matrix_rm_carch<-phyloseq::distance(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps, method = "l")

set.seed(2)
pairwise.adonis(landes_dist_matrix_rm_carch, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps)$carcharhinus_detected)

#eDNAindex
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 no_carch_detected vs carch_detected  1  97.45902 1.103158 0.01898619   0.312      0.312 
#hellinger
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 no_carch_detected vs carch_detected  1  97.45902 1.103158 0.01898619   0.312      0.312    
set.seed(2)
pairwise.adonis(landes_dist_matrix_rm_carch, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps)$blacktip_site)

#eDNAindex
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 carch_site vs non-carch_site  1  207.4864 2.401041 0.04042086   0.019      0.019   .

#hellinger
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 carch_site vs non-carch_site  1  207.4864 2.401041 0.04042086   0.028      0.028   .

set.seed(2)
btsite_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix_rm_carch, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps)$blacktip_site)
write_delim(btsite_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/btsite_landes_pairwise_adonis_results.eDNAindex.txt")

cca_ps_rm_carch <- ordinate( physeq = Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps,  method = "CCA", formula = ~ Microhabitat+Tide+hammerhead_site+hammerheads_detected+scaled_sst+scaled_depth+scaled_primaryprod+scaled_oxygen+scaled_phosphate + Condition(Site_Name), scale=TRUE)

summary(cca_ps_rm_carch)
anova(cca_ps_rm_carch, by = "margin")

#eDNAindex

# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          14.646     1.0000
# Conditioned     1.666     0.1138
# Constrained     2.725     0.1860
# Unconstrained  10.255     0.7002

# # Model: cca(formula = OTU ~ Microhabitat + Tide + hammerhead_site + hammerheads_detected + scaled_sst + scaled_depth + scaled_primaryprod + scaled_oxygen + scaled_phosphate + Condition(Site_Name), data = data)
# Df ChiSquare      F Pr(>F)  
# Microhabitat          2    0.5866 1.2298  0.018 *
#   Tide                  1    0.2817 1.1813  0.111  
# hammerhead_site       0    0.0000    Inf  1.000  
# hammerheads_detected  1    0.3122 1.3091  0.100 .
# scaled_sst            1    0.4116 1.7258  0.011 *
#   scaled_depth          1    0.2561 1.0740  0.360  
# scaled_primaryprod    1    0.2675 1.1217  0.304  
# scaled_oxygen         1    0.2669 1.1190  0.272  
# scaled_phosphate      1    0.2670 1.1195  0.260  
# Residual             43   10.2550                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##########################

#hellinger

#Partitioning of scaled Chi-square:
# Inertia Proportion
# Total          4.8968     1.0000
# Conditioned    0.6283     0.1283
# Constrained    0.9833     0.2008
# Unconstrained  3.2853     0.6709

# Model: cca(formula = OTU ~ Microhabitat + Tide + hammerhead_site + hammerheads_detected + scaled_sst + scaled_depth + scaled_primaryprod + scaled_oxygen + scaled_phosphate + Condition(Site_Name), data = data)
# Df ChiSquare      F Pr(>F)   
# Microhabitat          2    0.1934 1.2660  0.015 * 
#   Tide                  1    0.0796 1.0415  0.248   
# hammerhead_site       0    0.0000   -Inf  1.000   
# hammerheads_detected  1    0.0957 1.2523  0.184   
# scaled_sst            1    0.1673 2.1898  0.009 **
#   scaled_depth          1    0.0591 0.7736  0.624   
# scaled_primaryprod    1    0.0700 0.9160  0.451   
# scaled_oxygen         1    0.0700 0.9159  0.513   
# scaled_phosphate      1    0.0700 0.9163  0.520   
# Residual             43    3.2853                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



Shark_Menu_taxa_merge_species_named_traits.hell.ps@sam_data

#Make regression plot to see if there's a relationship between fewer reads post metabar (just indentifiable target eukaryotes) and water temperature per sample.  Need to check if there are actually fewer species, or just less DNA surviving at higher temperatures

fit = lm(Shark_Menu_ASVs_taxon_named.ps@sam_data$nb_reads ~ Shark_Menu_ASVs_taxon_named.ps@sam_data$Water_temp) # Run a regression analysis
plot(fit)

library(MASS)
fit1<-rlm(Shark_Menu_ASVs_taxon_named.ps@sam_data$nb_reads ~ Shark_Menu_ASVs_taxon_named.ps@sam_data$Water_temp)

plot(fit1)
summary(fit)

Shark_Menu_ASVs_taxon_named.ps@sam_data$nb_motus
readsbytemp <- plot_regression( nb_reads~ Water_temp, meta(Shark_Menu_ASVs_taxon_named.ps),shade = T,spag = F,show.CI =F,show.median = T,shade.alpha = 0.3, show.lm = T)
readsbytemp
summary(readsbytemp)
summary(fit)

#Three lower temp outliers with lots of reads, but they don't cross the Cook's distance lines in plot(fit).
#12, 17, 24 outliers 
Shark_Menu_ASVs_taxon_named.ps@sam_data[24]$Water_temp

# Basic scatter plot.
p1 <- ggplot(Shark_Menu_ASVs_taxon_named.ps@sam_data, aes(x=Water_temp, y=nb_reads)) + 
  geom_point( color="#69b3a2")
p1

# with linear trend
p2 <- ggplot(Shark_Menu_ASVs_taxon_named.ps@sam_data, aes(x=Water_temp, y=nb_reads)) +
  geom_point(color="#69b3a2") +
  geom_smooth(method=lm , color="black", se=TRUE,na.rm = TRUE)
p2

ASV_sample_data<-as.data.frame(cbind(sam_data(Shark_Menu_ASVs_taxon_named.ps)))

readsbywatertemp<-lm(formula = nb_reads~Water_temp,data = ASV_sample_data)
summary(readsbywatertemp)

# linear trend + confidence interval
p3 <- ggplot(Shark_Menu_ASVs_taxon_named.ps@sam_data, aes(x=Water_temp, y=nb_reads)) +
  geom_point(color="royalblue") +
  geom_smooth(method = "loess" , color="royalblue4", fill="#69b3a2", se=TRUE)+
  geom_smooth(method = "lm",color="red",se=FALSE)+
  ylab("Number of ASV Reads per Sample")+
  xlab("Water temperature")
p3#Multiple R-squared:  0.02597,	Adjusted R-squared:  0.008578 
#F-statistic: 1.493 on 1 and 56 DF,  p-value: 0.2268


p4<-ggplot(Shark_Menu_ASVs_taxon_named.ps@sam_data, aes(x=Water_temp, y=nb_motus)) +
  geom_point(color="royalblue") +
  geom_smooth(method = "loess" , color="royalblue4", fill="#69b3a2", se=TRUE)+
  geom_smooth(method = "lm",color="red",se=FALSE)+
  ylab("Number of ASVs per Sample")+
  xlab("Water temperature")
p4
# Multiple R-squared:  0.001356,	Adjusted R-squared:  -0.01648 
# F-statistic: 0.07605 on 1 and 56 DF,  p-value: 0.7837
  
fitotutemp<-lm(Shark_Menu_ASVs_taxon_named.ps@sam_data$nb_motus ~ Shark_Menu_ASVs_taxon_named.ps@sam_data$Water_temp) # Run a regression analysis
#plot(fitotutemp)
summary(fitotutemp)

#make a network graph----
trace(plot_net,edit = TRUE)
# link_layout = function(LinksData, vertexDT) {
#   linkstart = copy(vertexDT[LinksData$v1, 1:2])
#   linkend = copy(vertexDT[LinksData$v2, 1:2])
#   setnames(linkend, old = c("y", "x"), new = c("yend", 
#                                                "xend"))
#   LinksData <- copy(cbind(LinksData, linkstart, linkend))
#   return(LinksData)
# }
## total dataset----




plot_net(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, distance = "l", type ="taxa", maxdist = 0.07,
         laymeth = "circle", color = "class" , shape = NULL,
         rescale = TRUE, point_size = 5, point_alpha = .9, point_label = "species",
         hjust = 1.35, title = NULL)



ig <- make_network(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, max.dist=0.9)
plot_network(ig, Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps,label = "Site_Name", color = "Microhabitat")


#full dataset with bray etc. like the sharks below
menu_network<-plot_net(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, distance = "bray", type ="taxa", maxdist = 0.7,
         laymeth = "auto", color = "troph_type" , shape = NULL,
         rescale = TRUE, point_size = 5, point_alpha = .9, point_label = "species",
         hjust = 1, title = NULL)

menu_network


set.seed(199)
gtmenu<- graph_perm_test(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps,distance = "bray",type = "threshold.value",max.dist = 0.7,sampletype = "scaled_sst")

gtmenu$pval

plot_test_network(gtmenu)
plot_permutations(gtmenu)


set.seed(199)
gtmenu1<- graph_perm_test(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps,distance = "bray",type = "threshold.value",max.dist = 0.7,sampletype = "Microhabitat")

gtmenu1$pval

plot_test_network(gtmenu1)
plot_permutations(gtmenu1)

#calculate some summary statistics
menu_network_data<-as.data.frame(menu_network[["data"]])
#find the number of interactions for each taxon
menusiteinteractingspp<-pivot_longer(menu_network_data,cols = v1:v2,names_to = "v",values_to = "taxa")
#get rid of that " 1" at the end of everything just because it's annoying
interactingtaxanames<-gsub(pattern = " 1",replacement = "",x = menusiteinteractingspp$taxa)
menusiteinteractingspp$taxa<-interactingtaxanames
#count number of interactions per taxon
menusiteinteractingspp$taxa
menuinteractionsperspp<-menusiteinteractingspp%>%
  count(taxa,sort = TRUE,name = "n_menusiteinteractions")%>%
  mutate(menumorethan3=if_else(n_menusiteinteractions>3,
                                   TRUE,FALSE))%>%
  mutate(average_n_edges=mean(n_menusiteinteractions))%>%
  mutate(sd_edges=sd(n_menusiteinteractions))
menuinteractionsperspp

n_menuinteracting_spp<-length(menuinteractionsperspp$taxa)
n_menuinteracting_spp
n_morethanthreeinteractions_menusite<-menuinteractionsperspp%>%
  count(menumorethan3)
n_morethanthreeinteractions_menusite

#percent taxa with more than 4 edges
77/n_menuinteracting_spp


stronginteractions_menusite<-menusiteinteractingspp%>%
  filter(Distance<0.3)

n_strong_menusite_interacting_spp<-length(unique(stronginteractions_menusite$taxa))
n_strong_menusite_interacting_spp
#31
#percent of taxa with strong edges
n_strong_menusite_interacting_spp/n_menuinteracting_spp


#number of edges, strong and weak 0.3 cutoff
menu_network_data
menuedges<-menu_network_data%>%
  mutate(averagemenusitestrength=1-mean(Distance),
         sdmenusitedistance=sd(Distance))
menuedges
totalmenuedges<-length(menu_network[["data"]][["Distance"]])
totalmenuedges
strongmenuedges<-length(subset(menu_network[["data"]][["Distance"]],subset = menu_network[["data"]][["Distance"]]<0.3))
strongmenuedges
percentstrongmenuedges<-(strongmenuedges/totalmenuedges)*100
percentstrongmenuedges

weakmenuedges<-length(subset(menu_network[["data"]][["Distance"]],subset = menu_network[["data"]][["Distance"]]>0.3))
weakmenuedges

percentweakmenuedges<-(weakmenuedges/totalmenuedges)*100
percentweakmenuedges


## hammerhead sites vs. non-hammerhead sites----
#first wisconsin transform the abundances (eDNA index) Kelly et al. 2019 "Understanding PCR Processes to Draw Meaningful Conclusions from Environmental DNA Studies"

Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps<-phyloseq_standardize_otu_abundance(Shark_Menu_taxa_merge_species_named_traits.ps, method = "wisconsin")

#saveRDS(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, "../10_Phyloseq/Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps.RDS")

Hammerhead_sites.ps<-subset_samples(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, Site_Name %in% Hammerhead_sites)
Hammerhead_sites.ps<-prune_taxa(taxa_sums(Hammerhead_sites.ps)>0,Hammerhead_sites.ps)
Hammerhead_sites.ps<-prune_taxa(taxa_sums(Hammerhead_sites.ps)>0,Hammerhead_sites.ps)

Hammerhead_sites.ps@sam_data$Site_Name

nonHammerhead_sites.ps<-subset_samples(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, Site_Name %ni% Hammerhead_sites)
nonHammerhead_sites.ps@sam_data$Site_Name
nonHammerhead_sites.ps<-prune_taxa(taxa_sums(nonHammerhead_sites.ps)>0,nonHammerhead_sites.ps)

Hammerhead_sites.ps@tax_table


hh_network<-plot_net(Hammerhead_sites.ps, distance = 'bray', type ="taxa", maxdist = 0.7,
         laymeth = "circle", color = "troph_type" , shape = NULL,
         point_size = 5, point_alpha = .9, point_label = "species",
         hjust = 1, title = "Network of Taxa in Hammerhead Nursery Bays")
hh_network
hh_network[["data"]]

#calculate some summary statistics
hh_network_data<-as.data.frame(hh_network[["data"]])
#find the number of interactions for each taxon
hhsiteinteractingspp<-pivot_longer(hh_network_data,cols = v1:v2,names_to = "v",values_to = "taxa")
#get rid of that " 1" at the end of everything just because it's annoying
interactingtaxanames<-gsub(pattern = " 1",replacement = "",x = hhsiteinteractingspp$taxa)
hhsiteinteractingspp$taxa<-interactingtaxanames
#count number of interactions per taxon
hhsiteinteractingspp$taxa
hhinteractionsperspp<-hhsiteinteractingspp%>%
  count(taxa,sort = TRUE,name = "n_hhsiteinteractions")%>%
  mutate(hhmorethan3=if_else(n_hhsiteinteractions>3,
                          TRUE,FALSE))%>%
  mutate(average_n_edges=mean(n_hhsiteinteractions))%>%
  mutate(sd_edges=sd(n_hhsiteinteractions))
hhinteractionsperspp

n_hhinteracting_spp<-length(hhinteractionsperspp$taxa)
n_hhinteracting_spp
n_morethanthreeinteractions_hhsite<-hhinteractionsperspp%>%
  count(hhmorethan3)
n_morethanthreeinteractions_hhsite
#A tibble: 2 × 2
# hhmorethan3     n
# <lgl>       <int>
#   1 FALSE          32
# 2 TRUE           43
#percent taxa with more than 4 edges
43/75

stronginteractions_hhsite<-hhsiteinteractingspp%>%
  filter(Distance<0.3)
  
n_strong_hhsite_interacting_spp<-length(unique(stronginteractions_hhsite$taxa))
n_strong_hhsite_interacting_spp
#31
#percent of taxa with strong edges
n_strong_hhsite_interacting_spp/n_hhinteracting_spp


#number of edges, strong and weak 0.3 cutoff
hh_network_data
hhedges<-hh_network_data%>%
  mutate(averagehhsitestrength=1-mean(Distance),
         sdhhsitedistance=sd(Distance))
hhedges
totalhhedges<-length(hh_network[["data"]][["Distance"]])
totalhhedges
stronghhedges<-length(subset(hh_network[["data"]][["Distance"]],subset = hh_network[["data"]][["Distance"]]<0.3))
stronghhedges
percentstronghhedges<-(stronghhedges/totalhhedges)*100
percentstronghhedges

weakhhedges<-length(subset(hh_network[["data"]][["Distance"]],subset = hh_network[["data"]][["Distance"]]>0.3))
weakhhedges

percentweakhhedges<-(weakhhedges/totalhhedges)*100
percentweakhhedges

ggsave("../10_Phyloseq/Network of Taxa in Hammerhead Nursery Bays.jpg")

library("phyloseqGraphTest")

set.seed(199)
gt<- graph_perm_test(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps,sampletype = "hammerhead_site",distance = "bray",type = "threshold.value",max.dist = 0.7)
gt$pval
#0.044 mst seed 199 landes max.dist-0.17
plot_test_network(gt)
plot_permutations(gt)

nonhh_network<-plot_net(nonHammerhead_sites.ps, distance = "bray", type ="taxa", maxdist = 0.7,
         laymeth = "circle", color = "troph_type" , shape = NULL,
         point_size = 5, point_alpha = .9, point_label = "species",
         hjust = 1, title = "Network of Taxa in Shark Nursery Bays with no Hammerheads Detected")
nonhh_network

#calculate some summary statistics
nonhh_network_data<-as.data.frame(nonhh_network[["data"]])
#find the number of interactions for each taxon
nonhhsiteinteractingspp<-pivot_longer(nonhh_network_data,cols = v1:v2,names_to = "v",values_to = "taxa")
#get rid of that " 1" at the end of everything just because it's annoying
interactingtaxanames<-gsub(pattern = " 1",replacement = "",x = nonhhsiteinteractingspp$taxa)
nonhhsiteinteractingspp$taxa<-interactingtaxanames
#count number of interactions per taxon
nonhhsiteinteractingspp$taxa
nonhhinteractionsperspp<-nonhhsiteinteractingspp%>%
  count(taxa,sort = TRUE,name = "n_nonhhsiteinteractions")%>%
  mutate(nonhhmorethan3=if_else(n_nonhhsiteinteractions>3,
                             TRUE,FALSE))%>%
  mutate(average_n_edges=mean(n_nonhhsiteinteractions))%>%
  mutate(sd_edges=sd(n_nonhhsiteinteractions))

nonhhinteractionsperspp
n_nonhhinteracting_spp<-length(nonhhinteractionsperspp$taxa)
n_nonhhinteracting_spp
n_morethanthreeinteractions_nonhhsite<-nonhhinteractionsperspp%>%
  count(nonhhmorethan3)
n_morethanthreeinteractions_nonhhsite
# # A tibble: 2 × 2
# nonhhmorethan3     n
# <lgl>          <int>
#   1 FALSE             23
# 2 TRUE              81
#percent taxa with more than 4 edges
81/n_nonhhinteracting_spp

stronginteractions_nonhhsite<-nonhhsiteinteractingspp%>%
  filter(Distance<0.3)

n_strong_nonhhsite_interacting_spp<-length(unique(stronginteractions_nonhhsite$taxa))
n_strong_nonhhsite_interacting_spp


#percent of taxa with strong edges
n_strong_nonhhsite_interacting_spp/n_nonhhinteracting_spp


#number of edges, strong and weak 0.3 cutoff
nonhh_network_data
nonhhedges<-nonhh_network_data%>%
  mutate(averagenonhhsitestrength=1-mean(Distance),
         sdnonhhsitedistance=sd(Distance))
nonhhedges


totalnonhhedges<-length(nonhh_network[["data"]][["Distance"]])
totalnonhhedges
strongnonhhedges<-length(subset(nonhh_network[["data"]][["Distance"]],subset = nonhh_network[["data"]][["Distance"]]<0.3))
strongnonhhedges
percentstrongnonhhedges<-(strongnonhhedges/totalnonhhedges)*100
percentstrongnonhhedges

weaknonhhedges<-length(subset(nonhh_network[["data"]][["Distance"]],subset = nonhh_network[["data"]][["Distance"]]>0.3))
weaknonhhedges

percentweaknonhhedges<-(weaknonhhedges/totalnonhhedges)*100
percentweaknonhhedges



ggsave("../10_Phyloseq/Network of Taxa in Shark Nursery Bays with no Hammerheads Detected.jpg")




Carcharhinus_sites.ps<-subset_samples(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, Site_Name %in% Blacktip_sites)
Carcharhinus_sites.ps<-prune_taxa(taxa_sums(Carcharhinus_sites.ps)>0,Carcharhinus_sites.ps)
Carcharhinus_sites.ps<-prune_taxa(taxa_sums(Carcharhinus_sites.ps)>0,Carcharhinus_sites.ps)

Carcharhinus_sites.ps@sam_data$Site_Name

nonCarcharhinus_sites.ps<-subset_samples(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, Site_Name %ni% Blacktip_sites)
nonCarcharhinus_sites.ps@sam_data$Site_Name
nonCarcharhinus_sites.ps<-prune_taxa(taxa_sums(nonCarcharhinus_sites.ps)>0,nonCarcharhinus_sites.ps)

Carcharhinus_sites.ps@tax_table

carch_network<-plot_net(Carcharhinus_sites.ps, distance = 'bray', type ="taxa", maxdist = 0.7,
         laymeth = "auto", color = "troph_type" , shape = NULL,
         point_size = 5, point_alpha = .9, point_label = "species",
         hjust = 1, title = "Network of Taxa in Carcharhinus Nursery Bays")

carch_network

#calculate some summary statistics
carch_network_data<-as.data.frame(carch_network[["data"]])
#find the number of interactions for each taxon
carchsiteinteractingspp<-pivot_longer(carch_network_data,cols = v1:v2,names_to = "v",values_to = "taxa")
#get rid of that " 1" at the end of everything just because it's annoying
interactingtaxanames<-gsub(pattern = " 1",replacement = "",x = carchsiteinteractingspp$taxa)
carchsiteinteractingspp$taxa<-interactingtaxanames
#count number of interactions per taxon
carchsiteinteractingspp$taxa
carchinteractionsperspp<-carchsiteinteractingspp%>%
  count(taxa,sort = TRUE,name = "n_carchsiteinteractions")%>%
  mutate(carchmorethan3=if_else(n_carchsiteinteractions>3,
                             TRUE,FALSE))%>%
  mutate(average_n_edges=mean(n_carchsiteinteractions))%>%
  mutate(sd_edges=sd(n_carchsiteinteractions))
carchinteractionsperspp

n_carchinteracting_spp<-length(carchinteractionsperspp$taxa)
n_carchinteracting_spp
n_morethanthreeinteractions_carchsite<-carchinteractionsperspp%>%
  count(carchmorethan3)
n_morethanthreeinteractions_carchsite
# # A tibble: 2 × 2
# carchmorethan3     n
# <lgl>          <int>
#   1 FALSE             29
# 2 TRUE              73
#percent taxa with more than 4 edges
73/n_carchinteracting_spp


stronginteractions_carchsite<-carchsiteinteractingspp%>%
  filter(Distance<0.3)

n_strong_carchsite_interacting_spp<-length(unique(stronginteractions_carchsite$taxa))
n_strong_carchsite_interacting_spp
#31
#percent of taxa with strong edges
n_strong_carchsite_interacting_spp/n_carchinteracting_spp


#number of edges, strong and weak 0.3 cutoff
carch_network_data
carchedges<-carch_network_data%>%
  mutate(averagecarchsitestrength=1-mean(Distance),
         sdcarchsitedistance=sd(Distance))
carchedges
totalcarchedges<-length(carch_network[["data"]][["Distance"]])
totalcarchedges
strongcarchedges<-length(subset(carch_network[["data"]][["Distance"]],subset = carch_network[["data"]][["Distance"]]<0.3))
strongcarchedges
percentstrongcarchedges<-(strongcarchedges/totalcarchedges)*100
percentstrongcarchedges

weakcarchedges<-length(subset(carch_network[["data"]][["Distance"]],subset = carch_network[["data"]][["Distance"]]>0.3))
weakcarchedges

percentweakcarchedges<-(weakcarchedges/totalcarchedges)*100
percentweakcarchedges

ggsave("../10_Phyloseq/Network of Taxa in Carcharhinus Nursery Bays.jpg")


set.seed(1)
gtcarch<- graph_perm_test(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps,sampletype = "blacktip_site",distance = "bray",type = "threshold.value",max.dist = 0.7)
gtcarch$pval

plot_test_network(gtcarch)
plot_permutations(gtcarch)

noncarch_network<-plot_net(nonCarcharhinus_sites.ps, distance = "bray", type ="taxa", maxdist = 0.7,
         laymeth = "auto", color = "troph_type" , shape = NULL,
         point_size = 5, point_alpha = .9, point_label = "species",
         hjust = 1, title = "Network of Taxa in Shark Nursery Bays with no Carcharhinus Detected")
noncarch_network$data
noncarch_network

#calculate some summary statistics
noncarch_network_data<-as.data.frame(noncarch_network[["data"]])
#find the number of interactions for each taxon
noncarchsiteinteractingspp<-pivot_longer(noncarch_network_data,cols = v1:v2,names_to = "v",values_to = "taxa")
#get rid of that " 1" at the end of everything just because it's annoying
interactingtaxanames<-gsub(pattern = " 1",replacement = "",x = noncarchsiteinteractingspp$taxa)
noncarchsiteinteractingspp$taxa<-interactingtaxanames
#count number of interactions per taxon
noncarchsiteinteractingspp$taxa
noncarchinteractionsperspp<-noncarchsiteinteractingspp%>%
  count(taxa,sort = TRUE,name = "n_noncarchsiteinteractions")%>%
  mutate(noncarchmorethan3=if_else(n_noncarchsiteinteractions>3,
                                TRUE,FALSE))%>%
  mutate(average_n_edges=mean(n_noncarchsiteinteractions))%>%
  mutate(sd_edges=sd(n_noncarchsiteinteractions))
noncarchinteractionsperspp

n_noncarchinteracting_spp<-length(noncarchinteractionsperspp$taxa)
n_noncarchinteracting_spp
n_morethanthreeinteractions_noncarchsite<-noncarchinteractionsperspp%>%
  count(noncarchmorethan3)
n_morethanthreeinteractions_noncarchsite

#percent taxa with more than 4 edges
74/n_noncarchinteracting_spp


stronginteractions_noncarchsite<-noncarchsiteinteractingspp%>%
  filter(Distance<0.3)

n_strong_noncarchsite_interacting_spp<-length(unique(stronginteractions_noncarchsite$taxa))
n_strong_noncarchsite_interacting_spp
#31
#percent of taxa with strong edges
n_strong_noncarchsite_interacting_spp/n_noncarchinteracting_spp


#number of edges, strong and weak 0.3 cutoff
noncarch_network_data
noncarchedges<-noncarch_network_data%>%
  mutate(averagenoncarchsitestrength=1-mean(Distance),
         sdnoncarchsitedistance=sd(Distance))
noncarchedges
totalnoncarchedges<-length(noncarch_network[["data"]][["Distance"]])
totalnoncarchedges
strongnoncarchedges<-length(subset(noncarch_network[["data"]][["Distance"]],subset = noncarch_network[["data"]][["Distance"]]<0.3))
strongnoncarchedges
percentstrongnoncarchedges<-(strongnoncarchedges/totalnoncarchedges)*100
percentstrongnoncarchedges

weaknoncarchedges<-length(subset(noncarch_network[["data"]][["Distance"]],subset = noncarch_network[["data"]][["Distance"]]>0.3))
weaknoncarchedges

percentweaknoncarchedges<-(weaknoncarchedges/totalnoncarchedges)*100
percentweaknoncarchedges

ggsave("../10_Phyloseq/Network of Taxa in Shark Nursery Bays with no Carcharhinus Detected.jpg")









#Plot richness----
plot_richness(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, x="Site_Name","Microhabitat", measures=c("Shannon", "Simpson"), color="Site_Name")+
  scale_color_ucscgb()


ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Alpha_Diversity_untransformed_ASVs.jpg"))


# Transform data to proportions as appropriate for Bray-Curtis distances proportion of total reads per sample
Menu.ps.prop <- transform_sample_counts(Menu.ps, function(otu) otu/sum(otu))
plot_bar(Menu.ps.prop, x="Site_Name",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("Read Abundance")



### Percent composition barplot----
library(microbiome)

pseq <- Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps %>%
  aggregate_taxa(level = "order") %>%
  transform(transform = "compositional")

p <- plot_composition(pseq, plot.type = "barplot", group_by = "Microhabitat",verbose = TRUE,sample.sort = "neatmap")+
  ylab("Relative Proportion of Occurrence") + 
  xlab("Sample")+
  scale_fill_discrete(limits = c("Amphipoda", "Decapoda", "Euphausiacea", "Stomatopoda", "Balanomorpha", "Acanthuriformes", "Actinopteri class", "Anguilliformes", "Beloniformes", "Blenniformes", "Carangiformes", "Centrarchiformes", "Chaetodontiformes", "Clupeiformes", "Gerreiformes", "Gonorynchiformes", "Holocentriformes", "Istiophoriformes", "Labriformes", "Lutjaniformes","Mugiliformes","Perciformes","Pleuronectiformes","Scombriformes","Spariformes","Stomiiformes","Syngnathiformes","Tetraodontiformes","Carcharhiniformes")) +
  scale_fill_manual(values = c("Amphipoda"="darkorange", "Decapoda"="salmon", "Euphausiacea"="hotpink", "Scombridae"="darkolivegreen2", "Balanomorpha"="lightpink", "Centrarchiformes"="aquamarine4", "Chaetodontiformes"="green3", "Actinopteri class"="lightskyblue", "Clupeiformes"="lightseagreen", "Gerreiformes"="skyblue3", "Gonorynchiformes"="royalblue", "Anguilliformes"="seagreen3", "Beloniformes"="blue", "Holocentriformes"="slateblue1", "Blenniformes"="purple3", "Carangiformes"="royalblue3", "Centrarchiformes"="slategray2", "Istiophoriformes"="navyblue", "Labriformes"="yellowgreen", "Lutjaniformes"="darkolivegreen1","Mugiliformes"="darkolivegreen","Perciformes"="darkgreen","Pleuronectiformes"="palegreen","Scombriformes"="seagreen2","Spariformes"="aquamarine3","Stomiiformes"="grey","Syngnathiformes"="seagreen","Tetraodontiformes"="slateblue","Carcharhiniformes"="aquamarine")) 

p

p1 <- plot_composition(pseq, plot.type = "barplot", average_by = "Microhabitat",verbose = TRUE,sample.sort = "neatmap")+
  ylab("Relative Proportion of eDNA Transformed Reads") + 
  xlab("Sample")+
  


  scale_fill_discrete(limits = c("Amphipoda", "Decapoda", "Euphausiacea", "Stomatopoda", "Balanomorpha", "Acanthuriformes", "Actinopteri class", "Anguilliformes", "Beloniformes", "Blenniformes", "Carangiformes", "Centrarchiformes", "Chaetodontiformes", "Clupeiformes", "Gerreiformes", "Gonorynchiformes", "Holocentriformes", "Istiophoriformes", "Labriformes", "Lutjaniformes","Mugiliformes","Perciformes","Pleuronectiformes","Scombriformes","Spariformes","Stomiiformes","Syngnathiformes","Tetraodontiformes","Carcharhiniformes")) +
  scale_fill_manual(values = c("Amphipoda"="darkorange", "Decapoda"="salmon", "Euphausiacea"="hotpink", "Stomatopoda"="brown1","Scombridae"="darkolivegreen2", "Balanomorpha"="lightpink","Acanthuriformes"="green", "Centrarchiformes"="aquamarine4", "Chaetodontiformes"="green3", "Actinopteri class"="lightskyblue", "Clupeiformes"="lightseagreen", "Gerreiformes"="skyblue3", "Gonorynchiformes"="royalblue", "Anguilliformes"="seagreen3", "Beloniformes"="blue", "Blenniformes"="green3","Holocentriformes"="slateblue1", "Blenniformes"="purple3", "Carangiformes"="royalblue3", "Centrarchiformes"="slategray2", "Istiophoriformes"="navyblue", "Labriformes"="yellowgreen", "Lutjaniformes"="darkolivegreen1","Mugiliformes"="darkolivegreen","Perciformes"="darkgreen","Pleuronectiformes"="palegreen","Scombriformes"="seagreen2","Spariformes"="aquamarine3","Stomiiformes"="grey","Syngnathiformes"="seagreen","Tetraodontiformes"="slateblue","Carcharhiniformes"="aquamarine"))

p1

p2 <- plot_composition(pseq, plot.type = "barplot", average_by = "Site_Name",sample.sort = "neatmap",
                      otu.sort = "abundance", verbose = TRUE)+
  ylab("Relative Proportion of eDNA Transformed Reads") + 
  xlab("Site")+
  
  # # apply the purple-orange gradient to group 1
  # geom_col(
  #   aes(x = Abundance, group = xlabel, colour = Tax), 
  #   filter(pseq, pseq@tax_table$class!=Actinopteri), 
  #   size = 8, alpha = 0.5) +
  # scale_colour_gradientn(colours = c("purple", "orange")) +
  # labs(colour = "Purple-Orange") +
  # 
  # # start a new scale
  # new_scale_colour() +
  # 
  # # apply the black-grey gradient to group 2
  # geom_jitter(
  #   aes(x, group, colour = x), 
  #   filter(df, group == "group2"), 
  #   size = 8, alpha = 0.5) +
  # scale_colour_gradientn(colours = c("black", "grey80")) +
  # labs(
  #   colour = "Black-Grey",
  #   y = "Group")


 scale_fill_discrete(limits = c("Amphipoda", "Decapoda", "Euphausiacea", "Stomatopoda", "Balanomorpha", "Acanthuriformes", "Actinopteri class", "Anguilliformes", "Beloniformes", "Blenniformes", "Carangiformes", "Centrarchiformes", "Chaetodontiformes", "Clupeiformes", "Gerreiformes", "Gonorynchiformes", "Holocentriformes", "Istiophoriformes", "Labriformes", "Lutjaniformes","Mugiliformes","Perciformes","Pleuronectiformes","Scombriformes","Spariformes","Stomiiformes","Syngnathiformes","Tetraodontiformes","Carcharhiniformes")) +
  scale_fill_manual(values = c("Amphipoda"="darkorange", "Decapoda"="salmon", "Euphausiacea"="hotpink", "Scombridae"="darkolivegreen2", "Balanomorpha"="lightpink", "Centrarchiformes"="aquamarine4", "Chaetodontiformes"="green3", "Actinopteri class"="lightskyblue", "Clupeiformes"="lightseagreen", "Gerreiformes"="skyblue3", "Gonorynchiformes"="royalblue", "Anguilliformes"="seagreen3", "Beloniformes"="blue", "Holocentriformes"="slateblue1", "Blenniformes"="purple3", "Carangiformes"="royalblue3", "Centrarchiformes"="slategray2", "Istiophoriformes"="navyblue", "Labriformes"="yellowgreen", "Lutjaniformes"="darkolivegreen1","Mugiliformes"="darkolivegreen","Perciformes"="darkgreen","Pleuronectiformes"="palegreen","Scombriformes"="seagreen2","Spariformes"="aquamarine3","Stomiiformes"="steelblue","Syngnathiformes"="seagreen","Tetraodontiformes"="dodgerblue","Carcharhiniformes"="aquamarine")) 


p2

table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@tax_table[,"trophic_role"])

table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@tax_table[,"troph_type"])


#add shark site status
Shark_site_samples<-cbind(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps))

Shark_site_samples<-Shark_site_samples%>%
  mutate(Shark_site=if_else(hammerhead_site=="hh_site"&
                              blacktip_site=="carch_site",
                            "S.lewini_and_Carcharhinus",
                            if_else(hammerhead_site=="hh_site"&
                                      blacktip_site=="non-carch_site",
                            "S.lewini",
                            if_else(hammerhead_site=="non-hh_site"&
                                      blacktip_site=="carch_site",
                            "Carcharhinus","neither"))))



sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)<-Shark_site_samples
Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Shark_site

#plot eDNA index by shark site status and troph type----
hhstatusmerge.ps <- merge_samples(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, "hammerhead_site")
plot_bar(hhstatusmerge.ps, fill = "troph_type") + 
  geom_bar(aes(color=troph_type, fill=troph_type), stat="identity", position="stack")

#plot eDNA index by shark site status and troph type----
carchstatusmerge.ps <- merge_samples(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, "blacktip_site")
plot_bar(carchstatusmerge.ps, fill = "troph_type") + 
  geom_bar(aes(color=troph_type, fill=troph_type), stat = "identity", position="stack")

#I think I should do these as percent of occurrence since there were fewer samples without sharks than with, so the uneven sampling might be an issue

#http://gradientdescending.com/how-to-use-multiple-color-scales-in-ggplot-with-ggnewscale/
#library(ggnewscale)

# df |> 
#   ggplot() +
#   p2 <- plot_composition(pseq, plot.type = "barplot", average_by = "Site_Name",sample.sort = "neatmap",
#                          otu.sort = "abundance", verbose = TRUE)+
#   ylab("Relative Percent of Occurrence") + 
#   xlab("Site")+
#   # apply the purple-orange gradient to group 1
#   geom_jitter(
#     aes(x, group, colour = x), 
#     filter(df, group == "group1"), 
#     size = 8, alpha = 0.5) +
#   scale_colour_gradientn(colours = c("purple", "orange")) +
#   labs(colour = "Purple-Orange") +
#   
#   # start a new scale
#   new_scale_colour() +
#   
#   # apply the black-grey gradient to group 2
#   geom_jitter(
#     aes(x, group, colour = x), 
#     filter(df, group == "group2"), 
#     size = 8, alpha = 0.5) +
#   scale_colour_gradientn(colours = c("black", "grey80")) +
#   labs(
#     colour = "Black-Grey",
#     y = "Group")
# #


plot_bar(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, "Microhabitat", fill = "order", facet_grid = ~Site_Name)+
  #geom_bar(aes(color=order, fill=order), stat="identity", position="stack")+
  ylab("eDNA Index Transformed Abundance")+
scale_fill_discrete(limits = c("Amphipoda", "Decapoda", "Euphausiacea", "Stomatopoda", "Balanomorpha", "Acanthuriformes", "Actinopteri class", "Anguilliformes", "Beloniformes", "Blenniformes", "Carangiformes", "Centrarchiformes", "Chaetodontiformes", "Clupeiformes", "Gerreiformes", "Gonorynchiformes", "Holocentriformes", "Istiophoriformes", "Labriformes", "Lutjaniformes","Mugiliformes","Perciformes","Pleuronectiformes","Scombriformes","Spariformes","Stomiiformes","Syngnathiformes","Tetraodontiformes","Carcharhiniformes")) +
  scale_fill_manual(values = c("Amphipoda"="darkorange", "Decapoda"="salmon", "Euphausiacea"="hotpink", "Scombridae"="darkolivegreen2", "Balanomorpha"="lightpink", "Centrarchiformes"="aquamarine4", "Chaetodontiformes"="green3", "Actinopteri class"="lightskyblue", "Clupeiformes"="lightseagreen", "Gerreiformes"="skyblue3", "Gonorynchiformes"="royalblue", "Anguilliformes"="seagreen3", "Beloniformes"="blue", "Holocentriformes"="slateblue1", "Blenniformes"="purple3", "Carangiformes"="royalblue3", "Centrarchiformes"="slategray2", "Istiophoriformes"="navyblue", "Labriformes"="yellowgreen", "Lutjaniformes"="darkolivegreen1","Mugiliformes"="darkolivegreen","Perciformes"="darkgreen","Pleuronectiformes"="palegreen","Scombriformes"="seagreen2","Spariformes"="aquamarine3","Stomiiformes"="steelblue","Syngnathiformes"="seagreen","Tetraodontiformes"="dodgerblue","Carcharhiniformes"="aquamarine")) 


#plot just how many samples fit into each category but facet like above----
sample_table<-as.data.frame(cbind(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))
sample_table<-sample_table%>%
  #select(Site_Name,Microhabitat,hammerhead_site,blacktip_site)%>%
  group_by(Site_Name,Microhabitat)%>%
  count()%>%
  ungroup()

library(ggsci)
ggplot(sample_table, aes(x=Microhabitat, y=n, fill=Microhabitat))+
  geom_bar(stat='identity')+
  facet_grid(~Site_Name)+
  #scale_fill_aaas()
  ylab("Number of Samples")+
  ggtitle("Count of Samples in Each Site and Microhabitat Category")+
  scale_y_continuous(n.breaks = 13)+
  theme(axis.text.x = element_text(angle = 90))
  #theme(axis.text.x = element_blank(), axis.ticks = element_blank())

#save as PDF!

#hammerhead site samples vs non-hammerhead site samples
sample_table<-as.data.frame(cbind(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))
sample_table<-sample_table%>%
  group_by(hammerhead_site,Site_Name)%>%
  mutate(hammerhead_site=if_else(hammerhead_site=="hh_site",
                                 "Hammerhead Site",
                                 "Non-Hammerhead Site"))%>%
  count()%>%
  ungroup()

library(ggsci)
ggplot(sample_table, aes(x=hammerhead_site, y=n, fill=Site_Name))+
  geom_bar(stat='identity', position = "stack")+
  #facet_grid(~Site_Name)+
  scale_fill_viridis_d()+
  ylab("Number of Samples")+
  xlab("Hammerhead Site Status")+
  ggtitle("Count of Samples in Sites where Scalloped Hammerhead DNA was Detected or Not Detected")+
  scale_y_continuous(n.breaks = 20)

#theme(axis.text.x = element_blank(), axis.ticks = element_blank())

#save as PDF

#carcharhinus site samples vs non-carcharhinus site samples
sample_table<-as.data.frame(cbind(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))
sample_table<-sample_table%>%
  group_by(blacktip_site,Site_Name)%>%
  mutate(blacktip_site=if_else(blacktip_site=="carch_site",
                                 "Carcharhinus Site",
                                 "Non-Carcharhinus Site"))%>%
  count()%>%
  ungroup()

library(ggsci)
ggplot(sample_table, aes(x=blacktip_site, y=n, fill=Site_Name))+
  geom_bar(stat='identity', position = "stack")+
  #facet_grid(~Site_Name)+
  scale_fill_viridis_d()+
  ylab("Number of Samples")+
  xlab("Carcharhinus Site Status")+
  ggtitle("Count of Samples in Sites where Carcharhinus DNA was Detected or Not Detected")+
  scale_y_continuous(n.breaks = 20)
#theme(axis.text.x = element_blank(), axis.ticks = element_blank())

#save as PDF and then combine the two 

#transform to percent of sample
#Menu.r.ps = transform_sample_counts(Shark_Menu_taxa_merge_species_named_traits.ps, function(x) x/sum(x))

#select only otus with variance of their transformed abundance greater than 0.00001
#Menu.p.ps = filter_taxa(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, function(x) var(x) > 1e-03, TRUE)

#plot just the ones that vary more than the cutoff in the line above.  
# plot_bar(Menu.r.ps, "Microhabitat", fill = "order",facet_grid = ~Site_Name)+
#   ylab("eDNA Index Transformed Abundance")


##View the Top 20 most abundant taxa in the dataset----

#take just the top 20 taxa
top25 = names(sort(taxa_sums(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps), decreasing = TRUE)[1:25])
top25_taxa.ps = prune_taxa(top25, Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)
plot_bar(top25_taxa.ps, x="Site_Name",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")+
  ylab("eDNA Index Transformed Read Abundance")

top25taxa<-as.data.frame(tax_table(top25_taxa.ps))%>%
  select(species,order,class)
top25taxa


p3<-plot_bar(top25_taxa.ps, x="Microhabitat",fill="species",facet_grid = ~Site_Name)+ 
  geom_bar(aes(color=NULL, fill=species), stat="identity", position="stack")+
  ylab("eDNA Index Transformed Read Abundance")+
  xlab("Site Name")+
  ggtitle("Top 25 Most Detected Taxa by Site and Microhabitat")+
scale_fill_discrete(limits = c("Pachygrapsus socius","Hippa pacifica","Crinotonia attenuatus","Mugil cephalus", "Mugil thoburni", "Mugil genus","Carcharhinus genus",  "Gerres cinereus", "Halichoeres nicholsi", "Anisotremus scapularis", "Bodianus diplotaenia", "Archosargus genus", "Sphoeroides annulatus","Diodon genus", "Scombridae family","Thunnus genus","Pontinus genus", "Stegastes beebei", "Stegastes arcifrons", "Caranx genus", "Abudefduf genus", "Tylosurus crocodilus", "Hyporhamphus genus", "Epinephelus genus", "Eucinostomus genus")) +
  scale_fill_manual(values = c("Pachygrapsus socius"="salmon1", "Hippa pacifica"="tomato","Crinotonia attenuatus"="red","Mugil cephalus"="lightskyblue", "Mugil thoburni"="cornflowerblue", "Mugil genus"="lightskyblue1","Carcharhinus genus"="greenyellow",  "Gerres cinereus"="skyblue3", "Halichoeres nicholsi"="aquamarine", "Anisotremus scapularis"="springgreen3", "Bodianus diplotaenia"="turquoise", "Archosargus genus"="dodgerblue3", "Sphoeroides annulatus"="mediumpurple","Diodon genus"="lightslateblue", "Scombridae family"="cadetblue1", "Thunnus genus"="cadetblue3","Pontinus genus"="darkseagreen", "Stegastes beebei"="blue", "Stegastes arcifrons"="dodgerblue", "Caranx genus"="darkslategray4", "Abudefduf genus"="darkseagreen1", "Tylosurus crocodilus"="deepskyblue", "Hyporhamphus genus"="deepskyblue3", "Epinephelus genus"="cyan", "Eucinostomus genus"="chartreuse4"))

p3

#count how many samples of each category there are:

sample_table<-sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)

#count how many taxa are in each category

subset.ps <- ps_filter(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, Microhabitat =="middle of bay")
subset.ps <- prune_samples(sample_sums(subset.ps)>0, subset.ps)
ntaxa(subset.ps)


#ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Top50_ASVs_by_site_and_genus.jpg"))

# top10 = names(sort(taxa_sums(Menu.ps.prop), decreasing = TRUE)[1:10])
# top10_Menu.ps.prop = prune_taxa(top10, Menu.ps.prop)
# plot_bar(Menu.ps.prop, x="Site_Name",fill="genus")+ 
#   geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
#   ylab("Proportional Read Abundance")
# 
# ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Top50_ASVs_by_site_and_genus_proportions.jpg"))
# 

#Merge_species but keep taxa that aren't ID'ed to species----
Menu_taxa_merge.ps <- tax_glom(Menu.ps, taxrank = "species", NArm = FALSE)
#Merge_species and remove taxa that aren't ID'ed to species----
Menu_species_merge1.ps<-tax_glom(Menu.ps, taxrank = "species", NArm=TRUE)


#convert to presence absence----

Menu_species_merge1.pa.ps<-phyloseq_standardize_otu_abundance(Menu_species_merge1.ps, method = "pa")

#plot frequency of occurrence(FOO) by genus and site
plot_bar(Menu_species_merge1.pa.ps, x="Site_Name",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")+
  ylab("Frequency of Occurrence")

#ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Frequency_of_species_occurrence_by_Site.jpg"))

##find the top 20 most frequently occurring species----

top20occ = names(sort(taxa_sums(Menu_species_merge1.pa.ps), decreasing = TRUE)[1:20])
top20occ_Menu_species_merge1.pa.ps = prune_taxa(top20occ, Menu_species_merge1.pa.ps)

plot_bar(top20occ_Menu_species_merge1.pa.ps, x="Site_Name",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")+
  ylab("Frequency of Occurrence")

# #Merge taxa by trophic_role
# Menu_trophic_merge.ps <- merge_taxa(Shark_Menu_taxa_merge_species_named_traits.ps,eqtaxa = "trophic_role")
# 
# Menu_trophic_merge_eDNAindex.ps<-phyloseq_standardize_otu_abundance( Menu_trophic_merge.ps, method = "wisconsin")
# 
# 
# plot_bar(Menu_trophic_merge_eDNAindex.ps, fill = "trophic_role") + 
#   geom_bar(aes(color=trophic_role, fill=trophic_role), stat="identity", position="stack")+
#   #facet_wrap(facets = c("hammerhead_site","blacktip_site"))+
#   ylab("eDNA Index of taxa merged by Trophic Role")




# 
#plot percent of occurrence (POO) by genus and site----
#so that ASVs found in more samples in a given site have a higher POO value than ones found in fewer samples.
#the 100% value will be number of positive detections for each site, and the values for each taxon in the list will be the rowsum of the presence-absence transformed otu table grouped by site.
#remove sites with less than (cutoff) samples?



#Plot trophic level percent of occurrence (relative to occurrences per hammerhead or carcharhinus "treatments")----


##(Menu_species_merge1.pa.ps is just the ASVs identified to species level)

###Use Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps
###Use Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps for hammerheads

Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps

#Take the phyloseq tables into dataframes and visualize
plot.new()
##Hammerhead site trophic percentages----
otu_df<-cbind(data.frame(otu_table( Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps)))
otu_df$otu<-rownames(otu_df)

otu_df_long<-otu_df%>%
  pivot_longer(cols = -otu,names_to = "sample",values_to = "presence")
otu_df_long$sample <-gsub(pattern = ".",
                          "-",
                          x=otu_df_long$sample,
                          fixed = TRUE)



sample_df<-cbind(data.frame(sample_data( Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps)))
sample_df<-rownames_to_column(sample_df, var = "sample")
sample_site_key<-sample_df%>%
  dplyr::select(sample,hammerhead_site)

sample_site_key$hammerhead_site <-gsub(pattern = "hh_site",
                            "Hammerhead Site",
                            "non-hh_site","Non-Hammerhead Site",
                            x=sample_site_key$hammerhead_site,
                            fixed = TRUE)

taxa_df<-cbind(data.frame(tax_table( Shark_Menu_taxa_merge_species_named_traits_rm_hh.ps)))
taxa_df<-rownames_to_column(taxa_df,var="otu")
otu_species_key<-taxa_df%>%
  dplyr::select(otu,trophic_role)

otu_df_long<-otu_df_long%>%
  left_join(sample_site_key)%>%
  left_join(otu_species_key)

merge_site_otu_counts<-otu_df_long%>%
  group_by(hammerhead_site,trophic_role)%>%
  summarise(total_count=sum(presence))%>%
  ungroup()

merge_site_counts<-merge_site_otu_counts%>%
  group_by(hammerhead_site)%>%
  mutate(total_site=sum(total_count))%>%
  ungroup()%>%
  mutate(rel_prct=(total_count/total_site)*100)%>%
  mutate(trophic_role=factor(trophic_role,
                             levels=c("Predator","Mesopredator","Prey")))

merge_site_counts<-merge_site_counts%>%
  filter(trophic_role!="NA")


#plot 
ggplot(merge_site_counts,fill="trophic_role")+
  geom_bar(aes(x=hammerhead_site, y=rel_prct, fill=trophic_role), stat="identity", position="dodge")+
  title("Relative Percent of Occurrence of trophic roles by hammerhead nursery status")+
  ylab("Relative Percent of Occurrence")+
  xlab("Hammerhead Site Status")+
  scale_fill_discrete(limits =c("Predator","Mesopredator","Prey"))+
  scale_fill_manual(values = c("Predator"="red2","Mesopredator"="darkorange","Prey"="gold"))
  


ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Relative_percent_of_trophic_role_occurrence_by_hhstatus.jpg"))

##Carcharhinus status trophic percentages----

###Use Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps for Carcharhinus sites

sample_df<-cbind(data.frame(sample_data( Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps)))
sample_df<-rownames_to_column(sample_df, var = "sample")
sample_site_key<-sample_df%>%
  dplyr::select(sample,blacktip_site)

sample_site_key$blacktip_site <-gsub(pattern = "carch_site",
                                       "Carcharhinus Site",
                                       "non-carch_site","Non-Carcharhinus Site",
                                       x=sample_site_key$blacktip_site,
                                       fixed = TRUE)

taxa_df<-cbind(data.frame(tax_table( Shark_Menu_taxa_merge_species_named_traits_rm_carch.ps)))
taxa_df<-rownames_to_column(taxa_df,var="otu")
otu_species_key<-taxa_df%>%
  dplyr::select(otu,trophic_role)

otu_df_long<-otu_df_long%>%
  left_join(sample_site_key)%>%
  left_join(otu_species_key)

merge_site_otu_counts<-otu_df_long%>%
  group_by(blacktip_site,trophic_role)%>%
  summarise(total_count=sum(presence))%>%
  ungroup()

merge_site_counts<-merge_site_otu_counts%>%
  group_by(blacktip_site)%>%
  mutate(total_site=sum(total_count))%>%
  ungroup()%>%
  mutate(rel_prct=(total_count/total_site)*100)%>%
  mutate(trophic_role=factor(trophic_role,
                             levels=c("Predator","Mesopredator","Prey")))
merge_site_counts<-merge_site_counts%>%
  filter(trophic_role!="NA")

#plot 
ggplot(merge_site_counts,fill="trophic_role")+
  geom_bar(aes(x=blacktip_site, y=rel_prct, fill=trophic_role), stat="identity", position="dodge")+
  title("Relative Percent of Occurrence of trophic roles by Carcharhinus status")+
  ylab("Relative Percent of Occurrence")+
  xlab("Carcharhinus Site Status")+
  scale_fill_discrete(limits =c("Predator","Mesopredator","Prey"))+
  scale_fill_manual(values = c("Predator"="red2","Mesopredator"="darkorange","Prey"="gold"))



ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Relative_percent_of_trophic_role_occurrence_by_carchstatus.jpg"))


#multiple plots for various levels of taxonomy?

#two color schemes in one plot (fish v crustaceans)

##make a column in the data for Fish or crustaceans, make fill and color that.
#make shades of orange and shades of blue

#use different color palettes in the same plot for different categories.

#beeswarm plot?



######done plotting relative percent of occurrence by species (no NAs)

#Plot relative percent of occurrence by site for families----

#Menu_taxa_merge.ps is merging all taxa including those not ID'ed to the species level
#convert to presence absence:
Menu_taxa_merge.pa.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge.ps, method = "pa")

#Put it into dataframes
otu_df<-cbind(data.frame(otu_table(Menu_taxa_merge.pa.ps)))
otu_df$otu<-rownames(otu_df)

otu_df_long<-otu_df%>%
  pivot_longer(cols = -otu,names_to = "sample",values_to = "presence")
otu_df_long$sample <-gsub(pattern = ".",
                          "-",
                          x=otu_df_long$sample,
                          fixed = TRUE)

sample_df<-cbind(data.frame(sample_data(Menu_taxa_merge.pa.ps)))
sample_df<-rownames_to_column(sample_df, var = "sample")
sample_site_key<-sample_df%>%
  dplyr::select(sample,Site_Name)

taxa_df<-cbind(data.frame(tax_table(Menu_taxa_merge.pa.ps)))
taxa_df<-rownames_to_column(taxa_df,var="otu")
otu_species_key<-taxa_df%>%
  dplyr::select(otu,family,genus,species)

otu_df_long<-otu_df_long%>%
  left_join(sample_site_key)%>%
  left_join(otu_species_key)

#group by different taxonomic levels (family)

merge_site_otu_counts<-otu_df_long%>%
  group_by(Site_Name,family)%>%
  summarise(total_count=sum(presence))%>%
  ungroup()

merge_site_counts<-merge_site_otu_counts%>%
  group_by(Site_Name)%>%
  mutate(total_site=sum(total_count))%>%
  ungroup()%>%
  mutate(rel_prct=total_count/total_site)


#plot 
ggplot(merge_site_counts,fill="family")+
  geom_bar(aes(x=Site_Name, y=rel_prct, color=family, fill=family), color="black",stat="identity", position="stack")+
  ylab("Relative Percent of Occurrence by Site")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Relative_percent_of_family_occurrence_by_Site.jpg"))






#calculate otu occurrence (POO?) with the metagMisc package----

POO<-phyloseq_otu_occurrence(
  Menu_taxa_merge_species_named.ps,
  variable = "Site_Name",
  taxa_frequency = "percentage",
  drop_zeroes = TRUE,
  justdf = FALSE,
  long = FALSE)



plot_bar(POO,fill="species")+
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")+
  ylab("Percent of Occurrence")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Percent_of_species_occurrence_by_Site_and_genus.jpg"))

#by family
plot_bar(POO,fill="family")+
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")+
  ylab("Percent of Occurrence")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Percent_of_species_occurrence_by_Site_and_family.jpg"))

#by order
plot_bar(POO,fill="order")+
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack")+
  ylab("Percent of Occurrence")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Percent_of_species_occurrence_by_Site_and_order.jpg"))



##20 rarest ASVs ----
lowest20 = names(sort(taxa_sums(Menu.ps), decreasing = FALSE)[1:20])
lowest20_Menu.ps.prop = prune_taxa(lowest20, Menu.ps)
plot_bar(lowest20_Menu.ps.prop, x="Site_Name",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat = "identity", position="stack")+
  ylab("Proportional Read Abundance")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Lowest20_ASVs_by_site_and_genus.jpg"))


###############

library(RColorBrewer)
#top20 ASVs in the dataset
top20 = names(sort(taxa_sums(Menu.ps), decreasing = TRUE)[1:20])
top20_Menu.ps.prop = prune_taxa(top20, Menu.ps.prop)

plot_bar(top20_Menu.ps.prop, x="Microhabitat",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("Proportional Read Abundance")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Top20_ASVs_by_site_and_genus.jpg"))

#Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(Menu.ps))
standf = function(x, t=total) round(t * (x / sum(x)))
Menu.ps.med.norm = transform_sample_counts(Menu.ps, standf)

plot_bar(Menu.ps.med.norm, x="Microhabitat",fill="family")+ 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")+
  ylab("Read Abundance")
ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"median_normalized_read_abundance_combined_microhabitat.jpg"))



plot_bar(Menu.ps.med.norm, fill = "genus", facet_grid = "Site_Name")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Median_normalized_genus_faceted_by_site.jpg"))


plot_bar(Menu.ps.med.norm, fill = "genus", facet_grid = "Microhabitat")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Median_normalized_genus_faceted_by_microhabitat.jpg"))


plot_bar(Menu.ps.prop, fill = "genus", facet_grid = "Site_Name")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("proportional read abundance")
ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_proportions_genus_faceted_by_site.jpg"))



plot_bar(Menu.ps.prop, fill = "family", facet_grid = "Site_Name")+ 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")+
  ylab("proportional read abundance")




#different organization of the graph
plot_bar(Menu.ps.prop, x="order", fill = "order", facet_grid = Microhabitat~Site_Name) +
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack")


#just the top 20
plot_bar(top20_Menu.ps.prop, x="order", fill = "order", facet_grid = Microhabitat~Site_Name) +
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack")


Menu.ps.ord <- ordinate(Menu.ps.prop, "NMDS", "bray")
plot_ordination(Menu.ps,Menu.ps.ord , type="Site_Name", color="Site_Name", 
                shape="Microhabitat", title="biplot", label = "Site_Name") +  
  geom_point(size=3)


#This is the most important part to me to make the otu table readable.
#merge motus by prey species
Menu_taxa_merge.ps <- tax_glom(Menu.ps, taxrank = "species", NArm = FALSE)



Menu_taxa_merge.ps.prop <- transform_sample_counts(Menu_taxa_merge.ps, function(otu) otu/sum(otu))

plot_bar(Menu_taxa_merge.ps.prop, fill = "family") + 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")


# Inspect and Save the taxa-merged phyloseq object
Menu_taxa_merge.ps
sample_variables(Menu_taxa_merge.ps)
rank_names(Menu_taxa_merge.ps)
tax_table(Menu_taxa_merge.ps)

saveRDS(Menu_taxa_merge.ps, file = paste0("../10_Phyloseq/",Primer,"/",Primer,"_taxa_merged_phyloseq.rds"))

#Heatmaps----
plot_heatmap(Menu_taxa_merge.ps.prop, method = "NMDS", distance = "bray", taxa.label = "genus" )

sorted_Menu <- prune_taxa(names(sort(taxa_sums(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps ),TRUE)[1:171]), Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)
samples<-as.data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps))
samples<-samples%>%
  mutate(site_microhabitat=paste(Site_Name, Microhabitat))
sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)<-samples

plot_heatmap(sorted_Menu, taxa.label= "species", distance="l")

# Plot richness----
plot_richness(Menu_taxa_merge.ps, x="Water_temp", color="Microhabitat", measures=c("Simpson"))+
  scale_color_locuszoom()+
  scale_shape_girafe_filled() +
  ggplot2::stat_smooth(
    ggplot2::aes(colour = Microhabitat))


ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_alpha_diversity_by water_temp_and_Microhabitat.jpg"))

plot_richness(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps , x="Site_Name", color="Site_Name", measures = c("Observed"))+ 
  geom_boxplot()+
  scale_color_aaas()

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_alpha_diversity_by_site.jpg"))

plot_richness(Shark_Menu_taxa_merge_species_named_traits.ps, x="Microhabitat", color="Microhabitat", measures = c("Observed","invSimpson"))+ 
  geom_boxplot()+
  scale_color_locuszoom()

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_alpha_diversity_by_microhabitat.jpg"))



#Ordinations----#Ordinations----#Ordinations----

#Following the microviz tutorial https://david-barnett.github.io/microViz/index.html
#awesome stuff!

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)

#install.packages(
#  "microViz",
#  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#)

#install.packages("ggraph") # for taxatree_plots()
#install.packages("DT") # for tax_fix_interactive()
#install.packages("corncob") # for example datasets and beta binomial models
library(microViz)
library(ggraph)
library(DT)
library(corncob)


#fix taxa in Menu_taxa_merge.ps
#tax_fix_interactive(Menu_taxa_merge.ps)

pseq<-Menu_taxa_merge.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )


# # play with ordinations interactively----


pseq
pseq <- Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps %>%
 phyloseq_validate()
#
#tax_fix_interactive(pseq)
microViz::ord_explore(pseq)

#code from the above program...
# 
# 
# #PCoA compositional transformation dist=bray with ellipses
# pseq %>%
#   tax_transform(rank = "genus", trans = "compositional") %>%
#   dist_calc(dist = "bray") %>%
#   ord_calc(
#     method = "auto"
#   ) %>% 
#   ord_plot(
#     axes = c(1, 2),
#     colour = "Site_Name", fill = "Site_Name",
#     shape = "Microhabitat", 
#     alpha = 0.5,
#     size = 2
#   ) + 
#   scale_shape_girafe_filled() +
#   ggplot2::stat_ellipse(
#     ggplot2::aes(colour = Site_Name)
#   )+
#   scale_color_flatui()
# 
# ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_compositional PcoA with ellipses.jpg"))
# 
# 
# 
# #PCA of compositional proportion of prey with arrows
# 
# pseq %>%
#   tax_transform(rank = "species", trans = "compositional") %>%
#   ord_calc(
#     method = "PCA"
#   ) %>% 
#   ord_plot(
#     axes = c(1, 2),
#     plot_taxa = 1:5,
#     colour = "Microhabitat", fill = "Microhabitat",
#     shape = "circle", alpha = 0.7,
#     size = 2
#   )+
#   scale_color_locuszoom()
# 
# ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_PCA_with_arrows_by_microhabitat.jpg"))
# 
# #PCA of binary transformed (presence/absence) prey species.
# pseq %>%
#   tax_transform(rank = "species", trans = "binary") %>%
#   ord_calc(
#     method = "auto"
#   ) %>% 
#   ord_plot(
#     axes = c(1, 2),
#     colour = "Site_Name", fill = "Site_Name",
#     shape = "circle", alpha = 0.5,
#     size = 2
#   ) +
#   ggplot2::stat_ellipse(
#     ggplot2::aes(colour = Site_Name)
#   )+
#   scale_color_flatui()
# 
# ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_presence_absence_PCA_with_ellipses.jpg"))
# 
# 
# 
# print(paste0("Phyloseq exploratory visualizations for this run can be found in ../10_Phyloseq/",Primer,"/"))
# 





