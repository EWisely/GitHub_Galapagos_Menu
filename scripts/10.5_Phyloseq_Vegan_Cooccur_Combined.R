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

new_sample_data<-cbind(data.frame(sample_data(Menu1.ps)))

new_sample_data1<-new_sample_data%>%
  mutate(hammerheads_detected=
           if_else(sample_id %in% Hammerhead_samples, 
                   1,
                   0))%>%
  mutate(carcharhinus_detected=
           if_else(sample_id %in% Blacktip_samples, 
                   1,
                   0))%>%
  mutate(blacktip_site=
           if_else(Site_Name %in% Blacktip_sites, 
                   1,
                   0))%>%
  mutate(hammerhead_site=
           if_else(Site_Name %in% Hammerhead_sites,
                   1,
                   0))


new_sample_data1
Shark_Menu.ps<-Menu1.ps

sample_data(Shark_Menu.ps)<-new_sample_data1

Shark_Menu.ps@sam_data$carcharhinus_detected


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
  select("species")
#50 ASVs were not identifiable to species level. (because there are 121 species plus other ASVs)
#71 unique species found!

sealife_table<-fb_tbl("species","sealifebase") %>% 
  mutate(sci_name = paste(Genus, Species)) %>%
  filter(sci_name %in% species_names$species) %>% 
  select(sci_name, SpecCode, FBname,DemersPelag,DepthRangeShallow, DepthRangeDeep, Vulnerability, Length)

fish_table<-fb_tbl("species","fishbase")%>% 
  mutate(sci_name = paste(Genus, Species)) %>%
  filter(sci_name %in% species_names$species) %>% 
  select(sci_name, SpecCode, FBname,DemersPelag,DepthRangeShallow, DepthRangeDeep, Vulnerability, Length)

Menu_species_table<-rbind(fish_table, sealife_table)

crustacean_ecology_traits<-ecology(species_list = species_names$species, server = "sealifebase")
crustacean_ecology_traits


fish_ecology_traits<-ecology(species_list = species_names$species, server = "fishbase")
fish_ecology_traits

#get the FoodTroph, FoodSeTroph, Neritic, Intertidal, Oceanic, Epipelagic, Neritic, Mesopelagic, Estuaries, Mangroves, FeedingType, Solitary, and AddRems columns (and add them to the Menu_species_table)
fish_ecology_traits1<-fish_ecology_traits%>%
  select(Species,FoodTroph, FoodSeTroph, Neritic, Intertidal, Oceanic, Epipelagic, Neritic, Mesopelagic, Estuaries, Mangroves, FeedingType, Solitary, AddRems)

crustacean_ecology_traits<-crustacean_ecology_traits%>%
  dplyr::rename(Mesopelagic = mesopelagic)

crustacean_ecology_traits1<-crustacean_ecology_traits%>%
  select(Species, FoodTroph, FoodSeTroph, Neritic, Intertidal, Oceanic, Epipelagic, Neritic, Mesopelagic, Estuaries, Mangroves, FeedingType, Solitary, AddRems)

Menu_species_table<-Menu_species_table%>%
  dplyr::rename(Species = sci_name)


##Make a traits table for everything that could be ID'ed to species.----
Menu_traits_table<-left_join(Menu_species_table, fish_ecology_traits1, by = "Species")

Menu_traits_table<-left_join(Menu_traits_table, crustacean_ecology_traits1)


#Add the traits to the taxonomy table----

taxonomy_table<-cbind(data.frame(tax_table(Shark_Menu_taxa_merge_species_named.hell.ps)))

taxonomy_table<-rownames_to_column(taxonomy_table, var="ID")

Menu_traits_table<-Menu_traits_table%>%
  dplyr::rename(species = Species)

taxo_traits_table<-left_join(taxonomy_table,Menu_traits_table, by ="species")

taxo_traits_table<-column_to_rownames(taxo_traits_table, var = "ID")

taxo_traits_matrix<-as.matrix(taxo_traits_table)

##put traits back into the phyloseq object----

Shark_Menu_taxa_merge_species_named_traits.hell.ps<-Shark_Menu_taxa_merge_species_named.hell.ps

tax_table(Shark_Menu_taxa_merge_species_named_traits.hell.ps)<-taxo_traits_matrix

Shark_Menu_taxa_merge_species_named_traits.hell.ps@tax_table

saveRDS(Shark_Menu_taxa_merge_species_named_traits.hell.ps, file = paste0("../10_Phyloseq/",Primer,"/",Primer,"_hellinger_transformed_traits_added_taxa_named_Phyloseq.rds"))

###Also for the ASV phyloseq object----
ASV_taxonomy_table<-cbind(data.frame(tax_table(Shark_Menu_ASVs_taxon_named.hell.ps)))

ASV_taxonomy_table<-rownames_to_column(ASV_taxonomy_table, var="ID")

Menu_traits_table<-Menu_traits_table%>%
  dplyr::rename(species = Species)

ASV_traits_table<-left_join(ASV_taxonomy_table,Menu_traits_table, by ="species")

ASV_traits_table<-column_to_rownames(ASV_traits_table, var = "ID")

ASV_traits_matrix<-as.matrix(ASV_traits_table)

Shark_Menu_ASVs_taxon_named_traits.hell.ps<-Shark_Menu_ASVs_taxon_named.hell.ps

tax_table(Shark_Menu_ASVs_taxon_named_traits.hell.ps)<-ASV_traits_matrix

Shark_Menu_ASVs_taxon_named_traits.hell.ps@tax_table


Shark_Menu_ASVs_taxon_named_traits.hell.ps@sam_data$Microhabitat

Shark_Menu_ASVs_taxon_named_traits.hell.ps@sam_data$Tide


##Make a phyloseq object of taxa that has DemersPelag info----
Shark_Menu_taxa_merge_species_named_DemersPelag.hell.ps<- subset_taxa(Shark_Menu_taxa_merge_species_named_traits.hell.ps, DemersPelag != "NA")

Shark_Menu_taxa_merge_species_named_DemersPelag.hell.ps@tax_table

##Make a phyloseq object of taxa that has no missing traits----
Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps<- subset_taxa(Shark_Menu_taxa_merge_species_named_DemersPelag.hell.ps, FoodTroph != "NA")

Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps@tax_table

#has all traits complete except DepthRangeShallow and DepthRangeDeep, I'll just delete those columns for the gllvm

##Remove samples that are missing temperature data----

Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps <- Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps %>% ps_filter(Water_temp !="NA")

Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps

#try a gllvm using the traits phyloseq----
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

#the traits were only available for species. The list was built off of this phyloseq object: Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps

menu_gllvm_data$trait
menu_gllvm_trait<-cbind(data.frame(tax_table(Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps)))

menu_gllvm_trait<-rownames_to_column(menu_gllvm_trait, var="ID")
colnames(menu_gllvm_trait)
menu_gllvm_trait<-menu_gllvm_trait%>%
  select(ID, FoodTroph, Epipelagic, Neritic, Mesopelagic, Solitary, DemersPelag, Intertidal, Estuaries, Oceanic, Mangroves)
menu_gllvm_trait

menu_gllvm_data$trait<-menu_gllvm_trait

menu_gllvm_data$abund
ALLtraits_otu_matrix<-as(otu_table(Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps), "matrix")
ALLtraits_otu_table<-as.data.frame(ALLtraits_otu_matrix)
menu_gllvm_data$abund<-t(ALLtraits_otu_table)
menu_gllvm_data$abund


menu_gllvm_data$x
menu_gllvm_env_vars<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_ALL_traits.hell.ps)))
menu_gllvm_env_vars<-menu_gllvm_env_vars%>%
  select(Site_Name,Microhabitat,Tide, Water_temp)
menu_gllvm_env_vars

menu_gllvm_data$x<-menu_gllvm_env_vars

fitp <- gllvm(y = menu_gllvm_data$abund, family = poisson(), num.lv = 2)
fitnb <- gllvm(y = menu_gllvm_data$abund, family = "negative.binomial", num.lv = 2)
AIC(fitp)
## [1] 704.4989
AIC(fitnb)
## [1] 772.4951

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

plot.new()
ordiplot(fitpx, biplot = TRUE)
abline(h = 0, v = 0, lty=2)

cr <- getResidualCor(fitpx)
#install.packages("corrplot")
library(corrplot)
plot.new()
corrplot(cr, diag = FALSE, type = "lower", method = "square", tl.srt = 25)


X <- menu_gllvm_data$x[,2:3]
fitpx1 <- gllvm(menu_gllvm_data$abund, X, family = poisson(), num.lv = 1)
coefplot(fitpx1, mfrow = c(1,2), cex.ylab = 0.8)
ggsave(filename = "../11_Vegan/Menu_combined2/Tide_gllvm_coeff.jpg")
coefplot(fitpx1, mfrow = c(2,3), cex.ylab = 0.8)
ggsave(filename = "../11_Vegan/Menu_combined2/Microhabitat_and_Tide_gllvm_coeff.jpg")

#the coefplots all crossed 0 so are not significant

##Include trait information in the gllvm ----

#https://jenniniku.github.io/gllvm/articles/vignette1.html#incorporating-functional-traits-into-fourth-corner-models
criteria <- NULL
for(i in 1:5){
  fiti <- gllvm(menu_gllvm_data$abund, menu_gllvm_data$x, family = poisson(), num.lv = i, sd.errors = FALSE,
                formula = ~ Site_Name + Microhabitat + Water_temp, seed = 1234)
  criteria[i] <- summary(fiti)$AICc
  names(criteria)[i] = i
}

criteria
#1 latent variable had the lowest AIC criteria number, so we'll use that.

fit_env <- gllvm(menu_gllvm_data$abund, menu_gllvm_data$x, family = poisson(), num.lv = 1,
                 formula = ~Site_Name+Tide+Microhabitat, seed = 1234)
dev.off()
coefplot(fit_env, cex.ylab = 0.7, mfrow=c(2,3),mar = c(4,6,2,1))
ggsave(filename = "../11_Vegan/GLLVM_species_microhabitat_coeffs.jpg")

# Fit GLLVM without environmental variables and 1 latent variable:
fit1lv <- gllvm(menu_gllvm_data$abund, family = poisson(), num.lv = 1, seed = 1234)

# Correlation matrix
library(gclus)
cr0 <- getResidualCor(fit1lv)
dev.off()
corrplot(cr0[order.single(cr0), order.single(cr0)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")

##correlation plots from gllvm----
# Residual correlation matrix:
cr <- getResidualCor(fit_env)
library(corrplot); library(gclus)
dev.off()
corrplot(cr[order.single(cr), order.single(cr)], diag = FALSE, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "red")



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

pairwise_landes_distances<-as.matrix(print(phyloseq::distance(Shark_Menu_taxa_merge_species_named.hell.ps, method="l", binary=F, type= "samples", upper=TRUE, diag=TRUE)))

## Site_Name*Microhabitat----
set.seed(200)
landes_permanova<- adonis2(as.dist(pairwise_landes_distances)~Shark_Menu_taxa_merge_species_named.hell.ps@sam_data$Site_Name*Shark_Menu_taxa_merge_species_named.hell.ps@sam_data$Microhabitat)

landes_permanova
write_delim(landes_permanova, file = "../11_Vegan/Menu_combined2/landes_permanova_results.txt")

#Site is significant p=0.001, Microhabitat p=0.022 less so.


#following https://github.com/pmartinezarbizu/pairwiseAdonis

landes_dist_matrix <- phyloseq::distance(Shark_Menu_taxa_merge_species_named.hell.ps, method ="l")

set.seed(200)
vegan::adonis2(landes_dist_matrix ~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named.hell.ps)$Site_Name+Shark_Menu_taxa_merge_species_named.hell.ps@sam_data$Microhabitat+Shark_Menu_taxa_merge_species_named.hell.ps@sam_data$Tide)
#Site p=0.001, Microhabitat p=0.004, Tide 0.001

###Just site_name----
set.seed(2)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named.hell.ps)$Site_Name)

set.seed(2)
site_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named.hell.ps)$Site_Name)
write_delim(site_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/site_landes_pairwise_adonis_results.txt")


###just microhabitat----
set.seed(2)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named.hell.ps)$Microhabitat)

set.seed(2)
microhabitat_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named.hell.ps)$Microhabitat)
write_delim(landes_pairwise_adonis, "../11_Vegan/Menu_combined2/microhabitat_landes_pairwise_adonis_results.txt")

#mangroves vs. beach were significant at and middle of bay vs beach was significant

##Check for differences between hammerhead sites and non-hammerhead sites.----
#First remove hammerheads so they don't confound the calculation.
Shark_Menu_taxa_merge_species_named_traits_rm_hh.hell.ps <- Shark_Menu_taxa_merge_species_named_traits.hell.ps %>% subset_taxa(species != "Sphyrna lewini")

set.seed(2)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.hell.ps)$hammerhead_site)

set.seed(2)
hhsite_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.hell.ps)$hammerhead_site)
write_delim(hhsite_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/hhsite_landes_pairwise_adonis_results.txt")
#significant 0.001

set.seed(2)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.hell.ps)$hammerheads_detected)

set.seed(2)
hhsite_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.hell.ps)$hammerheads_detected)
write_delim(hhsite_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/hhdetected_landes_pairwise_adonis_results.txt")
#also significant

cca_ps_rm_hh <- ordinate( physeq = Shark_Menu_taxa_merge_species_named_traits_rm_hh.hell.ps,  method = "CCA", formula = ~ Microhabitat+Tide+hammerhead_site+hammerheads_detected + Condition(Site_Name), scale=TRUE)

summary(cca_ps_rm_hh)
anova(cca_ps_rm_hh, by = "margin")

# Model: cca(formula = OTU ~ Microhabitat + Tide + hammerhead_site + hammerheads_detected + Condition(Site_Name), data = data)
# Df ChiSquare      F Pr(>F)    
# Microhabitat          2    0.2420 1.5032  0.001 ***
#   Tide                  2    0.1660 1.0312  0.262    
# hammerhead_site       0    0.0000    Inf  1.000    
# hammerheads_detected  1    0.0954 1.1850  0.204    
# Residual             47    3.7830  

#When controlling for Site_Name, Microhabitat was significant, hammerheads_detected was significant, and hammerhead_site couldn't be tested, it seems.  

##Check for differences between Carcharhinus sites and non-carcharhinus sites.----
#First remove Carcharhinus reads so they don't confound the calculation.
Shark_Menu_taxa_merge_species_named_traits_rm_carch.hell.ps <- Shark_Menu_taxa_merge_species_named_traits.hell.ps %>% subset_taxa(genus != "Carcharhinus")

set.seed(2)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.hell.ps)$carcharhinus_detected)

# pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 0 vs 1  1  157.3002 1.767298 0.0305934   0.077      0.077 

set.seed(2)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.hell.ps)$blacktip_site)
# pairs Df SumsOfSqs F.Model         R2 p.value p.adjusted sig
# 1 1 vs 0  1  206.2803  2.3406 0.04011957   0.028      0.028   .

set.seed(2)
btsite_landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.hell.ps)$blacktip_site)
write_delim(btsite_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/btsite_landes_pairwise_adonis_results.txt")

cca_ps_rm_carch <- ordinate( physeq = Shark_Menu_taxa_merge_species_named_traits_rm_carch.hell.ps,  method = "CCA", formula = ~ Microhabitat+Tide+blacktip_site+carcharhinus_detected + Condition(Site_Name), scale=TRUE)

summary(cca_ps_rm_carch)
anova(cca_ps_rm_carch, by = "margin")

# Model: cca(formula = OTU ~ Microhabitat + Tide + blacktip_site + carcharhinus_detected + Condition(Site_Name), data = data)
# Df ChiSquare      F Pr(>F)    
# Microhabitat           2    0.2422 1.5100  0.001 ***
#   Tide                   2    0.1689 1.0533  0.239    
# blacktip_site          0    0.0000    Inf  1.000    
# carcharhinus_detected  1    0.0852 1.0623  0.293    
# Residual              47    3.7689                  
# ---


##########################

#make a network graph----

## total dataset----

plot_net(Shark_Menu_taxa_merge_species_named_traits.hell.ps, distance = "l", type ="taxa", maxdist = 0.09,
         laymeth = "circle", color = "class" , shape = NULL,
         rescale = TRUE, point_size = 5, point_alpha = .9, point_label = "species",
         hjust = 1.35, title = NULL)

#figure out how to control for site and environmental variables.



ig <- make_network(Shark_Menu_taxa_merge_species_named_traits.hell.ps, max.dist=0.9)
plot_network(ig, Shark_Menu_taxa_merge_species_named_traits.hell.ps,label = "Site_Name")

Menu.ps<-ps_filter(Menu.ps, Microhabitat != "passive mid bay")



#Plot richness----
plot_richness(Menu.ps, x="Site_Name","Microhabitat", measures=c("Shannon", "Simpson"), color="Site_Name")+
  scale_color_ucscgb()


ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Alpha_Diversity_untransformed_ASVs.jpg"))


# Transform data to proportions as appropriate for Bray-Curtis distances proportion of total reads per sample
Menu.ps.prop <- transform_sample_counts(Menu.ps, function(otu) otu/sum(otu))

##View the Top 50 most abundant ASVs in the dataset----

#take just the top 50 ASVs
top50 = names(sort(taxa_sums(Menu.ps), decreasing = TRUE)[1:50])
top50_Menu.ps = prune_taxa(top50, Menu.ps)
plot_bar(top50_Menu.ps, x="Site_Name",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("Read Abundance")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Top50_ASVs_by_site_and_genus.jpg"))

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

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Frequency_of_species_occurrence_by_Site.jpg"))

##find the top 20 most frequently occurring species----

top20occ = names(sort(taxa_sums(Menu_species_merge1.pa.ps), decreasing = TRUE)[1:20])
top20occ_Menu_species_merge1.pa.ps = prune_taxa(top20occ, Menu_species_merge1.pa.ps)

plot_bar(top20occ_Menu_species_merge1.pa.ps, x="Site_Name",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")+
  ylab("Frequency of Occurrence")

#Merge samples by site ----
Menu_site_merge_taxa_merge.ps <- merge_samples(Menu_taxa_merge.ps, "Site_Name")

plot_bar(Menu_site_merge_taxa_merge.ps, fill = "family") + 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")+
  ylab("Relative Read Abundance")




# 
#plot percent of occurrence (POO) by genus and site----
#so that ASVs found in more replicates in a given site have a higher POO value than ones found in fewer replicates.
#the 100% value will be number of positive detections for each site, and the values for each taxon in the list will be the rowsum of the presence-absence transformed otu table grouped by site.
#remove sites with less than (cutoff) sample replicates?



#Plot percent of occurrence (relative to occurrences per site)----


##(Menu_species_merge1.pa.ps is just the ASVs identified to species level)

###Use Shark_Menu_taxa_merge_species_named_traits.hell.ps

#Take the phyloseq tables into dataframes and visualize

otu_df<-cbind(data.frame(otu_table(Shark_Menu_taxa_merge_species_named_traits.hell.ps)))
otu_df$otu<-rownames(otu_df)

otu_df_long<-otu_df%>%
  pivot_longer(cols = -otu,names_to = "sample",values_to = "presence")
otu_df_long$sample <-gsub(pattern = ".",
                          "-",
                          x=otu_df_long$sample,
                          fixed = TRUE)



sample_df<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.hell.ps)))
sample_df<-rownames_to_column(sample_df, var = "sample")
sample_site_key<-sample_df%>%
  select(sample,Site_Name)

taxa_df<-cbind(data.frame(tax_table(Shark_Menu_taxa_merge_species_named_traits.hell.ps)))
taxa_df<-rownames_to_column(taxa_df,var="otu")
otu_species_key<-taxa_df%>%
  select(otu,species)

otu_df_long<-otu_df_long%>%
  left_join(sample_site_key)%>%
  left_join(otu_species_key)

merge_site_otu_counts<-otu_df_long%>%
  group_by(Site_Name,species)%>%
  summarise(total_count=sum(presence))%>%
  ungroup()

merge_site_counts<-merge_site_otu_counts%>%
  group_by(Site_Name)%>%
  mutate(total_site=sum(total_count))%>%
  ungroup()%>%
  mutate(rel_prct=(total_count/total_site)*100)


#plot 
ggplot(merge_site_counts,fill="species")+
  geom_bar(aes(x=Site_Name, y=rel_prct, color=species, fill=species), color="black",stat="identity", position="stack")+
  ylab("Relative Percent of Occurrence by Site")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Relative_percent_of_species_occurrence_by_Site.jpg"))

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
  select(sample,Site_Name)

taxa_df<-cbind(data.frame(tax_table(Menu_taxa_merge.pa.ps)))
taxa_df<-rownames_to_column(taxa_df,var="otu")
otu_species_key<-taxa_df%>%
  select(otu,family,genus,species)

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
  Menu_species_merge1.pa.ps,
  variable = "Site_Name",
  taxa_frequency = TRUE,
  drop_zeroes = TRUE,
  justdf = FALSE,
  long = FALSE)



plot_bar(POO,fill="genus")+
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
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
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
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


# Plot richness----
plot_richness(Menu_taxa_merge.ps, x="Water_temp", color="Microhabitat", measures=c("Simpson"))+
  scale_color_locuszoom()+
  scale_shape_girafe_filled() +
  ggplot2::stat_smooth(
    ggplot2::aes(colour = Microhabitat))


ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_alpha_diversity_by water_temp_and_Microhabitat.jpg"))

plot_richness(Menu_taxa_merge.ps, x="Site_Name", color="Site_Name", measures = c("Observed","Shannon","Simpson","InvSimpson","Chao1"))+ 
  geom_boxplot()+
  scale_color_aaas()
ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_alpha_diversity_by_site.jpg"))

plot_richness(Menu_taxa_merge.ps, x="Microhabitat", color="Microhabitat", measures = c("Observed","Shannon","Simpson","InvSimpson","Chao1"))+ 
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
# pseq <- pseq %>%
#   phyloseq_validate()
# 
# ord_explore(pseq) 
# 
# #code from the above program... 
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





