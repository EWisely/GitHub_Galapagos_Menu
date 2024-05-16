#Eldridge Wisely
#Date: 5-13-24
#Do permanova analysis of combined BerryCrust and MiFish phyloseq objects 

library(vegan)
library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(phyloseq)
library(metagMisc)
library(tidyverse)
#remotes::install_github("mikemc/speedyseq")
library(speedyseq)
#install_github("kpmainali/CooccurrenceAffinity")
library(CooccurrenceAffinity)
library(cooccur)
library(iNEXT)
#remotes::install_github("ropensci/rfishbase")
library(rfishbase)
library(microViz)
library(ggplot2)



Primer<-"Menu_combined2" #options Menu_combined has empty PCRs, Menu_combined1 has had the empty PCRs removed, and Menu_combined2 is aggregated by sample (PCRs from the same sample combined), and empty samples removed, phylogenetic tree added.

Menu.ps<-readRDS(file= paste("../09_Metabar_to_Phyloseq/",Primer,"_Phyloseq.rds", sep = ""))

Menu.ps

#filter out passive mid-bay samples because they didn't seem to work well
Menu1.ps<- ps_filter(Menu.ps, Microhabitat != "passive mid bay")

Menu1.ps

Passive.ps<-ps_filter(Menu.ps, Microhabitat=="passive mid bay")

#Merge_species but keep taxa that aren't ID'ed to species----
Menu_taxa_merge.ps <- phyloseq::tax_glom(Menu1.ps, taxrank = "species", NArm = FALSE)
#Merge_species and remove taxa that aren't ID'ed to species----
Menu_species_merge.ps<-phyloseq::tax_glom(Menu1.ps, taxrank = "species", NArm=TRUE)


#convert to presence absence----


Menu_taxa_merge.pa.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge.ps, method = "pa")

Menu_species_merge.pa.ps<-phyloseq_standardize_otu_abundance(Menu_species_merge.ps, method = "pa")

#distance measures of taxa----
#calculate all the distance measures

dist_methods <- unlist(distanceMethodList)
print(dist_methods)
#Remove the distance methods that require a phylogenetic tree (for now)
#dist_methods <- dist_methods[-(1:3)]
# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods=="ANY")]
#loop over each distance method offered in Phyloseq
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- phyloseq::distance(Menu1.ps, method=i)
  # Calculate ordination
  iMDS  <- ordinate(Menu1.ps, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(Menu1.ps, iMDS, color="Site_Name", shape="Microhabitat")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

library("plyr")
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Site_Name, shape=Microhabitat))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for MiFish Menu dataset")
p

#sim looks good, binary, and jaccard and bray... l is Lande's index (Lande,1996 citation in Zotero) 
#betadiver indices come from Measuring beta diversity for presence–absence data, Patricia Koleff, Kevin J. Gaston, Jack J. Lennon 2003

print(plist[["sim"]])#axis1 26.1%, axis2 19.1%
ggsave(paste0("../11_Vegan/",Primer,"/MDS using distance method sim.jpg"))
print(plist[["jsd"]])#axis1 28.2, axis2 15.5%
ggsave(paste0("../11_Vegan/",Primer,"/MDS using distance method jsd.jpg"))
print(plist[["bray"]])#axis1 17.8%,axis2 10.7%
ggsave(paste0("../11_Vegan/",Primer,"/MDS using distance method bray.jpg"))
print(plist[["jaccard"]])#axis1 8.8%,axis2 5.6%
ggsave(paste0("../11_Vegan/",Primer,"/MDS using distance method jaccard.jpg"))
print(plist[["binary"]])#axis1 10.2%, axis2 8.1%
ggsave(paste0("../11_Vegan/",Primer,"/MDS using distance method binary.jpg"))
print(plist[["unifrac"]])#axis1 30.4%, axis2 12.2%
ggsave(paste0("../11_Vegan/",Primer,"/MDS using distance method unifrac.jpg"))
print(plist[["wunifrac"]])#axis1 60.7%, axis2 14.3%
ggsave(paste0("../11_Vegan/",Primer,"/MDS using distance method wunifrac.jpg"))
print(plist[["raup"]])#axis1 87.9%, axis2 71.9%
ggsave(paste0("../11_Vegan/",Primer,"/MDS using distance method raup.jpg"))


#combine results and shade according to Microhabitat
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Microhabitat, shape=Site_Name))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for MiFish Menu dataset")
p


##calculate Raup-Crick distances and put them in a symmetrical matrix (presence/absence)----
pairwise_raupcrick_distances<-as.matrix(print(phyloseq::distance(Menu_taxa_merge.ps, method="raup", binary=TRUE, type= "samples", upper=TRUE, diag=TRUE)))


#Yay, the matrix is symmetrical!
set.seed(200)
raupcrick_permanova<- adonis2(pairwise_raupcrick_distances~Menu_taxa_merge.ps@sam_data$Site_Name)

#following https://github.com/pmartinezarbizu/pairwiseAdonis

raupcrick_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="raup")

set.seed(200)
vegan::adonis2(raupcrick_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

#Different p-values, but none significant

#Calculate sim distance (Simpson?)
sim_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="sim")

set.seed(200)
sim_permanova<-vegan::adonis2(sim_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

write_delim(sim_permanova, file = "../11_Vegan/Menu_combined2/sim_permanova_results.txt")
#with sim distances, both site p=0.012 and microhabitat p=0.002 are significant 

pairwise_landes_distances<-as.matrix(print(phyloseq::distance(Menu_taxa_merge.ps, method="l", binary=T, type= "samples", upper=TRUE, diag=TRUE)))

set.seed(200)
landes_permanova<- adonis2(as.dist(pairwise_landes_distances)~Menu_taxa_merge.ps@sam_data$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

write_delim(landes_permanova, file = "../11_Vegan/Menu_combined2/landes_permanova_results.txt")

#Site is significant p=0.001, Microhabitat p=0.022 less so.


#following https://github.com/pmartinezarbizu/pairwiseAdonis

landes_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="l")

set.seed(200)
vegan::adonis2(landes_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)
#Site p=0.001, Microhabitat p=0.022

set.seed(200)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name)
set.seed(200)
landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name)
write_delim(landes_pairwise_adonis, "../11_Vegan/Menu_combined2/landes_pairwise_adonis_results.txt")

#interestingly, Puerto Grande vs. RB1&2, and Cerro Brujo are the only ones that are significant and even then, barely

#Check if water temperature effects presence absence of all taxa in each PCR replicate
#remove PCRs that don't have temperature recorded
Menu_taxa_merge.temp.ps <- ps_filter(Menu_taxa_merge.ps, Water_temp != "NA")

#make 7 bins for temps 

Menu_taxa_merge.temp.ps<- mutate_sample_data(Menu_taxa_merge.temp.ps, temp_bin=cut(Water_temp, breaks=6))

#redo distance matrix... do Raup-Crick indices
#Raup-Crick indices for presence–absence data should be able to handle unknown (and variable) sample sizes. Most of these indices are discussed by Krebs (1999) and Legendre & Legendre (2012) -vegan package vegdist help page

pairwise_raupcrick_distances<-as.matrix(print(phyloseq::distance(Menu_taxa_merge.temp.ps, method="raup", binary=TRUE, type= "samples", upper=TRUE, diag=TRUE)))

#do permanova
set.seed(200)
temp_raup_permanova<-adonis2(as.dist(pairwise_raupcrick_distances)~Menu_taxa_merge.temp.ps@sam_data$temp_bin)

write_delim(temp_raup_permanova, "../11_Vegan/Menu_combined2/raup-crick_permanova_temp_bins.txt")

#p-value is 0.165 with 5 breaks, p=0.03 with 6 breaks

#following https://github.com/pmartinezarbizu/pairwiseAdonis

unifrac_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="unifrac")
set.seed(200)
vegan::adonis2(unifrac_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

#the unifrac is presence/absence and takes phylogeny into account
set.seed(200)
pairwise.adonis(unifrac_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Tide)
#tide insignificant in all pairwise comparisons
pairwise.adonis(unifrac_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Microhabitat)
#microhabitat insignificant in all pairwise comparisons
pairwise.adonis(unifrac_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name)
#site insignificant in all pairwise comparisons

##calculate Jsd distances and put them in a symmetrical matrix taking read abundance into account----
pairwise_jsd_distances<-as.matrix(print(phyloseq::distance(Menu_taxa_merge.ps, method="jsd", binary=FALSE, type= "samples", upper=TRUE, diag=TRUE)))

set.seed(200)
jsd_permanova<- adonis2(as.dist(pairwise_jsd_distances)~Menu_taxa_merge.ps@sam_data$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat*Menu_taxa_merge.ps@sam_data$Tide)

write_delim(jsd_permanova, "../11_Vegan/Menu_combined2/jsd_permanova_results.txt")

#when taking read abundances into account using jsd , site p-value is 0.001, Microhabitat p-value is 0.001, Tide is 0.348(not significant)

set.seed(200)
pairwise.adonis(pairwise_jsd_distances, phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name)
set.seed(200)
pairwise_jsd_results_site<-pairwise.adonis(pairwise_jsd_distances, phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name)
#Puerto Grande vs Rosa Blanca 1 p=0.021, all others not significant
write_delim(pairwise_jsd_results_site, "../11_Vegan/Menu_combined2/jsd_pairwise_results_by_site.txt")
##############CoocurrenceAffinity package #################

#Use CoocurrenceAffinity package to look for groups of taxa that occur together more than could be expected by chance----

# compute the affinity between elements in rows (= Taxa)----

taxa_affinity <- affinity(data = Menu_taxa_merge.pa.ps@otu_table, row.or.col = "row", squarematrix = c("all"))
plotgg(data = taxa_affinity, variable = "alpha_mle", legendlimit = "datarange")

##Merge_families and remove taxa that aren't ID'ed to family----
Menu_family_merge.ps<-phyloseq::tax_glom(Menu.ps, taxrank = "family", NArm=TRUE)

###convert to presence absence----
Menu_family_merge.pa.ps<-phyloseq_standardize_otu_abundance(Menu_family_merge.ps, method = "pa")


family_matrix<-as(otu_table(Menu_family_merge.pa.ps),"matrix")
family_names<-as.data.frame(tax_table(Menu_family_merge.pa.ps))%>%
  select("family")
family_table<- as.data.frame(family_matrix)
rownames(family_table)<-family_names$family



family_affinity<-affinity(data = family_table, row.or.col = "row", squarematrix = c("all"))
plotgg(data = family_affinity, variable = "alpha_mle", legendlimit = "datarange", text.size = 3.5)

ggsave(paste0("../11_Vegan/",Primer,"/",Primer,"_family_coocurrence_affinity_plot.jpg"))


##Then just the ones ID'ed to species----


species_matrix = as(otu_table(Menu_species_merge.pa.ps), "matrix")
species_names<-as.data.frame(tax_table(Menu_species_merge.pa.ps))%>%
  select("species")
species_table = as.data.frame(species_matrix)
rownames(species_table)<-species_names$species



species_affinity<-affinity(data = species_table, row.or.col = "row", squarematrix = c("all"))
plotgg(data = species_affinity, variable = "alpha_mle", legendlimit = "datarange")

ggsave(paste0("../11_Vegan/",Primer,"/",Primer,"_species_coocurrence_affinity_plot.jpg"))



#Cooccur package----

#make a species x site matrix
#merge samples by Site_Name

Menu_taxa_merge_site_merge.pa.ps<-merge_samples(Menu_taxa_merge.pa.ps,group = "Site_Name")

#Take the phyloseq tables into dataframes and visualize

otu_df<-cbind(data.frame(otu_table(Menu_taxa_merge_site_merge.pa.ps)))
#otu_df$site<-rownames(otu_df)


taxa_df<-cbind(data.frame(tax_table(Menu_taxa_merge_site_merge.pa.ps)))
taxa_df<-rownames_to_column(taxa_df,var="otu")
otu_species_key<-taxa_df%>%
  select(otu,family,genus,species)

otu_df_t<-t(otu_df)
otu_df_t<-as.data.frame(otu_df_t)
otu_df_t$otu<-rownames(otu_df_t)
otu_df_t<-left_join(otu_df_t,otu_species_key)

species_by_site<-otu_df_t%>%
  select(-c("genus","family","otu"))%>%
  filter(species!="NA")%>%
  column_to_rownames(var="species")

species_by_site_mat<-as.matrix(species_by_site)
species_cooccur<-cooccur(species_by_site_mat, spp_names = TRUE)

class(species_cooccur)
summary(species_cooccur)

species_cooccur_prob_table<-prob.table(species_cooccur)

plot(species_cooccur) # add "plotrand = TRUE" to include completely random species
ggsave("../11_Vegan/Menu_combined2/species_cooccurrence_matrix")

hammerhead_cooccurring_sp<-pair(mod = species_cooccur, spp = "Sphyrna lewini")

pair.profile(mod = species_cooccur)

#do the same with genus, starting with tax_glom or phyloseq equivalent

#Hill Numbers----
#phyloseq_inext(Menu.ps)

#rFishbase----

sp_valid<-validate_names(
  species_names$species, 
  server = "fishbase")

# cr_valid<-validate_names(
#   species_names$species,
#   server ="sealifebase"
# )

print(sp_valid)

ecology_traits<-ecology(species_list = species_names)

#Find which samples had hammerheads and which had blacktips----

Blacktips.ps <- subset_taxa(Menu_taxa_merge.ps, genus =="Carcharhinus")
Blacktips.ps <- prune_samples(sample_sums(Blacktips.ps)>0, Blacktips.ps)
Blacktip_sites<-Blacktips.ps@sam_data$Site_Name
Blacktip_samples<-Blacktips.ps@sam_data$sample_id
Blacktip_samples.ps<-prune_samples(x = Menu_taxa_merge.ps,samples = Blacktip_samples)
Blacktip_samples.ps<-prune_taxa(taxa_sums(Blacktip_samples.ps)>0,Blacktip_samples.ps)

saveRDS(Blacktip_samples.ps, file=paste0("../11_Vegan/Menu_combined2/Blacktip_samples_Phyloseq.rds"))


Hammerheads.ps<-subset_taxa(Menu_taxa_merge.ps,genus=="Sphyrna")
Hammerheads.ps <- prune_samples(sample_sums(Hammerheads.ps)>0, Hammerheads.ps)
Hammerhead_sites<-Hammerheads.ps@sam_data$Site_Name
Hammerhead_samples<-Hammerheads.ps@sam_data$sample_id
Hammerhead_samples.ps<-prune_samples(x = Menu_taxa_merge.ps,samples = Hammerhead_samples)
Hammerhead_samples.ps<-prune_taxa(taxa_sums(Hammerhead_samples.ps)>0,Hammerhead_samples.ps)

saveRDS(Hammerhead_samples.ps, file=paste0("../11_Vegan/Menu_combined2/Hammerhead_samples_Phyloseq.rds"))


