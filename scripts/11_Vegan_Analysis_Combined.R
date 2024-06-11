#Eldridge Wisely
#Date: 5-13-24
#Do permanova analysis of combined BerryCrust and MiFish phyloseq objects 

library(plyr)
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
library(ggvegan)





Primer<-"Menu_combined2" #options Menu_combined has empty PCRs, Menu_combined1 has had the empty PCRs removed, and Menu_combined2 is aggregated by sample (PCRs from the same sample combined), and empty samples removed, phylogenetic tree added.

Menu.ps<-readRDS(file= paste("../09_Metabar_to_Phyloseq/",Primer,"_Phyloseq.rds", sep = ""))

Menu.ps

#filter out passive mid-bay samples because they didn't seem to work well
Menu1.ps<- ps_filter(Menu.ps, Microhabitat != "passive mid bay")

Menu1.ps

Passive.ps<-ps_filter(Menu.ps, Microhabitat=="passive mid bay")

#Merge_species but keep taxa that aren't ID'ed to species----
Menu_taxa_merge.ps <- phyloseq::tax_glom(Menu1.ps, taxrank = "species", NArm = FALSE)

#rename with species names if available
Menu_taxa_merge_tax_fix.ps<-Menu_taxa_merge.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

Menu_taxa_merge_species_named.ps<-tax_rename(Menu_taxa_merge_tax_fix.ps,rank = "species")
Menu_taxa_merge_family_named.ps<-tax_rename(Menu_taxa_merge_tax_fix.ps,rank = "family",)




#Merge_species and remove taxa that aren't ID'ed to species----
Menu_species_merge.ps<-phyloseq::tax_glom(Menu1.ps, taxrank = "species", NArm=TRUE)
#tax_fix_interactive(Menu_species_merge.ps)
Menu_species_merge_tax_fix.ps<-Menu_species_merge.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

Menu_species_merge_species_named.ps<-tax_rename(Menu_species_merge_tax_fix.ps,rank = "species")

#convert to presence absence----

Menu_taxa_merge.pa.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge.ps, method = "pa")

#and with the species named objects
Menu_taxa_merge_species_named.pa.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named.ps, method = "pa")

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

#library("plyr")
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Site_Name, shape=Microhabitat))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Menu dataset")
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
print(plist[["-3"]])#axis1 87.9%, axis2 71.9%
ggsave(paste0("../11_Vegan/",Primer,"/MDS using distance method -3.jpg"))

#combine results and shade according to Microhabitat
df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Microhabitat, shape=Site_Name))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Menu dataset")
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

##Calculate unifrac distance matrix----
unifrac_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="unifrac")

set.seed(200)
vegan::adonis2(unifrac_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

#not significant with unifrac either.

##Calculate sim distance (Simpson?)----
sim_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="sim")

set.seed(200)
sim_permanova<-vegan::adonis2(sim_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

write_delim(sim_permanova, file = "../11_Vegan/Menu_combined2/sim_permanova_results.txt")
#with sim distances, both site p=0.012 and microhabitat p=0.002 are significant 

##calculate Landes distances----

pairwise_landes_distances<-as.matrix(print(phyloseq::distance(Menu_taxa_merge.ps, method="l", binary=T, type= "samples", upper=TRUE, diag=TRUE)))

set.seed(200)
landes_permanova<- adonis2(as.dist(pairwise_landes_distances)~Menu_taxa_merge.ps@sam_data$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

landes_permanova
write_delim(landes_permanova, file = "../11_Vegan/Menu_combined2/landes_permanova_results.txt")



#Site is significant p=0.001, Microhabitat p=0.022 less so.


#following https://github.com/pmartinezarbizu/pairwiseAdonis

landes_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="l")

set.seed(200)
vegan::adonis2(landes_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)
#Site p=0.001, Microhabitat p=0.022

set.seed(2)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name)

set.seed(2)
landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name)
write_delim(landes_pairwise_adonis, "../11_Vegan/Menu_combined2/landes_pairwise_adonis_results.txt")

#interestingly, Puerto Grande vs. RB1&2, and Puerto Grande vs. Cartago Bay are the only ones that are significant and even then, barely

set.seed(2)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Microhabitat)

set.seed(2)
landes_pairwise_adonis<-pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Microhabitat)
write_delim(landes_pairwise_adonis, "../11_Vegan/Menu_combined2/landes_pairwise_adonis_results_microhabitat.txt")

#mangroves vs. beach were significant at p-adj 0.024 but the effect was fairly small at R2=0.09468624

###landes distances with presence absence data----
pa_landes_dist_matrix <- phyloseq::distance(Menu_taxa_merge.pa.ps, method ="l")

set.seed(200)
vegan::adonis2(pa_landes_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.pa.ps)$Site_Name*Menu_taxa_merge.pa.ps@sam_data$Microhabitat)
#Site p=0.001, Microhabitat p=0.037

set.seed(2)
pairwise.adonis(pa_landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.pa.ps)$Site_Name)

set.seed(2)
pa_landes_pairwise_adonis<-pairwise.adonis(pa_landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.pa.ps)$Site_Name)

write_delim(pa_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/pa_landes_pairwise_adonis_results.txt")

#Puerto Grande vs. RB1&2, and Puerto Grande vs. Cartago Bay are the only ones that are significant and the effect of Site remains the same between presence/absence and read abundance data with the landes distances!

set.seed(2)
pairwise.adonis(pa_landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.pa.ps)$Microhabitat)

set.seed(2)
pa_landes_pairwise_adonis<-pairwise.adonis(pa_landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.pa.ps)$Microhabitat)
write_delim(pa_landes_pairwise_adonis, "../11_Vegan/Menu_combined2/landes_pairwise_adonis_results_microhabitat.txt")

#mangroves vs. beach were significant at p-adj 0.024 but the effect was fairly small at R2=0.09468624 #this effect remained the same between read abundance and presence/absence data for landes distances.





#Check if water temperature effects presence absence of all taxa in each PCR replicate----
#remove PCRs that don't have temperature recorded
Menu_taxa_merge.temp.ps <- ps_filter(Menu_taxa_merge.ps, Water_temp != "NA")

#make 6 bins for temps 

Menu_taxa_merge.temp.ps<- mutate_sample_data(Menu_taxa_merge.temp.ps, temp_bin=cut(Water_temp, breaks=6))

#redo distance matrix... do Raup-Crick indices
#Raup-Crick indices for presence–absence data should be able to handle unknown (and variable) sample sizes. Most of these indices are discussed by Krebs (1999) and Legendre & Legendre (2012) -vegan package vegdist help page

pairwise_raupcrick_distances<-as.matrix(print(phyloseq::distance(Menu_taxa_merge.temp.ps, method="raup", binary=TRUE, type= "samples", upper=TRUE, diag=TRUE)))

#do permanova
set.seed(200)
tempbin_raup_permanova<-adonis2(as.dist(pairwise_raupcrick_distances)~Menu_taxa_merge.temp.ps@sam_data$temp_bin)

tempbin_raup_permanova

write_delim(tempbin_raup_permanova, "../11_Vegan/Menu_combined2/raup-crick_permanova_temp_bins.txt")

set.seed(200)
temp_raup_permanova<-adonis2(as.dist(pairwise_raupcrick_distances)~Menu_taxa_merge.temp.ps@sam_data$Water_temp)

temp_raup_permanova

write_delim(temp_raup_permanova, "../11_Vegan/Menu_combined2/raup-crick_permanova_Water_temp.txt")


#p-value is 0.165 with 5 breaks, p=0.03 with 6 breaks, 0.074 with no breaks

#following https://github.com/pmartinezarbizu/pairwiseAdonis

unifrac_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="unifrac")
set.seed(200)
vegan::adonis2(unifrac_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

#the unifrac is presence/absence and takes phylogeny into account----
set.seed(200)
pairwise.adonis(unifrac_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Tide)
#tide insignificant in all pairwise comparisons
pairwise.adonis(unifrac_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Microhabitat)
#microhabitat insignificant in all pairwise comparisons
pairwise.adonis(unifrac_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name)
#site insignificant in all pairwise comparisons

#check temp
temp_unifrac_dist_matrix <- phyloseq::distance(Menu_taxa_merge.temp.ps, method ="unifrac")
pairwise.adonis(temp_unifrac_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.temp.ps)$temp_bin)
#temperature bin not significant in unifrac pairwise comparisons


landes_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="l")
set.seed(200)
vegan::adonis2(landes_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

#the landes is robust to presence absence and factors rarer taxa more than some other distance matrices.----

landes_dist_matrix <- phyloseq::distance(Menu_taxa_merge.ps, method ="l")
set.seed(200)
vegan::adonis2(landes_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name*Menu_taxa_merge.ps@sam_data$Microhabitat)

set.seed(200)
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Tide)
#tide significant incoming/low vs incoming/high R2=0.26753111 padj=0.021
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Microhabitat)
#microhabitat significant mangroves vs beach 0.033
pairwise.adonis(landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.ps)$Site_Name)
#same as above site landes comparison

#check temp
temp_landes_dist_matrix <- phyloseq::distance(Menu_taxa_merge.temp.ps, method ="l")
pairwise.adonis(temp_landes_dist_matrix, phyloseq::sample_data(Menu_taxa_merge.temp.ps)$temp_bin)
#temperature bin not significant in landes pairwise comparisons

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

#Just get the ones that don't show up in every sample----
Nonshared_Menu_taxa_merge_species_named.ps<-phyloseq_extract_non_shared_otus(Menu_taxa_merge_species_named.ps)

nonshared_raupcrick_dist_matrix <- phyloseq::distance(Nonshared_Menu_taxa_merge_species_named.ps, method ="raup")

set.seed(200)
vegan::adonis2(nonshared_raupcrick_dist_matrix ~ phyloseq::sample_data(Nonshared_Menu_taxa_merge_species_named.ps)$Site_Name*Nonshared_Menu_taxa_merge_species_named.ps@sam_data$Microhabitat)

##Calculate landes distance matrix----
nonshared_pairwise_landes_distances<-as.matrix(print(phyloseq::distance(Nonshared_Menu_taxa_merge_species_named.ps, method="l", binary=T, type= "samples", upper=TRUE, diag=TRUE)))

set.seed(200)
nonshared_landes_permanova<- adonis2(as.dist(nonshared_pairwise_landes_distances)~Nonshared_Menu_taxa_merge_species_named.ps@sam_data$Site_Name*Nonshared_Menu_taxa_merge_species_named.ps@sam_data$Microhabitat)

nonshared_landes_permanova
write_delim(nonshared_landes_permanova, file = "../11_Vegan/Menu_combined2/nonshared_species_landes_permanova_results.txt")
#Interaction term of Site Name and Microhabitat significant 0.011, Site and Microhabitat not significant on their own when looking at non-shared taxa

#Vegan class stuff----


#use the Hellinger transformation for the ASV counts (square-root of relative abundances)

Menu_taxa_merge.hell.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named.ps, method = "hellinger")

#Try the Wisconsin transformation: species first standardized by maxima then sites by site totals.
Menu_taxa_merge.wisc.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named.ps, method = "wisconsin")

Menu_taxa_merge_species_named.pa.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named.ps, method = "pa")



#get rid of samples with no temp data
Menu_taxa_merge_temp.hell.ps<-ps_filter(Menu_taxa_merge.hell.ps, Water_temp != "NA")


cca_ps_site_name<-ordinate(physeq = Menu_taxa_merge.hell.ps, method = "CCA", formula = ~Site_Name)
summary(cca_ps_site_name)
anova(cca_ps_site_name)


cca_ps_taxa_merge <- ordinate( physeq = Menu_taxa_merge_temp.hell.ps,  method = "CCA", formula = ~ Microhabitat+Tide+Water_temp + Condition(Site_Name), scale=TRUE)

summary(cca_ps_taxa_merge)
anova(cca_ps_taxa_merge, by = "margin")

#microhabitat and water temp were significant when controlling for site, but not tide.

# Model: cca(formula = OTU ~ Microhabitat + Tide + Water_temp + Condition(Site_Name), data = data)
# Df ChiSquare      F Pr(>F)    
# Microhabitat  2   0.20839 1.3777  0.002 ** 
#   Tide          6   0.49539 1.0917  0.113    
# Water_temp    1   0.14455 1.9112  0.001 ***
#   Residual     41   3.10086                  
# ---

cca_ps_taxa_merge1 <- ordinate( physeq = Menu_taxa_merge_temp.hell.ps,  method = "CCA", formula = ~ Microhabitat+Water_temp + Condition(Site_Name+Tide), scale=TRUE)

summary(cca_ps_taxa_merge1)


# Call:
#   cca(formula = OTU ~ Microhabitat + Water_temp + Condition(Site_Name +      Tide), data = data) 
# 
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          4.6764    1.00000
# Conditioned    1.1950    0.25555
# Constrained    0.3805    0.08136
# Unconstrained  3.1009    0.66309
# 
# Eigenvalues, and their contribution to the scaled Chi-square 
# after removing the contribution of conditiniong variables
# 
# Importance of components:
#   CCA1    CCA2    CCA3    CA1     CA2     CA3
# Eigenvalue            0.21007 0.11410 0.05628 0.4422 0.23239 0.21670
# Proportion Explained  0.06034 0.03277 

anova(cca_ps_taxa_merge1, by = "margin")

# Model: cca(formula = OTU ~ Microhabitat + Water_temp + Condition(Site_Name + Tide), data = data)
# Df ChiSquare      F Pr(>F)    
# Microhabitat  2   0.20839 1.3777  0.003 ** 
#   Water_temp    1   0.14455 1.9112  0.001 ***
#   Residual     41   3.10086   



p_taxa_merge<-plot(cca_ps_taxa_merge1, scaling = "site", hill = TRUE, display= c("species","bp"))
identify(p_taxa_merge, what = "species")
#these are good ones to identify [1]   1   7  23  27  44  59  60  64 102 117 119
#ooh, especially these: [1]   7  10  30  49  59  60  64  90 102 112 117 119

#saved a screen shot of this as Best CCA of species-vs-temp_and_microhabitat.png in 11_Vegan folder. 

#how can I automate this part?


#make a prettier plot with ggvegan?
ggvegan::ordiggplot(cca_ps_taxa_merge)+
  ggvegan::geom_ordi_point(data = cca_ps_taxa_merge)+
  ggvegan::geom_ordi_arrow(score = "ChiSquare")




#check if all of this holds true with the pa transformed object----
#get rid of samples with no temp data
Menu_taxa_merge_temp.pa.ps<-ps_filter(Menu_taxa_merge_species_named.pa.ps, Water_temp != "NA")


cca_ps_site_name_pa<-ordinate(physeq = Menu_taxa_merge_species_named.pa.ps, method = "CCA", formula = ~Site_Name+Microhabitat+Tide)
summary(cca_ps_site_name_pa)
anova(cca_ps_site_name_pa, by="margin")


cca_ps_taxa_merge_pa <- ordinate( physeq = Menu_taxa_merge_temp.pa.ps,  method = "CCA", formula = ~ Microhabitat+Water_temp+Tide + Condition(Site_Name), scale=TRUE)

summary(cca_ps_taxa_merge_pa)
anova(cca_ps_taxa_merge_pa, by = "margin")
###tide and water temp are significant for pa but not Microhabitat


cca_ps_taxa_merge_pa1 <- ordinate( physeq = Menu_taxa_merge_temp.pa.ps,  method = "CCA", formula = ~ Water_temp +Tide +Condition(Site_Name+Microhabitat), scale=TRUE)

summary(cca_ps_taxa_merge_pa1)

anova(cca_ps_taxa_merge_pa1, by = "margin")

p_taxa_merge_pa<-plot(cca_ps_taxa_merge_pa1, scaling = "site", hill = TRUE, display= c("species","bp"))
#identify(p_taxa_merge_pa, what = "species")

#Tide became significant with presence absence and I realized it's not making much sense.  I feel it should either be high/low or incoming/outgoing not various combinations of both.



#PCoA----
Pcoa<-ordinate(physeq = Menu_taxa_merge.hell.ps, method = , formula = ~Site_Name)
summary(cca_ps_site_name)
anova(cca_ps_site_name)

dpcoa_taxa_merge<-ordinate(physeq = Menu_taxa_merge_species_named.ps, method = "DPCoA")

summary(dpcoa_taxa_merge)
dpcoa_taxa_merge<-DPCoA(Menu_taxa_merge_species_named.ps,correction = "caillez")
#not sure how to plot this

##############CoocurrenceAffinity package #################

#Use CoocurrenceAffinity package to look for groups of taxa that occur together more than could be expected by chance----

# compute the affinity between elements in rows (= Taxa)----

taxa_affinity <- affinity(data = Menu_taxa_merge_species_named.pa.ps@otu_table, row.or.col = "row", squarematrix = c("all"))
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



plotgg(data = species_affinity, variable = "simpson", legendlimit = "datarange",text.size = 5)

plotgg(data = species_affinity, variable = "sorensen", legendlimit = "datarange",text.size = 5)

plotgg(data = species_affinity, variable = "jaccard", legendlimit = "datarange")

#plotgg(data = species_affinity, variable="p_value", legendlimit = "datarange")

#find which affinities have signficant p-values:
species_affinity_summary<-summary(species_affinity$all)

species_affinity_summary


#Create a site species table with the number of PCRs in which each species was present----
## Start with the original un-trimmed Menu_combined


Menu_combined.ps<-readRDS(file= paste("../09_Metabar_to_Phyloseq/Menu_Combined_Phyloseq.rds", sep = ""))

Menu_combined.ps

#compute by site after binary transformation of the ASVs assigned to species level leaving NAs
Menu_combined_taxa_merge.ps<-phyloseq::tax_glom(Menu_combined.ps, taxrank = "species", NArm = FALSE)
#transform to presence absence
Menu_combined_taxa_merge.pa.ps<-phyloseq_standardize_otu_abundance(Menu_combined_taxa_merge.ps, method = "pa")
#merge samples and sum
Menu_taxa_merge_site_merge.pa.ps<-merge_samples(Menu_combined_taxa_merge.pa.ps, group = "Site_Name",fun = "sum")
#the numbers of the otu table are now the number of PCRs that had that species or taxa present.
#make a matrix and dataframe of it.
taxa_site_matrix = as(otu_table(Menu_taxa_merge_site_merge.pa.ps), "matrix")
taxa_site_table = as.data.frame(taxa_site_matrix)

species_site_binary<-dataprep(data = taxa_site_matrix, row.or.col = "row", datatype = "abundance",
         threshold = 1, class0.rule = "less")

site_affinity<-affinity(data = species_site_binary, row.or.col = "col", squarematrix = "all",sigPval = "0.05")

plotgg(data = site_affinity, variable = "alpha_mle",legendlimit = "datarange")

species_site_affinity<-affinity(data = species_site_binary, row.or.col = "row",sigPval = "0.05",squarematrix = "all")

plotgg(data = species_site_affinity, variable = "alpha_mle",legendlimit = "datarange")

#that's all of the taxa.  How can I only show the ones with significant p-values?

#
species_affinity1 <- affinity(data = species_site_binary, row.or.col = "row", lev=0.95, pvalType="Blaker")

plotgg(species_affinity1, variable = "alpha_mle",legendlimit = "datarange")
       
species_aff_sig05<-subset(species_affinity1$all, subset = p_value<0.05)

species_aff_sig05
nrow(species_aff_sig05)



species_affinity1 <- affinity(data = species_site_binary, row.or.col = "row", which.row.or.col = sigrows, lev=0.95, pvalType="Blaker")


species_affinity1$all$p_value

species_aff_sig05$p_value

species_affinity1

#Maybe a better heatmap?----
#devtools::install_github("rlbarter/superheat")
library(superheat)


superheat(mtcars,
          # normalize variables
          scale = T,
          # order rows/cols based on heirarchical clustering
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          # plot miles per gallon next to the rows
          yr = mtcars$mpg,
          yr.axis.name = "miles per gallon",
          # plot correlation with mpg above columns
          yt = cor(mtcars)[, "mpg"],
          yt.plot.type = "bar",
          yt.axis.name = "correlation with mpg",
          # increase size of left labels
          left.label.size = 0.45)


mtcars
species_aff_sig05_select<-species_aff_sig05%>%
  select(where(is.numeric))
species_aff_sig05_select

heatmap<-superheat(species_aff_sig05_select,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
           # plot p-values next to the rows
            yr = species_aff_sig05_select$alpha_mle,
           yr.axis.name = "p-values",
          # plot alpha mle above columns
          yt = cor(species_aff_sig05_select)[, "p_value"], ,
          yt.plot.type = "bar",
          yt.axis.name = "alpha_mle",
          # increase size of left labels
          left.label.size = 0.45)


#Maybe give up on heatmaps and do a volcano plot instead?----


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
#This takes a while

otu_df_inext<-cbind(data.frame(otu_table(Menu.ps)))

#This takes a while
iNext_Menu<-iNEXT(otu_df_inext, q=0, datatype="abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)

ggiNEXT(iNext_Menu, type = 1, se = TRUE)+
  xlab("ASV Read Count by Sample")+
  ylab("ASV Richness")+
  theme(legend.position = 'none')

ggsave(filename = "../11_Vegan/Menu_combined2/iNext_ASV_richness_v_coverage_rarefaction_curve.jpg")

Menu_taxa_merge_species_named_site_sums.ps<-merge_samples(Menu_taxa_merge_species_named.ps, group = "Site_Name",fun = "sum")

otu_df_inext_by_site<-cbind(data.frame(otu_table(Menu_taxa_merge_species_named_site_sums.ps)))
otu_df_inext_by_site_trans<-t(otu_df_inext_by_site)

#create it once, it takes a long time.
iNext_Menu_by_site<-iNEXT(otu_df_inext_by_site_trans, q=c(0,1,2), datatype = "abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)

ggiNEXT(iNext_Menu_by_site, type = 2, se = TRUE,facet.var = "Order.q")+
  xlab("Taxon Read Count by Site")+
  ylab("Proportion of Total Taxonomic Richness")+
  xlim(0,30000)
ggsave("../11_Vegan/Menu_combined2/iNext_taxonomic_richness_by read_coverage_per_site.jpg")


#rFishbase----



sp_valid<-validate_names(
  species_names$species, 
  server = "fishbase")
print(sp_valid)

# cr_valid<-validate_names(
#   species_names$species,
#   server ="sealifebase"
# )


fb_tbl("species","sealifebase")

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
  rename(Mesopelagic =mesopelagic)

crustacean_ecology_traits1<-crustacean_ecology_traits%>%
  select(Species, FoodTroph, FoodSeTroph, Neritic, Intertidal, Oceanic, Epipelagic, Neritic, Mesopelagic, Estuaries, Mangroves, FeedingType, Solitary, AddRems)

Menu_species_table<-Menu_species_table%>%
  rename(Species = sci_name)


##Make a traits table for everything that could be ID'ed to species.----
Menu_traits_table<-left_join(Menu_species_table, fish_ecology_traits1, by = "Species")

Menu_traits_table<-left_join(Menu_traits_table, crustacean_ecology_traits1)

fish_common_names<-rfishbase::common_names(species_list = species_names$species)
#crustacean_common_names<-rfishbase::common_names(species_list = species_names$species, server = "sealifebase")

print(fish_common_names, n = 384)


taxa_matrix = as(otu_table(Menu_taxa_merge_species_named.pa.ps), "matrix")
taxa_names<-as.data.frame(tax_table(Menu_taxa_merge_species_named.pa.ps))%>%
  select("species")
taxon_table = as.data.frame(taxa_matrix)
rownames(taxon_table)<-taxa_names$species

taxa_names$species
nrow(taxa_names)
nrow(species_names)
fish_ecology_traits<-ecology(species_list = taxa_names$species, server = "fishbase")
fish_ecology_traits

crustacean_ecology_traits<-ecology(species_list = taxa_names$species, server = "sealifebase")
crustacean_ecology_traits


rfishbase::fb_tables()
#fish_food<-rfishbase::diet_items(species_list = species_names$species)
#get FoodI, FoodII, FoodIII, only available for fish, apparently
#fish_food<-fish_food%>%
  select(FoodI,FoodII,FoodIII)

#crust_food<-rfishbase::diet_items(species_list = species_names$species, server="sealifebase")

#it's not able to find the species of my species list.

fish_morphology<-rfishbase::morphology(species_list = species_names$species)

colnames(fish_morphology)

fish_morphology_select<-fish_morphology%>%
  select(Species, BodyShapeI, BodyShapeII, TypeofMouth, PosofMouth, Forklength)
fish_morphology_select


rfishbase::morphometrics(species_list = species_names$species)
rfishbase::predators(species_list = species_names$species)


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

##Make new columns in the sample data for blacktips detected and hammerheads detected with presence absence.----  

taxa_merge_sample_data<-cbind(data.frame(sample_data(Menu_taxa_merge_species_named.ps)))

taxa_merge_sample_data1<-taxa_merge_sample_data%>%
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


taxa_merge_sample_data1
Menu_taxa_merge_species_named1.ps<-Menu_taxa_merge_species_named.ps

sample_data(Menu_taxa_merge_species_named1.ps)<-taxa_merge_sample_data1

Menu_taxa_merge_species_named1.ps@sam_data$carcharhinus_detected


##add in fishbase traits to the taxonomy table so that the species can be merged/grouped by traits.  ----

#transformation options: https://rdrr.io/github/vmikk/metagMisc/man/phyloseq_standardize_otu_abundance.html

##do a new cca with the new variables----
#use the Hellinger transformation for the ASV counts (square-root of relative abundances)

Menu_taxa_merge1.hell.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named1.ps, method = "hellinger")

#Try the Wisconsin transformation: species first standardized by maxima then sites by site totals.
Menu_taxa_merge1.wisc.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named1.ps, method = "wisconsin")

Menu_taxa_merge_species_named1.pa.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named1.ps, method = "pa")

#get rid of samples with no temp data
Menu_taxa_merge_temp1.hell.ps<-ps_filter(Menu_taxa_merge1.hell.ps, Water_temp != "NA")


cca_blacktip_site<-ordinate(physeq = Menu_taxa_merge1.hell.ps, method = "CCA", formula = ~ carcharhinus_detected +blacktip_site)
summary(cca_blacktip_site)
anova(cca_blacktip_site, by ="margin")

cca_hammerhead_site<-ordinate(physeq = Menu_taxa_merge1.hell.ps, method = "CCA", formula = ~ hammerheads_detected +hammerhead_site)
summary(cca_hammerhead_site)
anova(cca_hammerhead_site, by ="margin")

cca_shark_site<-ordinate(physeq = Menu_taxa_merge1.hell.ps, method = "CCA", formula = ~ hammerhead_site+ blacktip_site + hammerheads_detected + carcharhinus_detected)
summary(cca_shark_site)
anova(cca_shark_site, by ="margin")


p_cca_shark_site<-plot(cca_shark_site, scaling = "symmetric", hill = TRUE, display= c("species","site","bp"))
identify(p_cca_shark_site, what = "species")
#identify(p_taxa_merge, what = "sites")



#indicator species analysis of blacktip samples and hammerhead samples (obviously ignoring those sharks)----
#https://emf-creaf.github.io/indicspecies/articles/IndicatorSpeciesAnalysis.html#additional-functions-to-estimate-and-test-the-association-between-species-and-groups-of-sites

library(indicspecies)

#let's work with the temp_bins
Menu_taxa_merge_species_named1.temp.ps <- ps_filter(Menu_taxa_merge_species_named1.ps, Water_temp != "NA")

#make matrix of species as columns and samples as rows
otu_df_indicspecies<-cbind(data.frame(otu_table(Menu_taxa_merge_species_named1.temp.ps)))
otu_df_indicspecies_trans<-as.data.frame(t(otu_df_indicspecies))

#create groupings


#make 6 bins for temps 

  Menu_taxa_merge_species_named1.temp.ps<- mutate_sample_data(Menu_taxa_merge_species_named1.temp.ps, temp_bin=cut(Water_temp, breaks=6))

add_temp_groupings<-cbind(data.frame(sample_data(Menu_taxa_merge_species_named1.temp.ps)))

add_temp_groupings<-add_temp_groupings%>%
  mutate(temp_group=
           if_else(temp_bin=="(22.2,23.5]",
                   1,
                   if_else(temp_bin=="(23.5,24.8]",
                           2,
                           if_else(temp_bin=="(24.8,26]",
                                   3,
                                   if_else(temp_bin=="(26,27.3]",
                                           4,
                                           if_else(temp_bin=="(27.3,28.5]",
                                                   5,
                                                   if_else(temp_bin=="(28.5,29.8]",
                                                           6,
                                                           NA)))))))

temp_groups<-add_temp_groupings$temp_group
indval <- multipatt(otu_df_indicspecies_trans, temp_groups, 
                    control = how(nperm=999)) 
#class(otu_df_indicspecies_trans)
#nrow(otu_df_indicspecies_trans)
#class(temp_groups)
#length(temp_groups)

summary(indval) 
#Thunnus genus significantly associated with temp_group 4 (temp bin (26,27.3]),p-value=0.042

summary(indval, indvalcomp=TRUE)
indval$sign

otu_df_indicspecies_trans.pa <- ifelse(otu_df_indicspecies_trans>0,1,0)
temp_phi <- multipatt(otu_df_indicspecies_trans.pa, temp_groups, func = "r.g", 
                 control = how(nperm=999)) 
summary(temp_phi)
# List of species associated to each combination: 
#   
#   Group 4  #sps.  1 
# stat p.value  
# Thunnus genus 1 0.78   0.013 *
#   
#   Group 3+5  #sps.  1 
# stat p.value    
# Prionurus genus 1 0.746   0.001 ***
#   
#   Group 1+2+4  #sps.  1 
# stat p.value  
# Sphoeroides annulatus 1 0.635    0.04 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#hammerhead sites:
hammerheads_detected<-add_temp_groupings$hammerheads_detected
hammerhead_phi <- multipatt(otu_df_indicspecies_trans.pa, hammerheads_detected, func = "r.g", 
                      control = how(nperm=999)) 
summary(hammerhead_phi)
# List of species associated to each combination: 
#   
#   Group 1  #sps.  2 
# stat p.value    
# Sphyrna lewini 1  1.000   0.001 ***
#   Gerres cinereus 1 0.658   0.034 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
hammerhead_site<-add_temp_groupings$hammerhead_site
hammerhead_site_phi <- multipatt(otu_df_indicspecies_trans.pa, hammerhead_site, func = "r.g", 
                            control = how(nperm=999)) 
summary(hammerhead_site_phi)

# Selected number of species: 16 
# Number of species associated to 1 group: 16 
# 
# List of species associated to each combination: 
#   
#   Group 0  #sps.  13 
# stat p.value   
# Daira americana 1                0.445   0.004 **
#   Diodon genus 1                   0.445   0.003 **
#   Synalpheus genus 1               0.427   0.002 **
#   Cronius ruber 1                  0.400   0.008 **
#   Prionurus genus 1                0.368   0.014 * 
#   Palaemon ritteri 1               0.346   0.013 * 
#   Scarus ghobban 1                 0.319   0.029 * 
#   Microcassiope genus 1            0.317   0.035 * 
#   Eriphides hispida 1              0.316   0.042 * 
#   Neogonodactylus bahiahondensis 1 0.316   0.043 * 
#   Grapsus grapsus 1                0.311   0.031 * 
#   Abudefduf genus 1                0.296   0.038 * 
#   Calcinus explorator 1            0.292   0.038 * 
#   
#   Group 1  #sps.  3 
# stat p.value  
# Callinectes arcuatus 1 0.302   0.033 *
#   Sphyrna lewini 1       0.302   0.028 *
#   Scomberomorus genus 1  0.302   0.027 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

blacktip_site<-add_temp_groupings$blacktip_site
blacktip_site_phi <- multipatt(otu_df_indicspecies_trans.pa, blacktip_site, func = "r.g", 
                                 control = how(nperm=999)) 
summary(blacktip_site_phi)

# List of species associated to each combination: 
#   
#   Group 0  #sps.  4 
# stat p.value   
# Myripristis genus 1   0.343   0.010 **
#   Calcinus explorator 1 0.295   0.040 * 
#   Palaemon ritteri 1    0.290   0.047 * 
#   Synalpheus genus 1    0.290   0.050 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

carcharhinus_detected<-add_temp_groupings$carcharhinus_detected
blacktips_phi <- multipatt(otu_df_indicspecies_trans.pa, carcharhinus_detected, func = "r.g", 
                               control = how(nperm=999)) 
summary(blacktips_phi)
# Total number of species: 121
# Selected number of species: 3 
# Number of species associated to 1 group: 3 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  3 
# stat p.value    
# Carcharhinus genus 1   1.000   0.001 ***
#   Balanus trigonus 1     0.468   0.021 *  
#   Callinectes arcuatus 1 0.437   0.039 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

round((blacktips_phi$str),3)
round((hammerhead_phi$str),3)

#try a gllvm----
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

#the traits were only available for species. The list was built off of this phyloseq object: Menu_species_merge.pa.ps

menu_gllvm_data$trait
Menu_traits_table
menu_gllvm_data$trait<-Menu_traits_table

menu_gllvm_data$abund
species_table
menu_gllvm_data$abund<-t(species_table)

menu_gllvm_data$x
menu_gllvm_env_vars<-cbind(data.frame(sample_data(Menu_species_merge.pa.ps)))
menu_gllvm_env_vars<-menu_gllvm_env_vars%>%
  select(Site_Name,Microhabitat,Tide)
menu_gllvm_env_vars

menu_gllvm_data$x<-menu_gllvm_env_vars

fitp <- gllvm(y = menu_gllvm_data$abund, family = poisson(), num.lv = 2)
fitnb <- gllvm(y = menu_gllvm_data$abund, family = "negative.binomial", num.lv = 2)
AIC(fitp)
## [1] 3343.571
AIC(fitnb)
## [1] 3817.39

par(mfrow = c(1,2))
plot(fitnb, which = 1:2)

fitLAp <- gllvm(y=menu_gllvm_data$abund, family = poisson(), method = "LA", num.lv = 2)
fitLAnb <- gllvm(y=menu_gllvm_data$abund, family = "negative.binomial", method = "LA", num.lv = 2)
fitLAzip <- gllvm(y=menu_gllvm_data$abund, family = "ZIP", method = "LA", num.lv = 2)
AIC(fitLAp)
AIC(fitLAnb)
AIC(fitLAzip)

#fitLAp has the lowest AIC at 3335.085

plot(fitLAp)
plot.new()
ordiplot(fitLAp, biplot = TRUE)
         abline(h = 0, v = 0, lty=2)
dev.off()


# Parameters:
coef(fitLAp)

getLV(fitLAp)

# Standard errors for parameters:
fitLAp$sd

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


X <- menu_gllvm_data$x[,1:2]
fitpx1 <- gllvm(menu_gllvm_data$abund, X, family = poisson(), num.lv = 1)
coefplot(fitpx1, mfrow = c(1,2), cex.ylab = 0.8)
ggsave(filename = "../11_Vegan/Menu_combined2/Microhabitat_gllvm_coeff.jpg")
coefplot(fitpx1, mfrow = c(2,4), cex.ylab = 0.8)
ggsave(filename = "../11_Vegan/Menu_combined2/Site_and_Microhabitat_gllvm_coeff.jpg")


##Include trait information in the gllvm ----

#https://jenniniku.github.io/gllvm/articles/vignette1.html#incorporating-functional-traits-into-fourth-corner-models
criteria <- NULL
for(i in 1:5){
  fiti <- gllvm(menu_gllvm_data$abund, menu_gllvm_data$x, family = poisson(), num.lv = i, sd.errors = FALSE,
                formula = ~ Site_Name + Microhabitat + Tide, seed = 1234)
  criteria[i] <- summary(fiti)$AICc
  names(criteria)[i] = i
}

criteria
#1 latent variable had the lowest AIC criteria number, so we'll use that.

fit_env <- gllvm(menu_gllvm_data$abund, menu_gllvm_data$x, family = poisson(), num.lv = 1,
                 formula = ~Site_Name+Tide+Microhabitat, seed = 1234)
dev.off()
coefplot(fit_env, cex.ylab = 0.7, mfrow=c(1,2))


# Fit GLLVM without environmental variables and 1 latent variable:
fit1lv <- gllvm(menu_gllvm_data$abund, family = poisson(), num.lv = 1, seed = 1234)

# Correlation matrix
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

#These plots are contradictory... clearly the environmental variables that I tested are having a big impact.  


###add traits finally----
#have to do this for only fish, because the crustaceans have too many NA values for traits.

menu_no_na_gllvm_data<-menu_gllvm_data

menu_no_na_gllvm_data$trait<-menu_no_na_gllvm_data$trait[1:34,]
menu_no_na_gllvm_data$trait<-menu_no_na_gllvm_data$trait[,-6]
menu_no_na_gllvm_data$trait<-menu_no_na_gllvm_data$trait[,-19]
menu_no_na_gllvm_data$trait

menu_no_na_gllvm_data$abund

#still not working with the NA's removed from trait.  Maybe I have to filter the abund as well to match the remaining species.

# fit_4th <- gllvm(menu_no_na_gllvm_data$abund, menu_no_na_gllvm_data$x, menu_no_na_gllvm_data$trait, family = poisson(), num.lv = 1, 
#                  formula = menu_no_na_gllvm_data$abund ~ (Site_Name + Microhabitat + Tide) +
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

#try EcoPhyloMapper----
#https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13914
#https://github.com/ptitle/epm/wiki

#install.packages("epm")
library(epm)
library(sf)



#inspect the squirrel data and see how the Menu data can be formatted properly

library(rgbif)

squirrelGenera <- c('Spermophilus', 'Citellus', 'Tamias', 'Sciurus', 'Glaucomys', 'Marmota', 'Cynomys')
spOccList <- vector('list', length(squirrelGenera))
names(spOccList) <- squirrelGenera

for (i in 1:length(squirrelGenera)) {
  
  message('\t', squirrelGenera[i])
  
  key <- name_suggest(q = squirrelGenera[i], rank='genus')$data$key[1]
  xx <- occ_data(taxonKey = key, hasCoordinate = TRUE, hasGeospatialIssue = FALSE, limit = 10000)
  spOccList[[i]] <- as.data.frame(xx$data)
}
sharedHeaders <- Reduce(intersect, lapply(spOccList, colnames))
spOccList <- lapply(spOccList, function(x) x[, sharedHeaders])
occTable <- do.call(rbind, spOccList)

nrow(occTable) # 40206

# a little cleanup: Let's strip out the genus and species from the scientific name field, and drop any that are genus only
newNames <- occTable$acceptedScientificName
newNames <- lapply(strsplit(newNames, '\\s+'), function(x) x[1:2])
newNames <- sapply(newNames, function(x) paste0(x, collapse = '_'))
nrow(occTable) == length(newNames)

occTable <- cbind(occTable, genSp = newNames)

occTable <- occTable[!grepl('^BOLD:|,', occTable$genSp),]

sort(unique(occTable$genSp))

# split into a list of species tables
spOccList2 <- split(occTable, occTable$genSp)

sort(sapply(spOccList2, nrow))

#make the sf points list
sfPtsList <- lapply(spOccList2, function(x) st_as_sf(x[, c('genSp', 'decimalLongitude', 'decimalLatitude', 'basisOfRecord', 'year')], coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326))

EAproj <- '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

spPtsListEA <- lapply(sfPtsList, function(x) sf::st_transform(x, crs =EAproj))

extentPoly <- "POLYGON ((-2906193 4690015, -3110223 4376122, -3377032 4172092, -3769399 3873894, -4083291 3340276, -4067597 2728184, -2702163 -1509370, 782048.5 -4036208, 1425529 -3973429, 1802200 -3785093, 2163177 -2639384, 3622779 -2027293, 2900826 5396274, 578018.1 5207938, -2906193 4690015))"

squirrelPtsEPM100 <- createEPMgrid(spPtsListEA, resolution = 100000, extent = extentPoly, cellType = 'hex')

plot(squirrelPtsEPM100)


