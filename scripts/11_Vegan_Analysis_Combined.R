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

#convert to eDNA index (wisconsin transformation)

Menu_taxa_merge_species_named.eDNAindex.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named.ps, method = "wisconsin")



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
  iDist <- phyloseq::distance(Menu_taxa_merge_species_named.eDNAindex.ps, method=i)
  # Calculate ordination
  iMDS  <- ordinate(Menu_taxa_merge_species_named.eDNAindex.ps, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(Menu_taxa_merge_species_named.eDNAindex.ps, iMDS, color="Site_Name", shape="Microhabitat")
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


##calculate jsd distances and put them in a symmetrical matrix (presence/absence)----
# pairwise_jsd_distances<-as.matrix(print(phyloseq::distance(Menu_taxa_merge.ps, method="jsd", binary=TRUE, type= "samples", upper=TRUE, diag=TRUE)))
# 
# 
# #Yay, the matrix is symmetrical!
# set.seed(200)
# jsd_permanova<- adonis2(pairwise_jsd_distances~Menu_taxa_merge.ps@sam_data$Site_Name)

#following https://github.com/pmartinezarbizu/pairwiseAdonis

jsd_dist_matrix <- phyloseq::distance(Menu_taxa_merge_species_named.eDNAindex.ps, method ="jsd")

set.seed(200)
vegan::adonis2(jsd_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge_species_named.eDNAindex.ps)$Site_Name*Menu_taxa_merge_species_named.eDNAindex.ps@sam_data$Microhabitat)

#Both significant jsd on eDNA index values

##Calculate unifrac distance matrix----
unifrac_dist_matrix <- phyloseq::distance(Menu_taxa_merge_species_named.eDNAindex.ps, method ="unifrac")

set.seed(200)
vegan::adonis2(unifrac_dist_matrix ~ phyloseq::sample_data(Menu_taxa_merge_species_named.eDNAindex.ps)$Site_Name*Menu_taxa_merge_species_named.eDNAindex.ps@sam_data$Microhabitat)

#not significant with unifrac though.

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
#identify(p_taxa_merge, what = "species")
#these are good ones to identify [1]   1   7  23  27  44  59  60  64 102 117 119
#ooh, especially these: [1]   7  10  30  49  59  60  64  90 102 112 117 119

#saved a screen shot of this as Best CCA of species-vs-temp_and_microhabitat.png in 11_Vegan folder. 

#how can I automate this part?


#make a prettier plot with ggvegan?
ggvegan::ordiggplot(cca_ps_taxa_merge)+
  ggvegan::geom_ordi_point(data = cca_ps_taxa_merge)+
  ggvegan::geom_ordi_arrow(score = "ChiSquare")

#Vegan: try simper to see which are the species most affecting the distributions but use the wisconsin transformed numbers----
simper_table_wisc<-t(as.data.frame(cbind(otu_table(Menu_taxa_merge.wisc.ps))))

simper(simper_table_wisc, group = Menu_taxa_merge.wisc.ps@sam_data$Microhabitat, permutations = (999))


#Vegan: Use wisconsin transformed data (eDNA index) for CCA tests----
#Try the Wisconsin transformation: species first standardized by maxima then sites by site totals.
Menu_taxa_merge.wisc.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named.ps, method = "wisconsin")

Menu_taxa_merge_species_named.pa.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named.ps, method = "pa")



#get rid of samples with no temp data
Menu_taxa_merge_temp.wisc.ps<-ps_filter(Menu_taxa_merge.wisc.ps, Water_temp != "NA")


cca_ps_site_name2<-ordinate(physeq = Menu_taxa_merge.wisc.ps, method = "CCA", formula = ~Site_Name)
summary(cca_ps_site_name2)
# Call:
#   cca(formula = OTU ~ Site_Name, data = data) 
# 
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          14.723     1.0000
# Constrained     1.702     0.1156
# Unconstrained  13.021     0.8844

anova(cca_ps_site_name2,permutations = 9999)
# Permutation test for cca under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: cca(formula = OTU ~ Site_Name, data = data)
# Df ChiSquare     F Pr(>F)  
# Model     6    1.7022 1.133 0.0516 .
# Residual 52   13.0211               
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

cca_ps_taxa_merge2 <- ordinate( physeq = Menu_taxa_merge_temp.wisc.ps,  method = "CCA", formula = ~ Microhabitat+Tide+Water_temp + Condition(Site_Name), scale=TRUE)

summary(cca_ps_taxa_merge2)
# Call:
#   cca(formula = OTU ~ Microhabitat + Tide + Water_temp + Condition(Site_Name),      data = data) 
# 
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          14.587    1.00000
# Conditioned     1.741    0.11935
# Constrained     1.295    0.08877
# Unconstrained  11.551    0.79188

anova(cca_ps_taxa_merge2, by = "margin")
# Permutation test for cca under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# Model: cca(formula = OTU ~ Microhabitat + Tide + Water_temp + Condition(Site_Name), data = data)
# Df ChiSquare      F Pr(>F)    
# Microhabitat  2    0.6338 1.2894  0.001 ***
#   Tide          1    0.2752 1.1198  0.140    
# Water_temp    1    0.3726 1.5161  0.001 ***
#   Residual     47   11.5515                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

cca_ps_taxa_merge2 <- ordinate( physeq = Menu_taxa_merge_temp.wisc.ps,  method = "CCA", formula = ~ Microhabitat+Water_temp + Condition(Site_Name+Tide), scale=TRUE)

summary(cca_ps_taxa_merge2)
# Call:
#   cca(formula = OTU ~ Microhabitat + Water_temp + Condition(Site_Name +      Tide), data = data) 
# 
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          14.587    1.00000
# Conditioned     2.025    0.13880
# Constrained     1.011    0.06932
# Unconstrained  11.551    0.79188

anova(cca_ps_taxa_merge2, by = "margin")

# Permutation test for cca under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# Model: cca(formula = OTU ~ Microhabitat + Water_temp + Condition(Site_Name + Tide), data = data)
# Df ChiSquare      F Pr(>F)    
# Microhabitat  2    0.6338 1.2894  0.001 ***
#   Water_temp    1    0.3726 1.5161  0.001 ***
#   Residual     47   11.5515                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


p_taxa_merge2<-plot(cca_ps_taxa_merge2, scaling = "site", hill = TRUE, display= c("species","bp"))
#identify(p_taxa_merge2, what = "species")


#and finally the wisconsin transformed (eDNA index) with environmental data added.

Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps<-readRDS("~/Desktop/Galapagos projects/Menu/10_Phyloseq/Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps.RDS")

cca_ps_sst_eDNAindex<-ordinate(physeq=Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, method = "CCA", formula = ~Site_Name+Microhabitat+Tide+scaled_sst)

summary(cca_ps_sst_eDNAindex)
# Call:
#   cca(formula = OTU ~ Site_Name + Microhabitat + Tide + scaled_sst,      data = data) 
# 
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          14.723     1.0000
# Constrained     3.003     0.2039
# Unconstrained  11.721     0.7961
anova(cca_ps_sst_eDNAindex, by="margin") 

# Model: cca(formula = OTU ~ Site_Name + Microhabitat + Tide + scaled_sst, data = data)
# Df ChiSquare      F Pr(>F)    
# Site_Name     6    1.6593 1.1326  0.041 *  
#   Microhabitat  2    0.6278 1.2856  0.001 ***
#   Tide          1    0.2766 1.1328  0.137    
# scaled_sst    1    0.3773 1.5452  0.001 ***
#   Residual     48   11.7205                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)

#Condition out Site effects and tide
cca_ps_sst_eDNAindex<-ordinate(physeq=Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps, method = "CCA", formula = ~Microhabitat+scaled_sst+scaled_oxygen+scaled_phosphate+Condition(Site_Name+Tide))

summary(cca_ps_sst_eDNAindex)
# Call:
#   cca(formula = OTU ~ Microhabitat + scaled_sst + scaled_nitrate + Condition(Site_Name + Tide), data = data) 
# 
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          14.723     1.0000
# Conditioned     1.984     0.1348
# Constrained     1.626     0.1105
# Unconstrained  11.113     0.7548

anova(cca_ps_sst_eDNAindex, by="margin") 

# Model: cca(formula = OTU ~ Microhabitat + scaled_sst + scaled_oxygen + scaled_phosphate + Condition(Site_Name + Tide), data = data)
# Df ChiSquare      F Pr(>F)    
# Microhabitat      2    0.5973 1.2363  0.006 ** 
#   scaled_sst        1    0.4265 1.7656  0.001 ***
#   scaled_oxygen     1    0.2146 0.8884  0.744    
# scaled_phosphate  1    0.2355 0.9748  0.558    
# Residual         46   11.1129                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



p_cca_ps_sst_eDNAindex<-plot(cca_ps_sst_eDNAindex, scaling = "site", hill = TRUE, display= c("species","bp"))
#identify(p_cca_ps_sst_eDNAindex, what = "species")

#Vegan: try simper to see which are the species most affecting the distributions but use the wisconsin transformed numbers----
simper_table_wisc1<-t(as.data.frame(cbind(otu_table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps))))

simper(simper_table_wisc1, group = Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$hammerheads_detected, permutations = (999))

simper_hhdetected<-simper(simper_table_wisc1, group = Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$hammerheads_detected, permutations = (999))


simper(simper_table_wisc1, group = Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$carcharhinus_detected, permutations = (999))

simper_carchdetected<-simper(simper_table_wisc1, group = Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$carcharhinus_detected, permutations = (999))




#check if all of this holds true with the pa transformed object----
#get rid of samples with no temp data
Menu_taxa_merge_temp.pa.ps<-ps_filter(Menu_taxa_merge_species_named.pa.ps, Water_temp != "NA")


cca_ps_site_name_pa<-ordinate(physeq = Menu_taxa_merge_species_named.pa.ps, method = "CCA", formula = ~Site_Name+Microhabitat+Tide) #These explain 22% of the variance which is 4.38
summary(cca_ps_site_name_pa)
anova(cca_ps_site_name_pa, by="margin") 

#Permutation test for cca under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# Model: cca(formula = OTU ~ Site_Name + Microhabitat + Tide, data = data)
# Df ChiSquare      F Pr(>F)    
# Site_Name     6    0.6463 1.5460  0.001 ***
#   Microhabitat  2    0.2187 1.5693  0.002 ** 
#   Tide          1    0.0949 1.3617  0.071 .  
# Residual     49    3.4138     


cca_ps_taxa_merge_pa <- ordinate( physeq = Menu_taxa_merge_temp.pa.ps,  method = "CCA", formula = ~ Microhabitat+Water_temp+Tide + Condition(Site_Name), scale=TRUE)

summary(cca_ps_taxa_merge_pa)

#
# Call:
#   cca(formula = OTU ~ Microhabitat + Water_temp + Tide + Condition(Site_Name),      data = data) 
# 
# Partitioning of scaled Chi-square:
#   Inertia Proportion
# Total          4.3401     1.0000
# Conditioned    0.6545     0.1508
# Constrained    0.4766     0.1098
# Unconstrained  3.2090     0.7394
anova(cca_ps_taxa_merge_pa, by = "margin")

# Model: cca(formula = OTU ~ Microhabitat + Water_temp + Tide + Condition(Site_Name), data = data)
# Df ChiSquare      F Pr(>F)    
# Microhabitat  2    0.2062 1.5097  0.003 ** 
#   Water_temp    1    0.1557 2.2811  0.001 ***
#   Tide          1    0.0997 1.4598  0.013 *  
#   Residual     47    3.2090                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

cca_ps_taxa_merge_pa1 <- ordinate( physeq = Menu_taxa_merge_temp.pa.ps,  method = "CCA", formula = ~ Water_temp +Tide +Condition(Site_Name+Microhabitat), scale=TRUE)

summary(cca_ps_taxa_merge_pa1)

anova(cca_ps_taxa_merge_pa1, by = "margin")

p_taxa_merge_pa<-plot(cca_ps_taxa_merge_pa1, scaling = "site", hill = TRUE, display= c("species","bp"))
#identify(p_taxa_merge_pa, what = "species")





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

#"The statistical parameter α quantifies the degree of preferential selection for balls colored B to occupy boxes already occupied by balls colored A."



# compute the affinity between elements in rows (= Taxa)----

taxa_affinity <- affinity(data = Shark_Menu_taxa_merge_species_named_traits.pa.ps@otu_table, row.or.col = "row", squarematrix = c("all"))
plotgg(data = taxa_affinity, variable = "alpha_mle", legendlimit = "datarange")+
  theme(axis.text.x= element_text(angle=90,vjust=0.5))


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


##Then the species----


species_matrix = as(otu_table(Shark_Menu_taxa_merge_species_named_traits.pa.ps), "matrix")
species_names<-as.data.frame(tax_table(Shark_Menu_taxa_merge_species_named_traits.pa.ps))%>%
  dplyr::select("species")
species_table = as.data.frame(species_matrix)
rownames(species_table)<-species_names$species


species_affinity<-affinity(data = species_table, row.or.col = "row", sigPval = "0.05",lev=0.95, pvalType="Blaker", squarematrix=c("alpha_mle_sig","alpha_mle","p_value"))



plotgg(data = species_affinity, variable = "alpha_mle", legendlimit = "datarange", text.size=.3)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  theme(axis.text=element_text(size=5))

species_affinity$alpha_mle_sig

species_aff_sig05<-subset(species_affinity$all, subset = p_value<0.05)

species_aff_sig05
species_affinity
nrow(species_aff_sig05)


#ggsave(paste0("../11_Vegan/",Primer,"/",Primer,"_species_coocurrence_affinity_plot.jpg"))

library(superheat)
species_aff_sig05_names1<-gsub(x = species_aff_sig05$entity_1, pattern = " 1",replacement = "")
species_aff_sig05_names2<-gsub(x = species_aff_sig05$entity_2, pattern = " 1",replacement = "")

species_aff_sig05_names<-species_aff_sig05%>%
  mutate(entity1=species_aff_sig05_names1,
         entity2=species_aff_sig05_names2)


species_aff_sig05_wide<-species_aff_sig05_names%>%
  pivot_wider(names_from = entity1,id_cols = entity2, values_from = alpha_mle)

class(species_aff_sig05_wide)
species_aff_sig05_wide1<-column_to_rownames(species_aff_sig05_wide, var = "entity2")
species_aff_sig05_wide2<-replace(species_aff_sig05_wide1,is.na(species_aff_sig05_wide1),0)


heatmap<-superheat(species_aff_sig05_wide1,
                   # normalize variables
                   scale = F,
                   # order rows/cols based on heirarchical clustering
                   pretty.order.rows = FALSE,
                   pretty.order.cols = FALSE,
                   # increase size of left labels
                   left.label.size = 0.18,
                   bottom.label.size = 0.2,
                   bottom.label.text.angle = 90,
                   left.label.text.size = 2.5,
                   bottom.label.text.size = 2.5,
                   force.grid.hline = TRUE,
                   force.grid.vline = TRUE,
                   extreme.values.na = TRUE,
                   heat.na.col = "grey",
                   order.cols=samplecols,
                   order.rows=samplerows)

#order rows and columns to make it easier to see the sharks
samplecols<-heatmap$order.cols
samplerows<-heatmap$order.rows
samplerows<-c(11,1:10,12:98)
samplecols<-c(51,1:50,52:87)

heatmap<-superheat(species_aff_sig05_wide1,
                   # normalize variables
                   scale = F,
                   # order rows/cols based on heirarchical clustering
                   pretty.order.rows = FALSE,
                   pretty.order.cols = FALSE,
                   # increase size of left labels
                   left.label.size = 0.18,
                   bottom.label.size = 0.2,
                   bottom.label.text.angle = 90,
                   left.label.text.size = 2.5,
                   bottom.label.text.size = 2.5,
                   force.grid.hline = TRUE,
                   force.grid.vline = TRUE,
                   extreme.values.na = TRUE,
                   heat.na.col = "grey",
                   order.cols=samplecols,
                   order.rows=samplerows)

plotgg(data = species_affinity, variable = "simpson", legendlimit = "datarange",text.size = 5)

plotgg(data = species_affinity, variable = "sorensen", legendlimit = "datarange",text.size = 5)

plotgg(data = species_affinity, variable = "jaccard", legendlimit = "datarange")

#plotgg(data = species_affinity, variable="p_value", legendlimit = "datarange")

#find which affinities have signficant p-values:
species_affinity_summary<-summary(species_affinity$all)

species_affinity_summary




#Create a site species table with the number of PCRs in which each species was present----
# ## Start with the original un-trimmed Menu_combined
# 
# 
# Menu_combined.ps<-readRDS(file= paste("../09_Metabar_to_Phyloseq/Menu_Combined_Phyloseq.rds", sep = ""))
# 
# Menu_combined.ps
# 
# #compute by site after binary transformation of the ASVs assigned to species level leaving NAs
# Menu_combined_taxa_merge.ps<-phyloseq::tax_glom(Menu_combined.ps, taxrank = "species", NArm = FALSE)
# #transform to presence absence
# Menu_combined_taxa_merge.pa.ps<-phyloseq_standardize_otu_abundance(Menu_combined_taxa_merge.ps, method = "pa")

#change from the above code using the pa transformed Shark Menu
#site affinity matrix----
#merge samples and sum
Shark_Menu_taxa_merge_species_named_traits_sitemerge.pa.ps<-merge_samples(Shark_Menu_taxa_merge_species_named_traits.pa.ps , group = "Site_Name",fun = "sum")
#the numbers of the otu table are now the number of PCRs that had that species or taxa present.
#make a matrix and dataframe of it.
taxa_site_matrix = as(otu_table(Shark_Menu_taxa_merge_species_named_traits_sitemerge.pa.ps), "matrix")
taxa_site_table = as.data.frame(taxa_site_matrix)

species_site_binary<-dataprep(data = taxa_site_matrix, row.or.col = "row", datatype = "abundance",
         threshold = 1, class0.rule = "less")

site_affinity<-affinity(data = species_site_binary, row.or.col = "col", squarematrix = "all",sigPval = "0.05")

plotgg(data = site_affinity, variable = "alpha_mle",legendlimit = "datarange")
#save this to 11_Vegan/vis as PDF

#species by site-wise affinity----

species_site_affinity<-affinity(data = species_site_binary, row.or.col = "row",sigPval = "0.05",squarematrix = "all")

plotgg(data = species_site_affinity, variable = "alpha_mle",legendlimit = "datarange")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  theme(axis.text=element_text(size=5))

#that's all of the taxa.  How can I only show the ones with significant p-values?

#
species_site_affinity1 <- affinity(data = species_site_binary, row.or.col = "row", lev=0.95, pvalType="Blaker")

plotgg(species_site_affinity1, variable = "alpha_mle",legendlimit = "datarange")
       
species_site_aff_sig05<-subset(species_site_affinity1$all, subset = p_value<0.05)

species_site_aff_sig05
nrow(species_site_aff_sig05)
write.csv(species_site_aff_sig05, file = "../11_Vegan/species_site_cooccurrence_affinity_sig05.csv")


#Maybe a better heatmap?----
#devtools::install_github("rlbarter/superheat")
library(superheat)


# superheat(mtcars,
#           # normalize variables
#           scale = T,
#           # order rows/cols based on heirarchical clustering
#           pretty.order.rows = TRUE,
#           pretty.order.cols = TRUE,
#           # plot miles per gallon next to the rows
#           yr = mtcars$mpg,
#           yr.axis.name = "miles per gallon",
#           # plot correlation with mpg above columns
#           yt = cor(mtcars)[, "mpg"],
#           yt.plot.type = "bar",
#           yt.axis.name = "correlation with mpg",
#           # increase size of left labels
#           left.label.size = 0.45)
# 
# 
# mtcars
# species_aff_sig05_select<-species_aff_sig05%>%
#   select(where(is.numeric))
# species_aff_sig05_select

species_aff_sig05_names3<-gsub(x = species_site_aff_sig05$entity_1, pattern = " 1",replacement = "")
species_aff_sig05_names4<-gsub(x = species_site_aff_sig05$entity_2, pattern = " 1",replacement = "")

species_site_aff_sig05_names<-species_site_aff_sig05%>%
  mutate(entity1=species_aff_sig05_names3,
         entity2=species_aff_sig05_names4)


species_site_aff_sig05_wide<-species_site_aff_sig05_names%>%
  pivot_wider(names_from = entity1,id_cols = entity2, values_from = alpha_mle)

class(species_site_aff_sig05_wide)
species_site_aff_sig05_wide1<-column_to_rownames(species_site_aff_sig05_wide, var = "entity2")
species_site_aff_sig05_wide2<-replace(species_site_aff_sig05_wide1,is.na(species_site_aff_sig05_wide1),0)


heatmap<-superheat(species_site_aff_sig05_wide2,
          # normalize variables
          scale = F,
          # order rows/cols based on heirarchical clustering
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          # increase size of left labels
          #left.label.size = 0.45,
          bottom.label.text.angle = 90,
          left.label.size = 0.2,
          bottom.label.size = 0.41,
          left.label.text.size = 3,
          bottom.label.text.size = 3,
          force.grid.hline = TRUE,
          force.grid.vline = TRUE, 
          heat.na.col = "grey")

#first make it with pretty order rows and cols with species_site_aff_sig05_wide2, then run this code and re-run the above with species_site_aff_sig05_wide1
sitecols<-heatmap$order.cols
siterows<-heatmap$order.rows

heatmap<-superheat(species_site_aff_sig05_wide1,
                   # normalize variables
                   scale = F,
                   # order rows/cols based on heirarchical clustering
                   pretty.order.rows = FALSE,
                   pretty.order.cols = FALSE,
                   # increase size of left labels
                   #left.label.size = 0.45,
                   bottom.label.text.angle = 90,
                   left.label.size = 0.2,
                   bottom.label.size = 0.41,
                   left.label.text.size = 3,
                   bottom.label.text.size = 3,
                   force.grid.hline = TRUE,
                   force.grid.vline = TRUE, 
                   heat.na.col = "grey",
                   order.cols = sitecols,
                   order.rows = siterows)

#Cooccur package----

#make a species x site matrix
#merge samples by Site_Name

Menu_taxa_merge_site_merge.pa.ps<-merge_samples(Menu_taxa_merge.pa.ps,group = "Site_Name")
#Menu_taxa_merge_site_merge.pa.ps<-merge_samples(Menu_taxa_merge_species_named.pa.ps,group = "Site_Name")

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
ggsave("../11_Vegan/Menu_combined2/species_cooccurrence_matrix.jpg")

hammerhead_cooccurring_sp<-pair(mod = species_cooccur, spp = "Sphyrna lewini")
write.csv(hammerhead_cooccurring_sp, file = "../hammerhead_cooccurr_species_sitemerge.csv")

carcharhinus_coocurring_sp<-pair(mod=species_cooccur, spp="Carcharhinus genus")

pair.profile(mod = species_cooccur)
ggsave("../11_Vegan/pos_neg_rand_cooccur_bysite.jpg")



#######do the same with Final phyloseq object (taxa-merged_species-named)----

#Take the phyloseq tables into dataframes and visualize

species_by_sample_mat<-cbind(as.matrix(otu_table(Menu_taxa_merge_species_named.pa.ps)))

species_by_sample_names<-rownames(species_by_sample_mat)
species_by_sample_names<-gsub(pattern = " 1",replacement = "",x = species_by_sample_names)
rownames(species_by_sample_mat)<-species_by_sample_names


species_cooccur<-cooccur(species_by_sample_mat, spp_names = TRUE, thresh = TRUE, true_rand_classifier = 0.1, prob = "hyper")

class(species_cooccur)
summary(species_cooccur)

species_cooccur_prob_table<-prob.table(species_cooccur)

plot(species_cooccur) # add "plotrand = TRUE" to include completely random species
ggsave("../11_Vegan/Menu_combined2/species_cooccurrence_matrix.jpg")

hammerhead_cooccurring_sp<-pair(mod = species_cooccur, spp = "Sphyrna lewini")

carcharhinus_coocurring_sp<-pair(mod=species_cooccur, spp="Carcharhinus genus")

Mugil_thoburni_cooccurring_sp<-pair(mod = species_cooccur, spp="Mugil thoburni")

Gerres_cinereus_coocurring_sp<-pair(mod = species_cooccur, spp="Gerres cinereus")

significant_pair_coocurr<-as.data.frame(print(species_cooccur))
pair.profile(mod = species_cooccur)
ggsave("../11_Vegan/significant_pairings_cooccur.jpg")

#Then cooccur taxa merged species named by site----

Menu_taxa_merge_species_named_site_merge.pa.ps<-merge_samples(Menu_taxa_merge_species_named.pa.ps,group = "Site_Name")

species_by_site_mat<-t(cbind(as.matrix(otu_table(Menu_taxa_merge_species_named_site_merge.pa.ps))))

species_by_site_names<-rownames(species_by_site_mat)
species_by_site_names<-gsub(pattern = " 1",replacement = "",x = species_by_site_names)
rownames(species_by_site_mat)<-species_by_site_names


species_cooccur_site<-cooccur(species_by_site_mat, spp_names = TRUE, thresh = TRUE, true_rand_classifier = 0.1, prob = "hyper")

class(species_cooccur_site)
summary(species_cooccur_site)

species_cooccur_prob_table<-prob.table(species_cooccur_site)

plot(species_cooccur_site) # add "plotrand = TRUE" to include completely random species
ggsave("../11_Vegan/Menu_combined2/species_cooccurrence_matrix_by_site.jpg")

hammerhead_cooccurring_sp_site<-pair(mod = species_cooccur_site, spp = "Sphyrna lewini")

carcharhinus_coocurring_sp_site<-pair(mod = species_cooccur_site, spp = "Carcharhinus genus")

Mugil_thoburni_cooccurring_sp<-pair(mod = species_cooccur_site, spp="Mugil thoburni")

Gerres_cinereus_coocurring_sp<-pair(mod = species_cooccur_site, spp="Gerres cinereus")

significant_pair_coocurr_site<-as.data.frame(print(species_cooccur_site))
pair.profile(mod = species_cooccur_site)
ggsave("../11_Vegan/significant_pairings_cooccur_by_site.jpg")

#Then cooccur taxa merged species named by microhabitat----

Menu_taxa_merge_species_named_microhabitat_merge.pa.ps<-merge_samples(Menu_taxa_merge_species_named.pa.ps,group = "Microhabitat")

species_by_microhabitat_mat<-t(cbind(as.matrix(otu_table(Menu_taxa_merge_species_named_microhabitat_merge.pa.ps))))

species_by_microhabitat_names<-rownames(species_by_microhabitat_mat)
species_by_microhabitat_names<-gsub(pattern = " 1",replacement = "",x = species_by_microhabitat_names)
rownames(species_by_microhabitat_mat)<-species_by_microhabitat_names


species_cooccur_microhabitat<-cooccur(species_by_microhabitat_mat, spp_names = TRUE, thresh = TRUE, true_rand_classifier = 0.1, prob = "hyper")

class(species_cooccur_microhabitat)
summary(species_cooccur_microhabitat)

species_cooccur_prob_table<-prob.table(species_cooccur_microhabitat)

plot(species_cooccur_microhabitat) # add "plotrand = TRUE" to include completely random species
ggsave("../11_Vegan/Menu_combined2/species_cooccurrence_matrix_by_microhabitat.jpg")

hammerhead_cooccurring_sp<-pair(mod = species_cooccur_microhabitat, spp = "Sphyrna lewini")

Mugil_thoburni_cooccurring_sp<-pair(mod = species_cooccur_microhabitat, spp="Mugil thoburni")

Gerres_cinereus_coocurring_sp<-pair(mod = species_cooccur_microhabitat, spp="Gerres cinereus")

significant_pair_coocurr_microhabitat<-as.data.frame(print(species_cooccur_microhabitat))
pair.profile(mod = species_cooccur_microhabitat)
ggsave("../11_Vegan/significant_pairings_cooccur_by_microhabitat.jpg")




#iNext Hill Numbers----
#This takes a while



otu_df_inext<-cbind(data.frame(otu_table(Menu1.ps)))
taxa_df_iNext<-cbind(data.frame(otu_table(Menu_taxa_merge_species_named.ps)))

#This takes a while
iNext_Menu0<-iNEXT(otu_df_inext, q=0, datatype="abundance", size=NULL, knots=40, se=TRUE, conf=0.95, nboot=100, endpoint=100000)

iNext_Menu1<-iNEXT(otu_df_inext, q=1, datatype="abundance", size=NULL, knots=40, se=TRUE, conf=0.95, nboot=100, endpoint=100000)

iNext_Menu2<-iNEXT(otu_df_inext, q=2, datatype="abundance", size=NULL, knots=40, se=TRUE, conf=0.95, nboot=100, endpoint=100000)

#Sample-size-based R/E curve = type=1
ggiNEXT(iNext_Menu, type = 1, se = TRUE)+
  title("ASV Diversity Rarefaction and Extrapolation Curves by Sample")+
  xlab("Assigned ASV Reads")+
  ylab("ASV Diversity Hill q=2")+
  theme(legend.position = 'none')

#ggsave(filename = "../11_Vegan/Menu_combined2/iNext_ASV_richness_v_coverage_rarefaction_curve.jpg")

taxa_df_iNext<-cbind(data.frame(otu_table(Menu_taxa_merge_species_named.ps)))
taxa_df_iNext_trans<-t(taxa_df_iNext)

#create it once, it takes a long time.
iNext_menu_taxa_by_sample<-iNEXT(taxa_df_iNext, q=c(0,1,2), datatype = "abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)

#Sample-size-based R/E curve = type=1
#sample completeness curve= type=2
ggiNEXT(iNext_menu_taxa_by_sample, type = 1, se = TRUE,facet.var = "Order.q")+
  xlab("Assigned ASV Reads")+
  ylab("Diversity Measures of Samples")+
  xlim(0,15000)+
  theme(legend.position = 'none')+
  title("Rarefaction and Extrapolation Curves of Taxa per Sample by ASV Read Counts")

ggsave("../11_Vegan/Menu_combined2/iNext_taxonomic_richness_rarefaction_and_extrapolation_curves.jpg")


Menu_taxa_merge_species_named_site_sums.ps<-merge_samples(Menu_taxa_merge_species_named.ps, group = "Site_Name",fun = "sum")

otu_df_inext_by_site<-cbind(data.frame(otu_table(Menu_taxa_merge_species_named_site_sums.ps)))
otu_df_inext_by_site_trans<-t(otu_df_inext_by_site)

#create it once, it takes a long time.
iNext_Menu_by_site<-iNEXT(otu_df_inext_by_site_trans, q=c(0,1,2), datatype = "abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)

#Sample-size-based R/E curve = type=1
#sample completeness curve= type=2
ggiNEXT(iNext_Menu_by_site, type = 1, se = TRUE,facet.var = "Order.q")+
  xlab("Assigned ASV Reads")+
  ylab("Diversity Measures of Sites")+
  xlim(0,30000)+
  scale_color_viridis_d(option = "A")
ggsave("../11_Vegan/Menu_combined2/iNext_taxonomic_diversity_by_read_coverage_per_site.jpg")

Menu_taxa_merge_species_named_microhabitat_sums.ps<-merge_samples(Menu_taxa_merge_species_named.ps, group = "Microhabitat",fun = "sum")

otu_df_inext_by_microhabitat<-cbind(data.frame(otu_table(Menu_taxa_merge_species_named_microhabitat_sums.ps)))
otu_df_inext_by_microhabitat_trans<-t(otu_df_inext_by_microhabitat)

#create it once, it takes a long time.
iNext_Menu_by_microhabitat<-iNEXT(otu_df_inext_by_microhabitat_trans, q=c(0,1,2), datatype = "abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)

#Sample-size-based R/E curve = type=1
#sample completeness curve= type=2
ggiNEXT(iNext_Menu_by_microhabitat, type = 1, se = TRUE,facet.var = "Order.q", color.var = "Assemblage")+
  title("Rarefaction Curves of Diversity Measures of Microhabitats")+
  xlab("Assigned ASV Reads")+
  ylab("Diversity Measures of Microhabitats")+
  xlim(0,30000)+
  scale_color_viridis_d(option = "D")
ggsave("../11_Vegan/Menu_combined2/iNext_taxonomic_diversity_by read_coverage_per_microhabitat.jpg")



Menu_taxa_merge_species_named_hhsite_sums.ps<-merge_samples(Shark_Menu_taxa_merge_species_named_traits.ps, group = "hammerhead_site",fun = "sum")

otu_df_inext_by_hhstatus<-cbind(data.frame(otu_table(Menu_taxa_merge_species_named_hhsite_sums.ps)))
otu_df_inext_by_hhstatus_trans<-t(otu_df_inext_by_hhstatus)

#this part takes a long time
iNext_Menu_by_hhstatus<-iNEXT(otu_df_inext_by_hhstatus_trans, q=c(0,1,2), datatype = "abundance", size=NULL, endpoint=NULL, knots=40, se=TRUE, conf=0.95, nboot=50)

plot.new()
#Sample-size-based R/E curve = type=1
#sample completeness curve= type=2
ggiNEXT(iNext_Menu_by_hhstatus, type = 1, se = TRUE,facet.var = "Assemblage", color.var = "Assemblage")+
  title("Rarefaction Curves of Diversity Measures of Hammerhead Sites and Non-Hammerhead Sites")+
  xlab("Assigned ASV Reads")+
  ylab("Diversity Measures of Hammerhead Sites and Non-Hammerhead Sites")+
  xlim(0,30000)
ggsave("../11_Vegan/Menu_combined2/iNext_taxonomic_diversity_by read_coverage_by_hhstatus.jpg")



#iNextPD visualization of trees
## install iNEXT package 
## install the latest version from github
library(devtools)
install_github('JohnsonHsieh/iNextPD')

## import packages
library(iNextPD)
library(ggplot2)
library(ade4)

data(bird)
str(bird)

bird.lab <- rownames(bird$abun)
bird.phy <- ade4::newick2phylog(bird$tre)
plot(bird.phy)
table.phylog(bird$abun, bird.phy, csize=4, f.phylog=0.7)

Menu_taxa_merge_species_named_eDNAindex_site_sums.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named_site_sums.ps, method = "wisconsin")

otu_df_inextpd_by_site<-cbind(data.frame(otu_table(Menu_taxa_merge_species_named_eDNAindex_site_sums.ps)))
otu_df_inextpd_by_site_trans<-t(otu_df_inextpd_by_site)

#gsub out the .1 with nothing and then the . with space in the colnames

names<-rownames(otu_df_inextpd_by_site_trans)
names1<-gsub(".1","",names)
names2<-gsub("."," ",names1, fixed = TRUE)

#re-name to match tree
names3<-rownames(otu_df_inextpd_by_site_trans)
names4<-gsub(".","_",names3, fixed = TRUE)


menu.lab<-names2
menu.abun<-otu_df_inextpd_by_site_trans
rownames(menu.abun)<-names4

menu.phy1<- (Menu_taxa_merge_species_named_eDNAindex_site_sums.ps@phy_tree)
class(menu.phy1)
# convert the tree to a character string
library(ape)
tree.character <- paste(write.tree(menu.phy1), collapse = "\n")

# change to phylog
menu.phy <- ade4::newick2phylog(tree.character)
#plot 
plot(menu.phy)
table.phylog(menu.abun, menu.phy, csize=1, f.phylog=0.7,labels.row = menu.lab,clabel.row = .6,clabel.col = 1,clegend = 0)


#think about doing this with site/microhabitat as a grouping variable instead of just site



#can't do more with iNextPD because of rcpp errors.

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
#  select(FoodI,FoodII,FoodIII)

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



#indicator species analysis of temperature groupings 
#https://emf-creaf.github.io/indicspecies/articles/IndicatorSpeciesAnalysis.html#additional-functions-to-estimate-and-test-the-association-between-species-and-groups-of-sites

library(indicspecies)

#let's work with the temp_bins
Menu_taxa_merge_species_named1.temp.ps <- ps_filter(Menu_taxa_merge_species_named1.ps, Water_temp != "NA")

#eDNA index
Menu_taxa_merge_species_named1.temp.ps<-phyloseq_standardize_otu_abundance(Menu_taxa_merge_species_named1.temp.ps, method="wisconsin")

#make matrix of species as columns and samples as rows
otu_df_indicspecies<-cbind(data.frame(otu_table(Menu_taxa_merge_species_named1.temp.ps)))
otu_df_indicspecies_trans<-as.data.frame(t(otu_df_indicspecies))

#create groupings


#make 6 bins for temps 

  Menu_taxa_merge_species_named1.temp.ps<- mutate_sample_data(Menu_taxa_merge_species_named1.temp.ps, temp_bin=cut(Water_temp, breaks=6))

add_temp_groupings<-cbind(data.frame(sample_data(Menu_taxa_merge_species_named1.temp.ps)))

add_temp_groupings<-add_temp_groupings%>%
  mutate(temp_group=
           if_else(temp_bin=="(21.8,23.2]",
                   1,
                   if_else(temp_bin=="(23.2,24.5]",
                           2,
                           if_else(temp_bin=="(24.5,25.8]",
                                   3,
                                   if_else(temp_bin=="(25.8,27.1]",
                                           4,
                                           if_else(temp_bin=="(27.1,28.5]",
                                                   5,
                                                   if_else(temp_bin=="(28.5,29.8]",
                                                           6,
                                                           NA)))))))

temp_groups<-add_temp_groupings$temp_group
indval <- multipatt(otu_df_indicspecies_trans, temp_groups,func = "IndVal.g", 
                    control = how(nperm=999)) 
#class(otu_df_indicspecies_trans)
#nrow(otu_df_indicspecies_trans)
#class(temp_groups)
#length(temp_groups)

summary(indval) 



summary(indval, indvalcomp=TRUE)
indval$sign

#indicator species analysis----
#start from eDNA index traits phyloseq object
#make matrix of species as columns and samples as rows
otu_df_indicspecies<-cbind(data.frame(otu_table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))
otu_df_indicspecies_trans<-as.data.frame(t(otu_df_indicspecies))



#fix the names of the taxa
indicspecies_names<-colnames(otu_df_indicspecies_trans)
indicspecies_names1<-gsub(pattern = " 1",replacement = "",indicspecies_names)
colnames(otu_df_indicspecies_trans)<-indicspecies_names1


otu_df_indicspecies_trans.pa <- ifelse(otu_df_indicspecies_trans>0,1,0)
# temp_phi <- multipatt(otu_df_indicspecies_trans.pa, temp_groups, func = "IndVal.g", 
#                  control = how(nperm=999)) 
# summary(temp_phi)
sample_indic<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))


#hammerhead sites:
hammerheads_detected<-sample_indic$hammerheads_detected
hammerhead_phi <- multipatt(otu_df_indicspecies_trans, hammerheads_detected, func = "IndVal.g", 
                      control = how(nperm=999)) 
summary(hammerhead_phi)


# Component ‘A’ is sample estimate of the probability that the surveyed site belongs to the target site group given the fact that the species has been found. This conditional probability is called the specificity or positive predictive value of the species as indicator of the site group.
# Component ‘B’ is sample estimate of the probability of finding the species in sites belonging to the site group. This second conditional probability is called the fidelity or sensitivity of the species as indicator of the target site group.
# 
# To display the indicator value components ‘A’ and ‘B’ one simply uses:
#   
#   summary(indval, indvalcomp=TRUE)

summary(hammerhead_phi, indvalcomp = TRUE)

# Multilevel pattern analysis
# ---------------------------
#   
#   Association function: IndVal.g
# Significance level (alpha): 0.05
# 
# Total number of species: 121
# Selected number of species: 4 
# Number of species associated to 1 group: 4 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  4 
# stat p.value    
# Sphyrna lewini         1.000   0.001 ***
#   Sphoeroides annulatus  0.873   0.016 *  
#   Caranx hippos          0.703   0.018 *  
#   Selar crumenophthalmus 0.684   0.023 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 





hammerhead_site<-sample_indic$hammerhead_site
hammerhead_site_phi <- multipatt(otu_df_indicspecies_trans, hammerhead_site, func = "IndVal.g", 
                            control = how(nperm=999)) 
summary(hammerhead_site_phi)

# Multilevel pattern analysis
# ---------------------------
#   
#   Association function: IndVal.g
# Significance level (alpha): 0.05
# 
# Total number of species: 121
# Selected number of species: 19 
# Number of species associated to 1 group: 19 
# 
# List of species associated to each combination: 
#   
#   Group 0  #sps.  13 
# stat p.value    
# Teleophrys cristulipes 0.824   0.002 ** 
#   Calcinus explorator    0.793   0.005 ** 
#   Palaemon ritteri       0.760   0.002 ** 
#   Prionurus genus        0.707   0.003 ** 
#   Daira americana        0.681   0.001 ***
#   Synalpheus genus       0.623   0.001 ***
#   Diodon genus           0.613   0.011 *  
#   Grapsus grapsus        0.610   0.019 *  
#   Cronius ruber          0.603   0.003 ** 
#   Microcassiope genus    0.515   0.016 *  
#   Thalassoma genus       0.493   0.035 *  
#   Arothron hispidus      0.473   0.035 *  
#   Eriphides hispida      0.426   0.039 *  
#   
#   Group 1  #sps.  6 
# stat p.value  
# Thunnus genus          0.515   0.023 *
#   Selar crumenophthalmus 0.465   0.031 *
#   Hippa pacifica         0.444   0.019 *
#   Callinectes arcuatus   0.400   0.026 *
#   Sphyrna lewini         0.400   0.027 *
#   Scomberomorus genus    0.400   0.031 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

blacktip_site<-sample_indic$blacktip_site
blacktip_site_phi <- multipatt(otu_df_indicspecies_trans, blacktip_site, func = "IndVal.g", 
                                 control = how(nperm=999)) 
summary(blacktip_site_phi)


# Multilevel pattern analysis
# ---------------------------
#   
#   Association function: IndVal.g
# Significance level (alpha): 0.05
# 
# Total number of species: 121
# Selected number of species: 6 
# Number of species associated to 1 group: 6 
# 
# List of species associated to each combination: 
#   
#   Group non-carch_site  #sps.  6 
# stat p.value   
# Prionurus genus        0.732   0.004 **
#   Halichoeres dispilus   0.591   0.049 * 
#   Diodon genus           0.582   0.025 * 
#   Synalpheus genus       0.552   0.037 * 
#   Myripristis genus      0.459   0.013 * 
#   Petrolisthes glasselli 0.429   0.026 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

carcharhinus_detected<-sample_indic$carcharhinus_detected
blacktips_phi <- multipatt(otu_df_indicspecies_trans, carcharhinus_detected, func = "IndVal.g", 
                               control = how(nperm=999)) 
summary(blacktips_phi)
# Multilevel pattern analysis
# ---------------------------
#   
#   Association function: IndVal.g
# Significance level (alpha): 0.05
# 
# Total number of species: 121
# Selected number of species: 5 
# Number of species associated to 1 group: 5 
# 
# List of species associated to each combination: 
#   
#   Group 1  #sps.  5 
# stat p.value    
# Carcharhinus genus     1.000   0.001 ***
#   Scarus ghobban         0.710   0.046 *  
#   Callinectes arcuatus   0.577   0.029 *  
#   Balanus trigonus       0.537   0.013 *  
#   Petrolisthes edwardsii 0.501   0.049 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
round((blacktips_phi$str),3)
round((hammerhead_phi$str),3)

#tide incoming or outgoing
tide_incoming<-sample_indic$Tide
tides_phi <- multipatt(otu_df_indicspecies_trans, tide_incoming, func = "IndVal.g", 
                           control = how(nperm=999)) 
summary(tides_phi)

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
coefplot(fitpx1, mfrow = c(2,2), cex.ylab = 0.8)

coefplot(fitpx1, mfrow = c(1,4), cex.ylab = 0.8)


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


#make a melted phyloseq object with psmelt to get it into the format necessary.  Use taxa merged species named object since we want to map taxa, not OTUs.


melted_Shark_Menu_taxa_merge_species_named_traits.eDNAindex = psmelt(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)

menu_occ<-melted_Shark_Menu_taxa_merge_species_named_traits.eDNAindex%>%
  select(-c(Associated_Blood_Samples,Associated_Fecal_Samples))%>%
  filter(Abundance>0)
menu_occ$OTU<-gsub(" 1","",x = menu_occ$OTU)

# split into a list of species tables
spOccList2 <- split(occTable, occTable$genSp)

menu_occ2<-split(menu_occ, menu_occ$OTU)

sort(sapply(spOccList2, nrow))

sort(sapply(menu_occ2,nrow))

#make the sf points list
sfPtsList <- lapply(spOccList2, function(x) st_as_sf(x[, c('genSp', 'decimalLongitude', 'decimalLatitude', 'basisOfRecord', 'year')], coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326))

sfPtsList <- lapply(menu_occ2, function(x) st_as_sf(x[, c('OTU', 'DD_long', 'DD_lat', 'DemersPelag','FoodTroph', 'Microhabitat')], coords = c('DD_long', 'DD_lat'), crs = 4326))

EAproj <- '+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

EAproj<-'EPSG:31986'

spPtsListEA <- lapply(sfPtsList, function(x) sf::st_transform(x, crs =EAproj))

extentPoly <- "POLYGON ((-2906193 4690015, -3110223 4376122, -3377032 4172092, -3769399 3873894, -4083291 3340276, -4067597 2728184, -2702163 -1509370, 782048.5 -4036208, 1425529 -3973429, 1802200 -3785093, 2163177 -2639384, 3622779 -2027293, 2900826 5396274, 578018.1 5207938, -2906193 4690015))"

interactiveExtent(spPtsListEA, cellType = "hex", bb = NULL)

extentPoly<- "POLYGON ((-615271.7 -140697, -615271.7 -141212.4, -615787 5659.045, -418927.8 12873.78, -419443.1 -153580.5, -614756.4 -141212.4, -615271.7 -140697))"

MenuPtsEPM <- createEPMgrid(spPtsListEA, resolution = 100000, extent = extentPoly, cellType = 'hex')

plot(MenuPtsEPM)

#try hilldiv ----

#https://www.biorxiv.org/content/10.1101/545665v1.full.pdf

#for phylogenetic (and non-phylogenetic) diversity based on Hill numbers.

# Compute neutral or phylogenetic Hill numbers from a single sample (vector) or count table (matrix). Hill numbers or numbers equivalents of diversity indices are diversity measures that compute diversity in effective number of OTUs, i.e. the number of equally abundant OTUs that would be needed to give the same value of diversity.

#The Hill number of q=0 yields a richness value, the Hill number of q=1 is the exponential of the Shannon index, and the Hill number of q=2 is the multiplicative inverse of the Simpson index ​(Jost, 2006)​.

#Hill numbers can also be computed while taking into account the phylogenetic or functional relationships among OTUs. When these so-called ‘phylogenetic Hill numbers’ are computed, the diversity is measured in effective number of equally abundant and equally distinct lineages ​(Chao et al., 2010)​. For two systems with identical number of types and relative abundances, the one with the largest phylogenetic differences across OTUs will be the one with the highest phylogenetic diversity or phylodiversity. Similar to neutral Hill numbers, phylogenetic Hill numbers are also closely related to popular phylogenetic diversity indices: Faith’s PD (when q=0), Allen’s H (when q=1) and Rao’s Q (when q=2) ​(Rao, 1982; Faith, 1992; Allen, Kon, & Bar-Yam, 2009)​.

#although the Hill numbers framework was originally developed for abundance data, it can also be applied to incidence data, a type of information broadly employed when dealing with, for example, ecological niche-related issues. In abundance-based approaches the DNA sequence is the unit upon which diversity is computed, while in incidence-based approaches the sample is the unit upon which diversity is measured ​(Chao et al., 2014)​. When computing Hill numbers from incidence data, it is important to note that the interpretation of both the measure and the measurement unit is slightly different to that of abundance data. Abundance-based Hill numbers measure the effective number of equally abundant OTUs in the system, while incidence-based Hill numbers measure the effective number of equally frequent (across samples) OTUs in the system ​(Chao et al., 2014)​.

#install.packages("hilldiv")
library(hilldiv)
##depth cov
data(bat.diet.otutable)
Menu1.ps<- ps_filter(Menu.ps, Microhabitat != "passive mid bay")

otu_df_hilldiv<-cbind(data.frame(otu_table(Menu1.ps)))
Menu_eDNAindex.ps<-phyloseq_standardize_otu_abundance(Menu1.ps, method = "wisconsin")
otu_df_eDNAindex_hilldiv<-cbind(data.frame(otu_table(Menu_eDNAindex.ps)))



depth_cov(otu_df_hilldiv,0)#richness (presence/absence)
depth_cov(otu_df_hilldiv,qvalue=1)# exponential of the Shannon index
depth_cov(otu_df_hilldiv,qvalue=2)#multiplicative inverse of the Simpson index

otu_df_hilldiv_80<-depth_filt(otu_df_hilldiv,80)
##hill div----

data(bat.diet.otutable)
data(bat.diet.tree)
str(bat.diet.tree)
str(Menu.ps@phy_tree)

#Make ultrametric tree for full Menu.ps object 
menu_tree_hilldiv<-Menu1.ps@phy_tree


menu_tree_hilldiv_ultra<-chronos(menu_tree_hilldiv,lambda = 1,model = "correlated")
class(menu_tree_hilldiv_ultra)
write.tree(menu_tree_hilldiv_ultra, 'my_ultra.tre')
menu_tree_hilldiv_ultra <- ape::read.tree('my_ultra.tre')
class(menu_tree_hilldiv_ultra)


data(bat.diet.hierarchy)

#One sample
bat.diet.sample <- bat.diet.otutable[,1]
hill_div(bat.diet.sample,0)
hill_div(bat.diet.sample,qvalue=1)

#One sample (phylogenetic)
names(bat.diet.sample) <- rownames(bat.diet.otutable)
hill_div(bat.diet.sample,1,bat.diet.tree)


#Multiple samples
hill_div(bat.diet.otutable,0)

###Menu: Get Hill diversity q=1 considering phylogenetic diversity for each sample----
# #Allen’s H using eDNAindex transformed read abundances----
# 
# AllensH_otu_sample_eDNAindex<-hill_div(otu_df_eDNAindex_hilldiv,1,menu_tree_hilldiv_ultra)
# AllensH_otu_sample_eDNAindex
# class(AllensH_otu_sample_eDNAindex)
# 
# Allens.df<-rownames_to_column(as.data.frame(AllensH_otu_sample_eDNAindex),var="Sample")
# 
# 
# ###Menu: Get Hill diversity q=0 considering phylogenetic diversity for each sample----
# #Faith's PD (presence/absence)----
# 
# FaithsPD_otu_sample<-hill_div(otu_df_hilldiv,0,menu_tree_hilldiv_ultra)
# FaithsPD_otu_sample
# 
# Faiths.df<-rownames_to_column(as.data.frame(FaithsPD_otu_sample),var="Sample")

#Using taxa-merged object instead of just straight ASVs

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

################### taxonomy merged analysis #######
###Menu: Get Hill diversity q=0 considering phylogenetic diversity for each taxa-merged sample----
#Faith's PD (presence/absence)----

FaithsPD_taxa_sample<-hill_div(Menu_taxmerge_eDNAindex_otudf,0,Menu_taxmerge_eDNAindex_ultra_tree)
FaithsPD_taxa_sample

Faiths.df<-rownames_to_column(as.data.frame(FaithsPD_otu_sample),var="Sample")

#Allen’s H using eDNAindex transformed taxa-merged read abundances----
AllensH_taxa_sample_eDNAindex<-hill_div(Menu_taxmerge_eDNAindex_otudf,1,Menu_taxmerge_eDNAindex_ultra_tree)
AllensH_otu_sample_eDNAindex

Allens.df<-rownames_to_column(as.data.frame(AllensH_otu_sample_eDNAindex),var="Sample")



PDotusample<-full_join(Faiths.df,Allens.df)


Richness_otu_sample<-hill_div(Menu_taxmerge_eDNAindex_otudf,0)
Richness_otu_sample

Richness.df<-rownames_to_column(as.data.frame(Richness_otu_sample),var="Sample")

richdiv_otusample.df<-full_join(PDotusample,Richness.df)

Shannonexp_otu_sample<-hill_div(Menu_taxmerge_eDNAindex_otudf,1)
Shannonexp_otu_sample

Shannonexp.df<-rownames_to_column(as.data.frame(Shannonexp_otu_sample),var="Sample")

richdiv_otusample.df<-full_join(richdiv_otusample.df,Shannonexp.df)

sample_df_eDNAindex_hilldiv<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))

sample_df_eDNAindex_hilldiv<-rownames_to_column(sample_df_eDNAindex_hilldiv,var = "Sample")%>%
  select(Sample, DD_lat, DD_long, Site_Name, Microhabitat)

Menu_otu_hilldiv_by_sample<-full_join(sample_df_eDNAindex_hilldiv,richdiv_otusample.df)

#write to csv and map in Tableau
write.csv(Menu_otu_hilldiv_by_sample, "../11_Vegan/menu_otu_diversity_by_sample.csv")

###Average these sample-wise numbers by site----

Menu_taxmerge_eDNAindex_envvars<-rownames_to_column(Menu_taxmerge_eDNAindex_envvars, var="Sample")

Menu_taxmerge_eDNAindex.df<-left_join(Menu_otu_hilldiv_by_sample,Menu_taxmerge_eDNAindex_envvars)

Menu_diversity_averagesample_bysite<-Menu_taxmerge_eDNAindex.df%>%
  group_by(Site_Name)%>%
  summarize(mean_FaithsPD = mean(FaithsPD_otu_sample),
            sd_FaithsPD=sd(FaithsPD_otu_sample),
            mean_AllensH=mean(AllensH_otu_sample_eDNAindex),
            sd_AllensH=sd(AllensH_otu_sample_eDNAindex),
            mean_Richness=mean(Richness_otu_sample),
            sd_Richness=sd(Richness_otu_sample),
            mean_Shannonexp=mean(Shannonexp_otu_sample),
            sd_Shannonexp=sd(Shannonexp_otu_sample),
            var_FaithsPD = var(FaithsPD_otu_sample),
            var_AllensH=var(AllensH_otu_sample_eDNAindex),
            var_Richness=var(Richness_otu_sample),
            var_Shannonexp=var(Shannonexp_otu_sample),
            mean_site_water_temp=mean(Water_temp,na.rm=TRUE),
            var_site_water_temp=var(Water_temp,na.rm=TRUE),
            mean_primary_prod=mean(primary_productivity),
            var_primary_prod=var(primary_productivity),
            hammerhead_site=first(hammerhead_site),
            carcharhinus_site=first(blacktip_site))%>%
  mutate(hammerhead_site=if_else(hammerhead_site=="hh_site",
                                 1,0),
         carcharhinus_site=if_else(carcharhinus_site=="carch_site",
                                   1,0))




###Add some sample information, write to csv and map in Tableau, analyze with corrplot----
Menu_site.df<-Menu_taxmerge_eDNAindex_envvars%>%
  dplyr::select(Site_Name,DD_lat,DD_long)%>%
  group_by(Site_Name)%>%
  summarize(meanLAT=mean(DD_lat),
            meanLON=mean(DD_long))

Menu_diversity_averagesamplebysite.df<-inner_join(Menu_site.df,Menu_diversity_averagesample_bysite, multiple = "any",keep = FALSE)
 



Menu_diversity_averagesamplebysite.df<-column_to_rownames(Menu_diversity_averagesamplebysite.df, var="Site_Name")


Menu_diversity_averagesamplebysite.df1<-Menu_diversity_averagesamplebysite.df%>%
  dplyr::select(mean_FaithsPD,mean_AllensH,mean_Richness,mean_Shannonexp,mean_site_water_temp,var_site_water_temp,mean_primary_prod,var_primary_prod,meanLAT,meanLON,hammerhead_site,carcharhinus_site)


Menu_cor<-cor(Menu_diversity_averagesamplebysite.df1)
head(round(Menu_cor,2))


corrplot(Menu_cor, method="circle",type="lower")

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
p.mat <- cor.mtest(Menu_cor)
head(p.mat[, 1:5])

# mark out the ones with p-values greater than 0.05
corrplot(Menu_cor, type="lower", 
         p.mat = p.mat, sig.level = 0.05)

# Leave blank on no significant coefficient
corrplot(Menu_cor, type="lower", 
         p.mat = p.mat, sig.level = 0.05, insig = "blank")

#Best correlation plot----
corrplot(Menu_cor, type="lower",# order = "hclust",#addCoefasPercent = TRUE, 
         addCoef.col = "black",number.cex = .5, # Add coefficient of correlation
         tl.col="black", tl.srt=45,tl.cex = .5, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=TRUE,
         #title = "Correlation of Crustacean and Fish Diversity Measures with Temperature, Primary Productivity, and Shark Nursery Status by Site" 
)




#Make dataframes for mapping average taxa/species richness per sample per site, average otu Faith's PD per sample per site, same for Allen's H and species/taxa exponent of Shannon's diversity for both fish and crustaceans separately.   ----


#MiFish hilldiv----
#get just the MiFish results phyloseq
MiFish.ps<-readRDS("../09_Metabar_to_Phyloseq/MiFish_Phyloseq.rds")
#filter out passive mid-bay samples because they didn't seem to work well
MiFish1.ps<- ps_filter(MiFish.ps, Microhabitat != "passive mid bay")

MiFish1.ps

MiFish_eDNAindex.ps<-phyloseq_standardize_otu_abundance(MiFish1.ps,method = "wisconsin")

#make otu table and ultrametric tree for MiFish ASVs
MiFish_otu_df<-cbind(data.frame(otu_table(MiFish_eDNAindex.ps)))
#Make ultrametric tree for otu object 
MiFish_tree<-MiFish_eDNAindex.ps@phy_tree
MiFish_ultra_tree<-chronos(MiFish_tree,lambda = 1,model = "relaxed")
write.tree(MiFish_ultra_tree, 'MiFish_ultra.tre')
MiFish_ultra_tree <- ape::read.tree('MiFish_ultra.tre')



MiFish_eDNAindex.ps<-phyloseq_standardize_otu_abundance(MiFish1.ps, method = "wisconsin")
MiFish_otu_df_eDNAindex<-cbind(data.frame(otu_table(MiFish_eDNAindex.ps)))



#Merge_species but keep taxa that aren't ID'ed to species
MiFish_taxa_merge.ps <- phyloseq::tax_glom(MiFish1.ps, taxrank = "species", NArm = FALSE)


#rename with species names if available
MiFish_taxa_merge_tax_fix.ps<-MiFish_taxa_merge.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

MiFish_taxa_merge_species_named.ps<-tax_rename(MiFish_taxa_merge_tax_fix.ps,rank = "species")

MiFish_taxa_merge_eDNAindex.ps<-phyloseq_standardize_otu_abundance(MiFish_taxa_merge_species_named.ps,method = "wisconsin")

MiFish_taxa_df_eDNAindex<-cbind(data.frame(otu_table(MiFish_taxa_merge_eDNAindex.ps)))
#Make ultrametric tree for taxa merged object 
MiFish_taxamerge_tree<-MiFish_taxa_merge_eDNAindex.ps@phy_tree
MiFish_taxamerge_ultra_tree<-chronos(MiFish_taxamerge_tree,lambda = 1,model = "relaxed")
write.tree(MiFish_taxamerge_ultra_tree, 'MiFish_taxamerge_ultra.tre')
MiFish_taxamerge_ultra_tree <- ape::read.tree('MiFish_taxamerge_ultra.tre')
#make taxa_df names match tips of the tree
mfnames<-rownames(MiFish_taxa_df_eDNAindex)
mfnames<-gsub(pattern = " ","_",mfnames)
rownames(MiFish_taxa_df_eDNAindex)<-mfnames

##do Hill diversity indices for MiFish----
###Allen’s H using eDNAindex transformed taxamerged otu read abundances----

MiFish_AllensH_otu_sample_eDNAindex<-hill_div(MiFish_taxa_df_eDNAindex,1,MiFish_taxamerge_ultra_tree)
MiFish_AllensH_otu_sample_eDNAindex

MiFish_Allens.df<-rownames_to_column(as.data.frame(MiFish_AllensH_otu_sample_eDNAindex),var="Sample")


###MiFish: Get Hill diversity q=0 considering phylogenetic diversity for each sample----
####Faith's PD (presence/absence)----

MiFish_FaithsPD_otu_sample_eDNAindex<-hill_div(MiFish_taxa_df_eDNAindex,0,MiFish_taxamerge_ultra_tree)
MiFish_FaithsPD_otu_sample_eDNAindex

MiFish_FaithsPD.df<-rownames_to_column(as.data.frame(MiFish_FaithsPD_otu_sample_eDNAindex),var="Sample")

MiFish_phylohilldiv_otus<-full_join(MiFish_FaithsPD.df,MiFish_Allens.df)

#Now do the richness measures of taxa, not ASVs
MiFish_Richness_taxa_sample<-hill_div(MiFish_taxa_df_eDNAindex,0)
MiFish_Richness_taxa_sample

MiFish_Richness.df<-rownames_to_column(as.data.frame(MiFish_Richness_taxa_sample),var="Sample")

MiFish_richdiv_otusample.df<-full_join(MiFish_phylohilldiv_otus,MiFish_Richness.df)

MiFish_Shannonexp_taxa_sample<-hill_div(MiFish_taxa_df_eDNAindex,1)
MiFish_Shannonexp_taxa_sample

MiFish_Shannonexp.df<-rownames_to_column(as.data.frame(MiFish_Shannonexp_taxa_sample),var="Sample")

MiFish_richdiv_otusample.df<-full_join(MiFish_richdiv_otusample.df,MiFish_Shannonexp.df)


MiFish_sample_df<-cbind(data.frame(sample_data(MiFish_taxa_merge_eDNAindex.ps)))

MiFish_sample_df<-rownames_to_column(MiFish_sample_df,var = "Sample")%>%
  select(Sample, DD_lat, DD_long, Site_Name, Microhabitat)

MiFish_diversity_by_sample.df<-full_join(MiFish_sample_df,MiFish_richdiv_otusample.df)




#write to csv and map in Tableau
write.csv(MiFish_diversity_by_sample.df, "../11_Vegan/MiFish_otu_and_taxa_diversity_by_sample.csv")
###Average these sample-wise measures by site----
MiFish_diversity_averagesample_bysite<-MiFish_diversity_by_sample.df%>%
  group_by(Site_Name)%>%
  summarize(MiFish_mean_FaithsPD = mean(MiFish_FaithsPD_otu_sample_eDNAindex),
            MiFish_sd_FaithsPD=sd(MiFish_FaithsPD_otu_sample_eDNAindex),
            MiFish_mean_AllensH=mean(MiFish_AllensH_otu_sample_eDNAindex),
            MiFish_sd_AllensH=sd(MiFish_AllensH_otu_sample_eDNAindex),
            MiFish_mean_Richness=mean(MiFish_Richness_taxa_sample),
            MiFish_sd_Richness=sd(MiFish_Richness_taxa_sample),
            MiFish_mean_Shannonexp=mean(MiFish_Shannonexp_taxa_sample),
            MiFish_sd_Shannonexp=sd(MiFish_Shannonexp_taxa_sample))

#add some sample information, write to csv and map in Tableau
MiFish_site.df<-MiFish_sample_df%>%
  select(Site_Name,DD_lat,DD_long)%>%
  group_by(Site_Name)%>%
  summarize(meanLAT=mean(DD_lat),
            meanLON=mean(DD_long))

MiFish_diversity_averagesamplebysite.df<-inner_join(MiFish_site.df,MiFish_diversity_averagesample_bysite, multiple = "any",keep = FALSE)

write.csv(MiFish_diversity_averagesamplebysite.df, "../11_Vegan/MiFish_average_sample_diversity_by_site.csv")



#BerryCrust hilldiv----
#get just the BerryCrust results phyloseq
BerryCrust.ps<-readRDS("../09_Metabar_to_Phyloseq/BerryCrust_Phyloseq.rds")
#filter out passive mid-bay samples because they didn't seem to work well
BerryCrust1.ps<- ps_filter(BerryCrust.ps, Microhabitat != "passive mid bay")

BerryCrust_eDNAindex.ps<-phyloseq_standardize_otu_abundance(BerryCrust1.ps,method = "wisconsin")


#make otu table and ultrametric tree for BerryCrust ASVs
BerryCrust_otu_df<-cbind(data.frame(otu_table(BerryCrust_eDNAindex.ps)))
#Make ultrametric tree for otu object 
BerryCrust_tree<-BerryCrust_eDNAindex.ps@phy_tree
BerryCrust_ultra_tree<-chronos(BerryCrust_tree,lambda = 1,model = "relaxed")
write.tree(BerryCrust_ultra_tree, 'BerryCrust_ultra.tre')
BerryCrust_ultra_tree <- ape::read.tree('BerryCrust_ultra.tre')

#Merge_species but keep taxa that aren't ID'ed to species
BerryCrust_taxa_merge.ps <- phyloseq::tax_glom(BerryCrust1.ps, taxrank = "species", NArm = FALSE)
#rename with species names if available
BerryCrust_taxa_merge_tax_fix.ps<-BerryCrust_taxa_merge.ps %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

BerryCrust_taxa_merge_species_named.ps<-tax_rename(BerryCrust_taxa_merge_tax_fix.ps,rank = "species")

BerryCrust_taxa_merge_eDNAindex.ps<-phyloseq_standardize_otu_abundance(BerryCrust_taxa_merge_species_named.ps,method = "wisconsin")

#make otu table and for BerryCrust taxa
BerryCrust_taxa_otu_df<-cbind(data.frame(otu_table(BerryCrust_taxa_merge_eDNAindex.ps)))

#Make ultrametric tree for otu object 
BerryCrust_taxamerge_tree<-BerryCrust_taxa_merge_eDNAindex.ps@phy_tree
BerryCrust_taxamerge_ultra_tree<-chronos(BerryCrust_taxamerge_tree,lambda = 1,model = "relaxed")
write.tree(BerryCrust_taxamerge_ultra_tree, 'BerryCrust_taxamerge_ultra.tre')
BerryCrust_taxamerge_ultra_tree <- ape::read.tree('BerryCrust_taxamerge_ultra.tre')

#make names of BerryCrust_taxa_otu_df match the tip labels
bcnames<-rownames(BerryCrust_taxa_otu_df)
bcnames<-gsub(pattern = " ","_",bcnames)
rownames(BerryCrust_taxa_otu_df)<-bcnames

##do Hill diversity indices for BerryCrust----
###Allen’s H using eDNAindex transformed otu read abundances----

BerryCrust_AllensH_otu_sample_eDNAindex<-hill_div(BerryCrust_taxa_otu_df,1,BerryCrust_taxamerge_ultra_tree)
BerryCrust_AllensH_otu_sample_eDNAindex

BerryCrust_Allens.df<-rownames_to_column(as.data.frame(BerryCrust_AllensH_otu_sample_eDNAindex),var="Sample")


###BerryCrust: Get Hill diversity q=0 considering phylogenetic diversity for each sample----
####Faith's PD (presence/absence)----

BerryCrust_FaithsPD_otu_sample_eDNAindex<-hill_div(BerryCrust_taxa_otu_df,0,BerryCrust_taxamerge_ultra_tree)
BerryCrust_FaithsPD_otu_sample_eDNAindex

BerryCrust_FaithsPD.df<-rownames_to_column(as.data.frame(BerryCrust_FaithsPD_otu_sample_eDNAindex),var="Sample")

BerryCrust_phylohilldiv_otus<-full_join(BerryCrust_FaithsPD.df,BerryCrust_Allens.df)

#Now do the richness measures of taxa, not ASVs
BerryCrust_Richness_taxa_sample<-hill_div(BerryCrust_taxa_otu_df,0)
BerryCrust_Richness_taxa_sample

BerryCrust_Richness.df<-rownames_to_column(as.data.frame(BerryCrust_Richness_taxa_sample),var="Sample")

BerryCrust_richdiv_otusample.df<-full_join(BerryCrust_phylohilldiv_otus,BerryCrust_Richness.df)

BerryCrust_Shannonexp_taxa_sample<-hill_div(BerryCrust_taxa_otu_df,1)
BerryCrust_Shannonexp_taxa_sample

BerryCrust_Shannonexp.df<-rownames_to_column(as.data.frame(BerryCrust_Shannonexp_taxa_sample),var="Sample")

BerryCrust_richdiv_otusample.df<-full_join(BerryCrust_richdiv_otusample.df,BerryCrust_Shannonexp.df)


BerryCrust_sample.df<-cbind(data.frame(sample_data(BerryCrust_taxa_merge_eDNAindex.ps)))
Menu_sample.df<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))
Menu_sample.df<-rownames_to_column(Menu_sample.df, var="Sample")

BerryCrust_sample_df<-rownames_to_column(BerryCrust_sample.df,var = "Sample")%>%
  dplyr::select(Sample, DD_lat, DD_long, Site_Name, Microhabitat,Water_temp)

BerryCrust_sample1.df<-left_join(BerryCrust_sample_df,Menu_sample.df)
BerryCrust_sample2.df<-BerryCrust_sample1.df%>%
  dplyr::select(Sample, DD_lat, DD_long, Site_Name, Microhabitat,Water_temp,Tide,hammerheads_detected,carcharhinus_detected,hammerhead_site,blacktip_site, sst, depth, silicate, oxygen, phosphate, primary_productivity, nitrate,chloride)

BerryCrust_diversity_by_sample.df<-full_join(BerryCrust_sample2.df,BerryCrust_richdiv_otusample.df)

#write to csv and map in Tableau
write.csv(BerryCrust_diversity_by_sample.df, "../11_Vegan/BerryCrust_otu_and_taxa_diversity_by_sample.csv")

###Average these sample-wise numbers by site----

BerryCrust_diversity_averagesample_bysite<-BerryCrust_diversity_by_sample.df%>%
  group_by(Site_Name)%>%
  summarize(BerryCrust_mean_FaithsPD = mean(BerryCrust_FaithsPD_otu_sample_eDNAindex),
            BerryCrust_sd_FaithsPD=sd(BerryCrust_FaithsPD_otu_sample_eDNAindex),
            BerryCrust_mean_AllensH=mean(BerryCrust_AllensH_otu_sample_eDNAindex),
            BerryCrust_sd_AllensH=sd(BerryCrust_AllensH_otu_sample_eDNAindex),
            BerryCrust_mean_Richness=mean(BerryCrust_Richness_taxa_sample),
            BerryCrust_sd_Richness=sd(BerryCrust_Richness_taxa_sample),
            BerryCrust_mean_Shannonexp=mean(BerryCrust_Shannonexp_taxa_sample),
            BerryCrust_sd_Shannonexp=sd(BerryCrust_Shannonexp_taxa_sample),
            BerryCrust_var_FaithsPD = var(BerryCrust_FaithsPD_otu_sample_eDNAindex),
            BerryCrust_var_AllensH=var(BerryCrust_AllensH_otu_sample_eDNAindex),
            BerryCrust_var_Richness=var(BerryCrust_Richness_taxa_sample),
            BerryCrust_var_Shannonexp=var(BerryCrust_Shannonexp_taxa_sample),
            mean_site_water_temp=mean(Water_temp,na.rm=TRUE),
            var_site_water_temp=var(Water_temp,na.rm=TRUE),
            mean_primary_prod=mean(primary_productivity),
            var_primary_prod=var(primary_productivity),
            mean_nitrate=mean(nitrate),
            var_nitrate=var(nitrate),
            hammerhead_site=first(hammerhead_site),
            carcharhinus_site=first(blacktip_site))%>%
  mutate(hammerhead_site=if_else(hammerhead_site=="hh_site",
                                 1,0),
         carcharhinus_site=if_else(carcharhinus_site=="carch_site",
                                   1,0))
            
            


###Add some sample information, write to csv and map in Tableau, analyze with corrplot----
BerryCrust_site.df<-BerryCrust_sample_df%>%
  dplyr::select(Site_Name,DD_lat,DD_long)%>%
  group_by(Site_Name)%>%
  summarize(meanLAT=mean(DD_lat),
            meanLON=mean(DD_long))

BerryCrust_diversity_averagesamplebysite.df<-inner_join(BerryCrust_site.df,BerryCrust_diversity_averagesample_bysite, multiple = "any",keep = FALSE)

write.csv(BerryCrust_diversity_averagesamplebysite.df, "../11_Vegan/BerryCrust_average_sample_diversity_by_site.csv")


####BerryCrust div test by site----

#For phylogenetic analyses, use OTUs, BerryCrust_otu_df,BerryCrust_ultra_tree,
#For taxonomic diversity analyses use taxonomic groups BerryCrust_taxa_otu_df
div_test(BerryCrust_taxa_otu_df,qvalue=0,hierarchy=menu_sample_site_hierarchy)#richness
div_test(BerryCrust_otu_df,qvalue=1,hierarchy=menu_sample_site_hierarchy,tree=BerryCrust_ultra_tree) #Allen's H
div_test(BerryCrust_otu_df,qvalue=0,hierarchy=menu_sample_site_hierarchy,tree=BerryCrust_ultra_tree)#Faith's PD
div_test(BerryCrust_otu_df,2,menu_sample_site_hierarchy,BerryCrust_ultra_tree)#Rao's
div_test(BerryCrust_taxa_otu_df,qvalue=1,hierarchy=menu_sample_site_hierarchy,posthoc=TRUE)#Shannon's exponent

#no significant p-values by Site_Name


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

##Phylogenetic Hill q=0 (Faith's PD) distance matrix with taxa merged object ----

Menu_taxmerged_Hilldistmat<-pair_dis(countable = Menu_taxmerge_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1,tree = Menu_taxmerge_eDNAindex_ultra_tree)



Menu_taxmerged_Hilldistmat$L1_SqN
#make it symmetric.
Menu_taxmerged_Hill_Sor_dissim<-as.dist(Menu_taxmerged_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#Allen's H distance matrix (q=1)
Allensdistmat<-pair_dis(countable = Menu_taxmerge_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1,tree = Menu_taxmerge_eDNAindex_ultra_tree)

Allensdistmat$L1_SqN
#make it symmetric.
Allensdistmat_Sor_dissim<-as.dist(Allensdistmat$L1_VqN,diag = TRUE, upper = TRUE)
Allensdistmat_Jaccard_dissim<-as.dist(Allensdistmat$L1_SqN,diag = TRUE, upper = TRUE)



vegan::adonis2(Menu_taxmerged_Hill_Jaccard_dissim ~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)$Site_Name*Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Microhabitat)



set.seed(1234)
vegan::adonis2(Menu_taxmerged_Hill_Jaccard_dissim ~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)$Site_Name+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$hammerheads_detected+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$carcharhinus_detected+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Microhabitat+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Tide
               +Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$scaled_sst)

vegan::adonis2(Menu_taxmerged_Hill_Jaccard_dissim ~Microhabitat+scaled_Water_temp+Tide+scaled_oxygen+scaled_primaryprod,data = Menu_taxmerge_eDNAindex_envvars,by="margin",na.action = na.omit,sqrt.dist = FALSE)

#for q=0 Faith's PD, Site_Name, Microhabitat, and scaled sst significant using Sor or Jaccard type distance matrix
#only Microhabitat and scaled sst significant for Allen's H Jaccard-type dissim matrix

#q=1 Allen's H
# vegan::adonis2(formula = Menu_taxmerged_Hill_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_oxygen + scaled_primaryprod, data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, by = "margin", na.action = na.omit)
# Df SumOfSqs      R2      F Pr(>F)
# Microhabitat        2  0.08584 0.13532 5.2247  0.005 **
#   scaled_Water_temp   1  0.02588 0.04079 3.1499  0.060 .
# Tide                1  0.00903 0.01423 1.0987  0.349
# scaled_oxygen       1  0.07502 0.11826 9.1320  0.004 **
#   scaled_primaryprod  1  0.04029 0.06352 4.9048  0.019 *
#   Residual           51  0.41897 0.66044
# Total              57  0.63439 1.00000
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#q=0 Faith's PD
# vegan::adonis2(formula = Menu_taxmerged_Hill_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_oxygen + scaled_primaryprod, data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, by = "margin", na.action = na.omit)
# Df SumOfSqs      R2      F Pr(>F)
# Microhabitat        2   0.2999 0.04892 1.6865  0.094 .
# scaled_Water_temp   1   0.3326 0.05425 3.7405  0.013 *
#   Tide                1   0.1066 0.01739 1.1992  0.309
# scaled_oxygen       1   0.4139 0.06753 4.6559  0.010 **
#   scaled_primaryprod  1   0.3246 0.05295 3.6506  0.021 *
#   Residual           51   4.5343 0.73971
# Total              57   6.1299 1.00000
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#for q=1 Shannon's exp non-phylogenetic, Microhabitat, scaled o2, Tide were significant, water temp almost, but when replaced with sst its significant...

#Check for differences between hammerhead sites and non-hammerhead sites.----


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
Menu_taxmerged_rmhh_Shannons_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rmhh_eDNAindex_otudf,qvalue = 1,hierarchy=menu_sample_site_hierarchy,level = 1)#,tree = Menu_taxmerge_rmhh_eDNAindex_ultra_tree)

#get Jaccard dissimilarity matrix
Menu_taxmerged_rmhh_Shannons_Hilldistmat$L1_SqN
#make it symmetric.
Menu_taxmerged_rmhh_Richness_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_rmhh_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)


set.seed(200)
pairwise.adonis(Menu_taxmerged_rmhh_Richness_Hill_Jaccard_dissim, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_hh.eDNAindex.ps)$hammerhead_site)

##Check for differences between carcharhinus sites and non-carcharhinus sites.----


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
namechange_rm_carch<-gsub(pattern = " ",replacement = "_",namechange_rmhh)
rownames(Menu_taxmerge_rm_carch_eDNAindex_otudf)<-namechange_rm_carch

#make sample table for later
rm_carch_eDNAindex_envvars<-as.data.frame(cbind(sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps)))

##Hill distance matrix with taxa merged rm_carch object ----
Menu_taxmerged_rm_carch_Richness_Hilldistmat<-pair_dis(countable = Menu_taxmerge_rm_carch_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1)#,tree = Menu_taxmerge_rm_carch_eDNAindex_ultra_tree)

#get Jaccard dissimilarity matrix
Menu_taxmerged_rm_carch_Hilldistmat$L1_SqN
#make it symmetric.
Menu_taxmerged_rm_carch_Richness_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)


set.seed(200)
pairwise.adonis(Menu_taxmerged_rm_carch_Richness_Hill_Jaccard_dissim, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps)$blacktip_site)






#cool network graph of Jaccard-type dissimilarity of sites' Allen's H
pair_dis_plot(Menu_taxmerged_Hill_Jaccard_dissim,hierarchy=menu_sample_site_hierarchy,type="qgraph") #cool network graph of Jaccard-type dissimilarity of sites' Allen's H

# tree_depth(tree = menu_tree_hilldiv_ultra, abund = otu_df_eDNAindex_hilldiv)
# tree_depth(tree= Menu_taxmerge_eDNAindex_ultra_tree, abund = Menu_taxmerge_eDNAindex_otudf)


Menu_taxmerged_Hilldistmat_invSimpson<-pair_dis(countable = Menu_taxmerge_eDNAindex_otudf,qvalue = 2,hierarchy=menu_sample_site_hierarchy,level = 1)

Menu_taxmerged_Hilldistmat_invSimpson$L1_SqN
#make it symmetric.
Menu_taxmerged_Hilldistmat_invSimpson_dissim<-as.dist(Menu_taxmerged_Hilldistmat_invSimpson$L1_SqN,diag = TRUE, upper = TRUE)

vegan::adonis2(Menu_taxmerged_Hilldistmat_invSimpson_dissim ~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)$Site_Name*Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Microhabitat)


set.seed(200)
vegan::adonis2(Menu_taxmerged_Hill_Shannonexp_Jaccard_dissim ~ phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)$Site_Name+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$hammerheads_detected+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$carcharhinus_detected+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Microhabitat+Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Tide
               +Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$scaled_sst)

#without taking phylogeny into account, keeping all other things the same(q=1, L1_SqN) now everything except hammerheads detected is highly significant.

# phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)$Site_Name 1.8319  0.001 ***
#   Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$hammerheads_detected    1.0590  0.381    
# Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$carcharhinus_detected   2.0432  0.001 ***
#   Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Microhabitat            2.2682  0.001 ***
#   Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$Tide                    1.5673  0.027 *  
#   Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@sam_data$scaled_sst              2.4031  0.001 ***

#Incidence-based----
bat.diet.otutable.incidence <- to.incidence(bat.diet.otutable,bat.diet.hierarchy)
hill_div(bat.diet.otutable.incidence,qvalue=1)
hill_div(to.incidence(bat.diet.otutable,bat.diet.hierarchy),1)

#make hierarchy df for menu samples by site and by microhabitat
menu_sample_data<-cbind(data.frame(sample_data(Menu1.ps)))
menu_sample_site_hierarchy<-menu_sample_data%>%
  select(sample_id,Site_Name)

menu_sample_microhabitat_hierarchy<-menu_sample_data%>%
  select(sample_id,Microhabitat)

menu_otu_site_incidence <- to.incidence(otu_df_hilldiv,menu_sample_site_hierarchy)
hill_div(menu_otu_site_incidence,qvalue=0)
hill_div(to.incidence(otu_df_hilldiv,menu_sample_site_hierarchy),1)

#add tree information

FaithsPD_otu_incidence_site<-hill_div(menu_otu_site_incidence,qvalue=0,tree =menu_tree_hilldiv_ultra )
FaithsPD_otu_incidence_site


##index div----
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)

#One sample
bat.diet.sample <- bat.diet.otutable[,1]
index_div(bat.diet.sample)
index_div(bat.diet.sample,index="shannon")

#Multiple samples
index_div(bat.diet.otutable)
index_div(bat.diet.otutable,tree=bat.diet.tree,index="faith")

# ###Menu: Make an otu by site Faith's PD calculation----
# #Merge by site
# Menu_site_merge.ps<-merge_samples(Menu1.ps, group = "Site_Name", fun=sum)
# #calculate eDNA index by site
# Menu_site_merge_eDNAindex.ps<-phyloseq_standardize_otu_abundance(Menu_site_merge.ps, method = "wisconsin")
# #create dataframe of otu table
# site_merge_hilldiv<-cbind(data.frame(otu_table(Menu_site_merge.ps)))
# #make ultrametric tree from phyloseq tree object
# menu_tree_hilldiv_site<-Menu_site_merge.ps@phy_tree
# menu_tree_hilldiv_site_ultra<-chronos(menu_tree_hilldiv_site,lambda = 1,model = "correlated")
# class(menu_tree_hilldiv_site_ultra)
# write.tree(menu_tree_hilldiv_site_ultra, 'site_ultra.tre')
# menu_tree_hilldiv_site_ultra <- ape::read.tree('site_ultra.tre')
# class(menu_tree_hilldiv_site_ultra)
# 
# #Calculate Faith's phylogenetic diversity by site
#  FaithsPD_otu_site<-index_div(site_merge_hilldiv,tree=menu_tree_hilldiv_site_ultra,index="faith")
#  FaithsPD_otu_site
# 

#Incidence-based
bat.diet.otutable.incidence <- to.incidence(bat.diet.otutable,bat.diet.hierarchy)
index_div(bat.diet.otutable.incidence)
index_div(bat.diet.otutable.incidence,index="simpson")
index_div(to.incidence(bat.diet.otutable,bat.diet.hierarchy),tree=bat.diet.tree)

##div test----
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)

div_test(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)
div_test(bat.diet.otutable,qvalue=1,hierarchy=bat.diet.hierarchy,tree=bat.diet.tree)
div_test(bat.diet.otutable,2,bat.diet.hierarchy,bat.diet.tree)
div_test(bat.diet.otutable,qvalue=1,hierarchy=bat.diet.hierarchy,posthoc=TRUE)


###Menu div test by site----
div_test(otu_df_hilldiv,qvalue=0,hierarchy=menu_sample_site_hierarchy)#richness
div_test(otu_df_hilldiv,qvalue=1,hierarchy=menu_sample_site_hierarchy,tree=menu_tree_hilldiv_ultra) #Allen's H
div_test(otu_df_hilldiv,qvalue=0,hierarchy=menu_sample_site_hierarchy,tree=menu_tree_hilldiv_ultra)#Faith's PD
div_test(otu_df_hilldiv,2,menu_sample_site_hierarchy,menu_tree_hilldiv_ultra)#Rao's
div_test(otu_df_hilldiv,qvalue=1,hierarchy=menu_sample_site_hierarchy,posthoc=TRUE)


#use eDNA index taxmerged species named final phyloseq object
#Menu_taxmerge_eDNAindex_otudf, Menu_taxmerge_eDNAindex_ultra_tree, Menu_taxmerge_eDNAindex_envvars


##div test plot----

#Menu: plot Faith's phylogenetic diversity by site and by microhabitat
#With no post-hoc analyses

divtestres <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=0,hierarchy=menu_sample_site_hierarchy,tree = Menu_taxmerge_eDNAindex_ultra_tree) #Faith's PD
div_test_plot(divtestres,chart="box")
div_test_plot(divtestres,chart="violin")
div_test_plot(divtestres,chart="jitter")


#With post-hoc analyses
divtest.res.ph <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=0,hierarchy=menu_sample_site_hierarchy,tree =Menu_taxmerge_eDNAindex_ultra_tree, posthoc=TRUE)#Faith'sPD by site
div_test_plot(divtest.res.ph,chart="box")#,posthoc=TRUE,threshold = 0.05)

divtest.res1.ph <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=1,hierarchy=menu_sample_site_hierarchy,posthoc=TRUE),tree = Menu_taxmerge_eDNAindex_ultra_tree)#Allen's H by site
div_test_plot(divtest.res1.ph,chart="box")#,posthoc=TRUE,threshold=0.05)

divtest.res2.ph <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=0,hierarchy=menu_sample_site_hierarchy,posthoc=TRUE)#Richness and Shannnons by site not sig
div_test_plot(divtest.res2.ph,chart="box")#,posthoc=TRUE,threshold=0.05)
#Menu: plot Faith's phylogenetic diversity by microhabitat
#With no post-hoc analyses

divtestres <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=0,hierarchy=menu_sample_microhabitat_hierarchy,tree = Menu_taxmerge_eDNAindex_ultra_tree)
div_test_plot(divtestres,chart="box")
div_test_plot(divtestres,chart="violin")
div_test_plot(divtestres,chart="jitter")

div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=0,hierarchy=menu_sample_microhabitat_hierarchy,tree = Menu_taxmerge_eDNAindex_ultra_tree)

divtestres <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=0,hierarchy=menu_sample_microhabitat_hierarchy,tree = Menu_taxmerge_eDNAindex_ultra_tree)
div_test_plot(divtestres,chart="box")

#With post-hoc analyses
divtest.res.ph <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=1,hierarchy=menu_sample_microhabitat_hierarchy,posthoc=TRUE,tree = Menu_taxmerge_eDNAindex_ultra_tree)#Allen's H by microhabitat
div_test_plot(divtest.res.ph,chart="box",posthoc=TRUE,threshold=0.05)


divtest.res3.ph <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=0,hierarchy=menu_sample_microhabitat_hierarchy,posthoc=TRUE)#Richness by microhabitat
div_test_plot(divtest.res3.ph,chart="box")#,posthoc=TRUE,threshold=0.05)

#hhsite vs non-hhsite
menu_sample_data1<-rownames_to_column(menu_sample_data,var="sample_id")
menu_sample_hhsite_hierarchy<-menu_sample_data1%>%
  dplyr::select(sample_id,hammerhead_site)
rownames(menu_sample_hhsite_hierarchy)<-menu_sample_hhsite_hierarchy$sample_id


divtestres <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=0,hierarchy=menu_sample_hhsite_hierarchy,tree = Menu_taxmerge_eDNAindex_ultra_tree) #Faith's PD
div_test_plot(divtestres,chart="box")
div_test_plot(divtestres,chart="violin")
div_test_plot(divtestres,chart="jitter")

#With post-hoc analyses
divtest.res.ph <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=0,hierarchy=menu_sample_hhsite_hierarchy,posthoc = TRUE)#,tree = Menu_taxmerge_eDNAindex_ultra_tree)#Allen's H by hammerhead site non-sig
div_test_plot(divtest.res.ph,chart="box")#,posthoc=TRUE,threshold=0.05)



#carchsite vs non-carchsite
menu_sample_data1<-rownames_to_column(menu_sample_data,var="sample_id")
menu_sample_carchsite_hierarchy<-menu_sample_data1%>%
  dplyr::select(sample_id,blacktip_site)
rownames(menu_sample_carchsite_hierarchy)<-menu_sample_carchsite_hierarchy$sample_id

divtestres <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=1,hierarchy=menu_sample_carchsite_hierarchy)#,tree = Menu_taxmerge_eDNAindex_ultra_tree) #Shannonexp
div_test_plot(divtestres,chart="box")
div_test_plot(divtestres,chart="violin")
div_test_plot(divtestres,chart="jitter")

#With post-hoc analyses
divtest.res.ph <- div_test(Menu_taxmerge_eDNAindex_otudf,qvalue=1,hierarchy=menu_sample_carchsite_hierarchy,posthoc = TRUE)#,tree = Menu_taxmerge_eDNAindex_ultra_tree)#Allen's H by carchcharhinus site non-sig
div_test_plot(divtest.res.ph,chart="jitter")#,posthoc=TRUE,threshold=0.05)



#following Stier et al----
##1000 rarefaction replicates to get every sample to the same coverage and calculate diversity indices----
multi_rare_div<-phyloseq_mult_raref_div(physeq = Shark_Menu_taxa_merge_species_named_traits.ps,iter = 1000,verbose = TRUE, parallel = TRUE,SampSize = min(sample_sums(Shark_Menu_taxa_merge_species_named_traits.ps)))
#this just did Observed richness and Shannon diversity

#This creates a new phyloseq object after averaging the 1000 rarefaction iterations
multi_rare_avg.ps<-phyloseq_mult_raref_avg(physeq = Shark_Menu_taxa_merge_species_named_traits.ps,iter = 1000,verbose = TRUE, parallel = TRUE,SampSize = min(sample_sums(Shark_Menu_taxa_merge_species_named_traits.ps)))

otu_table(Shark_Menu_taxa_merge_species_named_traits.ps)
otu_table(multi_rare_avg.ps)

###Plot richness----
# plot_richness(multi_rare_avg.ps, x="hammerhead_site","Microhabitat", measures=c("Shannon", "Simpson"), color="Site_Name")+
#   scale_color_ucscgb()
# 

###eDNA index transform rarefied phyloseq object----
multi_rare_avg.eDNAindex.ps<-phyloseq_standardize_otu_abundance(multi_rare_avg.ps, method = "wisconsin")
otu_table(multi_rare_avg.eDNAindex.ps)

#create input objects for hilldiv from rarefied eDNA indexed taxa-merged phyloseq
Menu_taxmerge_raref_eDNAindex_otudf<-cbind(data.frame(otu_table(multi_rare_avg.ps)))
#Make ultrametric tree for otu object 
Menu_taxmerge_raref_eDNAindex_tree<-multi_rare_avg.ps@phy_tree
Menu_taxmerge_raref_eDNAindex_ultra_tree<-chronos(Menu_taxmerge_raref_eDNAindex_tree,lambda = 1,model = "relaxed")
write.tree(Menu_taxmerge_raref_eDNAindex_ultra_tree, 'Menu_taxmerge_raref_eDNAindex_ultra.tre')
Menu_taxmerge_raref_eDNAindex_ultra_tree <- ape::read.tree('Menu_taxmerge_raref_eDNAindex_ultra.tre')

#make rownames of otu table match tip names of the tree
rarefnamechange<-rownames(Menu_taxmerge_raref_eDNAindex_otudf)
rarefnamechange<-gsub(pattern = " ",replacement = "_",rarefnamechange)
rownames(Menu_taxmerge_raref_eDNAindex_otudf)<-rarefnamechange

#make envvars for later
Menu_taxmerge_raref_eDNAindex_envvars<-cbind(data.frame(sample_data(multi_rare_avg.ps)))

###calculate Hill alpha and beta diversity metrics for hammerhead sites and carcharhinus sites (and non-)----

#alpha =sample richness, #beta =jaccard from Stier et al. (incidence-based) Can do the same thing for alpha diversity that I did for the full unrarefied menu object above...

###Menu: Get Hill diversity q=0 considering phylogenetic diversity for each taxa-merged sample----
#Faith's PD (presence/absence)----

raref_FaithsPD_taxa_sample<-hill_div(Menu_taxmerge_raref_eDNAindex_otudf,0,Menu_taxmerge_eDNAindex_ultra_tree)
raref_FaithsPD_taxa_sample

raref_Faiths.df<-rownames_to_column(as.data.frame(raref_FaithsPD_taxa_sample),var="Sample")

#Allen’s H using eDNAindex transformed taxa-merged read abundances----
raref_AllensH_taxa_sample_eDNAindex<-hill_div(Menu_taxmerge_raref_eDNAindex_otudf,1,Menu_taxmerge_raref_eDNAindex_ultra_tree)
raref_AllensH_taxa_sample_eDNAindex

raref_Allens.df<-rownames_to_column(as.data.frame(raref_AllensH_taxa_sample_eDNAindex),var="Sample")



raref_PDotusample<-full_join(raref_Faiths.df,raref_Allens.df)


raref_Richness_otu_sample<-hill_div(Menu_taxmerge_raref_eDNAindex_otudf,0)
raref_Richness_otu_sample

raref_Richness.df<-rownames_to_column(as.data.frame(raref_Richness_otu_sample),var="Sample")

raref_richdiv_otusample.df<-full_join(raref_PDotusample,raref_Richness.df)

raref_Shannonexp_otu_sample<-hill_div(Menu_taxmerge_raref_eDNAindex_otudf,1)
raref_Shannonexp_otu_sample

raref_Shannonexp.df<-rownames_to_column(as.data.frame(raref_Shannonexp_otu_sample),var="Sample")

raref_richdiv_otusample.df<-full_join(raref_richdiv_otusample.df,raref_Shannonexp.df)

raref_sample_df_eDNAindex_hilldiv<-cbind(data.frame(sample_data(multi_rare_avg.ps)))

raref_sample_df_eDNAindex_hilldiv<-rownames_to_column(raref_sample_df_eDNAindex_hilldiv,var = "Sample")%>%
  dplyr::select(Sample, DD_lat, DD_long, Site_Name, Microhabitat,hammerhead_site,blacktip_site)

raref_Menu_otu_hilldiv_by_sample<-full_join(raref_sample_df_eDNAindex_hilldiv,raref_richdiv_otusample.df)


#diversity partitioning
div_part(countable = Menu_taxmerge_raref_eDNAindex_otudf,qvalue = 1,hierarchy = menu_sample_hhsite_hierarchy)#,tree = Menu_taxmerge_raref_eDNAindex_ultra_tree)
#Faith's PD
# $Hierarchical_levels
# [1] 3
# 
# $Order_diversity
# [1] 0
# 
# $Hill_numbers
# L1       L2       L3 
# 2.412765 5.184751 5.184751 
# 
# $Sample_size
# N1 N2 N3 
# 59  2  1 
# 
# $Beta
# B1_2     B2_3 
# 2.148884 1.000000 

#Shannon
# $Hierarchical_levels
# [1] 3
# 
# $Order_diversity
# [1] 1
# 
# $Hill_numbers
# L1        L2        L3 
# 5.453123 21.203602 23.741874 
# 
# $Sample_size
# N1 N2 N3 
# 59  2  1 
# 
# $Beta
# B1_2     B2_3 
# 3.888341 1.119709 


####calculate beta dispersion with correction----

#create dissimilarity matrix for betadisper to use just like the distance matrix used above.
#remove each shark:


Menu_taxmerged_Faiths_Hilldistmat<-pair_dis(countable = Menu_taxmerge_eDNAindex_otudf,qvalue = 0,hierarchy=menu_sample_site_hierarchy,level = 1,tree = Menu_taxmerge_eDNAindex_ultra_tree) 



Menu_taxmerged_Faiths_Hilldistmat$L1_SqN
#make it symmetric.
Menu_taxmerged_Faiths_Hill_Sor_dissim<-as.dist(Menu_taxmerged_Faiths_Hilldistmat$L1_VqN,diag = TRUE, upper = TRUE)
Menu_taxmerged_Faiths_Hill_Jaccard_dissim<-as.dist(Menu_taxmerged_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

#calculate dispersion around the centroid using the distance matrix above----
Faith_hhbetadisp<-betadisper(d = Menu_taxmerged_Faiths_Hill_Jaccard_dissim,group = Menu_taxmerge_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Faith_hhbetadisp)
permutest(Faith_hhbetadisp)

Richness_hhbetadisp<-betadisper(d = Menu_taxmerged_Richness_Hill_Jaccard_dissim,group = Menu_taxmerge_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Richness_hhbetadisp)
vegan::permutest(Richness_hhbetadisp)

Shannon_hhbetadisp<-betadisper(d=Menu_taxmerged_Hill_Shannonexp_Jaccard_dissim,group=Menu_taxmerge_eDNAindex_envvars$hammerhead_site,bias.adjust = TRUE,sqrt.dist = FALSE)
TukeyHSD(Shannon_hhbetadisp)
permutest(Shannon_hhbetadisp)


##div_profile----

data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)

#One sample example
bat.diet.sample <- bat.diet.otutable[,1]
div_profile(count=bat.diet.sample,qvalues=seq(from = 0, to = 5, by = (0.1)))

#One sample example (phylogenetic Hill numbers)
names(bat.diet.sample) <- rownames(bat.diet.otutable)
div_profile(count=bat.diet.sample,qvalues=seq(from = 0, to = 5, by = (0.1)),tree=bat.diet.tree)

#Multiple samples
div_profile(bat.diet.otutable)

#Multiple groups (gamma diversity)
div_profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="gamma")

#Multiple groups (alpha diversity)
div_profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="alpha")

##div profile plot----

data(bat.diet.otutable)
data(bat.diet.hierarchy)

#One sample example
bat.diet.sample <- bat.diet.otutable[,1]
profile.onesample <- div_profile(count=bat.diet.sample,qvalues=seq(from = 0, to = 5, by = (0.1)))
div_profile_plot(profile.onesample)

#Multiple samples
profile.multiplesamples <- div_profile(bat.diet.otutable)
div_profile_plot(profile.multiplesamples)

#Multiple samples
profile.multiplesamples <- div_profile(otu_df_hilldiv)
div_profile_plot(profile.multiplesamples)

#Multiple groups (gamma diversity)
profile.multiplegroups <- div_profile(bat.diet.otutable,hierarchy=bat.diet.hierarchy,level="gamma")
div_profile_plot(profile.multiplegroups)

###Menu Multiple groups (gamma diversity)----
profile.multiplegroups <- div_profile(otu_df_hilldiv,hierarchy=menu_sample_site_hierarchy,tree = menu_tree_hilldiv_ultra,level="gamma")

div_profile_plot(profile.multiplegroups)



##div part ----
data(bat.diet.otutable)
data(bat.diet.tree)
data(bat.diet.hierarchy)

#Two level examples (L1=sample (alpha diversity), L2=whole system (gamma diversity))
div_part(bat.diet.otutable,qvalue=1)
div_part(bat.diet.otutable,qvalue=0,tree=bat.diet.tree)

#Three-level example (L1=sample, L2=species, L3=whole system)
div_part(bat.diet.otutable,qvalue=0,hierarchy=bat.diet.hierarchy)




###Menu diversity partitioning---
div_part(otu_df_hilldiv,qvalue=1,hierarchy = menu_sample_site_hierarchy,tree = menu_tree_hilldiv_ultra)

#Faith's PD by site
# $Hierarchical_levels
# [1] 3
# 
# $Order_diversity
# [1] 0
# 
# $Hill_numbers
# L1        L2        L3 
# 7.335224 30.666815 30.666815 
# 
# $Sample_size
# N1 N2 N3 
# 59  7  1 
# 
# $Beta
# B1_2    B2_3 
# 4.18076 1.00000 

#Allen's H
# $Hierarchical_levels
# [1] 3
# 
# $Order_diversity
# [1] 1
# 
# $Hill_numbers
# L1       L2       L3 
# 2.219375 4.093147 4.210879 
# 
# $Sample_size
# N1 N2 N3 
# 59  7  1 
# 
# $Beta
# B1_2     B2_3 
# 1.844279 1.028763 

#Beta diversity is higher between samples within a site than between sites



#Corrplot of Water_temp and primary productivity vs diversity measures from hill_div----

library(corrplot)

##Combine MiFish and BerryCrust diversity dataframes by site----
Fish_and_Crustacean_hilldiv_bysite<-full_join(MiFish_diversity_averagesample_bysite,BerryCrust_diversity_averagesamplebysite.df)

# BerryCrust_diversity_averagesamplebysite.df1<-column_to_rownames( BerryCrust_diversity_averagesamplebysite.df, var="Site_Name")

Fish_and_Crustacean_hilldiv_bysite.df<-column_to_rownames(Fish_and_Crustacean_hilldiv_bysite, var="Site_Name")

# BerryCrust_diversity_averagesamplebysite.df2<-BerryCrust_diversity_averagesamplebysite.df1%>%
#   select(BerryCrust_mean_FaithsPD,BerryCrust_mean_AllensH,BerryCrust_mean_Richness,BerryCrust_mean_Shannonexp,mean_site_water_temp,var_site_water_temp,mean_primary_prod,var_primary_prod,hammerhead_site,carcharhinus_site)

Fish_and_Crustacean_hilldiv_bysite.df1<-Fish_and_Crustacean_hilldiv_bysite.df%>%
     dplyr::select(BerryCrust_mean_FaithsPD,BerryCrust_mean_AllensH,BerryCrust_mean_Richness,BerryCrust_mean_Shannonexp,MiFish_mean_FaithsPD,MiFish_mean_AllensH,MiFish_mean_Richness,MiFish_mean_Shannonexp, mean_site_water_temp,var_site_water_temp,mean_primary_prod,var_primary_prod,mean_nitrate,var_nitrate,hammerhead_site,carcharhinus_site)%>%
  rename(MiFish_mean_Shannon=MiFish_mean_Shannonexp,
         BerryCrust_mean_Shannon=BerryCrust_mean_Shannonexp
         )


BC_cor<-cor(Fish_and_Crustacean_hilldiv_bysite.df1)
head(round(BC_cor,2))


corrplot(BC_cor, method="circle",type="lower")

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
p.mat <- cor.mtest(BC_cor)
head(p.mat[, 1:5])

# mark out the ones with p-values greater than 0.05
corrplot(BC_cor, type="lower", 
         p.mat = p.mat, sig.level = 0.05)

# Leave blank on no significant coefficient
corrplot(BC_cor, type="lower", 
         p.mat = p.mat, sig.level = 0.05, insig = "blank")

#Best correlation plot----
corrplot(BC_cor, type="lower",# order = "hclust",#addCoefasPercent = TRUE, 
         addCoef.col = "black",number.cex = .5, # Add coefficient of correlation
         tl.col="black", tl.srt=45,tl.cex = .5, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=TRUE,
         #title = "Correlation of Crustacean and Fish Diversity Measures with Temperature, Primary Productivity, and Shark Nursery Status by Site" 
)

view(p.mat)
write.csv(p.mat,"../11_Vegan/p.mat_sitewise_average_corrplot.csv")


#try phyloraster----
library(phyloraster)
library(terra)
library(ape)
library(phylobase)

data <- load.data.rosauer()
head(data$presab)


Menu1_pa.ps<-phyloseq_standardize_otu_abundance(Menu1.ps, method = "pa")

menu_pa<-as.data.frame(t(cbind(otu_table(Menu1_pa.ps))))


menu_pa<-rownames_to_column(menu_pa, var = "sample")
  

menu_loc<-as.data.frame(cbind(sample_data(Menu1_pa.ps)))

menu_loc<-menu_loc%>%
  select(DD_long,DD_lat)

menu_loc<-rownames_to_column(menu_loc, var = "sample")

menu_loc_pa<-full_join(menu_loc,menu_pa)
menu_loc_pa<-menu_loc_pa%>%
  select(-sample)

menu_tree<-Menu1_pa.ps@phy_tree
class(menu_tree)
plot(menu_tree, cex=0.3)

r1<-df2rast(x= menu_loc_pa, CRS = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#Error: [raster,matrix(xyz)] x cell sizes are not regular

data$tree
class(data$tree)

plot(data$tree, cex = 0.65)

data <- load.data.rosauer()
r <- df2rast(x = data$presab, 
             CRS = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

class(r)

plot(r)



shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp", 
                               package = "phyloraster"))

colors <- rainbow(length(unique(shp$BINOMIAL)),
                  alpha = 0.5)
position <- match(shp$BINOMIAL,
                  unique(shp$BINOMIAL))
colors <- colors[position]
plot(shp, col = colors, lty = 0,
     main = "Spatial polygons")
library(maps)
maps::map(add = TRUE)

r2 <- shp2rast(shp, sps.col = "BINOMIAL", ymask = FALSE, background = 0, 
               resolution = 0.5)
r2

plot(r2[[9]])

library(terra)

shp <- terra::vect(system.file("extdata", "shps_iucn_spps_rosauer.shp",
                               package="phyloraster"))

# create a polygon to use as mask with an extent
e <- terra::ext(113, 123, -43.64, -33.90)
p <- terra::as.polygons(e, crs="")
# cut by the total extension of the polygons
coun.crop <- terra::crop(p, 
                         terra::ext(shp)) 
coun.rast <- terra::rasterize(coun.crop,
                              terra::rast(terra::ext(shp), resolution = 0.5))

# rasterizing with the mask of the polygon
shp.t <- shp2rast(shp, y = coun.rast, sps.col = "BINOMIAL", ymask = TRUE)
plot(shp.t[[1]], col = c("grey", "green"))

data <- load.data.rosauer()
names(data$raster) == data$tree$tip.label

ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", 
                                   package = "phyloraster"))
dataprep <- phylo.pres(x = ras, tree = tree)

names(dataprep$x) == tree$tip.label

ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
sr <- rast.sr(x = ras)
sr

plot(sr, main = "Species richness")

ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
wer <- rast.we(x = ras)
wer

wer$WE

plot(wer$WE, main ="Weighted Endemism")

ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", 
                                   package = "phyloraster"))
dataprep <- phylo.pres(x = ras, tree = tree, pruning = "tree")

pdr <- rast.pd(x = dataprep$x, edge.path = dataprep$edge.path, 
               branch.length = dataprep$branch.length)

plot(pdr$PD, main = "Phylogenetic diversity")

ras <- terra::rast(system.file("extdata", "rast.presab.tif", 
                               package = "phyloraster"))
tree <- ape::read.tree(system.file("extdata", "tree.nex", 
                                   package = "phyloraster"))
per <- rast.pe(x = dataprep$x, tree)
per

plot(per$PE, main = "Phylogenetic Endemism")

x <- terra::rast(system.file("extdata", "rast.presab.tif", package="phyloraster"))
# phylogenetic tree
tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))
data <- phylo.pres(x, tree)
ed <- rast.ed(data$x, tree)
ed

terra::plot(ed, main = "Evolutionary Distinctiveness")

tree <- ape::read.tree(system.file("extdata", "tree.nex", package="phyloraster"))

ed <- phyloraster::species.ed(tree)
head(ed)


