library(vegan)
library(tidyverse)
library(hilldiv)
library(vegan)

#dbRDA----

## prep taxonomic (taxa-merged) species-named eDNA index for hilldiv----
Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps<-readRDS("../10_Phyloseq/Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps.RDS")

# #make otu table and ultrametric tree for Menu_taxmerge_eDNAindex ASVs(taxa)
# Menu_taxmerge_eDNAindex_otudf<-cbind(data.frame(otu_table(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))
# #Make ultrametric tree for otu object
# Menu_taxmerge_eDNAindex_tree<-Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps@phy_tree
# Menu_taxmerge_eDNAindex_ultra_tree<-chronos(Menu_taxmerge_eDNAindex_tree,lambda = 1,model = "relaxed")
# write.tree(Menu_taxmerge_eDNAindex_ultra_tree, 'Menu_taxmerge_eDNAindex_ultra.tre')
# Menu_taxmerge_eDNAindex_ultra_tree <- ape::read.tree('Menu_taxmerge_eDNAindex_ultra.tre')
# 
# #make rownames of otu table match tip names of the tree
# namechange<-rownames(Menu_taxmerge_eDNAindex_otudf)
# namechange<-gsub(pattern = " ",replacement = "_",namechange)
# rownames(Menu_taxmerge_eDNAindex_otudf)<-namechange

#make envvars
Menu_taxmerge_eDNAindex_envvars<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))


#use hilldiv distance matrix on eDNA transformed otu table to do a dbRDA----
#these objects came from hilldiv in script 11!
#Then proceed with dbrda function from vegan
dbRDA_Faith <- dbrda(Menu_taxmerged_Faiths_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars,na.action = na.omit,sqrt.dist=FALSE) #dbrda from vegan
summary(dbRDA_Faith)
#Jaccard
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total           6.130     1.0000
# Conditioned     1.235     0.2015
# Constrained     1.529     0.2494
# Unconstrained   3.366     0.5491

#Sorensen
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total          2.9458     1.0000
# Conditioned    0.6018     0.2043
# Constrained    0.8074     0.2741
# Unconstrained  1.5366     0.5216

anova(dbRDA_Faith, by ="margin")
#Jaccard
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


#Sorensen
# Model: dbrda(formula = Menu_taxmerged_Faiths_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)   
# Microhabitat           2  0.11669 1.6327  0.162   
# scaled_Water_temp      1  0.07463 2.0884  0.149   
# Tide                   1  0.04670 1.3069  0.263   
# scaled_nitrate         1  0.24573 6.8763  0.004 **
#   hammerheads_detected   1  0.14598 4.0850  0.023 * 
#   carcharhinus_detected  1  0.01991 0.5570  0.566   
# scaled_primaryprod     1  0.20817 5.8253  0.006 **
#   Residual              43  1.53664                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

plot(dbRDA_Faith)

Shark_Menu_taxa_merge_species_named_traits_temp.eDNAindex.ps<- Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps %>% ps_filter(Water_temp !="NA")


speDNAindex<-t(as.data.frame(cbind(otu_table(Shark_Menu_taxa_merge_species_named_traits_temp.eDNAindex.ps))))
sppscores(dbRDA_Faith)<-decostand(speDNAindex,"range")


plot(dbRDA_Faith)


dbRDA_Richness <- dbrda(Menu_taxmerged_Richness_Hill_Sor_dissim ~Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan
summary(dbRDA_Richness) #Richness
#Jaccard
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total          18.455     1.0000
# Conditioned     2.611     0.1415
# Constrained     3.801     0.2060
# Unconstrained  12.043     0.6525

#Sorensen
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total          13.897     1.0000
# Conditioned     2.140     0.1540
# Constrained     3.318     0.2387
# Unconstrained   8.439     0.6072
anova(dbRDA_Richness, by ="margin")
#Jaccard
# Model: dbrda(formula = Menu_taxmerged_Richness_Hill_Jaccard_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)   
# Microhabitat           2   0.9171 1.6373  0.005 **
#   scaled_Water_temp      1   0.4612 1.6469  0.027 * 
#   Tide                   1   0.3479 1.2423  0.177   
# scaled_nitrate         1   0.5414 1.9333  0.005 **
#   hammerheads_detected   1   0.3657 1.3057  0.140   
# carcharhinus_detected  1   0.4200 1.4996  0.050 * 
#   scaled_primaryprod     1   0.4186 1.4945  0.054 . 
# Residual              43  12.0426                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Sorensen
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

dbRDA_Shannon<-dbrda(Menu_taxmerged_Shannons_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(dbRDA_Shannon)
#Jaccard
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total          21.769     1.0000
# Conditioned     3.606     0.1656
# Constrained     4.564     0.2097
# Unconstrained  13.599     0.6247

#Sorensen
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total          18.033     1.0000
# Conditioned     3.525     0.1955
# Constrained     4.492     0.2491
# Unconstrained  10.016     0.5554

anova(dbRDA_Shannon, by ="margin") 
#Jaccard
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

#Sorensen

# Model: dbrda(formula = Menu_taxmerged_Shannons_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)    
# Microhabitat           2   1.2451 2.6727  0.001 ***
#   scaled_Water_temp      1   0.5537 2.3770  0.005 ** 
#   Tide                   1   0.4587 1.9692  0.018 *  
#   scaled_nitrate         1   0.6329 2.7170  0.003 ** 
#   hammerheads_detected   1   0.3343 1.4350  0.158    
# carcharhinus_detected  1   0.6950 2.9835  0.001 ***
#   scaled_primaryprod     1   0.4205 1.8052  0.041 *  
#   Residual              43  10.0161                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
plot(dbRDA_Shannon)

dbRDA_Allen<-dbrda(Menu_taxmerged_Allens_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), Menu_taxmerge_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(dbRDA_Allen)
#Sorensen
# dbrda(formula = Menu_taxmerged_Allens_Hill_Sor_dissim ~ Microhabitat +      scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected +      carcharhinus_detected + scaled_primaryprod + Condition(Site_Name),      data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE,      na.action = na.omit) 
# 
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total         0.20699    1.00000
# Conditioned   0.01005    0.04856
# Constrained   0.08688    0.41971
# Unconstrained 0.11006    0.53173
anova(dbRDA_Allen, by ="margin")  

#Sorensen
# Model: dbrda(formula = Menu_taxmerged_Allens_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = Menu_taxmerge_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs       F Pr(>F)   
# Microhabitat           2 0.027812  5.4329  0.008 **
#   scaled_Water_temp      1 0.014124  5.5180  0.029 * 
#   Tide                   1 0.004045  1.5802  0.247   
# scaled_nitrate         1 0.025723 10.0497  0.003 **
#   hammerheads_detected   1 0.002634  1.0290  0.313   
# carcharhinus_detected  1 0.005472  2.1377  0.158   
# scaled_primaryprod     1 0.013927  5.4409  0.029 * 
#   Residual              43 0.110063                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



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

## Remove Hammerheads and see if that makes a difference----


HH_rm_dbRDA_Faith <- dbrda(Menu_taxmerged_rmhh_Faiths_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), HH_rm_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)#,na.action = na.omit) #dbrda from vegan


summary(HH_rm_dbRDA_Faith)
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total           6.407     1.0000
# Conditioned     1.198     0.1870
# Constrained     1.618     0.2526
# Unconstrained   3.591     0.5604

anova(HH_rm_dbRDA_Faith,by ="margin")

# Model: dbrda(formula = Menu_taxmerged_rmhh_Faiths_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = HH_rm_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)   
# Microhabitat           2  0.15968 1.3723  0.288   
# scaled_Water_temp      1  0.13967 2.4006  0.087 . 
# Tide                   1  0.08007 1.3763  0.269   
# scaled_nitrate         1  0.32669 5.6151  0.006 **
#   hammerheads_detected   1  0.16595 2.8524  0.049 * 
#   carcharhinus_detected  1  0.03790 0.6514  0.571   
# scaled_primaryprod     1  0.24637 4.2346  0.009 **
#   Residual              43  2.50175                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

HH_rm_dbRDA_Allen <- dbrda(Menu_taxmerged_rmhh_Allens_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), HH_rm_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)#,na.action = na.omit) #dbrda from vegan


summary(HH_rm_dbRDA_Allen)
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total         0.724629   1.000000
# Conditioned   0.002224   0.003069
# Constrained   0.295480   0.407767
# Unconstrained 0.426925   0.589164

anova(HH_rm_dbRDA_Allen,by ="margin")
# Model: dbrda(formula = Menu_taxmerged_rmhh_Allens_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = HH_rm_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs       F Pr(>F)  
# Microhabitat           2  0.08960  4.5123  0.031 *
#   scaled_Water_temp      1  0.05618  5.6588  0.036 *
#   Tide                   1  0.01119  1.1267  0.323  
# scaled_nitrate         1  0.10034 10.1058  0.011 *
#   hammerheads_detected   1  0.00227  0.2288  0.644  
# carcharhinus_detected  1  0.01709  1.7212  0.203  
# scaled_primaryprod     1  0.05337  5.3752  0.043 *
#   Residual              43  0.42693                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

HH_rm_dbRDA_Shannon <- dbrda(Menu_taxmerged_rmhh_Shannons_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), HH_rm_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(HH_rm_dbRDA_Shannon)

# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total          18.004     1.0000
# Conditioned     3.520     0.1955
# Constrained     4.461     0.2478
# Unconstrained  10.023     0.5567

anova(HH_rm_dbRDA_Shannon, by = "margin")

# Model: dbrda(formula = Menu_taxmerged_rmhh_Shannons_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = HH_rm_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)    
# Microhabitat           2   1.2460 2.6728  0.001 ***
#   scaled_Water_temp      1   0.5555 2.3830  0.003 ** 
#   Tide                   1   0.4568 1.9598  0.021 *  
#   scaled_nitrate         1   0.6342 2.7209  0.003 ** 
#   hammerheads_detected   1   0.3051 1.3090  0.228    
# carcharhinus_detected  1   0.6952 2.9825  0.002 ** 
#   scaled_primaryprod     1   0.4215 1.8081  0.055 .  
# Residual              43  10.0232                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


HH_rm_dbRDA_Richness <- dbrda(Menu_taxmerged_rmhh_Richness_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate +hammerheads_detected+carcharhinus_detected+ scaled_primaryprod+Condition(Site_Name), HH_rm_eDNAindex_envvars, sqrt.dist=FALSE, na.action=na.omit)

summary(HH_rm_dbRDA_Richness)
# Partitioning of squared Unknown distance:
#   Inertia Proportion
# Total          13.863     1.0000
# Conditioned     2.139     0.1543
# Constrained     3.255     0.2348
# Unconstrained   8.469     0.6109

anova(HH_rm_dbRDA_Richness,by = "margin")
# Model: dbrda(formula = Menu_taxmerged_rmhh_Richness_Hill_Sor_dissim ~ Microhabitat + scaled_Water_temp + Tide + scaled_nitrate + hammerheads_detected + carcharhinus_detected + scaled_primaryprod + Condition(Site_Name), data = HH_rm_eDNAindex_envvars, sqrt.dist = FALSE, na.action = na.omit)
# Df SumOfSqs      F Pr(>F)   
# Microhabitat           2   0.7833 1.9885  0.005 **
#   scaled_Water_temp      1   0.4093 2.0783  0.033 * 
#   Tide                   1   0.2425 1.2314  0.258   
# scaled_nitrate         1   0.5000 2.5388  0.008 **
#   hammerheads_detected   1   0.2230 1.1324  0.313   
# carcharhinus_detected  1   0.3862 1.9608  0.033 * 
#   scaled_primaryprod     1   0.3449 1.7510  0.057 . 
# Residual              43   8.4691                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Check for differences between Carcharhinus sites and non-carcharhinus sites.----
#First remove Carcharhinus reads so they don't confound the calculation.
Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps <- Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps %>% subset_taxa(genus != "Carcharhinus")

landes_dist_matrix_rm_carch<-phyloseq::distance(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps, method = "l")

set.seed(2)
pairwise.adonis(landes_dist_matrix_rm_carch, phyloseq::sample_data(Shark_Menu_taxa_merge_species_named_traits_rm_carch.eDNAindex.ps)$carcharhinus_detected)

#dbRDAs on Fish and Crustaceans separately----
#MiFish hilldiv----
#get just the MiFish results phyloseq
MiFish.ps<-readRDS("../09_Metabar_to_Phyloseq/MiFish_Phyloseq.rds")
#filter out passive mid-bay samples because they didn't seem to work well
MiFish1.ps<- ps_filter(MiFish.ps, Microhabitat != "passive mid bay")

MiFish1.ps

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

MiFish_sample_data<-cbind(data.frame(sample_data(MiFish_taxa_merge_eDNAindex.ps)))
MiFish_samples<-rownames(MiFish_sample_data)

Menu_sample_data<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))%>%
  rownames_to_column(var="sample_id")
MiFish_sample_data<-Menu_sample_data%>%
  filter(sample_id %in% MiFish_samples)

length(MiFish_samples)
nrow(Menu_sample_data)
nrow(MiFish_sample_data)

MiFish_sample_site_hierarchy<-MiFish_sample_data%>%
  select(sample_id,Site_Name)

##do Hill diversity distance matrices for MiFish----

##Phylogenetic Hill q=0 (Faith's PD) distance matrix with taxa merged object ----

MiFish_taxmerged_Faiths_Hilldistmat<-pair_dis(countable = MiFish_taxa_df_eDNAindex,qvalue = 0,hierarchy=MiFish_sample_site_hierarchy,level = 1,tree = MiFish_taxamerge_ultra_tree)

##Phylogenetic Hill q=1 (Allen's H) distance matrix with taxa merged object ----

MiFish_taxmerged_Allens_Hilldistmat<-pair_dis(countable = MiFish_taxa_df_eDNAindex,qvalue = 1,hierarchy=MiFish_sample_site_hierarchy,level = 1,tree = MiFish_taxamerge_ultra_tree)

##Non-Phylogenetic Hill q=0 (Richness) distance matrix with taxa merged object ----

MiFish_taxmerged_Richness_Hilldistmat<-pair_dis(countable = MiFish_taxa_df_eDNAindex,qvalue = 0,hierarchy=MiFish_sample_site_hierarchy,level = 1)

##Non-Phylogenetic Hill q=1 (exponent of Shannon's) distance matrix with taxa merged object ----

MiFish_taxmerged_Shannons_Hilldistmat<-pair_dis(countable = MiFish_taxa_df_eDNAindex,qvalue = 1,hierarchy=MiFish_sample_site_hierarchy,level = 1)

#Get Jaccard-type turnover complement (with respect to gamma) as a symmetrical distance matrix for each----
MiFish_taxmerged_Faiths_Hill_Jaccard_dissim<-as.dist(MiFish_taxmerged_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

MiFish_taxmerged_Allens_Hill_Jaccard_dissim<-as.dist(MiFish_taxmerged_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

MiFish_taxmerged_Richness_Hill_Jaccard_dissim<-as.dist(MiFish_taxmerged_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

MiFish_taxmerged_Shannons_Hill_Jaccard_dissim<-as.dist(MiFish_taxmerged_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)




#BerryCrust hilldiv----
#get just the BerryCrust results phyloseq
BerryCrust.ps<-readRDS("../09_Metabar_to_Phyloseq/BerryCrust_Phyloseq.rds")
#filter out passive mid-bay samples because they didn't seem to work well
BerryCrust1.ps<- ps_filter(BerryCrust.ps, Microhabitat != "passive mid bay")

BerryCrust1.ps

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

BerryCrust_taxa_df_eDNAindex<-cbind(data.frame(otu_table(BerryCrust_taxa_merge_eDNAindex.ps)))
#Make ultrametric tree for taxa merged object 
BerryCrust_taxamerge_tree<-BerryCrust_taxa_merge_eDNAindex.ps@phy_tree
BerryCrust_taxamerge_ultra_tree<-chronos(BerryCrust_taxamerge_tree,lambda = 1,model = "relaxed")
write.tree(BerryCrust_taxamerge_ultra_tree, 'BerryCrust_taxamerge_ultra.tre')
BerryCrust_taxamerge_ultra_tree <- ape::read.tree('BerryCrust_taxamerge_ultra.tre')
#make taxa_df names match tips of the tree
mfnames<-rownames(BerryCrust_taxa_df_eDNAindex)
mfnames<-gsub(pattern = " ","_",mfnames)
rownames(BerryCrust_taxa_df_eDNAindex)<-mfnames

BerryCrust_sample_data<-cbind(data.frame(sample_data(BerryCrust_taxa_merge_eDNAindex.ps)))
BerryCrust_samples<-rownames(BerryCrust_sample_data)

Menu_sample_data<-cbind(data.frame(sample_data(Shark_Menu_taxa_merge_species_named_traits.eDNAindex.ps)))%>%
  rownames_to_column(var="sample_id")
BerryCrust_sample_data<-Menu_sample_data%>%
  filter(sample_id %in% BerryCrust_samples)

length(BerryCrust_samples)
nrow(Menu_sample_data)
nrow(BerryCrust_sample_data)

BerryCrust_sample_site_hierarchy<-BerryCrust_sample_data%>%
  select(sample_id,Site_Name)

##do Hill diversity distance matrices for BerryCrust----

##Phylogenetic Hill q=0 (Faith's PD) distance matrix with taxa merged object ----

BerryCrust_taxmerged_Faiths_Hilldistmat<-pair_dis(countable = BerryCrust_taxa_df_eDNAindex,qvalue = 0,hierarchy=BerryCrust_sample_site_hierarchy,level = 1,tree = BerryCrust_taxamerge_ultra_tree)

##Phylogenetic Hill q=1 (Allen's H) distance matrix with taxa merged object ----

BerryCrust_taxmerged_Allens_Hilldistmat<-pair_dis(countable = BerryCrust_taxa_df_eDNAindex,qvalue = 1,hierarchy=BerryCrust_sample_site_hierarchy,level = 1,tree = BerryCrust_taxamerge_ultra_tree)

##Non-Phylogenetic Hill q=0 (Richness) distance matrix with taxa merged object ----

BerryCrust_taxmerged_Richness_Hilldistmat<-pair_dis(countable = BerryCrust_taxa_df_eDNAindex,qvalue = 0,hierarchy=BerryCrust_sample_site_hierarchy,level = 1)

##Non-Phylogenetic Hill q=1 (exponent of Shannon's) distance matrix with taxa merged object ----

BerryCrust_taxmerged_Shannons_Hilldistmat<-pair_dis(countable = BerryCrust_taxa_df_eDNAindex,qvalue = 1,hierarchy=BerryCrust_sample_site_hierarchy,level = 1)

#Get Jaccard-type turnover complement (with respect to gamma) as a symmetrical distance matrix for each----
BerryCrust_taxmerged_Faiths_Hill_Jaccard_dissim<-as.dist(BerryCrust_taxmerged_Faiths_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

BerryCrust_taxmerged_Allens_Hill_Jaccard_dissim<-as.dist(BerryCrust_taxmerged_Allens_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

BerryCrust_taxmerged_Richness_Hill_Jaccard_dissim<-as.dist(BerryCrust_taxmerged_Richness_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)

BerryCrust_taxmerged_Shannons_Hill_Jaccard_dissim<-as.dist(BerryCrust_taxmerged_Shannons_Hilldistmat$L1_SqN,diag = TRUE, upper = TRUE)


#dbRDAs for MiFish----

##Faiths
dbRDA_MiFish_Faiths <- dbrda(MiFish_taxmerged_Faiths_Hill_Jaccard_dissim ~scaled_Water_temp + Microhabitat + Tide + scaled_nitrate +scaled_primaryprod+hammerheads_detected+carcharhinus_detected+ Condition(Site_Name), MiFish_sample_data,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan

summary(dbRDA_MiFish_Faiths) #Faiths

anova(dbRDA_MiFish_Faiths, by ="margin")

plot(dbRDA_MiFish_Faiths)

##Allens
dbRDA_MiFish_Allens <- dbrda(MiFish_taxmerged_Allens_Hill_Jaccard_dissim ~scaled_Water_temp + Microhabitat + Tide + scaled_nitrate +scaled_primaryprod+hammerheads_detected+carcharhinus_detected+ Condition(Site_Name), MiFish_sample_data,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan

summary(dbRDA_MiFish_Allens) #Allens

anova(dbRDA_MiFish_Allens, by ="margin")

plot(dbRDA_MiFish_Allens)

#Richness
dbRDA_MiFish_Richness <- dbrda(MiFish_taxmerged_Richness_Hill_Jaccard_dissim ~scaled_Water_temp + Microhabitat + Tide + scaled_nitrate +scaled_primaryprod+hammerheads_detected+carcharhinus_detected+ Condition(Site_Name), MiFish_sample_data,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan

summary(dbRDA_MiFish_Richness) #Richness

anova(dbRDA_MiFish_Richness, by ="margin")

plot(dbRDA_MiFish_Richness)

#Shannons
dbRDA_MiFish_Shannons <- dbrda(MiFish_taxmerged_Shannons_Hill_Jaccard_dissim ~scaled_Water_temp + Microhabitat + Tide + scaled_nitrate +scaled_primaryprod+hammerheads_detected+carcharhinus_detected+ Condition(Site_Name), MiFish_sample_data,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan

summary(dbRDA_MiFish_Shannons) #Shannons

anova(dbRDA_MiFish_Shannons, by ="margin")

plot(dbRDA_MiFish_Shannons)


#dbRDAs for BerryCrust----

##Faiths
dbRDA_BerryCrust_Faiths <- dbrda(BerryCrust_taxmerged_Faiths_Hill_Jaccard_dissim ~scaled_Water_temp + Microhabitat + Tide + scaled_nitrate +scaled_primaryprod+hammerheads_detected+carcharhinus_detected+ Condition(Site_Name), BerryCrust_sample_data,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan

summary(dbRDA_BerryCrust_Faiths) #Faiths

anova(dbRDA_BerryCrust_Faiths, by ="margin")

plot(dbRDA_BerryCrust_Faiths)

##Allens
dbRDA_BerryCrust_Allens <- dbrda(BerryCrust_taxmerged_Allens_Hill_Jaccard_dissim ~scaled_Water_temp + Microhabitat + Tide + scaled_nitrate +scaled_primaryprod+hammerheads_detected+carcharhinus_detected+ Condition(Site_Name), BerryCrust_sample_data,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan

summary(dbRDA_BerryCrust_Allens) #Allens

anova(dbRDA_BerryCrust_Allens, by ="margin")

plot(dbRDA_BerryCrust_Allens)

#Richness
dbRDA_BerryCrust_Richness <- dbrda(BerryCrust_taxmerged_Richness_Hill_Jaccard_dissim ~scaled_Water_temp + Microhabitat + Tide + scaled_nitrate +scaled_primaryprod+hammerheads_detected+carcharhinus_detected+ Condition(Site_Name), BerryCrust_sample_data,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan

summary(dbRDA_BerryCrust_Richness) #Richness

anova(dbRDA_BerryCrust_Richness, by ="margin")

plot(dbRDA_BerryCrust_Richness)

#Shannons
dbRDA_BerryCrust_Shannons <- dbrda(BerryCrust_taxmerged_Shannons_Hill_Jaccard_dissim ~scaled_Water_temp + Microhabitat + Tide + scaled_nitrate +scaled_primaryprod+hammerheads_detected+carcharhinus_detected+ Condition(Site_Name), BerryCrust_sample_data,sqrt.dist = FALSE,na.action = na.omit) #dbrda from vegan

summary(dbRDA_BerryCrust_Shannons) #Shannons

anova(dbRDA_BerryCrust_Shannons, by ="margin")

plot(dbRDA_BerryCrust_Shannons)






