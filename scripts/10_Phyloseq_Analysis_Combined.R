#Eldridge Wisely
#Date: 4-30-24
#Do analysis and visualization of combined BerryCrust and MiFish phyloseq objects 

library("ggplot2")      # graphics
#library("readxl")       # necessary to import the data from Excel file
#library("tibble")       # Needed for converting column to row names
# Load requested package for plotting
#library("reshape2")
library("tidyverse")
library("dplyr")        # filter and reformat data frames
library(phyloseq)
#library(adegenet)
#library(ape)
#library(phangorn)
library(ggsci)
library(metagMisc)
library(ggiraph)
library(microViz)

'%ni%' <- Negate("%in%")

#set variables
Primer<-"Menu_combined1" #Menu_combined =includes empty PCRs, Menu_combined1 =empty PCRs have been removed


#Import Phyloseq object----
Menu.ps<-readRDS(file= paste("../09_Metabar_to_Phyloseq/",Primer,"_Phyloseq.rds", sep = ""))


#make a network graph----
ig <- make_network(Menu.ps, max.dist=0.9)
plot_network(ig, Menu.ps,label = "Microhabitat")

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

#Take the phyloseq tables into dataframes and visualize

otu_df<-cbind(data.frame(otu_table(Menu_species_merge1.pa.ps)))
otu_df$otu<-rownames(otu_df)

otu_df_long<-otu_df%>%
  pivot_longer(cols = -otu,names_to = "sample",values_to = "presence")
otu_df_long$sample <-gsub(pattern = ".",
                          "-",
                          x=otu_df_long$sample,
                          fixed = TRUE)



sample_df<-cbind(data.frame(sample_data(Menu_species_merge1.pa.ps)))
sample_df<-rownames_to_column(sample_df, var = "sample")
sample_site_key<-sample_df%>%
  select(sample,Site_Name)

taxa_df<-cbind(data.frame(tax_table(Menu_species_merge1.pa.ps)))
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





