#Eldridge Wisely
#Date: 4-30-24
#Do analysis and visualization of BerryCrust and MiFish phyloseq objects separately

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

'%ni%' <- Negate("%in%")

#set variables
Primer<-"MiFish" #choices: MiFish, BerryCrust, Menu_Combined


#Import Metabarlist----
Menu.ps<-readRDS(file= paste("../09_Metabar_to_Phyloseq/",Primer,"_Phyloseq.rds", sep = ""))


#Make phylogenetic tree visualization of cleaned ID'ed ASVs
head(phy_tree(Menu.ps)$node.label, 10)

ntaxa(Menu.ps)

plot_tree(Menu.ps, nodelabf=nodeplotboot(), ladderize="left", color="Microhabitat")+
  scale_color_locuszoom()

plot_tree(Menu.ps, nodelabf=nodeplotboot(), color="Microhabitat", label.tips="species", ladderize="left")+
  scale_color_locuszoom()
ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Phylogenetic_Tree_Species_by_Microhabitat.jpg"))


plot_tree(Menu.ps, nodelabf=nodeplotboot(), ladderize="left", color="order", shape="Microhabitat")+ scale_color_igv()
ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Phylogenetic_Tree_Order_by_Microhabitat.jpg"))


plot_tree(Menu.ps, nodelabf=nodeplotboot(), color="Microhabitat", label.tips="species", shape="class",ladderize="left")+
    scale_color_lancet()
ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Phylogenetic_Tree_Classes_by_Microhabitat.jpg"))


#take just the top 50 ASVs
top50 = names(sort(taxa_sums(Menu.ps), decreasing = TRUE)[1:50])
top50_Menu.ps = prune_taxa(top50, Menu.ps)

plot(phy_tree(top50_Menu.ps), show.node.label = TRUE)

plot_tree(top50_Menu.ps, nodelabf=nodeplotboot(), color="order", ladderize="left",label.tips="species") + coord_polar(theta="y")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Circular_Phylogenetic_Tree_Top50ASVs.jpg"))

#make a network graph----
ig <- make_network(Menu.ps, max.dist=0.9)
plot_network(ig, Menu.ps,label = "Site_Name")


#Plot richness----
plot_richness(Menu.ps, x="Site_Name","Microhabitat", measures=c("Shannon", "Simpson"), color="Site_Name")+
  scale_color_ucscgb()


ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Alpha_Diversity_untransformed_ASVs.jpg"))


# Transform data to proportions as appropriate for Bray-Curtis distances proportion of total reads per sample
Menu.ps.prop <- transform_sample_counts(Menu.ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(Menu.ps.prop, method="NMDS", distance="bray")


plot_ordination(Menu.ps.prop, ord.nmds.bray, color="Microhabitat"
                , title="Bray NMDS")

#The above seems like it's better for the non_aggregated dataset
##View the Top 50 most abundant ASVs in the dataset----

#take just the top 50 ASVs
top50 = names(sort(taxa_sums(Menu.ps), decreasing = TRUE)[1:50])
top50_Menu.ps = prune_taxa(top50, Menu.ps)
plot_bar(top50_Menu.ps, x="Site_Name",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("Read Abundance")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Top50_ASVs_by_site_and_genus.jpg"))

top50 = names(sort(taxa_sums(Menu.ps.prop), decreasing = TRUE)[1:50])
top50_Menu.ps.prop = prune_taxa(top50, Menu.ps.prop)
plot_bar(top50_Menu.ps.prop, x="Site_Name",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("Read Abundance")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Top50_ASVs_by_site_and_genus_proportions.jpg"))


#Merge_species but keep taxa that aren't ID'ed to species
Menu_taxa_merge.ps <- tax_glom(Menu.ps, taxrank = "species", NArm = FALSE)

Menu_species_merge1.ps<-tax_glom(Menu.ps, taxrank = "species", NArm=TRUE)


#convert to presence absence


Menu_species_merge1.pa.ps<-phyloseq_standardize_otu_abundance(Menu_species_merge1.ps, method = "pa")
#plot frequency of occurrence(FOO) by genus and site
plot_bar(Menu_species_merge1.pa.ps, x="Site_Name",fill="species")+ 
  geom_bar(aes(color=species, fill=species), stat="identity", position="stack")+
  ylab("Frequency of Occurrence")

#find the top 20 most frequently occurring species

top20occ = names(sort(taxa_sums(Menu_species_merge1.pa.ps), decreasing = TRUE)[1:20])
top20occ_Menu_species_merge1.pa.ps = prune_taxa(top20occ, Menu_species_merge1.pa.ps)

plot_bar(top20occ_Menu_species_merge1.pa.ps, x="Site_Name",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("Frequency of Occurrence")

# 
# #plot percent of occurrence (POO) by genus and site
# 
# POO<-phyloseq_otu_occurrence(
#   top20occ_Menu_species_merge1.pa.ps,
#   variable = "Site_Name",
#   taxa_frequency = "percentage",
#   drop_zeroes = FALSE,
#   justdf = FALSE,
#   long = TRUE)
# 
# 
# 
# plot_bar(POO, x="Site_Name",fill="genus")+ 
#   geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
#   ylab("Percent of Occurrence")
# 



##50 rarest ASVs ----
lowest50 = names(sort(taxa_sums(Menu.ps.prop), decreasing = FALSE)[1:50])
lowest50_Menu.ps.prop = prune_taxa(lowest50, Menu.ps.prop)
plot_bar(lowest50_Menu.ps.prop, x="Site_Name",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("Read Abundance")

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Lowest50_ASVs_by_site_and_genus_proportions.jpg"))


###############

library(RColorBrewer)
#top20 ASVs in the dataset
top20 = names(sort(taxa_sums(Menu.ps), decreasing = TRUE)[1:20])
top20_Menu.ps.prop = prune_taxa(top20, Menu.ps.prop)

plot_bar(top20_Menu.ps.prop, x="Microhabitat",fill="genus")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("Proportional Read Abundance")

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

#ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_Median_normalized_genus_faceted_by_site.jpg"))


plot_bar(Menu.ps.prop, fill = "genus", facet_grid = "Site_Name")+ 
  geom_bar(aes(color=genus, fill=genus), stat="identity", position="stack")+
  ylab("proportional read abundance")
ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_proportions_genus_faceted_by_site.jpg"))



plot_bar(Menu.ps.prop, fill = "family", facet_grid = "Site_Name")+ 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")+
  ylab("proportional read abundance")



# #Merge samples by donor species
# Menu_donor_merge_MF_GENUS.ps <- merge_samples(Menu_MF_GENUS.ps, "Species..scientific.name.")
# plot_bar(Menu_donor_merge_MF_GENUS.ps, fill = "family") + 
#   geom_bar(aes(color=family, fill=family), stat="identity", position="stack")

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


# Inspect and Save the Species-merged phyloseq object
Menu_taxa_merge.ps
sample_variables(Menu_taxa_merge.ps)
rank_names(Menu_taxa_merge.ps)
tax_table(Menu_taxa_merge.ps)

saveRDS(Menu_taxa_merge.ps, file = paste0("../10_Phyloseq/",Primer,"/",Primer,"_specied_merged_phyloseq.rds"))

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
  scale_color_locuszoom()
ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_alpha_diversity_by_site.jpg"))

plot_richness(Menu_taxa_merge.ps, x="Microhabitat", color="Microhabitat", measures = c("Observed","Shannon","Simpson","InvSimpson","Chao1"))+ 
  geom_boxplot()+
  scale_color_locuszoom()
ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_alpha_diversity_by_microhabitat.jpg"))



#Ordinations----

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
    suffix_rank = "current"
  )

# Menu_taxa_merge.ps %>%
#   tax_fix(
#     min_length = 4,
#     unknowns = c(""),
#     sep = " ", anon_unique = TRUE,
#     suffix_rank = "classified"
#   )

# play with ordinations interactively----
#pseq <- pseq %>%
#  phyloseq_validate()

#ord_explore(pseq) 

#code from the above program... 

#PCoA compositional transformation dist=bray with ellipses
pseq %>%
  tax_transform(rank = "genus", trans = "compositional") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Site_Name", fill = "Site_Name",
    shape = "Microhabitat", 
    alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Site_Name)
  )+
  scale_color_flatui()

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_compositional PcoA with ellipses.jpg"))



#PCA of compositional proportion of prey with arrows

pseq %>%
  tax_transform(rank = "species", trans = "compositional") %>%
  ord_calc(
    method = "PCA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:5,
    colour = "Microhabitat", fill = "Microhabitat",
    shape = "circle", alpha = 0.7,
    size = 2
  )+
  scale_color_locuszoom()

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_PCA_with_arrows_by_microhabitat.jpg"))

#PCA of binary transformed (presence/absence) prey species.
pseq %>%
  tax_transform(rank = "species", trans = "binary") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Site_Name", fill = "Site_Name",
    shape = "circle", alpha = 0.5,
    size = 2
  ) +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Site_Name)
  )+
  scale_color_flatui()

ggsave(paste0("../10_Phyloseq/",Primer,"/",Primer,"_presence_absence_PCA_with_ellipses.jpg"))



print(paste0("Phyloseq exploratory visualizations for this run can be found in ../10_Phyloseq/",Primer,"/"))



