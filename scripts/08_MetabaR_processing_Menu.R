#Eldridge Wisely
#Date: 4-29-24
#MetabaR filtering of Menu project for either/both primers


#Metabar and Phyloseq section
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("tibble")       # Needed for converting column to row names
# Load requested package for plotting
library("reshape2")
library("tidyverse")
library("dplyr")        # filter and reformat data frames
library(devtools)
#install_github("tobiasgf/lulu")
require(lulu)
library(metabaR)
library(phyloseq)
library(colorspace)
require(paletteer)
library(ggsci)
library(RColorBrewer)
library(colorRamps)



'%ni%' <- Negate("%in%")

#set variables
Primer<-"BerryCrust"


#Import Metabarlist----
Menu<-readRDS(file= paste("../07_lulu_metabar/",Primer,"_Menu_Metabarlist.rds", sep = ""))

# Compute the number of reads per pcr
Menu$pcrs$nb_reads <- rowSums(Menu$reads)

# Compute the number of motus per pcr
Menu$pcrs$nb_motus <- rowSums(Menu$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Menu <- subset_metabarlist(Menu, table = "pcrs",
                                 indices = Menu$pcrs$nb_reads>0)
summary_metabarlist(Menu)


#Rarefy with Hill numbers and look at rarefaction curves----
Menu.raref = hill_rarefaction(Menu, nboot = 20, nsteps = 10)
head(Menu.raref$hill_table)
gghill_rarefaction(Menu.raref) 


# Define a vector containing the Site_name info for each pcrs 
control_type <- Menu$pcrs$control_type

# Use of gghill_rarefaction requires a vector with named pcrs
control_type <- setNames(control_type,rownames(Menu$pcrs))

# Plot
p <- gghill_rarefaction(Menu.raref, group=control_type)
p + scale_fill_manual(values = c("goldenrod", "brown", "grey")) +
  scale_color_manual(values = c("goldenrod", "brown", "grey")) +
  labs(color="control type")

ggsave(paste0(Primer,"_output/",Primer,"_Hill_rarefaction_and_coverage.jpg"))

# Identify contaminants----

## Identifying extraction contaminants(separately for field,extraction, pcr, and sequencing)----


#extraction negatives list:
MiFish_ext_neg<-c("ExtCon2-1and2","ExtCon2-3and4","GW13neg-1and2","GW13neg-3and4","GW8A-1and2","GW8A-3and4")

BerryCrust_ext_neg<-c("ExtCon2-1and2","ExtCon2-3and4","GW13neg-1and2","GW13neg-3and4","GW7A-1and2","GW7A-3and4")

if(Primer =="MiFish") ext_neg<-MiFish_ext_neg else ext_neg<-BerryCrust_ext_neg
  
Menu <- contaslayer(Menu, 
                          controls = ext_neg,
                          output_col = "not_an_extraction_conta", 
                    method = "all"
)

table(Menu$motus$not_an_extraction_conta)
#summary_metabarlist(Menu)

##Field Negatives----

MiFish_field_neg_pcr<-c("FN2-1and2","FN2-3and4","FN3-1and2","FN3-3and4","FN4-1and2","FN4-3and4","FN5-1and2","FN5-3and4","FN6-1and2","FN6-3and4","GW7A-1and2","GW7A-3and4")

BerryCrust_field_neg_pcr<-c("FN2-3and4","FN3-1and2","FN3-3and4","FN4-1and2","FN4-3and4","FN5-1and2","FN5-3and4","FN6-1and2","FN6-3and4","GW7A-1and2","GW7A-3and4")

if(Primer =="MiFish") field_neg_pcr<-MiFish_field_neg_pcr else field_neg_pcr<-BerryCrust_field_neg_pcr
Menu <- contaslayer(Menu, 
                    controls = field_neg_pcr,
                    output_col = "not_a_field_conta", 
                    method = "max")

table(Menu$motus$not_a_field_conta)
summary_metabarlist(Menu)

# Identifying pcr contaminants
Menu <- contaslayer(Menu, 
                          control_types = "pcr",
                          output_col = "not_a_pcr_conta",
                          method = "max"
)
table(Menu$motus$not_a_pcr_conta)


#Sequencing
Menu <- contaslayer(Menu, 
                          control_types = "sequencing",
                          output_col = "not_a_sequencing_conta", 
                    method = "max"
)

table(Menu$motus$not_a_sequencing_conta)
summary_metabarlist(Menu)



Menu$motus$not_a_conta<-if_else(rowSums(Menu$motus[,c("not_a_field_conta","not_an_extraction_conta","not_a_pcr_conta", "not_a_sequencing_conta")])==4,
                                TRUE,
                                FALSE)
table(Menu$motus$not_a_conta)


# Compute relative abundance of all pcr contaminants together 
a <- data.frame(conta.relab = rowSums(Menu$reads[,!Menu$motus$not_a_conta]) / 
                  rowSums(Menu$reads))
# Add information on control types
a$control_type <- Menu$pcrs$control_type[match(rownames(a), rownames(Menu$pcrs))]

ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) + 
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") + 
  theme_bw() + 
  scale_y_log10()


# flag pcrs with total contaminant relative abundance > 20% of reads)
Menu$pcrs$low_contamination_level <- 
  ifelse(a$conta.relab[match(rownames(Menu$pcrs), rownames(a))]>2e-1,  F, T)

# Proportion of potentially functional (TRUE) vs. failed (FALSE) pcrs
# (controls included) based on this criterion
table(Menu$pcrs$low_contamination_level) / nrow(Menu$pcrs)

##Flag non-target motus ----
#next the tutorial excludes motus without a family ID and plots the similarity scores and removes the ones with low scores.  

#Flag MOTUs corresponding to target (TRUE) vs. non-target (FALSE) taxa 
Menu$motus$target_taxon <- !is.na(Menu$motus$family)


# Proportion of each of these over total number of MOTUs
table(Menu$motus$target_taxon) / nrow(Menu$motus)


# Intersection with extraction contaminant flags (not contaminant = T)
table(Menu$motus$target_taxon, 
      Menu$motus$not_a_conta)


#Next remove contaminants and off-target reads and look at the remaining patterns----
Menu1 <- subset_metabarlist(Menu, "motus", 
                                  indices = rowSums(Menu$motus[,c("not_a_conta", "target_taxon")]) == 2)


# Compute the number of reads per pcr
Menu1$pcrs$nb_reads <- rowSums(Menu1$reads)

# Compute the number of motus per pcr
Menu1$pcrs$nb_motus <- rowSums(Menu1$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Menu1 <- subset_metabarlist(Menu1, table = "pcrs",
                                  indices = Menu1$pcrs$nb_reads>0)
summary_metabarlist(Menu1)

#161 motus remaining average 19.79 -MiFish Menu

# Plot the unweighted distribution of MOTUs similarity scores 
a <- 
  ggplot(Menu1$motus, aes(x=pct_id)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 97, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# MOTUs")

# Same for the weighted distribution
b <- 
  ggplot(Menu1$motus, 
         aes(x=pct_id, y = after_stat(count), weight = COUNT)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 97, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# Reads")

# Combine plots into one
library(cowplot)
ggdraw() + 
  draw_plot(a, x=0, y=0, width = 0.5) + 
  draw_plot(b, x=0.5, y=0, width = 0.5)
#make a note of the cutoff for best match %.  


########This is where I can play with different values for ID#####

# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
Menu1$motus$not_degraded <-
  ifelse(Menu1$motus$pct_id < 96, F, T)

# Proportion of each of these over total number of MOTUs
table(Menu1$motus$not_degraded) / nrow(Menu1$motus)

# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
Menu$motus$not_degraded <-
  ifelse(Menu$motus$pct_id < 96 &!is.na(Menu$motus$pct_id), F, T)

# Proportion of each of these over total number of MOTUs
table(Menu$motus$not_degraded) / nrow(Menu$motus)




# Proportion of each of these over total number of MOTUs
table(Menu1$motus$not_degraded) / nrow(Menu1$motus)
#> For Menu1
#
#FALSE      TRUE 
#0.02484472 0.97515528

# Intersection with other flags
table(Menu$motus$target_taxon, 
      Menu$motus$not_a_conta, 
      Menu$motus$not_degraded)



#Moving on to detecting PCR outliers by sequencing depth
ggplot(Menu$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") + 
  geom_vline(xintercept = 70, lty=2, color="orange") + # threshold
  scale_x_log10() + 
  labs(x="# Reads (with all MOTUs and PCRs)", 
       y="# PCRs") +
  theme_bw() + 
  theme(panel.grid = element_blank())


#play with the xintercept above to visualize your cutoff level

# Flag pcrs with an acceptable sequencing depth (TRUE) or inacceptable one (FALSE)
Menu$pcrs$seqdepth_ok <- ifelse(Menu$pcrs$nb_reads < 70, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(Menu$pcrs$seqdepth_ok[Menu$pcrs$type=="sample"]) /
  nrow(Menu$pcrs[Menu$pcrs$type=="sample",])
#and for the semi-cleaned dataset
Menu1$pcrs$seqdepth_ok <- ifelse(Menu1$pcrs$nb_reads < 70, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(Menu1$pcrs$seqdepth_ok[Menu1$pcrs$type=="sample"]) /
  nrow(Menu1$pcrs[Menu1$pcrs$type=="sample",])

#     FALSE       TRUE 
#      0.032 0.968

#Compare replication of PCR replicates ----
#They use this next function to compare the similarity of biological controls with the expectation they will be more similar to each other than to other samples.  

# # Subsetting the metabarlist
# Menu_sub <- subset_metabarlist(Menu, 
#                                      table = "pcrs", 
#                                      indices = Menu$pcrs$nb_reads>0 & (
#                                        is.na(Menu$pcrs$control_type) ))
# 
#                                           #|Menu$pcrs$low_contamination_level=="TRUE"))
# 
# # First visualization
# comp1 = pcr_within_between(Menu_sub)
# check_pcr_thresh(comp1)
# 
# #cool, this looks good.  Now... for flagging not-well-replicated pcrs. 
# # Subsetting the metabarlist
# Menu1_sub <- subset_metabarlist(Menu1, 
#                                       table = "pcrs", 
#                                       indices = Menu1$pcrs$nb_reads>0 & (
#                                         is.na(Menu1$pcrs$control_type) ))
# 
# 
# 
# # First visualization
# comp2 = pcr_within_between(Menu1_sub)
# check_pcr_thresh(comp2)
# 
# #The distance between samples was much higher than within samples!

#However, when flagging for replication, it's best to use the unfiltered dataset. You lose fewer pcrs because there's more to compare to when including non-target species like fish in the case of BerryCrust primers.

# # Flagging
# Menu_sub <- pcrslayer(Menu_sub, output_col = "replicating_pcr", plot = F, thresh.method = "intersect")#,method = "pairwise")
# 
# 
# # Proportion of replicating pcrs (TRUE)
# table(Menu_sub$pcrs$replicating_pcr) /
#   nrow(Menu_sub$pcrs)
# #98.461538 % of Menu before removal of extraction contaminants etc.
# 
# # Intersection with the sequencing depth criterion
# table(Menu_sub$pcrs$seqdepth_ok,
#       Menu_sub$pcrs$replicating_pcr)
# 
# # Distinguish between pcrs obtained from samples from controls
# mds = check_pcr_repl(Menu_sub,
#                      groups = Menu_sub$pcrs$type,
#                      funcpcr = Menu_sub$pcrs$replicating_pcr)
# mds + labs(color = "pcr type") + scale_color_manual(values = c("gray", "cyan4"))
# 
# #add results of this to the original Menu (and to Menu1)
# Menu$pcrs$replicating_pcr <- NA
# Menu$pcrs[rownames(Menu_sub$pcrs),"replicating_pcr"] <- Menu_sub$pcrs$replicating_pcr
# 
# 
# Menu1$pcrs$replicating_pcr <- NA
# Menu1$pcrs[rownames(Menu_sub$pcrs),"replicating_pcr"] <- Menu_sub$pcrs$replicating_pcr
# 
# 
# mds <- check_pcr_repl(Menu_sub, groups = Menu_sub$pcrs$rep)
# mds + labs(color = "PCR replicate")
# 
# mds <- check_pcr_repl(Menu1_sub, groups = Menu1_sub$pcrs$rep)
# mds + labs(color = "PCR replicate")

#Summarize noise in the dataset----
##Motus artefacts pie chart----

# Create a table of MOTUs quality criteria
# noise is identified as FALSE in Menu, the "!" transforms it to TRUE
motus.qual <- !Menu$motus[,c("not_a_conta", "target_taxon", "not_degraded")]
colnames(motus.qual) <- c("contaminant", "untargeted_taxon", "degraded_seq")


# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(Menu$motus$COUNT~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)
apply(motus.qual, 2, function(x) sum(Menu$motus$COUNT[x])/sum(Menu$motus$COUNT))

tmp.motus <-
  apply(sapply(1:ncol(motus.qual), function(x) {
    ifelse(motus.qual[,x]==T, colnames(motus.qual)[x], NA)}), 1, function(x) {
      paste(sort(unique(x)), collapse = "|")
    })
tmp.motus <- as.data.frame(gsub("^$", "not_artefactual", tmp.motus))
colnames(tmp.motus) <-  "artefact_type"

ggplot(tmp.motus, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") +
  coord_polar(theta="y") + theme_void() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.direction = "vertical") +
  ggtitle(paste(Primer,"Menu ASVs artefacts overview"))

ggsave(filename = paste0(Primer,"_output/",Primer,"_Menu_ASVs_artifacts_overview.jpg"))

##PCR artefacts pie chart ----

# Create a table of pcrs quality criteria
# noise is identified as FALSE in Menu, the "!" transforms it to TRUE
pcrs.qual <- !Menu$pcrs[,c("low_contamination_level", "seqdepth_ok")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used
# excluding controls
prop.table(table(apply(pcrs.qual[Menu$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[Menu$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[Menu$pcrs$type=="sample",])

tmp.pcrs <-
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[Menu$pcrs$type=="sample",x]==T,
           colnames(pcrs.qual)[x], NA)}), 1, function(x) {
             paste(sort(unique(x)), collapse = "|")
           })
tmp.pcrs <- as.data.frame(gsub("^$", "not_artefactual", tmp.pcrs))

colnames(tmp.pcrs) <- "artefact_type"

ggplot(tmp.pcrs, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") +
  coord_polar(theta="y") + theme_void() +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.direction = "vertical") +
  ggtitle(paste(Primer,"Menu PCRs artefacts overview"))

ggsave(filename = paste0(Primer,"_output/",Primer,"_Menu_PCRs_artifacts_overview.jpg"))


#Tag jumps filtering ----

##In unfiltered dataset ----

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 3e-2, 5e-2)

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(Menu,x))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- Menu$pcrs$control_type[match(tmp$sample, rownames(Menu$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

ggplot(tmp2, aes(x=as.factor(threshold), y=value)) +
  geom_boxplot(color="grey40") +
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.01"), col="orange", lty=2) +
  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=5) +
  theme_bw() +
  scale_y_log10() +
  labs(x="ASV pcr : total abundance filtering threshold", y="# Reads/ASVs") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=40, h=1),
        legend.position = "none")

ggsave(paste0(Primer,"_output/",Primer,"_threshold_filtering_for_tag_jumps.jpg"))

#Remove spurious signal ----
# Use tag-jump corrected metabarlist with the threshold identified above----
MiFish_threshold<-"t_0.01"
BerryCrust_threshold<-"t_1e-04"

if (Primer =="MiFish") t_val<-MiFish_threshold else t_val<-BerryCrust_threshold

tmp <- tests[[t_val]]

## Subset on MOTUs: ----
#we keep motus that are defined as TRUE following the
# three criteria below (sum of three TRUE is equal to 3 with the rowSums function)
tmp <- subset_metabarlist(tmp, "motus",
                          indices = rowSums(tmp$motus[,c("not_a_conta", "target_taxon",
                                                         "not_degraded")]) == 3)
summary_metabarlist(tmp)

##subset on PCRs table criteria----

# Subset on pcrs and keep only controls
#Menu_clean!----
Menu_clean <- subset_metabarlist(tmp, "pcrs",
                                     indices = rowSums(tmp$pcrs[,c("low_contamination_level",
                                                                   "seqdepth_ok")]) == 2 &
                                       tmp$pcrs$type == "sample")
summary_metabarlist(Menu_clean)


#check for empty motus
if(sum(colSums(Menu_clean$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(Menu_clean$reads)==0)>0){print("empty pcrs present")}

#recount reads now that they're cleaned
Menu_clean$motus$count = colSums(Menu_clean$reads)
Menu_clean$pcrs$nb_reads_postmetabaR = rowSums(Menu_clean$reads)
Menu_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(Menu_clean$reads>0, T, F))


#subset the metabarlist to all PCRs with more than 0 reads
Menu_clean1 <- subset_metabarlist(Menu_clean, table = "pcrs",
                           indices = rowSums(Menu_clean$reads)>0)

if(sum(colSums(Menu_clean1$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(Menu_clean1$reads)==0)>0){print("empty pcrs present")}



summary_metabarlist(Menu_clean1)
summary_metabarlist(Menu_clean)




#Compare before and after basic stats ----
check <- melt(Menu_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR",
                                     "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "richness", "abundance")

ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot( color = "darkgrey") +
  geom_jitter(alpha=0.1, color = "darkgrey") +
  theme_bw() +
  facet_wrap(~type, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle=45, h=1))
ggsave(paste0(Primer,"_output/",Primer,"ASVs_Before_and_After_MetabaR.jpg"))


#Compare before and after beta diversity----

# Get row data only for samples
tmp <- subset_metabarlist(Menu, table = "pcrs",
                          indices = Menu$pcrs$type == "sample")

# Add sample biological information for checks

##put Site_Name into pcrs for tmp and clean (and clean1 in case there were empty PCRs)----

tmp$pcrs$Site_Name <- tmp$samples$Site_Name[match(tmp$pcrs$sample_id, rownames(tmp$samples))]


Menu_clean$pcrs$Site_Name <-
  Menu_clean$samples$Site_Name[match(Menu_clean$pcrs$sample_id,
                                        rownames(Menu_clean$samples))]

Menu_clean1$pcrs$Site_Name <-
  Menu_clean1$samples$Site_Name[match(Menu_clean1$pcrs$sample_id,
                                     rownames(Menu_clean1$samples))]


##put Site_Name into pcrs for tmp and clean (and clean1 in case there were empty PCRs)----

tmp$pcrs$Microhabitat <- tmp$samples$Microhabitat[match(tmp$pcrs$sample_id, rownames(tmp$samples))]


Menu_clean$pcrs$Microhabitat <-
  Menu_clean$samples$Microhabitat[match(Menu_clean$pcrs$sample_id,
                                     rownames(Menu_clean$samples))]

Menu_clean1$pcrs$Microhabitat <-
  Menu_clean1$samples$Microhabitat[match(Menu_clean1$pcrs$sample_id,
                                      rownames(Menu_clean1$samples))]


# Build PCoA ordinations

##PCOA ordination of raw vs clean data by Site----
# Build PCoA ordinations
mds1 <- check_pcr_repl(tmp,
                       groups = tmp$pcrs$Site_Name)
mds2 <- check_pcr_repl(Menu_clean1,
                       groups = 
                         Menu_clean1$pcrs$Site_Name)
                         

a <- mds1 + labs(color = "Site_name") +
  scale_color_ucscgb()+
  theme(legend.position = "none") +
  ggtitle("Raw data")
b <- mds2 + labs(color = "Site_name") +
  scale_color_ucscgb()+
  ggtitle("Clean data")

# Assemble plots
leg <- get_legend(b + guides(shape=F) +
                    theme(legend.position = "right",
                          legend.direction = "vertical"))
ggdraw() +
  draw_plot(a, x=0, y=0, width = 0.4, height = 1) +
  draw_plot(b + guides(color=F, shape=F), x=0.42, y=0, width = 0.4, height = 1) +
  draw_grob(leg, x=0.4, y=0)

ggsave(paste0(Primer,"_output/",Primer,"_raw_vs_clean_beta_diversity_of_ASVs_by_Site.jpg"))

##PCOA ordination of raw vs clean data by Microhabitat----
# Build PCoA ordinations
mds1 <- check_pcr_repl(tmp,
                       groups = tmp$pcrs$Microhabitat)
mds2 <- check_pcr_repl(Menu_clean1,
                       groups = 
                         Menu_clean1$pcrs$Microhabitat)


a <- mds1 + labs(color = "Microhabitat") +
  scale_color_ucscgb()+
  theme(legend.position = "none") +
  ggtitle("Raw data")
b <- mds2 + labs(color = "Microhabitat") +
  scale_color_ucscgb()+
  ggtitle("Clean data")

# Assemble plots
leg <- get_legend(b + guides(shape=F) +
                    theme(legend.position = "right",
                          legend.direction = "vertical"))
ggdraw() +
  draw_plot(a, x=0, y=0, width = 0.4, height = 1) +
  draw_plot(b + guides(color=F, shape=F), x=0.42, y=0, width = 0.4, height = 1) +
  draw_grob(leg, x=0.4, y=0)

ggsave(paste0(Primer,"_output/",Primer,"_raw_vs_clean_beta_diversity_of_ASVs_by_Microhabitat.jpg"))


#Aggregate PCRs by sample----

Menu_clean_agg <- aggregate_pcrs(Menu_clean1)
summary_metabarlist(Menu_clean_agg)

#Look at rarefaction curves of cleaned data by sample----

Menu_clean_agg.raref = hill_rarefaction(Menu_clean_agg, nboot = 20, nsteps = 10)
#Location <- paste(Menu_clean_agg$samples$Site_Name, Menu_clean_agg$samples$Microhabitat)
#Location <- setNames(Location,rownames(Menu_clean_agg$samples))


# Define a vector containing the Site_name info for each pcrs 
Site <- Menu_clean_agg$samples$Site_Name

# Use of gghill_rarefaction requires a vector with named pcrs
Site <- setNames(Site,rownames(Menu_clean_agg$samples))


#plot
p <- gghill_rarefaction(Menu_clean_agg.raref, group=Site)
p + scale_fill_ucscgb() +
  scale_color_ucscgb() +
  labs(color="Site")

ggsave(paste0(Primer,"_output/",Primer,"_Hill_rarefaction_cleaned_samples_by_Site.jpg"))

#Make phylogenetic plots!!!----


#Menu_clean_agg$motus<-Menu_clean_agg$motus%>%
#  mutate(Species=str_split_i(pattern = " ",species,i = -1))

menu_clean_agg_nonpassive<-subset_metabarlist(Menu_clean_agg,
                                              table = "samples",
                                              indices=!grepl("passive", Menu_clean_agg$samples$Microhabitat))

#check for empty motus
if(sum(colSums(menu_clean_agg_nonpassive$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(menu_clean_agg_nonpassive$reads)==0)>0){print("empty pcrs present")}

#recount reads now that they're cleaned
menu_clean_agg_nonpassive$motus$count = colSums(menu_clean_agg_nonpassive$reads)
menu_clean_agg_nonpassive$pcrs$nb_reads_postmetabaR = rowSums(menu_clean_agg_nonpassive$reads)
menu_clean_agg_nonpassive$pcrs$nb_motus_postmetabaR = rowSums(ifelse(menu_clean_agg_nonpassive$reads>0, T, F))


#subset the metabarlist to all PCRs with more than 0 reads
menu_clean_agg_nonpassive1 <- subset_metabarlist(menu_clean_agg_nonpassive, table = "pcrs",
                                  indices = rowSums(menu_clean_agg_nonpassive$reads)>0)

if(sum(colSums(menu_clean_agg_nonpassive1$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(menu_clean_agg_nonpassive1$reads)==0)>0){print("empty pcrs present")}
menu_clean_agg_nonpassive1$motus$count = colSums(menu_clean_agg_nonpassive1$reads)


summary_metabarlist(Menu_clean_agg)


summary_metabarlist(menu_clean_agg_nonpassive1)





taxo.col <- c(
  "superkingdom", "phylum", "class",
  "order", "family", "genus","species")

eukaryota <- subset_metabarlist(menu_clean_agg_nonpassive1,
                                   table = "motus",
                                   indices = grepl("Eukaryota", 
                                                   menu_clean_agg_nonpassive1$motus$superkingdom))

taxo.col1 <- c("phylum", "class",
  "order", "family", "genus","species")


####BerryCrust section----
if (Primer=="BerryCrust") {

  malacostraca <- subset_metabarlist(menu_clean_agg_nonpassive1,
                                 table = "motus",
                                 indices = grepl("Malacostraca", 
                                                 menu_clean_agg_nonpassive1$motus$class))

  taxo.col2 <- c("class",
                 "order", "family", "genus","species")
  ggtaxplot(malacostraca, taxo.col2)

  ggsave("BerryCrust_output/BerryCrust_phylogenetic_diversity_Menu_Malacostraca.jpg")
  
  decapoda <- subset_metabarlist(menu_clean_agg_nonpassive1,
                                     table = "motus",
                                     indices = grepl("Decapoda", 
                                                     menu_clean_agg_nonpassive1$motus$order))
  
  ggtaxplot(decapoda, taxo.col)
  
  ggsave("BerryCrust_output/BerryCrust_phylogenetic_diversity_Menu_Decapoda.jpg")
  
  arthropoda <- subset_metabarlist(menu_clean_agg_nonpassive1,
                                 table = "motus",
                                 indices = grepl("Arthropoda", 
                                                 menu_clean_agg_nonpassive1$motus$phylum))
  
  ggtaxplot(arthropoda, taxo.col1)
  
  pgrapsussocius<-subset_metabarlist(menu_clean_agg_nonpassive1,
                                     table = "motus",
                                     indices = grepl("Pachygrapsus socius", 
                                                     menu_clean_agg_nonpassive1$motus$species))
  pgrapsussocius_readcounts<-sum(pgrapsussocius$reads)
  pgrapsussocius_readcounts
  BerryCrust_readcounts<-sum(menu_clean_agg_nonpassive1$reads)
  BerryCrust_readcounts
  
  #percent reads assigned to Pachygrapsus socius
  pgrapsussocius_readcounts/BerryCrust_readcounts

  arthropoda_hammerhead_nurseries <- subset_metabarlist(arthropoda,
                                                        table = "samples",
                                                        indices = (arthropoda$samples$Site_Name == "Puerto Grande" |arthropoda$samples$Site_Name == "Cartago Bay, Isabela"))
  
  ggtaxplot(arthropoda_hammerhead_nurseries, taxo.col)
  
  ggsave("BerryCrust_output/Arthropoda_phylogeny_hammerhead_nurseries.jpg")
  arthropoda_non_hammerhead_nurseries <- subset_metabarlist(arthropoda,
                                                            table = "samples",
                                                            indices = (arthropoda$samples$Site_Name != "Puerto Grande" &
                                                                         arthropoda$samples$Site_Name != "Cartago Bay, Isabela"))
  
  ggtaxplot(arthropoda_non_hammerhead_nurseries, taxo.col)
  
  ggsave("BerryCrust_output/Arthropoda_phylogeny_non_hammerhead_nurseries.jpg")
  

}else {print("not BerryCrust")}
  
###MiFish section ----
if (Primer=="MiFish"){

  MiFish_read_counts<-sum(menu_clean_agg_nonpassive1$reads)
  MiFish_read_counts
  Mugil_metabarlist<-subset_metabarlist(menu_clean_agg_nonpassive1,table = "motus",indices = grepl("Mugil",menu_clean_agg_nonpassive1$motus$genus))
  Mugil_read_counts<-sum(Mugil_metabarlist$reads)
  Mugil_read_counts
  
  Mugil_read_counts/MiFish_read_counts
  
  chordata<-subset_metabarlist(menu_clean_agg_nonpassive1,
                               table = "motus",
                               indices = grepl("Chordata",
                                               menu_clean_agg_nonpassive1$motus$phylum))
  
  ggtaxplot(chordata, taxo.col1)
  
  ggsave(paste0(Primer,"_output/",Primer,"_phylogenetic_diversity_Menu_Eukaryota.jpg"))
  
  
  summary_metabarlist(chordata)
  
  chondrichthyes <- subset_metabarlist(menu_clean_agg_nonpassive1,
                                 table = "motus",
                                 indices = grepl("Chondrichthyes", 
                                                 menu_clean_agg_nonpassive1$motus$class))
  
  ggtaxplot(chondrichthyes, taxo.col)
  
  ggsave("MiFish_output/MiFish_phylogenetic_diversity_Menu_Chondrichthyes.jpg")
  
  hammerhead_nurseries <- subset_metabarlist(menu_clean_agg_nonpassive1,
                                                        table = "samples",
                                                        indices = (menu_clean_agg_nonpassive1$samples$Site_Name == "Puerto Grande" |menu_clean_agg_nonpassive1$samples$Site_Name == "Cartago Bay, Isabela"))
 
   # Compute the number of reads per pcr
  hammerhead_nurseries$pcrs$nb_reads <- rowSums(hammerhead_nurseries$reads)
  
  # Compute the number of motus per pcr
  hammerhead_nurseries$pcrs$nb_motus <- rowSums(hammerhead_nurseries$reads>0)
  
  #subset the metabarlist to all PCRs with more than 0 reads
  hammerhead_nurseries <- subset_metabarlist(hammerhead_nurseries, table = "pcrs",
                             indices = hammerhead_nurseries$pcrs$nb_reads>0)
  summary_metabarlist(hammerhead_nurseries)
  
  ggtaxplot(hammerhead_nurseries, taxo.col)
  
  ggsave("MiFish_output/Menu_hammerhead_nurseries.jpg")
  
  
  non_hammerhead_nurseries <- subset_metabarlist(menu_clean_agg_nonpassive1,
                                                            table = "samples",
                                                            indices = (menu_clean_agg_nonpassive1$samples$Site_Name != "Puerto Grande" &
                                                                         menu_clean_agg_nonpassive1$samples$Site_Name != "Cartago Bay, Isabela"))
  
  ggtaxplot(non_hammerhead_nurseries, taxo.col)
  
  ggsave("MiFish_output/Menu_non_hammerhead_nurseries.jpg")
  

  
  
  
} else print("not MiFish")
  
  
  
  

#save intermediate files for each primer----

write.csv(Menu_clean$reads, paste("../07_lulu_metabar/",Primer,"_Metabar_cleaned_reads.csv",sep=""))
write.csv(Menu_clean$motus, paste("../07_lulu_metabar/",Primer,"_Metabar_cleaned_motus.csv",sep=""))
write.csv(Menu_clean$pcrs, paste("../07_lulu_metabar/",Primer,"_Metabar_cleaned_pcrs.csv",sep=""))
write.csv(Menu_clean$samples, paste("../07_lulu_metabar/",Primer,"_Metabar_cleaned_samples.csv",sep=""))

summary_metabarlist(Menu_clean_agg)
summary_metabarlist(Menu_clean1)
summary_metabarlist(Menu_clean)
summary_metabarlist((menu_clean_agg_nonpassive1))

saveRDS(Menu_clean_agg, file=paste("../07_lulu_metabar/",Primer,"_Menu_Clean_Metabarlist.rds",sep=""))

print(paste0("Done processing ",Primer," with MetabaR, visualizations can be found in 07_lulu_metabar/ and phylogenetic diversity plots are in ",Primer,"_output/"))

