#Eldridge Wisely
#4-26-24
#Compare percent ID columns of local and global assignments, and pick the lowest taxonomic rank shared between any winning global assignments and the local species checklist.

library(tidyverse)
library(prettyunits)
library(taxonomizr)

#Load global obitools database results and local usearch database results ----

#global EMBL database obitools results
obi_result95<-readr::read_delim("../03_obitools_results/BerryCrust_Menu_95_named.tab")

colnames(obi_result95)
obi_result95<-obi_result95%>%
  select(c(ID,TAXID,SCIENTIFIC_NAME,BEST_IDENTITY))

#Load local vsearch results using the in-silico PCR'ed mito file added to lsu and 16S without insilico PCR-----

lca_vsearch<-readr::read_delim("../05_vsearch_obitools_results/lca_crustaceans_Galapagos_top5_comprehensive_galapagos_results.txt", col_names = c("ID","sintax"), delim = '\t')

#parse the sintax taxonomy column into 7 new columns.

lca_vsearch<-lca_vsearch%>%
  separate_wider_delim(sintax, delim = ',', names = c("Domain","Phylum","Class","Order","Family","Genus","Species"),too_few = "align_start")

lca_vsearch$Domain<-str_remove_all(lca_vsearch$Domain,pattern = "d:")
lca_vsearch$Phylum<-str_remove_all(lca_vsearch$Phylum, pattern = "p:")
lca_vsearch$Class<-str_remove_all(lca_vsearch$Class, pattern = "c:")
lca_vsearch$Order<-str_remove_all(lca_vsearch$Order, pattern = "o:")
lca_vsearch$Family<-str_remove_all(lca_vsearch$Family, pattern = "f:")
lca_vsearch$Genus<-str_remove_all(lca_vsearch$Genus, pattern = "g:")
lca_vsearch$Species<-str_remove_all(lca_vsearch$Species, pattern = "s:")
lca_vsearch$Species<-str_replace_all(lca_vsearch$Species, pattern = "_", replacement = " ")

#read in userout for vsearch  

pct_vsearch<-readr::read_delim("../05_vsearch_obitools_results/userout_crustaceans_Galapagos_top5_comprehensive_galapagos_results.txt", col_names =FALSE, delim = '\t')

pct_vsearch<-pct_vsearch%>%
  separate(col = X2, sep = ";",into = c("ACC","sintax"))%>%
  dplyr::rename(ID=X1, vsearch_pctid=X3, vsearch_alnlen=X4, vsearch_mism=X5, vsearch_opens=X6)


#Import the obitools ASVs file so that the sequences of each can be used in making the MOTUS table for MetabaR later.

query_ASVs<-read_delim("../05_vsearch_obitools_results/upper_BerryCrust_Menu_95_named_cleared_tags.tab",col_names = FALSE)

query_ASVs<-query_ASVs%>%
  separate(X1, into = c("ID","Count"), sep=" ")%>%
  dplyr::rename(Sequence=X2)%>%
  select(c("ID","Sequence"))



#combine the LCA taxonomy with the percent_id column (because in vsearch I specified only the top matches so the multiples all share the same pctid with each other, so we can pick any of them)
vsearch_results<- left_join(lca_vsearch, pct_vsearch, by ="ID", multiple ="any")

vsearch_results_seq<-full_join(vsearch_results, query_ASVs, by="ID")

vsearch_results_all_matches<-full_join(lca_vsearch,pct_vsearch, by ="ID")

phylum_only<-vsearch_results_all_matches%>%
  dplyr::filter(is.na(Class))%>%
  dplyr::filter(Domain =="Eukaryota")
nrow(phylum_only)
#41 observations LCA'ed to just Phylum

unique_sp<-vsearch_results_all_matches%>%
  select(Species)%>%
  unique()
nrow(unique_sp)
#46 including NA

unique_genuses<-vsearch_results_all_matches%>%
  select(Genus)%>%
  unique()
nrow(unique_genuses)
#44 including NA


#VSEARCH combine results of obitools and vsearch_global with lca------

#join by ID column

lca_obi_combined<- full_join(vsearch_results_seq, obi_result95, by ="ID")

lca_obi_combined<-lca_obi_combined%>%
  mutate(lca_name = coalesce(Species,Genus,Family,Order,Class,Phylum,Domain))


#Compare global and local database assignments and pick the best assignment for each ASV----

#if BEST_IDENTITY*100 is greater than vsearch_pctid then mutate fish_combined$preferred_name is SCIENTIFIC_NAME.
#if vsearch_pctid is greater than or equal to BEST_IDENTITY*100, then mutate fish_combined$preferred_name is lca_name


best_ID_combined<-lca_obi_combined%>%
  mutate(local_pctid=
           if_else(is.na(vsearch_pctid),
                   0,
                   vsearch_pctid))%>%
  mutate(global_pctid=
           if_else(BEST_IDENTITY=="NA",
                   0,
                   signif(BEST_IDENTITY,3)*100))%>%
  mutate(preferred_pctid =
           if_else(local_pctid >= global_pctid, 
                   local_pctid, 
                   global_pctid))%>%
  mutate(preferred_name=
           if_else(local_pctid >= global_pctid, 
                   lca_name, 
                   SCIENTIFIC_NAME))%>%
  mutate(database=
           if_else(local_pctid >= global_pctid,
                   "local",
                   "global"))


#summarize changes to taxonomic classifications----



#count how many got assigned to local vs. global and number of unassigned.

#total ASVs
total_ASVs=nrow(obi_result95)
total_ASVs
#to make sure we didn't lose any
final_ASVs=nrow(best_ID_combined)
final_ASVs
#1438

#Updated number of total ASVs assigned to a taxon
assigned<-nrow(best_ID_combined[is.na(best_ID_combined$preferred_name) ==FALSE, ])
assigned
#4263

#Percent of all ASVs assigned to a taxon after combining usearch and obitools
(assigned/nrow(best_ID_combined))*100
#37.72233

#compared to just obitools
obi_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
obi_assigned
#1005

#increase in assigned ASVs from global only
assigned - obi_assigned
#3258


lca_vsearch_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$lca_name) ==FALSE, ])
lca_vsearch_assigned
#3747

#After comparing the taxonomic assignments of global and local to get preferred names----

#count the number of times the global database was used for the preferred assignment
globally_assigned<-nrow(best_ID_combined[best_ID_combined$database == 'global'& is.na(best_ID_combined$preferred_name) ==FALSE, ])

globally_assigned
#525

#put these in a new dataframe
global_preferred_assignments<-best_ID_combined%>%
  filter(database=="global")

#print global_preferred_assignments to a csv file to maybe use as supplement.
write.csv(global_preferred_assignments, "BerryCrust_output/global_preferred_assignments_before_final_LCA.csv")

#percentage of all ASVs assigned to the global database
(globally_assigned/nrow(best_ID_combined))*100
#4.645607%
#percentage of identfied ASVs assigned to the global database
(globally_assigned/assigned)*100
#12.31527%

#count the number of times the local database was used for the preferred assignment
locally_assigned<-nrow(best_ID_combined[best_ID_combined$database == 'local'& is.na(best_ID_combined$preferred_name) ==FALSE, ])
locally_assigned
#3738
#percentage of all ASVs
(locally_assigned/nrow(best_ID_combined))*100
#33.07672%

#percentage of identified sequences
(locally_assigned/assigned)*100
#87.68473% 


#number of occurences where local assignment changed the existing global assignment
global_to_local_assigment<-nrow(best_ID_combined[best_ID_combined$database == 'local'& is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
global_to_local_assigment
#480 taxa updated from existing global assignments


#number of unique taxa after choosing the best assignment between global and local:

summary_best_ID_combined<-best_ID_combined%>%
  summarise(n_globally_IDed_taxa= n_distinct(SCIENTIFIC_NAME),
            n_locally_IDed_taxa= n_distinct(lca_name),
            n_combined_IDed_taxa= n_distinct(preferred_name)

)
  
#LCA of global and local ID'ed taxa-----
#need the taxonomizr database to be already prepared
meta_best_ID_combined<-best_ID_combined%>%
  separate(sintax, into=c("sintax-to-genus", "vsearch_species"),sep = "s:", remove = FALSE)

v_species<-gsub("_"," ",meta_best_ID_combined$vsearch_species)
meta_best_ID_combined$vsearch_species<-v_species


taxaId<-getId(v_species,'accessionTaxa.sql')
print(taxaId)

meta_best_ID_combined$v_taxID<-taxaId

#get the taxonomy for the vsearch TaxIDs (even though we have the sintax format already, this is getting it ready to do LCA in this program)

local_levels<-getTaxonomy(taxaId,'accessionTaxa.sql')
print(local_levels)

global_taxaId<-meta_best_ID_combined$TAXID
global_levels<-getTaxonomy(global_taxaId,'accessionTaxa.sql')
print(global_levels)

global_taxa<-meta_best_ID_combined$SCIENTIFIC_NAME


#condenseTaxa(taxa)
##   superkingdom phylum     class      order family genus species
## 1 "Eukaryota"  "Chordata" "Mammalia" NA    NA     NA    NA
#This function can also be fed a large number of grouped hits, e.g. BLAST hits for high throughput sequencing reads after filtering for the best hits for each read, and output a condensed taxonomy for each grouping:


IDs<-meta_best_ID_combined$ID

local_taxas<-as.data.frame(local_levels, row.names = paste(IDs,"_local"))
global_taxas<-as.data.frame(global_levels, row.names = paste(IDs, "_global"))
taxas<-rbind(local_taxas, global_taxas)

doubleIDs<-append(IDs,IDs)
condensed_lca_global_vs_local<-condenseTaxa(taxas, groupings = doubleIDs)

condensed_lca_global_vs_local<-rownames_to_column(as.data.frame(condensed_lca_global_vs_local), var = "ID")

condensed_lca_global_vs_local<-condensed_lca_global_vs_local%>%
  mutate(global_v_local_lca_name = coalesce(species,genus,family,order,class,phylum,superkingdom))


#add global_v_local_lca_name column to best_ID_combined

lca_global_vs_local_assignments<-condensed_lca_global_vs_local%>%
  select(ID,global_v_local_lca_name)

best_ID_combined<-full_join(best_ID_combined, lca_global_vs_local_assignments, by ="ID")

#after visually looking at best_ID_combined for local database sequences that were assigned to very different taxa by global and local databases, and even within the local database assignments (phylum only dataframe), there were over 3,000 sequences that matched to Calanus sinicus with a few HQ619236 HQ619232 HQ619230 HQ619237 HQ619234 HQ619232 HQ619231 HQ619235 HQ619228 from the same study.  I'll check these sequences in the database file to see what's going on.  I also need to remove nans from the database file. #3214 Calanus sinicus sequences, #41 ASVs reduced ID to eukaryota after global_v_local lca step.  Also remove d:Bacteria, and check MT872704, MT672041

#re-ran this code with the new cleaned vsearch local database file results and nothing was phylum only!  the best_ID_combined file looks like they're converging on genus level agreements now too!

#updated database resulted in #440 taxa updated from existing global assignments (down from 480), 44.11474% locally assigned sequences (out of assigned sequences), 3.94% of all ASVs, 446 ASVs had greater or equal percent ID as the global database and were therefore preferred (down from 3738 (although 3214 of those were Calanus sinicus problematic sequences, likely with lots of Ns)), 55% of assigned ASVs are still from the global database, 4.9% of all ASVs.
#6 ASVs were assigned locally that had no assignment with the global obitools database.


#This can happen with a ver incomplete local database.  Next steps:


#Compare the global taxonomic assignments and the local species list----


#When global is still the preferred assignment after comparing pctid, and global_v_local LCA name is NA, find the genus or family in the local database and reduce the ID to that level.

#make a SCIENTIFIC_GENUS column from SCIENTIFIC_NAME, and if SCIENTIFIC_GENUS is found in the local database checklist ("../custom_db/comprehensive_galapagos_crustaceans_list.txt", or "../custom_db/comprehensive_galapagos_fish_list.txt") then, put the genus in local_relation column.  

#read in the comprehensive local checklist and isolate the genus column
local_checklist<-read.csv("../custom_db/comprehensive_galapagos_crustaceans_list.txt", header = FALSE)

local_genuses<-local_checklist%>%
  separate(V1, into = c("local_genus", "local_species"), sep=" ")










#When global is the database instead of local, use the LCA between global and local as preferred name, and the lower pctid (local), as the preferred_pctid----








