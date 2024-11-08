#Eldridge Wisely
#4-26-24
#Compare percent ID columns of local and global assignments, and pick the lowest taxonomic rank shared between any winning global assignments and the local species checklist.

library(tidyverse)
library(prettyunits)
library(taxonomizr)

#set variables
Primer ="MiFish" #Options MiFish or BerryCrust

 
Local_advantage=TRUE #default is TRUE

#Local advantage controls the behavior of the program when global_pctid = local_pctid.  If set to TRUE, then the local name becomes the preferred name when the percent ID values are equal.  If not, the global name remains the preferred name, and the script will check for the presence of the global species, genus, and family in the local checklist made by combining information from GBIF, OBIS, and the Darwin center, and the assignment will be downgraded to the lowest shared taxonomic rank between the global assignment and the local checklist.


##################### Begin script ################

#This takes a long time and plenty of space, but is necessary for this script
prepareDatabase('accessionTaxa.sql')
#accessed 4-26-24

# Local copies of my data 

# #the final crustacean files I used
# BerryCrust_vsearch_lca_file <- "../05_vsearch_obitools_results/lca_crustaceans_Galapagos_top5_comprehensive_galapagos_results.txt"
# 
# BerryCrust_vsearch_userout_file <- "../05_vsearch_obitools_results/userout_crustaceans_Galapagos_top5_comprehensive_galapagos_results.txt"
# 
# BerryCrust_obitools_results_file<-"../03_obitools_results/BerryCrust_Menu_95_named.tab"
# 
# BerryCrust_ASVs_file<-"../05_vsearch_obitools_results/upper_BerryCrust_Menu_95_named_cleared_tags.tab"
# 
# #the final fish files I used
# #MiFish_vsearch_lca_file<-"../05_vsearch_obitools_results/no_prob_v2_lca_top5_comprehensive_galapagos_results.txt"
# 
# #MiFish_vsearch_userout_file<-"../05_vsearch_obitools_results/no_prob_v2_userout_top5_comprehensive_galapagos_results.txt"
# 
# MiFish_obitools_results_file<-"../03_obitools_results/MiFish_Menu_95_named.tab"
# 
# MiFish_ASVs_file<-"../05_vsearch_obitools_results/upper_MiFish_Menu_95_named_cleared_tags.tab"
# 
# #new fish files 5-1 with the database cleaned the same as the crustaceans file
# MiFish_vsearch_lca_file<-"../05_vsearch_obitools_results/lca_MiFish5-1db_Galapagos_top5_comprehensive_galapagos_results.txt"
# 
# MiFish_vsearch_userout_file<-"../05_vsearch_obitools_results/userout_MiFish5-1db_Galapagos_top5_comprehensive_galapagos_results.txt"
# 

#GitHub copy of my data (for release upon publication)

#MiFish 

MiFish_vsearch_userout_file<-"data/script_06_inputs/userout_MiFish5-1db_Galapagos_top5_comprehensive_galapagos_results.txt"

MiFish_vsearch_lca_file<-"data/script_06_inputs/lca_MiFish5-1db_Galapagos_top5_comprehensive_galapagos_results.txt"

MiFish_obitools_results_file<-"data/script_06_inputs/MiFish_Menu_95_named.tab"

MiFish_ASVs_file<-"data/script_06_inputs/upper_MiFish_Menu_95_named_cleared_tags.tab"

#BerryCrust
BerryCrust_vsearch_lca_file <- "data/script_06_inputs/lca_crustaceans_Galapagos_top5_comprehensive_galapagos_results.txt"

BerryCrust_vsearch_userout_file <- "data/script_06_inputs/userout_crustaceans_Galapagos_top5_comprehensive_galapagos_results.txt"

BerryCrust_obitools_results_file<-"data/script_06_inputs/BerryCrust_Menu_95_named.tab"

BerryCrust_ASVs_file<-"data/script_06_inputs/upper_BerryCrust_Menu_95_named_cleared_tags.tab"

#Load global obitools database results and local usearch database results ----

if (Primer=="MiFish"){
#global EMBL database obitools results
obi_result95<-readr::read_delim(MiFish_obitools_results_file)
}else{
  obi_result95<-readr::read_delim(BerryCrust_obitools_results_file)}

colnames(obi_result95)
obi_result95<-obi_result95%>%
  select(c(ID,TAXID,SCIENTIFIC_NAME,BEST_IDENTITY,COUNT))

#Load local vsearch results using the in-silico PCR'ed mito file added to lsu and 16S without insilico PCR-----
if (Primer=="MiFish"){
lca_vsearch<-readr::read_delim(MiFish_vsearch_lca_file, col_names = c("ID","sintax"), delim = '\t')
}else{
  lca_vsearch<-readr::read_delim(BerryCrust_vsearch_lca_file, col_names = c("ID","sintax"), delim = '\t')}


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

if (Primer=="MiFish"){
pct_vsearch<-readr::read_delim(MiFish_vsearch_userout_file, col_names =FALSE, delim = '\t')
}else{
  pct_vsearch<-readr::read_delim(BerryCrust_vsearch_userout_file, col_names =FALSE, delim = '\t')
}
pct_vsearch<-pct_vsearch%>%
  separate(col = X2, sep = ";",into = c("ACC","sintax"))%>%
  dplyr::rename(ID=X1, vsearch_pctid=X3, vsearch_alnlen=X4, vsearch_mism=X5, vsearch_opens=X6)


#Import the obitools ASVs file so that the sequences of each can be used in making the MOTUS table for MetabaR later.

if (Primer=="MiFish"){
query_ASVs<-read_delim(MiFish_ASVs_file,col_names = FALSE)
}else{
  query_ASVs<-read_delim(BerryCrust_ASVs_file,col_names = FALSE)
}

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

if (Local_advantage ==TRUE){
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
}else{
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
             if_else(local_pctid > global_pctid, 
                     local_pctid, 
                     global_pctid))%>%
    mutate(preferred_name=
             if_else(local_pctid > global_pctid, 
                     lca_name, 
                     SCIENTIFIC_NAME))%>%
    mutate(database=
             if_else(local_pctid > global_pctid,
                     "local",
                     "global"))
  
  
}

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
#546

#Percent of all ASVs assigned to a taxon after combining usearch and obitools
(assigned/nrow(best_ID_combined))*100
#37.72233

#compared to just obitools
obi_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
obi_assigned
#542

#increase in assigned ASVs from global only
assigned - obi_assigned
#4


lca_vsearch_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$lca_name) ==FALSE, ])
lca_vsearch_assigned
#188

#After comparing the taxonomic assignments of global and local to get preferred names----

#count the number of times the global database was used for the preferred assignment
globally_assigned<-nrow(best_ID_combined[best_ID_combined$database == 'global'& is.na(best_ID_combined$preferred_name) ==FALSE, ])

globally_assigned
#403

#put these in a new dataframe
global_preferred_assignments<-best_ID_combined%>%
  filter(database=="global")

#print global_preferred_assignments to a csv file to maybe use as supplement.
write.csv(global_preferred_assignments, paste0(Primer,"_output/global_preferred_assignments_before_final_LCA.csv"))

#percentage of all ASVs assigned to the global database
(globally_assigned/nrow(best_ID_combined))*100
#28.02503%
#percentage of identified ASVs assigned to the global database
(globally_assigned/assigned)*100
#73.80952%

#count the number of times the local database was used for the preferred assignment
locally_assigned<-nrow(best_ID_combined[best_ID_combined$database == 'local'& is.na(best_ID_combined$preferred_name) ==FALSE, ])
locally_assigned
#23
#percentage of all ASVs
(locally_assigned/nrow(best_ID_combined))*100
#1.599444%

#percentage of identified sequences
(locally_assigned/assigned)*100
#4.212454% 


#number of occurences where local assignment changed the existing global assignment
global_to_local_assigment<-nrow(best_ID_combined[best_ID_combined$database == 'local'& is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
global_to_local_assigment
#19 taxa updated from existing global assignments


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
#print(taxaId)

meta_best_ID_combined$v_taxID<-taxaId

#get the taxonomy for the vsearch TaxIDs (even though we have the sintax format already, this is getting it ready to do LCA in this program)

local_levels<-getTaxonomy(taxaId,'accessionTaxa.sql')
#print(local_levels)

global_taxaId<-meta_best_ID_combined$TAXID
global_levels<-getTaxonomy(global_taxaId,'accessionTaxa.sql')
#print(global_levels)

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


##When global is the database instead of local, use the LCA between global and local as preferred name----
if (Local_advantage ==TRUE){
best_ID_combined<-best_ID_combined%>%
  mutate(preferred_name=
           if_else(global_pctid > local_pctid & is.na(lca_name)==FALSE,
                   global_v_local_lca_name, 
                   preferred_name))%>%
  mutate(database=
           if_else(global_pctid > local_pctid & is.na(lca_name)==FALSE,
                   "lca_global_v_local",
                   database))
}else{
  best_ID_combined<-best_ID_combined%>%
    mutate(preferred_name=
             if_else(global_pctid >= local_pctid & is.na(lca_name)==FALSE,
                     global_v_local_lca_name, 
                     preferred_name))%>%
    mutate(database=
             if_else(global_pctid >= local_pctid & is.na(lca_name)==FALSE,
                     "lca_global_v_local",
                     database))
}

#put the global names and lca reassignments in a new dataframe
global_preferred_assignments<-best_ID_combined%>%
  filter(database%in%c("lca_global_v_local","global"))

#print global_preferred_assignments to a csv file to maybe use as supplement.
write.csv(global_preferred_assignments, paste0(Primer,"_output/global_preferred_assignments_after_obi-vsearch_LCA.csv"))


#updated database resulted in #440 taxa updated from existing global assignments (down from 480), 44.11474% locally assigned sequences (out of assigned sequences), 3.94% of all ASVs, 446 ASVs had greater or equal percent ID as the global database and were therefore preferred (down from 3738 (although 3214 of those were Calanus sinicus problematic sequences, likely with lots of Ns)), 55% of assigned ASVs are still from the global database, 4.9% of all ASVs.
#6 ASVs were assigned locally that had no assignment with the global obitools database.


#This can happen with a fairly incomplete local database.  Next steps:


#Compare the global taxonomic assignments and the local species checklist----


##When global is still the database after comparing pctid, and global_v_local LCA name is NA, find the genus or family in the local checklist and reduce the ID to that level.----


#read in the comprehensive local checklist and isolate the genus column


if (Primer=="MiFish"){
local_checklist<-read.csv("../custom_db/comprehensive_galapagos_fish_list.txt", header = FALSE)

#clean it up
local_checklist<-local_checklist%>%
  dplyr::mutate(V1=gsub("Gen. ", "",V1))%>%
  dplyr::mutate(V1=gsub("indet. ", "",V1))%>%
  dplyr::mutate(V1=gsub("\"", "",V1))%>%
  dplyr::mutate(V1=gsub("sp. ", "",V1))%>%
  dplyr::mutate(V1=gsub("cf. ", "",V1))%>%
  dplyr::rename(Scientific_name=V1)

}else{
  local_checklist<-read.csv("../custom_db/comprehensive_galapagos_crustaceans_list.txt", header = FALSE)

#clean it up
local_checklist<-local_checklist%>%
  dplyr::mutate(V1=gsub("Gen. ", "",V1))%>%
  dplyr::mutate(V1=gsub("indet. ", "",V1))%>%
  dplyr::mutate(V1=gsub("\"", "",V1))%>%
  dplyr::mutate(V1=gsub("c.f. ", "",V1))%>%
  dplyr::mutate(V1=gsub("cf.", "",V1))%>%
  dplyr::rename(Scientific_name=V1)
}

#separate genus from species
local_genuses<-local_checklist%>%
  separate(Scientific_name, into = c("listed_genus", "listed_species"), sep=" ",remove = FALSE)

local_checklist<-unique(full_join(local_checklist, local_genuses, by="Scientific_name"))
##get taxonomy for everything in the comprehensive Galapagos species checklist from taxonomizr-----

#need the taxonomizr database to be already prepared
checklist_taxIDs<-getId(local_checklist$Scientific_name,'accessionTaxa.sql')
#print(checklist_taxIDs)

local_checklist$taxID<-checklist_taxIDs

#get the taxonomy for the checklist TaxIDs

checklist_levels<-getTaxonomy(checklist_taxIDs,'accessionTaxa.sql')
#print(checklist_levels)
checklist_levels<-as.data.frame(checklist_levels)%>%
  dplyr::rename(Scientific_name=species)
  

local_checklist<-full_join(local_checklist, checklist_levels, by="Scientific_name")

# If the database=="global", check if the preferred_name is in the local_checklist, if it is, then leave it.  If not, get the global_taxas entry for that ID and check if the genus is in the local_checklist$genus column, go up the columns until a match.


#put the global names and lca reassignments in a new dataframe
global_preferred<-best_ID_combined%>%
  filter(database=="global")

checklist_best_ID_combined<-global_preferred%>%
  select(ID,Sequence,TAXID,SCIENTIFIC_NAME,global_pctid,preferred_name,preferred_pctid)

global_preferred_levels<-getTaxonomy(global_preferred$TAXID,'accessionTaxa.sql')
global_preferred_IDs<-global_preferred$ID

global_preferred_levels<-as.data.frame(global_preferred_levels)%>%
  mutate(ID=global_preferred_IDs)

global_preferred_levels<-full_join(global_preferred_levels,checklist_best_ID_combined)

global_preferred_levels$local_relative<-NA
#stop at family and everything beyond that level becomes NA
global_preferred_levels<-global_preferred_levels%>%
  mutate(local_relative=
           if_else(species%in%local_checklist$Scientific_name & is.na(species)==FALSE,
                   species, 
                   if_else(genus%in%local_checklist$listed_genus& is.na(genus)==FALSE,
                           genus,
                           if_else(genus%in%local_checklist$genus& is.na(genus)==FALSE,
                                   genus,
                                   if_else(family%in%local_checklist$family& is.na(family)==FALSE,
                                           family,
                                           NA)))))


local_relative_df<-global_preferred_levels%>%
  select(ID,local_relative)


#put the local_relative into the best_ID_combined dataframe
best_ID_combined<-full_join(best_ID_combined,local_relative_df, by="ID")
best_ID_combined<-best_ID_combined%>%
  mutate(preferred_name=
           if_else(database=="global",
                   local_relative,
                   preferred_name))


off_target_global_preferred_list<-global_preferred_levels

off_target_global_preferred_list<-off_target_global_preferred_list%>%
  mutate(preferred_name=local_relative)%>%
  mutate(database="global")

#print global_preferred_assignments to a csv file to maybe use as supplement.
write.csv(off_target_global_preferred_list, paste0(Primer,"_output/global_preferred_assignments_after_local_db_and_checklist_LCA.csv"))


#re-summarize changes after LCA of global vs. local.  Now all should be local.-----

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
#545

#Percent of all ASVs assigned to a taxon after combining vsearch and obitools
(assigned/nrow(best_ID_combined))*100
#37.89986

#compared to just obitools
obi_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
obi_assigned
#542

#this is a decrease in assignment of 526 crustacean ASVs (because of off-target amplification of insects and bryozoans, and cnidarians)
assigned - obi_assigned
#3


lca_vsearch_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$lca_name) ==FALSE, ])
lca_vsearch_assigned
#188

#After comparing the taxonomic assignments of global and local to get preferred names----

#count the number of times the global database was used for the preferred assignment
lca_global_v_local_assigned<-nrow(best_ID_combined[best_ID_combined$database == 'lca_global_v_local'& is.na(best_ID_combined$preferred_name) ==FALSE, ])

lca_global_v_local_assigned
#45

globally_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$preferred_name) ==FALSE& best_ID_combined$preferred_name==best_ID_combined$SCIENTIFIC_NAME, ])

globally_assigned
#357


#count the number of times the local database was used for the preferred assignment
locally_assigned<-nrow(best_ID_combined[is.na(best_ID_combined$preferred_name) ==FALSE& best_ID_combined$preferred_name==best_ID_combined$lca_name, ])
locally_assigned

#143

#percentage of all ASVs
(locally_assigned/nrow(best_ID_combined))*100
#9.944367%

#percentage of identified sequences
(locally_assigned/assigned)*100
#26.23853%


#number of occurences where local assignment changed the existing global assignment
global_to_local_assigment<-nrow(best_ID_combined[best_ID_combined$database == 'local'& is.na(best_ID_combined$SCIENTIFIC_NAME) ==FALSE, ])
global_to_local_assigment
#139 taxa updated from existing global assignments


#number of global taxa that got updated with local checklist
global_to_local_relative<-nrow(best_ID_combined[best_ID_combined$preferred_name == best_ID_combined$local_relative & best_ID_combined$preferred_name!= best_ID_combined$SCIENTIFIC_NAME &is.na(best_ID_combined$preferred_name)==FALSE, ])
global_to_local_relative




#Make new motu (taxa) table for MetabaR ----

## MOTUs characteristics table
library(metabaR)

final_names<-best_ID_combined$preferred_name

final_taxaId<-getId(final_names,'accessionTaxa.sql')
#print(final_taxaId)

final_taxa<-getTaxonomy(final_taxaId,'accessionTaxa.sql')
#print(final_taxa)
class(final_taxa)

#make a new dataframe with the final taxa

final_taxa<-as.data.frame(final_taxa)
IDs<-best_ID_combined$ID

rownames(final_taxa)<-IDs

nrow(final_taxa)
#1438

#select only the columns I want for metabaR and rename them (lowercase sequence) taking out "preferred" for brevity.
final_taxa$pct_id<-best_ID_combined$preferred_pctid

final_taxa$Scientific_name<-final_names
final_taxa$taxID<-final_taxaId
final_taxa$ASV<-IDs
final_taxa$sequence<-best_ID_combined$Sequence
final_taxa$database<-best_ID_combined$database
final_taxa$COUNT<-best_ID_combined$COUNT


final_taxa<-final_taxa%>%
  mutate(pct_id=
           if_else(is.na(Scientific_name)==TRUE,
                   NA,
                   pct_id))%>%
  mutate(database=
           if_else(is.na(Scientific_name)==TRUE,
                   NA,
                   database))


class(final_taxa)
write.csv(final_taxa,paste0("../06_local_vs_global_results/",Primer,"_Menu_ready_for_MetabaR.csv"))

write.csv(best_ID_combined,paste0("../06_local_vs_global_results/",Primer,"_best_ID_combined.csv"))

print(paste0("Finished cross validating taxonomic assignments for ",Primer," with Local advantage set to ",Local_advantage))



post_summary_best_ID_combined<-best_ID_combined%>%
  summarise(n_global_v_checklist=n_distinct(local_relative),
            n_combined_IDed_taxa= n_distinct(preferred_name))
            
  
post_summary_best_ID_combined
globalvlocal<-best_ID_combined%>%
  filter(database=="lca_global_v_local")
summaryglobalvlocal<-globalvlocal%>%
  summarise(n_global_v_local= n_distinct(preferred_name))
summaryglobalvlocal

finallocal<-best_ID_combined%>%
  filter(preferred_name==lca_name)
summaryfinallocal<-finallocal%>%
  summarise(n_finallocal= n_distinct(preferred_name))
summaryfinallocal

finalglobal<-best_ID_combined%>%
  filter(preferred_name==SCIENTIFIC_NAME)
summaryfinalglobal<-finalglobal%>%
  summarise(n_finalglobal= n_distinct(preferred_name))
summaryfinalglobal

