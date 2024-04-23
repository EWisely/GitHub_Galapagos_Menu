#Eldridge Wisely
#4-13-24
#compare the percent ID columns of the global and local database hits and pick the highest.  Use global assignment if values are equal. Create a new table to be used in MetabaR as the MOTUS table.
#Then compare to local and global databases to see how many identified species from the samples have been found locally or regionally before.
#tag the new MOTUS table with new columns for Local and Regional

library(tidyverse)
library(prettyunits)
library(taxonomizr)

#Load global obitools database results and local usearch database results ----

#global EMBL database obitools results
obi_fish95<-readr::read_delim("../03_obitools_results/MiFish_Menu_95_named.tab")

colnames(obi_fish95)
obi_fish95<-obi_fish95%>%
  select(c(ID,TAXID,SCIENTIFIC_NAME,BEST_IDENTITY))

#local database usearch results
usearch_fish95<-readr::read_delim("../04_usearch_obitools_results/Menu_Galapagos_fish_usearch_results.txt", col_names = FALSE)

ufish_95<-separate(data = usearch_fish95, col = X1, sep = " ",into = c("ID","COUNT","UMI"))
ufish_95<-separate(data = ufish_95, col = X2, sep = " ",into = c("NCBI_ACC","GENUS","SPECIES","DEFINITION"))
ufish_95<-rename(ufish_95, usearch_pctid=X3, Sequence=X9)

ufish_95$Scientific_Name<-paste(ufish_95$GENUS, ufish_95$SPECIES, sep=" ")

colnames(ufish_95)
ufish_95<-ufish_95%>%
  select(c(ID,NCBI_ACC,Scientific_Name,usearch_pctid,Sequence))


#Vsearch with comprehensive OBIS_GBIF_Darwin list----
lca_vsearch_MiFish<-readr::read_delim("../05_vsearch_obitools_results/lca_top5_comprehensive_galapagos_results.txt", col_names = c("ID","sintax"), delim = '\t')

#parse the sintax taxonomy column into 7 new columns.

lca_vsearch_MiFish<-lca_vsearch_MiFish%>%
  separate_wider_delim(sintax, delim = ',', names = c("Domain","Phylum","Class","Order","Family","Genus","Species"),too_few = "align_start")

lca_vsearch_MiFish$Domain<-str_remove_all(lca_vsearch_MiFish$Domain,pattern = "d:")
lca_vsearch_MiFish$Phylum<-str_remove_all(lca_vsearch_MiFish$Phylum, pattern = "p:")
lca_vsearch_MiFish$Class<-str_remove_all(lca_vsearch_MiFish$Class, pattern = "c:")
lca_vsearch_MiFish$Order<-str_remove_all(lca_vsearch_MiFish$Order, pattern = "o:")
lca_vsearch_MiFish$Family<-str_remove_all(lca_vsearch_MiFish$Family, pattern = "f:")
lca_vsearch_MiFish$Genus<-str_remove_all(lca_vsearch_MiFish$Genus, pattern = "g:")
lca_vsearch_MiFish$Species<-str_remove_all(lca_vsearch_MiFish$Species, pattern = "s:")
lca_vsearch_MiFish$Species<-str_replace_all(lca_vsearch_MiFish$Species, pattern = "_", replacement = " ")
  
#read in userout for vsearch  
#pct_vsearch_MiFish<-readr::read_delim("../05_vsearch_obitools_results/userout_top5_comprehensive_galapagos_results.txt", col_names = c("ID","ACC","PCTID","ALNLEN","MISM","OPENS")) #This code works for the results when compared to the un-dereplicated, not-TAXID'ed file.  But when using the same sintax file as I used in the lca command, I'll have to do some more wrangling.

pct_vsearch_MiFish<-readr::read_delim("../05_vsearch_obitools_results/userout_top5_comprehensive_galapagos_results.txt", col_names =FALSE, delim = '\t')

pct_vsearch_MiFish<-pct_vsearch_MiFish%>%
  separate(col = X2, sep = ";",into = c("ACC","sintax"))%>%
  dplyr::rename(ID=X1, vsearch_pctid=X3, vsearch_alnlen=X4, vsearch_mism=X5, vsearch_opens=X6)

#Import the obitools ASVs file so that the sequences of each can be used in making the MOTUS table for MetabaR later.

query_ASVs<-read_delim("../05_vsearch_obitools_results/upper_MiFish_Menu_95_named_cleared_tags.tab",col_names = FALSE)

query_ASVs<-query_ASVs%>%
  separate(X1, into = c("ID","Count"), sep=" ")%>%
  dplyr::rename(Sequence=X2)%>%
  select(c("ID","Sequence"))


#combine the LCA taxonomy with the percent_id column (because in vsearch I specified only the top matches so the multiples all share the same pctid with each other, so we can pick any of them)
MiFish_vsearch_results<- left_join(lca_vsearch_MiFish, pct_vsearch_MiFish, by ="ID", multiple ="any")

MiFish_vsearch_results_seq<-full_join(MiFish_vsearch_results, query_ASVs, by="ID")


#v2_vsearch using the in-silico PCR'ed mito file added to ssu and 12S without insilico PCR-----

v2_lca_vsearch_MiFish<-readr::read_delim("../05_vsearch_obitools_results/v2_lca_top5_comprehensive_galapagos_results.txt", col_names = c("ID","sintax"), delim = '\t')

#parse the sintax taxonomy column into 7 new columns.

v2_lca_vsearch_MiFish<-v2_lca_vsearch_MiFish%>%
  separate_wider_delim(sintax, delim = ',', names = c("Domain","Phylum","Class","Order","Family","Genus","Species"),too_few = "align_start")

v2_lca_vsearch_MiFish$Domain<-str_remove_all(v2_lca_vsearch_MiFish$Domain,pattern = "d:")
v2_lca_vsearch_MiFish$Phylum<-str_remove_all(v2_lca_vsearch_MiFish$Phylum, pattern = "p:")
v2_lca_vsearch_MiFish$Class<-str_remove_all(v2_lca_vsearch_MiFish$Class, pattern = "c:")
v2_lca_vsearch_MiFish$Order<-str_remove_all(v2_lca_vsearch_MiFish$Order, pattern = "o:")
v2_lca_vsearch_MiFish$Family<-str_remove_all(v2_lca_vsearch_MiFish$Family, pattern = "f:")
v2_lca_vsearch_MiFish$Genus<-str_remove_all(v2_lca_vsearch_MiFish$Genus, pattern = "g:")
v2_lca_vsearch_MiFish$Species<-str_remove_all(v2_lca_vsearch_MiFish$Species, pattern = "s:")
v2_lca_vsearch_MiFish$Species<-str_replace_all(v2_lca_vsearch_MiFish$Species, pattern = "_", replacement = " ")

#read in userout for vsearch  

v2_pct_vsearch_MiFish<-readr::read_delim("../05_vsearch_obitools_results/v2_userout_top5_comprehensive_galapagos_results.txt", col_names =FALSE, delim = '\t')

v2_pct_vsearch_MiFish<-v2_pct_vsearch_MiFish%>%
  separate(col = X2, sep = ";",into = c("ACC","sintax"))%>%
  dplyr::rename(ID=X1, vsearch_pctid=X3, vsearch_alnlen=X4, vsearch_mism=X5, vsearch_opens=X6)


#combine the LCA taxonomy with the percent_id column (because in vsearch I specified only the top matches so the multiples all share the same pctid with each other, so we can pick any of them)
v2_MiFish_vsearch_results<- left_join(v2_lca_vsearch_MiFish, v2_pct_vsearch_MiFish, by ="ID", multiple ="any")

v2_MiFish_vsearch_results_seq<-full_join(v2_MiFish_vsearch_results, query_ASVs, by="ID")

v2_MiFish_vsearch_results_all_matches<-full_join(v2_lca_vsearch_MiFish,v2_pct_vsearch_MiFish, by ="ID")

v2_chordata_only<-v2_MiFish_vsearch_results_all_matches%>%
  dplyr::filter(is.na(Class))%>%
  dplyr::filter(Domain =="Eukaryota")
#184 observations, and they all are mixing fish with one of two ray sequences from the same study from 2014 using ion torrent data.

#v3_vsearch using the in-silico PCR'ed mito file added to ssu and 12S without insilico PCR-----

v3_lca_vsearch_MiFish<-readr::read_delim("../05_vsearch_obitools_results/no_prob_v2_lca_top5_comprehensive_galapagos_results.txt", col_names = c("ID","sintax"), delim = '\t')

#parse the sintax taxonomy column into 7 new columns.

v3_lca_vsearch_MiFish<-v3_lca_vsearch_MiFish%>%
  separate_wider_delim(sintax, delim = ',', names = c("Domain","Phylum","Class","Order","Family","Genus","Species"),too_few = "align_start")

v3_lca_vsearch_MiFish$Domain<-str_remove_all(v3_lca_vsearch_MiFish$Domain,pattern = "d:")
v3_lca_vsearch_MiFish$Phylum<-str_remove_all(v3_lca_vsearch_MiFish$Phylum, pattern = "p:")
v3_lca_vsearch_MiFish$Class<-str_remove_all(v3_lca_vsearch_MiFish$Class, pattern = "c:")
v3_lca_vsearch_MiFish$Order<-str_remove_all(v3_lca_vsearch_MiFish$Order, pattern = "o:")
v3_lca_vsearch_MiFish$Family<-str_remove_all(v3_lca_vsearch_MiFish$Family, pattern = "f:")
v3_lca_vsearch_MiFish$Genus<-str_remove_all(v3_lca_vsearch_MiFish$Genus, pattern = "g:")
v3_lca_vsearch_MiFish$Species<-str_remove_all(v3_lca_vsearch_MiFish$Species, pattern = "s:")
v3_lca_vsearch_MiFish$Species<-str_replace_all(v3_lca_vsearch_MiFish$Species, pattern = "_", replacement = " ")

#read in userout for vsearch  

v3_pct_vsearch_MiFish<-readr::read_delim("../05_vsearch_obitools_results/no_prob_v2_userout_top5_comprehensive_galapagos_results.txt", col_names =FALSE, delim = '\t')

v3_pct_vsearch_MiFish<-v3_pct_vsearch_MiFish%>%
  separate(col = X2, sep = ";",into = c("ACC","sintax"))%>%
  dplyr::rename(ID=X1, vsearch_pctid=X3, vsearch_alnlen=X4, vsearch_mism=X5, vsearch_opens=X6)


#combine the LCA taxonomy with the percent_id column (because in vsearch I specified only the top matches so the multiples all share the same pctid with each other, so we can pick any of them)
v3_MiFish_vsearch_results<- left_join(v3_lca_vsearch_MiFish, v3_pct_vsearch_MiFish, by ="ID", multiple ="any")

v3_MiFish_vsearch_results_seq<-full_join(v3_MiFish_vsearch_results, query_ASVs, by="ID")

v3_MiFish_vsearch_results_all_matches<-full_join(v3_lca_vsearch_MiFish,v3_pct_vsearch_MiFish, by ="ID")

v3_chordata_only<-v3_MiFish_vsearch_results_all_matches%>%
  dplyr::filter(is.na(Class))%>%
  dplyr::filter(Domain =="Eukaryota")
#0 observations LCA'ed to just Chordata!  YAY!







#USEARCH combine results of obitools and usearch_global no lca----

#join by ID column

fish_combined<- full_join(ufish_95, obi_fish95, by ="ID")


#Compare global and local database assignments and pick the best assignment for each ASV----

#if BEST_IDENTITY*100 is greater than or equal to usearch_pctid then mutate fish_combined$consensus_Sciname is SCIENTIFIC_NAME.  and consensus_taxID is TAXID.
#if usearch_pctid is greater than BEST_IDENTITY*100, then mutate fish_combined$consensus_Sciname is Scientific_Name


best_ID_fish_combined<-fish_combined%>%
  mutate(local_pctid=
           if_else(usearch_pctid=="NA",
                   0,
                   usearch_pctid))%>%
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
                   Scientific_Name, 
                   SCIENTIFIC_NAME))%>%
  mutate(database=
           if_else(local_pctid > global_pctid,
                   "local",
                   "global"))

  #select(ID,preferred_pctid,preferred_name,database)

#summarize changes to taxonomic classifications----


#count how many got assigned to local vs. global and number of unassigned.

#total ASVs
total_ASVs=nrow(obi_fish95)
total_ASVs
#to make sure we didn't lose any
final_ASVs=nrow(best_ID_fish_combined)
final_ASVs
#1438

#Updated number of total ASVs assigned to a taxon
assigned<-nrow(best_ID_fish_combined[is.na(best_ID_fish_combined$database) ==FALSE, ])
assigned
#688

#Percent of all ASVs assigned to a taxon after combining usearch and obitools
(assigned/nrow(best_ID_fish_combined))*100
#47.84423

#compared to just obitools
obi_assigned<-nrow(best_ID_fish_combined[is.na(best_ID_fish_combined$SCIENTIFIC_NAME) ==FALSE, ])
obi_assigned
#542

#this is an increase in assignment of 146 ASVs
assigned - obi_assigned
#146


usearch_assigned<-nrow(best_ID_fish_combined[is.na(best_ID_fish_combined$Scientific_Name) ==FALSE, ])
usearch_assigned
#688

#After combining the taxonomic assignments of global and local ----

#count the number of times the global database was used for the preferred assignment
globally_assigned<-nrow(best_ID_fish_combined[best_ID_fish_combined$database == 'global'& is.na(best_ID_fish_combined$database) ==FALSE, ])

globally_assigned
#333

#percentage of all ASVs assigned to the global database
(globally_assigned/nrow(best_ID_fish_combined))*100
#23.15716%
#percentage of identfied ASVs assigned to the global database
(globally_assigned/assigned)*100
#48.40116%

#count the number of times the local database was used for the preferred assignment
locally_assigned<-nrow(best_ID_fish_combined[best_ID_fish_combined$database == 'local'& is.na(best_ID_fish_combined$database) ==FALSE, ])
locally_assigned
#355
#percentage of all ASVs
(locally_assigned/nrow(best_ID_fish_combined))*100
#24.68707%

#percentage of identified sequences
(locally_assigned/assigned)*100
#51.59884%


#number of occurences where local assignment changed the existing global assignment
global_to_local_assigment<-nrow(best_ID_fish_combined[best_ID_fish_combined$database == 'local'& is.na(best_ID_fish_combined$SCIENTIFIC_NAME) ==FALSE, ])
global_to_local_assigment
#209 taxa updated from existing global assignments

# Make new motu (taxa) table for MetabaR ----

## MOTUs characteristics table
library(metabaR)
library(taxonomizr)
prepareDatabase('accessionTaxa.sql')

acc<-best_ID_fish_combined$NCBI_ACC[!is.na(best_ID_fish_combined$NCBI_ACC)]
print(acc)


taxaId<-accessionToTaxa(accessions = acc, sqlFile ="accessionTaxa.sql", version = "version")
getTaxonomy(taxaId,'accessionTaxa.sql')
#put this info into best_ID_fish_combined as more columns matching on acc 
#mutate a new column with preferred_TAXID if_else local>global etc.
#select only the columns I want for metabaR and rename them (lowercase sequence) taking out "preferred" for brevity.





#VSEARCH combine results of obitools and vsearch_global with lca------

#join by ID column

lca_fish_combined<- full_join(v3_MiFish_vsearch_results_seq, obi_fish95, by ="ID")

lca_fish_combined<-lca_fish_combined%>%
  mutate(lca_name = coalesce(Species,Genus,Family,Order,Class,Phylum,Domain))


#Compare global and local database assignments and pick the best assignment for each ASV----

#if BEST_IDENTITY*100 is greater than vsearch_pctid then mutate fish_combined$consensus_Sciname is SCIENTIFIC_NAME.  and consensus_taxID is TAXID.
#if vsearch_pctid is greater than or equal to BEST_IDENTITY*100, then mutate fish_combined$consensus_Sciname is lca_name


lca_best_ID_fish_combined<-lca_fish_combined%>%
  mutate(local_pctid=
           if_else(vsearch_pctid=="NA",
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
total_ASVs=nrow(obi_fish95)
total_ASVs
#to make sure we didn't lose any
final_ASVs=nrow(lca_best_ID_fish_combined)
final_ASVs
#1438

#Updated number of total ASVs assigned to a taxon
assigned<-nrow(lca_best_ID_fish_combined[is.na(lca_best_ID_fish_combined$database) ==FALSE, ])
assigned
#790

#Percent of all ASVs assigned to a taxon after combining usearch and obitools
(assigned/nrow(lca_best_ID_fish_combined))*100
#54.93741

#compared to just obitools
obi_assigned<-nrow(lca_best_ID_fish_combined[is.na(lca_best_ID_fish_combined$SCIENTIFIC_NAME) ==FALSE, ])
obi_assigned
#542

#this is an increase in assignment of 248 ASVs
assigned - obi_assigned
#248


lca_vsearch_assigned<-nrow(lca_best_ID_fish_combined[is.na(lca_best_ID_fish_combined$lca_name) ==FALSE, ])
lca_vsearch_assigned
#790

#After comparing the taxonomic assignments of global and local to get preferred names----

#count the number of times the global database was used for the preferred assignment
globally_assigned<-nrow(lca_best_ID_fish_combined[lca_best_ID_fish_combined$database == 'global'& is.na(lca_best_ID_fish_combined$database) ==FALSE, ])

globally_assigned
#26

#put these in a new dataframe
global_preferred_assignments<-lca_best_ID_fish_combined%>%
  filter(database=="global")

#print global_preferred_assignments to a csv file to maybe use as supplement.
write.csv(global_preferred_assignments, "MiFish_output/global_preferred_assignments_before_final_LCA.csv")

#percentage of all ASVs assigned to the global database
(globally_assigned/nrow(lca_best_ID_fish_combined))*100
#1.808067%
#percentage of identfied ASVs assigned to the global database
(globally_assigned/assigned)*100
#3.291139%

#count the number of times the local database was used for the preferred assignment
locally_assigned<-nrow(lca_best_ID_fish_combined[lca_best_ID_fish_combined$database == 'local'& is.na(lca_best_ID_fish_combined$database) ==FALSE, ])
locally_assigned
#764
#percentage of all ASVs
(locally_assigned/nrow(lca_best_ID_fish_combined))*100
#53.129%

#percentage of identified sequences
(locally_assigned/assigned)*100
#96.70886%


#number of occurences where local assignment changed the existing global assignment
global_to_local_assigment<-nrow(lca_best_ID_fish_combined[lca_best_ID_fish_combined$database == 'local'& is.na(lca_best_ID_fish_combined$SCIENTIFIC_NAME) ==FALSE, ])
global_to_local_assigment
#516 taxa updated from existing global assignments


#When global is the database instead of local, use the LCA between global and local as preferred name, and the lower pctid (local), as the preferred_pctid----

#LCA of global and local ID'ed taxa-----
#need the taxonomizr database to be already prepared
meta_lca_best_ID_fish_combined<-lca_best_ID_fish_combined%>%
  separate(sintax, into=c("sintax-to-genus", "vsearch_species"),sep = "s:", remove = FALSE)

v_species<-gsub("_"," ",meta_lca_best_ID_fish_combined$vsearch_species)
meta_lca_best_ID_fish_combined$vsearch_species<-v_species


taxaId<-getId(v_species,'accessionTaxa.sql')
print(taxaId)

meta_lca_best_ID_fish_combined$v_taxID<-taxaId

#get the taxonomy for the vsearch TaxIDs (even though we have the sintax format already, this is getting it ready to do LCA in this program)

local_levels<-getTaxonomy(taxaId,'accessionTaxa.sql')
print(local_levels)

global_taxaId<-meta_lca_best_ID_fish_combined$TAXID
global_levels<-getTaxonomy(global_taxaId,'accessionTaxa.sql')
print(global_levels)

global_taxa<-meta_lca_best_ID_fish_combined$SCIENTIFIC_NAME


#condenseTaxa(taxa)
##   superkingdom phylum     class      order family genus species
## 1 "Eukaryota"  "Chordata" "Mammalia" NA    NA     NA    NA
#This function can also be fed a large number of grouped hits, e.g. BLAST hits for high throughput sequencing reads after filtering for the best hits for each read, and output a condensed taxonomy for each grouping:


IDs<-meta_lca_best_ID_fish_combined$ID

local_taxas<-as.data.frame(local_levels, row.names = paste(IDs,"_local"))
global_taxas<-as.data.frame(global_levels, row.names = paste(IDs, "_global"))
taxas<-rbind(local_taxas, global_taxas)

doubleIDs<-append(IDs,IDs)
condensed_lca_global_vs_local<-condenseTaxa(taxas, groupings = doubleIDs)

condensed_lca_global_vs_local<-rownames_to_column(as.data.frame(condensed_lca_global_vs_local), var = "ID")

condensed_lca_global_vs_local<-condensed_lca_global_vs_local%>%
  mutate(global_v_local_lca_name = coalesce(species,genus,family,order,class,phylum,superkingdom))


#add global_v_local_lca_name column to lca_best_ID_fish_combined

#assign global_v_local_lca_name to preferred name when global scores higher than local



#re-summarize changes after LCA of global vs. local.  Now all should be local.-----


# Make new motu (taxa) table for MetabaR ----

## MOTUs characteristics table
library(metabaR)
library(taxonomizr)
prepareDatabase('accessionTaxa.sql')

acc<-best_ID_fish_combined$NCBI_ACC[!is.na(best_ID_fish_combined$NCBI_ACC)]
print(acc)


taxaId<-accessionToTaxa(accessions = acc, sqlFile ="accessionTaxa.sql", version = "version")
getTaxonomy(taxaId,'accessionTaxa.sql')
#put this info into best_ID_fish_combined as more columns matching on acc 
#mutate a new column with preferred_TAXID if_else local>global etc.
#select only the columns I want for metabaR and rename them (lowercase sequence) taking out "preferred" for brevity.








#create a TRUE or FALSE vector if the number two is in list A answer by julius-vainora on https://stackoverflow.com/questions/53086053/how-to-check-if-a-list-contains-a-certain-element-in-r

#sapply(A, `%in%`, x = 2)

#ooh!  neat!
#https://ropensci.org/blog/2017/01/25/obis/


##### This section under construction #####


# #### check occurence data of remaining identified taxa that aren't local or regional ####
# #https://docs.ropensci.org/rgbif/articles/downloading_a_long_species_list.html
# 
# 
# #import species lists
# species_by_line<-"../03_obitools_results/Menu_1-24_obicrust_results_95.txt"
# 
# long_checklist <- readr::read_tsv(species_by_line, col_names = "Species")
# 
# 
# # match the names 
# gbif_taxon_keys <- long_checklist %>% 
#   head(1000) %>% # only first 1000 names 
#   name_backbone_checklist() %>% # match to backbone 
#   filter(!matchType == "NONE") %>% # get matched names
#   pull(usageKey) 
# 
# # gbif_taxon_keys should be a long vector like this c(2977832,2977901,2977966,2977835,2977863)
# 
# # download the data
# occ_download(
#   pred_in("taxonKey", gbif_taxon_keys), # important to use pred_in
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SIMPLE_CSV"
# )
# 
# 
# 
# 
# 
# #MiFish Menu 95% ID 1-24 EMBL obitools assignment List Citation:
# # <<gbif download>>
# # Your download is being processed by GBIF:
# #   https://www.gbif.org/occurrence/download/0152701-240321170329656
# # Most downloads finish within 15 min.
# # Check status with
# # occ_download_wait('0152701-240321170329656')
# # After it finishes, use
# # d <- occ_download_get('0152701-240321170329656') %>%
# #   occ_download_import()
# # to retrieve your download.
# # Download Info:
# #   Username: eldridgewisely
# # E-mail: eldridge.wisely@gmail.com
# # Format: SIMPLE_CSV
# # Download key: 0152701-240321170329656
# # Created: 2024-04-12T22:11:41.739+00:00
# # Citation Info:  
# #   Please always cite the download DOI when using this data.
# # https://www.gbif.org/citation-guidelines
# # DOI: 10.15468/dl.jfptf8
# # Citation:
# #   GBIF Occurrence Download https://doi.org/10.15468/dl.jfptf8 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-04-12
# 
# obi95fish <- occ_download_get('0152701-240321170329656') %>%
#   occ_download_import()
# 
# 
# #BerryCrust Menu 95% ID 1-24 EMBL obitools assignment List Citation:
# #<<gbif download>>
# # Your download is being processed by GBIF:
# #   https://www.gbif.org/occurrence/download/0152722-240321170329656
# # Most downloads finish within 15 min.
# # Check status with
# # occ_download_wait('0152722-240321170329656')
# # After it finishes, use
# # d <- occ_download_get('0152722-240321170329656') %>%
# #   occ_download_import()
# # to retrieve your download.
# # Download Info:
# #   Username: eldridgewisely
# # E-mail: eldridge.wisely@gmail.com
# # Format: SIMPLE_CSV
# # Download key: 0152722-240321170329656
# # Created: 2024-04-12T22:14:53.674+00:00
# # Citation Info:  
# #   Please always cite the download DOI when using this data.
# # https://www.gbif.org/citation-guidelines
# # DOI: 10.15468/dl.u38qqq
# # Citation:
# #   GBIF Occurrence Download https://doi.org/10.15468/dl.u38qqq Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-04-12
# 
# #obi95crustaceans <- occ_download_get('0152722-240321170329656') %>%
# #     occ_download_import()
# 
# #https://search.r-project.org/CRAN/refmans/rgbif/html/occ_data.html
# 
# #Galapagos Fish List Citation:
# # #Download Info:
# # Username: eldridgewisely
# # E-mail: eldridge.wisely@gmail.com
# # Format: SIMPLE_CSV
# # Download key: 0145620-240321170329656
# # Created: 2024-04-11T22:21:13.971+00:00
# # Citation Info:  
# #   Please always cite the download DOI when using this data.
# # https://www.gbif.org/citation-guidelines
# # DOI: 10.15468/dl.aw6e6v
# # Citation:
# #   GBIF Occurrence Download https://doi.org/10.15468/dl.aw6e6v Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-04-11
# 
# #Galapagos_fish_gbif<-occ_download_get('0145620-240321170329656') %>%
# #  occ_download_import()
# 
# #Galapagos Crustaceans List Citation:
# # Username: eldridgewisely
# # E-mail: eldridge.wisely@gmail.com
# # Format: SIMPLE_CSV
# # Download key: 0145629-240321170329656
# # Created: 2024-04-11T22:24:26.762+00:00
# # Citation Info:  
# #   Please always cite the download DOI when using this data.
# # https://www.gbif.org/citation-guidelines
# # DOI: 10.15468/dl.n55tv9
# # Citation:
# #   GBIF Occurrence Download https://doi.org/10.15468/dl.n55tv9 Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-04-11
# 
# 
# #Galapagos_crustaceans_gbif<-occ_download_get('0145629-240321170329656') %>%
# #  occ_download_import()
# 
# 
# #https://search.r-project.org/CRAN/refmans/rgbif/html/occ_data.html
# 
# # explore the Galapagos Fish List Dataset
# occ_count(datasetKey='d8cd16ba-bb74-4420-821e-083f2bac17c2') # from the Galapagos Fish List Dataset
# occ_count(datasetKey='d8cd16ba-bb74-4420-821e-083f2bac17c2',hasCoordinate=TRUE) # all with coordinates
# occ_count(datasetKey='d8cd16ba-bb74-4420-821e-083f2bac17c2',basisOfRecord='OBSERVATION') # with this basis of record
# occ_count(taxonKey=2435099, hasCoordinate=TRUE) # with coordinates
# occ_count(datasetKey='d8cd16ba-bb74-4420-821e-083f2bac17c2',country="EC") # from Ecuador 
# occ_count(datasetKey='d8cd16ba-bb74-4420-821e-083f2bac17c2') # from the Galapagos Fish List Dataset
# occ_count_country("EC")
# 
# # Search on latitidue and longitude
# occ_data(decimalLatitude=0, decimalLongitude=-90, limit = 10)
# ## or using bounding box, converted to WKT internally
# 
# occ_data(geometry=c(-95,-85,-1,1), datasetKey='d8cd16ba-bb74-4420-821e-083f2bac17c2')
# 
# 
# #by administrative boundary of Galapagos ECU.9_1
# occ_data(gadmGid = "ECU.9_1", datasetKey='d8cd16ba-bb74-4420-821e-083f2bac17c2')
# 
# 
# occ_data(geometry=c(-95,-85,-1,1), gadmGid = "ECU.9_1", datasetKey='d8cd16ba-bb74-4420-821e-083f2bac17c2')
# 
# 
# 
# 
