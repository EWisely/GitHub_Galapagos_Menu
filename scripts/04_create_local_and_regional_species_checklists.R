#Eldridge Wisely
#
#This script searches for local (Galapagos) and regional (TEP) checklists for fish and crustacean species from GBIF (Global Biodiversity Information Facility ) and OBIS (Ocean Biodiversity Information System https://obis.org), combines them with curated species checklists downloaded from the Charles Darwin Research Station Natural History Collections database (https://datazone.darwinfoundation.org/en/checklist/checklists-archive) in preparation for downloading with CRABS program

#devtools::install_github("james-thorson/FishLife")
#devtools::install_github("cfree14/freeR")

library(usethis)
#usethis::edit_r_environ()
#input GBIF login information
library(rgbif)
library(rfishbase)
library(dplyr)
library(readr)  
library(worrms)
library(taxize)
library(robis)
library(tidyverse)
library(taxonomizr)
library(readr)
#library(freeR)
#library(FishLife)


#https://docs.ropensci.org/rgbif/articles/gbif_credentials.html
#usethis::edit_r_environ()


#OBIS----
####Galapagos fish checklist####
obis_galfish<-checklist(c("Agnatha", "Chondrichthyes", "Osteichthyes"),geometry = "POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))")

#### Pull out only the marine fish species list for the local database
obis_galfish_species<-obis_galfish%>%
  filter(taxonRank=="Species" , is_terrestrial==FALSE)%>%
  select(scientificName)%>%
  rename(Species=scientificName)

#### Get the local marine fish taxa list for checking the results against
obis_galfish_taxa<-obis_galfish%>%
  filter(is_terrestrial==FALSE)%>%
  select(scientificName)%>%
  rename(Taxa=scientificName)

####TEP fish checklist####
obis_tepfish<-checklist(c("Agnatha", "Chondrichthyes", "Osteichthyes"),geometry = "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))")

#### Get the regional marine fish taxa list for checking the results against
obis_tepfish_taxa<-obis_tepfish%>%
  filter(is_terrestrial==FALSE)%>%
  select(scientificName)%>%
  rename(Taxa=scientificName)

####Galapagos crustacean checklist####
obis_galcrust<-checklist("Crustacea",geometry = "POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))")

#### Pull out only the marine crustacean species list for the local database
obis_galcrust_species<-obis_galcrust%>%
  filter(taxonRank=="Species" , is_marine==TRUE)%>%
  select(scientificName)%>%
  rename(Species=scientificName)

#### Get the local marine crustacean taxa list for checking the results against
obis_galcrust_taxa<-obis_galcrust%>%
  filter(is_terrestrial==FALSE)%>%
  select(scientificName)%>%
  rename(Taxa=scientificName)

####TEP crustacean checklist####
obis_tepcrust<-checklist("Crustacea",geometry = "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))")

#### Get the regional marine crustacean taxa list for checking the results against
obis_tepcrust_taxa<-obis_tepcrust%>%
  filter(is_terrestrial==FALSE)%>%
  select(scientificName)%>%
  rename(Taxa=scientificName)



#GBIF----

## Download GBIF dataset of occurences within the Galapagos region by geometry ####

##for wkt geometry draw it on this webapp https://wktmap.com----

#for just the Galapagos: POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))
#for the TEP: POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))


#All Galapagos Species (not just fish and crustaceans)

#<<gbif download>>
# Your download is being processed by GBIF:
#   https://www.gbif.org/occurrence/download/0157668-240321170329656
# Most downloads finish within 15 min.
# Check status with
# occ_download_wait('0157668-240321170329656')
# After it finishes, use
# d <- occ_download_get('0157668-240321170329656') %>%
#   occ_download_import()
# to retrieve your download.
# Download Info:
#   Username: eldridgewisely
# E-mail: eldridge.wisely@gmail.com
# Format: SPECIES_LIST
# Download key: 0157668-240321170329656
# Created: 2024-04-13T17:23:24.975+00:00
# Citation Info:  
#   Please always cite the download DOI when using this data.
# https://www.gbif.org/citation-guidelines
# DOI: 10.15468/dl.rqpmhf
# Citation:
#   GBIF Occurrence Download https://doi.org/10.15468/dl.rqpmhf Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-04-13


#Add just backbone key for fish or crustaceans!


#### Find GBIF keys for crustaceans ####
#Find list of Classses within Crustacea subphylum
#Worms (World Register of Marine Species)
worrms::wm_name2id(name = "Crustacea")
#1066
worrms::wm_external(id = 1066, type = "ncbi")
#6657

taxize::worms_downstream(id = 1066, downto = "class")
#id          name  rank
#1    1070   Maxillopoda class
#2    1069  Branchiopoda class
#3  150314 Cephalocarida class
#4    1278      Hexapoda class
#5    1067     Remipedia class
#6    1080      Copepoda class
#7  889925   Hexanauplia class
#8    1071  Malacostraca class
#9    1083 Tantulocarida class
#10  22388   Thecostraca class
#11 845958 Ichthyostraca class
#12   1078     Ostracoda class
crustacean_classes<-taxize::worms_downstream(id = 1066, downto = "class")
crustacean_class_names<-crustacean_classes$name

Crustacean_backbone_keys <- crustacean_class_names %>% 
  name_backbone_checklist() %>% # match to backbone 
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) 

 #### Download Galapagos Crustacean list ####

# occ_download(
#   pred_within("POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))"),
#   pred_in("taxonKey", Crustacean_backbone_keys), # important to use pred_in
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SPECIES_LIST"
# )

# <<gbif download>>
# Your download is being processed by GBIF:
#   https://www.gbif.org/occurrence/download/0172650-240321170329656
# Most downloads finish within 15 min.
# Check status with
# occ_download_wait('0172650-240321170329656')
# After it finishes, use
# d <- occ_download_get('0172650-240321170329656') %>%
#   occ_download_import()
# to retrieve your download.
# Download Info:
#   Username: eldridgewisely
# E-mail: eldridge.wisely@gmail.com
# Format: SPECIES_LIST
# Download key: 0172650-240321170329656
# Created: 2024-04-15T18:25:32.739+00:00
# Citation Info:  
#   Please always cite the download DOI when using this data.
# https://www.gbif.org/citation-guidelines
# DOI: 10.15468/dl.ap833f
# Citation:
#   GBIF Occurrence Download https://doi.org/10.15468/dl.ap833f Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-04-15



##### Download TEP Crustacean list #####


# occ_download(
#   pred_within("POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))"),
#   pred_in("taxonKey", Crustacean_backbone_keys), # important to use pred_in
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SPECIES_LIST"
# )

# <<gbif download>>
# Your download is being processed by GBIF:
#   https://www.gbif.org/occurrence/download/0172684-240321170329656
# Most downloads finish within 15 min.
# Check status with
# occ_download_wait('0172684-240321170329656')
# After it finishes, use
# d <- occ_download_get('0172684-240321170329656') %>%
#   occ_download_import()
# to retrieve your download.
# Download Info:
#   Username: eldridgewisely
# E-mail: eldridge.wisely@gmail.com
# Format: SPECIES_LIST
# Download key: 0172684-240321170329656
# Created: 2024-04-15T18:30:18.884+00:00
# Citation Info:  
#   Please always cite the download DOI when using this data.
# https://www.gbif.org/citation-guidelines
# DOI: 10.15468/dl.d8qknd
# Citation:
#   GBIF Occurrence Download https://doi.org/10.15468/dl.d8qknd Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-04-15

#### Find GBIF keys for fish ####
#"Agnatha", "Chondrichthyes", "Osteichthyes"
#"Osteichtheyes" has a wormsID but no gbif id, so using Actinopterygii instead.

#Find list of Classes within 
#Worms (World Register of Marine Species)
worrms::wm_name2id(name = "Chondrichthyes")
#1517375
taxize::worms_downstream(id = 1517375, downto = "class")
#id           name  rank
#1 10193 Elasmobranchii class
#2 10196    Holocephali class

worrms::wm_name2id(name = "Actinopterygii")
#10194 
taxize::worms_downstream(id = 10194, downto = "class")
#       id        name  rank
#1 1517379 Chondrostei class
#2  293496   Teleostei class

#Osteichthyes = 152352

worrms::wm_name2id(name = "Agnatha")
#1829
taxize::worms_downstream(id = 1829, downto = "class")
# id               name  rank
# 1  10189             Myxini class
# 2 843662       Petromyzonti class
# 3  10190 Cephalaspidomorphi class


fish_classes<-taxize::worms_downstream(id = c(1829,10194,1517375), downto = "class")
fish_class_names<-fish_classes$name

Fish_backbone_keys <- fish_class_names %>% 
  name_backbone_checklist() %>% # match to backbone 
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) 

#### Download Galapagos Fish Checklist from GBIF ####
# occ_download(
#   pred_within("POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))"),
#   pred_in("taxonKey", Fish_backbone_keys), # important to use pred_in
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SPECIES_LIST"
# )
# <<gbif download>>
# Your download is being processed by GBIF:
#   https://www.gbif.org/occurrence/download/0172899-240321170329656
# Most downloads finish within 15 min.
# Check status with
# occ_download_wait('0172899-240321170329656')
# After it finishes, use
# d <- occ_download_get('0172899-240321170329656') %>%
#   occ_download_import()
# to retrieve your download.
# Download Info:
#   Username: eldridgewisely
# E-mail: eldridge.wisely@gmail.com
# Format: SPECIES_LIST
# Download key: 0172899-240321170329656
# Created: 2024-04-15T19:05:40.808+00:00
# Citation Info:  
#   Please always cite the download DOI when using this data.
# https://www.gbif.org/citation-guidelines
# DOI: 10.15468/dl.ktkgat
# Citation:
#   GBIF Occurrence Download https://doi.org/10.15468/dl.ktkgat Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-04-15



##### Download TEP Fish list #####


# occ_download(
#   pred_within("POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))"),
#   pred_in("taxonKey", Fish_backbone_keys), # important to use pred_in
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   format = "SPECIES_LIST"
# )

# <<gbif download>>
# Your download is being processed by GBIF:
#   https://www.gbif.org/occurrence/download/0172907-240321170329656
# Most downloads finish within 15 min.
# Check status with
# occ_download_wait('0172907-240321170329656')
# After it finishes, use
# d <- occ_download_get('0172907-240321170329656') %>%
#   occ_download_import()
# to retrieve your download.
# Download Info:
#   Username: eldridgewisely
# E-mail: eldridge.wisely@gmail.com
# Format: SPECIES_LIST
# Download key: 0172907-240321170329656
# Created: 2024-04-15T19:06:54.789+00:00
# Citation Info:  
#   Please always cite the download DOI when using this data.
# https://www.gbif.org/citation-guidelines
# DOI: 10.15468/dl.qjtewn
# Citation:
#   GBIF Occurrence Download https://doi.org/10.15468/dl.qjtewn Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-04-15


##### Import Galapagos Crustaceans and get the species list ####

GBIF_Gal_crustaceans_list <- occ_download_get('0172650-240321170329656') %>%
   occ_download_import()

GBIF_Gal_crustaceans_species <- GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(Species=species)


#get list including higher taxonomic levels to compare the results file with
GBIF_Gal_crustaceans_taxa1<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(taxa=species)
GBIF_Gal_crustaceans_taxa2<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="GENUS")%>%
  select(genus)%>%
  rename(taxa=genus)
GBIF_Gal_crustaceans_taxa3<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="FAMILY")%>%
  select(family)%>%
  rename(taxa=family)
GBIF_Gal_crustaceans_taxa4<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="ORDER")%>%
  select(order)%>%
  rename(taxa=order)
GBIF_Gal_crustaceans_taxa5<-GBIF_Gal_crustaceans_list%>% 
  filter(taxonRank=="CLASS")%>%
  select(class)%>%
  rename(taxa=class)
GBIF_Gal_crustaceans_taxa<-unique(rbind(GBIF_Gal_crustaceans_taxa1,GBIF_Gal_crustaceans_taxa2,GBIF_Gal_crustaceans_taxa3,GBIF_Gal_crustaceans_taxa4, GBIF_Gal_crustaceans_taxa5))


##### Import TEP Crustaceans and get the taxa list ####

GBIF_TEP_crustaceans_list <- occ_download_get('0172684-240321170329656') %>%
   occ_download_import()

GBIF_TEP_crustaceans_taxa1<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(taxa=species)
GBIF_TEP_crustaceans_taxa2<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="GENUS")%>%
  select(genus)%>%
  rename(taxa=genus)
GBIF_TEP_crustaceans_taxa3<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="FAMILY")%>%
  select(family)%>%
  rename(taxa=family)
GBIF_TEP_crustaceans_taxa4<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="ORDER")%>%
  select(order)%>%
  rename(taxa=order)
GBIF_TEP_crustaceans_taxa5<-GBIF_TEP_crustaceans_list%>% 
  filter(taxonRank=="CLASS")%>%
  select(class)%>%
  rename(taxa=class)
GBIF_TEP_crustaceans_taxa<-unique(rbind(GBIF_TEP_crustaceans_taxa1,GBIF_TEP_crustaceans_taxa2,GBIF_TEP_crustaceans_taxa3,GBIF_TEP_crustaceans_taxa4, GBIF_TEP_crustaceans_taxa5))
  

##### Import Galapagos Fish and get the species list ####

GBIF_Gal_fish_list <- occ_download_get('0172899-240321170329656') %>%
     occ_download_import()

GBIF_Gal_fish_species<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(Species=species)

#get list including higher taxonomic levels to compare the results file with
GBIF_Gal_fish_taxa1<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(taxa=species)
GBIF_Gal_fish_taxa2<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="GENUS")%>%
  select(genus)%>%
  rename(taxa=genus)
GBIF_Gal_fish_taxa3<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="FAMILY")%>%
  select(family)%>%
  rename(taxa=family)
GBIF_Gal_fish_taxa4<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="ORDER")%>%
  select(order)%>%
  rename(taxa=order)
GBIF_Gal_fish_taxa5<-GBIF_Gal_fish_list%>% 
  filter(taxonRank=="CLASS")%>%
  select(class)%>%
  rename(taxa=class)
GBIF_Gal_fish_taxa<-unique(rbind(GBIF_Gal_fish_taxa1,GBIF_Gal_fish_taxa2,GBIF_Gal_fish_taxa3,GBIF_Gal_fish_taxa4, GBIF_Gal_fish_taxa5))


##### Import TEP Fish and get the taxa list ####
GBIF_TEP_fish_list<-occ_download_get('0172907-240321170329656') %>%
     occ_download_import()

GBIF_TEP_fish_taxa1<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="SPECIES")%>%
  select(species)%>%
  rename(taxa=species)
GBIF_TEP_fish_taxa2<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="GENUS")%>%
  select(genus)%>%
  rename(taxa=genus)
GBIF_TEP_fish_taxa3<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="FAMILY")%>%
  select(family)%>%
  rename(taxa=family)
GBIF_TEP_fish_taxa4<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="ORDER")%>%
  select(order)%>%
  rename(taxa=order)
GBIF_TEP_fish_taxa5<-GBIF_TEP_fish_list%>% 
  filter(taxonRank=="CLASS")%>%
  select(class)%>%
  rename(taxa=class)
GBIF_TEP_fish_taxa<-unique(rbind(GBIF_TEP_fish_taxa1,GBIF_TEP_fish_taxa2,GBIF_TEP_fish_taxa3,GBIF_TEP_fish_taxa4, GBIF_TEP_fish_taxa5))

#Darwin foundation ----

#### Galapagos Marine Fish List 2016 ####
curl::curl_download(url = "https://datazone.darwinfoundation.org/media/pdf/checklist/2016Aug24_Tirado-Sanchez_et_al_Galapagos_Pisces_Checklist.csv", destfile = "data/2016Aug24_Tirado-Sanchez_et_al_Galapagos_Pisces_Checklist.csv")
Darwin_fish<-read.csv("data/2016Aug24_Tirado-Sanchez_et_al_Galapagos_Pisces_Checklist.csv", fileEncoding = "latin1", check.names = FALSE)

#Combine Genus and species columns when species is not just sp. and filter everything else
Darwin_fish_species<- Darwin_fish%>%
  filter(`Specific Epithtet`!="sp.")%>%
  unite(Genus_species, c("Genus", "Specific Epithtet"), sep =" ")%>%
  select(Genus_species)%>%
  rename(Species=Genus_species)

#Get names at all taxonomic levels
Darwin_fish_taxa1<-Darwin_fish_species%>% 
  rename(taxa=Species)
Darwin_fish_taxa2<-Darwin_fish%>% 
  select(Genus)%>%
  rename(taxa=Genus)
Darwin_fish_taxa3<-Darwin_fish%>% 
  select(Family)%>%
  rename(taxa=Family)
Darwin_fish_taxa4<-Darwin_fish%>% 
  select(Order)%>%
  rename(taxa=Order)
Darwin_fish_taxa5<-Darwin_fish%>% 
  select(Class)%>%
  rename(taxa=Class)
Darwin_fish_taxa<-unique(rbind(Darwin_fish_taxa1,Darwin_fish_taxa2,Darwin_fish_taxa3,Darwin_fish_taxa4, Darwin_fish_taxa5))

####Galapagos Marine Crustaceans List 2016 ####

curl::curl_download(url = "https://datazone.darwinfoundation.org/media/pdf/checklist/2016Sep30_Tirado-Sanchez_et_al_Galapagos_Marine_crustaceans_Checklist.csv", destfile = "data/2016Sep30_Tirado-Sanchez_et_al_Galapagos_Marine_crustaceans_Checklist.csv")
Darwin_crustaceans<-read.csv("data/2016Sep30_Tirado-Sanchez_et_al_Galapagos_Marine_crustaceans_Checklist.csv", fileEncoding = "latin1", check.names = FALSE)

#Combine Genus and species columns when species is not just sp. and filter everything else
Darwin_crustacean_species<- Darwin_crustaceans%>%
  filter(`Specific Epithtet`!="sp.")%>%
  unite(Genus_species, c("Genus", "Specific Epithtet"), sep =" ")%>%
  select(Genus_species)%>%
  rename(Species=Genus_species)


#Get names at all taxonomic levels
Darwin_crustacean_taxa1<-Darwin_crustacean_species%>% 
  rename(taxa=Species)
Darwin_crustacean_taxa2<-Darwin_crustaceans%>% 
  select(Genus)%>%
  rename(taxa=Genus)
Darwin_crustacean_taxa3<-Darwin_crustaceans%>% 
  select(Family)%>%
  rename(taxa=Family)
Darwin_crustacean_taxa4<-Darwin_crustaceans%>% 
  select(Order)%>%
  rename(taxa=Order)
Darwin_crustacean_taxa5<-Darwin_crustaceans%>% 
  select(Class)%>%
  rename(taxa=Class)
Darwin_crustacean_taxa<-unique(rbind(Darwin_crustacean_taxa1,Darwin_crustacean_taxa2,Darwin_crustacean_taxa3,Darwin_crustacean_taxa4, Darwin_crustacean_taxa5))

#Combine OBIS, GBIF and Darwin lists----

Species_obis_gbif_darwin_galcrust<-unique(rbind(obis_galcrust_species,GBIF_Gal_crustaceans_species, Darwin_crustacean_species))
Species_obis_gbif_darwin_galcrust<-Species_obis_gbif_darwin_galcrust%>%
  filter(Species!="NA")

#order alphabetically
Species_obis_gbif_darwin_galcrust$Species<- Species_obis_gbif_darwin_galcrust$Species[order(Species_obis_gbif_darwin_galcrust$Species)]





Species_obis_gbif_darwin_galfish<-unique(rbind(obis_galfish_species,GBIF_Gal_fish_species, Darwin_fish_species))
Species_obis_gbif_darwin_galfish<-Species_obis_gbif_darwin_galfish%>%
  filter(Species!="NA")

#order alphabetically
Species_obis_gbif_darwin_galfish$Species<- Species_obis_gbif_darwin_galfish$Species[order(Species_obis_gbif_darwin_galfish$Species)]

##Write Fish list to a txt file to use as CRABS input----

write_delim(Species_obis_gbif_darwin_galfish, "../custom_db/comprehensive_galapagos_fish_list.txt", delim = '\t', col_names = FALSE)
nrow(Species_obis_gbif_darwin_galfish)
#1147


##Write Crustaceans list to a txt file to use as CRABS input----

write_delim(Species_obis_gbif_darwin_galcrust, "../custom_db/comprehensive_galapagos_crustaceans_list.txt", delim = '\t', col_names = FALSE)
nrow(Species_obis_gbif_darwin_galcrust)
#4498

