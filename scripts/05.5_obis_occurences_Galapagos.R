#Eldridge Wisely  
#4-12-24
#make a checklist of Galapagos species according to (Ocean Biodiversity Information System) OBIS database 
library(robis)
library(tidyverse)
#library(worrms)
#library(curl)
library(taxize)
library(taxonomizr)


#for wkt geometry draw it on this webapp https://wktmap.com

#for just the Galapagos: POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))
#for the TEP: POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))

#Galapagos fish checklist
obis_galfish<-checklist(c("Agnatha", "Chondrichthyes", "Osteichthyes"),geometry = "POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))")

#TEP fish checklist
obis_tepfish<-checklist(c("Agnatha", "Chondrichthyes", "Osteichthyes"),geometry = "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))")


#Galapagos crustacean checklist
obis_galcrust<-checklist("Crustacea",geometry = "POLYGON ((-93.339844 -3.162456, -93.339844 2.547988, -87.1875 2.547988, -87.1875 -3.162456, -93.339844 -3.162456))")

#TEP crustacean checklist
obis_tepcrust<-checklist("Crustacea",geometry = "POLYGON ((-117.421875 31.952162, -91.933594 -6.315299, -81.386719 -6.315299, -76.113281 7.710992, -82.089844 8.581021, -87.011719 13.581921, -104.238281 20.303418, -112.5 32.249974, -117.421875 31.952162))")

#next check the results of obitools assignment with this list, and if obitools assigned taxon is in this list, put TRUE in present_local(GAL) or present_regional(TEP) column of obitools results to end up in the MOTUs table.  At first it should be formatted ASV, obi_name, obi_TAXID, obi_ident, present_local, present_regional

obi_fish95<-readr::read_delim("../03_obitools_results/MiFish_Menu_95_named.tab")

colnames(obi_fish95)
obi_fish95<-obi_fish95%>%
  select(c(ID,TAXID,SCIENTIFIC_NAME,BEST_IDENTITY))%>%
  filter(BEST_IDENTITY>=0.95)

obi_fish_wormsIDs<-taxize::get_wormsid(obi_fish95$SCIENTIFIC_NAME)
#• Total: 542 
#• Found: 532 
#• Not Found: 10



obi_fish95 <- obi_fish95 %>%  mutate(obi_fish_wormsIDs)

colnames(obis_galfish)

local_obi_fish<-cbind(obi_fish95,sapply(obi_fish95$SCIENTIFIC_NAME, `%in%`, obis_tepfish$scientificName))
local_obi_fish<-rename(local_obi_fish,Present_regional="sapply(obi_fish95$SCIENTIFIC_NAME, `%in%`, obis_tepfish$scientificName)")

percent_regionally_present<-sum(local_obi_fish["Present_regional"])/nrow(local_obi_fish)
percent_regionally_present
#0.8173432 when using scientific names
#0.8154982 when using wormsIDs
#0.8118081 when using NCBI ids

local_obi_fish<-cbind(local_obi_fish,sapply(obi_fish95$SCIENTIFIC_NAME, `%in%`, obis_galfish$scientificName))
local_obi_fish<-rename(local_obi_fish,Present_local="sapply(obi_fish95$SCIENTIFIC_NAME, \`%in%\`, obis_galfish$scientificName)")

percent_locally_present<-sum(local_obi_fish["Present_local"])/nrow(local_obi_fish)
percent_locally_present
#0.7121771

#local usearch database results
#taxonomizr::prepareDatabase()

#load table and split columns appropriately
usearch_fish95<-readr::read_delim("../04_usearch_obitools_results/Menu_Galapagos_fish_usearch_results.txt", col_names = FALSE)

ufish_95<-separate(data = usearch_fish95, col = X1, sep = " ",into = c("ID","COUNT","UMI"))
ufish_95<-separate(data = ufish_95, col = X2, sep = " ",into = c("NCBI_ACC","GENUS","SPECIES","DEFINITION"))
ufish_95<-rename(ufish_95, usearch_pctid=X3)

ufish_95$Scientific_Name<-paste(ufish_95$GENUS, ufish_95$SPECIES, sep=" ")

colnames(ufish_95)
ufish_95<-ufish_95%>%
  select(c(ID,NCBI_ACC,Scientific_Name,usearch_pctid))%>%
  filter(usearch_pctid>=0.95)


#ufish_taxIDs<- taxize::genbank2uid(list(ufish_95$NCBI_ACC), batch_size=100)

#taxonomizr::

colnames(obis_galfish)

local_usearch_fish<-cbind(ufish_95,sapply(ufish_95$Scientific_Name, `%in%`, obis_tepfish$scientificName))
local_usearch_fish<-rename(local_usearch_fish,Present_regional="sapply(ufish_95$Scientific_Name, `%in%`, obis_tepfish$scientificName)")

percent_regionally_present<-sum(local_usearch_fish[,5])/nrow(local_usearch_fish)
percent_regionally_present
#0.997093

local_usearch_fish<-cbind(local_usearch_fish,sapply(ufish_95$Scientific_Name, `%in%`, obis_galfish$scientificName))
local_usearch_fish<-rename(local_usearch_fish,Present_local="sapply(ufish_95$Scientific_Name, `%in%`, obis_galfish$scientificName)")

percent_locally_present<-sum(local_usearch_fish[,6])/nrow(local_usearch_fish)
percent_locally_present
#0.9927326






       
       
       
       #select mo.*, ma.*,  -- better enumerate the columns you want here
#exists (select 1 from modeles_scrape ms where ms.marque = ma.display_name ) as isactive
#from modeles  mo
#inner join marques ma using(marque_id)
