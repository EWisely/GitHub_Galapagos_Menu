#Eldridge Wisely
#April 20, 2024
#input: results list from small CRABS download (with keep originals --no) created by du -s *>filename.txt then sed to remove file endings resulting in filesize and Species name in files for each search term.
#output: lists of species for large CRABS download (with keep originals --yes) to get the CRABS_ncbi_download.fasta for each command with all species in the list and their full title in the FASTA file.

library(readr)
library(tidyverse)


#Read in Galapagos Fish CRABS results----

#download matches by name and mitochondrion or mitochondrial with CRABS
#read in CRABS download results
CRABS_ncbi_mito<-read_delim(file = "../custom_db/combo_db_ncbi_mito.txt", col_names = c("bytes","Species"))
nrow(CRABS_ncbi_mito)
#1147

#number that had mito sequences in NCBI
ncbi_galfish_mito<-CRABS_ncbi_mito%>%
  filter(bytes>0)

nrow(ncbi_galfish_mito)

#874

#order alphabetically
ncbi_galfish_mito$Species<- ncbi_galfish_mito$Species[order(ncbi_galfish_mito$Species)]

ncbi_galfish_mito_list<-ncbi_galfish_mito%>%
  select(Species)

#Write Fish list of only ones with NCBI results to a txt file to use as CRABS input with keep original flag =yes

write_delim(ncbi_galfish_mito_list, "../custom_db/comprehensive_galapagos_fish_list_ncbi.txt", delim = '\t', col_names = FALSE)

# continue summary statistics

#percent of names in the list with mito results
(nrow(ncbi_galfish_mito)/nrow(CRABS_ncbi_mito))*100
#76.19878 %

#get list of names not present in NCBI database
no_mito_match<-CRABS_ncbi_mito%>%
  filter(bytes==0)

nrow(no_mito_match)
#273


#Read in Galapagos Crustaceans CRABS results----

#download matches by name and mitochondrion or mitochondrial with CRABS
##read in CRABS download results for "mitochondrion"----
crust_CRABS_ncbi_mito<-read_delim(file = "../custom_db/crust_combo_db_ncbi_mito.txt", col_names = c("bytes","Species"))
nrow(crust_CRABS_ncbi_mito)
#4489

#number that had mito sequences in NCBI
ncbi_galcrust_mito<-crust_CRABS_ncbi_mito%>%
  filter(bytes>0)

nrow(ncbi_galcrust_mito)

#302


#order alphabetically
ncbi_galcrust_mito$Species<- ncbi_galcrust_mito$Species[order(ncbi_galcrust_mito$Species)]

ncbi_galcrust_mito_list<-ncbi_galcrust_mito%>%
  select(Species)

#Write Crustaceans list of only ones with NCBI results to a txt file to use as CRABS input with keep original flag =yes

write_delim(ncbi_galcrust_mito_list, "../custom_db/mito_galapagos_crustacean_list_ncbi.txt", delim = '\t', col_names = FALSE)

# continue summary statistics

#percent of names in the list with mito results
(nrow(ncbi_galcrust_mito)/nrow(crust_CRABS_ncbi_mito))*100
#6.727 %

#get list of names not present in NCBI database
crust_no_mito_match<-crust_CRABS_ncbi_mito%>%
  filter(bytes==0)

nrow(crust_no_mito_match)
#4187

##read in CRABS download results for "large subunit"----
crust_CRABS_ncbi_lsu<-read_delim(file = "../custom_db/crust_combo_db_ncbi_lsu.txt", col_names = c("bytes","Species"))
nrow(crust_CRABS_ncbi_lsu)
#4498

#number that had mito sequences in NCBI
ncbi_galcrust_lsu<-crust_CRABS_ncbi_lsu%>%
  filter(bytes>0)

nrow(ncbi_galcrust_lsu)

#806

#order alphabetically
ncbi_galcrust_lsu$Species<- ncbi_galcrust_lsu$Species[order(ncbi_galcrust_lsu$Species)]

ncbi_galcrust_lsu_list<-ncbi_galcrust_lsu%>%
  select(Species)

#Write Crustaceans list of only ones with NCBI results to a txt file to use as CRABS input with keep original flag =yes

write_delim(ncbi_galcrust_lsu_list, "../custom_db/large_subunit_galapagos_crustacean_list_ncbi.txt", delim = '\t', col_names = FALSE)

# continue summary statistics

#percent of names in the list with lsu results
(nrow(ncbi_galcrust_lsu)/nrow(crust_CRABS_ncbi_lsu))*100
#17.91908 %

#get list of names not present in NCBI database
crust_no_lsu_match<-crust_CRABS_ncbi_lsu%>%
  filter(bytes==0)

nrow(crust_no_lsu_match)
#3692

##read in CRABS download results for "16S"----
crust_CRABS_ncbi_16S<-read_delim(file = "../custom_db/crust_combo_db_ncbi_16S.txt", col_names = c("bytes","Species"))
nrow(crust_CRABS_ncbi_16S)
#4496

#number that had 16S sequences in NCBI
ncbi_galcrust_16S<-crust_CRABS_ncbi_16S%>%
  filter(bytes>0)

nrow(ncbi_galcrust_16S)

#1441

#order alphabetically
ncbi_galcrust_16S$Species<- ncbi_galcrust_16S$Species[order(ncbi_galcrust_16S$Species)]

ncbi_galcrust_16S_list<-ncbi_galcrust_16S%>%
  select(Species)

#Write Crustaceans list of only ones with NCBI results to a txt file to use as CRABS input with keep original flag =yes

write_delim(ncbi_galcrust_16S_list, "../custom_db/16S_galapagos_crustacean_list_ncbi.txt", delim = '\t', col_names = FALSE)

# continue summary statistics

#percent of names in the list with 16S results
(nrow(ncbi_galcrust_16S)/nrow(crust_CRABS_ncbi_16S))*100
#32.05071 %

#get list of names not present in NCBI database
crust_no_16S_match<-crust_CRABS_ncbi_16S%>%
  filter(bytes==0)

nrow(crust_no_16S_match)
#3055


# total crustacean list statistics


