#Eldridge Wisely
#Date: 4-30-24
#Export fasta for alignment, Create Phyloseq objects from individual Metabar objects.

#library("ggplot2")      # graphics
#library("readxl")       # necessary to import the data from Excel file
#library("tibble")       # Needed for converting column to row names
# Load requested package for plotting
#library("reshape2")
library("tidyverse")
library("dplyr")        # filter and reformat data frames
library(metabaR)
library(phyloseq)
library(adegenet)
library(ape)
#library(phangorn)

'%ni%' <- Negate("%in%")

#set variables
Primer<-"BerryCrust"


#Import Metabarlist----
Menu<-readRDS(file= paste("../07_lulu_metabar/",Primer,"_Menu_Clean_Metabarlist.rds", sep = ""))

# Compute the number of reads per pcr
Menu$pcrs$nb_reads <- rowSums(Menu$reads)

# Compute the number of motus per pcr
Menu$pcrs$nb_motus <- rowSums(Menu$reads>0)

#subset the metabarlist to all PCRs with more than 0 reads
Menu <- subset_metabarlist(Menu, table = "pcrs",
                           indices = Menu$pcrs$nb_reads>0)
summary_metabarlist(Menu)

#Export Fasta for alignment to make phylogeny----

 fasta_generator(
   Menu,
   id = rownames(Menu$motus),
   output_file = paste0("../09_Metabar_to_Phyloseq/",Primer,"_Cleaned.fasta"),
   annotation = NULL,
   annot_sep = ";"
 )

#Read fasta into Geneious Prime and align and export tree file----


#Put Metabar object into Phyloseq format----


#Phyloseq!
#otu_table - Works on any numeric matrix. You must also specify if the species are rows or columns
#sample_data - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object
#tax_table - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
#phyloseq - Takes as argument an otu_table and any unordered list of valid phyloseq components: sample_data, tax_table, phylo, or XStringSet. The tip labels of a phylo-object (tree) must match the OTU names of the otu_table, and similarly, the sequence names of an XStringSet object must match the OTU names of the otu_table.
#merge_phyloseq - Can take any number of phyloseq objects and/or phyloseq components, and attempts to combine them into one larger phyloseq object. This is most-useful for adding separately-imported components to an already-created phyloseq object.


##otu table----

clean.otu.df<-as.data.frame(t(Menu$reads))
clean.otu.df<-rownames_to_column(clean.otu.df, var="id")


##taxa table----

clean.taxa.df<-as.data.frame(Menu$motus)
clean.taxa.df<-rownames_to_column(clean.taxa.df, var="id")
taxa.df<-clean.taxa.df%>%
  select(id, superkingdom, phylum, class, order,family,genus,species)


##samples table----

clean.samples.df<-as.data.frame(Menu$samples)
clean.samples.df<-rownames_to_column(clean.samples.df, var ="sample_id")

clean.pcrs.df<-Menu$pcrs

samples.df<-full_join(clean.pcrs.df,clean.samples.df, by="sample_id")


samples.df <- samples.df %>% 
  tibble::column_to_rownames("sample_id")%>%
  select(-c(Site_Name.x,Microhabitat.x))%>%
  dplyr::rename(Site_Name ="Site_Name.y")%>%
  dplyr::rename(Microhabitat ="Microhabitat.y")


##make otu matrix----

clean.otu.df<-column_to_rownames(clean.otu.df, var="id")
clean.otu.mat<-as.matrix(clean.otu.df)

#make taxa matrix
taxa.df<-column_to_rownames(taxa.df, var="id")
clean.taxa.mat<-as.matrix(taxa.df)

##phy_tree ----

newick<-read_tree(paste0("../09_Metabar_to_Phyloseq/",Primer,"_Cleaned Alignment Clustal Omega FastTree Tree.newick"))


##ASV sequences ----
library(msa)
Seqs <- readDNAStringSet(paste0("../09_Metabar_to_Phyloseq/",Primer,"_Cleaned.fasta"))


#Put them together into a phyloseq object
Menu.ps<- phyloseq(otu_table(clean.otu.mat, taxa_are_rows = TRUE), 
                         sample_data(samples.df), 
                         tax_table(clean.taxa.mat),
                          phy_tree(newick),
                          refseq(Seqs))



Menu.ps


saveRDS(Menu.ps, file=paste("../09_Metabar_to_Phyloseq/",Primer,"_Phyloseq.rds",sep=""))

print(paste0("Finished creating a Phyloseq object out of ",Primer," Metabarlist!"))

