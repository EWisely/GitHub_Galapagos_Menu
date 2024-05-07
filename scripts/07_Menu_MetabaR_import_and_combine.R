#Eldridge Wisely
#Date: 4-29-24
#LULU and make Metabarlists for Menu project MiFish and BerryCrust primers 

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

'%ni%' <- Negate("%in%")

#set variables
Primer<-"MiFish"


#Read in named tab file of ASV counts from obitools and cleaned of all extraneous columns
otutab <- read.csv(paste("../03_obitools_results/",Primer,"_Menu_named_tab_LULU.txt",sep = ""), sep='\t', header=TRUE, as.is=TRUE, row.names = "ID")

matchlist<- read.table(paste("../03_obitools_results/",Primer,"_Menu_95_named_matchlist.txt",sep = ""), header=FALSE, as.is=TRUE, stringsAsFactors=FALSE)

#run LULU ----
#curated_result <- lulu(otutab, matchlist)

#curated_result$curated_table # Curated OTU table
#curated_result$curated_count # Number of OTUs retained
#curated_result$curated_otus # IDs of curated OTUs

#curated_result$discarded_count # OTUs discarded
#curated_result$otu_map # total - total read count, spread - the number of samples the OTU is present in
# parent_id - ID of OTU with which this OTU was merged (or self)
# curated - ("parent" or "merged"), was this OTU kept as a valid OTU (parent) or merged with another
# rank - The rank of the OTU in terms of decreasing spread and read count

#curated_result$original_table # Original OTU table
#write.csv(curated_result$curated_table,paste("../07_lulu_metabar/lulu_curated_",Primer,"_Menu_otutable.tab",sep = ""))

lulubackup<-read.csv(paste("../07_lulu_metabar/lulu_curated_",Primer,"_Menu_otutable.tab",sep = ""))

#Metabarlist requirements---- 
#Arguments
#reads	
#MOTU abundance table. Rows and rownames of the table should correspond to PCRs and their names respectively. Columns and colnames should correspond to MOTUs and their names. Rownames in this table should correspond to PCR names respectively.

#motus	
#MOTU characteristics table (e.g. taxonomy, sequence, etc.). Rows and rownames of the table should correspond to MOTUs and their names respectively, and the columns to their characteristics. Mandatory fields: 'sequence', i.e. the sequence representative of the MOTU.

#pcrs	
#PCR characteristics table (e.g. tags, primers, plate wells, etc.). Rows and rownames of the table should correspond to PCRs and their names respectively, and the columns to their characteristics. Mandatory fields: (i) 'sample_id', i.e. the name of each biological sample. (ii) 'type', i.e. the type of PCR; can be 'sample' or 'control'. (iii) 'control_type', i.e. the type of control if applicable. Should be either: 'NA' for samples, 'extraction' for extraction negative controls, 'pcr' for PCR negative controls, 'sequencing' for sequencing negative controls (e.g. unused tag combinations), or 'positive' for positive controls.

#samples	
#Sample characteristics table. Rows and rownames of the table should correspond to biological samples and their names respectively, and the columns to their environnemental characteristics.

#Format for MetabaR----

##reads----
#if doing a fresh LULU run:
#lulu_table<-otu_table(curated_result$curated_table, taxa_are_rows = TRUE) 


lulubackup<-column_to_rownames(lulubackup, var = "X")
lulu_table<-otu_table(lulubackup, taxa_are_rows = TRUE) #lulubackup


lulu_table.df<-as.data.frame(lulu_table)
lulu_table.df<-rownames_to_column(lulu_table.df, var ="id")

lulu_table.df<-as.data.frame(lulu_table)

Metabar_formatted_reads<-t(lulu_table.df)

Metabar_formatted_reads<-as.data.frame(Metabar_formatted_reads)

#reads all looked good and just needed . replaced with - in the sample names

Metabar_formatted_reads<-rownames_to_column(Metabar_formatted_reads, var = "sample")
Metabar_formatted_reads$sample
Metabar_formatted_reads$sample<-gsub(".","-",Metabar_formatted_reads$sample, fixed = TRUE)
Metabar_formatted_reads$sample
Metabar_formatted_reads<-column_to_rownames(Metabar_formatted_reads, var = "sample")


Metabar_formatted_reads<-as.matrix(Metabar_formatted_reads)



#reads (abundance table)
metabar_reads_table<-read.csv(paste("../07_lulu_metabar/lulu_curated_",Primer,"_Menu_otutable.tab",sep = ""), header = TRUE, strip.white = TRUE)
metabar_reads_table<-column_to_rownames(metabar_reads_table, var="X")

metabar_reads_table<-t(metabar_reads_table)
reads_colnames<-colnames(metabar_reads_table)
reads_rownames<-rownames(metabar_reads_table)

#fix the accidental replacement of - with . 
reads_colnames <- sapply(reads_colnames, gsub, pattern = ".", replacement = "-", fixed = TRUE)
reads_rownames <- sapply(reads_rownames, gsub, pattern = ".", replacement = "-", fixed =TRUE)

#put it in a numerical matrix then put the rownames and colnames back on
Metabar_formatted_reads<-base::as.matrix(metabar_reads_table, 
                                         ncol = ncol(metabar_reads_table))

Metabar_formatted_reads<-as.numeric(Metabar_formatted_reads)
Metabar_formatted_reads<-base::matrix(Metabar_formatted_reads, 
                                      ncol = ncol(metabar_reads_table))
colnames(Metabar_formatted_reads)<-reads_colnames
rownames(Metabar_formatted_reads)<-reads_rownames


##motus----
taxonomy_motus<-read.csv(paste("../06_local_vs_global_results/",Primer,"_Menu_ready_for_MetabaR.csv",sep = ""), row.names = "ASV")

lulu_table.df<-rownames_to_column(lulu_table.df, var ="id")
taxonomy_motus<-rownames_to_column(taxonomy_motus, var="id")
Metabar_formatted_motus<- inner_join(lulu_table.df,taxonomy_motus,by="id")
Metabar_formatted_motus<-column_to_rownames(Metabar_formatted_motus, var="id")
Metabar_formatted_motus<-Metabar_formatted_motus%>%
  select(superkingdom,phylum,class,order,family,genus,species,taxID,pct_id,Scientific_name,database,COUNT,sequence)



##pcrs----

pcrs <- data.frame(
  sample_id = sapply(strsplit(rownames(Metabar_formatted_reads), "-"), "[[", 1),
  rep = sapply(strsplit(rownames(Metabar_formatted_reads), "-"), "[[", 2),
  type = "sample",
  control_type = "NA",
  row.names = rownames(Metabar_formatted_reads)
)
pcrs_to_edit<-rownames_to_column(pcrs, var = "id")

write_csv(pcrs_to_edit, paste("../07_lulu_metabar/",Primer,"_pcrs_to_edit.csv",sep=""))

#edit in excel to make sure that type and control type are set correctly! save as tab-delimited text file with the following name: {Primer}_Metabar_formatted_PCRs_edited.txt

pcrs<- read_tsv(paste("../07_lulu_metabar/",Primer,"_Metabar_formatted_PCRs_edited.txt", sep=""), col_names = TRUE)
Metabar_formatted_pcrs<-column_to_rownames(pcrs, var="id")
Metabar_formatted_pcrs<- as.data.frame(Metabar_formatted_pcrs)

##samples----
sample_metadata<-read.delim("../07_lulu_metabar/Menu_sample_metadata.txt")

sample_df<-data.frame(Sample_ID=pcrs$sample_id)

Metabar_formatted_samples<-semi_join(sample_metadata,sample_df)

Metabar_formatted_samples<-as.data.frame(Metabar_formatted_samples)


Metabar_formatted_samples<-column_to_rownames(Metabar_formatted_samples, var= "Sample_ID")

#make a metabarlist!----



GAL_Menu<-metabarlist_generator(reads=Metabar_formatted_reads, Metabar_formatted_motus, pcrs=Metabar_formatted_pcrs, samples=Metabar_formatted_samples)
summary_metabarlist(GAL_Menu)


#save intermediate files for each primer----

write.csv(Metabar_formatted_reads, paste("../07_lulu_metabar/",Primer,"_Metabar_formatted_reads.csv",sep=""))
write.csv(Metabar_formatted_motus, paste("../07_lulu_metabar/",Primer,"_Metabar_formatted_motus.csv",sep=""))
write.csv(Metabar_formatted_pcrs, paste("../07_lulu_metabar/",Primer,"_Metabar_formatted_pcrs.csv",sep=""))
write.csv(Metabar_formatted_samples, paste("../07_lulu_metabar/",Primer,"_Metabar_formatted_samples.csv",sep=""))



saveRDS(GAL_Menu, file=paste("../07_lulu_metabar/",Primer,"_Menu_Metabarlist.rds",sep=""))
print(paste0("Script 07 for ",Primer," Done!"))
