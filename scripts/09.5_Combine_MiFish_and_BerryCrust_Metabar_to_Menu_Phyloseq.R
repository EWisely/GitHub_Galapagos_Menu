#Eldridge Wisely
#Date: 4-30-24
#Export fasta for alignment, Create Phyloseq objects from individual Metabar objects

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
library(msa)

'%ni%' <- Negate("%in%")

#set variables
Primer1<-"MiFish"
Primer2<-"BerryCrust"


#Import Metabarlist intermediary files----

##Primer1 = MiFish----

df.MiFish_Menu_reads<-read.csv(paste("../07_lulu_metabar/",Primer1,"_Metabar_cleaned_reads.csv",sep=""))
df.MiFish_Menu_motus<-read.csv(paste("../07_lulu_metabar/",Primer1,"_Metabar_cleaned_motus.csv",sep=""))
df.MiFish_Menu_pcrs<-read.csv(paste("../07_lulu_metabar/",Primer1,"_Metabar_cleaned_pcrs.csv",sep=""))
df.MiFish_Menu_samples<-read.csv(paste("../07_lulu_metabar/",Primer1,"_Metabar_cleaned_samples.csv",sep=""))

##Primer2 = BerryCrust----

df.BerryCrust_Menu_reads<-read.csv(paste("../07_lulu_metabar/",Primer2,"_Metabar_cleaned_reads.csv",sep=""))
df.BerryCrust_Menu_motus<-read.csv(paste("../07_lulu_metabar/",Primer2,"_Metabar_cleaned_motus.csv",sep=""))
df.BerryCrust_Menu_pcrs<-read.csv(paste("../07_lulu_metabar/",Primer2,"_Metabar_cleaned_pcrs.csv",sep=""))
df.BerryCrust_Menu_samples<-read.csv(paste("../07_lulu_metabar/",Primer2,"_Metabar_cleaned_samples.csv",sep=""))

### Merge reads tables----
df.Menu_merged.reads<-full_join(df.MiFish_Menu_reads,df.BerryCrust_Menu_reads)
### Export to a new .csv
df.Menu_merged.reads<-column_to_rownames(df.Menu_merged.reads, var = "X")
write.csv(df.Menu_merged.reads, "../09_Metabar_to_Phyloseq/Menu_merged_reads.csv")


### Merge motus tables----
df.Menu_merged.motus<-full_join(df.MiFish_Menu_motus,df.BerryCrust_Menu_motus)
### Export to a new .csv
df.Menu_merged.motus<-column_to_rownames(df.Menu_merged.motus, var = "X")
write.csv(df.Menu_merged.motus, "../09_Metabar_to_Phyloseq/Menu_merged_motus.csv")

### Merge pcrs tables----

df.Menu_merged.pcrs<-full_join(df.MiFish_Menu_pcrs,df.BerryCrust_Menu_pcrs, by = c("X","sample_id","rep","type","control_type","Site_Name","Microhabitat"))
### Export to a new .csv
df.Menu_merged.pcrs<-column_to_rownames(df.Menu_merged.pcrs, var = "X")
write.csv(df.Menu_merged.pcrs, "../09_Metabar_to_Phyloseq/Menu_merged_pcrs.csv")


### Merge samples tables ----



df.Menu_merged.samples<-union_all(df.MiFish_Menu_samples,df.BerryCrust_Menu_samples)

df.Menu_merged.samples<-df.Menu_merged.samples %>% distinct(X, .keep_all = TRUE)

df.Menu_merged.samples<-column_to_rownames(df.Menu_merged.samples, var="X")

write.csv(df.Menu_merged.samples, "../09_Metabar_to_Phyloseq/Menu_merged_samples.csv")



#bring the full fish and crustacean dataset back into MetabaR.
Menu_Combined<-tabfiles_to_metabarlist(
  file_reads= "../09_Metabar_to_Phyloseq/Menu_merged_reads.csv",
  file_motus="../09_Metabar_to_Phyloseq/Menu_merged_motus.csv",
  file_pcrs ="../09_Metabar_to_Phyloseq/Menu_merged_pcrs.csv",
  file_samples ="../09_Metabar_to_Phyloseq/Menu_merged_samples.csv",
  files_sep = ","
)

### make a combined fasta file ----

fasta_generator(
  Menu_Combined,
  id = rownames(Menu_Combined$motus),
  output_file = paste0("../09_Metabar_to_Phyloseq/Menu_Combined.fasta"),
  annotation = NULL,
  annot_sep = ";"
)


#Put Merged Metabar object into Phyloseq format----


#Phyloseq!
#otu_table - Works on any numeric matrix. You must also specify if the species are rows or columns
#sample_data - Works on any data.frame. The rownames must match the sample names in the otu_table if you plan to combine them as a phyloseq-object
#tax_table - Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
#phyloseq - Takes as argument an otu_table and any unordered list of valid phyloseq components: sample_data, tax_table, phylo, or XStringSet. The tip labels of a phylo-object (tree) must match the OTU names of the otu_table, and similarly, the sequence names of an XStringSet object must match the OTU names of the otu_table.
#merge_phyloseq - Can take any number of phyloseq objects and/or phyloseq components, and attempts to combine them into one larger phyloseq object. This is most-useful for adding separately-imported components to an already-created phyloseq object.


##otu table----

clean.otu.df<-as.data.frame(t(Menu_Combined$reads))
clean.otu.df<-rownames_to_column(clean.otu.df, var="id")


##taxa table----

clean.taxa.df<-as.data.frame(Menu_Combined$motus)
clean.taxa.df<-rownames_to_column(clean.taxa.df, var="id")
taxa.df<-clean.taxa.df%>%
  select(id, superkingdom, phylum, class, order,family,genus,species)


##samples table----

clean.samples.df<-as.data.frame(Menu_Combined$samples)
clean.samples.df<-rownames_to_column(clean.samples.df, var ="sample_id")

clean.pcrs.df<-Menu_Combined$pcrs
clean.pcrs.df<-rownames_to_column(clean.pcrs.df, var="sample")


#Before I can do this I need to either aggregate the PCRs again... or combine by PCR name instead of sample name


samples.df<-full_join(clean.pcrs.df,clean.samples.df, by=c("sample_id","Site_Name","Microhabitat"))


samples.df <- samples.df %>% 
  tibble::column_to_rownames("sample")


##make otu matrix----

clean.otu.df<-column_to_rownames(clean.otu.df, var="id")
clean.otu.mat<-as.matrix(clean.otu.df)

#make taxa matrix
taxa.df<-column_to_rownames(taxa.df, var="id")
clean.taxa.mat<-as.matrix(taxa.df)


##ASV sequences ----
library(msa)
Seqs <- readDNAStringSet(paste0("../09_Metabar_to_Phyloseq/Menu_Combined.fasta"))


#Put them together into a phyloseq object
Menu_Combined.ps<- phyloseq(otu_table(clean.otu.mat, taxa_are_rows = TRUE), 
                         sample_data(samples.df), 
                         tax_table(clean.taxa.mat),
                          refseq(Seqs))



Menu_Combined.ps


saveRDS(Menu_Combined.ps, file=paste("../09_Metabar_to_Phyloseq/Menu_Combined_Phyloseq.rds",sep=""))

print(paste0("Finished creating a Phyloseq object out of Combined Metabarlist!"))

