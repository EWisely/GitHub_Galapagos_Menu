# GitHub_Galapagos_Menu

#Hi Kristina (and Eric and Renata!), 

Here are my slides I made for our meeting: https://docs.google.com/presentation/d/1SH_PuTqrIKmlbS94F_BS7VBjIYiq1QZ9V5mEtxXeRzE/edit?usp=sharing

I look forward to working together to make my script for validating eDNA taxonomic assignments a publishable product that people could use!  

The problem: incomplete (and uncurated) reference DNA databases used for taxonomic assignment of query sequences from eDNA (water, soil, air) samples

;)
#Global database (all vertebrates when assigning fish or all invertebrates when assigning crustacean sequences for example):

Global_Pros <- c("likely to at least match relative when the local isn't sequenced", "most completeness possible", "unbiased (least bias anyway)")

Global_Cons <- c("likely to get a matching relative when the local isn't sequenced", "that's a lot of invasive species!","local species added to the database are filtered out by a common database cleaning step")

#Local database (sequences downloaded from the global database (NCBI in this case) matching the local species checklist -from script 04):

Local_Pros <- c("no invasive species", "the problematic local database cleaning step can be left out of the CRABS workflow")

Local_Cons <- c("when the local species has not been sequenced the result is "NA", "different results for different species checklists", "local assignments may be overreaching the limitations of the method -low primer specificity could result in all members of a genus matching equally to the query sequence")

#By combining the results of both methods and comparing them methodically to choose the best assignment database for each sequence, we can have a cross-validated taxonomic assignment of the eDNA sequences which can take advantage of the pros of each method and mitigate the cons.  

I think that script 06 is the one that would be really useful for people doing metabarcoding studies. The ones after that are just me working through the Galapagos marine eDNA dataset.  
Script 04 is potentially also broadly useful by downloading species checklists for a given region (user-defined coordinates) from OBIS and GBIF (and I combined it with a species checklist csv downloaded from the Darwin center website) to create a list of local species for:
1) making a list of species to download as the local database (the actual downloading, cleaning, and formatting of the local database is done with a linux program called CRABS https://github.com/gjeunen/reference_database_creator). And
2) it's used in the last step of my script 6 where if the local and global database can't be reconciled (for example, when the local database had no match), to find a match in the local checklist (sometimes at a higher taxonomic level).




 
