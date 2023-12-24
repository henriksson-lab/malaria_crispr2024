 

# Two data set:
#   
#   gRNA seq: Tags (2023)
# Try to understand what it means and if we should include
# 
# gRNA seq: KO
# Priority
# Main Fig: Line graph all genes with both guides starting from day 5 (exclude day 4)
# Based on above calculate growth rate using global normalization, both gRNA separate
# Supplemental FigLine graph all genes with both guides starting from day 4



### groups: dispensable, slow growing, essential; 3x4 genes x 2 grna
## must be from earlier? 
#/corgi/otherdataset/ellenbushell/crispr_pools/2023march_pools/counts.csv                          

# #dispensable
# 0515000
# 0933700
# 1013300 
# 1037800
# 
# #slow
# 1034400 
# 1101400
# 1104200
# 1401600
# 
# #essential
# 0211000
# 0706400
# 1039700
# 1214100


# /corgi/otherdataset/ellenbushell/crispr_pools/2023march_pools/gene_annot.csv

geneannot <- read.csv("/corgi/otherdataset/ellenbushell/crispr_pools/2023march_pools/gene_annot.csv")



################################################################################
### "Tags" analysis, /corgi/otherdataset/ellenbushell/crispr_pools/EB_barseq_slowpool_CRISPR
################################################################################

#these are normal gene + grna ID

#just color each gene, then different symbols for each grna

#need metadata for all;
#which tag is what?

#"PBANKA_1224200_G1"

cnt <- read.csv("/corgi/otherdataset/ellenbushell/crispr_pools/EB_barseq_slowpool_CRISPR/counts.csv")

#rownames(cnt) %in% 
cnt




