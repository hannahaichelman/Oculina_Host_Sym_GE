# This script analyzes the ITS2 data from the Oculina arbuscula holobiont experiment. 
# Author: Hannah E Aichelman
# Last Updated: January 1, 2023


#### Load packages and read in data ####

#packages
#install.packages("decontam")
library(decontam)
packageVersion("decontam") #1.16.0 HA's version
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(viridis)
library(plyr)
#library(microbiomeutilities)

# SymPortal ITS2 DIV Analysis
# cleaned file up to remove extraneous info in the header in excel, but original file here:
# Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/ITS2_Analyses/20210826_aichelman/its2_type_profiles/166_20210826_DBV_20210826T064400.profiles.absolute.abund_and_meta.txt
setwd("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/09_Supplemental_Info/Symbiont_Identity/ITS2_Analyses")
its2_divs = read.csv("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/09_Supplemental_Info/Symbiont_Identity/ITS2_Analyses/20210826_aichelman/its2_type_profiles/SymPortal_RawDIVs.csv")
head(its2_divs)

# remove samples not included in the dataset and clones
its2_divs = its2_divs %>%
  filter(frag != "DW")

str(its2_divs)

# add some identifying info
treatment = as.factor(sapply(strsplit(its2_divs$frag, split = ""), "[[", 1))

genotype = as.factor(sapply(strsplit(its2_divs$frag, split = ""), "[[", 2))

coral_id = substring(its2_divs$frag, 2)

metadata = data.frame(its2_divs$frag, treatment, genotype, coral_id)

metadata = metadata %>%
  dplyr::rename(frag = its2_divs.frag)

# need to figure out symbiotic state, which we have in the PAM data
pam = read.csv('~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/06_Sym_Physiology/Oculina_PAM.csv')

pam2 = pam %>%
  select(coral_id, sym_state)

pam3 = pam2 %>%
  distinct()

# merge symbiotic state info with its2 data
samdf = left_join(metadata, pam3, by = "coral_id")
str(samdf)

samdf$sym_state = as.factor(samdf$sym_state)
levels(samdf$sym_state) = c("W","B")

samdf$full_id <- paste(samdf$coral_id,samdf$treatment,samdf$sym_state, sep = "_")

#rownames have to match between counts table & sample data table or else phyloseq will throw a fit
rownames(samdf) = samdf$full_id

#making a taxa table for phyloseq
its2_divs_ps = left_join(its2_divs, samdf, by = "frag")

its2_divs_ps = its2_divs_ps %>%
  select(full_id, B2.B2e, B2, B5a) %>%
  column_to_rownames("full_id")

taxa <- data.frame(colnames(its2_divs_ps)) #extract sym data

colnames(taxa) <- c("DIV") #changing the column name to be more user-friendly
taxa$genus = str_sub(taxa$DIV, 1, 1)
str(taxa)

taxa$DIV = as.factor(taxa$DIV)
taxa$genus = as.factor(taxa$genus)

#rownmaes also have to match between the columns of the counts table & the taxa table or ps freaks out
rownames(taxa) <- taxa$DIV

taxa.m <- as.matrix(taxa) #also has to be a matrix

ps.its2 <- phyloseq(sample_data(samdf),
                    otu_table(its2_divs_ps,taxa_are_rows=FALSE),
                    tax_table(taxa.m))
ps.its2
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 3 taxa and 39 samples ]
# sample_data() Sample Data:       [ 39 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 3 taxa by 2 taxonomic ranks ]

#### Bar plot - raw table pre-processing ####
# Skipping decontam for now since we are ignoring negative controls per ben hume
ps.rel <- transform_sample_counts(ps.its2, function(OTU) OTU/sum(OTU))

plot_bar(ps.rel, x="frag",fill="DIV")+
  theme_classic()

#keeps samples with summed counts greater than 0
ps.its2.no0 <- prune_samples(sample_sums(ps.its2)!=0, ps.its2)
ps.its2.no0 # don't lose any, can use ps.rel above

# save phyloseq object 
saveRDS(ps.its2.no0, "ps.its2.RDS")

#### Bar plot - post-processing ####

its2_cols_purples = c("#8c96c6", "#810f7c", "#edf8fb")

barplot = plot_bar(ps.rel, fill="DIV") +
  theme_bw() +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "DIV", values = its2_cols_purples) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.4, hjust= 0)) 
barplot

ggsave(barplot, filename = "its2_barplot.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)

