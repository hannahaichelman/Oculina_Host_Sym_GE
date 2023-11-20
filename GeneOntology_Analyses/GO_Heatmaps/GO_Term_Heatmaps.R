## This script plots heat maps of the genes contained within interesting GO terms for the Breviolum psygmophilum symbionts, 
## both in host and in culture

#### Libraries and files ####
library(tidyverse)
library(plotly)
library(viridis)
library(scales)
library(plyr)

setwd("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/")


#### Heat Maps of Interesting GOs - Syms in Host ####
## Symbionts in Host ##
# read in go.obo database to figure out the go category id's that we want to pull.

go.obo = read.delim("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/go.obo")

# these terms are visually pulled from the Cold GO tree outputs, all terms related to photosynthesis 
inds = c(which(go.obo$format.version..1.2 == "name: tetrapyrrole binding"), # MF
         which(go.obo$format.version..1.2 == "name: chlorophyll binding"), # MF
         which(go.obo$format.version..1.2 == "name: photosystem"), # CC
         which(go.obo$format.version..1.2 == "name: thylakoid membrane"), # CC
         which(go.obo$format.version..1.2 == "name: photosynthesis, light harvesting"), # BP
         which(go.obo$format.version..1.2 == "name: protein-chromophore linkage") # BP
) 

want.go.obo = go.obo %>% 
  dplyr::slice(sort(unique(c(inds - 1)))) %>%
  mutate_at("format.version..1.2", str_replace, "id: ", "")

##pulling all of the photosynthesis terms and making a data frame
rldpval = read.csv("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/tables/SymInHost_RLDandPVALS.csv", header = TRUE) %>%
  dplyr::rename("gene" = "X")

# read in iso2go file for the algae
iso2go_sym = read.delim("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/B_psygmophilum_isogroup_to_GOterm.tab", sep = "\t", header = FALSE) %>%
  dplyr::rename("gene" = "V1") %>%
  dplyr::rename("GO_ID" = "V2")
 
#pulling all of the photosynthesis terms from the iso2go 
list = as.list(want.go.obo)
print(list)

#pulling all of the photosynthesis terms from the iso2go 
photo_genes = iso2go_sym %>%
  filter(grepl("GO:0009521|GO:0009765|GO:0016168|GO:0018298|GO:0042651|GO:0046906", GO_ID)) %>%
  left_join(rldpval)
head(photo_genes)
str(photo_genes)

# add gene name to this df - trim down names and add to photo_genes
iso2gene = read.delim("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_TagSeq/References/B_psygmophilum_transcriptome/B_psygmophilum_isogroup_to_genename.tab", sep = "\t", header = FALSE) %>%
  dplyr::rename("gene" = "V1") %>%
  dplyr::rename("gene_name" = "V2")
head(iso2gene)

# shorten gene name for easier reading on heatmap
iso2gene$gene_name = gsub(" OS.*","",as.character(iso2gene$gene_name))
iso2gene$gene_name = gsub(", .*","",as.character(iso2gene$gene_name))

head(iso2gene)
View(iso2gene)

photo_genes_anno = photo_genes %>%
  left_join(iso2gene)
head(photo_genes_anno)

# set raw p-value for GO enriched
p.val = 0.10 

# filter based on p-values from deseq results
conds=photo_genes_anno[photo_genes_anno$pval.cold.syminhost<=p.val & !is.na(photo_genes_anno$pval.cold.syminhost),]
length(conds[,1])
#6
head(conds)
rownames(conds)<-NULL

# make gene names row names
conds2 = conds %>%   
  column_to_rownames(var = "gene_name")
head(conds2)

# make unique row names - only need for p=0.2
#rownames(conds) <- make.unique(conds$gene_name)

#row.names(conds)=conds$gene_name
exp = conds2[, c(3:23)]
head(exp)


means=apply(exp,1,mean) # calculate means of rows
explc=exp-means # subtracting them
head(explc)

# order column names of explc so we can group columns by treatment
explc = explc %>%
  select("OC4_F_B","OD5_F_B","OF8_F_B","OI2_F_B","OJ14_F_B","OM2_F_B","OR8_F_B","OC9_C_B","OD4_C_B","OF7_C_B","OI1_C_B","OJ13_C_B","OM1_C_B", 
         "OR7_C_B","OC5_H_B","OD6_H_B","OF9_H_B","OI3_H_B","OJ15_H_B","OM3_H_B","OR9_H_B")
# this explc object is what we can use to make our heatmap

# now make heatmap
# set color palette
col0=colorRampPalette(rev(c("#c51b7d","#e9a3c9","#ffffff", "#5ab4ac","#01665e")))(100)

# make treatment data frame
treatment = as.factor(sapply(strsplit(colnames(explc), split = "_"), "[[", 2)) %>%
  revalue(c("C" = "control", "F" = "cold", "H" = "heat"))
expDesign = data.frame(colnames(explc), treatment)
expDesign = expDesign %>%
  column_to_rownames(var = "colnames.explc.")

expDesign$treatment = ordered(expDesign$treatment, levels = c("cold", "control", "heat"))
expDesign = expDesign[order(expDesign$treatment), , drop = FALSE]

# check that our input dataframes match before making heatmap
colnames(explc) == row.names(expDesign)

my_colour = list(treatment = c(cold = "#74c476", heat = "#fd8d3c", control = "#a6611a"))

#heat map of all photosynthesis genes and save
library(pheatmap)

photo.heatmap = pheatmap(explc, cluster_cols = F, scale = "row", color = col0, annotation_col = expDesign, annotation_colors = my_colour, show_rownames = TRUE, show_colnames = FALSE, border_color = "NA")
ggsave(photo.heatmap, file = "/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/plots/syminhost_heatmap_photosynthesis_p=.1.pdf", width=8, height=2.5, units=c("in"), useDingbats=FALSE)

#### Heat Maps of Interesting Photosynthesis GOs - Syms in Culture ####
## Symbionts in Culture ##
# Photosynthesis genes first
# read in go.obo database to figure out the go category id's that we want to pull.

go.obo = read.delim("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/go.obo")

# these terms are visually pulled from the GO tree outputs, all terms related to photosynthesis 
inds = c(which(go.obo$format.version..1.2 == "name: photosynthesis, light harvesting"), # BP
         which(go.obo$format.version..1.2 == "name: photosynthesis"), # BP
         which(go.obo$format.version..1.2 == "name: protein-chromophore linkage"), # BP
         which(go.obo$format.version..1.2 == "name: chloroplast-nucleus signaling pathway"), # BP
         which(go.obo$format.version..1.2 == "name: photosystem"), # CC
         which(go.obo$format.version..1.2 == "name: thylakoid membrane"), # CC
         which(go.obo$format.version..1.2 == "name: tetrapyrrole binding"), # MF
         which(go.obo$format.version..1.2 == "name: chlorophyll binding") # MF
) 

want.go.obo = go.obo %>% 
  dplyr::slice(sort(unique(c(inds - 1)))) %>%
  mutate_at("format.version..1.2", str_replace, "id: ", "")

##pulling all of the photosynthesis terms and making a data frame
rldpval = read.csv("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/Sym_TagSeq/Culture_RLDandPVALS.csv", header = TRUE) %>%
  dplyr::rename("gene" = "X")

# remove isogroup from rldpval gene
rldpval$gene = gsub("isogroup", "", rldpval$gene)

# read in iso2go file for the algae
iso2go_sym = read.delim("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/B_psygmophilum_isogroup_to_GOterm.tab", sep = "\t", header = FALSE) %>%
  dplyr::rename("gene" = "V1") %>%
  dplyr::rename("GO_ID" = "V2")

#pulling all of the photosynthesis terms from the iso2go 
list = as.list(want.go.obo)
print(list)

photo_genes = iso2go_sym %>%
  filter(grepl("GO:0009521|GO:0009765|GO:0010019|GO:0015979|GO:0016168|GO:0018298|GO:0042651|GO:0046906", GO_ID)) %>%
  left_join(rldpval)
head(photo_genes)
str(photo_genes)

# add gene name to this df - trim down names and add to photo_genes
iso2gene = read.delim("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_TagSeq/References/B_psygmophilum_transcriptome/B_psygmophilum_isogroup_to_genename.tab", sep = "\t", header = FALSE) %>%
  dplyr::rename("gene" = "V1") %>%
  dplyr::rename("gene_name" = "V2")
head(iso2gene)

iso2gene$gene_name = gsub("OS=.*", "", iso2gene$gene_name)
iso2gene$gene_name = gsub(", .*","",as.character(iso2gene$gene_name))

head(iso2gene)
str(iso2gene)

photo_genes_anno = photo_genes %>%
  left_join(iso2gene)
head(photo_genes_anno)

# set raw p-value for GO enriched
p.val = 0.01

# filter based on p-values from deseq results
conds=photo_genes_anno[photo_genes_anno$pval.cold.culture<=p.val & !is.na(photo_genes_anno$pval.cold.culture),]
length(conds[,1])
#85, p = 0.05
#59, p = 0.01
head(conds)

# add annotation as row names and remove p-values and identifying info from dataframe
# make unique rows for p=0.2
rownames(conds) <- make.unique(conds$gene_name)
#row.names(conds)=conds$gene_name

exp = conds[, c(3:22)]
head(exp)

means=apply(exp,1,mean) # calculate means of rows
explc=exp-means # subtracting them
head(explc)

# reorganize column names to make heatmap clustered by treatment
explc = explc %>%
  select("Cool_1","Cool_2","Cool_3","Cool_4","Control_1.1","Control_1","Control_2.1","Control_2","Control_3.1","Control_3","Control_4.1","Control_4",
         "Heat_1.1","Heat_1","Heat_2.1","Heat_2","Heat_3.1","Heat_3","Heat_4.1","Heat_4")
# this explc object is what we can use to make our heatmap

# now make heatmap
# set color palette
col0=colorRampPalette(rev(c("#c51b7d","#e9a3c9","#ffffff", "#5ab4ac","#01665e")))(100)

# make treatment data frame
treatment = as.factor(sapply(strsplit(colnames(explc), split = "_"), "[[", 1)) %>%
  revalue(c("Control" = "control", "Cool" = "cold", "Heat" = "heat"))
expDesign = data.frame(colnames(explc), treatment)
expDesign = expDesign %>%
  column_to_rownames(var = "colnames.explc.")

expDesign$treatment = ordered(expDesign$treatment, levels = c("cold", "control", "heat"))
expDesign = expDesign[order(expDesign$treatment), , drop = FALSE]

# check that our input dataframes match before making heatmap
colnames(explc) == row.names(expDesign)

my_colour = list(treatment = c(cold = "#74c476", heat = "#fd8d3c", control = "#a6611a"))

#heat map of all photosynthesis genes and save
library(pheatmap)

photo.heatmap = pheatmap(explc, cluster_cols = FALSE, scale = "row", color = col0, annotation_col = expDesign, annotation_colors = my_colour, show_rownames = TRUE, show_colnames = FALSE, border_color = "NA")
ggsave(photo.heatmap, file = "/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/plots/culture_heatmap_photosynthesis_p=.01.pdf", width=9, height=7, units=c("in"), useDingbats=FALSE)


 #### Heat Maps of Interesting Stress GOs - Syms in Culture ####
# Stress genes
# read in go.obo database to figure out the go category id's that we want to pull.

go.obo = read.delim("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/go.obo")

# these terms are visually pulled from the GO tree outputs, all terms related to stress 
inds = c(which(go.obo$format.version..1.2 == "name: cellular response to chemical stress"), # BP
         which(go.obo$format.version..1.2 == "name: cellular response to oxidative stress"), # BP
         which(go.obo$format.version..1.2 == "name: hydrogen peroxide metabolic process"), # BP
         which(go.obo$format.version..1.2 == "name: proteasomal ubiquitin−independent protein catabolic process"), # BP
         which(go.obo$format.version..1.2 == "name: protein folding"), # BP
         which(go.obo$format.version..1.2 == "name: unfolded protein binding"), # MF
         which(go.obo$format.version..1.2 == "name: antioxidant"), # MF
         which(go.obo$format.version..1.2 == "name: oxidoreductase, acting on peroxide as acceptor") # MF
) 

want.go.obo = go.obo %>% 
  dplyr::slice(sort(unique(c(inds - 1)))) %>%
  mutate_at("format.version..1.2", str_replace, "id: ", "")

##pulling all of the stress terms and making a data frame
rldpval = read.csv("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/Sym_TagSeq/Culture_RLDandPVALS.csv", header = TRUE) %>%
  dplyr::rename("gene" = "X")

# remove isogroup from rldpval gene
rldpval$gene = gsub("isogroup", "", rldpval$gene)

# read in iso2go file for the algae
iso2go_sym = read.delim("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/Oculina_GO_Analyses/B_psygmophilum_isogroup_to_GOterm.tab", sep = "\t", header = FALSE) %>%
  dplyr::rename("gene" = "V1") %>%
  dplyr::rename("GO_ID" = "V2")

#pulling all of the photosynthesis terms from the iso2go 
list = as.list(want.go.obo)
print(list)

stress_genes = iso2go_sym %>%
  filter(grepl("GO:0006457|GO:0034599|GO:0042743|GO:0051082|GO:0062197", GO_ID)) %>%
  left_join(rldpval)
head(stress_genes)
str(stress_genes)

# add gene name to this df - trim down names and add to photo_genes
iso2gene = read.delim("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_TagSeq/References/B_psygmophilum_transcriptome/B_psygmophilum_isogroup_to_genename.tab", sep = "\t", header = FALSE) %>%
  dplyr::rename("gene" = "V1") %>%
  dplyr::rename("gene_name" = "V2")
head(iso2gene)

iso2gene$gene_name = gsub("OS=.*", "", iso2gene$gene_name)
iso2gene$gene_name = gsub(", .*","",as.character(iso2gene$gene_name))
head(iso2gene)
str(iso2gene)

stress_genes_anno = stress_genes %>%
  left_join(iso2gene)
head(stress_genes_anno)

# set raw p-value for GO enriched
p.val = 0.01

# filter based on p-values from deseq results
conds=stress_genes_anno[stress_genes_anno$pval.cold.culture<=p.val & !is.na(stress_genes_anno$pval.cold.culture),]
length(conds[,1])
#109, p = 0.1
#92, p = 0.05
#69, p = 0.01
head(conds)

# add annotation as row names and remove p-values and identifying info from dataframe
# remove NA gene names first
conds2 = conds %>%
  drop_na(gene_name)

# make unique rows for p=0.2
rownames(conds2) <- make.unique(conds2$gene_name)
#row.names(conds)=conds$gene_name

exp = conds2[, c(3:22)]
head(exp)

means=apply(exp,1,mean) # calculate means of rows
explc=exp-means # subtracting them
head(explc)

# reorganize column names to make heatmap clustered by treatment
explc = explc %>%
  select("Cool_1","Cool_2","Cool_3","Cool_4","Control_1.1","Control_1","Control_2.1","Control_2","Control_3.1","Control_3","Control_4.1","Control_4",
         "Heat_1.1","Heat_1","Heat_2.1","Heat_2","Heat_3.1","Heat_3","Heat_4.1","Heat_4")
# this explc object is what we can use to make our heatmap

# now make heatmap
# set color palette
#col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
col0=colorRampPalette(rev(c("#c51b7d","#e9a3c9","#ffffff", "#5ab4ac","#01665e")))(100)

# make treatment data frame
treatment = as.factor(sapply(strsplit(colnames(explc), split = "_"), "[[", 1)) %>%
  revalue(c("Control" = "control", "Cool" = "cold", "Heat" = "heat"))
expDesign = data.frame(colnames(explc), treatment)
expDesign = expDesign %>%
  column_to_rownames(var = "colnames.explc.")

expDesign$treatment = ordered(expDesign$treatment, levels = c("cold", "control", "heat"))

my_colour = list(treatment = c(cold = "#74c476", heat = "#fd8d3c", control = "#a6611a"))

#heat map of all photosynthesis genes and save
library(pheatmap)

stress.heatmap = pheatmap(explc, cluster_cols = FALSE, scale = "row", color = col0, annotation_col = expDesign, annotation_colors = my_colour, show_rownames = TRUE, show_colnames = FALSE, border_color = "NA")
ggsave(stress.heatmap, file = "/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/plots/culture_heatmap_stress_p=.01.pdf", width=9, height=8, units=c("in"), useDingbats=FALSE)
