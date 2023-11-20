library(adegenet)
library(tidyverse)
library(plyr)
library(DESeq2)


#### Host Data ####
countDataHost <- read.table("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/CountsFiles/Oculina_counts_newref_host.txt")

names(countDataHost) = c("OA4_C_W", "OA5_F_W", "OA6_H_W", 
                         "OC4_F_B",	"OC5_H_B",	"OC9_C_B",	
                         "OD4_C_B",	"OD5_F_B",	"OD6_H_B",
                         "OE10_F_W", "OE11_C_W", "OE3_H_W",
                         "OF7_C_B",	"OF8_F_B",	"OF9_H_B",
                         "OH11_F_W", "OH15_H_W", "OH1_C_W",
                         "OI1_C_B",	"OI2_F_B",	"OI3_H_B",	
                         "OJ13_C_B",	"OJ14_F_B",	"OJ15_H_B",	
                         "OK1_C_W", "OK2_F_W", "OK3_H_W",
                         "OL6_C_B",	"OL7_F_B",	"OL8_H_B",	
                         "OM1_C_B",	"OM2_F_B",	"OM3_H_B",
                         "ON4_C_W", "ON5_F_W", "ON6_H_W",
                         "OO1_C_W", "OO2_F_W",
                         "OP1_C_W", "OP2_F_W", "OP3_H_W",
                         "OQ11_H_W", "OQ1_C_W", "OQ4_F_W",
                         "OR7_C_B",	"OR8_F_B",	"OR9_H_B")

treatmentHost = as.factor(sapply(strsplit(colnames(countDataHost), split = "_"), "[[", 2)) %>%
  revalue(c("C" = "control", "F" = "cold", "H" = "heat"))

genotypeHost  = as.factor(sapply(strsplit(colnames(countDataHost), split = ""), "[[", 2))

symState  = as.factor(sapply(strsplit(colnames(countDataHost), split = "_"), "[[", 3)) %>%
  revalue(c("B" = "sym", "W" = "apo"))

expDesign_Host = data.frame(colnames(countDataHost), treatmentHost, genotypeHost, symState)
expDesign_Host$type = "host"
expDesign_Host$treat_type = as.factor(paste(expDesign_Host$treatmentHost,expDesign_Host$symState, sep = "_"))
names(expDesign_Host) = c("sample", "treatment", "genotype", "sym_state", "type", "treat_type")
str(expDesign_Host)

colSums(countDataHost)

# remove clones and outlier genotype A from experimental design and counts
expDesign_Host_clones_removed = expDesign_Host %>%
  filter(genotype!="P") %>%
  filter(genotype!="O") %>%
  filter(genotype!="L") %>%
  filter(genotype!="K") %>%
  filter(genotype!="A") #strange on cluster dendrogram

countDataHost_clones_removed = countDataHost %>%
  select(-c(OP1_C_W, OP2_F_W, OP3_H_W, OO1_C_W, OO2_F_W, OL6_C_B, OL7_F_B, OL8_H_B, OK1_C_W, OK2_F_W, OK3_H_W, OA4_C_W, OA5_F_W, OA6_H_W))

# Run DESeq
dds<-DESeqDataSetFromMatrix(countData=countDataHost_clones_removed, colData=expDesign_Host_clones_removed, design=~treat_type) #can only test for the main effects of treatment
dds = DESeq(dds)

# rlog transformation
rlogged = rlog(dds, blind = TRUE)

dat_rlog = as.data.frame(assay(rlogged))

# Now run DAPC - look for optimal number of clusters
clus = find.clusters(t(dat_rlog), max.n.clust=30) # 20 pc's and 6 clusters
names(clus)
clus$grp = expDesign_Host_clones_removed$treat_type
#head(clus$grp, 6)
#clus$size
# 5 3 5 6 5 9

# use a-score to find the optimal number of PCs
dapc1 = dapc(t(dat_rlog), n.da = 100, n.pca = 20, clus$grp)
temp = optim.a.score(dapc1) # optimal number of PCs = 10


# colours 

myCol = c("#9ecae1", "#3182bd","lightgrey","grey","#fc9272", "#de2d26")

shapes = c(21, 19, 21, 19, 21, 19)

# now lets build a discriminant function for these six groups:
dp=dapc(t(dat_rlog),clus$grp) # 2pcs, 6 discriminant functions

pdf("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/GE_Plasticity/dapc_2PC_6DF_host.pdf")
scatter(dp,1,1, bg = "white", legend = TRUE, posi.leg = "topleft", col = myCol, pch = shapes, solid = 0.8)
dev.off()


pdf("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/GE_Plasticity/dapc_pca_2PC_6DF_host.pdf")
scatter(dp, bg = "white", scree.da = TRUE, legend = TRUE, posi.leg = "topleft",
        solid = 0.8, col = myCol, scree.pca = TRUE,
        clabel = 0)
dev.off()

#### Symbiont in Host Data ####

countDataSym <- read.table("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/CountsFiles/Oculina_Counts_newref_sym.txt")

names(countDataSym) = c("OA4", "OA5", "OA6", 
                        "OC4_F_B",	"OC5_H_B",	"OC9_C_B",	
                        "OD4_C_B",	"OD5_F_B",	"OD6_H_B",
                        "OE10", "OE11", "OE3",
                        "OF7_C_B",	"OF8_F_B",	"OF9_H_B",
                        "OH11", "OH15", "OH1",
                        "OI1_C_B",	"OI2_F_B",	"OI3_H_B",	
                        "OJ13_C_B",	"OJ14_F_B",	"OJ15_H_B",	
                        "OK1", "OK2", "OK3",
                        "OL6_C_B",	"OL7_F_B",	"OL8_H_B",	
                        "OM1_C_B",	"OM2_F_B",	"OM3_H_B",
                        "ON4", "ON5", "ON6",
                        "OO1", "OO2",
                        "OP1", "OP2", "OP3",
                        "OQ11", "OQ1", "OQ4",
                        "OR7_C_B",	"OR8_F_B",	"OR9_H_B")

countDataSym_brown = countDataSym %>%
  select(-c(OA4, OA5, OA6, OE10, OE11, OE3, OH11, OH15, OH1, OK1, OK2, OK3, ON4, ON5, ON6, OO1, OO2, OP1, OP2, OP3, OQ11, OQ1, OQ4))

# We found a pair of clones in the symbiotic Oculina dataset, M's and L's, so removing same clone pair from the host and other analyses (L's)
countDataSym_brown_noclone = countDataSym_brown %>%
  select(-c(OL6_C_B,  OL7_F_B,  OL8_H_B))

treatmentSym = as.factor(sapply(strsplit(colnames(countDataSym_brown_noclone), split = "_"), "[[", 2)) %>%
  revalue(c("C" = "control", "F" = "cold", "H" = "heat"))

genotypeSym  = as.factor(sapply(strsplit(colnames(countDataSym_brown_noclone), split = ""), "[[", 2))

expDesign_Sym = data.frame(colnames(countDataSym_brown_noclone), treatmentSym, genotypeSym)
expDesign_Sym$type = "sym_inhost"
expDesign_Sym$treat_type = as.factor(paste(expDesign_Sym$treatmentSym,expDesign_Sym$type, sep = "_"))
names(expDesign_Sym) = c("sample", "treatment", "genotype", "type", "treat_type")
str(expDesign_Sym)

dds_Sym<-DESeqDataSetFromMatrix(countData=countDataSym_brown_noclone, colData=expDesign_Sym, design=~treatment) #can only test for the main effects of treatment
dds_Sym = DESeq(dds_Sym)

# rlog transformation
rlog_syminhost = rlog(dds_Sym, blind = TRUE)

dat_rlog_syminhost = as.data.frame(assay(rlog_syminhost))

# Now run DAPC - look for optimal number of clusters
clus_syminhost = find.clusters(t(dat_rlog_syminhost), max.n.clust=20) # 15 pc's and 3 clusters
names(clus_syminhost)
clus_syminhost$grp = expDesign_Sym$treatment

# use a-score to find the optimal number of PCs
dapc1_syminhost = dapc(t(dat_rlog_syminhost), n.da = 100, n.pca = 20, clus_syminhost$grp)
temp_syminhost = optim.a.score(dapc1_syminhost) # optimal number of PCs = 1


# colours 

cols_sym = c("control" = "#a6611a", "cold" = "#74c476", "heat" = "#fd8d3c")

# now lets build a discriminant function for these six groups:
dp_syminhost=dapc(t(dat_rlog_syminhost),clus_syminhost$grp) # 2pcs, 3 discriminant functions

pdf("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/GE_Plasticity/dapc_2PC_3DF_syminhost.pdf")
scatter(dp_syminhost,1,1, bg = "white", legend = TRUE, posi.leg = "topright", col = cols_sym, solid = 0.8)
dev.off()


pdf("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/GE_Plasticity/dapc_pca_2PC_3DF_syminhost.pdf")
scatter(dp_syminhost, bg = "white", scree.da = TRUE, legend = TRUE, posi.leg = "topleft",
        solid = 0.8, col = cols_sym, scree.pca = TRUE,
        clabel = 0)
dev.off()


#### Symbiont in Culture Data ####

countDataCulture <- read.table("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/CountsFiles/B_psygmophilum_counts.txt")

names(countDataCulture)=c( "Control_1.1", "Control_1", "Control_2.1", "Control_2", "Control_3.1", "Control_3", "Control_4.1", "Control_4", "Cool_1", "Cool_2", "Cool_3", "Cool_4","Heat_1.1", "Heat_1", "Heat_2.1", "Heat_2", "Heat_3.1", "Heat_3", "Heat_4.1", "Heat_4")
row.names(countDataCulture)=sub("", "isogroup", rownames(countDataCulture))
head(countDataCulture)

treatment = as.factor(sapply(strsplit(colnames(countDataCulture), split = "_"), "[[", 1)) %>%
  revalue(c("Control" = "control", "Cool" = "cold", "Heat" = "heat"))

expDesignCulture = data.frame(treatment) 
expDesignCulture$sample = colnames(countDataCulture)

str(expDesignCulture)
dds_Culture<-DESeqDataSetFromMatrix(countData=countDataCulture, colData=expDesignCulture, design=~treatment) #can only test for the main effects of treatment
dds_Culture = DESeq(dds_Culture)

rlog_culture = rlogTransformation(dds_Culture, blind = TRUE)

dat_rlog_culture = as.data.frame(assay(rlog_culture))

# Now run DAPC - look for optimal number of clusters
clus_culture = find.clusters(t(dat_rlog_culture), max.n.clust=15) # 15 pc's and 3 clusters
names(clus_culture)
clus_culture$grp = expDesignCulture$treatment

# use a-score to find the optimal number of PCs
dapc1_culture = dapc(t(dat_rlog_culture), n.da = 100, n.pca = 15, clus_culture$grp)
temp_culture = optim.a.score(dapc1_culture) # optimal number of PCs = 1


# colours 

cols_culture = c("control" = "#a6611a", "cold" = "#74c476", "heat" = "#fd8d3c")

# now lets build a discriminant function for these six groups:
dp_culture=dapc(t(dat_rlog_culture),clus_culture$grp) # 2pcs, 3 discriminant functions

pdf("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/GE_Plasticity/dapc_2PC_3DF_syminculture.pdf")
scatter(dp_culture,1,1, bg = "white", legend = TRUE, posi.leg = "top", col = cols_culture, solid = 0.8)
dev.off()


pdf("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/GE_Plasticity/dapc_pca_2PC_3DF_syminculture.pdf")
scatter(dp_culture, bg = "white", scree.da = TRUE, legend = TRUE, posi.leg = "topleft",
        solid = 0.8, col = cols_culture, scree.pca = TRUE,
        clabel = 0)
dev.off()


