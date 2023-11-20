#assemble_ortholog_counts.R
# this script takes the counts file output from the ortholog pipeline and compiles based on the shared ortholog id to create the counts file that we
# will input into deseq for comparative gene expression
rm(list=ls())
library(tidyverse)
library(data.table)

# This updated analysis includes all comparisons - (brown host, white host, syms in culture, syms in symbiosis)
setwd('~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/03_DESeq_Scripts/CountsFiles')


# READ IN THE ORTHOLOGS ---------------------------------------------------

#read in single copy orthologs
sOrthos = read_tsv('singleCopyOrthos.txt',
                   col_names=c('orthoGroup', 'contig'))

#strip away the .p suffixes
sOrthos$contig = sapply(sOrthos$contig, function(x) return(strsplit(x, '.p', fixed=TRUE)[[1]][1]))



# READ IN THE COUNTS ------------------------------------------------------

# read in full counts file
all_counts = read_tsv('all_ortho_counts.txt')
colnames(all_counts)[1]='contig'

# separate out host counts
hostCounts = all_counts[, c(1, grep("Oculina", colnames(all_counts)))] 
hostCounts = hostCounts[hostCounts$contig %like% "OARB", ]  

# re-name host columns based on treatment and sym state
names(hostCounts) = c("contig", "OA4_C_W", "OA5_F_W", "OA6_H_W", 
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

colnames(hostCounts)[2:48] <- paste(colnames(hostCounts)[2:48], "host", sep = ".")

hostCounts_clones_removed = hostCounts %>%
  select(-c(OP1_C_W.host, OP2_F_W.host, OP3_H_W.host, 
            OO1_C_W.host, OO2_F_W.host, 
            OL6_C_B.host, OL7_F_B.host, OL8_H_B.host, 
            OK1_C_W.host, OK2_F_W.host, OK3_H_W.host, 
            OA4_C_W.host, OA5_F_W.host, OA6_H_W.host))

# separate out brown and white host counts
#brownhostCounts = hostCounts[, c(1, grep("_B", colnames(hostCounts)))]
#head(brownhostCounts)

#whitehostCounts = hostCounts[, c(1, grep("_W", colnames(hostCounts)))]
#head(whitehostCounts)

# separate out symbiont in symbiosis counts from big counts dataset  
symCounts = all_counts[, c(1, grep("Oculina", colnames(all_counts)))] 
symCounts = symCounts[symCounts$contig %like% "BPSG", ]  

# rename columns
names(symCounts) = c("contig", "OA4_C_W", "OA5_F_W", "OA6_H_W", 
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

# select only samples from brown hosts, then add in .sym suffix to column names
symCounts = symCounts[, c(1, grep("_B", colnames(symCounts)))]
colnames(symCounts)[2:25] <- paste(colnames(symCounts)[2:25], "sym", sep = ".")
head(symCounts)

# remove clones that are removed in the host analyses.
symCounts_clones_removed = symCounts %>%
  select(-c(OL6_C_B.sym, OL7_F_B.sym, OL8_H_B.sym))

# separate out culture counts from bit counts dataset
#cultureCounts = all_counts[, c(1, grep("Culture", colnames(all_counts)))] 
#cultureCounts = cultureCounts[cultureCounts$contig %like% "BPSG", ]  

# remove annoying parts of the column names
#for ( col in 1:ncol(cultureCounts)){
#  colnames(cultureCounts)[col] <-  sub(".fastq.trim.sam.counts", "", colnames(cultureCounts)[col])
#}


dim(hostCounts_clones_removed)
# [1] 1917   34
dim(symCounts_clones_removed)
# [1] 1919   22
#dim(cultureCounts)
# [1] 1919   21

#filter contigs for mean read count greater than CUTOFF
CUTOFF=2
filter_by_mean = function(countdf){
  c = countdf %>% 
    dplyr::select(-contig)
  rmns = apply(c, 1, function(x) mean(x, na.rm=TRUE))
  keep = rmns > CUTOFF
  reduced = countdf[keep,]
  print('removed genes with mean count less than 2')
  print(paste('before =', nrow(countdf)))
  print(paste('after =', nrow(reduced)))
  return(reduced)
}

rhostCounts = filter_by_mean(hostCounts_clones_removed)
# [1] "removed genes with mean count less than 2"
# [1] "before = 1917"
# [1] "after = 1381"

rsymCounts = filter_by_mean(symCounts_clones_removed)
# [1] "removed genes with mean count less than 2"
# [1] "before = 1919"
# [1] "after = 250"

#rcultureCounts = filter_by_mean(cultureCounts)
# [1] "removed genes with mean count less than 2"
# [1] "before = 1919"
# [1] "after = 1306"


# MERGE WITH ORTHOLOGS ----------------------------------------------------

mhostCounts = rhostCounts %>% 
  inner_join(sOrthos, by='contig') %>%  #inner_join returns all rows from x where there are matching values in y, and all columns from x and y. 
  dplyr::select(-contig)

msymCounts = rsymCounts %>% 
  inner_join(sOrthos, by='contig') %>% 
  dplyr::select(-contig)

#mcultureCounts = rcultureCounts %>% 
#  inner_join(sOrthos, by='contig') %>% 
#  dplyr::select(-contig)


dim(mhostCounts) #[1] 1381   34
dim(msymCounts) #[1] 250  22
#dim(mcultureCounts) #[1] 1306   21

#now merge the datasets together
counts0 = mhostCounts %>% 
  inner_join(msymCounts, by = 'orthoGroup')
  #inner_join(mcultureCounts, by = 'orthoGroup')

nrow(counts0) #185

#put columns in order and convert to dataframe for DESEQ
counts = counts0[,order(colnames(counts0))] %>% 
  data.frame()
rownames(counts)=counts$orthoGroup
counts$orthoGroup<-NULL
write.table(counts, "OrthoCounts_deseq.txt", sep="\t") #this is the counts file we will want to come back to


#build coldata - ended up doing this by hand in excel based on the column names from OrthoCounts_deseq.txt


# RUN DESEQ TO CHECK EXPRESSION AGREEMENT ---------------------------------
library(DESeq2)

expDesign = read.csv("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/08_tables/ExpDesign_orthos_noculture.csv")
expDesign$treat_type = as.factor(paste(expDesign$temp,expDesign$type, sep = "_"))

#run DESeq
dds<-DESeqDataSetFromMatrix(countData=counts,
                            colData = expDesign,
                            design =~ treat_type)
#dds <- DESeq(dds,
#             fitType='local')
dds <- DESeq(dds)
resultsNames(dds)
colData(dds)
head(dds)
#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot Oculina")

#rlog transformation
rlogged = rlogTransformation(dds, blind = TRUE)

# DIFFERENTIAL EXPRESSION ------------
#compare host heat to host control
resHostH.HostC=results(dds, contrast=c("treat_type","heat_brown_host","control_brown_host"), independentFiltering = F)
resHostH.HostC
table(resHostH.HostC$pvalue<0.05)
#FALSE  TRUE 
#151    34
table(resHostH.HostC$padj<0.05)
#FALSE  TRUE 
#177    8 
summary(resHostH.HostC)

#can make a volcano plot with this code:
resHostH.HostC %>% 
  data.frame() %>% 
  mutate(sig=padj<0.1) %>% 
  ggplot(aes(x=log2FoldChange, y=-log(pvalue), color=sig)) +
  geom_point(alpha=0.3) +
  scale_color_manual(values=c('black', 'red'))

#MA plot
plotMA(resHostH.HostC, main="Host Heat v. Control")

resHostH.sig = subset(resHostH.HostC, resHostH.HostC$padj<0.05)
#write.table(resHostH.sig, "HostHeatSigOrthos.txt", sep="\t")

#compare host freezing to host control
resHostF.HostC=results(dds, contrast=c("treat_type","cold_brown_host","control_brown_host"), independentFiltering = F)
resHostF.HostC
table(resHostF.HostC$pvalue<0.05)
#FALSE  TRUE 
#114   71 
table(resHostF.HostC$padj<0.05)
#FALSE  TRUE 
#129    56 
summary(resHostF.HostC)

resHostF.sig = subset(resHostF.HostC, resHostF.HostC$padj<0.05)
#write.table(resHostF.sig, "HostFreezingSigOrthos.txt", sep="\t")

#compare symbiont hot to symbiont control
resSymH.SymC=results(dds, contrast=c("treat_type","heat_sym_inhost","control_sym_inhost"), independentFiltering = F)
resSymH.SymC
table(resSymH.SymC$pvalue<0.05)
#FALSE  TRUE 
#162    23 
table(resSymH.SymC$padj<0.05)
#FALSE  TRUE 
#179     6 
summary(resSymH.SymC)

resSymH.sig = subset(resSymH.SymC, resSymH.SymC$padj<0.05)
#write.table(resSymH.sig, "SymHeatSigOrthos.txt", sep="\t")

#compare symbiont freezing to symbiont control
resSymF.SymC=results(dds, contrast=c("treat_type","cold_sym_inhost","control_sym_inhost"), independentFiltering = F)
resSymF.SymC
table(resSymF.SymC$pvalue<0.05)
#FALSE  TRUE 
#173    12 
table(resSymF.SymC$padj<0.05)
#FALSE  TRUE 
#181     4 
summary(resSymF.SymC)

resSymF.sig = subset(resSymF.SymC, resSymF.SymC$padj<0.05)
write.table(resSymF.SymC, "SymFreezingSigOrthos.txt", sep="\t")

#compare culture hot to culture control
#resCultureH.CultureC=results(dds, contrast=c("treat_type","heat_culture","control_culture"), independentFiltering = F)
#resCultureH.CultureC
#table(resCultureH.CultureC$pvalue<0.05)
#FALSE  TRUE 
#117    58 
#table(resCultureH.CultureC$padj<0.05)
#FALSE  TRUE 
#136     39 
#summary(resCultureH.CultureC)

#resCultureH.sig = subset(resCultureH.CultureC, resCultureH.CultureC$padj<0.05)
#write.table(resCultureH.CultureC, "SymHeatSigOrthos.txt", sep="\t")

#compare symbiont freezing to symbiont control
#resCultureF.CultureC=results(dds, contrast=c("treat_type","cold_culture","control_culture"), independentFiltering = F)
#resCultureF.CultureC
#table(resCultureF.CultureC$pvalue<0.05)
#FALSE  TRUE 
#120    55 
#table(resCultureF.CultureC$padj<0.05)
#FALSE  TRUE 
#149     26 
#summary(resCultureF.CultureC)

#resCultureF.sig = subset(resCultureF.CultureC, resCultureF.CultureC$padj<0.05)
#write.table(resCultureF.CultureC, "SymFreezingSigOrthos.txt", sep="\t")

# GET PVALS AND MAKE TABLE ------------
head(resHostH.HostC) #host heat vs. control
valsHostH=cbind(resHostH.HostC$pvalue, resHostH.HostC$padj)
head(valsHostH)
colnames(valsHostH)=c("pval.hostH", "padj.hostH")
length(valsHostH[,1])
table(complete.cases(valsHostH))

head(resHostF.HostC) #host freezing v. control
valsHostF=cbind(resHostF.HostC$pvalue, resHostF.HostC$padj)
head(valsHostF)
colnames(valsHostF)=c("pval.hostF", "padj.hostF")
length(valsHostF[,1])
table(complete.cases(valsHostF))

head(resSymH.SymC) #sym heat vs. control
valsSymH=cbind(resSymH.SymC$pvalue, resSymH.SymC$padj)
head(valsSymH)
colnames(valsSymH)=c("pval.symH", "padj.symH")
length(valsSymH[,1])
table(complete.cases(valsSymH))

head(resSymF.SymC) #sym freezing vs. control
valsSymF=cbind(resSymF.SymC$pvalue, resSymF.SymC$padj)
head(valsSymF)
colnames(valsSymF)=c("pval.symF", "padj.symF")
length(valsSymF[,1])
table(complete.cases(valsSymF))

#head(resCultureH.CultureC) #culture heat vs. control
#valsCultureH=cbind(resCultureH.CultureC$pvalue, resCultureH.CultureC$padj)
#head(valsCultureH)
#colnames(valsCultureH)=c("pval.symH", "padj.symH")
#length(valsCultureH[,1])
#table(complete.cases(valsCultureH))

#head(resCultureF.CultureC) #culture freezing vs. control
#valsCultureF=cbind(resCultureF.CultureC$pvalue, resCultureF.CultureC$padj)
#head(valsCultureF)
#colnames(valsCultureF)=c("pval.symF", "padj.symF")
#length(valsCultureF[,1])
#table(complete.cases(valsCultureF))

#rlog transform data
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
length(rld[,1])


# GET VARIANCE STABILIZED COUNTS AND CHECK ORTHLOG CORRELATION ------------
#idea here is that if they are really the 'same' genes, their expression values
#should be reasonably well-correlated

vsd = varianceStabilizingTransformation(dds,
          fitType='local')

vsdf = data.frame(assay(vsd))
hostdat = vsdf %>% 
  dplyr::select(contains('.host'))
symdat = vsdf %>% 
  dplyr::select(contains('.sym'))
#culturedat = vsdf  
#  dplyr::select(contains('Culture'))

mnhost = apply(hostdat, 1, function(x) mean(x, na.rm=TRUE))
mnsym = apply(symdat, 1, function(x) mean(x, na.rm=TRUE))
#mnculture = apply(culturedat, 1, function(x) mean(x, na.rm=TRUE))

mdat = data.frame(mnhost, mnsym)
mdat %>% 
  ggplot(aes(x=mnhost,y=mnsym)) +
  geom_point(alpha=0.3)
lm1=lm(mnhost~mnsym)
summary(lm1)
cor(x=mnhost, y=mnsym)

