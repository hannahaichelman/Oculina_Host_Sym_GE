# Set working directory
setwd("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/09_Supplemental_Info/SNPCalling/Oculina")
#clustering / PCoA based on identity by state (IBS) 
#(for low and/or uneven coverage)


#First we start by reading in the bams files.
#Change what files you are working with to plot the cluster dendrogram with and without genotype A

#bams=read.table("bams_brown")[,1] # list of bam files for symbiotic/brown oculina individuals only
bams=read.table("bams")[,1] # list of bam files for all oculina individuals
bams.noA=read.table("bams_noA")[,1] # list of bam files for all oculina individuals except genotype A that looks weird in dendrogram

goods=c(1:length(bams))
goods.noA=c(1:length(bams.noA))

#ma = as.matrix(read.table("Oarb_browns.ibsMat")) # identity by state matrix for symbiotic/brown oculina individuals only
ma = as.matrix(read.table("Oarb_newref.ibsMat")) # all oculina individuals
ma.noA = as.matrix(read.table("Oarb_newref_noA.ibsMat")) # all oculina individuals except genotype A

ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])

ma.noA=ma.noA[goods.noA,goods.noA]
dimnames(ma.noA)=list(bams.noA[goods.noA],bams.noA[goods.noA])

# trim column and row names to make the cluster dendrogram nicer to read
colnames(ma) = c("OA4_C_W", "OA5_F_W", "OA6_H_W", 
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

row.names(ma) = c("OA4_C_W", "OA5_F_W", "OA6_H_W", 
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

# now plot dendrograms
# first with all samples
hc=hclust(as.dist(ma),"ave")
hc.plot = plot(hc,cex=0.7) 

# now plot without genotype A
hc.noA=hclust(as.dist(ma.noA),"ave")
plot(hc.noA,cex=0.7) 
abline(h=0.19, col='darkgray', lty = 2, lwd=2) # add a vertical line to indicate cut-off for clones
# This cut-off aligns well with what James Fifer has used in the past (0.15), 
# this means that any samples that cluster together below that threshold are clones

# Sample names are read as follows: 
# OI3_H_B = Oculina, genet I, fragment 3, heat treatment, brown/symbiotic
# This analysis indicates that for brown/symbiotic individuals, we only have one unexpected clonal group (M and L)
# When we look at all individuals, we see that genotype A arranges strangely across the dendrogram, we remove from analysis for this reason
# For white/aposymbiotic individuals, we see that there are two more clonal groups: N+O+P and H+K.
