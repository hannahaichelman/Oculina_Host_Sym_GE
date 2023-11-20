# This script analyzes the cell count data from the Breviolum psygmophilum culture experiment. 
# Author: Hannah E Aichelman
# Last Updated: January 1, 2023

#### Load in and format data ####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Rmisc)

setwd("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/09_Supplemental_Info/Culture_Growth")

# We have cell count data from two separate experiments, the test experiment to determine the timing of stationary growth phase (test_counts), 
# and the final experiment where we sampled cells for gene expression (exp_counts). 

test_counts = read.csv("SemiContinuousTestExperiment_CellCounts.csv")
str(test_counts)
test_counts$Day = as.numeric(test_counts$Day)
test_counts$Flask_ID = as.factor(test_counts$Flask_ID)

exp_counts = read.csv("FinalExperiment_CellCounts.csv")
str(exp_counts)
exp_counts$Day = as.numeric(exp_counts$Day)
exp_counts$Flask_ID = as.factor(exp_counts$Flask_ID)



#### Plot Semi Continuous Test Experiment Cell Counts ####
head(test_counts)
test_summary = summarySE(data = test_counts, measurevar = "Cells_Per_mL", groupvars = c("Day", "Flask_ID"))

flask_cols = c("#addd8e", "#31a354", "#006837")

test.plot = ggplot(test_summary, aes(x = Day, fill = Flask_ID))+
  theme_bw()+
  geom_errorbar(aes(x = Day, ymax = Cells_Per_mL+se, ymin = Cells_Per_mL-se), width = .3, position = position_dodge(width=0.5))+
  geom_point(aes(y = Cells_Per_mL), color = "black", shape = 21, size = 3, position = position_dodge(width = 0.5))+
  scale_y_continuous(name = "Cell Density (cells mL-1)")+
  scale_fill_manual(values = flask_cols)+
  scale_x_continuous(breaks = c(3,5,7,10,12,14,17,19,21,25)) +
  ggtitle("A. Test Experiment") +
  theme(legend.position = c(.8, .15), legend.background = element_rect(color = "black"))
test.plot

ggsave(test.plot, file = "test.culture.growth.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)



#### Plot Final Experiment Cell Counts ####
head(exp_counts)
exp_summary = summarySE(data = exp_counts, measurevar = "Cell_Per_mL", groupvars = c("Day", "Temp_C", "Flask_ID"))

exp_summary$Flask_ID = substring(exp_summary$Flask_ID, 5)

exp.plot = ggplot(exp_summary, aes(x = Day, fill = Flask_ID))+
  theme_bw()+
  geom_errorbar(aes(x = Day, ymax = Cell_Per_mL+se, ymin = Cell_Per_mL-se), width = .3, position = position_dodge(width=0.9))+
  geom_point(aes(y = Cell_Per_mL, fill = Flask_ID), color = "black", size = 3, shape = 21, position = position_dodge(width = 0.9))+
  scale_y_continuous(name = "Cell Density (cells mL-1)")+
  #scale_color_manual(values = flask_cols)+
  scale_x_continuous(breaks = c(4,6,8,11,13,15)) +
  ggtitle("C. Full Experiment - All Flasks") 
exp.plot

ggsave(exp.plot, file = "exp.culture.growth.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)


# now average over flask id too
head(exp_counts)
exp_counts$treatment = as.factor(sapply(strsplit(as.character(exp_counts$Flask_ID), split = "_"), "[[", 2)) %>%
  revalue(c("Control" = "control", "Cool" = "cold", "Heat" = "heat"))

exp_summary2 = summarySE(data = exp_counts, measurevar = "Cell_Per_mL", groupvars = c("Day", "treatment"))

# set colors
cols_sym = c("control" = "#a6611a", "cold" = "#74c476", "heat" = "#fd8d3c")

exp.plot2 = ggplot(exp_summary2, aes(x = Day, fill = treatment))+
  theme_bw()+
  geom_errorbar(aes(x = Day, ymax = Cell_Per_mL+se, ymin = Cell_Per_mL-se), width = .3, position = position_dodge(width=0.9))+
  geom_point(aes(y = Cell_Per_mL, fill = treatment), color = "black", size = 3, shape = 21, position = position_dodge(width = 0.9))+
  scale_y_continuous(name = "Cell Density (cells mL-1)")+
  scale_fill_manual(values = cols_sym)+
  scale_x_continuous(breaks = c(4,6,8,11,13,15)) +
  ggtitle("B. Full Experiment") +
  theme(legend.position = c(.2, .15), legend.background = element_rect(color = "black"))
exp.plot2

ggsave(exp.plot2, file = "exp.culture.growth.collapsed.pdf", width=4, height=4, units=c("in"), useDingbats=FALSE)

