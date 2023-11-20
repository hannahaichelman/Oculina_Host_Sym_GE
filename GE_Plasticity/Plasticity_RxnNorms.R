library(tidyverse)
library(ggplot2)
library(stringr)
library(ggpubr)


# This script takes in the plasticity scores output from Colleen Bove's PCA Plasticity function 
# and plots reaction norms of PCA plasticity

setwd("~/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/GE_Plasticity")

# read in host plasticity, make columns for genotype and fragment
host_plast = read_csv("host_plasticity.csv") %>%
  select(-...1) %>%
  select(name, treat_type, sym_state, dist) %>%
  filter(sym_state == "sym") %>%
  mutate(genotype = str_extract(.$name, "^.{2}")) %>%
  mutate(fragment = str_extract(.$name, "^.{3}")) %>%
  mutate(treatment = str_sub(treat_type, 1, 4)) %>%
  mutate_at(c('treat_type','sym_state','treatment'), as.factor) %>%
  filter(treat_type != "control_sym") %>%
  filter(treat_type != "control_apo") %>%
  select(-treat_type, -sym_state) %>%
  select(name, treatment, genotype, fragment, dist) %>%
  mutate(type = "host")

head(host_plast)

# read in symbiont in host plasticity, make columns for genotype and fragment
syminhost_plast = read_csv("symbiont_inhost_plasticity.csv") %>%
  select(-...1) %>%
  select(name, treatment, dist) %>%
  mutate(genotype = str_extract(.$name, "^.{2}")) %>%
  mutate(fragment = str_extract(.$name, "^.{3}")) %>%
  mutate_at(c('treatment'), as.factor) %>%
  filter(treatment != "control") %>%
  select(name, treatment, genotype, fragment, dist) %>%
  mutate(type = "sym_in_host")

head(syminhost_plast)

# read in symbiont in culture and symbiont in host combined data plasticity, make columns for genotype and fragment
culturesym_plast = read_csv("culturesymcombined_plasticity.csv") %>%
  select(-...1) %>%
  mutate(genotype = str_extract(.$name, "^.{2}")) %>%
  mutate(fragment = str_extract(.$name, "^.{3}")) %>%
  mutate(treatment = str_sub(treat_type, 1, 4)) %>%
  mutate_at(c('treat_type','type','treatment'), as.factor)


# remove plasticity values in control treatments
culturesym_plast = culturesym_plast %>%
  filter(treat_type != "control_sym") %>%
  filter(treat_type != "control_apo")


# plot host and sym in host data together on reaction norm of gene expression plasticity
host_sym_plast = rbind(host_plast, syminhost_plast)

host_plast_heat = host_plast %>%
  filter(treatment == "heat")

host_plast_cold = host_plast %>%
  filter(treatment == "cold")

syminhost_plast_heat = syminhost_plast %>%
  filter(treatment == "heat")

syminhost_plast_cold = syminhost_plast %>%
  filter(treatment == "cold")

host_sym_plast_heat = host_sym_plast %>%
  filter(treatment == "heat")

host_sym_plast_cold = host_sym_plast %>%
  filter(treatment == "cold")

plast_plot_heat = ggplot(host_sym_plast_heat, aes(x = type, y = dist, group = genotype)) +
  geom_point(aes(color = genotype), size = 2) +
  geom_line(aes(group = genotype)) +
  ylab("Gene Expression Plasticity") +
  xlab("Data Type") +
  ggtitle("Heat Challenge")
plast_plot_heat

plast_plot_cold = ggplot(host_sym_plast_cold, aes(x = type, y = dist, group = genotype)) +
  geom_point(aes(color = genotype), size = 2) +
  geom_line(aes(group = genotype)) +
  ylab("Gene Expression Plasticity") +
  xlab("Data Type") +
  ggtitle("Cold Challenge")
plast_plot_cold


plast_combined = ggarrange(plast_plot_heat, plast_plot_cold, 
                           nrow = 1,
                           common.legend = TRUE, legend = "right")
plast_combined


plot(x = host_plast_heat$dist, y = syminhost_plast_heat$dist,
     xlab = "Host Gene Expression Plasticity",
     ylab = "Symbiont Gene Expression Plasticity",
     main = "Oculina Heat Challenge")

plot(x = host_plast_cold$dist, y = syminhost_plast_cold$dist,
     xlab = "Host Gene Expression Plasticity",
     ylab = "Symbiont Gene Expression Plasticity",
     main = "Oculina Cold Challenge")



    