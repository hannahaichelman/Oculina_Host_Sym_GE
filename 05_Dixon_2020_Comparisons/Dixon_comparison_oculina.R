### Comparing the GO Delta Ranks to the stress module by Dixon et al. 2020
library(tidyverse)
library(patchwork)
library(cowplot)
library(ggpubr)
library(plotly)


# Read in Dixon's 'type A' cluster for each GO category
# downloaded from here: https://github.com/grovesdixon/Acropora_gene_expression_meta/tree/master/go_mwu
Dixon_BP = read.table("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/05_Dixon_2020_Comparisons/MWU_BP_clusterAstress_For_MWU.csv", header = T)
Dixon_MF = read.table("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/05_Dixon_2020_Comparisons/MWU_MF_clusterAstress_For_MWU.csv", header = T)
Dixon_CC = read.table("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/05_Dixon_2020_Comparisons/MWU_CC_clusterAstress_For_MWU.csv", header = T)

#### Biological Processes ####
# Read in MWU GO files
Ocu_SymHost_Heat_BP = read.table("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/07_GeneOntology_Analyses/GO_MWU_output_files/MWU_BP_oculina_hot_brown_host_GO.csv", header = T)
Ocu_SymHost_Cold_BP = read.table("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/07_GeneOntology_Analyses/GO_MWU_output_files/MWU_BP_oculina_cold_brown_host_GO.csv", header = T)
Ocu_ApoHost_Heat_BP = read.table("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/07_GeneOntology_Analyses/GO_MWU_output_files/MWU_BP_oculina_hot_white_host_GO.csv", header = T)
Ocu_ApoHost_Cold_BP = read.table("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/07_GeneOntology_Analyses/GO_MWU_output_files/MWU_BP_oculina_cold_white_host_GO.csv", header = T)


## heat symbiotic host
# intersect() makes comparisons are made row-wise, so that in the data frame case, intersect(A,B) is a data frame with those rows that are both in A and in B.
heat_symhost_goods = intersect(Ocu_SymHost_Heat_BP$term, Dixon_BP$term)
heat_sym = Ocu_SymHost_Heat_BP[Ocu_SymHost_Heat_BP$term %in% heat_symhost_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% heat_symhost_goods,]

# Combine them
heat_symhost_BP_data = merge(heat_sym, Dixon_BP_set, by="term")

heat_symhost_BP_plot = ggplot(heat_symhost_BP_data, aes(delta.rank.x, delta.rank.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,15)) +
  geom_hex() +
  labs( x = "Symbiotic Host",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  xlim(c(-6000,6000)) +
  theme_bw()
heat_symhost_BP_plot


# look for most correlated categories
# can id delta ranks using this plotly:
heat_symhost_plotly = ggplot(heat_symhost_BP_data, aes(x = delta.rank.x, y = delta.rank.y)) +
  geom_point()
ggplotly(heat_symhost_plotly)

# and then filter data to figure out who is who from the plotly:
heat_symhost_BP_neat_up = heat_symhost_BP_data %>%
  filter(delta.rank.x > 3000 & delta.rank.y > 1000)

heat_symhost_BP_neat_down = heat_symhost_BP_data %>%
  filter(delta.rank.x < -3000) %>%
  filter(delta.rank.y < -800)


## heat aposymbiotic host
heat_apohost_goods = intersect(Ocu_ApoHost_Heat_BP$term, Dixon_BP$term)
heat_apo = Ocu_ApoHost_Heat_BP[Ocu_ApoHost_Heat_BP$term %in% heat_apohost_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% heat_apohost_goods,]

# Combine them
heat_apohost_BP_data = merge(heat_apo, Dixon_BP_set, by="term")

heat_apohost_BP_plot = ggplot(heat_apohost_BP_data, aes(delta.rank.x, delta.rank.y)) +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,15)) +
  geom_hex() +
  labs( x = "Aposymbiotic Host",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, lty = 2, color="black", se =FALSE) +
  xlim(c(-6000,6000)) +
  theme_bw()
heat_apohost_BP_plot

# look for most correlated categories
# can id delta ranks using this plotly:
heat_apohost_plotly = ggplot(heat_apohost_BP_data, aes(x = delta.rank.x, y = delta.rank.y)) +
  geom_point()
ggplotly(heat_apohost_plotly)

# and then filter data to figure out who is who:
heat_apohost_BP_neat_up = heat_apohost_BP_data %>%
  filter(delta.rank.x > 3000 & delta.rank.y > 1000)

heat_apohost_BP_neat_down = heat_apohost_BP_data %>%
  filter(delta.rank.x < -3000) %>%
  filter(delta.rank.y < -1000)

## cold symbiotic host
cold_symhost_goods = intersect(Ocu_SymHost_Cold_BP$term, Dixon_BP$term)
cold_sym = Ocu_SymHost_Cold_BP[Ocu_SymHost_Cold_BP$term %in% cold_symhost_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% cold_symhost_goods,]

# Combine them
cold_symhost_BP_data = merge(cold_sym, Dixon_BP_set, by="term")

cold_symhost_BP_plot = ggplot(cold_symhost_BP_data, aes(delta.rank.x, delta.rank.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,15)) +
  geom_hex() +
  labs( x = "Symbiotic Host",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, color="black", se =FALSE) +
  xlim(c(-6000,6000)) +
  theme_bw()
cold_symhost_BP_plot

# look for most correlated categories
# can id delta ranks using this plotly:
cold_symhost_plotly = ggplot(cold_symhost_BP_data, aes(x = delta.rank.x, y = delta.rank.y)) +
  geom_point()
ggplotly(cold_symhost_plotly)

# and then filter data to figure out who is who:
cold_symhost_BP_neat_up = cold_symhost_BP_data %>%
  filter(delta.rank.x > 4000 & delta.rank.y > 900)

cold_symhost_BP_neat_down = cold_symhost_BP_data %>%
  filter(delta.rank.x < -3000) %>%
  filter(delta.rank.y < -1000)

## cold aposymbiotic host
cold_apohost_goods = intersect(Ocu_ApoHost_Cold_BP$term, Dixon_BP$term)
cold_apo = Ocu_ApoHost_Cold_BP[Ocu_ApoHost_Cold_BP$term %in% cold_apohost_goods,]
Dixon_BP_set = Dixon_BP[Dixon_BP$term %in% cold_apohost_goods,]

# Combine them
cold_apohost_BP_data = merge(cold_apo, Dixon_BP_set, by="term")

cold_apohost_BP_plot = ggplot(cold_apohost_BP_data, aes(delta.rank.x, delta.rank.y)) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0,15)) +
  geom_hex() +
  labs( x = "Aposymbiotic Host",
        y = "Dixon") +
  geom_vline(xintercept = 0, linetype = 2, alpha = 0.75) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75) +
  geom_smooth(method=lm, lty = 2, color="black", se =FALSE) +
  xlim(c(-6000,6000)) +
  theme_bw()
cold_apohost_BP_plot

# look for most correlated categories
# can id delta ranks using this plotly:
cold_apohost_plotly = ggplot(cold_apohost_BP_data, aes(x = delta.rank.x, y = delta.rank.y)) +
  geom_point()
ggplotly(cold_apohost_plotly)

# and then filter data to figure out who is who:
cold_apohost_BP_neat_up = cold_apohost_BP_data %>%
  filter(delta.rank.x > 3000 & delta.rank.y > 1000)

cold_apohost_BP_neat_down = cold_apohost_BP_data %>%
  filter(delta.rank.x < -3000) %>%
  filter(delta.rank.y < -1000)

## merge all plots
all_BP = ggarrange(heat_symhost_BP_plot, heat_apohost_BP_plot, 
          cold_symhost_BP_plot, cold_apohost_BP_plot,
          ncol = 2, nrow = 2)
all_BP
ggsave(all_BP, file = "Dixon_Comparisons_Host_BP.pdf", width=6, height=4, units=c("in"), useDingbats=FALSE)
