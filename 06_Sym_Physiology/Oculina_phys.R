# Analyze Oculina arbuscula Fv/Fm data

library(tidyverse)
library(Rmisc)
library(lme4)
library(car)
library(emmeans)

#### Read in & organize data ####
# Read in Fv/Fm data and separate by symbiotic and aposybiotic

pam = read.csv('~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/06_Sym_Physiology/Oculina_PAM.csv')

str(pam)
head(pam)
#pam$day = as.factor(pam$day)
pam$treatment = factor(pam$treatment, levels = c("control","cold","heat"))
pam$sym_state = as.factor(pam$sym_state)
pam$genet = as.factor(str_sub(pam$coral_id, start = 1, end = 1))

pam2 = pam %>%
  mutate(avgfvfm = rowMeans(select(., fvfm1, fvfm2, fvfm3))) %>%
  filter(complete.cases(avgfvfm))

pam.sym = pam2 %>%
  subset(sym_state == 'Sym') %>%
  filter(genet == "C" | genet == "D" | genet == "F" | genet == "I" | genet == "J" | genet == "M" | genet == "R")

pam.apo = pam2 %>%
  subset(sym_state == 'Apo')

# remove POLKA genotypes from dataset


# Read in temperature data
host_exp_temp = read.csv("~/Dropbox/BU/Host_Buffering/FinalGithub/Oculina_Host_Sym_GE/02_ExperimentalDesign/Temperature_data/Oculina_temperature.csv") %>%
  select(-X) %>%
  mutate(treatment = factor(treatment))
head(host_exp_temp)
str(host_exp_temp)
levels(host_exp_temp$treatment) = list(cold = "cold", control = "control", heat = "hot") 

# find average daily temperature and filter out data so we can merge with pam
day_temps = summarySE(host_exp_temp, measurevar = "temp_C", groupvars = c("day","treatment"))
day_temps_to_merge = day_temps %>%
  select(day, treatment, temp_C) %>%
  filter(day == "2"| day == "5" | day == "8" | day == "11" | day == "14")


# merge temp and pam
pam_temps = merge(pam.sym, day_temps_to_merge, by = c("day","treatment"))

#### Plot Data ####
# plot fv/fm data
pam.sym.summary = summarySE(data = pam_temps, measurevar = "avgfvfm", groupvars = c("day", "treatment", "temp_C"))

cols_treatment = c("control" = "#a6611a", "cold" = "#74c476", "heat" = "#fd8d3c")
#coef = 30
# point and whiskers plot
pam.sym.plot = ggplot(pam.sym.summary, aes(x = day, color = treatment))+
  theme_bw()+
  geom_point(aes(y = avgfvfm), size = 3, position = position_dodge(width = 0.3))+
  geom_errorbar(aes(x = day, ymax = avgfvfm+se, ymin = avgfvfm-se), width = .3, position = position_dodge(width=0.3))+
  #geom_point(aes(y = temp_C/coef), size = 3, alpha = 0.5) +
  scale_y_continuous(name = "Average Photochemical Efficiency (Fv/Fm)")+
  #scale_y_continuous(name = "Average Photochemical Efficiency (Fv/Fm)", sec.axis = sec_axis(~.*coef, name = "Temperature (°C)"))+
  scale_color_manual(values = cols_treatment)+
  scale_x_continuous(breaks = c(2,5,8,11,14)) +
  theme(legend.position = 'none')
pam.sym.plot

ggsave(pam.sym.plot, file = "sym.host.fvfm.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# raw fv/fm data and means + se
pam.sym.plot.2 <- ggplot(pam_temps,aes(x = day, y = avgfvfm))+
  theme_bw()+
  geom_jitter(aes(color = treatment, fill = treatment), 
              position=position_dodge(width=1), 
              alpha=0.2, pch = 21,
              color = "black") +
  geom_errorbar(data = pam.sym.summary, aes(x = day, ymax = avgfvfm+se, ymin = avgfvfm-se, color = treatment), width = .2, position = position_dodge(width=1)) +
  geom_point(data = pam.sym.summary, mapping = aes(x=day, y=avgfvfm, color = treatment, fill = treatment), size = 2.5, pch = 21, color = "black", position = position_dodge(width=1))+ 
  scale_fill_manual(name = "Treatment",
                    breaks = c("cold","control","heat"),
                    values = cols_treatment)+
  scale_color_manual(name = "Treatment",
                     breaks = c("cold","control","heat"),
                     values = cols_treatment)+
  xlab("Day")+
  ylab("Photochemical Efficiency (Fv/Fm)")+
  scale_x_continuous(breaks = c(2,5,8,11,14))
#geom_hline(yintercept=0, linetype='dotted', color = 'gray')+
#theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) 
pam.sym.plot.2  
ggsave(pam.sym.plot.2, file = "sym.host.fvfm_alldata.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)



# plot temperature data
temp.pam.plot = day_temps_to_merge %>%
  ggplot(aes(x = day, y = temp_C, color = treatment))+
  theme_bw() +
  geom_point(aes(color = treatment), size=2, alpha = 0.9)+
  geom_line() +
  scale_color_manual(name = "Treatment", values = cols_treatment)+
  xlab("Day") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(2,32,4)) +
  scale_x_continuous(breaks = c(2,5,8,11,14))
#theme(legend.position = 'none')
temp.pam.plot
ggsave(temp.pam.plot, file = "temps.fvfm.pdf", width=5, height=2, units=c("in"), useDingbats=FALSE)



#### Stats ####
head(pam_temps)
str(pam_temps)
pam_temps$day = as.factor(pam_temps$day)

lm.pam = lmer(avgfvfm ~ day*treatment + (1|genet), data = pam_temps, REML = TRUE)
summary(lm.pam)
Anova(lm.pam)

# Response: avgfvfm
# Chisq Df Pr(>Chisq)    
# day            75.074  4  1.922e-15 ***
# treatment     163.371  2  < 2.2e-16 ***
# day:treatment 151.534  8  < 2.2e-16 ***
  
emms<-emmeans(lm.pam, ~treatment|day) #, adjust="Bonferoni"
pairs(emms, interaction = "pairwise") %>% rbind(adjust="fdr")

# day treatment_pairwise estimate     SE  df t.ratio p.value
# 2   control - cold     -0.00167 0.0198 282  -0.084  0.9328
# 2   control - heat      0.01689 0.0198 282   0.854  0.4544
# 2   cold - heat         0.01855 0.0203 282   0.916  0.4508
# 5   control - cold      0.03647 0.0195 282   1.869  0.1342
# 5   control - heat      0.03078 0.0195 282   1.578  0.1928
# 5   cold - heat        -0.00568 0.0197 282  -0.288  0.8289
# 8   control - cold      0.05378 0.0195 282   2.757  0.0155
# 8   control - heat      0.02771 0.0195 282   1.421  0.2348
# 8   cold - heat        -0.02607 0.0197 282  -1.321  0.2560
# 11  control - cold      0.22160 0.0195 282  11.359  <.0001
# 11  control - heat      0.03163 0.0195 282   1.621  0.1928
# 11  cold - heat        -0.18997 0.0197 282  -9.624  <.0001
# 14  control - cold      0.23822 0.0195 282  12.211  <.0001
# 14  control - heat      0.08213 0.0195 282   4.210  0.0001
# 14  cold - heat        -0.15608 0.0197 282  -7.908  <.0001



  