# This script analyzes the temperature data from the Oculina arbuscula holobiont experiment and the Breviolum psygmophilum culture experiment
# as well as environmental data from the coral collection site (Radio Island NC). 
# Author: Hannah E Aichelman

#### Load in and format data ####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Rmisc)


host_exp = read.csv("02_ExperimentalDesign/Temperature_data/Oculina_temperature.csv") %>%
  select(-X) %>%
  mutate(treatment = factor(treatment))
head(host_exp)
str(host_exp)

culture_exp = read.csv("02_ExperimentalDesign/Temperature_data/All_Temps.csv") %>%
  mutate(Incubator = factor(Incubator))
head(culture_exp)
str(culture_exp)

# work with date time format to subset date
culture_exp$DateTime<-strptime(culture_exp$DateTime, format="%m/%d/%y %H:%M")
culture_exp$datetime_ct <- as.POSIXct(culture_exp$DateTime, format="%Y-%m-%dT%H:%M:%S")

culture_exp$Date<-format(culture_exp$datetime_ct,"%D")
culture_exp$Date<-as.POSIXlt(culture_exp$Day, format="%m/%d/%y")
culture_exp$Day = as.factor(culture_exp$Day)

#### Make plots ####

#cols = c("control" = "#a6611a", "cold" = "#74c476", "hot" = "#fd8d3c")
cols_host = c("control" = "grey", "cold" = "#3182bd", "hot" = "#de2d26")

host.plot = host_exp %>%
  ggplot(aes(x = day, y = temp_C, color = treatment))+
  theme_bw() +
  geom_point(aes(color = treatment), size=2, alpha = 0.9)+
  #geom_line() +
  scale_color_manual(name = "Treatment", values = cols_host)+
  xlab("Day") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(2,32,4)) +
  scale_x_continuous(breaks = seq(1,15,1))+
  geom_hline(yintercept = 32, linetype = "solid", color = "#de2d26", size = 1) +
  geom_hline(yintercept = 6, linetype = "solid", color = "#3182bd", size = 1) +
  ggtitle("A. Holobiont Experiment") +
  theme(legend.position = c(.1, .2), legend.background = element_rect(color = "black"))
host.plot
ggsave(host.plot, filename = "/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/plots/Host_Temp_Plot.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)

# now culture plot
cols_sym = c("control" = "#a6611a", "cool" = "#74c476", "heat" = "#fd8d3c")

culture.plot = culture_exp %>%
  ggplot(aes(x = Date, y = TempC, color = Incubator))+
  theme_bw() +
  geom_point(aes(color = Incubator), size=2, alpha = 0.9, shape = 21)+
  scale_color_manual(name = "Treatment", values = cols_sym)+
  xlab("Day") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(2,32,4)) +
  scale_x_discrete(labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))+
  geom_hline(yintercept = 32, linetype = "dashed", color = "#fd8d3c", size = 1) +
  geom_hline(yintercept = 6, linetype = "dashed", color = "#74c476", size = 1) +
  #geom_hline(yintercept = 18, linetype = "dashed", color = "#a6611a", size = 1) +
  ggtitle("B. Culture Experiment") +
  theme(legend.position = c(.1, .2), legend.background = element_rect(color = "black"))
culture.plot
ggsave(culture.plot, filename = "Culture_Temp_Plot.pdf", width=5, height=4, units=c("in"), useDingbats=FALSE)


# combine host and symbiont plots
figs.combined = ggarrange(host.plot, culture.plot, ncol = 2, nrow = 1)
ggsave(figs.combined, filename = "Combined_Temp_Plot.pdf", width=8, height=4, units=c("in"), useDingbats=FALSE)


#### Ribbon plot experiments together ####
culture_exp_tomerge = culture_exp %>%
  select(Date, Incubator, TempC) %>%
  mutate(Date = factor(Date)) %>%
  mutate(Day = case_when(
    Date == "10/30/20" ~ "1",
    Date == "10/31/20" ~ "2",
    Date == "11/01/20" ~ "3",
    Date == "11/02/20" ~ "4",
    Date == "11/03/20" ~ "5",
    Date == "11/04/20" ~ "6",
    Date == "11/05/20" ~ "7",
    Date == "11/06/20" ~ "8",
    Date == "11/07/20" ~ "9",
    Date == "11/08/20" ~ "10",
    Date == "11/09/20" ~ "11",
    Date == "11/10/20" ~ "12",
    Date == "11/11/20" ~ "13",
    Date == "11/12/20" ~ "14",
    Date == "11/13/20" ~ "15")) %>%
  mutate(Day = as.numeric(Day)) %>%
  group_by(Day, Incubator) %>% 
  dplyr::summarise(across(where(is.numeric), list(mean = mean, sd = sd, min = min, max = max), na.rm = TRUE)) %>%
  mutate(Location = "Culture") %>%
  mutate(Location = factor(Location)) %>%
  dplyr::rename("treatment"="Incubator") %>%
  mutate(Day = factor(Day)) %>%
  mutate(treatment = fct_recode(treatment,
                                "control" = "control",
                                "cold" = "cool",
                                "heat" = "heat"))

str(culture_exp_tomerge)

host_exp_tomerge = host_exp %>%
  select(day, treatment, temp_C) %>%
  mutate(day = factor(day)) %>%
  dplyr::rename("Day"="day") %>%
  dplyr::rename("TempC" = "temp_C") %>%
  group_by(Day, treatment) %>% 
  dplyr::summarise(across(where(is.numeric), list(mean = mean, sd = sd, min = min, max = max), na.rm = TRUE)) %>%
  mutate(Location = "InHost") %>%
  mutate(Location = factor(Location)) %>%
  mutate(treatment = fct_recode(treatment,
                                "control" = "control",
                                "cold" = "cool",
                                "heat" = "heat"))

str(host_exp_tomerge)

temps_combined = rbind(culture_exp_tomerge, host_exp_tomerge)
temps_combined$treat_type = as.factor(paste(temps_combined$Location, temps_combined$treatment, sep = "_"))
str(temps_combined)

# now plot
combined.colors = c("Culture_cold" = "#74c476","Culture_control" = "#a6611a","Culture_heat" = "#fd8d3c","InHost_cold"= "#3182bd", "InHost_control" = "grey", "InHost_hot" = "#de2d26")

ribbon.plot = ggplot(data = temps_combined, aes(x = Day, y = TempC_mean, color = treat_type, group = treat_type))+
  geom_ribbon(aes(ymin = TempC_mean - TempC_sd, ymax = TempC_mean + TempC_sd, fill = treat_type), alpha = 0.3, colour = NA)+
  geom_line()+
  theme_bw() +
  scale_color_manual(name = "Treatment", values = combined.colors) +
  xlab("Day") +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(2,32,4)) +
  scale_x_discrete(labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))+
  geom_hline(yintercept = 32, linetype = "dashed", color = "#fd8d3c", size = 1) +
  geom_hline(yintercept = 6, linetype = "dashed", color = "#74c476", size = 1) +
  ggtitle("B. Culture Experiment") +
  theme(legend.position = c(.1, .2), legend.background = element_rect(color = "black"))
culture.plot

#### Summarize temperatures ####
str(host_exp)
host_temps = summarySE(data = host_exp, measurevar = "temp_C", groupvars = c("treatment","day"))
host_temps_fullavg = summarySE(data = host_exp, measurevar = "temp_C", groupvars = c("treatment"))

# treatment  N   temp_C        sd         se       ci
# 1      cold 81 11.26914 3.4784207 0.38649119 0.769142
# 2   control 81 17.86914 0.6653274 0.07392526 0.147116
# 3       hot 84 24.48690 4.1919826 0.45738280 0.909716

max(host_exp$temp_C)
#[1] 31.2
min(host_exp$temp_C)
#[1] 6.3

# take a look at salinity
host_exp_sal = host_exp %>%
  drop_na(salinity)
host_sal = summarySE(data = host_exp_sal, measurevar = "salinity", groupvars = "treatment")
# .id   N salinity        sd         se         ci
# 1 <NA> 243 33.98724 0.2279123 0.01462058 0.02879984

# treatment  N salinity        sd         se         ci
# 1      cold 81 33.86667 0.1313393 0.01459325 0.02904149
# 2   control 81 33.94321 0.2246877 0.02496530 0.04968253
# 3       hot 81 34.15185 0.2127858 0.02364286 0.04705079

str(culture_exp)
culture_temps = summarySE(data = culture_exp, measurevar = "TempC", groupvars = c("Incubator","Date"))
culture_temps_fullavg = summarySE(data = culture_exp, measurevar = "TempC", groupvars = c("Incubator"))

max(culture_exp$TempC)
#[1] 31.14
min(culture_exp$TempC)
#[1] 5.4


#### Plot Radio Island Temps for context with treatments ####

## NOAA buoy data for Radio Island
library(lubridate)


boing<-read.delim("02_ExperimentalDesign/Temperature_Data/RadioIsland_Temps/RadioIslandTempNOAABuoy_2017full.txt", na.strings='999')
head(boing)
boing<-cbind("datetime"=paste(paste(boing$X.YY.yr, boing$MM.mo, boing$DD.dy, sep="-")," ",paste(boing$hh.hr,boing$mm.mn, sep=":"),sep=""),boing)
boing$datetime<-strptime(boing$datetime, format="%Y-%m-%d %H:%M")
boing<-boing[order(boing$datetime),]

boing$datetime = as.POSIXct(boing$datetime, format = "%Y-%m-%d %H:%M:%S")

mean(boing$WTMP.degC, na.rm=TRUE)
#[1] 18.78945
min <- min(boing$WTMP.degC, na.rm=TRUE)
#[1] 3
max <- max(boing$WTMP.degC, na.rm=TRUE)
#[1] 28.5


# plot buoy temp data
ritemp_plot = ggplot(data = boing, aes(x = datetime, y = WTMP.degC)) +
  geom_line() +
  theme_bw() +
  scale_y_continuous(name = "Temperature (°C)", breaks = seq(2,32,4)) +
  scale_x_datetime(breaks = seq(as.POSIXct("2017-01-01 00:00:00"),
                                as.POSIXct("2017-12-31 23:54:00 EST"),
                                "1 month"),
                   date_labels = "%b-%Y") +
  geom_hline(yintercept = 31.8, linetype = "solid", color = "#de2d26", size = .5, alpha = 0.8) +
  geom_hline(yintercept = 32.2, linetype = "dashed", color = "#fd8d3c", size = .5, alpha = 0.8) +
  geom_hline(yintercept = 17.8, linetype = "solid", color = "grey", size = .5, alpha = 0.8) +
  geom_hline(yintercept = 18.2, linetype = "dashed", color = "#a6611a", size = .5, alpha = 0.8) +
  geom_hline(yintercept = 5.8, linetype = "solid", color = "#3182bd", size = .5, alpha = 0.8) +
  geom_hline(yintercept = 6.2, linetype = "dashed", color = "#74c476", size = .5, alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.title.x = element_blank())

ritemp_plot 
ggsave(ritemp_plot, filename = "RadioIsland_Temp_Plot.pdf", width=4, height=3.5, units=c("in"), useDingbats=FALSE)
