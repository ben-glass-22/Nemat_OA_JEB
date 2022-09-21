#Statistical analysis and plots for Nematostella ocean acidification experiment

#Load packages
library(ggplot2)
library(ggpubr)
library(plotrix)
library(tidyverse)
library(tidyr)
library(dplyr)
library(car)
library(emmeans)
library(Rmisc)
library(oce)
library(seacarb)
library(lubridate)
library(mgcv)
library(MuMIn)

#Carbonate chemistry
carbchemdata <- read.csv(file.choose(), header = TRUE)
carbchemdata$Date <- as.POSIXct(carbchemdata$Date, format = "%Y-%m-%d", origin = "1970-01-01")
carbchemdata$Treatment <- as.factor(carbchemdata$Treatment)

carbchemdata_SC <- carb(flag=8, var1=carbchemdata$pH., var2=carbchemdata$ALK_seacarb, S=carbchemdata$Salinity, T=carbchemdata$Temp)
carbchemdata_SC$Date <- carbchemdata$Date
carbchemdata_SC$Treatment <- carbchemdata$Treatment

carbchemsummary <- carbchemdata_SC %>%
  group_by(Treatment) %>%
  summarize_at(vars(pH, pCO2insitu, HCO3, CO3, T, DIC, S, ALK, OmegaAragonite, OmegaCalcite),
               funs(mean, std.error, sd, min, max), na.rm=T)

cc.pCO2.lm <- lm(pCO2insitu ~ Treatment, data = carbchemdata_SC, na.action = na.exclude)
Anova(cc.pCO2.lm)
#Anova shows significant differences in pCO2insitu between treatments

cc.pCO2.tukey <- emmeans(cc.pCO2.lm, list(pairwise ~ Treatment), adjust = "tukey")
cc.pCO2.tukey
#pCO2 for acidic treatment (1375) significantly higher than Ambient (494); df = 48, t-ratio = 6.295, p < .0001

cc.pH.lm <- lm(pH ~ Treatment, data = carbchemdata_SC, na.action = na.exclude)
Anova(cc.pH.lm)
#Anova shows significant differences in pH between treatments

cc.pH.tukey <- emmeans(cc.pH.lm, list(pairwise ~ Treatment), adjust = "tukey")
cc.pH.tukey
#pH for acidic treatment (7.402) significantly lower than Ambient (7.717); df = 102, t-ratio = -18.823, p < .0001

cc.alk.lm <- lm(ALK ~ Treatment, data = carbchemdata_SC, na.action = na.exclude)
Anova(cc.alk.lm)
#Anova shows significant differences in total alkalinity between treatments

cc.alk.tukey <- emmeans(cc.alk.lm, list(pairwise ~ Treatment), adjust = "tukey")
cc.alk.tukey
#TA for acidic treatment (0.001217) significantly higher than Ambient (0.000982); df = 48, t-ratio = 2.206, p = 0.032

cc.arag.lm <- lm(OmegaAragonite ~ Treatment, data = carbchemdata_SC, na.action = na.exclude)
Anova(cc.arag.lm)
#Anova shows significant differences in aragonite saturation state between treatments

cc.arag.tukey <- emmeans(cc.arag.lm, list(pairwise ~ Treatment), adjust = "tukey")
cc.arag.tukey
#aragonite saturation state for acidic treatment (0.213) significantly lower than Ambient (0.362); df = 48, t-ratio = -5.491, p-value < 0.001

ggplot(data = carbchemdata_SC, aes(x = Date, y = pH, color = Treatment, fill = Treatment)) +
  geom_smooth(size = 1, alpha = 0.4) +
  geom_point() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.89, 0.89), legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  ylim(7, 8) +
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) + 
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_x_datetime(name = "Week", breaks = "2 weeks", labels = c("Week 0","Week 1","Week 3", "Week 5", "Week 7", 
                                                                 "Week 9", "Week 11", "Week 13", "Week 15"))

ggsave("pH.png", width = 3, height = 3, units = c("in"))

#Male fecundity
mfecunditydata <- read.csv(file.choose(), header = TRUE)
mfecunditydata$Treatment <- factor(mfecunditydata$Treatment, levels = c("Ambient", "Acidic"))
mfecunditydata$Date <- factor(mfecunditydata$Date, levels = c("Week 2", "Week 4", "Week 7",
                                                              "Week 9", "Week 11", "Week 13"))

mfecunditymodel <- lm(Fecundity ~ Treatment + Date, data = mfecunditydata)
Anova(mfecunditymodel, type = "II")
#no significant differences detected

mfecundityposthoc <- emmeans(mfecunditymodel, list(pairwise~Treatment))
mfecundityposthoc

ggplot(data = mfecunditydata, aes_string(x = "Date", y = "Fecundity_div", 
                                         fill = "Treatment")) + 
  geom_point(position = position_jitterdodge(), aes_string(color = "Treatment")) +
  geom_boxplot(alpha = 0.40, outlier.shape = NA) + 
  labs(x = "Week") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.18, 0.80), legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9), axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  scale_fill_manual(values = c("lightblue", "tomato")) + 
  scale_color_manual(values = c("lightblue", "tomato")) +
  scale_y_continuous(expression(Male~Fecundity~(10^{5}~sperm~mL^{-1})))

ggsave("Male_fecundity.png", width = 3.5, height = 2.5, units = c("in"))

#Sperm mitochondrial membrane potential (MMP)
mmpdata <- read.csv(file.choose(), header = TRUE)
mmpdata$Treatment <- factor(mmpdata$Treatment, levels = c("Ambient", "Acidic"))
mmpdata$Date <- factor(mmpdata$Date, levels = c("Week 11", "Week 13"))

mmpmodel <- lm(MMP ~ Treatment*Date, data = mmpdata)
Anova(mmpmodel, type = "III")
mmpposthoc <- emmeans(mmpmodel, list(pairwise~Treatment*Date), adjust = "tukey", simple = "Treatment")
mmpposthoc
#Significant for both weeks; Week 11: df = 8, t-ratio = -16.144, p-value < .0001; Week 13: df = 8, t-ratio = -4.604, p-value = 0.0017

ggplot(data = mmpdata, aes_string(x = "Date", y = "MMP", fill = "Treatment")) + 
  geom_point(position = position_jitterdodge(), aes_string(color = "Treatment")) +
  geom_boxplot(alpha = 0.40, outlier.shape = NA) +
  labs(x = 'Week', y = 'Sperm with High MMP (%)') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none", axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_color_manual(values = c("lightblue", "tomato")) +
  ylim(71, 85)

ggsave("Sperm_MMP.png", width = 3, height = 2.5, units = c("in"))

#Fertilization
fertdata <- read.csv(file.choose(), header = TRUE)
fertdata$Date <- factor(fertdata$Date, levels = c("Week 2", "Week 4", "Week 9"))
fertdata$Treatment <- factor(fertdata$Treatment, levels = c("Ambient", "Acidic"))

fertdatalm <- lm(data = fertdata, formula = Per_Fert ~ Sperm_Con_Smaller + Treatment)
Anova(fertdatalm, type = "II")
#Both concentration (F value = 5.8903, p-value = 0.018) and treatment (F value = 4.5351, p-value = 0.037) are significant

fertlmposthoc <- emmeans(fertdatalm, list(pairwise~Treatment))
fertlmposthoc

simspermcon <- seq(1.2875, 16.0050, 0.01)
acidic <- rep("Acidic", 1472)
Ambient <- rep("Ambient", 1472)

acidicdf <- data.frame(Sperm_Con_Smaller = simspermcon, Treatment = acidic)
Ambientdf <- data.frame(Sperm_Con_Smaller = simspermcon, Treatment = Ambient)
df_for_predict <- rbind(Ambientdf, acidicdf)

predicted <- predict(fertdatalm, newdata = df_for_predict, se.fit = TRUE)

predictions_df <- data.frame(Sperm_Con_Smaller = df_for_predict$Sperm_Con_Smaller, 
                             Treatment = df_for_predict$Treatment, 
                             low = predicted$fit - 2*predicted$se.fit, 
                             fit = predicted$fit, 
                             high = predicted$fit + 2*predicted$se.fit)

ggplot(data = predictions_df, aes(x = Sperm_Con_Smaller, y = fit)) +
  geom_path(aes(color = Treatment), size = 1) + 
  geom_ribbon(aes(x = Sperm_Con_Smaller, ymin = low, ymax = high, fill = Treatment), alpha = 0.4) +
  geom_point(data = fertdata, aes(x = Sperm_Con_Smaller, y = Per_Fert, color = Treatment)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        legend.position = "none") +
  scale_x_continuous(expand = c(0.015, 0.015), 
                     expression(Sperm~Concentration~(10^{5}~sperm~mL^{-1}))) + 
  ylab("Fertilization (%)") +
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato"))

ggsave("Fertilization.png", width = 6.5, height = 3, units = c("in"))

#Female fecundity
ffecunditydata <- read.csv(file.choose(), header = TRUE)
ffecunditydata$Treatment <- factor(ffecunditydata$Treatment, levels = c("Ambient", "Acidic"))
ffecunditydata$Date <- factor(ffecunditydata$Date, levels = c("Week 2", "Week 4", "Week 7",
                                                              "Week 9", "Week 11", "Week 13"))

ffecunditymodel <- lm(Eggs ~ Treatment*Date, data = ffecunditydata)
Anova(ffecunditymodel, type = "III")
ffecundityposthoc <- emmeans(ffecunditymodel, list(pairwise~Treatment*Date), 
                             adjust = "tukey", simple = "Treatment")
ffecundityposthoc
#significant difference detected for Week 4 (df = 146, t-ratio = 3.535, p-value = 0.0005) between Ambient and acidic

ggplot(data = ffecunditydata, aes_string(x = "Date", y = "Eggs", fill = "Treatment")) + 
  geom_point(position = position_jitterdodge(), aes_string(color = "Treatment")) +
  geom_boxplot(alpha = 0.40, outlier.shape = NA) + 
  labs(x = 'Week', y = 'Female Fecundity (eggs per female)') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none",
        legend.title = element_text(size = 8), axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 0), 
        axis.title.y = element_text(size = 8)) + 
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato"))

ggsave("Female_fecundity.png", width = 3.5, height = 2.5, units = c("in"))

#Egg size
eggsizedata <- read.csv(file.choose(), header = TRUE)
eggsizedata$Treatment <- factor(eggsizedata$Treatment, levels = c("Ambient", "Acidic"))
eggsizedata$Date <- factor(eggsizedata$Date, levels = c("Week 2", "Week 4", "Week 7",
                                                              "Week 9", "Week 11", "Week 13"))

eggsizemodel <- lm(Egg_Diameter ~ Treatment*Date, data = eggsizedata)
Anova(eggsizemodel, type = "III")
eggsizeposthoc <- emmeans(eggsizemodel, list(pairwise~Treatment*Date), adjust = "tukey", simple = "Date")
eggsizeposthoc
#significant differences detected between treatments for week 4 (df = 3067, t-ratio = -3.139, p-value = 0.0017), 
#week 9 (df = 3067, t-ratio = -4.233, p-value < .0001), and 11 (df = 3067, t-ratio = -7.056, p-value < .0001)

ggplot(data = eggsizedata, aes_string(x = "Date", y = "Egg_Diameter", fill = "Treatment")) +
  geom_point(position = position_jitterdodge(), aes_string(color = "Treatment")) +
  geom_boxplot(alpha = 0.40, outlier.shape = NA) +
  labs(x = 'Week', y = 'Egg Diameter (mm)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none",
        legend.title = element_text(size = 8), axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato"))

ggsave("Egg_size.png", width = 3.5, height = 2.5, units = c("in"))

#Egg number and size by bundle
eggnumsizedata <- read.csv(file.choose(), header = TRUE)
eggnumsizedata$Treatment <- factor(eggnumsizedata$Treatment, levels = c("Ambient", "Acidic"))
eggnumsizedata$Date <- factor(eggnumsizedata$Date, levels = c("Week 7", "Week 9", 
                                                              "Week 11", "Week 13"))

eggnumsizemodel <- lm(Eggs ~ Avg_egg_size, data = eggnumsizedata)
Anova(eggnumsizemodel, type = "II")
summary(eggnumsizemodel)

ggplot(data = eggnumsizedata, aes(x = Avg_egg_size, y = Eggs, color = Treatment)) + 
  geom_point(shape = 19) + 
  stat_cor(aes(group = 1), p.accuracy = 0.001, label.x = 0.185, label.y = 60,
           show.legend = FALSE, size = 3) + 
  geom_smooth(aes(x = Avg_egg_size, y = Eggs), method = "lm", se = TRUE, inherit.aes = FALSE,
              na.rm = TRUE, fullrange = TRUE, color = "black", fill = "lightgray", size = 1, alpha = 0.4) + 
  xlab("Average Egg Diameter (mm)") + 
  ylab("Eggs in Bundle") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.88, 0.92), legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato"))

ggsave("Egg_number_and_size.png", width = 3, height = 5, units = c("in"))

#Larval respiration
larvrespdata <- read.csv(file.choose(), header = TRUE)
larvrespdata$Treatment <- factor(larvrespdata$Treatment, levels = c("Ambient", "Acidic"))
larvrespdata$Date <- factor(larvrespdata$Date, levels = c("Week 11", "Week 13"))

larvrespmodel <- lm(Resp_rate ~ Treatment + Date, data = larvrespdata)
Anova(larvrespmodel, type = "II")
#No significant differences detected

ggplot(data = larvrespdata, aes_string(x = "Date", y = "Resp_rate", fill = "Treatment")) +
  geom_point(position = position_jitterdodge(), aes_string(color = "Treatment")) +
  geom_boxplot(alpha = 0.40, outlier.shape = NA) +
  labs(x = 'Week') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none",
        legend.title = element_text(size = 6), axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8, angle = 90, hjust = 0.5), 
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8), legend.key.size = unit(.15, "in")) + 
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_y_continuous(expression(Larval~Respiraiton~Rate~(nmol~O[2]~min^{-1}~larva^{-1})))

ggsave("Larval_respiration.png", width = 3, height = 3, units = c("in"))

#Planula development
planuladata <- read.csv(file.choose(), header = TRUE)
planuladata$Treatment <- factor(planuladata$Treatment, levels = c("Ambient", "Acidic"))
planuladata$Date <- factor(planuladata$Date, levels = c("Week 4", "Week 9"))

planulamodel <- lm(Prct_planula ~ Treatment + Date, data = planuladata)
Anova(planulamodel, type = "II")
#Only date (F value = 94.6403, p-value < .0001) is significant, not treatment (F value = 2.1288, p-value = 0.1515)
planulameans <- emmeans(planulamodel, list(pairwise~Date), adjust = "tukey")
planulameans

ggplot(data = planuladata, aes_string(x = "Date", y = "Prct_planula", fill = "Treatment")) +
  geom_point(position = position_jitterdodge(), aes_string(color = "Treatment")) +
  geom_boxplot(alpha = 0.40, outlier.shape = NA) +
  labs(x = 'Week', y = 'Larvae in Planula Stage at 3 DPF (%)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.2, 0.85), legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6), axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato"))

ggsave("Planula_development.png", width = 3, height = 3, units = c("in"))

#Larval settlement
settlementdata <- read.csv(file.choose(), header = TRUE)
settlementdata$Treatment <- factor(settlementdata$Treatment, levels = c("Ambient", "Acidic"))
settlementdata$Date <- factor(settlementdata$Date, levels = c("Week 4", "Week 9"))

settlementmodel <- lm(Prct_settled ~ Treatment + Date, data = settlementdata)
Anova(settlementmodel, type = "II")
#Only date (F value = 48.2134, p-value < .0001) is significant
settlemeans <- emmeans(settlementmodel, list(pairwise~Date), adjust = "tukey")
settlemeans

ggplot(data = settlementdata, aes_string(x = "Date", y = "Prct_settled", 
                                         fill = "Treatment")) +
  geom_point(position = position_jitterdodge(), aes_string(color = "Treatment")) +
  geom_boxplot(alpha = 0.40, outlier.shape = NA) +
  labs(x = 'Week', y = 'Larvae Settled at 7 DPF (%)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none",
        legend.title = element_text(size = 6), axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato"))

ggsave("Larval_settlement.png", width = 3, height = 3, units = c("in"))

#Larval heat tolerance
larvheatbinarydata <- read.csv(file.choose(), header = TRUE)
larvheatbinarydata$Treatment <- factor(larvheatbinarydata$Treatment, 
                                       levels = c("Ambient", "Acidic"))
larvheatmodel <- glm(Survival ~ Treatment + Temperature, data = larvheatbinarydata, 
                     family = "binomial")
Anova(larvheatmodel, type = "II")
summary(larvheatmodel)
#No significant differences between treatments (p-value = 0.1056)

larvheatdata <- read.csv(file.choose(), header = TRUE)
larvheatdata$Treatment <- factor(larvheatdata$Treatment, levels = c("Ambient", "Acidic"))

ggplot(data = larvheatbinarydata, aes_string(x = "Temperature", y = "Survival", 
                                       color = "Treatment", fill = "Treatment")) +
  stat_smooth(method = "glm", fullrange = TRUE, method.args = list(family = "binomial"), 
              alpha = 0.4, color = "black", fill = "lightgray", size = 1) +
  stat_cor(aes(group = 1), p.accuracy = 0.001, label.x = 40.75, label.y = .50, show.legend = FALSE,
           size = 3) + 
  geom_point(data = larvheatdata, aes_string(x = "Temperature", y = "Prop_Surv",
                                             color = "Treatment")) +
  labs(x = 'Temperature (°C)', y = 'Larvae Surviving (%)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none",
        legend.title = element_text(size = 6), axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_y_continuous(labels = function(x) paste(x * 100))

ggsave("Larval_heat_tolerance.png", width = 3, height = 3, units = c("in"))

#Adult respiration
adultrespdata <- read.csv(file.choose(), header = TRUE)
adultrespdata$Treatment <- factor(adultrespdata$Treatment, levels = c("Ambient", "Acidic"))
adultrespdata$Sex <- factor(adultrespdata$Sex, levels = c("Male", "Female"))

adultrespmodel <- lm(Resp_rate ~ Treatment + Sex, data = adultrespdata)
Anova(adultrespmodel, type = "II")
#No significant difference between treatments (p-value = 0.05112) or by sex (p-value = 0.21795)

ggplot(data = adultrespdata, aes_string(x = "Treatment", y = "Resp_rate", 
                                        fill = "Treatment")) +
  facet_wrap(~ Sex) +
  geom_point(position = position_jitterdodge(), aes_string(color = "Treatment")) +
  geom_boxplot(alpha = 0.40, outlier.shape = NA) +
  labs(x = 'Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none",
        legend.title = element_text(size = 6), axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8)) + 
  scale_fill_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato")) +
  scale_y_continuous(expression(Adult~Respiraiton~Rate~(nmol~O[2]~min^{-1}~anemone^{-1}))) +
  scale_color_manual(values = c("Ambient" = "lightblue", "Acidic" = "tomato"))

ggsave("Adult_respiration.png", width = 4, height = 4, units = c("in"))
