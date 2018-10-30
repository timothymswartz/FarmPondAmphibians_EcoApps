#AT WORK
setwd("C:/Users/tuj71018/Dropbox/Fall_2018/manuscripts/Chapter 2_Farm Pond Amphibians_Ecol Apps/Code and Data for Deposit")
#AT HOME
setwd("~/Dropbox/Fall_2017/Analysis_Fall_2017//Age Comparison Analysis")

data=read.csv("Age_Comparison_Analysis_PondHabitatVars_withAge_51_1617.csv")

#Packages to load---------------------------------
library(AICcmodavg)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(multcomp)
library(sjstats)
library(lsr)

#Add reference lines to plot
#http://www.sthda.com/english/wiki/ggplot2-add-straight-lines-to-a-plot-horizontal-vertical-and-regression-lines

#ANOVA effects sizes: http://three-mode.leidenuniv.nl/mtl/mtl_smc3_mancova.html

#For each variable, the year of the measurement was taken is noted as either a "1617" to indicate that measurements from 2016 and 2017 were averaged for this comparison. pH, was only measured in 2017 so no year is listed.
#The following abbreviations are used throughout:
  #FL= Floating vegetation percent cover
  #EM= Emergent vegetation percent cover
  #AQ= submerged/aquatic vegetation percent cover
  #LogArea= log-transformed area of each pond (in meteres^2)
  #cattle= presence (1) or absence (0) of cattle around the pond

#Dependent Variables (Habitat)
sl1617=data$SL1617
em1617=data$EM1617
fl1617=data$FL1617
aq1617=data$AQ1617
ph=data$pH

#Independent Variables
Age=data$AgeCategory
Cattle1617=data$Cattle1617
treesedge=data$TreesEdge

#Emergent Vegetation Cover (%) -----------------------------------------------------

#one-way ANCOVA
EM.mod <- lm(data$EM1617 ~ data$AgeCategory+data$TreesEdge+data$Cattle1617)
summary(EM.mod)
(anEm=anova(EM.mod))
#write.csv(anEm,"EMANCOVA.csv")

#derive ANCOVA effects size as eta-squared
(ETA.sq.EM <- aov(EM1617 ~ Cattle1617+AgeCategory+TreesEdge,data=data))
etaSquared(ETA.sq.EM,type=2)

#Run post hoc Tukey's test to determine significant differences among age categories 
a1 <- aov(EM1617 ~ AgeCategory,data=data)
tukEm=glht(a1,linfct=mcp(AgeCategory="Tukey"))
summary(tukEm)
tuk.cld<- cld(tukEm)
par(mar=c(1,1,5,1))
plot(tuk.cld) #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)

#Check for significant interactions
EM.INT.mod <- lm(data$EM1617 ~ data$AgeCategory*data$TreesEdge+data$AgeCategory*data$Cattle1617+data$TreesEdge+data$Cattle1617+data$AgeCategory)
summary(EM.INT.mod)
(anEM.INT=anova(EM.INT.mod))

#Generate plots
Emsummary=group_by(data, AgeCategory) %>%
  summarise(
    count = n(),
    mean = mean(EM1617, na.rm = TRUE),
    # se=count/(sqrt(EM1617)),
    sd = sd(EM1617, na.rm = TRUE),
    median = median(EM1617, na.rm = TRUE),
    IQR = IQR(EM1617, na.rm = TRUE)
  )
#write.csv(Emsummary,"EmSummary.csv")
(Emplot=ggplot(Emsummary,aes(AgeCategory,mean))+
    geom_line(aes(group=1))+
    geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
    geom_point(data=data, aes(AgeCategory,EM1617, size="1"),color="darkgreen",position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,
                                                                                                            dodge.width = 0.4))+
    geom_line(aes(x=AgeCategory,y=mean)))
(Emnice= Emplot + labs(y="Emergent Vegetation Cover (%)", x = "Pond Age (years)") + theme_pubr()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +geom_point())

#Floating Vegetation Cover (%)--------------------------------

#one-way ANCOVA
FL.mod <- lm(data$FL1617 ~ data$AgeCategory+data$TreesEdge+data$Cattle1617)
summary(FL.mod)
(anfl=anova(FL.mod))
#write.csv(anfl,"FLANCOVA.csv")

#ANCOVA effects partial eta squared 
(ETA.sq.FL <- aov(FL1617 ~ Cattle1617+AgeCategory+TreesEdge,data=data))
etaSquared(ETA.sq.FL,type=2)

#Run post hoc Tukey's test to determine significant differences among age categories 
a5 <- aov(FL1617 ~ AgeCategory,data=data)
tukFL=glht(a5,linfct=mcp(AgeCategory="Tukey"))
summary(tukFL)
tuk.cld<- cld(tukFL)
par(mar=c(1,1,5,1))
plot(tuk.cld)  #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)

#Check for significant interactions
FL.INT.mod <- lm(data$FL1617 ~ data$AgeCategory*data$TreesEdge+data$AgeCategory*data$Cattle1617+data$TreesEdge+data$Cattle1617+data$AgeCategory)
summary(FL.INT.mod)
(anFL.INT=anova(FL.INT.mod))

# Generate Plots
FLsummary=group_by(data, AgeCategory) %>%
  summarise(
    count = n(),
    mean = mean(FL1617, na.rm = TRUE),
    # se=count/(sqrt(AQ1617)),
    sd = sd(FL1617, na.rm = TRUE),
    median = median(FL1617, na.rm = TRUE),
    IQR = IQR(FL1617, na.rm = TRUE)
  )
#write.csv(FLsummary,"FLsummary.csv")
(FLplot=ggplot(FLsummary,aes(AgeCategory,mean))+
    geom_line(aes(group=1))+
    geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
    geom_point(data=data, aes(AgeCategory,FL1617, size="1"),color="darkgreen",position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,
                                                                                                            dodge.width = 0.4))+
    geom_line(aes(x=AgeCategory,y=mean)))
(FLnice=  FLplot + labs(y="Floating Vegetation Cover", x = "Pond Age") + theme_pubr() +geom_point()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))

# Submerged Vegetation Cover (%)------------------------------------------------

#one-way ANCOVA
SUB.mod <- lm(data$AQ1617 ~ data$AgeCategory+data$TreesEdge+data$Cattle1617)
summary(SUB.mod)
(ansub=anova(SUB.mod))
write.csv(ansub,"SUBANCOVA.csv")

#ANCOVA effects partial eta squared 
(ETA.sq.SUB <- aov(AQ1617 ~ Cattle1617+AgeCategory+TreesEdge,data=data))
etaSquared(ETA.sq.SUB,type=2)

#Run post hoc Tukey's test to determine significant differences among age categories 
a3 <- aov(AQ1617 ~ AgeCategory,data=data)
tukSub=glht(a3,linfct=mcp(AgeCategory="Tukey"))
summary(tukSub)
tuk.cld<- cld(tukSub)
par(mar=c(1,1,5,1))
plot(tuk.cld)  #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)

#Check for significant interactions
SUB.INT.mod <- lm(data$AQ1617 ~ data$AgeCategory*data$TreesEdge+data$AgeCategory*data$Cattle1617+data$TreesEdge+data$Cattle1617+data$AgeCategory)
summary(SUB.INT.mod)
(anSUB.INT=anova(SUB.INT.mod))

#Generate plots for the data
Subsummary=group_by(data, AgeCategory) %>%
  summarise(
    count = n(),
    mean = mean(AQ1617, na.rm = TRUE),
    # se=count/(sqrt(AQ1617)),
    sd = sd(AQ1617, na.rm = TRUE),
    median = median(AQ1617, na.rm = TRUE),
    IQR = IQR(AQ1617, na.rm = TRUE)
  )
#write.csv(Subsummary,"Subsummary.csv")
(Subplot=ggplot(Subsummary,aes(AgeCategory,mean))+
    geom_line(aes(group=1))+
    geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
    geom_point(data=data, aes(AgeCategory,AQ1617, size="1"),color="darkgreen",position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,
                                                                                                            dodge.width = 0.4))+
    geom_line(aes(x=AgeCategory,y=mean)))
(Subnice= Subplot + labs(y="Submerged Vegetation Cover (%)", x = "Pond Age") + theme_pubr() +geom_point()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))

#pH ------------------------------------------------

#one-way ANCOVA
PH.mod <- lm(data$pH ~ data$AgeCategory+data$TreesEdge+data$Cattle1617)
summary(PH.mod)
(anpH=anova(PH.mod))
#write.csv(anpH,"PHANCOVA.csv")

#ANCOVA effects partial eta squared 
(ETA.sq.pH <- aov(pH ~ Cattle1617+AgeCategory+TreesEdge,data=data))
etaSquared(ETA.sq.pH,type=2)

#Run post hoc Tukey's test to determine significant differences among age categories 
a2 <- aov(pH ~ AgeCategory,data=data)
tukpH=glht(a2,linfct=mcp(AgeCategory="Tukey"))
summary(tukpH)
tuk.cld<- cld(tukpH)
par(mar=c(1,1,5,1))   
plot(tuk.cld)  #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)

#Check for significant interactions
PH.INT.mod <- lm(data$pH ~ data$AgeCategory*data$TreesEdge+data$AgeCategory*data$Cattle1617+data$TreesEdge+data$Cattle1617+data$AgeCategory)
summary(PH.INT.mod)
(anPH.INT=anova(PH.INT.mod))

#Generate plots
pHsummary=group_by(data, AgeCategory) %>%
  summarise(
    count = n(),
    mean = mean(pH, na.rm = TRUE),
    # se=count/(sqrt(pH)),
    sd = sd(pH, na.rm = TRUE),
    median = median(pH, na.rm = TRUE),
    IQR = IQR(pH, na.rm = TRUE)
  )
#write.csv(pHsummary,"pHsummary.csv")
(pHplot=ggplot(pHsummary,aes(AgeCategory,mean))+
    geom_line(aes(group=1))+
    geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
    geom_point(data=data, aes(AgeCategory,pH, size="1"),color="darkgreen",position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,
                                                                                                        dodge.width = 0.4))+
    geom_line(aes(x=AgeCategory,y=mean)))
(pHnice= pHplot + labs(y="pH", x = "Pond Age") + theme_pubr() +geom_point() +theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))

#Slope (rise:run)--------------------------------------------

#one-way ANCOVA
SL.mod <- lm(data$SL1617 ~ data$AgeCategory+data$TreesEdge+data$Cattle1617)
summary(SL.mod)
(ansl=anova(SL.mod))
#write.csv(ansl,"SLANCOVA.csv")
#ANCOVA effects partial eta squared 
(ETA.sq.SL <- aov(SL1617 ~ Cattle1617+AgeCategory+TreesEdge,data=data))
etaSquared(ETA.sq.SL,type=2)

#Run post hoc Tukey's test to determine significant differences among age categories 
a4 <- aov(SL1617 ~ AgeCategory,data=data)
tukSL=glht(a4,linfct=mcp(AgeCategory="Tukey"))
summary(tukSL)
tuk.cld<- cld(tukSL)
par(mar=c(1,1,5,1))
plot(tuk.cld)  #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)

#Check for significant interactions
SL.INT.mod <- lm(data$SL1617 ~ data$AgeCategory*data$TreesEdge+data$AgeCategory*data$Cattle1617+data$TreesEdge+data$Cattle1617+data$AgeCategory)
summary(SL.INT.mod)
(anSL.INT=anova(SL.INT.mod))

#Generate plots
   SLsummary=group_by(data, AgeCategory) %>%
     summarise(
       count = n(),
       mean = mean(SL1617, na.rm = TRUE),
       # se=count/(sqrt(AQ1617)),
       sd = sd(SL1617, na.rm = TRUE),
       median = median(SL1617, na.rm = TRUE),
       IQR = IQR(SL1617, na.rm = TRUE)
     )
   write.csv(SLsummary,"SLsummary.csv")
   (SLplot=ggplot(SLsummary,aes(AgeCategory,mean))+
       geom_line(aes(group=1))+
       geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
       geom_point(data=data, aes(AgeCategory,SL1617, size="1"),color="darkgreen",position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,
                                                                                                               dodge.width = 0.4))+
       geom_line(aes(x=AgeCategory,y=mean)))
   (SLnice= SLplot + labs(y="Slope", x = "Pond Age") + theme_pubr() +geom_point()+theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))

# Generate final multi-paneled plot of all 5 variables -------------------
grid.arrange(Emnice, FLnice, Subnice, pHnice,SLnice, nrow = 3)
