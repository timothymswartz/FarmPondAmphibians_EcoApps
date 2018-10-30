library(unmarked)
library(AICcmodavg)
library(ggplot2)
library(ggpubr)
#https://sites.google.com/site/asrworkshop/home/schedule/r-occupancy-1
#https://cran.r-project.org/web/packages/AICcmodavg/AICcmodavg.pdf
#https://rdrr.io/cran/AICcmodavg/#vignettes
#https://rpubs.com/adamsmith_fws/CSNWRocc
setwd("~/Dropbox/Fall_2018/manuscripts/Chapter 2_Farm Pond Amphibians_Ecol Apps/Code and Data for Deposit")
data=read.csv("SSOM_1617_51_FullDataSet_6_7_18.csv")
options(scipen=999)  # turn off scientific notation like 1e+06


#All four frog species are analyzed separately using a single species, single-season occupancy model for each year of the study (2016 + 2017).
#Species names are abbreviated as follows:
  #HYLA: Hyla chrysoscelis/versicolor (Cope's/Eastern gray treefrog)
  #ACBL: Acris blanchardi (Blanchard's cricket frog)
  #LIBL: Lithobates blairi (Plains leopard frog)
  #LICA: Lithobates catesbeianus (American bullfrog)
  #In all cases, the year being analyzed is denoted as either "16" or "17", and the number of sites included is denoted with "51"

#HYLA16_51_ssom--------------------------------------------------------
hyla16y<-data[,6:9]
n<-nrow(data)
#site level individual covariates
hyla16.site<-data[,55:64]
#put everything together in unmarked data frame
#note that covariate can come from separate files
#BUT IF THERE ARE MULTIPLE INPUTS YOU HAVE TO USE AN "R list"
hyla16 <- unmarkedFrameOccu(y = hyla16y, siteCovs = data[55:64],obsCovs=list(Julian=data[35:38],Temp=data[39:42],Time=data[51:54]))
summary(hyla16) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
hylad61=occu(~1~1,hyla16)
hylad62=occu(~Time~1,hyla16)
hylad63=occu(~Temp~1,hyla16)
hylad64=occu(~Julian~1,hyla16)
hylad65=occu(~Julian+Time~1,hyla16)
hylad66=occu(~Temp*Time~1,hyla16)
hylad67=occu(~Temp+Time~1,hyla16)
hylad68=occu(~Temp+Julian~1,hyla16)
hylad69=occu(~Temp*Julian~1,hyla16)
#Goodness of Fit Test Bootstrap
#(gof_boot <- mb.gof.test(hylad61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
hyla16.Det.Cand.models=list(hylad61,hylad62,hylad63,hylad64,hylad66,hylad67,hylad68,hylad69)
hyla16.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = hyla16.Det.Cand.models,modnames = hyla16.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(hylad61,'det') # for detection
backTransform(hylad63,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
hyla16OccuNull=occu(~Temp~1,hyla16)
#Abiotic
hyla16OccuSL=occu(~Temp~SL16,hyla16)
hyla16OccuLogArea=occu(~Temp~LogArea,hyla16)
hyla16OccuCattle=occu(~Temp~Cattle16,hyla16)
hyla16.Abiotic.Occu.Cand.models=list(hyla16OccuSL,hyla16OccuLogArea,hyla16OccuCattle,hyla16OccuNull)
hyla16.Abiotic.Occu.Modnames=c("SL","LogArea","Cattle","Null")
print(aictab(cand.set = hyla16.Abiotic.Occu.Cand.models,modnames = hyla16.Abiotic.Occu.Modnames,second.ord = TRUE), digits = 4)
#Biotic
hyla16OccuEM=occu(~Temp~EM16,hyla16)
hyla16OccuAQ=occu(~Temp~AQ16,hyla16)
hyla16OccuFL=occu(~Temp~FL16,hyla16)
hyla16OccuFish=occu(~Temp~Fish16,hyla16)
#Biotic Candidate Set
hyla16.Biotic.Occu.Cand.models=list(hyla16OccuEM,hyla16OccuAQ,hyla16OccuFL,hyla16OccuFish,hyla16OccuNull)
hyla16.Biotic.ccu.Modnames=c("EM","AQ","FL","Fish","Null")
print(aictab(cand.set = hyla16.Biotic.Occu.Cand.models,modnames = hyla16.Biotic.ccu.Modnames,second.ord = TRUE), digits = 4)
#_______________________________
#Stage 2 - Interactions
hyla16OccuEM.Area=occu(~Temp~EM16+LogArea,hyla16)
#FINAL CANDIDATE SET
hyla16.FINAL.Occu.Cand.models=list(hyla16OccuEM,hyla16OccuLogArea,hyla16OccuEM.Area,hyla16OccuNull)
hyla16.FINAL.Occu.Modnames=c("EM","Area","EM.Area","Null")
#FINAL AICc Table
print(aictab(cand.set = hyla16.FINAL.Occu.Cand.models,modnames = hyla16.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#MODEL AVERAGING _____________________________
#Model averaging for each variable   ----> must exclude INTERACTION TERMS
modavg(cand.set = hyla16.FINAL.Occu.Cand.models, modnames = hyla16.FINAL.Occu.Modnames, parm = "EM16",parm.type = "psi",exclude = list("EM16:LogArea"),conf.level = 0.85)
modavg(cand.set = hyla16.FINAL.Occu.Cand.models, modnames = hyla16.FINAL.Occu.Modnames, parm = "LogArea",parm.type = "psi",exclude = list("EM16:LogArea"),conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predEM <- data.frame(EM16 = seq(from = min(data$EM16),to = max(data$EM16), by = 5),LogArea = mean(data$LogArea))
dat.predArea <- data.frame(LogArea = seq(from = min(data$LogArea),to = max(data$LogArea), by = 0.1),EM16 = mean(data$EM16))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionsHyla16EM=modavgPred(cand.set = hyla16.FINAL.Occu.Cand.models, modnames = hyla16.FINAL.Occu.Modnames,newdata = dat.predEM, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsHyla16EM,'predictionsHyla16EM.csv')
predictionsHyla16LogArea=modavgPred(cand.set = hyla16.FINAL.Occu.Cand.models, modnames = hyla16.FINAL.Occu.Modnames,newdata = dat.predArea, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsHyla16LogArea,'predictionsHyla16LogArea.csv')
#PLOTTING IN GGPLOT_____________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT!!!
    Hyla16EMpreds=as.data.frame(predictionsHyla16EM)
    Hyla16LogAreapreds=as.data.frame(predictionsHyla16LogArea)
#Write CSVs
#    write.csv(Hyla16EMpreds, file ="Hyla16EMpreds.csv")
#    write.csv(Hyla16LogAreapreds, file ="Hyla16LogAreapreds.csv")
#...HYLA 16 plots -------------------------------
(hyla16em <- ggplot(Hyla16EMpreds, aes(x=seq(from = 0, to = 65, by = 5), y=mod.avg.pred))+
   geom_line(data=Hyla16EMpreds)+
   geom_ribbon(data=Hyla16EMpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Emergent Vegetation Cover (%)"))
(hyla16logarea <- ggplot(Hyla16LogAreapreds, aes(x=seq(from = 2.3, to = 4.3, by = 0.1), y=mod.avg.pred))+
   geom_line(data=Hyla16LogAreapreds)+
   geom_ribbon(data=Hyla16LogAreapreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Log Area"))
  
    
# ACBL16_51_ssom------------------------------------------------------
acbl16y<-data[,2:5]
n<-nrow(data)
#site level individual covariates
acbl16.site<-data[,55:64]
#put everything together in unmarked data frame
#note that covariate can come from separate files
#BUT IF THERE ARE MULTIPLE INPUTS YOU HAVE TO USE AN "R list"
acbl16 <- unmarkedFrameOccu(y = acbl16y, siteCovs = data[55:64],obsCovs=list(Julian=data[35:38],Temp=data[39:42],Time=data[51:54]))
summary(acbl16) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
acbld61=occu(~1~1,acbl16)
acbld62=occu(~Time~1,acbl16)
acbld63=occu(~Temp~1,acbl16)
acbld64=occu(~Julian~1,acbl16)
acbld65=occu(~Julian+Time~1,acbl16)
acbld66=occu(~Temp*Time~1,acbl16)
acbld67=occu(~Temp+Time~1,acbl16)
acbld68=occu(~Temp+Julian~1,acbl16)
acbld69=occu(~Temp*Julian~1,acbl16)
#Goodness of Fit Test Bootstrap
#(gof_boot <- mb.gof.test(acbl16.global, nsim = 500, plot.hist=TRUE))
acbl16.global=occu(~Julian+Time+Temp~EM16+AQ16+FL16+LogArea+Fish16+SL16,acbl16)
#AICc Model Selection Table
acbl16.Det.Cand.models=list(acbld61,acbld62,acbld63,acbld64,acbld66,acbld67,acbld68,acbld69)
acbl16.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = acbl16.Det.Cand.models,modnames = acbl16.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(acbld61,'det') # for detection
backTransform(acbld68,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
acbl16OccuNull=occu(~Temp+Julian~1,acbl16)
#Abiotic
acbl16OccuSL=occu(~Temp+Julian~SL16,acbl16)
acbl16OccuLogArea=occu(~Temp+Julian~LogArea,acbl16)
acbl16OccuCattle=occu(~Temp+Julian~Cattle16,acbl16)
acbl16.Abiotic.Occu.Cand.models=list(acbl16OccuSL,acbl16OccuLogArea,acbl16OccuCattle,acbl16OccuNull)
acbl16.Abiotic.Occu.Modnames=c("SL","LogArea","Cattle","Null")
print(aictab(cand.set = acbl16.Abiotic.Occu.Cand.models,modnames = acbl16.Abiotic.Occu.Modnames,second.ord = TRUE), digits = 4)
#Biotic
acbl16OccuEM=occu(~Temp+Julian~EM16,acbl16)
acbl16OccuAQ=occu(~Temp+Julian~AQ16,acbl16)
acbl16OccuFL=occu(~Temp+Julian~FL16,acbl16)
acbl16OccuFish=occu(~Temp+Julian~Fish16,acbl16)
acbl16.Biotic.Occu.Cand.models=list(acbl16OccuEM,acbl16OccuAQ,acbl16OccuFL,acbl16OccuFish,acbl16OccuNull)
acbl16.Biotic.ccu.Modnames=c("EM","AQ","FL","Fish","Null")
print(aictab(cand.set = acbl16.Biotic.Occu.Cand.models,modnames = acbl16.Biotic.ccu.Modnames,second.ord = TRUE), digits = 4)
#_______________________________
#Stage 2 - Interactions
acbl16OccuEM.AQ=occu(~Temp+Julian~EM16+AQ16,acbl16)
#FINAL CANDIDATE SET
acbl16.FINAL.Occu.Cand.models=list(acbl16OccuEM,acbl16OccuAQ,acbl16OccuEM.AQ,acbl16OccuNull)
acbl16.FINAL.Occu.Modnames=c("EM","AQ","EM.AQ","Null")
#FINAL AICc Table
print(aictab(cand.set = acbl16.FINAL.Occu.Cand.models,modnames = acbl16.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#MODEL AVERAGING ______________________
#Model averaging for each variable   ----> must exclude INTERACTION TERMS
modavg(cand.set = acbl16.FINAL.Occu.Cand.models, modnames = acbl16.FINAL.Occu.Modnames, parm = "EM16",parm.type = "psi",exclude = list("EM16:AQ16"),conf.level = 0.85)
modavg(cand.set = acbl16.FINAL.Occu.Cand.models, modnames = acbl16.FINAL.Occu.Modnames, parm = "AQ16",parm.type = "psi",exclude = list("EM16:AQ16"),conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predEM <- data.frame(EM16 = seq(from = min(data$EM16),to = max(data$EM16), by = 5),AQ16 = mean(data$AQ16))
dat.predAQ <- data.frame(AQ16 = seq(from = min(data$AQ16),to = max(data$AQ16), by = 5),EM16 = mean(data$EM16))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionsacbl16EM=modavgPred(cand.set = acbl16.FINAL.Occu.Cand.models, modnames = acbl16.FINAL.Occu.Modnames,newdata = dat.predEM, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsacbl16EM,'predictionsacbl16EM.csv')
predictionsacbl16AQ=modavgPred(cand.set = acbl16.FINAL.Occu.Cand.models, modnames = acbl16.FINAL.Occu.Modnames,newdata = dat.predAQ, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsacbl16LogArea,'predictionsacbl16LogArea.csv')
summary(data$AQ16)
#PLOTTING IN GGPLOT________________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT!
acbl16EMpreds=as.data.frame(predictionsacbl16EM)
acbl16AQpreds=as.data.frame(predictionsacbl16AQ)
#write to csv
#write.csv(acbl16EMpreds, file ="acbl16EMpreds.csv")
#write.csv(acbl16AQpreds, file ="acbl16AQpreds.csv")

#...ACBL 16 plots -------------------------------
(acbl16em <- ggplot(acbl16EMpreds, aes(x=seq(from = 0, to = 65, by = 5), y=mod.avg.pred))+
    geom_line(data=acbl16EMpreds)+
    geom_ribbon(data=acbl16EMpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Emergent Vegetation Cover (%)"))
(acbl16aq <- ggplot(acbl16AQpreds, aes(x=seq(from = 0, to = 100, by = 5), y=mod.avg.pred))+
    geom_line(data=acbl16AQpreds)+
    geom_ribbon(data=acbl16AQpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Submerged Vegetation Cover (%)"))

#LIBL16_51_ssom ------------------------------------
libl16y<-data[,10:13]
n<-nrow(data)
#site level individual covariates
libl16.site<-data[,55:64]
#put everything together in unmarked data frame
#note that covariate can come from separate files
#BUT IF THERE ARE MULTIPLE INPUTS YOU HAVE TO USE AN "R list"
libl16 <- unmarkedFrameOccu(y = libl16y, siteCovs = data[55:64],obsCovs=list(Julian=data[35:38],Temp=data[39:42],Time=data[51:54]))
summary(libl16) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
libld61=occu(~1~1,libl16)
libld62=occu(~Time~1,libl16)
libld63=occu(~Temp~1,libl16)
libld64=occu(~Julian~1,libl16)
libld65=occu(~Julian+Time~1,libl16)
libld66=occu(~Temp*Time~1,libl16)
libld67=occu(~Temp+Time~1,libl16)
libld68=occu(~Temp+Julian~1,libl16)
libld69=occu(~Temp*Julian~1,libl16)
#Goodness of Fit Test Bootstrap

#(gof_boot <- mb.gof.test(libld61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
libl16.Det.Cand.models=list(libld61,libld62,libld63,libld64,libld66,libld67,libld68,libld69)
libl16.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = libl16.Det.Cand.models,modnames = libl16.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(libld61,'det') # for detection
backTransform(libld61,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
libl16OccuNull=occu(~1~1,libl16)
#Abiotic
libl16OccuSL=occu(~1~SL16,libl16)
libl16OccuLogArea=occu(~1~LogArea,libl16)
libl16OccuCattle=occu(~1~Cattle16,libl16)
libl16.Abiotic.Occu.Cand.models=list(libl16OccuSL,libl16OccuLogArea,libl16OccuCattle,libl16OccuNull)
libl16.Abiotic.Occu.Modnames=c("SL","LogArea","Cattle","Null")
print(aictab(cand.set = libl16.Abiotic.Occu.Cand.models,modnames = libl16.Abiotic.Occu.Modnames,second.ord = TRUE), digits = 4)
#Biotic
libl16OccuEM=occu(~1~EM16,libl16)
libl16OccuAQ=occu(~1~AQ16,libl16)
libl16OccuFL=occu(~1~FL16,libl16)
libl16OccuFish=occu(~1~Fish16,libl16)
libl16.Biotic.Occu.Cand.models=list(libl16OccuEM,libl16OccuAQ,libl16OccuFL,libl16OccuFish,libl16OccuNull)
libl16.Biotic.ccu.Modnames=c("EM","AQ","FL","Fish","Null")
print(aictab(cand.set = libl16.Biotic.Occu.Cand.models,modnames = libl16.Biotic.ccu.Modnames,second.ord = TRUE), digits = 4)
#_______________________________
#Stage 2 - Interactions
libl16OccuSL.Cattle=occu(~1~SL16+Cattle16,libl16)

#FINAL CANDIDATE SET
libl16.FINAL.Occu.Cand.models=list(libl16OccuSL,libl16OccuCattle,libl16OccuSL.Cattle,libl16OccuNull)
libl16.FINAL.Occu.Modnames=c("SL","Cattle","Sl.Cattle","Null")
#FINAL AICc Table
print(aictab(cand.set = libl16.FINAL.Occu.Cand.models,modnames=libl16.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#MODEL AVERAGING_____________________________
#Model averaging for each variable   ----> must exclude INTERACTION TERMS
modavg(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames, parm = "SL16",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames, parm = "Cattle16",parm.type = "psi",conf.level = 0.85)

#Create New data frames to store model averaged Predictions in 
dat.predSL <- data.frame(SL16 = seq(from = min(data$SL16),to = max(data$SL16), by = 0.01),Cattle16=mean(data$Cattle16))
dat.predCattle <- data.frame(Cattle16 = factor(c("0", "1")),SL16 = mean(data$SL16))
summary(data$SL16)
#Model-averaged predictions of psi across range of values
predictionslibl16SL=modavgPred(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames,newdata = dat.predSL, parm.type = "psi",conf.level = 0.85)
predictionslibl16Cattle=modavgPred(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames,newdata = dat.predCattle, parm.type = "psi",conf.level = 0.85)

#PLOTTING IN GGPLOT_____________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT!!!
libl16SLpreds=as.data.frame(predictionslibl16SL)
libl16Cattlepreds=as.data.frame(predictionslibl16Cattle)
#write to csv
#write.csv(libl16SLpreds, file ="libl16SLpreds.csv")
#write.csv(libl16Cattlepreds, file ="libl16Cattlepreds.csv")

#...LIBL 16 plots -------------------------------
(libl16sl <- ggplot(libl16SLpreds, aes(x=seq(from = 0, to = 0.35, by = 0.01), y=mod.avg.pred))+
    geom_line(data=libl16SLpreds)+
    geom_ribbon(data=libl16SLpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Pond Slope (rise:run)"))
(libl16Cattle <-ggplot(libl16Cattlepreds, aes(x=c(0,1), y=mod.avg.pred)) + 
    geom_bar(position=position_dodge(), stat="identity",alpha=0.3) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),width=.2,# Width of the error bars
                  position=position_dodge(.9))
    +theme_pubr()+
    scale_x_discrete(name="Cattle")+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1)))


#LICA16_51_ssom------------------------------------------------
lica16y<-data[,14:17]
n<-nrow(data)
#site level individual covariates
lica16.site<-data[,55:64]
#put everything together in unmarked data frame
#note that covariate can come from separate files
#BUT IF THERE ARE MULTIPLE INPUTS YOU HAVE TO USE AN "R list"
lica16 <- unmarkedFrameOccu(y = lica16y, siteCovs = data[55:64],obsCovs=list(Julian=data[35:38],Temp=data[39:42],Time=data[51:54]))
summary(lica16) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
licad61=occu(~1~1,lica16)
licad62=occu(~Time~1,lica16)
licad63=occu(~Temp~1,lica16)
licad64=occu(~Julian~1,lica16)
licad65=occu(~Julian+Time~1,lica16)
licad66=occu(~Temp*Time~1,lica16)
licad67=occu(~Temp+Time~1,lica16)
licad68=occu(~Temp+Julian~1,lica16)
licad69=occu(~Temp*Julian~1,lica16)
#Goodness of Fit Test Bootstrap

#(gof_boot <- mb.gof.test(licad61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
lica16.Det.Cand.models=list(licad61,licad62,licad63,licad64,licad66,licad67,licad68,licad69)
lica16.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = lica16.Det.Cand.models,modnames = lica16.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(licad61,'det') # for detection
backTransform(licad61,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
lica16OccuNull=occu(~1~1,lica16)
#Abiotic
lica16OccuSL=occu(~1~SL16,lica16)
lica16OccuLogArea=occu(~1~LogArea,lica16)
lica16OccuCattle=occu(~1~Cattle16,lica16)
lica16.Abiotic.Occu.Cand.models=list(lica16OccuSL,lica16OccuLogArea,lica16OccuCattle,lica16OccuNull)
lica16.Abiotic.Occu.Modnames=c("SL","LogArea","Cattle","Null")
print(aictab(cand.set = lica16.Abiotic.Occu.Cand.models,modnames = lica16.Abiotic.Occu.Modnames,second.ord = TRUE), digits = 4)
#Biotic
lica16OccuEM=occu(~1~EM16,lica16)
lica16OccuAQ=occu(~1~AQ16,lica16)
lica16OccuFL=occu(~1~FL16,lica16)
lica16OccuFish=occu(~1~Fish16,lica16)
lica16.Biotic.Occu.Cand.models=list(lica16OccuEM,lica16OccuAQ,lica16OccuFL,lica16OccuFish,lica16OccuNull)
lica16.Biotic.ccu.Modnames=c("EM","AQ","FL","Fish","Null")
print(aictab(cand.set = lica16.Biotic.Occu.Cand.models,modnames = lica16.Biotic.ccu.Modnames,second.ord = TRUE), digits = 4)
#_______________________________
#Stage 2 - Interactions
lica16OccuFL.Area=occu(~1~FL16+LogArea,lica16)
lica16OccuFL.by.Area=occu(~Temp+Julian~FL16*LogArea,lica16)
#FINAL CANDIDATE SET
lica16.FINAL.Occu.Cand.models=list(lica16OccuFL,lica16OccuLogArea,lica16OccuFL.Area,lica16OccuNull)
lica16.FINAL.Occu.Modnames=c("FL","Area","FL.Area","Null")
#FINAL AICc Table
print(aictab(cand.set = lica16.FINAL.Occu.Cand.models,modnames = lica16.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#MODEL AVERAGING __________________________
#Model averaging for each variable   ----> must exclude INTERACTION TERMS
modavg(cand.set = lica16.FINAL.Occu.Cand.models, modnames = lica16.FINAL.Occu.Modnames, parm = "FL16",parm.type = "psi",exclude = list("FL16:LogArea"),conf.level = 0.85)
modavg(cand.set = lica16.FINAL.Occu.Cand.models, modnames = lica16.FINAL.Occu.Modnames, parm = "LogArea",parm.type = "psi",exclude = list("FL16:LogArea"),conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predFL<- data.frame(FL16 = seq(from = min(data$FL16),to = max(data$FL16), by = 5),LogArea = mean(data$LogArea))
dat.predArea <- data.frame(LogArea = seq(from = min(data$LogArea),to = max(data$LogArea), by = 0.1),FL16 = mean(data$FL16))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionslica16FL=modavgPred(cand.set = lica16.FINAL.Occu.Cand.models, modnames = lica16.FINAL.Occu.Modnames,newdata = dat.predFL, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslica16EM,'predictionslica16EM.csv')
predictionslica16Area=modavgPred(cand.set = lica16.FINAL.Occu.Cand.models, modnames = lica16.FINAL.Occu.Modnames,newdata = dat.predArea, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslica16LogArea,'predictionslica16LogArea.csv')

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT!!!
lica16FLpreds=as.data.frame(predictionslica16FL)
lica16LogAreapreds=as.data.frame(predictionslica16Area)
summary(data$FL16)
#write to csv
#write.csv(lica16FLpreds, file ="lica16FLpreds.csv")
#write.csv(lica16LogAreapreds, file ="lica16LogAreapreds.csv")

#...LICA 16 plots -------------------------------
(lica16fl <- ggplot(lica16FLpreds, aes(x=seq(from = 0, to = 100, by = 5), y=mod.avg.pred))+
    geom_line(data=lica16FLpreds)+
    geom_ribbon(data=lica16FLpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Floating Vegetation Cover (%)"))
(lica16logarea <- ggplot(lica16LogAreapreds, aes(x=seq(from = 2.3, to = 4.3, by = 0.1), y=mod.avg.pred))+
    geom_line(data=lica16LogAreapreds)+
    geom_ribbon(data=lica16LogAreapreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Log Area"))


#HYLA17_51_SSOM---------------------------------------

hyla17y<-data[,18:21]
n<-nrow(data)
#site level individual covariates
hyla17.site<-data[,61:70]
#put everything together in unmarked data frame
#note that covariate can come from separate files
#BUT IF THERE ARE MULTIPLE INPUTS YOU HAVE TO USE AN "R list"
hyla17 <- unmarkedFrameOccu(y = hyla17y, siteCovs = data[61:70],obsCovs=list(Julian=data[43:46],Temp=data[47:50],Time=data[51:54]))
summary(hyla17) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
hylad71=occu(~1~1,hyla17)
hylad72=occu(~Time~1,hyla17)
hylad73=occu(~Temp~1,hyla17)
hylad74=occu(~Julian~1,hyla17)
hylad75=occu(~Julian+Time~1,hyla17)
hylad76=occu(~Temp*Time~1,hyla17)
hylad77=occu(~Temp+Time~1,hyla17)
hylad78=occu(~Temp+Julian~1,hyla17)
hylad79=occu(~Temp*Julian~1,hyla17)
#Goodness of Fit Test Bootstrap

#(gof_boot <- mb.gof.test(hylad61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
hyla17.Det.Cand.models=list(hylad61,hylad62,hylad63,hylad64,hylad66,hylad67,hylad68,hylad69)
hyla17.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = hyla17.Det.Cand.models,modnames = hyla17.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(hylad71,'det') # for detection
backTransform(hylad71,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
hyla17OccuNull=occu(~1~1,hyla17)
#Abiotic
hyla17OccuSL=occu(~1~SL17,hyla17)
hyla17OccuLogArea=occu(~1~LogArea,hyla17)
hyla17OccupH=occu(~1~pH,hyla17)
hyla17OccuCattle=occu(~1~Cattle17,hyla17)
hyla17.Abiotic.Occu.Cand.models=list(hyla17OccuSL,hyla17OccuLogArea,hyla17OccupH,hyla17OccuCattle,hyla17OccuNull)
hyla17.Abiotic.Occu.Modnames=c("SL","LogArea","pH","Cattle","Null")
print(aictab(cand.set = hyla17.Abiotic.Occu.Cand.models,modnames = hyla17.Abiotic.Occu.Modnames,second.ord = TRUE), digits = 4)
#Biotic
hyla17OccuEM=occu(~1~EM17,hyla17)
hyla17OccuAQ=occu(~1~AQ17,hyla17)
hyla17OccuFL=occu(~1~FL17,hyla17)
hyla17OccuFish=occu(~1~Fish17,hyla17)
hyla17.Biotic.Occu.Cand.models=list(hyla17OccuEM,hyla17OccuAQ,hyla17OccuFL,hyla17OccuFish,hyla17OccuNull)
hyla17.Biotic.ccu.Modnames=c("EM","AQ","FL","Fish","Null")
print(aictab(cand.set = hyla17.Biotic.Occu.Cand.models,modnames = hyla17.Biotic.ccu.Modnames,second.ord = TRUE), digits = 4)
#_______________________________
#Stage 2 - Interactions
hyla17OccuEM.pH=occu(~1~EM17+pH,hyla17)
#FINAL CANDIDATE SET
hyla17.FINAL.Occu.Cand.models=list(hyla17OccuEM,hyla17OccupH,hyla17OccuEM.pH,hyla17OccuNull)
hyla17.FINAL.Occu.Modnames=c("EM","pH","EM.pH","Null")
#FINAL AICc Table
print(aictab(cand.set = hyla17.FINAL.Occu.Cand.models,modnames = hyla17.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#MODEL AVERAGING __________________________
#Model averaging for each variable   ----> must exclude INTERACTION TERMS
modavg(cand.set = hyla17.FINAL.Occu.Cand.models, modnames = hyla17.FINAL.Occu.Modnames, parm = "EM17",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = hyla17.FINAL.Occu.Cand.models, modnames = hyla17.FINAL.Occu.Modnames, parm = "pH",parm.type = "psi",conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predEM<- data.frame(EM17 = seq(from = min(data$EM17),to = max(data$EM17), by = 5),pH = mean(data$pH))
dat.predpH <- data.frame(pH = seq(from = min(data$pH),to = max(data$pH), by = 0.1),EM17 = mean(data$EM17))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionshyla17EM=modavgPred(cand.set = hyla17.FINAL.Occu.Cand.models, modnames = hyla17.FINAL.Occu.Modnames,newdata = dat.predEM, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionshyla17EM,'predictionshyla17EM.csv')
predictionshyla17pH=modavgPred(cand.set = hyla17.FINAL.Occu.Cand.models, modnames = hyla17.FINAL.Occu.Modnames,newdata = dat.predpH, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionshyla17LogArea,'predictionshyla17LogArea.csv')

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT!!!
hyla17EMpreds=as.data.frame(predictionshyla17EM)
hyla17pHpreds=as.data.frame(predictionshyla17pH)
summary(data$EM17)
#write to csv
#write.csv(hyla17EMpreds, file ="hyla17EMpreds.csv")
#write.csv(hyla17pHpreds, file ="hyla17pHpreds.csv")

#...HYLA 17 plots -------------------------------
(hyla17em <- ggplot(hyla17EMpreds, aes(x=seq(from = 0, to = 70, by = 5), y=mod.avg.pred))+
   geom_line(data=hyla17EMpreds)+
   geom_ribbon(data=hyla17EMpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Emergent Vegetation Cover (%)"))
(hyla17pH <- ggplot(hyla17pHpreds, aes(x=seq(from = 6.94, to = 10.95, by = 0.1), y=mod.avg.pred))+
    geom_line(data=hyla17pHpreds)+
    geom_ribbon(data=hyla17pHpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="pH"))
summary(data$pH)
#ACBL17_51_SSOM---------------------------------------

acbl17y<-data[,22:25]
n<-nrow(data)
#site level individual covariates
acbl17.site<-data[,61:70]
#put everything together in unmarked data frame
#note that covariate can come from separate files
#BUT IF THERE ARE MULTIPLE INPUTS YOU HAVE TO USE AN "R list"
acbl17 <- unmarkedFrameOccu(y = acbl17y, siteCovs = data[61:70],obsCovs=list(Julian=data[43:46],Temp=data[47:50],Time=data[51:54]))
summary(acbl17) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
acbld71=occu(~1~1,acbl17)
acbld72=occu(~Time~1,acbl17)
acbld73=occu(~Temp~1,acbl17)
acbld74=occu(~Julian~1,acbl17)
acbld75=occu(~Julian+Time~1,acbl17)
acbld76=occu(~Temp*Time~1,acbl17)
acbld77=occu(~Temp+Time~1,acbl17)
acbld78=occu(~Temp+Julian~1,acbl17)
acbld79=occu(~Temp*Julian~1,acbl17)
#Goodness of Fit Test Bootstrap

#(gof_boot <- mb.gof.test(acbld61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
acbl17.Det.Cand.models=list(acbld71,acbld72,acbld73,acbld74,acbld76,acbld77,acbld78,acbld79)
acbl17.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = acbl17.Det.Cand.models,modnames = acbl17.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(acbld71,'det') # for detection
backTransform(acbld74,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
acbl17OccuNull=occu(~Julian~1,acbl17)
#Abiotic
acbl17OccuSL=occu(~Julian~SL17,acbl17)
acbl17OccuLogArea=occu(~Julian~LogArea,acbl17)
acbl17OccupH=occu(~Julian~pH,acbl17)
acbl17OccuCattle=occu(~Julian~Cattle17,acbl17)
acbl17.Abiotic.Occu.Cand.models=list(acbl17OccuSL,acbl17OccuLogArea,acbl17OccupH,acbl17OccuCattle,acbl17OccuNull)
acbl17.Abiotic.Occu.Modnames=c("SL","LogArea","pH","Cattle","Null")
print(aictab(cand.set = acbl17.Abiotic.Occu.Cand.models,modnames = acbl17.Abiotic.Occu.Modnames,second.ord = TRUE), digits = 4)
#Biotic
acbl17OccuEM=occu(~Julian~EM17,acbl17)
acbl17OccuAQ=occu(~Julian~AQ17,acbl17)
acbl17OccuFL=occu(~Julian~FL17,acbl17)
acbl17OccuFish=occu(~Julian~Fish17,acbl17)
acbl17.Biotic.Occu.Cand.models=list(acbl17OccuEM,acbl17OccuAQ,acbl17OccuFL,acbl17OccuFish,acbl17OccuNull)
acbl17.Biotic.ccu.Modnames=c("EM","AQ","FL","Fish","Null")
print(aictab(cand.set = acbl17.Biotic.Occu.Cand.models,modnames = acbl17.Biotic.ccu.Modnames,second.ord = TRUE), digits = 4)
#_______________________________
#Stage 2 - Interactions
acbl17OccuEM.pH=occu(~Julian~EM17+pH,acbl17)
acbl17OccuFL.by.Area=occu(~Temp+Julian~FL17*LogArea,acbl17)
#FINAL CANDIDATE SET
acbl17.FINAL.Occu.Cand.models=list(acbl17OccuEM,acbl17OccupH,acbl17OccuEM.pH,acbl17OccuNull)
acbl17.FINAL.Occu.Modnames=c("EM","pH","EM.pH","Null")
#FINAL AICc Table
print(aictab(cand.set = acbl17.FINAL.Occu.Cand.models,modnames = acbl17.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#MODEL AVERAGING __________________________
#Model averaging for each variable   ----> must exclude INTERACTION TERMS
modavg(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames, parm = "EM17",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames, parm = "pH",parm.type = "psi",conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predEM<- data.frame(EM17 = seq(from = min(data$EM17),to = max(data$EM17), by = 5),pH = mean(data$pH))
dat.predpH <- data.frame(pH = seq(from = min(data$pH),to = max(data$pH), by = 0.1),EM17 = mean(data$EM17))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionsacbl17EM=modavgPred(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames,newdata = dat.predEM, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsacbl17EM,'predictionsacbl17EM.csv')
predictionsacbl17pH=modavgPred(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames,newdata = dat.predpH, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsacbl17LogArea,'predictionsacbl17LogArea.csv')

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT!!!
acbl17EMpreds=as.data.frame(predictionsacbl17EM)
acbl17pHpreds=as.data.frame(predictionsacbl17pH)
summary(data$EM17)
#write to csv
#write.csv(acbl17EMpreds, file ="acbl17EMpreds.csv")
#write.csv(acbl17pHpreds, file ="acbl17pHpreds.csv")

#...ACBL 17 plots -------------------------------
(acbl17em <- ggplot(acbl17EMpreds, aes(x=seq(from = 0, to = 70, by = 5), y=mod.avg.pred))+
   geom_line(data=acbl17EMpreds)+
   geom_ribbon(data=acbl17EMpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Emergent Vegetation Cover (%)"))
(acbl17pH <- ggplot(acbl17pHpreds, aes(x=seq(from = 6.94, to = 10.95, by = 0.1), y=mod.avg.pred))+
    geom_line(data=acbl17pHpreds)+
    geom_ribbon(data=acbl17pHpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="pH"))


#LIBL17_51_SSOM---------------------------------------

libl17y<-data[,26:29]
n<-nrow(data)
#site level individual covariates
libl17.site<-data[,61:70]
#put everything together in unmarked data frame
#note that covariate can come from separate files
#BUT IF THERE ARE MULTIPLE INPUTS YOU HAVE TO USE AN "R list"
libl17 <- unmarkedFrameOccu(y = libl17y, siteCovs = data[61:70],obsCovs=list(Julian=data[43:46],Temp=data[47:50],Time=data[51:54]))
summary(libl17) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
libld71=occu(~1~1,libl17)
libld72=occu(~Time~1,libl17)
libld73=occu(~Temp~1,libl17)
libld74=occu(~Julian~1,libl17)
libld75=occu(~Julian+Time~1,libl17)
libld76=occu(~Temp*Time~1,libl17)
libld77=occu(~Temp+Time~1,libl17)
libld78=occu(~Temp+Julian~1,libl17)
libld79=occu(~Temp*Julian~1,libl17)
#Goodness of Fit Test Bootstrap

#(gof_boot <- mb.gof.test(libld61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
libl17.Det.Cand.models=list(libld71,libld72,libld73,libld74,libld76,libld77,libld78,libld79)
libl17.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = libl17.Det.Cand.models,modnames = libl17.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(libld71,'det') # for detection
backTransform(libld78,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
libl17OccuNull=occu(~Temp+Julian~1,libl17)
#Abiotic
libl17OccuSL=occu(~Temp+Julian~SL17,libl17)
libl17OccuLogArea=occu(~Temp+Julian~LogArea,libl17)
libl17OccupH=occu(~Temp+Julian~pH,libl17)
libl17OccuCattle=occu(~Temp+Julian~Cattle17,libl17)
libl17.Abiotic.Occu.Cand.models=list(libl17OccuSL,libl17OccuLogArea,libl17OccupH,libl17OccuCattle,libl17OccuNull)
libl17.Abiotic.Occu.Modnames=c("SL","LogArea","pH","Cattle","Null")
print(aictab(cand.set = libl17.Abiotic.Occu.Cand.models,modnames = libl17.Abiotic.Occu.Modnames,second.ord = TRUE), digits = 4)
#Biotic
libl17OccuEM=occu(~Temp+Julian~EM17,libl17)
libl17OccuAQ=occu(~Temp+Julian~AQ17,libl17)
libl17OccuFL=occu(~Temp+Julian~FL17,libl17)
libl17OccuFish=occu(~Temp+Julian~Fish17,libl17)
libl17.Biotic.Occu.Cand.models=list(libl17OccuEM,libl17OccuAQ,libl17OccuFL,libl17OccuFish,libl17OccuNull)
libl17.Biotic.ccu.Modnames=c("EM","AQ","FL","Fish","Null")
print(aictab(cand.set = libl17.Biotic.Occu.Cand.models,modnames = libl17.Biotic.ccu.Modnames,second.ord = TRUE), digits = 4)
#_______________________________
#Stage 2 - Interactions
libl17OccuAQ.Fish=occu(~Temp+Julian~AQ17+Fish17,libl17)
libl17OccuFL.by.Area=occu(~Temp+Julian~FL17*LogArea,libl17)
#FINAL CANDIDATE SET
libl17.FINAL.Occu.Cand.models=list(libl17OccuAQ,libl17OccuFish,libl17OccuAQ.Fish,libl17OccuNull)
libl17.FINAL.Occu.Modnames=c("AQ","Fish","AQ.Fish","Null")
#FINAL AICc Table
print(aictab(cand.set = libl17.FINAL.Occu.Cand.models,modnames = libl17.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#MODEL AVERAGING __________________________
#Model averaging for each variable   ----> must exclude INTERACTION TERMS
modavg(cand.set = libl17.FINAL.Occu.Cand.models, modnames = libl17.FINAL.Occu.Modnames, parm = "AQ17",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = libl17.FINAL.Occu.Cand.models, modnames = libl17.FINAL.Occu.Modnames, parm = "Fish17",parm.type = "psi",conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predAQ<- data.frame(AQ17 = seq(from = min(data$AQ17),to = max(data$AQ17), by = 5),Fish17 = mean(data$Fish17))
#dat.predFish <- data.frame(Fish17 = seq(from = min(data$Fish17),to = max(data$Fish17), by = 0.1),AQ17 = mean(data$AQ17))
dat.predFish<-data.frame(Fish17 = factor(c("0", "1")),AQ17 = mean(data$SL17))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionslibl17AQ=modavgPred(cand.set = libl17.FINAL.Occu.Cand.models, modnames = libl17.FINAL.Occu.Modnames,newdata = dat.predAQ, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslibl17EM,'predictionslibl17EM.csv')
predictionslibl17Fish=modavgPred(cand.set = libl17.FINAL.Occu.Cand.models, modnames = libl17.FINAL.Occu.Modnames,newdata = dat.predFish, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslibl17LogArea,'predictionslibl17LogArea.csv')

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT!!!
libl17AQpreds=as.data.frame(predictionslibl17AQ)
libl17Fishpreds=as.data.frame(predictionslibl17Fish)
summary(data$AQ17)
#write to csv
#write.csv(libl17AQpreds, file ="libl17AQpreds.csv")
#write.csv(libl17Fishpreds, file ="libl17Fishpreds.csv")

#...LIBL 17 plots -------------------------------
(libl17aq <- ggplot(libl17AQpreds, aes(x=seq(from = 0, to = 70, by = 5), y=mod.avg.pred))+
   geom_line(data=libl17AQpreds)+
   geom_ribbon(data=libl17AQpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Submerged Vegetation Cover (%)"))
(libl17Fish <-ggplot(libl17Fishpreds, aes(x=c(0,1), y=mod.avg.pred)) + 
    geom_bar(position=position_dodge(), stat="identity",alpha=0.3) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),width=.2,# Width of the error bars
                  position=position_dodge(.9))
    +theme_pubr()+
    scale_x_discrete(name="Fish")+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1)))
  
#LICA17_51_SSOM---------------------------------------

lica17y<-data[,30:33]
n<-nrow(data)
#site level individual covariates
lica17.site<-data[,61:70]
#put everything together in unmarked data frame
#note that covariate can come from separate files
#BUT IF THERE ARE MULTIPLE INPUTS YOU HAVE TO USE AN "R list"
lica17 <- unmarkedFrameOccu(y = lica17y, siteCovs = data[61:70],obsCovs=list(Julian=data[43:46],Temp=data[47:50],Time=data[51:54]))
summary(lica17) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
licad71=occu(~1~1,lica17)
licad72=occu(~Time~1,lica17)
licad73=occu(~Temp~1,lica17)
licad74=occu(~Julian~1,lica17)
licad75=occu(~Julian+Time~1,lica17)
licad76=occu(~Temp*Time~1,lica17)
licad77=occu(~Temp+Time~1,lica17)
licad78=occu(~Temp+Julian~1,lica17)
licad79=occu(~Temp*Julian~1,lica17)
#Goodness of Fit Test Bootstrap

#(gof_boot <- mb.gof.test(licad61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
lica17.Det.Cand.models=list(licad71,licad72,licad73,licad74,licad76,licad77,licad78,licad79)
lica17.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = lica17.Det.Cand.models,modnames = lica17.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(licad71,'det') # for detection
backTransform(licad71,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
lica17OccuNull=occu(~1~1,lica17)
#Abiotic
lica17OccuSL=occu(~1~SL17,lica17)
lica17OccuLogArea=occu(~1~LogArea,lica17)
lica17OccupH=occu(~1~pH,lica17)
lica17OccuCattle=occu(~1~Cattle17,lica17)
lica17.Abiotic.Occu.Cand.models=list(lica17OccuSL,lica17OccuLogArea,lica17OccupH,lica17OccuCattle,lica17OccuNull)
lica17.Abiotic.Occu.Modnames=c("SL","LogArea","pH","Cattle","Null")
print(aictab(cand.set = lica17.Abiotic.Occu.Cand.models,modnames = lica17.Abiotic.Occu.Modnames,second.ord = TRUE), digits = 4)
#Biotic
lica17OccuEM=occu(~1~EM17,lica17)
lica17OccuAQ=occu(~1~AQ17,lica17)
lica17OccuFL=occu(~1~FL17,lica17)
lica17OccuFish=occu(~1~Fish17,lica17)
lica17.Biotic.Occu.Cand.models=list(lica17OccuEM,lica17OccuAQ,lica17OccuFL,lica17OccuFish,lica17OccuNull)
lica17.Biotic.ccu.Modnames=c("EM","AQ","FL","Fish","Null")
print(aictab(cand.set = lica17.Biotic.Occu.Cand.models,modnames = lica17.Biotic.ccu.Modnames,second.ord = TRUE), digits = 4)
#_______________________________
#Stage 2 - Interactions
lica17OccuFL.Area=occu(~1~LogArea+FL17,lica17)
#lica17OccuFL.Area=occu(~1~FL17*LogArea,lica17)
#FINAL CANDIDATE SET
lica17.FINAL.Occu.Cand.models=list(lica17OccuFL,lica17OccuLogArea,lica17OccuFL.Area,lica17OccuNull)
lica17.FINAL.Occu.Modnames=c("FL","Area","FL.Area","Null")
#FINAL AICc Table
print(aictab(cand.set = lica17.FINAL.Occu.Cand.models,modnames = lica17.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#MODEL AVERAGING __________________________
#Model averaging for each variable   ----> must exclude INTERACTION TERMS
modavg(cand.set = lica17.FINAL.Occu.Cand.models, modnames = lica17.FINAL.Occu.Modnames, parm = "FL17",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = lica17.FINAL.Occu.Cand.models, modnames = lica17.FINAL.Occu.Modnames, parm = "LogArea",parm.type = "psi",conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predFL<- data.frame(FL17 = seq(from = min(data$FL17),to = max(data$FL17), by = 5),LogArea = mean(data$LogArea))
dat.predLogArea <- data.frame(LogArea = seq(from = min(data$LogArea),to = max(data$LogArea), by = 0.1),FL17 = mean(data$FL17))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionslica17Area=modavgPred(cand.set = lica17.FINAL.Occu.Cand.models, modnames = lica17.FINAL.Occu.Modnames,newdata = dat.predLogArea, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslica17EM,'predictionslica17EM.csv')
predictionslica17FL=modavgPred(cand.set = lica17.FINAL.Occu.Cand.models, modnames = lica17.FINAL.Occu.Modnames,newdata = dat.predFL, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslica17LogArea,'predictionslica17LogArea.csv')

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT!!!
lica17Areapreds=as.data.frame(predictionslica17Area)
lica17FLpreds=as.data.frame(predictionslica17FL)
summary(data$LogArea)
summary(data$FL17)
#write to csv
#write.csv(lica17Areapreds, file ="lica17Areapreds.csv")
#write.csv(lica17FLpreds, file ="lica17FLpreds.csv")

#...LICA 17 plots -------------------------------
(lica17fl <- ggplot(lica17FLpreds, aes(x=seq(from = 0, to = 95, by = 5), y=mod.avg.pred))+
   geom_line(data=lica17FLpreds)+
   geom_ribbon(data=lica17FLpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Floating Vegetation Cover (%)"))
(lica17area <- ggplot(lica17Areapreds, aes(x=seq(from = 2.3, to = 4.3, by = 0.1), y=mod.avg.pred))+
    geom_line(data=lica17Areapreds)+
    geom_ribbon(data=lica17Areapreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Log Area"))


#Final Figures--------------------------------------
ggarrange(acbl16em,acbl16aq,acbl17em,acbl17pH,
          hyla16em,hyla16logarea,hyla17em,hyla17pH,
          libl16sl,libl16Cattle,libl17aq,libl17Fish,
          lica16fl,lica16logarea,lica17fl,lica17area,ncol = 4, nrow = 4)
#FACET WRAP FIGURES-------------------------------------------------------
plot.multi=read.csv("ResultsforMultipanel.csv")
ggplot(plot.multi, aes(x = Value, y = Psi, fill = SpeciesVariable)) +
  geom_line() + 
  geom_ribbon(data=plot.multi,aes(ymin=LCL,ymax=UCL),fill="grey50",alpha=0.3,show.legend = NA)+
  facet_wrap(~ SpeciesVariable, ncol = 4,scales="free_x",
             #            labeller=as_labeller(c(ACBL16AQ=("ACBL (16) Submerged"),
             #                                   ACBL16EM="ACBL (16) Emergent (%)",
             #                                   ACBL17EM="ACBL (17) Emergent (%)", 
             #                                   ACBL17FL="ACBL (17) Floating (%)", 
             #                                   ACBL17PH="ACBL (17) pH",
             #                                   HYLA16Area="HYLA (16) Area (m^2)",
             #                                   HYLA16EM="HYLA (16) Emergent (%)",
             #                                   HYLA16Fish="HYLA (16) Fish",
             #                                   HYLA17AREA=" HYLA (17) Area (m^2)",
             #                                   HYLA17EM=" HYLA (17) Emergent (%)",
             #                                   HYLA17PH="HYLA (17) pH",
             #                                   HYLA17TREES="HYLA (17) Trees",
             #                                   LIBL16Cattle="LIBL (16) Cattle",
             #                                   LIBL16EM="LIBL (16) Emergent (%)",
             #                                   LIBL16SL="LIBL (16) Slope",
             #                                   LIBL17AQ="LIBL (17) Submerged (%)",
             #                                   LICA16AREA="LICA (16) Area (m^2)",
             #                                   LICA16FL="LICA (16) Floating (%)",
             #                                   LICA16TREES="LICA (16) Trees",
             #                                   LICA17AREA="LICA (17) Area (m^2)",
             #                                   LICA17FL="LICA (17) Floating (%)")),
             strip.position="top") +
  scale_x_continuous(name = "Parameter Value") + 
  scale_y_continuous(name = "Occupancy Probability") +
  theme(legend.position = "bottom") + 
  #guides(fill = F) +
  theme_bw()+ylab(NULL)+
  theme(panel.spacing = unit(1, "lines"))+
  #theme(strip.background = element_rect(fill="white", colour="black",size=.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))+
  theme(strip.text.x = element_text(size = 9.5, colour = "black"))+
  theme(axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black"))+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
