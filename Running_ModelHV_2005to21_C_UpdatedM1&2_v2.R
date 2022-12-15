rm(list=ls())
###################################################################################
library(jagsUI)
library(MCMCvis)

myscriptname="Running_ModelHV_2005to21_C_UpdatedM1&2_v2.R"

#####################
### SIMULATING SOME DATA - TO BE USED AS INITIAL VALUES HERE; 

source("Appendix 3 Simulated_Data.R")
Simulated_data <- SimPop(PHI3=0.9, mean_F=0.5, p1=0.5, p2=0.5, bias1=1, T=40)                                 # Run script to simulate data 

#####################
### READING AND PREPARING DATA FROM HARDANGERVIDDA
d <- read.csv("Data_Hardangervidda_2005_2021v2.csv", sep=";", header = T)
#2016: >300 reindeer killed by lightening was included as harvested
#2019: fallen stock (killed due to footrot symptoms) were included in harvested 

#Including data from 2005-2021
T_max <- nrow(d) #17


bugs.data.HVs <- list(SU = c(d[1:T_max,"SU"]),
                      J = c(d[1:T_max,"J"]),
                      H0f = d[1:T_max,"H0f"],
                      H0m = d[1:T_max,"H0m"],
                      H1f = d[1:T_max,"H1f"],
                      H1m = d[1:T_max,"H1m"],
                      Hadf = d[1:T_max,"Hadf"],
                      Hadm = d[1:T_max,"Hadm"],
                      T = T_max,
                      y_w=c(d[2:T_max, "N_v"],NA),
                      y_s=c(d[1:T_max, "N_s"]),
                      C0=d[1:T_max,"C0"],
                      Cm=d[1:T_max,"Cm"],
                      Cf=d[1:T_max,"Cf"],
                      B23=d[1:T_max,"B2"]+d[1:T_max,"B3."])

#vague prior for coverage of summer survey 2021 (p1[T])
bugs.data.HVs$alphapT<- 1 
bugs.data.HVs$betapT<- 1  


##################################
#Dataset including 2020 subset (counted and grouped similar to other years)
HV20subset_adf<-288 + 134 + 69
HV20subset_calf<-129 + 60 + 31
HV20subset_yearling<-48 + 13 + 7
HV20subset_adm<-89 + 28 + 6

sum20subset<-HV20subset_adf+HV20subset_calf+HV20subset_yearling+HV20subset_adm
sum20S<-sum(d[d$year==2020,c("C0", "Cf","Cm")])
#sum20subset/sum20S #16%


bugs.data.HV4s <- list(SU = c(d[1:T_max,"SU"]),
                      J = c(d[1:T_max,"J"]),
                      H0f = d[1:T_max,"H0f"],
                      H0m = d[1:T_max,"H0m"],
                      H1f = d[1:T_max,"H1f"],
                      H1m = d[1:T_max,"H1m"],
                      Hadf = d[1:T_max,"Hadf"],
                      Hadm = d[1:T_max,"Hadm"],
                      T = T_max,
                      y_w=c(d[2:T_max, "N_v"],NA),
                      y_s=c(d[1:T_max, "N_s"]),
                      C0=d[1:T_max,"C0"],
                      Cm=d[1:T_max,"Cm"],
                      Cf=d[1:T_max,"Cf"],
                      B23=d[1:T_max,"B2"]+d[1:T_max,"B3."],
                      C0sub20=HV20subset_calf,
                      C1sub20=HV20subset_yearling,
                      Cfsub20=HV20subset_adf,
                      B23sub20=HV20subset_adm
)

#vague prior for coverage of summer survey 2021 (p1[T])
bugs.data.HV4s$alphapT<- 1 
bugs.data.HV4s$betapT<- 1  


####### Folder for output results ###########
stopifnot(file.exists(myscriptname)) 
folder.out = "Output_Results/"
dir.create(folder.out)

path<-getwd()
pathR<-paste0(path,"/",folder.out)

###############################################################################
### PREPARING JAGS MODEL RUNS; 

# MCMC settings
niter <- 300000
nthin <- 3
nburn <- 50000
nchains <- 3

#From 2005 including model prediction pre-harvest 2022
length(Simulated_data$N[1, 20:37])

t=0
# Initial values
inits <- function() { list( phi3 = runif(1, 0.9, 0.99), 
                            mean.f = runif(1, 2, 3), 
                            mean.phi1 = runif(1, 2, 3), 
                            mean.phi2 = runif(1, 2, 3), 
                            sd_phi1 = runif(1, 0, 1),
                            sd_phi2 = runif(1, 0, 1),
                            sd_f = runif(1, 0, 1),
                            p2 = runif(T_max, 0.3, 0.9),
                            p1 = runif(T_max, 0.5, 0.99),
                            N0f = Simulated_data$N[1, (20+t):37]*350, 
                            N0m = Simulated_data$N[2, (20+t):37]*350, 
                            N1f = Simulated_data$N[3, (20+t):37]*350,
                            N1m = Simulated_data$N[4, (20+t):37]*350,
                            Nadf = Simulated_data$N[5, (20+t):37]*350,
                            Nadm = Simulated_data$N[6, (20+t):37]*350) }  


#######Run the models #########################################################
source("M1&2_variants_Hvidda_fall2021subset2020.R")

#Baseline model, not including the summer minimum counts 
#Alternative, if also including demographic counts for subset 2020: use data=bugs.data.HV4s
HV_M1C <- jags(data=bugs.data.HVs, inits=inits, 
               parameters.to.save=c(
                 "phi3", "phi1", "f", "p1", "p2", "Ntot","Xtot",
                 "Sjuv", "fert", "sd_f", "sd_phi1","Sjuv18","lambda"), 
               model.file="jagsMod_M1C_HV_fall2021p1.bug",n.chain=nchains, 
               n.iter=niter, n.burnin=nburn, parallel=TRUE)

HV_M_var <-HV_M1C
sumdat<-MCMCsummary(HV_M_var, round = 3)
sumdat
write.csv2(sumdat,file=paste0(pathR,"HvdatM1C_05_r1.csv"))


#Baseline model including summer counts of 2019 and 2021
#Alternative, if also including demographic counts for subset 2020: use data=bugs.data.HV4s
HV_M2C <- jags(data=bugs.data.HVs, inits=inits, 
               parameters.to.save=c( 
                 "phi3", "phi1", "f", "p1", "p2","Ntot","Xtot", 
                 "fert", "sd_f","Sjuv", "sd_phi1",
                 #"N0f", "N0m", "N1f", "N1m", "Nadf", "Nadm","Xadf","Xadm",
                 "Sjuv18","lambda",#"Mlambda18",
                 "DeadCalves[15]"), 
               model.file="jagsMod_M2C_HV_fall2021p1.bug",n.chain=nchains, 
               n.iter=niter, n.burnin=nburn, parallel=TRUE)

HV_M_var <-HV_M2C
sumdat<-MCMCsummary(HV_M_var, round = 3)
write.csv2(sumdat,file=paste0(pathR,"HvdatM2C_05_r1.csv"))


#Annual model, including summer counts of 2019 and 2021
#Including demographic counts for subset 2020
HV_M2CB23phi2 <- jags(data=bugs.data.HV4s, inits=inits, 
                      parameters.to.save=c("phi3", "phi2","phi1", "f", 
                                           "p1", "p2","ps", "Ntot", "Xtot", 
                                           "Sjuv", "fert", "sd_f", "sd_phi1", 
                                           "Sjuv18","Sjuvphi2_18","lambda"),
                      model.file="jagsMod_M2CB23phi2_HV_fall2021p1_E.bug",n.chain=nchains, 
                      n.iter=niter, n.burnin=nburn, parallel=TRUE)

HV_M_var2 <-HV_M2CB23phi2
sumdat<-MCMCsummary(HV_M_var2, round = 3)
write.csv2(sumdat,file=paste0(pathR,"HvdatM2CB23phi2_05_r1.csv"))

################################################################################
#Generate Figure 2A and B
################################################################################
blue = "#085DCB"

png(paste0(pathR,"Fig2A_EffPlot1_phi1_M2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)
MCMCplot(object = HV_M2C, 
         params = 'phi1', 
         col=blue,
         offset = 0.1,xlim=c(0.4,1.0),
         labels=c("2005",
                  "2006","2007","2008","2009","2010","2011",
                  "2012","2013","2014","2015","2016","2017",
                  "2018","2019","2020","2021"),xlab="Summer survival of calves (phi1)")
text(0.54,3.8,"2019",col=blue)

dev.off()



png(paste0(pathR,"Fig2B_EffPlot1_lambda_M2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)

MCMCplot(object = HV_M2C,
         params = 'lambda', 
         col=blue,
         offset = 0.1,
         labels=c(
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab=expression(paste("Population growth rate (",lambda,")")))
text(1.09,3.8,"2019",col=blue)
dev.off()


