library(grid)
library(lattice) #utilized to plot multiple plots on a single page

##Specifying the colors
red    = "#E27D60"
blue   = "#085DCB"


###############################################################################
##Comparing baseline model (HV_M2C) with model excluding summer minimum counts (HV_M1C)
###############################################################################

modA<-HV_M2C

png(paste0(pathR,"SuppPlot1_v2.png"),width = 7.8, height = 12, units = "in", pointsize = 12,
    res=1200,
    bg = "white", family = "", restoreConsole = TRUE)

par(mfrow=c(3,2))

MCMCplot(object = HV_M1C, 
         object2 = modA, 
         params = 'phi1', 
         col="black",
         col2=blue,
         offset = 0.1,xlim=c(0.4,1.0),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Summer survival of calves (phi1)")

MCMCplot(object = HV_M1C, 
         object2 = modA, 
         params = 'f', 
         col2=blue,
         offset = 0.1,xlim=c(0.4,1.0),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Early recruitment rate (f)")

MCMCplot(object = HV_M1C,
         object2 = modA, 
         params = 'lambda', 
         col2=blue,
         offset = 0.1,
         labels=c(
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Population growth rate (λ)")

MCMCplot(object = HV_M1C,
         object2 = modA, 
         params = 'Xtot', 
         col2=blue,
         offset = 0.1,
         labels=c("2005",
                  "2006","2007","2008","2009","2010","2011",
                  "2012","2013","2014","2015","2016","2017",
                  "2018","2019","2020","2021"),xlab="Total population size after harvest (Xtot)")

MCMCplot(object = HV_M1C, 
         object2 = modA, 
         params = 'p1', 
         col2=blue,
         offset = 0.1,xlim=c(0.4,1),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Survey coverage in summer (p1)")

MCMCplot(object = HV_M1C,
         object2 = modA, 
         params = 'p2', 
         offset = 0.1,xlim=c(0,1),
         col2=blue,
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Survey coverage in fall (p2)")


dev.off()


###############################################################################
##Comparing baseline model (HV_M2C) with annual model (HV_M2CB23phi2)
###############################################################################
modB<-HV_M2C
modA<-HV_M2CB23phi2

png(paste0(pathR,"SuppPlot2_v2.png"),width = 7.8, height = 12, units = "in", pointsize = 12,
    res=1200,
    bg = "white", family = "", restoreConsole = TRUE)

par(mfrow=c(3,2))

MCMCplot(object = modB, 
         object2 = modA, 
         params = 'phi1',
         col=blue,
         col2=red, 
         offset = 0.1,xlim=c(0.4,1.0),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Summer survival of calves (phi1)")

MCMCplot(object = modA,
         params = 'phi2', xlim=c(0.4,1),
         offset = 0.1,col=red,
         labels=c(
           "2005/06",
           "2006/07","2007/08","2008/09","2009/10","2010/11","2011/12",
           "2012/13","2013/14","2014/15","2015/16","2016/17","2017/18",
           "2018/19","2019/20","2020/21"),xlab="Annual survival of pre-harvest calves to yearlings (phi2)")

MCMCplot(object = modB,
         object2 = modA, 
         params = 'lambda', 
         col=blue,
         col2=red,
         offset = 0.1,
         labels=c(
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Population growth rate (λ)")

MCMCplot(object = modB,
         object2 = modA, 
         params = 'Xtot', 
         offset = 0.1,
         col=blue,
         col2=red,
         labels=c("2005",
                  "2006","2007","2008","2009","2010","2011",
                  "2012","2013","2014","2015","2016","2017",
                  "2018","2019","2020","2021"),xlab="Total population size after harvest (Xtot)")

MCMCplot(object = modB, 
         object2 = modA, 
         params = 'p1', 
         col=blue,
         col2=red,
         offset = 0.1,xlim=c(0.4,1),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Survey coverage in summer (p1)")

MCMCplot(object = modB,
         object2 = modA, 
         params = 'p2', 
         offset = 0.1,xlim=c(0,1),
         col=blue,
         col2=red,
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Survey coverage in fall (p2)")


dev.off()

###############################################################################
