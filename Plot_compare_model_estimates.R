
##Specifying the colors
red    = "#E27D60"
blue   = "#085DCB"
orange = "#E8A87C"

###############################################################################
##Comparing baseline model (HV_M2C) with model excluding summer minimum counts (HV_M1C)
###############################################################################

modA<-HV_M2C

png(paste0(pathR,"EffPlot1_phi1_M1M2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)
MCMCplot(object = HV_M1C, 
         object2 = modA, 
         params = 'phi1', 
         col="black",
         col2=blue,
         offset = 0.1,xlim=c(0.4,1.0),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate phi1")
dev.off()

png(paste0(pathR,"EffPlot1_f_M1M2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)
MCMCplot(object = HV_M1C, 
         object2 = modA, 
         params = 'f', 
         col2=blue,
         offset = 0.1,xlim=c(0.4,1.0),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate f")
dev.off()

png(paste0(pathR,"EffPlot1_p1_M1M2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)
MCMCplot(object = HV_M1C, 
         object2 = modA, 
         params = 'p1', 
         col2=blue,
         offset = 0.1,xlim=c(0.4,1),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate p1")
dev.off()

png(paste0(pathR,"EffPlot1_p2_M1M2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)

MCMCplot(object = HV_M1C,
         object2 = modA, 
         params = 'p2', 
         offset = 0.1,xlim=c(0,1),
         col2=blue,
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate p2")
dev.off()

png(paste0(pathR,"EffPlot1_Xtot_M1M2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)

MCMCplot(object = HV_M1C,
         object2 = modA, 
         params = 'Xtot', 
         col2=blue,
         offset = 0.1,
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate Xtot")
dev.off()


png(paste0(pathR,"EffPlot1_lambda_M1M2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)

MCMCplot(object = HV_M1C,
         object2 = modA, 
         params = 'lambda', 
         col2=blue,
         offset = 0.1,
         labels=c(
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate lambda")
dev.off()



###############################################################################
##Comparing baseline model (HV_M2C) with annual model (HV_M2CB23phi2)
###############################################################################
modB<-HV_M2C
modA<-HV_M2CB23phi2

png(paste0(pathR,"EffPlot1_phi1_M2M2phi2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)
MCMCplot(object = modB, 
         object2 = modA, 
         params = 'phi1',
         col=blue,
         col2=red, 
         offset = 0.1,xlim=c(0.4,1.0),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate phi1")
dev.off()

png(paste0(pathR,"EffPlot1_p1_M2M2phi2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)
MCMCplot(object = modB, 
         object2 = modA, 
         params = 'p1', 
         col=blue,
         col2=red,
         offset = 0.1,xlim=c(0.4,1),
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate p1")
dev.off()

png(paste0(pathR,"EffPlot1_p2_M2M2phi2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)

MCMCplot(object = modB,
         object2 = modA, 
         params = 'p2', 
         offset = 0.1,xlim=c(0,1),
         col=blue,
         col2=red,
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate p2")
dev.off()

png(paste0(pathR,"EffPlot1_Xtot_M2M2phi2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)

MCMCplot(object = modB,
         object2 = modA, 
         params = 'Xtot', 
         offset = 0.1,
         col=blue,
         col2=red,
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate Xtot")
dev.off()


png(paste0(pathR,"EffPlot1_lambda_M2M2phi2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)

MCMCplot(object = modB,
         object2 = modA, 
         params = 'lambda', 
         col=blue,
         col2=red,
         offset = 0.1,
         labels=c(
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate lambda")
dev.off()

png(paste0(pathR,"EffPlot1_f_M2M2phi2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)

MCMCplot(object = modB,
         object2 = modA, 
         params = 'f', 
         col=blue,
         col2=red,
         offset = 0.1,
         labels=c("2005",
           "2006","2007","2008","2009","2010","2011",
           "2012","2013","2014","2015","2016","2017",
           "2018","2019","2020","2021"),xlab="Parameter estimate f")
dev.off()

png(paste0(pathR,"EffPlot1_phi2_M2phi2.png"),width = 6, height = 6, units = "in", pointsize = 12,
    res=320,
    bg = "white", family = "", restoreConsole = TRUE)

MCMCplot(object = modA,
         params = 'phi2', xlim=c(0.4,1),
         offset = 0.1,col=red,
         labels=c(
           "2005/06",
           "2006/07","2007/08","2008/09","2009/10","2010/11","2011/12",
           "2012/13","2013/14","2014/15","2015/16","2016/17","2017/18",
           "2018/19","2019/20","2020/21"),xlab="Parameter estimate phi2")
dev.off()

###############################################################################
