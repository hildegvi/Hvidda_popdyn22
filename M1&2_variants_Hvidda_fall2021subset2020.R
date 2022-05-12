#Model variants - using informative or vague prior for p1[T]
#M1- model variants do not include summer minimum counts for last 3 years
#jagsMod_M1C_HV_fall2021p1.bug

#M2- model variants including summer minimum counts for 2019 and 2021
#jagsMod_M2C_HV_fall2021p1.bug - "Baseline model"
#jagsMod_M2CB23phi2_HV_fall2021p1.bug -"Annual model"

##############################################################
##jagsMod_M1C_HV_fall2021p1.bug ------------------------------
##############################################################
#summer minimum counts not included

sink("jagsMod_M1C_HV_fall2021p1.bug")
cat("
    model {
    
    
    ############################################################
    # Define the priors for the parameters
    ############################################################ 
    
    # Initial pre-harvest population sizes
    N0f[1] ~ dnorm(500, 0.0001)I(0,)     # calves
    N0m[1] ~ dnorm(500, 0.0001)I(0,)     # calves
    N1f[1]~ dnorm(300, 0.0001)I(0,)      # yearling females
    N1m[1]~ dnorm(300, 0.0001)I(0,)      # yearling males
    Nadf[1]~dnorm(500, 0.0001)I(0,)      # adult females
    Nadm[1]~dnorm(500, 0.0001)I(0,)      # adult females
    
    # Survival probabilites and productivity-estimates
    
    for(t in 1:(T-1)){
    p1[t]~ dunif(0,1)
    p2[t]~ dunif(0,1)
    }
    
    #prior for p1[T] - can be specified as informed or vague 
    for(t in T){
    p1[t]~ dbeta(alphapT,betapT)
    p2[t]~ dunif(0,1)
    }
    
    #prior for proportion of animals counted in subset 2020:
    ps~ dunif(0.1,0.2)

    ## Vital rates
    phi3 ~ dunif(0, 1)      
    mean.f ~ dunif(-5,5)
    mean.phi1 ~ dunif(-5,5)
    
    tau_f <- 1/(sd_f*sd_f)
    sd_f ~ dunif(0,2)
    
    tau_phi1 <- 1/(sd_phi1*sd_phi1)
    sd_phi1 ~ dunif(0,2)
    
    ##############################
    # DERIVED PARAMETERS
    #############################
    
    fert <- exp(mean.f)/(1+exp(mean.f))
    Sjuv <- exp(mean.phi1)/(1+exp(mean.phi1))  
    
    #############################
    # SYSTEM PROCESS
    #############################
    
    for (t in 1:(T)){  

    logit(f[t]) <- mean.f+eps_f[t]
    eps_f[t]~dnorm(0, tau_f)
    
    logit(phi1[t]) <- mean.phi1+eps_phi1[t]
    eps_phi1[t]~dnorm(0, tau_phi1)

    }


    for (t in 1:(T-2)){  
    N0f[t+1] ~ dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1])) 
    N0m[t+1] ~ dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1]))

    N1f[t+1] ~ dbin(phi3, max(50, round(N0f[t]-H0f[t])))
    N1m[t+1] ~ dbin(phi3, max(50, round(N0m[t]-H0m[t])))
    
    Nadf[t+1] ~ dbin(phi3, max(50, round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t])))
    Nadm[t+1] ~ dbin(phi3, max(50, round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t])))
    }
    
 
    #Average over 2001-2018 year 
    Sjuv18<-mean(phi1[2:(T-3)])

    #Use average over 6 last year for next year prediction
    fert6m<-mean(f[(T-6):T])
    Sjuv6m<-mean(phi1[(T-6):T])

    for (t in (T-1)){  
    N0f[t+1] ~dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1]))
    N0m[t+1] ~dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1]))

    N1f[t+1] ~ dbin(phi3, max(50, round(N0f[t]-H0f[t])))
    N1m[t+1] ~ dbin(phi3, max(50, round(N0m[t]-H0m[t])))
    
    Nadf[t+1] ~ dbin(phi3, max(50, round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t])))
    Nadm[t+1] ~ dbin(phi3, max(50, round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t])))
    }

    for (t in T){  
    N0f[t+1] ~dbin(Sjuv6m*fert6m/2, max(50, Nadf[t+1]))
    N0m[t+1] ~dbin(Sjuv6m*fert6m/2, max(50, Nadf[t+1]))

    N1f[t+1] ~ dbin(phi3, max(50, round(N0f[t]-H0f[t])))
    N1m[t+1] ~ dbin(phi3, max(50, round(N0m[t]-H0m[t])))
    
    Nadf[t+1] ~ dbin(phi3, max(50, round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t])))
    Nadm[t+1] ~ dbin(phi3, max(50, round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t])))
    }

    for (t in 1:(T+1)){
    Ntot[t] <- round(N0f[t] + N0m[t] + N1f[t] + N1m[t] + Nadf[t] + Nadm[t])      # summing up population vector to population size
    } 
    
    for (t in 1:(T)){
    Xtot[t] <- round(N0f[t]+N0m[t]+N1f[t]+N1m[t]+Nadf[t]+Nadm[t]-H0f[t]-H0m[t]-H1f[t]-H1m[t]-Hadf[t]-Hadm[t])    	# summing up pop vector to post-harvest pop size	
    Xadf[t] <- round(Nadf[t]-Hadf[t])    	
    Xadm[t] <- round(Nadm[t]-Hadm[t])    	
    }
    
    for (t in 1:(T-1))
    {
    lambda[t] <- Ntot[t+1]/Xtot[t]
    }
    
    #Derived parameter - Average before Ntot_2019/Xtot_2018
    #Mlambda18<-mean(lambda[1:(T-4)])

    #################################################################
    # OBSERVATION PROCESS: COUNT DATA
    #################################################################
    
    for(t in 1:(T)){
    y_w[t]~dpois(Xtot[t])
    }
    
    ############################################################
    # OBSERVATION PROCESS: STRUCTURE SURVEY DATA
    ############################################################
    
    for (t in 1:(T-2)){

    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))

    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t]))
    
    }
    
    #Year 2020 = T-1
    for (t in (T-1)){ 
    
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))

    #When T=2020, yearlings were counted together with adult females 
    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1m[t]+N1f[t]+Nadf[t]-H1m[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(Nadm[t]-Hadm[t]))
    
    #When T=2020, subset counted
    C0sub20 ~ dbin(ps*p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cfsub20 ~ dbin(ps*p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    C1sub20 ~ dbin(ps*p2[t], round(N1m[t]-H1m[t]))
    B23sub20 ~ dbin(ps*p2[t], round(Nadm[t]-Hadm[t]))

    }

    #Year = 2021
    for (t in T){
    
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))
    
    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t]))

    } 
    
    ############################################################
    ### Ending model; 
    ###########################################################
    
    } # End Model
    ",fill = TRUE)
sink()


##############################################################
##jagsMod_M2C_HV_fall2021p1.bug ------------------------------
##############################################################
#Include summer minimum counts for some of the last years


sink("jagsMod_M2C_HV_fall2021p1.bug")
cat("
    model {
    
    
    ############################################################
    # Define the priors for the parameters
    ############################################################ 
    
    # Initial pre-harvest population sizes
    N0f[1] ~ dnorm(500, 0.0001)I(0,)     # calves
    N0m[1] ~ dnorm(500, 0.0001)I(0,)     # calves
    N1f[1]~ dnorm(300, 0.0001)I(0,)      # yearling females
    N1m[1]~ dnorm(300, 0.0001)I(0,)      # yearling males
    Nadf[1]~dnorm(500, 0.0001)I(0,)      # adult females
    Nadm[1]~dnorm(500, 0.0001)I(0,)      # adult females
    
    # Survival probabilites and productivity-estimates
    
    for(t in 1:(T-1)){
    p1[t]~ dunif(0,1)
    p2[t]~ dunif(0,1)
    }
    
    #prior for p1[T] - can be specified as informed or vague 
    for(t in T){
    p1[t]~ dbeta(alphapT,betapT)
    p2[t]~ dunif(0,1)
    }
    
    #prior for proportion of animals counted in subset 2020:
    ps~ dunif(0.1,0.2)
    #p1[T+1]~ dunif(0,1)

    ## Vital rates
    phi3 ~ dunif(0, 1)      
    mean.f ~ dunif(-5,5)
    mean.phi1 ~ dunif(-5,5)
    
    tau_f <- 1/(sd_f*sd_f)
    sd_f ~ dunif(0,2)
    
    tau_phi1 <- 1/(sd_phi1*sd_phi1)
    sd_phi1 ~ dunif(0,2)
    
    ##############################
    # DERIVED PARAMETERS
    #############################
    
    fert <- exp(mean.f)/(1+exp(mean.f))
    Sjuv <- exp(mean.phi1)/(1+exp(mean.phi1))  
    
    #############################
    # SYSTEM PROCESS
    #############################
    
    for (t in 1:(T)){  

    logit(f[t]) <- mean.f+eps_f[t]
    eps_f[t]~dnorm(0, tau_f)
    
    logit(phi1[t]) <- mean.phi1+eps_phi1[t]
    eps_phi1[t]~dnorm(0, tau_phi1)

    }



    for (t in 1:(T-2)){  
    N0f[t+1] ~ dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1])) 
    N0m[t+1] ~ dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1]))

    N1f[t+1] ~ dbin(phi3, max(50, round(N0f[t]-H0f[t])))
    N1m[t+1] ~ dbin(phi3, max(50, round(N0m[t]-H0m[t])))
    
    Nadf[t+1] ~ dbin(phi3, max(50, round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t])))
    Nadm[t+1] ~ dbin(phi3, max(50, round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t])))
    }
    
    for (t in 1:T){
    DeadCalves[t]<-round(Nadf[t]*f[t])-(N0f[t]+N0m[t])
    }
    
    #Average over 2001-2018 year 
    Sjuv18<-mean(phi1[2:(T-3)])

    #Use average over 6 last year for next prediction
    fert6m<-mean(f[(T-6):T])
    Sjuv6m<-mean(phi1[(T-6):T])


    for (t in (T-1)){  
    N0f[t+1] ~dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1]))
    N0m[t+1] ~dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1]))

    N1f[t+1] ~ dbin(phi3, max(50, round(N0f[t]-H0f[t])))
    N1m[t+1] ~ dbin(phi3, max(50, round(N0m[t]-H0m[t])))
    
    Nadf[t+1] ~ dbin(phi3, max(50, round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t])))
    Nadm[t+1] ~ dbin(phi3, max(50, round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t])))
    }

    for (t in T){  
    N0f[t+1] ~dbin(Sjuv6m*fert6m/2, max(50, Nadf[t+1]))
    N0m[t+1] ~dbin(Sjuv6m*fert6m/2, max(50, Nadf[t+1]))

    N1f[t+1] ~ dbin(phi3, max(50, round(N0f[t]-H0f[t])))
    N1m[t+1] ~ dbin(phi3, max(50, round(N0m[t]-H0m[t])))
    
    Nadf[t+1] ~ dbin(phi3, max(50, round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t])))
    Nadm[t+1] ~ dbin(phi3, max(50, round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t])))
    
    }

    for (t in 1:(T+1)){
    Ntot[t] <- round(N0f[t] + N0m[t] + N1f[t] + N1m[t] + Nadf[t] + Nadm[t])      # summing up population vector to population size	
    } 
    
    for (t in 1:(T)){
    Xtot[t] <- round(N0f[t]+N0m[t]+N1f[t]+N1m[t]+Nadf[t]+Nadm[t]-H0f[t]-H0m[t]-H1f[t]-H1m[t]-Hadf[t]-Hadm[t])    	# summing up pop vector to post-harvest pop size	
    }
    
    for (t in 1:(T-1))
    {
    lambda[t] <- Ntot[t+1]/Xtot[t]
    }
    
    #Derived parameter - Average before Ntot_2019/Xtot_2018
    Mlambda18<-mean(lambda[1:(T-4)])

    #################################################################
    # OBSERVATION PROCESS: COUNT DATA
    #################################################################
    
    for(t in 1:(T)){
    y_w[t]~dpois(Xtot[t])
    }
    
    y_s[T-2]~dpois(Ntot[T-2]) #2019
    y_s[T]~dpois(Ntot[T]) #2021


    ############################################################
    # OBSERVATION PROCESS: STRUCTURE SURVEY DATA
    ############################################################
    
    for (t in 1:(T-2)){
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))

    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t]))
    }
    
    #Year 2020 = T-1
    for (t in (T-1)){ 
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))

    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    
    #When T=2020, yearlings were counted together with adult females 
    Cf[t] ~ dbin(p2[t], round(N1m[t]+N1f[t]+Nadf[t]-H1m[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(Nadm[t]-Hadm[t]))
    
    #When T=2020, subset counted
    C0sub20 ~ dbin(ps*p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cfsub20 ~ dbin(ps*p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    C1sub20 ~ dbin(ps*p2[t], round(N1m[t]-H1m[t]))
    B23sub20 ~ dbin(ps*p2[t], round(Nadm[t]-Hadm[t]))
    }

    #Year = 2021
    for (t in T){
    ## Last year - with no phi1[t] estimate available - use 
    ##J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/Sjuv6m))
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))
    
    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t]))
    } 
    
    ############################################################
    ### Ending model; 
    ###########################################################
    
    } # End Model
    ",fill = TRUE)
sink()




####################################################
##jagsMod_M2CB23phi2_HV_fall2021p1.bug--------------
####################################################

sink("jagsMod_M2CB23phi2_HV_fall2021p1_E.bug")
cat("
    model {
    
    
    ############################################################
    # Define the priors for the parameters
    ############################################################ 
    
    # Initial pre-harvest population sizes
    #N0f[1] ~ dnorm(500, 0.0001)I(0,)     # calves
    #N0m[1] ~ dnorm(500, 0.0001)I(0,)     # calves
    N1f[1]~ dnorm(300, 0.0001)I(0,)      # yearling females
    N1m[1]~ dnorm(300, 0.0001)I(0,)      # yearling males
    Nadf[1]~dnorm(500, 0.0001)I(0,)      # adult females
    Nadm[1]~dnorm(500, 0.0001)I(0,)      # adult females
    
    # Survival probabilites and productivity-estimates
    
    for(t in 1:(T-1)){
    p1[t]~ dunif(0,1)
    p2[t]~ dunif(0,1)
    }
    
    #prior for p1[T] - can be specified as informed or vague 
    for(t in T){
    p1[t]~ dbeta(alphapT,betapT)
    p2[t]~ dunif(0,1)
    }

    #Prior for subset counted in 2020
    ps~ dunif(0.1,0.2)

    ## Vital rates
    phi3 ~ dunif(0, 1)      
    mean.f ~ dunif(-5,5)
    mean.phi1 ~ dunif(-5,5)
    mean.phi2 ~ dunif(-5,5)

    tau_f <- 1/(sd_f*sd_f)
    sd_f ~ dunif(0,2)
    
    tau_phi1 <- 1/(sd_phi1*sd_phi1)
    sd_phi1 ~ dunif(0,2)
    
    tau_phi2 <- 1/(sd_phi2*sd_phi2)
    sd_phi2 ~ dunif(0,2)

    ##############################
    # DERIVED PARAMETERS
    #############################
    
    fert <- exp(mean.f)/(1+exp(mean.f))
    Sjuv <- exp(mean.phi1)/(1+exp(mean.phi1))  

    #############################
    # SYSTEM PROCESS
    #############################
    
    for (t in 1:(T-1)){  

    logit(f[t]) <- mean.f+eps_f[t]
    eps_f[t]~dnorm(0, tau_f)
    
    logit(phi1[t]) <- mean.phi1+eps_phi1[t]
    eps_phi1[t]~dnorm(0, tau_phi1)
    
    logit(phi2[t]) <- mean.phi2+eps_phi2[t] 

    }
    
    for (t in 1:5)
    {eps_phi2[t]~dnorm(0, 100000)}
    for(t in 6:(T-1)){ 
    eps_phi2[t]~dnorm(0, tau_phi2)
    }


    for (t in T){   
    
    logit(f[t]) <- mean.f+eps_f[t]
    eps_f[t]~dnorm(0, tau_f)
    
    logit(phi1[t]) <- mean.phi1+eps_phi1[t]
    eps_phi1[t]~dnorm(0, tau_phi1)
    
    }


    N0f[1] ~dbin(phi1[1]*f[1]/2, round(max(50, Nadf[1])))
    N0m[1] ~dbin(phi1[1]*f[1]/2, round(max(50, Nadf[1])))

    for (t in 1:(T-1)){  
    N0f[t+1] ~dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1]))
    N0m[t+1] ~dbin(phi1[t+1]*f[t+1]/2, max(50, Nadf[t+1]))
    
    N1f[t+1] ~ dbin(phi2[t], max(50, round(N0f[t]-H0f[t])))
    N1m[t+1] ~ dbin(phi2[t], max(50, round(N0m[t]-H0m[t])))
 
    Nadf[t+1] ~ dbin(phi3, max(50, round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t])))
    Nadm[t+1] ~ dbin(phi3, max(50, round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t])))
    }

    #Average over 2001-2018 year 
    Sjuv18<-mean(phi1[2:(T-3)])
    Sjuvphi2_18<-mean(phi2[1:(T-4)])
    
    #Use average over 6 last year for next prediction
    fert6m<-mean(f[(T-6):T])
    Sjuv6m<-mean(phi1[(T-6):T])
    SjuvW6m<-mean(phi2[(T-7):(T-1)])

    for (t in T){  
    N0f[t+1] ~dbin(Sjuv6m*fert6m/2, max(50, Nadf[t+1]))
    N0m[t+1] ~dbin(Sjuv6m*fert6m/2, max(50, Nadf[t+1]))

    N1f[t+1] ~ dbin(SjuvW6m, max(50, round(N0f[t]-H0f[t]))) 
    N1m[t+1] ~ dbin(SjuvW6m, max(50, round(N0m[t]-H0m[t])))
    
    Nadf[t+1] ~ dbin(phi3, max(50, round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t])))
    Nadm[t+1] ~ dbin(phi3, max(50, round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t])))
    }

    for (t in 1:(T+1)){
    Ntot[t] <- round(N0f[t] + N0m[t] + N1f[t] + N1m[t] + Nadf[t] + Nadm[t])      # summing up population vector to population size	
    } 
    
    for (t in 1:(T)){
    Xtot[t] <- round(N0f[t]+N0m[t]+N1f[t]+N1m[t]+Nadf[t]+Nadm[t]-H0f[t]-H0m[t]-H1f[t]-H1m[t]-Hadf[t]-Hadm[t])    	# summing up to post-harvest pop size
    Xadf[t] <- round(Nadf[t]-Hadf[t])    		
    Xadm[t] <- round(Nadm[t]-Hadm[t])    		
    }
    
    for (t in 1:(T-1))
    {
    lambda[t] <- Ntot[t+1]/Xtot[t]
    }

    #Derived parameter - Average before Ntot_2019/Xtot_2018
    Mlambda18<-mean(lambda[1:(T-4)])

    #################################################################
    # OBSERVATION PROCESS: COUNT DATA
    #################################################################
    
    for(t in 1:(T)){
    y_w[t]~dpois(Xtot[t])
    }
    
    y_s[T-2]~dpois(Ntot[T-2]) #2019
    y_s[T]~dpois(Ntot[T]) #2021


    ############################################################
    # OBSERVATION PROCESS: STRUCTURE SURVEY DATA
    ############################################################
    
    for (t in 1:3){
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))

    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t]))
    }
    
    for (t in 4:(T-3)){
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))

    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t]))
    B23[t] ~ dbin(p2[t], round(Nadm[t]-Hadm[t]))
    }

    #Year 2019
    for (t in (T-2)){
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t])) 

    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t]))
    B23[t] ~ dbin(p2[t], round(Nadm[t]-Hadm[t]))
    }

    #Year 2020
    for (t in (T-1)){
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/phi1[t]))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))

    #When T=2020, yearlings were counted together with adult females 
    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1m[t]+N1f[t]+Nadf[t]-H1m[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(Nadm[t]-Hadm[t]))
    
    #When T=2020, subset counted
    C0sub20 ~ dbin(ps*p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cfsub20 ~ dbin(ps*p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    C1sub20 ~ dbin(ps*p2[t], round(N1m[t]-H1m[t]))
    B23sub20 ~ dbin(ps*p2[t], round(Nadm[t]-Hadm[t]))
    }
    
    #Year 2021
    for (t in (T)){
    J[t] ~ dbin(p1[t], round((N0f[t]+N0m[t])/Sjuv6m))
    SU[t] ~ dbin(p1[t], round(N1m[t]+N1f[t]+Nadf[t]))
    
    C0[t] ~ dbin(p2[t], round(N0f[t]+N0m[t]-H0f[t]-H0m[t])) 
    Cf[t] ~ dbin(p2[t], round(N1f[t]+Nadf[t]-H1f[t]-Hadf[t]))
    Cm[t] ~ dbin(p2[t], round(N1m[t]+Nadm[t]-H1m[t]-Hadm[t]))
    B23[t] ~ dbin(p2[t], round(Nadm[t]-Hadm[t]))
    }
    
    ############################################################
    ### Ending model; 
    ###########################################################
    
    } # End Model
    ",fill = TRUE)
sink()




