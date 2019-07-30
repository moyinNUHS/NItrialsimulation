################################################################################################################
###################Using causal inference to address non-compliance in non inferiority trials###################
################################################################################################################

######## Set up #########
setwd("/Users/moyin/Desktop/VAP studd/Causal inference simulation") #set working directory 
rm(list=ls()) # Clean working environment 

# Required libraries 
library(Hmisc); library(rms); library(gsDesign);library(ivpack)
library(data.table); library(dplyr); library(plotly); library(ggpubr); library(ggplot2) 
library(Matching); library(tableone); library(MatchIt); library(geepack); library(scales)

# The colour blind friendly palette with black:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#set seed for simulation 
set.seed(1234)

#################CASE 1.1 Non-compliance caused by random process (affect both groups) #########################
################################################################################################################
#BIAS
simdata.bias.case1<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # n:                  number of participants per group; 
    # i0:                 probability of intervention = 1 when R=0
    # i1:                 probability of intervention = 1 when R=1
    # p0:                 mortality and recurrence rate in intervention 0; 
    # p1:                 mortality and recurrence rate in intervention 1;
    # confounder.eff.o:   effect of confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.u.eff.o: effect of unknown confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.eff.i:   effect of confounder on intervention (odds ratio);
    # confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio);
    # nIterations:        number of iterations 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()       #average of output from each simulation              
    .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.mpp<-.eff.ps<-c() #output from each simulation
    prop.c<-.prop.c<-c()
    
    #number of data points in simulation 
    interval<-10 
    
    #simulate and derive treatment effect 
    sim<- function() { 
        for(i in 1:interval) { print(paste ("interval value",i))
            for(l in 1:nIterations) { 
                #simulate data frame 
                
                #RANDOMISATION ratio 1:1
                simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                simdata$id<- seq(1,(2*n), by=1) #create participant id 
                
                #CONFOUNDER normal distribution ranging 0-1 
                simdata$confounder<- rep(NA,2*n) 
                simdata$confounder<- rnorm(2*n) 
                simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                
                #CONFOUNDER unknown normal distribution ranging 0-1 
                simdata$confounder.u<- rep(NA,2*n) 
                simdata$confounder.u<- rnorm(2*n) 
                simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                
                #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                
                prob.outcome.1<-    p1+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o 
                simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                
                prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                
                #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                simdata$intervention <- rep(NA,2*n)
                
                compliance<- seq(0.6,1, length.out = interval)*2*n  # set compliance from 60-100%,number of people who complied to protocol
                compliance<- (sample(1:(2*n), compliance[i]))       # IDs of participants who complied to protocol
                for (d in 1: (2*n)) {        
                    if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } 
                }
                
                for (d in 1: (2*n)) {                               # otherwise intervention is opposite of randomisation in the compliant cases 
                    if ((is.na(simdata$intervention))[d]) { 
                        if (simdata$randomisation[d] == 1) {
                            simdata$intervention[d] = 0
                        } else {
                            simdata$intervention [d]= 1
                        }
                    }
                }
                    
                #ACTUAL OUTCOMES depend on intervention
                simdata$outcome<- rep(NA,2*n)
                for (y in 1:(2*n)) {
                    if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                        simdata$outcome[y]<-simdata$outcome0[y]}
                }
                
                #estimating treatment effect 
                ## intention to treat 
                pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                .eff.itt[l]<- pz1.value-pz0.value            
                
                ## per protocol 
                pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                .eff.pp[l] <- p11.value-p00.value   
                
                ## as treated
                pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                .eff.at[l] <- pd1.value-pd0.value               
                
                ## linear binomial generalized linear models with inverse probability weights
                E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                ps<-predict(E.out, type="response")
                sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                .eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                
                ## propensity score matching 
                m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                pscore<- m$fitted.values
                mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                y_trt<-matched$outcome[matched$intervention==1]
                y_con<-matched$outcome[matched$intervention==0]
                diffy<-y_trt-y_con
                .eff.ps[l]<-(t.test(diffy))$estimate
                
                # iv with 2 stage regression
                ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                .eff.iv[l]<-ivmodel$coef[2] 
                
            }
            
            #mean of iterated data 
            eff.itt[i]<- mean(.eff.itt, na.rm=TRUE)
            eff.pp[i] <- mean(.eff.pp, na.rm=TRUE) 
            eff.at[i] <- mean(.eff.at, na.rm=TRUE)
            eff.iv[i] <- mean(.eff.iv, na.rm = TRUE)
            eff.mpp[i]<- mean(.eff.mpp, na.rm = TRUE)
            eff.ps[i] <- mean(.eff.ps, na.rm = TRUE)
            
            # mean proportion of patients who were compliant to protocol 
            prop.c [i]<- mean(.prop.c)
            
        }
        return(data.frame(prop.c, eff.iv, eff.itt, eff.pp, eff.at, eff.mpp, eff.ps))
    }
    
    sim.df<-sim()
    
    bias.plot <- ggplot(sim.df, aes(x=sim.df[1]))+ 
        geom_point(aes(y=sim.df[2], colour="Instrumental variable",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.4)) + 
        geom_point(aes(y=sim.df[4], colour="Per protocol",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[5], colour="As treated",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) + 
        geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
        geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y=sim.df[interval,4], colour='True effect'),linetype="dotted", size=2) +
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="bottom")+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.02, 0.12))
    
    return(bias.plot)
} 

#TYPE 1 ERROR
simdata.type1.case1<- function(n, i0, i1, NImargin, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # NImargin: non inferiority margin 
    
    #lower critical value given p1 (one sided of 95% CI = alpha= 0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE) 
    lower.value<-lower$conf.int [1] #95% CI of NI margin 
    
    #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
    p0=p1-NImargin 
    
    #number of data points in simulation 
    interval<-10 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
    type1.error.at<-type1.error.itt<-type1.error.iv<-type1.error.mpp<-type1.error.pp<-type1.error.ps<-c() 
    prop.c<-.prop.c<-c()
    
    if (NImargin>0) {
        
        #simulate data
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    compliance<- seq(0.6,1, length.out = interval)*2*n  # set compliance from 60-100% 
                    compliance<- (sample(1:(2*n), compliance[i]))
                    for (d in 1: (2*n)) {        
                        if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } 
                    }
                    
                    for (d in 1: (2*n)) {                               # otherwise intervention is opposite of randomisation in the compliant cases 
                        if ((is.na(simdata$intervention))[d]) { 
                            if (simdata$randomisation[d] == 1) {
                                simdata$intervention[d] = 0
                            } else {
                                simdata$intervention [d]= 1
                            }
                        }
                    }
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                # mean of type 1 error from iterated data 
                type1.error.itt[i]<-mean(eff.itt<lower.value)
                type1.error.pp[i] <- mean(eff.pp<lower.value)
                type1.error.at[i] <- mean(eff.at<lower.value)
                type1.error.mpp[i]<- mean(eff.mpp<lower.value)
                type1.error.iv[i] <- mean(eff.iv<lower.value)
                type1.error.ps[i] <- mean(eff.ps<lower.value)
                
                # mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c[l])
                
            }
            
            return(data.frame(prop.c, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.mpp, type1.error.ps))
            
        }
        
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.4)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat",alpha=0.4)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.4)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.4)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Type 1 error")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.00,0.4))
        
        return(plot)
        
    } else print("NI margin should be positive in this stimulation")
}

#POWER
simdata.power.case1<- function(n, i0, i1, p1, p0, confounder.eff.o, confounder.eff.i, confounder.u.eff.i,confounder.u.eff.o, NImargin, nIterations){  
    
    #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE)
    lower.value<-lower$conf.int [1]
    
    if ((NImargin>0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
        
        #make up vectors for simulations
        power.iv<-power.itt<- power.at<-power.pp<-power.mpp<-power.ps<-c() 
        eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()    
        prop.c<-.prop.c<-c()
        
        #number of data points in simulation 
        interval<-10
        
        #simulate and derive treatment effect 
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    compliance<- seq(0.6,1, length.out = interval)*2*n  # set compliance from 60-100% 
                    compliance<- (sample(1:(2*n), compliance[i]))
                    for (d in 1: (2*n)) {        
                        if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } 
                    }
                    
                    for (d in 1: (2*n)) {                               # otherwise intervention is opposite of randomisation in the compliant cases 
                        if ((is.na(simdata$intervention))[d]) { 
                            if (simdata$randomisation[d] == 1) {
                                simdata$intervention[d] = 0
                            } else {
                                simdata$intervention [d]= 1
                            }
                        }
                    }
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                #mean power from iterated data
                power.iv[i]<- mean(eff.iv>lower.value)
                power.itt[i]<-mean(eff.itt>lower.value)
                power.at[i]<- mean(eff.at>lower.value)
                power.pp[i]<- mean(eff.pp>lower.value)
                power.mpp[i]<-mean(eff.mpp>lower.value)
                power.ps[i]<- mean(eff.ps>lower.value)
                
                #mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c)
            }
            return(data.frame(prop.c, power.iv, power.itt, power.pp, power.at, power.mpp, power.ps))
        }
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.4)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.4)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.4)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.4)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Power")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.2, 1))
        
        return(plot)
    }
    else (print ("NI margin must be positive, and p1-p0 (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

#Run simulations for case 1.1
b1.1<-simdata.bias.case1   (n=300, p1=0.4, p0=0.3, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
t1.1<-simdata.type1.case1  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p1.1<- simdata.power.case1 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case1<-ggarrange(b1.1, t1.1, p1.1, ncol = 3, nrow = 1,labels = c(paste(LETTERS[1:3], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(case1, top = text_grob("Non-compliance caused by random processes affecting both groups", face = "bold", size = 12))

#################CASE 1.2 Non-compliance caused by random process (affect intervention group) ##################
################################################################################################################
#BIAS
simdata.bias.case2<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # n:                  number of participants per group; 
    # i0:                 probability of intervention = 1 when R=0
    # i1:                 probability of intervention = 1 when R=1
    # p0:                 mortality and recurrence rate in intervention 0; 
    # p1:                 mortality and recurrence rate in intervention 1;
    # confounder.eff.o:   effect of confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.u.eff.o: effect of unknown confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.eff.i:   effect of confounder on intervention (odds ratio);
    # confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio);
    # nIterations:        number of iterations 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
    .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.mpp<-.eff.ps<-c() 
    prop.c<-.prop.c<-c()
    
    #number of data points in simulation 
    interval<-10 
    
    #simulate and derive treatment effect 
    sim<- function() { 
        for(i in 1:interval) { print(paste ("interval value",i))
            for(l in 1:nIterations) { 
                #simulate data frame 
                
                #RANDOMISATION ratio 1:1
                simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                simdata$id<- seq(1,(2*n), by=1) #create participant id 
                
                #CONFOUNDER normal distribution ranging 0-1 
                simdata$confounder<- rep(NA,2*n) 
                simdata$confounder<- rnorm(2*n) 
                simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                
                #CONFOUNDER nknown normal distribution ranging 0-1 
                simdata$confounder.u<- rep(NA,2*n) 
                simdata$confounder.u<- rnorm(2*n) 
                simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                
                #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                
                prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                
                prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                
                #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                simdata$intervention <- rep(NA,2*n)
                
                for (d in 1: (2*n)) {        
                    if (simdata$randomisation[d]==0) {            # intervention is the same as randomisation in the control cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } 
                }
                
                compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                compliance<- (sample(1:n, compliance[i]))
                
                for (d in 1:n) { 
                    if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } else {
                        simdata$intervention[d]=0                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                    }
                }
            
                #ACTUAL OUTCOMES depend on intervention
                simdata$outcome<- rep(NA,2*n)
                for (y in 1:(2*n)) {
                    if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                        simdata$outcome[y]<-simdata$outcome0[y]}
                }
                
                #estimating treatment effect 
                ## intention to treat 
                pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                .eff.itt[l]<- pz1.value-pz0.value            
                
                ## per protocol 
                pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                .eff.pp[l] <- p11.value-p00.value   
                
                ## as treated
                pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                .eff.at[l] <- pd1.value-pd0.value               
                
                ## linear binomial generalized linear models with inverse probability weights
                E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                ps<-predict(E.out, type="response")
                sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                .eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                
                ## propensity score matching 
                m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                pscore<- m$fitted.values
                mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                y_trt<-matched$outcome[matched$intervention==1]
                y_con<-matched$outcome[matched$intervention==0]
                diffy<-y_trt-y_con
                .eff.ps[l]<-(t.test(diffy))$estimate
                
                # iv with 2 stage regression
                ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                .eff.iv[l]<-ivmodel$coef[2] 
                
            }
            
            #mean of iterated data 
            eff.itt[i]<- mean(.eff.itt, na.rm=TRUE)
            eff.pp[i] <- mean(.eff.pp, na.rm=TRUE) 
            eff.at[i] <- mean(.eff.at, na.rm=TRUE)
            eff.iv[i] <- mean(.eff.iv, na.rm = TRUE)
            eff.mpp[i]<- mean(.eff.mpp, na.rm = TRUE)
            eff.ps[i] <- mean(.eff.ps, na.rm = TRUE)
            
            # mean proportion of patients who were compliant to protocol 
            prop.c [i]<- mean(.prop.c)
            
        }
        return(data.frame(prop.c, eff.iv, eff.itt, eff.pp, eff.at, eff.mpp, eff.ps))
    }
    
    sim.df<-sim()
    
    bias.plot <- ggplot(sim.df, aes(x=sim.df[1]))+ 
        geom_point(aes(y=sim.df[2], colour="Instrumental variable",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.4)) + 
        geom_point(aes(y=sim.df[4], colour="Per protocol",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[5], colour="As treated",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) + 
        geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
        geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y=sim.df[interval,4], colour='True effect'),linetype="dotted", size=2) +
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="bottom")+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.02, 0.12))
    
    return(bias.plot)
} 

#TYPE 1 ERROR
simdata.type1.case2<- function(n, i0, i1, NImargin, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # NImargin: non inferiority margin 
    
    #lower critical value given p1 (one sided of 95% CI = alpha= 0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE) 
    lower.value<-lower$conf.int [1] #95% CI of NI margin 
    
    #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
    p0=p1-NImargin 
    
    #number of data points in simulation 
    interval<-10 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
    type1.error.at<-type1.error.itt<-type1.error.iv<-type1.error.mpp<-type1.error.pp<-type1.error.ps<-c() 
    prop.c<-.prop.c<-c()
    
    if (NImargin>0) {
        
        #simulate data
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    for (d in 1: (2*n)) {        
                        if (simdata$randomisation[d]==0) {            # intervention is the same as randomisation in the control cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } 
                    }
                    
                    compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                    compliance<- (sample(1:n, compliance[i]))
                    
                    for (d in 1:n) { 
                        if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } else {
                            simdata$intervention[d]=0                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                        }
                    }
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                # mean of type 1 error from iterated data 
                type1.error.itt[i]<-mean(eff.itt<lower.value)
                type1.error.pp[i] <- mean(eff.pp<lower.value)
                type1.error.at[i] <- mean(eff.at<lower.value)
                type1.error.mpp[i]<- mean(eff.mpp<lower.value)
                type1.error.iv[i] <- mean(eff.iv<lower.value)
                type1.error.ps[i] <- mean(eff.ps<lower.value)
                
                # mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c[l])
                
            }
            
            return(data.frame(prop.c, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.mpp, type1.error.ps))
            
        }
        
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.4)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat",alpha=0.4)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.4)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.4)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Type 1 error")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.00,0.4))
        
        return(plot)
        
    } else print("NI margin should be positive in this stimulation")
}

#POWER
simdata.power.case2<- function(n, i0, i1, p1, p0, confounder.eff.o, confounder.eff.i, confounder.u.eff.i,confounder.u.eff.o, NImargin, nIterations){  
    
    #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE)
    lower.value<-lower$conf.int [1]
    
    if ((NImargin>0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
        
        #make up vectors for simulations
        power.iv<-power.itt<- power.at<-power.pp<-power.mpp<-power.ps<-c() 
        eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()    
        prop.c<-.prop.c<-c()
        
        #number of data points in simulation 
        interval<-10
        
        #simulate and derive treatment effect 
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    for (d in 1: (2*n)) {        
                        if (simdata$randomisation[d]==0) {            # intervention is the same as randomisation in the control cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } 
                    }
                    
                    compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                    compliance<- (sample(1:n, compliance[i]))
                    
                    for (d in 1:n) { 
                        if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } else {
                            simdata$intervention[d]=0                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                        }
                    }
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                #mean power from iterated data
                power.iv[i]<- mean(eff.iv>lower.value)
                power.itt[i]<-mean(eff.itt>lower.value)
                power.at[i]<- mean(eff.at>lower.value)
                power.pp[i]<- mean(eff.pp>lower.value)
                power.mpp[i]<-mean(eff.mpp>lower.value)
                power.ps[i]<- mean(eff.ps>lower.value)
                
                #mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c)
            }
            return(data.frame(prop.c, power.iv, power.itt, power.pp, power.at, power.mpp, power.ps))
        }
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.4)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.4)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.4)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.4)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Power")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.2, 1))
        
        return(plot)
    }
    else (print ("NI margin must be positive, and p1-p0 (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

#Run simulations for case 1.2
b1.2<-simdata.bias.case2   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
t1.2<-simdata.type1.case2  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p1.2<- simdata.power.case2 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case2<-ggarrange(b1.2, t1.2, p1.2, ncol = 3, nrow = 1,labels = c(paste(LETTERS[1:3], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(case2, top = text_grob("Non-compliance caused by random processes affecting intervention or control group", face = "bold", size = 12))

#################CASE 1.3 Non-compliance caused by random process (affect control group) #######################
##################################################Not shown#####################################################
#BIAS
simdata.bias.case3<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # n:                  number of participants per group; 
    # i0:                 probability of intervention = 1 when R=0
    # i1:                 probability of intervention = 1 when R=1
    # p0:                 mortality and recurrence rate in intervention 0; 
    # p1:                 mortality and recurrence rate in intervention 1;
    # confounder.eff.o:   effect of confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.u.eff.o: effect of unknown confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.eff.i:   effect of confounder on intervention (odds ratio);
    # confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio);
    # nIterations:        number of iterations 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
    .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.mpp<-.eff.ps<-c() 
    prop.c<-.prop.c<-c()
    
    #number of data points in simulation 
    interval<-10 
    
    #simulate and derive treatment effect 
    sim<- function() { 
        for(i in 1:interval) { print(paste ("interval value",i))
            for(l in 1:nIterations) { 
                #simulate data frame 
                
                #RANDOMISATION ratio 1:1
                simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                simdata$id<- seq(1,(2*n), by=1) #create participant id 
                
                #CONFOUNDER normal distribution ranging 0-1 
                simdata$confounder<- rep(NA,2*n) 
                simdata$confounder<- rnorm(2*n) 
                simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                
                #CONFOUNDER nknown normal distribution ranging 0-1 
                simdata$confounder.u<- rep(NA,2*n) 
                simdata$confounder.u<- rnorm(2*n) 
                simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                
                #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                
                prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                
                prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                
                #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                simdata$intervention <- rep(NA,2*n)
                
                for (d in 1: (2*n)) {        
                    if (simdata$randomisation[d]==1) {            # intervention is the same as randomisation in the intervention cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } 
                }
                
                compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                compliance<- (sample((n+1):(2*n), compliance[i]))
                
                for (d in (n+1):(2*n)) { 
                    if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                        simdata$intervention[d] = simdata$randomisation[d]
                    } else {
                        simdata$intervention[d]=1                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                    }
                }
                
                #ACTUAL OUTCOMES depend on intervention
                simdata$outcome<- rep(NA,2*n)
                for (y in 1:(2*n)) {
                    if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                        simdata$outcome[y]<-simdata$outcome0[y]}
                }
                
                #estimating treatment effect 
                ## intention to treat 
                pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                .eff.itt[l]<- pz1.value-pz0.value            
                
                ## per protocol 
                pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                .eff.pp[l] <- p11.value-p00.value   
                
                ## as treated
                pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                .eff.at[l] <- pd1.value-pd0.value               
                
                ## linear binomial generalized linear models with inverse probability weights
                E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                ps<-predict(E.out, type="response")
                sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                .eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                
                ## propensity score matching 
                m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                pscore<- m$fitted.values
                mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                y_trt<-matched$outcome[matched$intervention==1]
                y_con<-matched$outcome[matched$intervention==0]
                diffy<-y_trt-y_con
                .eff.ps[l]<-(t.test(diffy))$estimate
                
                # iv with 2 stage regression
                ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                .eff.iv[l]<-ivmodel$coef[2] 
                
            }
            
            #mean of iterated data 
            eff.itt[i]<- mean(.eff.itt, na.rm=TRUE)
            eff.pp[i] <- mean(.eff.pp, na.rm=TRUE) 
            eff.at[i] <- mean(.eff.at, na.rm=TRUE)
            eff.iv[i] <- mean(.eff.iv, na.rm = TRUE)
            eff.mpp[i]<- mean(.eff.mpp, na.rm = TRUE)
            eff.ps[i] <- mean(.eff.ps, na.rm = TRUE)
            
            # mean proportion of patients who were compliant to protocol 
            prop.c [i]<- mean(.prop.c)
            
        }
        return(data.frame(prop.c, eff.iv, eff.itt, eff.pp, eff.at, eff.mpp, eff.ps))
    }
    
    sim.df<-sim()
    
    bias.plot <- ggplot(sim.df, aes(x=sim.df[1]))+ 
        geom_point(aes(y=sim.df[2], colour="Instrumental variable",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.4)) + 
        geom_point(aes(y=sim.df[4], colour="Per protocol",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[5], colour="As treated",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) + 
        geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
        geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y=sim.df[interval,4], colour='True effect'),linetype="dotted", size=2) +
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="bottom")+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(-0.02, 0.12))
    
    return(bias.plot)
} 

#TYPE 1 ERROR
simdata.type1.case3<- function(n, i0, i1, NImargin, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # NImargin: non inferiority margin 
    
    #lower critical value given p1 (one sided of 95% CI = alpha= 0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE) 
    lower.value<-lower$conf.int [1] #95% CI of NI margin 
    
    #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
    p0=p1-NImargin 
    
    #number of data points in simulation 
    interval<-10 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
    type1.error.at<-type1.error.itt<-type1.error.iv<-type1.error.mpp<-type1.error.pp<-type1.error.ps<-c() 
    prop.c<-.prop.c<-c()
    
    if (NImargin>0) {
        
        #simulate data
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    for (d in 1: (2*n)) {        
                        if (simdata$randomisation[d]==1) {            # intervention is the same as randomisation in the intervention cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } 
                    }
                    
                    compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                    compliance<- (sample((n+1):(2*n), compliance[i]))
                    
                    for (d in (n+1):(2*n)) { 
                        if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } else {
                            simdata$intervention[d]=1                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                        }
                    }
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                # mean of type 1 error from iterated data 
                type1.error.itt[i]<-mean(eff.itt<lower.value)
                type1.error.pp[i] <- mean(eff.pp<lower.value)
                type1.error.at[i] <- mean(eff.at<lower.value)
                type1.error.mpp[i]<- mean(eff.mpp<lower.value)
                type1.error.iv[i] <- mean(eff.iv<lower.value)
                type1.error.ps[i] <- mean(eff.ps<lower.value)
                
                # mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c[l])
                
            }
            
            return(data.frame(prop.c, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.mpp, type1.error.ps))
            
        }
        
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.4)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat",alpha=0.4)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.4)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.4)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Type 1 error")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.00,0.4))
        
        return(plot)
        
    } else print("NI margin should be positive in this stimulation")
}

#POWER
simdata.power.case3<- function(n, i0, i1, p1, p0, confounder.eff.o, confounder.eff.i, confounder.u.eff.i,confounder.u.eff.o, NImargin, nIterations){  
    
    #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE)
    lower.value<-lower$conf.int [1]
    
    if ((NImargin>0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
        
        #make up vectors for simulations
        power.iv<-power.itt<- power.at<-power.pp<-power.mpp<-power.ps<-c() 
        eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()    
        prop.c<-.prop.c<-c()
        
        #number of data points in simulation 
        interval<-10
        
        #simulate and derive treatment effect 
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and a random process - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    for (d in 1: (2*n)) {        
                        if (simdata$randomisation[d]==1) {            # intervention is the same as randomisation in the intervention cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } 
                    }
                    
                    compliance<- seq(0.2,1, length.out = interval)*n  # set compliance from 60-100% 
                    compliance<- (sample((n+1):(2*n), compliance[i]))
                    
                    for (d in (n+1):(2*n)) { 
                        if (simdata$id[d] %in% compliance) {            # intervention is the same as randomisation in the compliant cases 
                            simdata$intervention[d] = simdata$randomisation[d]
                        } else {
                            simdata$intervention[d]=1                   # otherwise intervention is opposite of randomisation in the non-compliant cases 
                        }
                    }
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                #mean power from iterated data
                power.iv[i]<- mean(eff.iv>lower.value)
                power.itt[i]<-mean(eff.itt>lower.value)
                power.at[i]<- mean(eff.at>lower.value)
                power.pp[i]<- mean(eff.pp>lower.value)
                power.mpp[i]<-mean(eff.mpp>lower.value)
                power.ps[i]<- mean(eff.ps>lower.value)
                
                #mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c)
            }
            return(data.frame(prop.c, power.iv, power.itt, power.pp, power.at, power.mpp, power.ps))
        }
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.4)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.4)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.4)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.4)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Power")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.2, 1))
        
        return(plot)
    }
    else (print ("NI margin must be positive, and p1-p0 (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

#Run simulations for case 1.3
b1.3<-simdata.bias.case3   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
t1.3<-simdata.type1.case3  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p1.3<- simdata.power.case3 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1, confounder.eff.o = 0,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case3<-ggarrange(b1.3, t1.3, p1.3, ncol = 3, nrow = 1,labels = c(paste(LETTERS[1:3], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(case3, top = text_grob("Case 1.3: Non-compliance caused by random processes affecting control group", face = "bold", size = 12))

############################CASE 2 Non-compliance caused by non random process  ################################
################################################################################################################
#BIAS
simdata.bias.case4<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # n:                  number of participants per group; 
    # i0:                 probability of intervention = 1 when R=0
    # i1:                 probability of intervention = 1 when R=1
    # p0:                 mortality and recurrence rate in intervention 0; 
    # p1:                 mortality and recurrence rate in intervention 1;
    # confounder.eff.o:   effect of confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.u.eff.o: effect of unknown confounder on outcome (absolute difference in probability of outcome in presence of confounder);
    # confounder.eff.i:   effect of confounder on intervention (odds ratio);
    # confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio);
    # nIterations:        number of iterations 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
    .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.mpp<-.eff.ps<-c() 
    prop.c<-.prop.c<-c()
    
    #number of data points in simulation 
    interval<-10 
    
    #simulate and derive treatment effect 
    sim<- function() { 
        for(i in 1:interval) { print(paste ("interval value",i))
            for(l in 1:nIterations) { 
                #simulate data frame 
                
                #RANDOMISATION ratio 1:1
                simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                simdata$id<- seq(1,(2*n), by=1) #create participant id 
                
                #CONFOUNDER normal distribution ranging 0-1 
                simdata$confounder<- rep(NA,2*n) 
                simdata$confounder<- rnorm(2*n) 
                simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                
                #CONFOUNDER nknown normal distribution ranging 0-1 
                simdata$confounder.u<- rep(NA,2*n) 
                simdata$confounder.u<- rnorm(2*n) 
                simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                
                #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                
                prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                
                prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                
                #INTERVENTION dependent on randomisation and confounders by logistric regression - creating increasing compliance
                simdata$intervention <- rep(NA,2*n)
                
                .b0<-logit (i0)
                b0<- seq (.b0,logit(0.01), length.out = interval)
                .b1<-log (confounder.eff.i)
                b1<- seq (.b1,0, length.out = interval)
                .b2<-log (confounder.u.eff.i)
                b2<- seq (.b2,0, length.out = interval)
                .b3<-logit (i1)- logit (i0)
                b3<- seq(.b3,log((0.99/(1-0.99))/(0.01/(1-0.01))), length.out = interval)
                
                if (b1[i]==0) { pi<- simdata$randomisation} else {
                    logit.pi<- b0[i]+b1[i]*simdata$confounder+b2[i]*simdata$confounder.u+b3[i]*simdata$randomisation
                    pi<-exp(logit.pi)/(1+exp(logit.pi))
                }
                simdata$intervention <- rbinom(2*n,1,prob=pi)
                
                #ACTUAL OUTCOMES depend on intervention
                simdata$outcome<- rep(NA,2*n)
                for (y in 1:(2*n)) {
                    if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                        simdata$outcome[y]<-simdata$outcome0[y]}
                }
                
                #estimating treatment effect 
                ## intention to treat 
                pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                .eff.itt[l]<- pz1.value-pz0.value            
                
                ## per protocol 
                pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                .eff.pp[l] <- p11.value-p00.value   
                
                ## as treated
                pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                .eff.at[l] <- pd1.value-pd0.value               
                
                ## linear binomial generalized linear models with inverse probability weights
                E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                ps<-predict(E.out, type="response")
                sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                .eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                
                ## propensity score matching 
                m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                pscore<- m$fitted.values
                mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                y_trt<-matched$outcome[matched$intervention==1]
                y_con<-matched$outcome[matched$intervention==0]
                diffy<-y_trt-y_con
                .eff.ps[l]<-(t.test(diffy))$estimate
                
                # iv with 2 stage regression
                ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                .eff.iv[l]<-ivmodel$coef[2] 
                
            }
            
            #mean of iterated data 
            eff.itt[i]<- mean(.eff.itt, na.rm=TRUE)
            eff.pp[i] <- mean(.eff.pp, na.rm=TRUE) 
            eff.at[i] <- mean(.eff.at, na.rm=TRUE)
            eff.iv[i] <- mean(.eff.iv, na.rm = TRUE)
            eff.mpp[i]<- mean(.eff.mpp, na.rm = TRUE)
            eff.ps[i] <- mean(.eff.ps, na.rm = TRUE)
            
            # mean proportion of patients who were compliant to protocol 
            prop.c [i]<- mean(.prop.c)
            
        }
        return(data.frame(prop.c, eff.iv, eff.itt, eff.pp, eff.at, eff.mpp, eff.ps))
    }
    
    sim.df<-sim()
    
    bias.plot <- ggplot(sim.df, aes(x=sim.df[1]))+ 
        geom_point(aes(y=sim.df[2], colour="Instrumental variable",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.4)) + 
        geom_point(aes(y=sim.df[4], colour="Per protocol",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[5], colour="As treated",  alpha=0.4)) + 
        geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) + 
        geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
        geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y=sim.df[interval,4], colour='True effect'),linetype="dotted", size=2) +
        xlab("Proportion of compliant participants")+
        ylab("Effect")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.position="bottom")+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.60, 1))+
        scale_y_continuous(limits=c(0, 0.2))
    
    return(bias.plot)
} 

#TYPE 1 ERROR
simdata.type1.case4<- function(n, i0, i1, NImargin, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
    # NImargin: non inferiority margin 
    
    #lower critical value given p1 (one sided of 95% CI = alpha= 0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE) 
    lower.value<-lower$conf.int [1] #95% CI of NI margin 
    
    #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
    p0=p1-NImargin 
    
    #number of data points in simulation 
    interval<-10 
    
    #make up vectors for simulations 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
    type1.error.at<-type1.error.itt<-type1.error.iv<-type1.error.mpp<-type1.error.pp<-type1.error.ps<-c() 
    prop.c<-.prop.c<-c()
    
    if (NImargin>0) {
        
        #simulate data
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and confounders by logistric regression - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    .b0<-logit (i0)
                    b0<- seq (.b0,logit(0.01), length.out = interval)
                    .b1<-log (confounder.eff.i)
                    b1<- seq (.b1,0, length.out = interval)
                    .b2<-log (confounder.u.eff.i)
                    b2<- seq (.b2,0, length.out = interval)
                    .b3<-logit (i1)- logit (i0)
                    b3<- seq(.b3,log((0.99/(1-0.99))/(0.01/(1-0.01))), length.out = interval)
                    
                    if (b1[i]==0) { pi<- simdata$randomisation} else {
                        logit.pi<- b0[i]+b1[i]*simdata$confounder+b2[i]*simdata$confounder.u+b3[i]*simdata$randomisation
                        pi<-exp(logit.pi)/(1+exp(logit.pi))
                    }
                    simdata$intervention <- rbinom(2*n,1,prob=pi)
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                # mean of type 1 error from iterated data 
                type1.error.itt[i]<-mean(eff.itt<lower.value)
                type1.error.pp[i] <- mean(eff.pp<lower.value)
                type1.error.at[i] <- mean(eff.at<lower.value)
                type1.error.mpp[i]<- mean(eff.mpp<lower.value)
                type1.error.iv[i] <- mean(eff.iv<lower.value)
                type1.error.ps[i] <- mean(eff.ps<lower.value)
                
                # mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c[l])
                
            }
            
            return(data.frame(prop.c, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.mpp, type1.error.ps))
            
        }
        
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.4)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat",alpha=0.4)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.4)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.4)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Type 1 error")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.00,0.4))
        
        return(plot)
        
    } else print("NI margin should be positive in this stimulation")
}

#POWER
simdata.power.case4<- function(n, i0, i1, p1, p0, confounder.eff.o, confounder.eff.i, confounder.u.eff.i,confounder.u.eff.o, NImargin, nIterations){  
    
    #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
    lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE)
    lower.value<-lower$conf.int [1]
    
    if ((NImargin>0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
        
        #make up vectors for simulations
        power.iv<-power.itt<- power.at<-power.pp<-power.mpp<-power.ps<-c() 
        eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()    
        prop.c<-.prop.c<-c()
        
        #number of data points in simulation 
        interval<-10
        
        #simulate and derive treatment effect 
        sim<- function() { 
            for(i in 1:interval) { print(paste ("interval value",i))
                for(l in 1:nIterations) { 
                    #simulate data frame 
                    
                    #RANDOMISATION ratio 1:1
                    simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
                    simdata$id<- seq(1,(2*n), by=1) #create participant id 
                    
                    #CONFOUNDER normal distribution ranging 0-1 
                    simdata$confounder<- rep(NA,2*n) 
                    simdata$confounder<- rnorm(2*n) 
                    simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #CONFOUNDER nknown normal distribution ranging 0-1 
                    simdata$confounder.u<- rep(NA,2*n) 
                    simdata$confounder.u<- rnorm(2*n) 
                    simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
                    
                    #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
                    simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
                    simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
                    
                    prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
                    simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
                    
                    prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
                    simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
                    
                    #INTERVENTION dependent on randomisation and confounders by logistric regression - creating increasing compliance
                    simdata$intervention <- rep(NA,2*n)
                    
                    .b0<-logit (i0)
                    b0<- seq (.b0,logit(0.01), length.out = interval)
                    .b1<-log (confounder.eff.i)
                    b1<- seq (.b1,0, length.out = interval)
                    .b2<-log (confounder.u.eff.i)
                    b2<- seq (.b2,0, length.out = interval)
                    .b3<-logit (i1)- logit (i0)
                    b3<- seq(.b3,log((0.99/(1-0.99))/(0.01/(1-0.01))), length.out = interval)
                    
                    if (b1[i]==0) { pi<- simdata$randomisation} else {
                        logit.pi<- b0[i]+b1[i]*simdata$confounder+b2[i]*simdata$confounder.u+b3[i]*simdata$randomisation
                        pi<-exp(logit.pi)/(1+exp(logit.pi))
                    }
                    simdata$intervention <- rbinom(2*n,1,prob=pi)
                    
                    #ACTUAL OUTCOMES depend on intervention
                    simdata$outcome<- rep(NA,2*n)
                    for (y in 1:(2*n)) {
                        if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                            simdata$outcome[y]<-simdata$outcome0[y]}
                    }
                    
                    #estimating treatment effect 
                    ## intention to treat 
                    pz1.value<- mean((filter(simdata, randomisation==1))$outcome)                 
                    pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
                    eff.itt[l]<- pz1.value-pz0.value            
                    
                    ## per protocol 
                    pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
                    .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol 
                    p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome) 
                    p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
                    eff.pp[l] <- p11.value-p00.value   
                    
                    ## as treated
                    pd1.value<- mean((filter(simdata, intervention==1))$outcome)                   
                    pd0.value<- mean((filter(simdata, intervention==0))$outcome)
                    eff.at[l] <- pd1.value-pd0.value               
                    
                    ## linear binomial generalized linear models with inverse probability weights
                    E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
                    ps<-predict(E.out, type="response")
                    sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
                    eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
                    
                    ## propensity score matching 
                    m<- glm(intervention~confounder, family = binomial(), data=simdata) 
                    pscore<- m$fitted.values
                    mout<-matchit(intervention~confounder, data=simdata, method='nearest')
                    match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
                    matched<-simdata[unlist(match[c('index.treated','index.control')]),]
                    y_trt<-matched$outcome[matched$intervention==1]
                    y_con<-matched$outcome[matched$intervention==0]
                    diffy<-y_trt-y_con
                    eff.ps[l]<-(t.test(diffy))$estimate
                    
                    # iv with 2 stage regression
                    ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) 
                    eff.iv[l]<-ivmodel$coef[2] 
                    
                }
                
                #mean power from iterated data
                power.iv[i]<- mean(eff.iv>lower.value)
                power.itt[i]<-mean(eff.itt>lower.value)
                power.at[i]<- mean(eff.at>lower.value)
                power.pp[i]<- mean(eff.pp>lower.value)
                power.mpp[i]<-mean(eff.mpp>lower.value)
                power.ps[i]<- mean(eff.ps>lower.value)
                
                #mean proportion of patients who were compliant to protocol 
                prop.c [i]<- mean(.prop.c)
            }
            return(data.frame(prop.c, power.iv, power.itt, power.pp, power.at, power.mpp, power.ps))
        }
        sim.df<-sim()
        
        plot <- ggplot(sim.df, aes(sim.df[1]))+ 
            geom_point(aes(y=sim.df[2], colour="Instrumental variable", alpha=0.4)) + 
            geom_point(aes(y=sim.df[3], colour="Intention to treat", alpha=0.4)) + 
            geom_point(aes(y=sim.df[4], colour="Per protocol", alpha=0.4)) + 
            geom_point(aes(y=sim.df[5], colour="As treated", alpha=0.4)) + 
            geom_point(aes(y=sim.df[6], colour="Modified per protocol", alpha=0.4)) +
            geom_point(aes(y=sim.df[7], colour="Propensity score matching", alpha=0.4)) + 
            geom_smooth(aes(y=sim.df[2], colour="Instrumental variable"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[3], colour="Intention to treat"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[4], colour="Per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[5], colour="As treated"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[6], colour="Modified per protocol"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            geom_smooth(aes(y=sim.df[7], colour="Propensity score matching"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
            xlab("Proportion of compliant participants")+
            ylab("Power")+
            theme_bw()+
            scale_colour_manual(values=cbPalette)+
            theme(legend.position="bottom")+
            theme(legend.title=element_blank())+
            scale_x_continuous(limits=c(0.6, 1))+
            scale_y_continuous(limits=c(0.2, 1))
        
        return(plot)
    }
    else (print ("NI margin must be positive, and p1-p0 (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

#Run simulations for case 2.1: Non-compliance caused by non random process (affect both groups)
## higher value of confounder makes intervention less likely, outcome more likely
b2.1<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3, i0=0.45, i1=0.55, confounder.eff.i= 0.5,     confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t2.1<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.5, confounder.eff.o = 0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.1<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.5, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
## higher value of confounder makes intervention more likely, outcome more likely
b2.12<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3, i0=0.45, i1=0.55, confounder.eff.i= 2,     confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t2.12<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 2, confounder.eff.o = 0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.12<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 2, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
## higher value of confounder makes outcome less likely, outcome less likely
b2.13<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3, i0=0.45, i1=0.55, confounder.eff.i= 0.5,     confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t2.13<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.5, confounder.eff.o = -0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.13<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.5, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
## higher value of confounder makes outcome more likely, outcome less likely
b2.14<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3, i0=0.45, i1=0.55, confounder.eff.i= 2,     confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t2.14<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 2, confounder.eff.o = -0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.14<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 2, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)

case4<-ggarrange(b2.1, t2.1, p2.1, 
                 b2.12,t2.12,p2.12,
                 b2.13,t2.13,p2.13,
                 b2.14,t2.14,p2.14,
                 ncol = 3, nrow = 4,labels = c(paste(LETTERS[1:12], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(case4, top = text_grob("Non-compliance caused by non-random processes affecting both groups", face = "bold", size = 12))

#Run simulations for case 2.2: Non-compliance caused by non random process (affect intervention groups)
b2.2<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t2.2<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.2<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
b2.21<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t2.21<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = -0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.21<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case5<-ggarrange(b2.2, t2.2, p2.2, 
                 b2.21, t2.21, p2.21, 
                 ncol = 3, nrow = 2,labels = c(paste(LETTERS[1:6], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(case5, top = text_grob("Non-compliance caused by non-random processes affecting intervention group", face = "bold", size = 12))

#Run simulations for case 2.3: Non-compliance caused by non random process (affect control groups)
b2.3<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t2.3<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = 0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.3<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = 0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
b2.31<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t2.31<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = -0.2, confounder.u.eff.i=1, confounder.u.eff.o=0, nIterations=1000)
p2.31<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000, confounder.eff.o = -0.2,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
case6<-ggarrange(b2.3, t2.3, p2.3, 
                 b2.31, t2.31, p2.31, 
                 ncol = 3, nrow = 2,labels = c(paste(LETTERS[1:6], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(case6, top = text_grob("Non-compliance caused by non-random processes affecting control group", face = "bold", size = 12))

#Run simulations for case 2.4: Non-compliance caused by non random process (affect control groups) with undocumented confounders
b2.4<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 1000, confounder.eff.o = 0.1,  confounder.u.eff.i=1000,confounder.u.eff.o=0.1,nIterations=1000)
t2.4<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000, confounder.eff.o = 0.1, confounder.u.eff.i=1000, confounder.u.eff.o=0.1, nIterations=1000)
p2.4<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000, confounder.eff.o = 0.1,  confounder.u.eff.i=1000,confounder.u.eff.o=0.1,nIterations=1000)
b2.41<-simdata.bias.case4   (n=300, p1=0.4, p0=0.3,       i0=0.01, i1=0.99, confounder.eff.i= 100000000, confounder.eff.o = 0.2,  confounder.u.eff.i=0.01,confounder.u.eff.o=-0.1,nIterations=1000)
t2.41<-simdata.type1.case4  (n=300, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 100000000, confounder.eff.o = 0.2, confounder.u.eff.i=0.01, confounder.u.eff.o=-0.1, nIterations=1000)
p2.41<- simdata.power.case4 (n=300, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 100000000, confounder.eff.o = 0.2,  confounder.u.eff.i=0.01,confounder.u.eff.o=-0.1,nIterations=1000)
case7<-ggarrange(b2.4, t2.4, p2.4, 
                 b2.41, t2.41, p2.41, 
                 ncol = 3, nrow = 2,labels = c(paste(LETTERS[1:6], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(case7, top = text_grob("Non-compliance caused by non-random processes with undocumented confounders", face = "bold", size = 12))


######################################GROUP SEQUENTIAL#####################################
###########################################################################################
simdata.gs<- function(n, p1, p0, i0, i1, NImargin, confounder.eff.o, confounder.eff.i,confounder.u.eff.i,confounder.u.eff.o, nIterations){  
  
  if ((NImargin > 0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong < NI 
    if (nIterations>1) {
      
      interval<-10 #number of data points in simulation 
      
      #make empty vectors for iterations 
      eff.itt1<-eff.pp1<-eff.at1<-eff.iv1<-eff.mpp1<-eff.ps1<- c()
      eff.itt2<-eff.pp2<-eff.at2<-eff.iv2<-eff.mpp2<-eff.ps2<- c()
      eff.itt3<-eff.pp3<-eff.at3<-eff.iv3<-eff.mpp3<-eff.ps3<- c()
      eff.itt4<-eff.pp4<-eff.at4<-eff.iv4<-eff.mpp4<-eff.ps4<- c()
      z.itt1<-z.pp1<-z.at1<-z.iv1<-z.mpp1<-z.ps1<- c()
      z.itt2<-z.pp2<-z.at2<-z.iv2<-z.mpp2<-z.ps2<- c()
      z.itt3<-z.pp3<-z.at3<-z.iv3<-z.mpp3<-z.ps3<-c()
      z.itt4<-z.pp4<-z.at4<-z.iv4<-z.mpp4<-z.ps4<- c()
      .z.itt<-.z.pp<-.z.at<-.z.iv<-c()
      .prop.c<-prop.c<-c()
      
      #simulate and derive treatment effect 
      sim<- function() {
        for(i in 1:interval) {  print(paste ("interval value",i))
          for(l in 1:nIterations) { 
            #simulate data frame 
            
            #RANDOMISATION ratio 1:1
            simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))
            simdata$id<- seq(1,(2*n), by=1) #create participant id 
            
            #CONFOUNDER normal distribution ranging 0-1 
            simdata$confounder<- rep(NA,2*n) 
            simdata$confounder<- rnorm(2*n) 
            simdata$confounder<- rbeta(n=n,shape1=2,shape2=2)
            
            #CONFOUNDER nknown normal distribution ranging 0-1 
            simdata$confounder.u<- rep(NA,2*n) 
            simdata$confounder.u<- rnorm(2*n) 
            simdata$confounder.u<- rbeta(n=n,shape1=2,shape2=2)
            
            #COUNTERFACTUAL OUTCOMES with or without intervention (dependent on confounders and intervention)
            simdata$outcome1 <- rep(NA,2*n)     #outcome with intervention 
            simdata$outcome0 <- rep(NA,2*n)     #outcome without intervention 
            
            prob.outcome.1<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o + (p1-p0) 
            simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
            
            prob.outcome.0<-    p0+ simdata$confounder*confounder.eff.o+simdata$confounder.u*confounder.u.eff.o
            simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
            
            #INTERVENTION dependent on randomisation and confounders by logistric regression - creating increasing compliance
            simdata$intervention <- rep(NA,2*n)
            
            .b0<-logit (i0)
            b0<- seq (.b0,logit(0.01), length.out = interval)
            .b1<-log (confounder.eff.i)
            b1<- seq (.b1,0, length.out = interval)
            .b2<-log (confounder.u.eff.i)
            b2<- seq (.b2,0, length.out = interval)
            .b3<-logit (i1)- logit (i0)
            b3<- seq(.b3,log((0.99/(1-0.99))/(0.01/(1-0.01))), length.out = interval)
            
            if (b1[i]==0) {pi<- simdata$randomisation} else {
              logit.pi<- b0[i]+b1[i]*simdata$confounder+b2[i]*simdata$confounder.u+b3[i]*simdata$randomisation
              pi<-exp(logit.pi)/(1+exp(logit.pi))
            }
            simdata$intervention <- rbinom(2*n,1,prob=pi)
            
            #ACTUAL OUTCOMES depend on intervention
            simdata$outcome<- rep(NA,2*n)
            for (y in 1:(2*n)) {
              if (simdata$intervention[y]==1) {simdata$outcome[y]<-simdata$outcome1[y]} else { 
                simdata$outcome[y]<-simdata$outcome0[y]}
            }
            
            #shuffle rows in simulated data 
            simdata <- simdata[sample(nrow(simdata)),]
            
            #data for each interim analysis 
            s1<-simdata[1:(2*n/4),]
            s2<-simdata[1:(2*n/4*2),]
            s3<-simdata[1:(2*n/4*3),]
            s4<-simdata[1:(2*n),]
            
            #Compute z values at each interim analysis (3 interim and 1 final) 
            
            #######s1
      
            ###estimate treatment effects s1
            ##intention to treat
            pz1.vector<- (filter(s1,randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
            pz0.vector<- (filter(s1,randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
            eff.itt1[l]<- pz1.value-pz0.value            
            
            ##per protocol
            pp1<- rbind(filter(s1, randomisation==1 & intervention==1), (filter(s1, randomisation==0 & intervention==0))) # perprotocol population
            p11.vector<- (filter(s1,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
            p00.vector<- (filter(s1,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
            eff.pp1[l] <- p11.value-p00.value           
            
            ##as treated 
            pd1.vector<- (filter(s1,intervention==1))$outcome;pd1.value<- mean(pd1.vector)
            pd0.vector<- (filter(s1,intervention==0))$outcome;pd0.value<- mean(pd0.vector)
            eff.at1[l] <- pd1.value-pd0.value            
            
            ##linear binomial generalized linear models with inverse probability weights
            E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=s1, na.action=na.exclude) #calculate denominators used in inverse probability weights
            ps<-predict(E.out, type="response")
            sptw<-s1$intervention*mean(s1$intervention)/ps+(1-s1$intervention)*(1-mean(s1$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
            eff.mpp1[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=s1)))$coefficients['intervention']# fit linear binomial model to the weighted data
            
            ##propensity score matching
            m<- glm(intervention~confounder, family = binomial(), data=s1)  
            pscore<- m$fitted.values
            mout<-matchit(intervention~confounder, data=s1, method='nearest')
            match<-Match(Tr=s1$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
            matched<-simdata[unlist(match[c('index.treated','index.control')]),]
            y_trt<-matched$outcome[matched$intervention==1]
            y_con<-matched$outcome[matched$intervention==0]
            diffy<-y_trt-y_con
            eff.ps1[l]<-(t.test(diffy))$estimate
            
            ##iv with 2 stage regression
            ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s1)
            eff.iv1[l]<-ivmodel$coef[2]                 
            
            #######s2
            ###estimate treatment effects s2
            ##intention to treat
            pz1.vector<- (filter(s2, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
            pz0.vector<- (filter(s2, randomisation==0))$outcome; pz0.value<-mean(pz0.vector)
            eff.itt2[l]<- pz1.value-pz0.value             
            
            ##per protocol
            pp2<- rbind(filter(s2, randomisation==1 & intervention==1), (filter(s2, randomisation==0 & intervention==0))) # perprotocol population
            p11.vector<- (filter(s2,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
            p00.vector<- (filter(s2,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
            eff.pp2[l] <- p11.value-p00.value 
            
            ##as treated 
            pd1.vector<- (filter(s2, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
            pd0.vector<- (filter(s2, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
            eff.at2[l] <- pd1.value-pd0.value            
            
            ##linear binomial generalized linear models with inverse probability weights
            E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=s2, na.action=na.exclude) #calculate denominators used in inverse probability weights
            ps<-predict(E.out, type="response")
            sptw<-s2$intervention*mean(s2$intervention)/ps+(1-s2$intervention)*(1-mean(s2$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
            eff.mpp2[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=s2)))$coefficients['intervention']# fit linear binomial model to the weighted data
            
            ##propensity score matching
            m<- glm(intervention~confounder, family = binomial(), data=s2)  
            pscore<- m$fitted.values
            mout<-matchit(intervention~confounder, data=s2, method='nearest')
            match<-Match(Tr=s2$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
            matched<-simdata[unlist(match[c('index.treated','index.control')]),]
            y_trt<-matched$outcome[matched$intervention==1]
            y_con<-matched$outcome[matched$intervention==0]
            diffy<-y_trt-y_con
            eff.ps2[l]<-(t.test(diffy))$estimate
            
            ##iv with 2 stage regression
            ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s2)
            eff.iv2[l]<-ivmodel$coef[2]                 
            
            #######s3
            ###estimate treatment effects s3
            ##intention to treat
            pz1.vector<- (filter(s3, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
            pz0.vector<- (filter(s3, randomisation==0))$outcome; pz0.value<-mean(pz0.vector)
            eff.itt3[l]<- pz1.value-pz0.value             
            
            ##per protocol
            pp3<- rbind(filter(s3, randomisation==1 & intervention==1), (filter(s3, randomisation==0 & intervention==0))) # perprotocol population
            p11.vector<- (filter(s3,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
            p00.vector<- (filter(s3,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
            eff.pp3[l] <- p11.value-p00.value            
            
            ##as treated
            pd1.vector<- (filter(s3, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
            pd0.vector<- (filter(s3, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
            eff.at3[l] <- pd1.value-pd0.value             

            ##linear binomial generalized linear models with inverse probability weights
            E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=s3, na.action=na.exclude) #calculate denominators used in inverse probability weights
            ps<-predict(E.out, type="response")
            sptw<-s3$intervention*mean(s3$intervention)/ps+(1-s3$intervention)*(1-mean(s3$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
            eff.mpp3[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=s3)))$coefficients['intervention']# fit linear binomial model to the weighted data
            
            ##propensity score matching 
            m<- glm(intervention~confounder, family = binomial(), data=s3) 
            pscore<- m$fitted.values
            mout<-matchit(intervention~confounder, data=s3, method='nearest')
            match<-Match(Tr=s3$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
            matched<-simdata[unlist(match[c('index.treated','index.control')]),]
            y_trt<-matched$outcome[matched$intervention==1]
            y_con<-matched$outcome[matched$intervention==0]
            diffy<-y_trt-y_con
            eff.ps3[l]<-(t.test(diffy))$estimate
            
            ##iv with 2 stage regression
            ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s3)
            eff.iv3[l]<-ivmodel$coef[2]                 
            
            #######s4
            ###estimate treatment effects s4
            ##intention to treat
            pz1.vector<- (filter(s4, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
            pz0.vector<- (filter(s4, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
            eff.itt4[l]<- pz1.value-pz0.value           
            
            ##per protocol
            pp4<- rbind(filter(s4, randomisation==1 & intervention==1), (filter(s4, randomisation==0 & intervention==0))) # perprotocol population
            .prop.c[l]<- nrow(pp4)/(2*n) #proportion of patients who were compliant to protocol
            p11.vector<- (filter(s4,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
            p00.vector<- (filter(s4,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
            eff.pp4[l] <- p11.value-p00.value            
            
            ##as treated
            pd1.vector<- (filter(s4, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
            pd0.vector<- (filter(s4, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
            eff.at4[l] <- pd1.value-pd0.value             
            
            ##linear binomial generalized linear models with inverse probability weights
            E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=s4, na.action=na.exclude) #calculate denominators used in inverse probability weights
            ps<-predict(E.out, type="response")
            sptw<-s4$intervention*mean(s4$intervention)/ps+(1-s4$intervention)*(1-mean(s4$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
            eff.mpp4[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=s4)))$coefficients['intervention']# fit linear binomial model to the weighted data
            
            ##propensity score matching 
            m<- glm(intervention~confounder, family = binomial(), data=s4) 
            pscore<- m$fitted.values
            mout<-matchit(intervention~confounder, data=s4, method='nearest')
            match<-Match(Tr=s4$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
            matched<-simdata[unlist(match[c('index.treated','index.control')]),]
            y_trt<-matched$outcome[matched$intervention==1]
            y_con<-matched$outcome[matched$intervention==0]
            diffy<-y_trt-y_con
            eff.ps4[l]<-(t.test(diffy))$estimate
            
            ##iv with 2 stage regression
            ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s4)
            eff.iv4[l]<-ivmodel$coef[2]                 
            
          }
          
          ###sds of effects s1
          sd.eff.itt1<-sd(eff.itt1)
          sd.eff.pp1<- sd(eff.pp1)
          sd.eff.at1<- sd(eff.at1)
          sd.eff.ps1<- sd(eff.ps1)
          sd.eff.mpp1<-sd(eff.mpp1)
          sd.eff.iv1<- sd(eff.iv1)
          
          ###mean z values s1
          z.itt1[i]<-(mean(eff.itt1)-NImargin)/sd.eff.itt1
          z.iv1[i]<- (mean(eff.iv1)- NImargin)/sd.eff.iv1
          z.pp1[i]<- (mean(eff.pp1)- NImargin)/sd.eff.pp1
          z.mpp1[i]<-(mean(eff.mpp1)-NImargin)/sd.eff.mpp1
          z.at1[i]<- (mean(eff.at1)- NImargin)/sd.eff.at1
          z.ps1[i]<- (mean(eff.ps1)- NImargin)/sd.eff.ps1
          
          ###sds of effects s2
          sd.eff.itt2<-sd(eff.itt2)
          sd.eff.pp2<- sd(eff.pp2)
          sd.eff.at2<- sd(eff.at2)
          sd.eff.ps2<- sd(eff.ps2)
          sd.eff.mpp2<-sd(eff.mpp2)
          sd.eff.iv2<- sd(eff.iv2)
          
          ###mean z values s2
          z.itt2[i]<-(mean(eff.itt2)-NImargin)/sd.eff.itt2
          z.iv2[i]<- (mean(eff.iv2)- NImargin)/sd.eff.iv2
          z.pp2[i]<- (mean(eff.pp2)- NImargin)/sd.eff.pp2
          z.mpp2[i]<-(mean(eff.mpp2)-NImargin)/sd.eff.mpp2
          z.at2[i]<- (mean(eff.at2)- NImargin)/sd.eff.at2
          z.ps2[i]<- (mean(eff.ps2)- NImargin)/sd.eff.ps2
          
          ###sds of effects s3
          sd.eff.itt3<-sd(eff.itt3)
          sd.eff.pp3<- sd(eff.pp3)
          sd.eff.at3<- sd(eff.at3)
          sd.eff.ps3<- sd(eff.ps3)
          sd.eff.mpp3<-sd(eff.mpp3)
          sd.eff.iv3<- sd(eff.iv3)
          
          ###mean z values s3
          z.itt3[i]<-(mean(eff.itt3)-NImargin)/sd.eff.itt3
          z.iv3[i]<- (mean(eff.iv3)- NImargin)/sd.eff.iv3
          z.pp3[i]<- (mean(eff.pp3)- NImargin)/sd.eff.pp3
          z.mpp3[i]<-(mean(eff.mpp3)-NImargin)/sd.eff.mpp3
          z.at3[i]<- (mean(eff.at3)- NImargin)/sd.eff.at3
          z.ps3[i]<- (mean(eff.ps3)- NImargin)/sd.eff.ps3
          
          ###sds of effects s4
          sd.eff.itt4<-sd(eff.itt4)
          sd.eff.pp4<- sd(eff.pp4)
          sd.eff.at4<- sd(eff.at4)
          sd.eff.ps4<- sd(eff.ps4)
          sd.eff.mpp4<-sd(eff.mpp4)
          sd.eff.iv4<- sd(eff.iv4)
          
          ###mean z values s4
          z.itt4[i]<-(mean(eff.itt4)-NImargin)/sd.eff.itt4
          z.iv4[i]<- (mean(eff.iv4)- NImargin)/sd.eff.iv4
          z.pp4[i]<- (mean(eff.pp4)- NImargin)/sd.eff.pp4
          z.mpp4[i]<-(mean(eff.mpp4)-NImargin)/sd.eff.mpp4
          z.at4[i]<- (mean(eff.at4)- NImargin)/sd.eff.at4
          z.ps4[i]<- (mean(eff.ps4)- NImargin)/sd.eff.ps4
          
          .z.itt<-list(z.itt1, z.itt2,z.itt3,z.itt4)
          .z.at<- list(z.at1, z.at2,z.at3,z.at4)
          .z.pp<- list(z.pp1,z.pp2,z.pp3,z.pp4)
          .z.mpp<-list(z.mpp1,z.mpp2,z.mpp3,z.mpp4)
          .z.ps<- list(z.ps1,z.ps2,z.ps3,z.ps4)
          .z.iv<- list(z.iv1,z.iv2,z.iv3,z.iv4)
          
          #mean proportion of patients who were compliant to protocol 
          prop.c [i]<- mean(.prop.c)
        } 
        
        #group sequential boundaries at each analysis with 80% power 
        z<-gsDesign(k=4, n.fix=2*n, delta0 = NImargin, n.I=c(nrow(s1),nrow(s2) ,nrow(s3),nrow(s4)), maxn.IPlan = 2*n, beta=0.2) 
        z<-z$upper$bound #4 boundaries 
        
        z.itt<-data.frame(.z.itt)
        z.at<- data.frame(.z.at)
        z.pp<- data.frame(.z.pp)
        z.mpp<-data.frame(.z.mpp)
        z.ps<- data.frame(.z.ps)
        z.iv<- data.frame(.z.iv)
        z.c<-  abs(data.frame(cbind(prop.c, z.itt, z.at, z.pp, z.mpp,z.ps,z.iv, 
                              rep(z[1],length(interval)),
                              rep(z[2],length(interval)),
                              rep(z[3],length(interval)),
                              rep(z[4],length(interval)))))
        
        names(z.c)<-c('prop.c','itt1','itt2','itt3','itt4',
                      'at1','at2','at3','at4',
                      'pp1','pp2','pp3','pp4',
                      'mpp1','mpp2','mpp3','mpp4',
                      'ps1','ps2','ps3','ps4',
                      'iv1','iv2','iv3','iv4',
                      'z1','z2','z3','z4')
        
        return(z.c)
      }
      
      z.c<- sim()
      
      plot.z.itt<- ggplot(z.c, aes(z.c[1]))+ 
        geom_point(aes(y = itt1, colour="Z value 1" , alpha=0.4))+
        geom_point(aes(y = itt2, colour="Z value 2",  alpha=0.4))+
        geom_point(aes(y = itt3, colour="Z value 3",  alpha=0.4))+
        geom_point(aes(y = itt4, colour="Z value 4",  alpha=0.4))+
        geom_smooth(aes(y= itt1, colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=itt2, colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=itt3, colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=itt4, colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y = z1, colour="Z value 1"),linetype="dotted")+
        geom_line(aes(y = z2, colour="Z value 2"),linetype="dotted")+
        geom_line(aes(y = z3, colour="Z value 3"),linetype="dotted")+
        geom_line(aes(y = z4, colour="Z value 4"),linetype="dotted")+
        xlab("Proportion of compliant participants")+
        ylab("Z value")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.55, 1))+
        scale_y_continuous(limits=c(0.5, 4))
      
      plot.z.at<- ggplot(z.c, aes(z.c[1]))+ 
        geom_point(aes(y = at1, colour="Z value 1" , alpha=0.4))+
        geom_point(aes(y = at2, colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = at3, colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = at4, colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y=at1, colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=at2, colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=at3, colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=at4, colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y = z1, colour="Z value 1"),linetype="dotted")+
        geom_line(aes(y = z2, colour="Z value 2"),linetype="dotted")+
        geom_line(aes(y = z3, colour="Z value 3"),linetype="dotted")+
        geom_line(aes(y = z4, colour="Z value 4"),linetype="dotted")+
        xlab("Proportion of compliant participants")+
        ylab("Z value")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.55, 1))+
        scale_y_continuous(limits=c(0.5, 4))
      
      plot.z.pp<- ggplot(z.c, aes(z.c[1]))+ 
        geom_point(aes(y = pp1, colour="Z value 1", alpha=0.4 ))+
        geom_point(aes(y = pp2, colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = pp3, colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = pp4, colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y=pp1, colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=pp2, colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=pp3, colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=pp4, colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y = z1, colour="Z value 1"),linetype="dotted")+
        geom_line(aes(y = z2, colour="Z value 2"),linetype="dotted")+
        geom_line(aes(y = z3, colour="Z value 3"),linetype="dotted")+
        geom_line(aes(y = z4, colour="Z value 4"),linetype="dotted")+
        xlab("Proportion of compliant participants")+
        ylab("Z value")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.55, 1))+
        scale_y_continuous(limits=c(0.5, 4))
      
      plot.z.mpp<- ggplot(z.c, aes(z.c[1]))+ 
        geom_point(aes(y = mpp1, colour="Z value 1" , alpha=0.4))+
        geom_point(aes(y = mpp2, colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = mpp3, colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = mpp4, colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y=mpp1, colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=mpp2, colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=mpp3, colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=mpp4, colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y = z1, colour="Z value 1"),linetype="dotted")+
        geom_line(aes(y = z2, colour="Z value 2"),linetype="dotted")+
        geom_line(aes(y = z3, colour="Z value 3"),linetype="dotted")+
        geom_line(aes(y = z4, colour="Z value 4"),linetype="dotted")+
        xlab("Proportion of compliant participants")+
        ylab("Z value")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.55, 1))+
        scale_y_continuous(limits=c(0.5, 4))
      
      plot.z.ps<- ggplot(z.c, aes(z.c[1]))+ 
        geom_point(aes(y = ps1, colour="Z value 1" , alpha=0.4))+
        geom_point(aes(y = ps2, colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = ps3, colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = ps4, colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y=ps1, colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=ps2, colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=ps3, colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=ps4, colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y = z1, colour="Z value 1"),linetype="dotted")+
        geom_line(aes(y = z2, colour="Z value 2"),linetype="dotted")+
        geom_line(aes(y = z3, colour="Z value 3"),linetype="dotted")+
        geom_line(aes(y = z4, colour="Z value 4"),linetype="dotted")+
        xlab("Proportion of compliant participants")+
        ylab("Z value")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.55, 1))+
        scale_y_continuous(limits=c(0.5, 4))
      
      plot.z.iv<- ggplot(z.c, aes(z.c[1]))+ 
        geom_point(aes(y = iv1, colour="Z value 1", alpha=0.4 ))+
        geom_point(aes(y = iv2, colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = iv3, colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = iv4, colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y=iv1, colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=iv2, colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=iv3, colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=iv4, colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_line(aes(y = z1, colour="Z value 1"),linetype="dotted")+
        geom_line(aes(y = z2, colour="Z value 2"),linetype="dotted")+
        geom_line(aes(y = z3, colour="Z value 3"),linetype="dotted")+
        geom_line(aes(y = z4, colour="Z value 4"),linetype="dotted")+
        xlab("Proportion of compliant participants")+
        ylab("Z value")+
        theme_bw()+
        scale_colour_manual(values=cbPalette)+
        theme(legend.title=element_blank())+
        scale_x_continuous(limits=c(0.6, 1))+
        scale_y_continuous(limits=c(0.5, 7))
      
      plot<- ggarrange(plot.z.itt, plot.z.at, plot.z.pp, plot.z.mpp,plot.z.ps,plot.z.iv, ncol = 3, nrow = 2,
                       labels = c("Intention to treat", "As treated","Per protocol","Modified PP","Propensity score","Instrumental variable"), common.legend = TRUE, legend = "bottom")
      
      return(plot)
    }else (print ("Number of iterations must be more than 1."))
  } else (print ("NI margin must be positive, and p1-p0 (in terms of unfavorable outcomes) must be less than NI margin in this simulation"))
}

gs1<-simdata.gs(n=200, p1=0.4, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000, NImargin=0.12)
annotate_figure(gs1, top = text_grob("Figure 1. Varying effect of unknown confounder", face = "bold", size = 12))

gs2<-simdata.gs(n=200, p1=0.4, p0=0.3, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000, NImargin=0.12) 
annotate_figure(gs2, top = text_grob("Figure 2. Varying effect of unknown confounder", face = "bold", size = 12))

gs3<-simdata.gs(n=200, p1=0.3, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000, NImargin=0.12) 
annotate_figure(gs3, top = text_grob("Figure 3. Varying effect of unknown confounder", face = "bold", size = 12))

gs4<-simdata.gs(n=200, p1=0.4, p0=0.3, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000, NImargin=0.12) 
annotate_figure(gs4, top = text_grob("Figure 4. Varying effect of unknown confounder", face = "bold", size = 12))

gs5<-simdata.gs(n=200, p1=0.4, p0=0.3, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000, NImargin=0.12) 
annotate_figure(gs5, top = text_grob("Figure 5. Varying effect of unknown confounder", face = "bold", size = 12))

gs6<-simdata.gs(n=200, p1=0.4, p0=0.3, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000, NImargin=0.12) 
annotate_figure(gs6, top = text_grob("Figure 6. Varying effect of unknown confounder", face = "bold", size = 12))

gs7<-simdata.gs(n=200, p1=0.4, p0=0.3, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000, NImargin=0.12) 
annotate_figure(gs7, top = text_grob("Figure 7. Varying effect of unknown confounder", face = "bold", size = 12))