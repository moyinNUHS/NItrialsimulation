########################################################################
###################Simulation for REGARD-VAP analysis###################
########################################################################

setwd("/Users/moyin/Desktop/VAP studd/Causal inference simulation") #set working directory 
rm(list=ls()) # Clean working environment 
library(Hmisc); library(rms); library(gsDesign);library(ivpack); library(ggplot2); 
library(data.table); library(dplyr); library(plotly); library(ggpubr) # Required libraries 
library(Matching); library(tableone); library(MatchIt); library(geepack)

set.seed(1234)

##########BIAS##########
########################
simdata.bias<- function(n, p0, p1, confounder.eff.o, nIterations){  
  # n: number of participants per group; 
  # p0: mortality and recurrence rate in intervention 0; 
  # p1: mortality and recurrence rate in intervention 1;
  # confounder.eff.o: effect of counfounder on outcome 0:1
  # nIterations: number of iterations 
  
  alpha=0.025
  odd1<- seq((1/1)/(1/1),(10/1)/(1/10), by=9)
  
  #make up vectors for simulations 
  eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
  .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.mpp<-.eff.ps<-c() 
  prop.c<-.prop.c<-c()
  
  #simulate and derive treatment effect 
  sim<- function() { 
    for(i in 1:length(odd1)) { print(i)
      for(l in 1:nIterations) { 
        #simulate data frame 
        
        #RANDOMISATION ratio 1:1
        simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))     #half randomised to 1, half randomised to 0 
        simdata$id<- seq(1,(2*n), by=1) #create participant id 
        
        #confounder normal distribution ranging 0-1 
        simdata$confounder<- rep(NA,2*n) 
        simdata$confounder<- rnorm((2*n),6, 2) # assume counfounder to be SOFA(severity), normally distributed 
        simdata$confounder [which(simdata$confounder<0)] <- -1*(simdata$confounder[which(simdata$confounder<0)]) #make sure all sofa scores are positive 
        simdata$confounder<- simdata$confounder/max(simdata$confounder) # scale from 0-1
        
        #INTERVENTION dependent on randomisation and confounders
        simdata$intervention <- rep(NA,2*n)
        
        odd0<- 1/5                                            #Number of intervention=1 when C=0 R=0/ Number of intervention=0 when C=0 R=0
        odd1<- seq((1/1)/(1/1),(10/1)/(1/10), by=9)           #(Number of intervention=1 when C=1 / Number of intervention=0 when when C=1) OVER
                                                              # (Number of of intervention =1 when C=0 / Number of of intervention =0 when C=0 ) (controlled for R)
                                                              # if confounder has no negative effect on compliance - (1/1)/(1/1)
                                                              # if confounder has strong negative effect on compliance - (1/10)/(10/1)
        odd2<- (50/1)/(1/50)                                  #(Number of intervention=1 when R=1 / Number of intervention=0 when when R=1) OVER
                                                              # (Number of of intervention =1 when R=0 / Number of of intervention =0 when R=0 ) (controlled for C)
                                                              # randomisation is assumed to have more effect on compliance than confounder 
                                                              
        b0<-log(odd0)
        b1<-log (odd1)
        b2<-log (odd2)
        
        logit.pi<- b0+b1[i]*simdata$confounder+b2*simdata$randomisation
        pi<-exp(logit.pi)/(1+exp(logit.pi))
        
        simdata$intervention <- rbinom(2*n,1,prob=pi)
        
        #OUTCOME dependent on confounders and intervention
        simdata$outcome <- rep(NA,2*n)
        
        odd0<- 1/5                                            # Number of outcome=1 when C=0 I=0/ Number of outcome=0 when C=0 I=0
        odd1<- confounder.eff.o                               # (Number of outcome=1 when C=1 / Number of outcome=0 when when C=1) OVER
                                                              # (Number of of outcome =1 when C=0 / Number of of outcome =0 when C=0 ) (controlled for R)
                                                              # Strength of confounder on outcome
        odd2<- (p1/(1-p1))/(p0/(1-p0))                        # (Number of outcome=1 when I=1 / Number of outcome=0 when when I=1) OVER
                                                              # (Number of of outcome =1 when I=0 / Number of of outcome=0 when I=0 ) (controlled for C)
        
        b0<-log(odd0)
        b1<-log (odd1)
        b2<-log (odd2)
        
        logit.pi<- b0+b1*simdata$confounder+b2*simdata$intervention
        pi<-exp(logit.pi)/(1+exp(logit.pi))
        
        simdata$outcome <- rbinom(2*n,1,prob=pi)
        
        #generate proportions from simulated data  
        p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome)
        p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
  
        pz1.value<- mean((filter(simdata, randomisation==1))$outcome) 
        pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
        
        pd1.value<- mean((filter(simdata, intervention==1))$outcome)
        pd0.value<- mean((filter(simdata, intervention==0))$outcome)
        
        #treatment effect 
        .eff.itt[l]<- pz1.value-pz0.value              #intention to treat 
        .eff.pp[l] <- p11.value-p00.value              #per protocol 
        .eff.at[l] <- pd1.value-pd0.value              # as treated 
        
        pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
        # Linear binomial generalized linear models with inverse probability weights
        E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
        ps<-predict(E.out, type="response")
        sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
        .eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
        
        m<- glm(intervention~confounder, family = binomial(), data=simdata) #propensity score matching 
        pscore<- m$fitted.values
        mout<-matchit(intervention~confounder, data=simdata, method='nearest')
        match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
        matched<-simdata[unlist(match[c('index.treated','index.control')]),]
        y_trt<-matched$outcome[matched$intervention==1]
        y_con<-matched$outcome[matched$intervention==0]
        diffy<-y_trt-y_con
        .eff.ps[l]<-(t.test(diffy))$estimate
        
        ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) #iv with 2 stage regression
        .eff.iv[l]<-ivmodel$coef[2] 
        
        #proportion of patients who were compliant to protocol 
        .prop.c[l]<- nrow(pp)/(2*n)
        
        }
      #mean of iterated data 
      eff.itt[i]<- mean(.eff.itt, na.rm=TRUE)
      eff.pp[i] <- mean(.eff.pp, na.rm=TRUE) 
      eff.at[i] <- mean(.eff.at, na.rm=TRUE)
      eff.iv[i]<- mean(.eff.iv, na.rm = TRUE)
      eff.mpp[i]<- mean(.eff.mpp, na.rm = TRUE)
      eff.ps[i]<- mean(.eff.ps, na.rm = TRUE)
        
      #proportion of patients who were compliant to protocol 
      prop.c [i]<- mean(.prop.c)
      
      }
    return(data.frame(prop.c, eff.iv, eff.itt, eff.pp, eff.at, eff.mpp, eff.ps))
  }
  
  sim.df<-sim()
  
  bias.plot <- ggplot(sim.df, aes(x=sim.df[1]))+ 
    geom_point(aes(y=sim.df[2], colour="IV effect")) + 
    geom_point(aes(y=sim.df[3], colour="ITT effect")) + 
    geom_point(aes(y=sim.df[4], colour="PP effect")) + 
    geom_point(aes(y=sim.df[5], colour="AT effect")) + 
    geom_point(aes(y=sim.df[6], colour="Modified PP effect")) + 
    geom_point(aes(y=sim.df[7], colour="Propensity score effect")) + 
    geom_smooth(aes(y=sim.df[2], colour="IV effect"),se=FALSE)+
    geom_smooth(aes(y=sim.df[3], colour="ITT effect"),se=FALSE)+
    geom_smooth(aes(y=sim.df[4], colour="PP effect"),se=FALSE)+
    geom_smooth(aes(y=sim.df[5], colour="AT effect"),se=FALSE)+
    geom_smooth(aes(y=sim.df[6], colour="Modified PP effect"),se=FALSE)+
    geom_smooth(aes(y=sim.df[7], colour="Propensity score effect"),se=FALSE)+
    geom_line(aes(y=p1-p0, colour='True effect'),linetype="dotted") +
    xlab("Proportion of compliant participants")+
    ylab("Effect")
  
  return(bias.plot)
} 

bias1<-simdata.bias(n=230,p1=0.4, p0=0.4, confounder.eff.o = (8/1)/(1/5),  nIterations=1000)
bias2<-simdata.bias(n=230,p1=0.4, p0=0.4, confounder.eff.o = (4/1)/(1/5),  nIterations=1000)
bias3<-simdata.bias(n=230,p1=0.4, p0=0.4, confounder.eff.o = (2/1)/(1/5),  nIterations=1000)

ggarrange(bias1, bias2, bias3, ncol = 3, nrow = 1)

bias4<-simdata.bias(n=230, p1=0.6, p0=0.4,  confounder.eff.o = (4/1)/(1/5),nIterations=1000) 
bias5<-simdata.bias(n=230, p1=0.5, p0=0.4,  confounder.eff.o = (4/1)/(1/5),nIterations=1000)
bias6<-bias2
bias7<-simdata.bias(n=230, p1=0.3, p0=0.4,  confounder.eff.o = (4/1)/(1/5),nIterations=1000) 
bias8<-simdata.bias(n=230, p1=0.2, p0=0.4,  confounder.eff.o = (4/1)/(1/5),nIterations=1000) 

ggarrange(bias4, bias5, bias6,  bias7,bias8, 
          ncol = 3, nrow = 2)

bias9<- bias5
bias10<-simdata.bias(n=300,p1=0.5, p0=0.4,  confounder.eff.o = (4/1)/(1/5),nIterations=1000)
bias11<-simdata.bias(n=400,p1=0.5, p0=0.4,  confounder.eff.o = (4/1)/(1/5),nIterations=1000)

ggarrange(bias9, bias10, bias11,  
          ncol = 3, nrow = 1)

######TYPE 1 ERROR######
########################
simdata.type1<- function(n, p1, NImargin, confounder.eff.o, nIterations){  
  # n: number of participants per group; 
  # p0: mortality and recurrence rate  in intervention 0; 
  # p1: mortality and recurrence rate in intervention 1;
  # NImargin: non inferiority margin 
  # confounder.eff.o: effect of counfounder on outcome 0:1
  # nIterations: number of iterations 
  
  #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
  lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE) 
  lower.value<-lower$conf.int [1] #95% CI of NI margin 
  
  #build simulated data based on null hypothesis pshort-plong >= +NI (mortality or recurrences)
  p0=p1-NImargin 
  
  # for simulations to run over the length of odd1
  odd1<- seq((1/1)/(1/1),(10/1)/(1/10), by=9)
  
  #make up vectors for simulations 
  eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()                   
  type1.error.at<-type1.error.itt<-type1.error.iv<-type1.error.mpp<-type1.error.pp<-type1.error.ps<-c() 
  prop.c<-.prop.c<-c()
  
  if (NImargin>0) {
    
    #simulate data
    sim<- function() { 
      for(i in 1:length(odd1)) { print(i)
        for(l in 1:nIterations) { 
          #simulate data frame 
        
          #RANDOMISATION ratio 1:1
          simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))     #half randomised to 1, half randomised to 0 
          simdata$id<- seq(1, (2*n), by=1)
          
          #confounder normal distribution ranging 0-1 
          simdata$confounder<- rep(NA,2*n) 
          simdata$confounder<- rnorm((2*n),6, 2) # assume counfounder to be SOFA(severity), normally distributed 
          simdata$confounder[which(simdata$confounder<0)] <- -1*(simdata$confounder[which(simdata$confounder<0)]) #make sure all sofa scores are positive 
          simdata$confounder<- simdata$confounder/max(simdata$confounder) # scale from 0-1
          
          #INTERVENTION dependent on randomisation and confounders
          simdata$intervention <- rep(NA,2*n)
          
          odd0<- 1/5                                            #Number of intervention=1 when C=0 R=0/ Number of intervention=0 when C=0 R=0
          odd1<- seq((1/1)/(1/1),(10/1)/(1/10), by=9)           #(Number of intervention=1 when when C=1 / Number of intervention=0 when when C=1) OVER
                                                                # (Number of of intervention =1 when C=0 / Number of of intervention =0 when C=0 ) (controlled for R)
                                                                # if confounder has no negative effect on compliance - (1/1)/(1/1)
                                                                # if confounder has strong negative effect on compliance - (1/10)/(10/1)
          odd2<- (50/1)/(1/50)                                  #(Number of intervention=1 when when R=1 / Number of intervention=0 when when R=1) OVER
                                                                # (Number of of intervention =1 when R=0 / Number of of intervention =0 when R=0 ) (controlled for C)
                                                                # randomisation is assumed to have more effect on compliance than confounder 
          
          b0<-log(odd0)
          b1<-log (odd1)
          b2<-log (odd2)
          
          logit.pi<- b0+b1[i]*simdata$confounder+b2*simdata$randomisation
          pi<-exp(logit.pi)/(1+exp(logit.pi))
          
          simdata$intervention <- rbinom(2*n,1,prob=pi)
          
          #OUTCOME dependent on confounders and intervention
          simdata$outcome <- rep(NA,2*n)
          
          odd0<- 1/5                                            # Number of outcome=1 when C=0 I=0/ Number of outcome=0 when C=0 I=0
          odd1<- confounder.eff.o                               # (Number of outcome=1 when C=1 / Number of outcome=0 when when C=1) OVER
                                                                # (Number of of outcome =1 when C=0 / Number of of outcome =0 when C=0 ) (controlled for R)
                                                                # Strength of confounder on outcome
          odd2<- (p1/(1-p1))/(p0/(1-p0))                        # (Number of outcome=1 when I=1 / Number of outcome=0 when when I=1) OVER
                                                                # (Number of of outcome =1 when I=0 / Number of of outcome=0 when I=0 ) (controlled for C)
          
          b0<-log(odd0)
          b1<-log (odd1)
          b2<-log (odd2)
          
          logit.pi<- b0+b1*simdata$confounder+b2*simdata$intervention
          pi<-exp(logit.pi)/(1+exp(logit.pi))
          
          simdata$outcome <- rbinom(2*n,1,prob=pi)
          
          #generate proportions in simulated data 
          p00.vector<- (filter(simdata,intervention==0 & randomisation==0))$outcome; p00.value<-mean(p00.vector)
          p11.vector<- (filter(simdata,intervention==1 & randomisation==1))$outcome; p11.value<-mean(p11.vector)
          
          pz1.vector<- (filter(simdata, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
          pz0.vector<- (filter(simdata, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
          
          pd1.vector<- (filter(simdata, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
          pd0.vector<- (filter(simdata, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
          
          #estimate treatment effects 
          eff.itt[l]<- pz1.value-pz0.value            #intention to treat 
          eff.pp[l] <- p11.value-p00.value            #per protocol
          eff.at[l] <- pd1.value-pd0.value            #as treated 
          
          pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
          # Linear binomial generalized linear models with inverse probability weights
          E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
          ps<-predict(E.out, type="response")
          sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
          eff.mpp[l]<-(geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp))$coefficients['intervention']# fit linear binomial model to the weighted data
          
          m<- glm(intervention~confounder, family = binomial(), data=simdata) #propensity score matching 
          pscore<- m$fitted.values
          mout<-matchit(intervention~confounder, data=simdata, method='nearest')
          match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
          matched<-simdata[unlist(match[c('index.treated','index.control')]),]
          y_trt<-matched$outcome[matched$intervention==1]
          y_con<-matched$outcome[matched$intervention==0]
          diffy<-y_trt-y_con
          eff.ps[l]<-(t.test(diffy))$estimate
          
          ivmodel=ivreg(simdata$outcome ~ simdata$intervention, ~ simdata$randomisation, x=TRUE, data=simdata)
          eff.iv[l]<-ivmodel$coef[2]                  #iv with 2 stage regression
          
          }
        
          #proportion of patients who were compliant to protocol 
          .prop.c[l]<- nrow(pp)/(2*n)
          
          #type 1 error
          type1.error.itt[i]<-mean(eff.itt<lower.value)
          type1.error.pp[i]<- mean(eff.pp<lower.value)
          type1.error.at[i]<- mean(eff.at<lower.value)
          type1.error.mpp[i]<- mean(eff.mpp<lower.value)
          type1.error.iv[i]<- mean(eff.iv<lower.value)
          type1.error.ps[i]<- mean(eff.ps<lower.value)
          
          #proportion of patients who were compliant to protocol 
          prop.c [i]<- mean(.prop.c[l])
          
          }
      return(data.frame(prop.c, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.mpp, type1.error.ps))
      }
    sim.df<-sim()
  
    type1.plot <- ggplot(sim.df, aes(sim.df[1]))+ 
    geom_point(aes(y=sim.df[2], colour="IV Type 1 error")) + 
    geom_point(aes(y=sim.df[3], colour="ITT Type 1 error")) + 
    geom_point(aes(y=sim.df[4], colour="PP Type 1 error")) + 
    geom_point(aes(y=sim.df[5], colour="AT Type 1 error")) + 
    geom_point(aes(y=sim.df[6], colour="Modified PP Type 1 error")) +
    geom_point(aes(y=sim.df[7], colour="Propensity score Type 1 error")) + 
    geom_smooth(aes(y=sim.df[2], colour="IV Type 1 error"),se=FALSE)+
    geom_smooth(aes(y=sim.df[3], colour="ITT Type 1 error"),se=FALSE)+
    geom_smooth(aes(y=sim.df[4], colour="PP Type 1 error"),se=FALSE)+
    geom_smooth(aes(y=sim.df[5], colour="AT Type 1 error"),se=FALSE)+
    geom_smooth(aes(y=sim.df[6], colour="Modified PP Type 1 error"),se=FALSE)+
    geom_smooth(aes(y=sim.df[7], colour="Propensity score Type 1 error"),se=FALSE)+
    geom_line(aes(y=0.025, colour='True Type 1 error'),linetype="dotted") +
    xlab("Proportion of compliant participants")+
    ylab("Type 1 error")
    
    return(type1.plot)
    
    } else print("NI margin should be positive in this stimulation")
}

t1<-simdata.type1(n=230,p1=0.4, confounder.eff.o=(8/1)/(1/5), NImargin=0.12, nIterations=1000)
t2<-simdata.type1(n=230,p1=0.4, confounder.eff.o=(4/1)/(1/5), NImargin=0.12, nIterations=1000)
t3<-simdata.type1(n=230,p1=0.4, confounder.eff.o=(2/1)/(1/5), NImargin=0.12, nIterations=1000)

ggarrange(t1, t2, t3, ncol = 3, nrow = 1)

t5<-simdata.type1(n=230,p1=0.2,confounder.eff.o=(4/1)/(1/5), NImargin=0.12,nIterations=1000)
t6<-t2
t7<-simdata.type1(n=230,p1=0.6,confounder.eff.o=(4/1)/(1/5), NImargin=0.12,nIterations=1000)

ggarrange(t5, t6, t7, ncol = 3, nrow = 1)

t8<-t2
t9<- simdata.type1(n=400,p1=0.4,confounder.eff.o=(4/1)/(1/5), NImargin=0.12,nIterations=1000)
t10<-simdata.type1(n=600,p1=0.4,confounder.eff.o=(4/1)/(1/5), NImargin=0.12,nIterations=1000)

ggarrange(t8, t9, t10, ncol = 3, nrow = 1)

##########POWER#########
########################
simdata.power<- function(n, p1, p0, confounder.eff.o, NImargin, nIterations){  
  # n: number of participants per group; 
  # p0: mortality and recurrences rate in intervention 0; 
  # p1: mortality and recurrences rate in intervention 1;
  # confounder.eff.o: effect of counfounder on outcome 0:1
  # NImargin: non inferiority margin (negative);
  # nIterations: number of iterations 
  
  #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
  lower<-prop.test(x=c(p1*n,(p1-NImargin)*n),n=c(n,n), correct=FALSE)
  lower.value<-lower$conf.int [1]
  
  if ((NImargin>0) & (p1-p0 < NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
    alpha=0.025
    odd1<- seq((1/1)/(1/1),(10/1)/(1/10), by=9)
    
    #make up vectors for simulations
    power.iv<-power.itt<- power.at<-power.pp<-power.mpp<-power.ps<-c() 
    eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.mpp<-eff.ps<-c()    
    prop.c<-.prop.c<-c()
      
    #simulate and derive treatment effect 
    sim<- function() { 
      for(i in 1:length(odd1)) { print(i)
        for(l in 1:nIterations) { 
          #simulate data frame 
            
          #RANDOMISATION ratio 1:1
          simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))     #half randomised to 1, half randomised to 0 
          simdata$id<- seq(1,(2*n), by=1)
              
          #confounder normal distribution ranging 0-1 
          simdata$confounder<- rep(NA,2*n) 
          simdata$confounder<- rnorm((2*n),6, 2) # assume counfounder to be SOFA(severity), normally distributed 
          simdata$confounder [which(simdata$confounder<0)] <- -1*(simdata$confounder[which(simdata$confounder<0)]) #make sure all sofa scores are positive 
          simdata$confounder<- simdata$confounder/max(simdata$confounder) # scale from 0-1
            
          #INTERVENTION dependent on randomisation and confounders
          simdata$intervention <- rep(NA,2*n)
            
          odd0<- 1/5                                            #Number of intervention=1 when C=0 R=0/ Number of intervention=0 when C=0 R=0
          odd1<- seq((1/1)/(1/1),(10/1)/(1/10), by=9)           #(Number of intervention=1 when when C=1 / Number of intervention=0 when when C=1) OVER
                                                                  # (Number of of intervention =1 when C=0 / Number of of intervention =0 when C=0 ) (controlled for R)
                                                                  # if confounder has no negative effect on compliance - (1/1)/(1/1)
                                                                  # if confounder has strong negative effect on compliance - (1/10)/(10/1)
          odd2<- (50/1)/(1/50)                                  #(Number of intervention=1 when when R=1 / Number of intervention=0 when when R=1) OVER
                                                                  # (Number of of intervention =1 when R=0 / Number of of intervention =0 when R=0 ) (controlled for C)
                                                                  # randomisation is assumed to have more effect on compliance than confounder 
            
          b0<-log(odd0)
          b1<-log (odd1)
          b2<-log (odd2)
            
          logit.pi<- b0+b1[i]*simdata$confounder+b2*simdata$randomisation
          pi<-exp(logit.pi)/(1+exp(logit.pi))
            
          simdata$intervention <- rbinom(2*n,1,prob=pi)
            
          #OUTCOME dependent on confounders and intervention
          simdata$outcome <- rep(NA,2*n)
          
          odd0<- 1/5                                            # Number of outcome=1 when C=0 I=0/ Number of outcome=0 when C=0 I=0
          odd1<- confounder.eff.o                               # (Number of outcome=1 when C=1 / Number of outcome=0 when when C=1) OVER
                                                                  # (Number of of outcome =1 when C=0 / Number of of outcome =0 when C=0 ) (controlled for R)
                                                                  # Strength of confounder on outcome
          odd2<- (p1/(1-p1))/(p0/(1-p0))                        # (Number of outcome=1 when I=1 / Number of outcome=0 when when I=1) OVER
                                                                  # (Number of of outcome =1 when I=0 / Number of of outcome=0 when I=0 ) (controlled for C)
            
          b0<-log(odd0)
          b1<-log (odd1)
          b2<-log (odd2)
            
          logit.pi<- b0+b1*simdata$confounder+b2*simdata$intervention
          pi<-exp(logit.pi)/(1+exp(logit.pi))
            
          simdata$outcome <- rbinom(2*n,1,prob=pi)
            
          #generate proportions in simulated data 
          p00.vector<- (filter(simdata,intervention==0, simdata$randomisation==0))$outcome; p00.value<-mean(p00.vector)
          p11.vector<- (filter(simdata,intervention==1, simdata$randomisation==1))$outcome; p11.value<-mean(p11.vector)
            
          pz1.vector<- (filter(simdata, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
          pz0.vector<- (filter(simdata, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
            
          pd1.vector<- (filter(simdata, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
          pd0.vector<- (filter(simdata, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
          
          #estimate treatment effects 
          eff.itt[l]<- pz1.value-pz0.value            #intention to treat 
          eff.pp[l] <- p11.value-p00.value            #per protocol
          eff.at[l] <- pd1.value-pd0.value            #as treated 
            
          pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
          # Linear binomial generalized linear models with inverse probability weights
          E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=pp, na.action=na.exclude) #calculate denominators used in inverse probability weights
          ps<-predict(E.out, type="response")
          sptw<-pp$intervention*mean(pp$intervention)/ps+(1-pp$intervention)*(1-mean(pp$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
          eff.mpp[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=pp)))$coefficients['intervention']# fit linear binomial model to the weighted data
 
          m<- glm(intervention~confounder, family = binomial(), data=simdata) #propensity score matching 
          pscore<- m$fitted.values
          mout<-matchit(intervention~confounder, data=simdata, method='nearest')
          match<-Match(Tr=simdata$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
          matched<-simdata[unlist(match[c('index.treated','index.control')]),]
          y_trt<-matched$outcome[matched$intervention==1]
          y_con<-matched$outcome[matched$intervention==0]
          diffy<-y_trt-y_con
          eff.ps[l]<-(t.test(diffy))$estimate
            
          ivmodel=ivreg(simdata$outcome ~ simdata$intervention, ~ simdata$randomisation, x=TRUE, data=simdata)
          eff.iv[l]<-ivmodel$coef[2]                  #iv with 2 stage regression
            
          #proportion of patients who were compliant to protocol 
          .prop.c[l]<- nrow(pp)/(2*n)
          
        }
            
        #power calculations mean
        power.iv[i]<- mean(eff.iv<lower.value)
        power.itt[i]<- mean(eff.itt<lower.value)
        power.at[i]<- mean(eff.at<lower.value)
        power.pp[i]<- mean(eff.pp<lower.value)
        power.mpp[i]<- mean(eff.mpp<lower.value)
        power.ps[i]<- mean(eff.ps<lower.value)
        
        #proportion of patients who were compliant to protocol 
        prop.c [i]<- mean(.prop.c)
        }
      return(data.frame(prop.c, power.iv, power.itt, power.pp, power.at, power.mpp, power.ps))
      }
    sim.df<-sim()
      
    plot <- ggplot(sim.df, aes(sim.df[1]))+ 
        geom_point(aes(y=sim.df[2], colour="IV Power")) + 
        geom_point(aes(y=sim.df[3], colour="ITT Power")) + 
        geom_point(aes(y=sim.df[4], colour="PP Power")) + 
        geom_point(aes(y=sim.df[5], colour="AT Power")) + 
        geom_point(aes(y=sim.df[6], colour="Modified PP Power")) +
        geom_point(aes(y=sim.df[7], colour="Propensity score Power")) + 
        geom_smooth(aes(y=sim.df[2], colour="IV Power"),se=FALSE)+
        geom_smooth(aes(y=sim.df[3], colour="ITT Power"),se=FALSE)+
        geom_smooth(aes(y=sim.df[4], colour="PP Power"),se=FALSE)+
        geom_smooth(aes(y=sim.df[5], colour="AT Power"),se=FALSE)+
        geom_smooth(aes(y=sim.df[6], colour="Modified PP Power"),se=FALSE)+
        geom_smooth(aes(y=sim.df[7], colour="Propensity score Power"),se=FALSE)+
        geom_line(aes(y=.8, colour='Target power'),linetype="dotted") +
        xlab("Proportion of compliant participants")+
        ylab("Power")
      
      return(plot)
    }
    else (print ("NI margin must be positive, and p1-p0 (in terms of negative outcomes) must be less than NI margin in this simulation"))
}

p1<-simdata.power(n=230, p1=0.4, p0=0.4, NImargin=0.12, confounder.eff.o=(8/1)/(1/5) ,nIterations=1000) 
p2<-simdata.power(n=230, p1=0.4, p0=0.4, NImargin=0.12, confounder.eff.o=(4/1)/(1/5), nIterations=1000)
p3<-simdata.power(n=230, p1=0.4, p0=0.4, NImargin=0.12, confounder.eff.o=(2/1)/(1/5), nIterations=1000)

ggarrange(p1, p2, p3, ncol = 3, nrow = 1)

p4<-simdata.power(n=230,p1=0.2,p0=0.4,NImargin=0.12, confounder.eff.o=(4/1)/(1/5), nIterations=1000)                                             
p5<-simdata.power(n=230,p1=0.3,p0=0.4,NImargin=0.12, confounder.eff.o=(4/1)/(1/5), nIterations=1000) 
p6<-p2
p7<-simdata.power(n=230,p1=0.5,p0=0.4,NImargin=0.12, confounder.eff.o=(4/1)/(1/5), nIterations=1000) 

ggarrange(p4, p5, p6, p7,  
          ncol = 2, nrow = 2)

p9<-p2
p10<-simdata.power(n=400, p1=0.4, p0=0.4, NImargin=0.12, confounder.eff.o=(4/1)/(1/5), nIterations=1000)
p11<-simdata.power(n=600, p1=0.4, p0=0.4, NImargin=0.12, confounder.eff.o=(4/1)/(1/5), nIterations=1000)

ggarrange(p9, p10, p11, ncol = 3, nrow = 1)

##########GROUP SEQUENTIAL#########
###################################
simdata.gs<- function(n, p1, p0, NImargin, confounder.eff.o,nIterations){  
  # n: number of participants per group:
  # p0: survival rate in intervention 0; 
  # p1: survival rate in intervention 1;
  # confounder.eff.o: effect of counfounder on outcome 0:1
  # NImargin: non inferiority margin (negative);
  # nIterations: number of iterations 
  
  if ((NImargin <0) & (p1-p0 > NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
    if (nIterations>1) {
    alpha=0.025
    odd1<- seq((1/1)/(1/1),(10/1)/(1/10), by=9)
    
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
    
    #simulate and derive treatment effect 
    sim<- function() {
      for(i in 1:length(odd1)) {  print (i)
        for(l in 1:nIterations) { 
          #simulate data frame 
          
          #RANDOMISATION ratio 1:1
          simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))     #half randomised to 1, half randomised to 0 
          simdata$id<- seq(1,(2*n), by=1)
            
          #confounder normal distribution ranging 0-1 
          simdata$confounder<- rep(NA,2*n) 
          simdata$confounder<- rnorm((2*n),6, 2) # assume counfounder to be SOFA(severity), normally distributed 
          simdata$confounder [which(simdata$confounder<0)] <- -1*(simdata$confounder[which(simdata$confounder<0)]) #make sure all sofa scores are positive 
          simdata$confounder<- simdata$confounder/max(simdata$confounder) # scale from 0-1
          
          #INTERVENTION dependent on randomisation and confounders
          simdata$intervention <- rep(NA,2*n)
          
          odd0<- 1/5                                            #Number of intervention=1 when C=0 R=0/ Number of intervention=0 when C=0 R=0
          odd1<- seq((1/1)/(1/1),(10/1)/(1/10), by=9)           #(Number of intervention=1 when when C=1 / Number of intervention=0 when when C=1) OVER
                                                                # (Number of of intervention =1 when C=0 / Number of of intervention =0 when C=0 ) (controlled for R)
                                                                # if confounder has no negative effect on compliance - (1/1)/(1/1)
                                                                # if confounder has strong negative effect on compliance - (1/10)/(10/1)
          odd2<- (50/1)/(1/50)                                  #(Number of intervention=1 when when R=1 / Number of intervention=0 when when R=1) OVER
                                                                # (Number of of intervention =1 when R=0 / Number of of intervention =0 when R=0 ) (controlled for C)
                                                                # randomisation is assumed to have more effect on compliance than confounder 
          
          b0<-log(odd0)
          b1<-log (odd1)
          b2<-log (odd2)
          
          logit.pi<- b0+b1[i]*simdata$confounder+b2*simdata$randomisation
          pi<-exp(logit.pi)/(1+exp(logit.pi))
          
          simdata$intervention <- rbinom(2*n,1,prob=pi)
          
          #OUTCOME dependent on confounders and intervention
          simdata$outcome <- rep(NA,2*n)
          
          odd0<- 1/5                                            # Number of outcome=1 when C=0 I=0/ Number of outcome=0 when C=0 I=0
          odd1<- confounder.eff.o                               # (Number of outcome=1 when C=1 / Number of outcome=0 when when C=1) OVER
                                                                # (Number of of outcome =1 when C=0 / Number of of outcome =0 when C=0 ) (controlled for R)
                                                                # Strength of confounder on outcome
          odd2<- (p1/(1-p1))/(p0/(1-p0))                        # (Number of outcome=1 when I=1 / Number of outcome=0 when when I=1) OVER
                                                                # (Number of of outcome =1 when I=0 / Number of of outcome=0 when I=0 ) (controlled for C)
          
          b0<-log(odd0)
          b1<-log (odd1)
          b2<-log (odd2)
          
          logit.pi<- b0+b1*simdata$confounder+b2*simdata$intervention
          pi<-exp(logit.pi)/(1+exp(logit.pi))
          
          simdata$outcome <- rbinom(2*n,1,prob=pi)
          
          #Shuffle rows 
          simdata <- simdata[sample(nrow(simdata)),]
          
          #data for each interim analysis 
          s1<-simdata[1:(2*n/4),]
          s2<-simdata[1:(2*n/4*2),]
          s3<-simdata[1:(2*n/4*3),]
          s4<-simdata[1:(2*n),]
          
          #Compute z values at each interim analysis (3 interim and 1 final) 
          
          #######s1
          #generate proportions in simulated data s1
          p11.vector<- (filter(s1,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
          p00.vector<- (filter(s1,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
          pz1.vector<- (filter(s1,randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
          pz0.vector<- (filter(s1,randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
          pd1.vector<- (filter(s1,intervention==1))$outcome;pd1.value<- mean(pd1.vector)
          pd0.vector<- (filter(s1,intervention==0))$outcome;pd0.value<- mean(pd0.vector)
          
          ###estimate treatment effects s1
          eff.itt1[l]<- pz1.value-pz0.value            #intention to treat 
          eff.pp1[l] <- p11.value-p00.value            #per protocol
          eff.at1[l] <- pd1.value-pd0.value            #as treated 
          
          pp1<- rbind(filter(s1, randomisation==1 & intervention==1), (filter(s1, randomisation==0 & intervention==0))) # perprotocol population
          # Linear binomial generalized linear models with inverse probability weights
          E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=s1, na.action=na.exclude) #calculate denominators used in inverse probability weights
          ps<-predict(E.out, type="response")
          sptw<-s1$intervention*mean(s1$intervention)/ps+(1-s1$intervention)*(1-mean(s1$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
          eff.mpp1[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=s1)))$coefficients['intervention']# fit linear binomial model to the weighted data
          
          m<- glm(intervention~confounder, family = binomial(), data=s1) #propensity score matching 
          pscore<- m$fitted.values
          mout<-matchit(intervention~confounder, data=s1, method='nearest')
          match<-Match(Tr=s1$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
          matched<-simdata[unlist(match[c('index.treated','index.control')]),]
          y_trt<-matched$outcome[matched$intervention==1]
          y_con<-matched$outcome[matched$intervention==0]
          diffy<-y_trt-y_con
          eff.ps1[l]<-(t.test(diffy))$estimate
          
          ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s1)
          eff.iv1[l]<-ivmodel$coef[2]                 #iv with 2 stage regression
          
          #######s2
          #generate proportions in simulated data s2
          p11.vector<- (filter(s2,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
          p00.vector<- (filter(s2,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
          pz1.vector<- (filter(s2, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
          pz0.vector<- (filter(s2, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
          pd1.vector<- (filter(s2, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
          pd0.vector<- (filter(s2, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
          
          ###estimate treatment effects s2
          eff.itt2[l]<- pz1.value-pz0.value            #intention to treat 
          eff.pp2[l] <- p11.value-p00.value            #per protocol
          eff.at2[l] <- pd1.value-pd0.value            #as treated 
          
          pp2<- rbind(filter(s2, randomisation==1 & intervention==1), (filter(s2, randomisation==0 & intervention==0))) # perprotocol population
          # Linear binomial generalized linear models with inverse probability weights
          E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=s2, na.action=na.exclude) #calculate denominators used in inverse probability weights
          ps<-predict(E.out, type="response")
          sptw<-s2$intervention*mean(s2$intervention)/ps+(1-s2$intervention)*(1-mean(s2$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
          eff.mpp2[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=s2)))$coefficients['intervention']# fit linear binomial model to the weighted data
       
          m<- glm(intervention~confounder, family = binomial(), data=s2) #propensity score matching 
          pscore<- m$fitted.values
          mout<-matchit(intervention~confounder, data=s2, method='nearest')
          match<-Match(Tr=s2$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
          matched<-simdata[unlist(match[c('index.treated','index.control')]),]
          y_trt<-matched$outcome[matched$intervention==1]
          y_con<-matched$outcome[matched$intervention==0]
          diffy<-y_trt-y_con
          eff.ps2[l]<-(t.test(diffy))$estimate
          
          ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s2)
          eff.iv2[l]<-ivmodel$coef[2]                 #iv with 2 stage regression
          
          #######s3
          #generate proportions in simulated data s3
          p11.vector<- (filter(s3,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
          p00.vector<- (filter(s3,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
          pz1.vector<- (filter(s3, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
          pz0.vector<- (filter(s3, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
          pd1.vector<- (filter(s3, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
          pd0.vector<- (filter(s3, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
          
          ###estimate treatment effects s3
          eff.itt3[l]<- pz1.value-pz0.value            #intention to treat 
          eff.pp3[l] <- p11.value-p00.value            #per protocol
          eff.at3[l] <- pd1.value-pd0.value            #as treated 
          
          pp3<- rbind(filter(s3, randomisation==1 & intervention==1), (filter(s3, randomisation==0 & intervention==0))) # perprotocol population
          # Linear binomial generalized linear models with inverse probability weights
          E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=s3, na.action=na.exclude) #calculate denominators used in inverse probability weights
          ps<-predict(E.out, type="response")
          sptw<-s3$intervention*mean(s3$intervention)/ps+(1-s3$intervention)*(1-mean(s3$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
          eff.mpp3[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=s3)))$coefficients['intervention']# fit linear binomial model to the weighted data
          
          m<- glm(intervention~confounder, family = binomial(), data=s3) #propensity score matching 
          pscore<- m$fitted.values
          mout<-matchit(intervention~confounder, data=s3, method='nearest')
          match<-Match(Tr=s3$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
          matched<-simdata[unlist(match[c('index.treated','index.control')]),]
          y_trt<-matched$outcome[matched$intervention==1]
          y_con<-matched$outcome[matched$intervention==0]
          diffy<-y_trt-y_con
          eff.ps3[l]<-(t.test(diffy))$estimate
          
          ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s3)
          eff.iv3[l]<-ivmodel$coef[2]                 #iv with 2 stage regression
          
          #######s4
          #generate proportions in simulated data s4
          p11.vector<- (filter(s4,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
          p00.vector<- (filter(s4,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
          pz1.vector<- (filter(s4, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
          pz0.vector<- (filter(s4, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
          pd1.vector<- (filter(s4, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
          pd0.vector<- (filter(s4, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
          
          ###estimate treatment effects s4
          eff.itt4[l]<- pz1.value-pz0.value            #intention to treat 
          eff.pp4[l] <- p11.value-p00.value            #per protocol
          eff.at4[l] <- pd1.value-pd0.value            #as treated 
          
          pp4<- rbind(filter(s4, randomisation==1 & intervention==1), (filter(s4, randomisation==0 & intervention==0))) # perprotocol population
          # Linear binomial generalized linear models with inverse probability weights
          E.out<-glm(intervention~confounder,family=binomial(link="logit"), data=s4, na.action=na.exclude) #calculate denominators used in inverse probability weights
          ps<-predict(E.out, type="response")
          sptw<-s4$intervention*mean(s4$intervention)/ps+(1-s4$intervention)*(1-mean(s4$intervention))/(1-ps)#create stabilized weights, using a null model with intervention as the dependent variable
          eff.mpp4[l]<-((geeglm(outcome~intervention, family=binomial(link="identity"), weight=sptw, id=id, data=s4)))$coefficients['intervention']# fit linear binomial model to the weighted data
          
          m<- glm(intervention~confounder, family = binomial(), data=s4) #propensity score matching 
          pscore<- m$fitted.values
          mout<-matchit(intervention~confounder, data=s4, method='nearest')
          match<-Match(Tr=s4$intervention, M=1, X= logit(pscore), replace=FALSE, caliper=0.2)
          matched<-simdata[unlist(match[c('index.treated','index.control')]),]
          y_trt<-matched$outcome[matched$intervention==1]
          y_con<-matched$outcome[matched$intervention==0]
          diffy<-y_trt-y_con
          eff.ps4[l]<-(t.test(diffy))$estimate
          
          ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s4)
          eff.iv4[l]<-ivmodel$coef[2]                 #iv with 2 stage regression
          
          #proportion of patients who were compliant to protocol 
          .prop.c[l]<- nrow(pp4)/(2*n)
        }
        
        ###sds of effects s1
        sd.eff.itt1<-sd(eff.itt1)
        sd.eff.pp1<- sd(eff.pp1)
        sd.eff.at1<- sd(eff.at1)
        sd.eff.ps1<- sd(eff.ps1)
        sd.eff.mpp1<-sd(eff.mpp1)
        sd.eff.iv1<- sd(eff.iv1)
        
        ###z values s1
        z.itt1[i]<-(mean(eff.itt1)-NImargin)/sd.eff.itt1
        z.iv1[i]<-(mean(eff.iv1)-NImargin)/sd.eff.iv1
        z.pp1[i]<-(mean(eff.pp1)-NImargin)/sd.eff.pp1
        z.mpp1[i]<-(mean(eff.mpp1)-NImargin)/sd.eff.mpp1
        z.at1[i]<-(mean(eff.at1)-NImargin)/sd.eff.at1
        z.ps1[i]<-(mean(eff.ps1)-NImargin)/sd.eff.ps1
        
        ###sds of effects s2
        sd.eff.itt2<-sd(eff.itt2)
        sd.eff.pp2<- sd(eff.pp2)
        sd.eff.at2<- sd(eff.at2)
        sd.eff.ps2<- sd(eff.ps2)
        sd.eff.mpp2<-sd(eff.mpp2)
        sd.eff.iv2<- sd(eff.iv2)
        
        ###z values s2
        z.itt2[i]<-(mean(eff.itt2)-NImargin)/sd.eff.itt2
        z.iv2[i]<-(mean(eff.iv2)-NImargin)/sd.eff.iv2
        z.pp2[i]<-(mean(eff.pp2)-NImargin)/sd.eff.pp2
        z.mpp2[i]<-(mean(eff.mpp2)-NImargin)/sd.eff.mpp2
        z.at2[i]<-(mean(eff.at2)-NImargin)/sd.eff.at2
        z.ps2[i]<-(mean(eff.ps2)-NImargin)/sd.eff.ps2
        
        ###sds of effects s3
        sd.eff.itt3<-sd(eff.itt3)
        sd.eff.pp3<- sd(eff.pp3)
        sd.eff.at3<- sd(eff.at3)
        sd.eff.ps3<- sd(eff.ps3)
        sd.eff.mpp3<-sd(eff.mpp3)
        sd.eff.iv3<- sd(eff.iv3)
        
        ###z values s3
        z.itt3[i]<-(mean(eff.itt3)-NImargin)/sd.eff.itt3
        z.iv3[i]<-(mean(eff.iv3)-NImargin)/sd.eff.iv3
        z.pp3[i]<-(mean(eff.pp3)-NImargin)/sd.eff.pp3
        z.mpp3[i]<-(mean(eff.mpp3)-NImargin)/sd.eff.mpp3
        z.at3[i]<-(mean(eff.at3)-NImargin)/sd.eff.at3
        z.ps3[i]<-(mean(eff.ps3)-NImargin)/sd.eff.ps3
        
        ###sds of effects s4
        sd.eff.itt4<-sd(eff.itt4)
        sd.eff.pp4<- sd(eff.pp4)
        sd.eff.at4<- sd(eff.at4)
        sd.eff.ps4<- sd(eff.ps4)
        sd.eff.mpp4<-sd(eff.mpp4)
        sd.eff.iv4<- sd(eff.iv4)
        
        ###z values s4
        z.itt4[i]<-(mean(eff.itt4)-NImargin)/sd.eff.itt4
        z.iv4[i]<-(mean(eff.iv4)-NImargin)/sd.eff.iv4
        z.pp4[i]<-(mean(eff.pp4)-NImargin)/sd.eff.pp4
        z.mpp4[i]<-(mean(eff.mpp4)-NImargin)/sd.eff.mpp4
        z.at4[i]<-(mean(eff.at4)-NImargin)/sd.eff.at4
        z.ps4[i]<-(mean(eff.ps4)-NImargin)/sd.eff.ps4
        
        .z.itt<-list(z.itt1, z.itt2,z.itt3,z.itt4)
        .z.at<-list(z.at1, z.at2,z.at3,z.at4)
        .z.pp<-list(z.pp1,z.pp2,z.pp3,z.pp4)
        .z.mpp<-list(z.mpp1,z.mpp2,z.mpp3,z.mpp4)
        .z.ps<-list(z.ps1,z.ps2,z.ps3,z.ps4)
        .z.iv<-list(z.iv1,z.iv2,z.iv3,z.iv4)
        
        #proportion of patients who were compliant to protocol 
        prop.c [i]<- mean(.prop.c)
      } 
      
      z<-gsDesign(k=4, n.fix=2*n, delta0 = NImargin, n.I=c(nrow(s1),nrow(s2) ,nrow(s3),nrow(s4)), maxn.IPlan = 2*n, beta=0.2) #boundary 
      z<-z$upper$bound #4 boundaries 
      z.itt<-data.frame(.z.itt)
      z.at<-data.frame(.z.at)
      z.pp<-data.frame(.z.pp)
      z.mpp<-data.frame(.z.mpp)
      z.ps<-data.frame(.z.ps)
      z.iv<-data.frame(.z.iv)
      z.c<-data.frame(cbind(prop.c, z.itt, z.at, z.pp, z.mpp,z.ps,z.iv, 
                            rep(z[1],length(odd1)),
                            rep(z[2],length(odd1)),
                            rep(z[3],length(odd1)),
                            rep(z[4],length(odd1))))
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
    print(z.c)
    
    plot.z.itt<- ggplot(z.c, aes(z.c[1]))+ 
      geom_point(aes(y = z.c[2], colour="Z value ITT 1" ))+
      geom_point(aes(y = z.c[3], colour="Z value ITT 2"))+
      geom_point(aes(y = z.c[4], colour="Z value ITT 3"))+
      geom_point(aes(y = z.c[5], colour="Z value ITT 4"))+
      geom_smooth(aes(y=z.c[2], colour="Z value ITT 1"),se=FALSE)+
      geom_smooth(aes(y=z.c[3], colour="Z value ITT 2"),se=FALSE)+
      geom_smooth(aes(y=z.c[4], colour="Z value ITT 3"),se=FALSE)+
      geom_smooth(aes(y=z.c[5], colour="Z value ITT 4"),se=FALSE)+
      geom_line(aes(y = z1, colour="Z value ITT 1"),linetype="dotted")+
      geom_line(aes(y = z2, colour="Z value ITT 2"),linetype="dotted")+
      geom_line(aes(y = z3, colour="Z value ITT 3"),linetype="dotted")+
      geom_line(aes(y = z4, colour="Z value ITT 4"),linetype="dotted")+
      xlab("Proportion of compliant participants")+
      ylab("Z value")
    
    plot.z.at<- ggplot(z.c, aes(z.c[1]))+ 
      geom_point(aes(y = z.c[6], colour="Z value AT 1" ))+
      geom_point(aes(y = z.c[7], colour="Z value AT 2"))+
      geom_point(aes(y = z.c[8], colour="Z value AT 3"))+
      geom_point(aes(y = z.c[9], colour="Z value AT 4"))+
      geom_smooth(aes(y=z.c[6], colour="Z value AT 1"),se=FALSE)+
      geom_smooth(aes(y=z.c[7], colour="Z value AT 2"),se=FALSE)+
      geom_smooth(aes(y=z.c[8], colour="Z value AT 3"),se=FALSE)+
      geom_smooth(aes(y=z.c[9], colour="Z value AT 4"),se=FALSE)+
      geom_line(aes(y = z1, colour="Z value AT 1"),linetype="dotted")+
      geom_line(aes(y = z2, colour="Z value AT 2"),linetype="dotted")+
      geom_line(aes(y = z3, colour="Z value AT 3"),linetype="dotted")+
      geom_line(aes(y = z4, colour="Z value AT 4"),linetype="dotted")+
      xlab("Proportion of compliant participants")+
      ylab("Z value")
    
    plot.z.pp<- ggplot(z.c, aes(z.c[1]))+ 
      geom_point(aes(y = z.c[10], colour="Z value PP 1" ))+
      geom_point(aes(y = z.c[11], colour="Z value PP 2"))+
      geom_point(aes(y = z.c[12], colour="Z value PP 3"))+
      geom_point(aes(y = z.c[13], colour="Z value PP 4"))+
      geom_smooth(aes(y=z.c[10], colour="Z value PP 1"),se=FALSE)+
      geom_smooth(aes(y=z.c[11], colour="Z value PP 2"),se=FALSE)+
      geom_smooth(aes(y=z.c[12], colour="Z value PP 3"),se=FALSE)+
      geom_smooth(aes(y=z.c[13], colour="Z value PP 4"),se=FALSE)+
      geom_line(aes(y = z1, colour="Z value PP 1"),linetype="dotted")+
      geom_line(aes(y = z2, colour="Z value PP 2"),linetype="dotted")+
      geom_line(aes(y = z3, colour="Z value PP 3"),linetype="dotted")+
      geom_line(aes(y = z4, colour="Z value PP 4"),linetype="dotted")+
      xlab("Proportion of compliant participants")+
      ylab("Z value")
    
    plot.z.mpp<- ggplot(z.c, aes(z.c[1]))+ 
      geom_point(aes(y = z.c[14], colour="Z value MPP 1" ))+
      geom_point(aes(y = z.c[15], colour="Z value MPP 2"))+
      geom_point(aes(y = z.c[16], colour="Z value MPP 3"))+
      geom_point(aes(y = z.c[17], colour="Z value MPP 4"))+
      geom_smooth(aes(y=z.c[14], colour="Z value MPP 1"),se=FALSE)+
      geom_smooth(aes(y=z.c[15], colour="Z value MPP 2"),se=FALSE)+
      geom_smooth(aes(y=z.c[16], colour="Z value MPP 3"),se=FALSE)+
      geom_smooth(aes(y=z.c[17], colour="Z value MPP 4"),se=FALSE)+
      geom_line(aes(y = z1, colour="Z value MPP 1"),linetype="dotted")+
      geom_line(aes(y = z2, colour="Z value MPP 2"),linetype="dotted")+
      geom_line(aes(y = z3, colour="Z value MPP 3"),linetype="dotted")+
      geom_line(aes(y = z4, colour="Z value MPP 4"),linetype="dotted")+
      xlab("Proportion of compliant participants")+
      ylab("Z value")
    
    plot.z.ps<- ggplot(z.c, aes(z.c[1]))+ 
      geom_point(aes(y = z.c[18], colour="Z value PS 1" ))+
      geom_point(aes(y = z.c[19], colour="Z value PS 2"))+
      geom_point(aes(y = z.c[20], colour="Z value PS 3"))+
      geom_point(aes(y = z.c[21], colour="Z value PS 4"))+
      geom_smooth(aes(y=z.c[18], colour="Z value PS 1"),se=FALSE)+
      geom_smooth(aes(y=z.c[19], colour="Z value PS 2"),se=FALSE)+
      geom_smooth(aes(y=z.c[20], colour="Z value PS 3"),se=FALSE)+
      geom_smooth(aes(y=z.c[21], colour="Z value PS 4"),se=FALSE)+
      geom_line(aes(y = z1, colour="Z value PS 1"),linetype="dotted")+
      geom_line(aes(y = z2, colour="Z value PS 2"),linetype="dotted")+
      geom_line(aes(y = z3, colour="Z value PS 3"),linetype="dotted")+
      geom_line(aes(y = z4, colour="Z value PS 4"),linetype="dotted")+
      xlab("Proportion of compliant participants")+
      ylab("Z value")
    
    plot.z.iv<- ggplot(z.c, aes(z.c[1]))+ 
      geom_point(aes(y = iv1, colour="Z value IV 1" ))+
      geom_point(aes(y = iv2, colour="Z value IV 2"))+
      geom_point(aes(y = iv3, colour="Z value IV 3"))+
      geom_point(aes(y = iv4, colour="Z value IV 4"))+
      geom_smooth(aes(y=iv1, colour="Z value IV 1"),se=FALSE)+
      geom_smooth(aes(y=iv2, colour="Z value IV 2"),se=FALSE)+
      geom_smooth(aes(y=iv3, colour="Z value IV 3"),se=FALSE)+
      geom_smooth(aes(y=iv4, colour="Z value IV 4"),se=FALSE)+
      geom_line(aes(y = z1, colour="Z value IV 1"),linetype="dotted")+
      geom_line(aes(y = z2, colour="Z value IV 2"),linetype="dotted")+
      geom_line(aes(y = z3, colour="Z value IV 3"),linetype="dotted")+
      geom_line(aes(y = z4, colour="Z value IV 4"),linetype="dotted")+
      xlab("Proportion of compliant participants")+
      ylab("Z value")
    
    plot<- ggarrange(plot.z.itt, plot.z.at, plot.z.pp, plot.z.mpp,plot.z.ps,plot.z.iv, ncol = 3, nrow = 2)
    
    return(plot)
    }else (print ("Number of iterations must be more than 1."))
  } else (print ("NI margin must be negative, and p1-p0 (in terms of favorable outcomes) must be more than NI margin in this simulation"))
}

gs1<-simdata.gs(n=230,p0=0.6, p1=0.6,confounder.eff.o=40,nIterations=250,NImargin=-.12)  
gs2<-simdata.gs(n=300,p0=0.6, p1=0.6,confounder.eff.o=20,nIterations=250,NImargin=-.12) 
gs3<-simdata.gs(n=230,p0=0.6, p1=0.6,confounder.eff.o=10,nIterations=250,NImargin=-.12) 

 
gs5<-simdata.gs(n=230,p0=0.6, p1=0.6,confounder.eff.o=20,nIterations=250,NImargin=-.12) 
gs6<-simdata.gs(n=230,p0=0.4, p1=0.6,confounder.eff.o=20,nIterations=250,NImargin=-.12) 


gs8<-simdata.gs(n=300,p0=0.6, p1=0.6,confounder.eff.o=20,nIterations=250,NImargin=-.12) 
gs9<-simdata.gs(n=400,p0=0.4, p1=0.6,confounder.eff.o=20,nIterations=250,NImargin=-.12) 


#using gsDesign package to determine sample size for REGARD-VAP 
n.fix<- nBinomial (p1=0.6, p2=0.6, delta0=-0.12)
x<-gsDesign(n.fix=n.fix, k=4, beta=0.2, delta0=-0.12)
ceiling(x$n.I)
y<-gsProbability(theta=x$delta*seq(0,2,0.25),d=x)
plot(y, plottype=6, lty=2, lwd=3)


