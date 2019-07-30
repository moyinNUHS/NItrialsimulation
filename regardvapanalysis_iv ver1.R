########################################################################
###################Simulation for REGARD-VAP analysis###################
########################################################################

setwd("/Users/moyin/Desktop/VAP studd/Causal inference simulation") #set working directory 
rm(list=ls()) # Clean working environment 
library(Hmisc); library(rms); library(gsDesign);library(ivpack); library(ggplot2); 
library(data.table); library(dplyr); library(plotly); library(ggpubr) # Required libraries 

set.seed(1234)

##########BIAS##########
########################
simdata.bias<- function(n, wna, pc0, pc1, pa, pn, nIterations){  
  # n: number of participants per group; 
  # wna: proportion of never takers over proportion of always takers; 
  # pc0: mortality and recurrence rate in compliers in intervention 0; 
  # pc1: mortality and recurrence rate in compliers in intervention 1;
  # pa: mortality and recurrence rate in always takers; 
  # pn: mortality and recurrence rate in never takers;
  # nIterations: number of iterations 
  
  alpha=0.025
  
  #define proportions of patients
  wc=seq(0.1,1, by=0.05)    # proportion of compliers
  wa= (1-wc)/(wna+1)             # proportion of always takers 
  wn= wa*wna                     # proportion of never takers
  
  w00=wc+wn                      # Randomised to 0, intervention 0 (C+NT), proportion
  w01=rep(wn, length (wc))       # Randomised to 0, intervention 1 (NT), proportion
  w10= wa                        # Randomised to 1, intervention 0 (AT), proportion
  w11=wc+wa                      # Randomised to 1, intervention 1 (C+AT), proportion
  
  #define proportions of events 
  p00=(wc*pc0+wn*pn)/(wc+wn)     # Randomised to 0, intervention 0 (C+NT), outcome
  p01=rep(pn, length(wc))        # Randomised to 0, intervention 1 (NT), outcome
  p10=rep(pa, length(wc))        # Randomised to 1, intervention 0 (AT), outcome
  p11=(wc*pc1+wa*pa)/(wc+wa)     # Randomised to 1, intervention 1 (C+AT), outcome
  
  #make up vectors for simulations 
  eff.iv<-eff.itt<-eff.pp<-eff.at<-eff.cace<-c()                   
  .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-.eff.cace<-c()                    
  
  #simulate and derive treatment effect 
  sim<- function() {
    for(i in 1:length(wc)) { print(i)
      for(l in 1:nIterations) { 
        #simulate data frame 
        simdata<- data.frame (                           
          "randomisation" = c(rep(1,n), rep(0,n)))     #half randomised to 1, half randomised to 0 
        
        simdata$intervention <- c(rbinom(n,1,w11[i]),   # for those randomised to 1, intervention proportion w11 
                                  rbinom(n,1,w10[i]))   # for those randomised to 0, intervention proportion w10
        
        simdata$outcome <- rep(NA,2*n)
        for(j in 1:(2*n)){
          if(simdata$randomisation[j]==1 & simdata$intervention[j]==1){
            simdata$outcome[j] <- rbinom(1,1,prob=p11[i])# for those randomised to 1, intervention 1, outcome proportion p11
          }  
          else if(simdata$randomisation[j]==1 & simdata$intervention[j]==0){
            simdata$outcome[j] <-rbinom(1,1,prob=p01[i]) # for those randomised to 1, intervention 0, outcome proportion p01
          }
          else if(simdata$randomisation[j]==0 & simdata$intervention[j]==1){
            simdata$outcome[j] <-rbinom(1,1,prob=p10[i]) # for those randomised to 0, intervention 1, outcome proportion p10
          }
          else if(simdata$randomisation[j]==0 & simdata$intervention[j]==0){
            simdata$outcome[j] <-rbinom(1,1,prob=p00[i]) # for those randomised to 0, intervention 0, outcome proportion p00
          }
        }
        
        #generate proportions from simulated data  
        p00.value<- mean((filter(simdata, intervention==0, randomisation==0))$outcome)
        p11.value<- mean((filter(simdata, intervention==1, randomisation==1))$outcome)
        
        w1.value<- w11.value<-  mean((filter(simdata, randomisation==1))$intervention)
        w0.value<- w10.value<-  mean((filter(simdata, randomisation==0))$intervention)
        wc.value<- w1.value-w0.value
        if (wc.value==0) {wc.value=wc[i]}
        
        pz1.value<- mean((filter(simdata, randomisation==1))$outcome) 
        pz0.value<- mean((filter(simdata, randomisation==0))$outcome)
        
        pd1.value<- mean((filter(simdata, intervention==1))$outcome)
        pd0.value<- mean((filter(simdata, intervention==0))$outcome)
        
        #treatment effect 
        .eff.itt[l]<- pz1.value-pz0.value              #intention to treat 
        .eff.cace[l] <- .eff.itt[l]/wc.value           #complier average causal effect 
        .eff.pp[l] <- p11.value-p00.value              #per protocol 
        .eff.at[l] <- pd1.value-pd0.value              # as treated 
        ivmodel=ivreg(simdata$outcome ~ simdata$intervention, ~ simdata$randomisation, x=TRUE, data=simdata)
        .eff.iv[l]<-ivmodel$coef[2]                    #iv with 2 stage regression
        
        #mean of iterated data 
        eff.itt[i]<- mean(.eff.itt)
        eff.cace[i] <- mean(.eff.cace)
        eff.pp[i] <- mean(.eff.pp) 
        eff.at[i] <- mean(.eff.at)
        eff.iv[i]<-mean(.eff.iv, na.rm=TRUE)
       
      }
    }
    return(data.frame(wc, eff.cace, eff.iv, eff.itt, eff.pp, eff.at))
  }
  
  sim.df<-sim()
  
  sim.df$pred.eff.iv<- predict(lm(eff.iv~log(wc), data=sim.df))
  sim.df$pred.eff.itt<- predict(lm(eff.itt~log(wc), data=sim.df))
  sim.df$pred.eff.at<- predict(lm(eff.at~log(wc), data=sim.df))
  sim.df$pred.eff.pp<- predict(lm(eff.pp~log(wc), data=sim.df))
  sim.df$pred.eff.cace1<- predict(lm(eff.cace~log(wc), data=sim.df))
  sim.df$pred.eff.cace<- sim.df$pred.eff.cace1-sim.df$pred.eff.cace1/50
  
  bias.plot <- ggplot(sim.df, aes(sim.df[1]))+ 
    geom_point(aes(y=sim.df[2], colour="CACE effect")) + 
    geom_point(aes(y=sim.df[3], colour="IV effect")) + 
    geom_point(aes(y=sim.df[4], colour="ITT effect")) + 
    geom_point(aes(y=sim.df[5], colour="PP effect")) + 
    geom_point(aes(y=sim.df[6], colour="AT effect")) + 
    geom_line(aes(y = pred.eff.cace, colour="CACE effect"))+
    geom_line(aes(y = pred.eff.iv, colour="IV effect" ))+
    geom_line(aes(y = pred.eff.itt, colour="ITT effect"))+
    geom_line(aes(y = pred.eff.pp, colour="PP effect"))+
    geom_line(aes(y = pred.eff.at, colour="AT effect"))+
    geom_line(aes(y=pc1-pc0, colour='True effect'),linetype="dotted") +
    xlab("Proportion of compliers")+
    ylab("Effect")
  
  return(bias.plot)
} 

bias1<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.4, pn=0.4, nIterations=1000)
bias2<-simdata.bias(n=230,wna=2,pc1=0.4, pc0=0.4, pa=0.4, pn=0.4, nIterations=1000)
bias3<-simdata.bias(n=230,wna=3,pc1=0.4, pc0=0.4, pa=0.4, pn=0.4, nIterations=1000)

bias11<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.2, pa=0.4, pn=0.4, nIterations=1000)
bias21<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.4, pn=0.4, nIterations=1000)
bias31<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.6, pa=0.4, pn=0.4, nIterations=1000)
bias41<-simdata.bias(n=230,wna=1,pc1=0.5, pc0=0.3, pa=0.4, pn=0.4, nIterations=1000)
bias51<-simdata.bias(n=230,wna=1,pc1=0.5, pc0=0.5, pa=0.4, pn=0.4, nIterations=1000)
bias61<-simdata.bias(n=230,wna=1,pc1=0.5, pc0=0.7, pa=0.4, pn=0.4, nIterations=1000)

bias12<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.6, pn=0.4, nIterations=1000)
bias22<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.4, pn=0.4, nIterations=1000)
bias32<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.2, pn=0.4, nIterations=1000)
bias42<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.7, pn=0.5, nIterations=1000)
bias52<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.5, pn=0.5, nIterations=1000)
bias62<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.3, pn=0.5, nIterations=1000)


ggarrange(bias12, bias22, bias32, 
          bias42, bias52, bias62,
          ncol = 3, nrow = 2)

######TYPE 1 ERROR######
########################
simdata.type1<- function(n, wna, pc1, pa, pn, NImargin, nIterations){  
  # n: number of participants per group; 
  # wna: proportion of never takers over proportion of always takers; 
  # pc0: mortality and recurrence rate in compliers in intervention 0; 
  # pc1: mortality and recurrence rate in compliers in intervention 1;
  # pa: mortality and recurrence rate in always takers; 
  # pn: mortality and recurrence rate in never takers;
  # NImargin: non inferiority margin 
  # nIterations: number of iterations 
  
  alpha=0.025
  
  #define proportions of patients
  wc=seq(0.1,0.99, by=0.05)    # proportion of compliers
  wa= (1-wc)/(wna+1)             # proportion of always takers 
  wn= wa*wna                     # proportion of never takers
  
  w00=wc+wn                      # Randomised to 0, intervention 0 (C+NT), proportion
  w01=rep(wn, length (wc))       # Randomised to 0, intervention 1 (NT), proportion
  w10= wa                        # Randomised to 1, intervention 0 (AT), proportion
  w11=wc+wa                      # Randomised to 1, intervention 1 (C+AT), proportion
  
  #define proportions of events 
  pc0=pc1-NImargin               # data simulated with assumption of null hypothesis pc1-pc0=NImargin 
  p00=(wc*pc0+wn*pn)/(wc+wn)     # Randomised to 0, intervention 0 (C+NT), outcome
  p01=rep(pn, length(wc))        # Randomised to 0, intervention 1 (NT), outcome
  p10=rep(pa, length(wc))        # Randomised to 1, intervention 0 (AT), outcome
  p11=(wc*pc1+wa*pa)/(wc+wa)     # Randomised to 1, intervention 1 (C+AT), outcome
  
  #make up vectors for simulations 
  type1.error.iv<-type1.error.itt<-type1.error.pp<-type1.error.at<-type1.error.cace<- c()       
  .type1.error.iv<-.type1.error.itt<- .type1.error.pp<-.type1.error.at<-.type1.error.cace<-c()
  
  #simulate and derive type 1 error
  sim<- function() {
    for(i in 1:length(wc)) { print(i)
      for(l in 1:nIterations) { 
        #simulate data frame 
        simdata<- data.frame (                           
          "randomisation" = c(rep(1,n), rep(0,n)))     #half randomised to 1, half randomised to 0 
        
        simdata$intervention <- c(rbinom(n,1,w11[i]),   # for those randomised to 1, intervention proportion w11 
                                  rbinom(n,1,w10[i]))   # for those randomised to 0, intervention proportion w10
        
        simdata$outcome <- rep(NA,2*n)
        for(j in 1:(2*n)){
          if(simdata$randomisation[j]==1 & simdata$intervention[j]==1){
            simdata$outcome[j] <- rbinom(1,1,prob=p11[i])# for those randomised to 1, intervention 1, outcome proportion p11
          }  
          else if(simdata$randomisation[j]==1 & simdata$intervention[j]==0){
            simdata$outcome[j] <-rbinom(1,1,prob=p01[i]) # for those randomised to 1, intervention 0, outcome proportion p01
          }
          else if(simdata$randomisation[j]==0 & simdata$intervention[j]==1){
            simdata$outcome[j] <-rbinom(1,1,prob=p10[i]) # for those randomised to 0, intervention 1, outcome proportion p10
          }
          else if(simdata$randomisation[j]==0 & simdata$intervention[j]==0){
            simdata$outcome[j] <-rbinom(1,1,prob=p00[i]) # for those randomised to 0, intervention 0, outcome proportion p00
          }
        }
        
        #generate proportions in simulated data 
        w1.vector<- w11.vector<-(filter(simdata,simdata$randomisation==1))$intervention; w1.value<-mean(w1.vector)
        w0.vector<- w10.vector<-(filter(simdata,simdata$randomisation==0))$intervention; w0.value<-mean(w0.vector)
        wc.value<- w1.value-w0.value
        if (wc.value==0) {wc.value=wc[i]}
        
        p11.vector<- (filter(simdata,intervention==1, simdata$randomisation==1))$outcome; p11.value<-mean(p11.vector)
        p01.vector<- (filter(simdata,intervention==0, simdata$randomisation==1))$outcome; p01.value<-mean(p01.vector)
        p10.vector<- (filter(simdata,intervention==1, simdata$randomisation==0))$outcome; p10.value<-mean(p10.vector)
        p00.vector<- (filter(simdata,intervention==0, simdata$randomisation==0))$outcome; p00.value<-mean(p00.vector)
        
        pz1.vector<- (filter(simdata, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
        pz0.vector<- (filter(simdata, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
        pd1.vector<- (filter(simdata, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
        pd0.vector<- (filter(simdata, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
        
        #estimate treatment effects 
        eff.itt<- pz1.value-pz0.value            #intention to treat 
        eff.cace<-  eff.itt/wc.value             #complier average causal effect 
        eff.pp <- p11.value-p00.value            #per protocol
        eff.at <- pd1.value-pd0.value            #as treated 
        ivmodel=ivreg(simdata$outcome ~ simdata$intervention, ~ simdata$randomisation, x=TRUE, data=simdata)
        eff.iv<-ivmodel$coef[2]                  #iv with 2 stage regression
        
        #variances of proportions and effects  
        vw1<-var(w1.vector)
        vw0<-var(w0.vector)
        
        v11<-var(p11.vector)
        v01<-var(p01.vector)
        v10<-var(p10.vector)
        v00<-var(p00.vector)
        
        v.1<-vw1*(v11+v01+(p11.value-p01.value)^2)+w1.value^2*v11+(1-w1.value)^2*v01
        v.0<-vw0*(v10+v00+(p10.value-p00.value)^2)+w0.value^2*v10+(1-w0.value)^2*v00
        
        var.eff.cace<- (v.1+v.0)/wc[i]^2 + eff.iv^2*(vw1+vw0)/wc[i]^2 - (2*eff.iv/wc[i]^2)*((p11.value-p01.value)*vw1+(p10.value-p00.value)*vw0)
        var.eff.itt<- pz1.value*(1-pz1.value)/length(pz1.vector) + pz0.value*(1-pz0.value)/length(pz0.vector)
        var.eff.pp<-  p11.value*(1-p11.value)/length(p11.vector) + p00.value*(1-p00.value)/length(p00.vector)
        var.eff.at<-  pd1.value*(1-pd1.value)/length(pd1.vector) + pd0.value*(1-pd0.value)/length(pd0.vector)
        if (is.na(eff.iv)) {var.eff.iv<-var.eff.cace} else {var.eff.iv<-robust.se(ivmodel)[2,2]^2 }
        
        #type 1 error
        .type1.error.cace[l]<- pnorm(qnorm(alpha)+(NImargin-eff.cace)/(var.eff.cace^0.5))
        .type1.error.itt[l]<-pnorm(qnorm(alpha)+(NImargin-eff.itt)/(var.eff.itt^0.5))
        .type1.error.pp[l]<- pnorm(qnorm(alpha)+(NImargin-eff.pp)/(var.eff.pp^0.5))
        .type1.error.at[l]<- pnorm(qnorm(alpha)+(NImargin-eff.at)/(var.eff.at^0.5))
        .type1.error.iv[l]<- pnorm(qnorm(alpha)+(NImargin-eff.iv)/var.eff.iv^0.5)
        
        #mean of Type 1 error
        type1.error.iv[i]<-mean(.type1.error.iv, na.rm=TRUE)
        type1.error.itt[i]<-mean(.type1.error.itt, na.rm=TRUE)
        type1.error.pp[i]<-mean(.type1.error.pp, na.rm=TRUE)
        type1.error.at[i]<-mean(.type1.error.at, na.rm=TRUE)
        type1.error.cace[i]<-mean(.type1.error.cace, na.rm=TRUE)
      }
    }
    
    return(data.frame(wc, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at, type1.error.cace))
  }
  
  sim.df<-sim()
  
  sim.df$pred.t1.iv<- predict(lm(type1.error.iv~log(wc), data=sim.df))
  sim.df$pred.t1.itt<- predict(lm(type1.error.itt~log(wc), data=sim.df))
  sim.df$pred.t1.at<- predict(lm(type1.error.at~log(wc), data=sim.df))
  sim.df$pred.t1.pp<- predict(lm(type1.error.pp~log(wc), data=sim.df))
  sim.df$pred.t1.cace1<- predict(lm(type1.error.cace~log(wc), data=sim.df))
  sim.df$pred.t1.cace<- sim.df$pred.t1.cace1-sim.df$pred.t1.cace1/50
  
  bias.plot <- ggplot(sim.df, aes(sim.df[1]))+ 
    geom_point(aes(y=sim.df[2], colour="IV Type 1 error")) + 
    geom_point(aes(y=sim.df[3], colour="ITT Type 1 error")) + 
    geom_point(aes(y=sim.df[4], colour="PP Type 1 error")) + 
    geom_point(aes(y=sim.df[5], colour="AT Type 1 error")) + 
    geom_point(aes(y=sim.df[6], colour="CACE Type 1 error")) + 
    geom_line(aes(y = pred.t1.iv, colour="IV Type 1 error" ))+
    geom_line(aes(y = pred.t1.cace, colour="CACE Type 1 error"))+
    geom_line(aes(y = pred.t1.itt, colour="ITT Type 1 error"))+
    geom_line(aes(y = pred.t1.pp, colour="PP Type 1 error"))+
    geom_line(aes(y = pred.t1.at, colour="AT Type 1 error"))+
    geom_line(aes(y=alpha, colour='Alpha'),linetype="dotted") +
    xlab("Proportion of compliers")+
    ylab("Type 1 error")
  
  return(bias.plot)
}

t1<-simdata.type1(n=230,wna=1,pc1=0.4, pa=0.4, pn=0.4, NImargin=0.12, nIterations=1000)
t2<-simdata.type1(n=230,wna=2,pc1=0.4, pa=0.4, pn=0.4, NImargin=0.12, nIterations=1000)
t3<-simdata.type1(n=230,wna=3,pc1=0.4, pa=0.4, pn=0.4, NImargin=0.12, nIterations=1000)

t4<-simdata.type1(n=230,wna=1,pc1=0.2, pa=0.4, pn=0.4, NImargin=0.12, nIterations=1000)
t5<-t1
t6<-simdata.type1(n=230,wna=1,pc1=0.6,  pa=0.4, pn=0.4, NImargin=0.12, nIterations=1000)
t7<-simdata.type1(n=230,wna=1,pc1=0.3,  pa=0.4, pn=0.4, NImargin=0.12, nIterations=1000)
t8<-simdata.type1(n=230,wna=1,pc1=0.5,  pa=0.4,  pn=0.4, NImargin=0.12,nIterations=1000)
t9<-simdata.type1(n=230,wna=1,pc1=0.7,  pa=0.4,  pn=0.4, NImargin=0.12, nIterations=1000)
 
t10<-simdata.type1(n=230,wna=1,pc1=0.4, pa=0.6, pn=0.4, NImargin=0.12,nIterations=1000)
t11<-t1
t12<-simdata.type1(n=230,wna=1,pc1=0.4, pa=0.2, pn=0.4, NImargin=0.12, nIterations=1000)
t13<-simdata.type1(n=230,wna=1,pc1=0.4, pa=0.7, pn=0.5, NImargin=0.12, nIterations=1000)
t14<-simdata.type1(n=230,wna=1,pc1=0.4, pa=0.5, pn=0.5, NImargin=0.12, nIterations=1000)
t15<-simdata.type1(n=230,wna=1,pc1=0.4, pa=0.3, pn=0.5, NImargin=0.12, nIterations=1000)

ggarrange(t10, t11, t12, 
          t13, t14, t15,
          ncol = 3, nrow = 2)

##########POWER#########
########################
simdata.power<- function(n, wna, pc1, pc0, pa, pn, NImargin, nIterations){  
  # n: number of participants per group; 
  # wna: proportion of never takers; 
  # pc0: survival rate in compliers in intervention 0; 
  # pc1: survival rate in compliers in intervention 1;
  # pa: survival rate in always takers; 
  # pn: survival rate in never takers;
  # NImargin: non inferiority margin (negative);
  # nIterations: number of iterations 
  
  alpha=0.025
  
  #define proportions of patients
  wc=seq(0.1,1, by=0.05)         # proportion of compliers
  wa= (1-wc)/(wna+1)             # proportion of always takers 
  wn= wa*wna                     # proportion of never takers
  
  w00=wc+wn                      # Randomised to 0, intervention 0 (C+NT), proportion
  w01=rep(wn, length (wc))       # Randomised to 0, intervention 1 (NT), proportion
  w10= wa                        # Randomised to 1, intervention 0 (AT), proportion
  w11=wc+wa                      # Randomised to 1, intervention 1 (C+AT), proportion
  
  #define proportions of events  # simulation based on alternative hypothesis that pc1-pc0>NImargin (negative)
  p00=(wc*pc0+wn*pn)/(wc+wn)     # Randomised to 0, intervention 0 (C+NT), outcome
  p01=rep(pn, length(wc))        # Randomised to 0, intervention 1 (NT), outcome
  p10=rep(pa, length(wc))        # Randomised to 1, intervention 0 (AT), outcome
  p11=(wc*pc1+wa*pa)/(wc+wa)     # Randomised to 1, intervention 1 (C+AT), outcome
  
  #make up vectors for simulations 
  power.iv<-power.itt<- power.at<-power.pp<-power.cace<-c() 
  .power.iv<-.power.itt<- .power.at<-.power.pp<-.power.cace<-c() 
  
  #simulate and derive treatment effect 
  sim<- function() {
    for(i in 1:length(wc)) { print(i)
      for(l in 1:nIterations) { 
        #simulate data frame 
        simdata<- data.frame (                           
          "randomisation" = c(rep(1,n), rep(0,n)))     #half randomised to 1, half randomised to 0 
        
        simdata$intervention <- c(rbinom(n,1,w11[i]),   # for those randomised to 1, intervention proportion w11 
                                  rbinom(n,1,w10[i]))   # for those randomised to 0, intervention proportion w10
        
        simdata$outcome <- rep(NA,2*n)
        for(j in 1:(2*n)){
          if(simdata$randomisation[j]==1 & simdata$intervention[j]==1){
            simdata$outcome[j] <- rbinom(1,1,prob=p11[i])# for those randomised to 1, intervention 1, outcome proportion p11
          }  
          else if(simdata$randomisation[j]==1 & simdata$intervention[j]==0){
            simdata$outcome[j] <-rbinom(1,1,prob=p01[i]) # for those randomised to 1, intervention 0, outcome proportion p01
          }
          else if(simdata$randomisation[j]==0 & simdata$intervention[j]==1){
            simdata$outcome[j] <-rbinom(1,1,prob=p10[i]) # for those randomised to 0, intervention 1, outcome proportion p10
          }
          else if(simdata$randomisation[j]==0 & simdata$intervention[j]==0){
            simdata$outcome[j] <-rbinom(1,1,prob=p00[i]) # for those randomised to 0, intervention 0, outcome proportion p00
          }
        }
        
        #generate proportions in simulated data 
        w1.vector<-(filter(simdata,simdata$randomisation==1))$intervention; w1.value<-mean(w1.vector)
        w0.vector<-(filter(simdata,simdata$randomisation==0))$intervention; w0.value<-mean(w0.vector)
        wc.value<- w1.value-w0.value
        if (wc.value==0) {wc.value=wc[i]} 
        p11.vector<- (filter(simdata,intervention==1, simdata$randomisation==1))$outcome; p11.value<-mean(p11.vector)
        p00.vector<- (filter(simdata,intervention==0, simdata$randomisation==0))$outcome; p00.value<-mean(p00.vector)
        pz1.vector<- (filter(simdata, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
        pz0.vector<- (filter(simdata, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
        pd1.vector<- (filter(simdata, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
        pd0.vector<- (filter(simdata, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
        
        #estimate treatment effects 
        eff.itt<- pz0.value-pz1.value            #intention to treat 
        eff.cace<-  eff.itt/wc.value             #complier average causal effect 
        eff.pp <- p00.value-p11.value            #per protocol
        eff.at <- pd0.value-pd1.value            #as treated 
        ivmodel=ivreg(simdata$outcome ~ simdata$intervention, ~ simdata$randomisation, x=TRUE, data=simdata)
        eff.iv<--ivmodel$coef[2]                 #iv with 2 stage regression
        
        #variances of proportions and effects  
        var.eff.cace<-pc0*(1-pc0)/(w00[i]*n-w01[i]*n)+ pc1*(1-pc1)/(w11[i]*n-w10[i]*n)
        var.eff.itt<- pz1.value*(1-pz1.value)/length(pz1.vector) + pz0.value*(1-pz0.value)/length(pz0.vector)
        var.eff.pp<-  p11.value*(1-p11.value)/length(p11.vector) + p00.value*(1-p00.value)/length(p00.vector)
        var.eff.at<-  pd1.value*(1-pd1.value)/length(pd1.vector) + pd0.value*(1-pd0.value)/length(pd0.vector)
        if (is.na(eff.iv)) {var.eff.iv<-var.eff.cace} else {var.eff.iv<-robust.se(ivmodel)[2,2]^2 }
        
        #z values 
        z.itt<-(eff.itt-NImargin)/sqrt(var.eff.itt)
        z.iv<-(eff.iv-NImargin)/sqrt(var.eff.iv)
        z.pp<-(eff.pp-NImargin)/sqrt(var.eff.pp)
        z.at<-(eff.at-NImargin)/sqrt(var.eff.at)
        z.cace<-(eff.cace-NImargin)/sqrt(var.eff.cace)
        
        #power calculations 
        .power.iv[l]<- pnorm(z.iv-qnorm(1-alpha))+pnorm(-z.iv-qnorm(1-alpha))
        .power.itt[l]<- pnorm(z.itt-qnorm(1-alpha))+pnorm(-z.itt-qnorm(1-alpha))
        .power.at[l]<- pnorm(z.at-qnorm(1-alpha))+pnorm(-z.at-qnorm(1-alpha))
        .power.pp[l]<-pnorm(z.pp-qnorm(1-alpha))+pnorm(-z.pp-qnorm(1-alpha))
        .power.cace[l]<-pnorm(z.cace-qnorm(1-alpha))+pnorm(-z.cace-qnorm(1-alpha))
        
        power.iv[i]<-mean(.power.iv, na.rm=TRUE)
        power.itt[i]<-mean(.power.itt,na.rm=TRUE)
        power.pp[i]<-mean(.power.pp,na.rm=TRUE)
        power.at[i]<-mean(.power.at,na.rm=TRUE)
        power.cace[i]<- mean(.power.cace,na.rm=TRUE)
      }
    }
    return(data.frame(wc, power.iv, power.itt, power.pp, power.at, power.cace))
  }
  
  sim.df<-sim()
  
  sim.df$pred.p.iv<- predict(lm(power.iv~log(wc), data=sim.df))
  sim.df$pred.p.itt<- predict(lm(power.itt~log(wc), data=sim.df))
  sim.df$pred.p.at<- predict(lm(power.at~log(wc), data=sim.df))
  sim.df$pred.p.pp<- predict(lm(power.pp~log(wc), data=sim.df))
  sim.df$pred.p.cace1<- predict(lm(power.cace~log(wc), data=sim.df))
  sim.df$pred.p.cace<- sim.df$pred.p.cace1-sim.df$pred.p.cace1/50
  
  plot <- ggplot(sim.df, aes(sim.df[1]))+ 
    geom_point(aes(y=sim.df[2], colour="IV Power")) + 
    geom_point(aes(y=sim.df[3], colour="ITT Power")) + 
    geom_point(aes(y=sim.df[4], colour="PP Power")) + 
    geom_point(aes(y=sim.df[5], colour="AT Power")) + 
    geom_point(aes(y=sim.df[6], colour="CACE Power")) + 
    geom_line(aes(y = pred.p.iv, colour="IV Power" ))+
    geom_line(aes(y = pred.p.cace, colour="CACE Power"))+
    geom_line(aes(y = pred.p.itt, colour="ITT Power"))+
    geom_line(aes(y = pred.p.pp, colour="PP Power"))+
    geom_line(aes(y = pred.p.at, colour="AT Power"))+
    geom_line(aes(y=.8, colour='Target power'),linetype="dotted") +
    xlab("Proportion of compliers")+
    ylab("Power")
  
  return(plot)
}

p1<-simdata.power(n=230,wna=1, pc1=0.4, pa=0.4,pc0=0.4, pn=0.4, NImargin=-0.12, nIterations=1000) #pc1-pc0=0
p2<-simdata.power(n=230,wna=2,pc1=0.4, pa=0.4, pc0=0.4,pn=0.4, NImargin=-0.12, nIterations=1000)
p3<-simdata.power(n=230,wna=3,pc1=0.4, pa=0.4, pc0=0.4,pn=0.4, NImargin=-0.12, nIterations=1000)

p4<-simdata.power(n=230,wna=1,pc1=0.4, pa=0.4, pc0=0.3, pn=0.4, NImargin=-0.12, nIterations=1000) #pc1-pc0=0.1
p5<-p1                                                
p6<-simdata.power(n=230,wna=1,pc1=0.4,  pa=0.4, pc0=0.5,pn=0.4, NImargin=-0.12, nIterations=1000) #pc1-pc0= -0.1
p7<-simdata.power(n=230,wna=1,pc1=0.5,  pa=0.4, pc0=0.4,pn=0.4, NImargin=-0.12, nIterations=1000) #pc1-pc0=0.1
p8<-simdata.power(n=230,wna=1,pc1=0.5,  pa=0.4, pc0=0.5,pn=0.4, NImargin=-0.12, nIterations=1000) #pc1-pc0=0
p9<-simdata.power(n=230,wna=1,pc1=0.5,  pa=0.4, pc0=0.6,pn=0.4, NImargin=-0.12, nIterations=1000) #pc1-pc0=-0.1

p10<-simdata.power(n=230,wna=1,pc1=0.4, pa=0.6, pc0=0.4,pn=0.4, NImargin=-0.12,nIterations=1000) 
p11<-p1
p12<-simdata.power(n=230,wna=1,pc1=0.4, pa=0.2, pc0=0.4,pn=0.4, NImargin=-0.12, nIterations=1000) #pc1-pc0=0
p13<-simdata.power(n=230,wna=1,pc1=0.4, pa=0.7, pc0=0.4,pn=0.5, NImargin=-0.12, nIterations=1000)
p14<-simdata.power(n=230,wna=1,pc1=0.4, pa=0.5, pc0=0.4,pn=0.5, NImargin=-0.12, nIterations=1000)
p15<-simdata.power(n=230,wna=1,pc1=0.4, pa=0.3, pc0=0.4,pn=0.5, NImargin=-0.12, nIterations=1000)

ggarrange(p1, p2, p3, 
          #p7, p8, p9,
          ncol = 3, nrow = 1)

ggarrange(p4, p5, p6, 
          p7, p8, p9,
          ncol = 3, nrow = 2)

ggarrange(p10, p11, p12, 
          p13, p14, p15,
          ncol = 3, nrow = 2)
