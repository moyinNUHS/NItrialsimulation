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
  eff.iv<-eff.itt<-eff.pp<-eff.at<-c()                   
  .eff.iv<-.eff.itt<-.eff.pp<-.eff.at<-c()                    
  
  #simulate and derive treatment effect 
  sim<- function() { 
    for(i in 1:length(wc)) { 
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
        .eff.pp[l] <- p11.value-p00.value              #per protocol 
        .eff.at[l] <- pd1.value-pd0.value              # as treated 
        ivmodel<-ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=simdata) #iv with 2 stage regression
        .eff.iv[l]<-ivmodel$coef[2]                   
        
        #mean of iterated data 
        eff.itt[i]<- mean(.eff.itt, na.rm=TRUE)
        eff.pp[i] <- mean(.eff.pp, na.rm=TRUE) 
        eff.at[i] <- mean(.eff.at, na.rm=TRUE)
        eff.iv[i]<-mean(.eff.iv, na.rm = TRUE)
       
      }
    }
    return(data.frame(wc, eff.iv, eff.itt, eff.pp, eff.at))
  }
  
  sim.df<-sim()
  
  sim.df$pred.eff.iv<- predict(lm(eff.iv~log(wc), data=sim.df))
  sim.df$pred.eff.itt<- predict(lm(eff.itt~log(wc), data=sim.df))
  sim.df$pred.eff.at<- predict(lm(eff.at~log(wc), data=sim.df))
  sim.df$pred.eff.pp<- predict(lm(eff.pp~log(wc), data=sim.df))
  
  bias.plot <- ggplot(sim.df, aes(sim.df[1]))+ 
    geom_point(aes(y=sim.df[2], colour="IV effect")) + 
    geom_point(aes(y=sim.df[3], colour="ITT effect")) + 
    geom_point(aes(y=sim.df[4], colour="PP effect")) + 
    geom_point(aes(y=sim.df[5], colour="AT effect")) + 
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

ggarrange(bias1, bias2, bias3, ncol = 3, nrow = 1)

bias4<-bias1 #0
bias5<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.2, pa=0.4, pn=0.4, nIterations=1000) #-0.2
bias6<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.6, pa=0.4, pn=0.4, nIterations=1000) # 0.2
bias7<-simdata.bias(n=230,wna=1,pc1=0.3, pc0=0.4, pa=0.4, pn=0.4, nIterations=1000) # 0.1
bias8<-simdata.bias(n=230,wna=1,pc1=0.5, pc0=0.4, pa=0.4, pn=0.4, nIterations=1000) # -0.1
bias9<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.2, pn=0.4, nIterations=1000) # -0.2
bias10<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.6, pn=0.4, nIterations=1000) # 0.2
bias11<-simdata.bias(n=230,wna=1,pc1=0.4, pc0=0.4, pa=0.4, pn=0.3, nIterations=1000) # 0.1


ggarrange(bias4, bias5, bias6,  bias7,
          bias8, bias9, bias10, bias11,
          ncol = 4, nrow = 2)

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
  
  #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
  lower<-prop.test(x=c(pc1*n,(pc1-NImargin)*n),n=c(n,n), correct=FALSE) 
  lower.value<-lower$conf.int [1] #95% CI of NI margin 
  
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
  type1.error.iv<-type1.error.itt<-type1.error.pp<-type1.error.at<- c()       
  eff.iv<-eff.itt<- eff.pp<-eff.at<-c()
  
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
        p11.vector<- (filter(simdata,intervention==1, simdata$randomisation==1))$outcome; p11.value<-mean(p11.vector)
        p00.vector<- (filter(simdata,intervention==0, simdata$randomisation==0))$outcome; p00.value<-mean(p00.vector)
        pz1.vector<- (filter(simdata, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
        pz0.vector<- (filter(simdata, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
        pd1.vector<- (filter(simdata, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
        pd0.vector<- (filter(simdata, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
        
        #estimate treatment effects 
        eff.itt[l]<- pz1.value-pz0.value            #intention to treat 
        eff.pp [l]<- p11.value-p00.value            #per protocol
        eff.at [l]<- pd1.value-pd0.value            #as treated 
        ivmodel=ivreg(simdata$outcome ~ simdata$intervention, ~ simdata$randomisation, x=TRUE, data=simdata)
        eff.iv[l]<-ivmodel$coef[2]              #iv with 2 stage regression
        
      }
      
      #type 1 error
      type1.error.itt[i]<-mean(eff.itt<lower.value)
      type1.error.pp[i]<- mean(eff.pp<lower.value)
      type1.error.at[i]<- mean(eff.at<lower.value)
      type1.error.iv[i]<- mean(eff.iv<lower.value)

    }
    return(data.frame(wc, type1.error.iv, type1.error.itt, type1.error.pp, type1.error.at))
  }
  
  sim.df<-sim()
  
  sim.df$pred.t1.iv<- predict(lm(type1.error.iv~log(wc), data=sim.df))
  sim.df$pred.t1.itt<- predict(lm(type1.error.itt~log(wc), data=sim.df))
  sim.df$pred.t1.at<- predict(lm(type1.error.at~log(wc), data=sim.df))
  sim.df$pred.t1.pp<- predict(lm(type1.error.pp~log(wc), data=sim.df))
  
  bias.plot <- ggplot(sim.df, aes(sim.df[1]))+ 
    geom_point(aes(y=sim.df[2], colour="IV Type 1 error")) + 
    geom_point(aes(y=sim.df[3], colour="ITT Type 1 error")) + 
    geom_point(aes(y=sim.df[4], colour="PP Type 1 error")) + 
    geom_point(aes(y=sim.df[5], colour="AT Type 1 error")) + 
    geom_line(aes(y = pred.t1.iv, colour="IV Type 1 error" ))+
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

ggarrange(t1, t2, t3, ncol = 3, nrow = 1)

simdata.type1(n=500,wna=1,pc1=0.4, pa=0.4, pn=0.4, NImargin=0.12, nIterations=10)

t4<-t1 
t5<-simdata.type1(n=230,wna=1,pc1=0.2,pa=0.4,pn=0.4,NImargin=0.12,nIterations=1000)#-0.12
t6<-simdata.type1(n=230,wna=1,pc1=0.6,pa=0.4,pn=0.4,NImargin=0.12,nIterations=1000)#-0.12
t7<-simdata.type1(n=230,wna=1,pc1=0.4,pa=0.2,pn=0.4,NImargin=0.12,nIterations=1000)#-0.32
t8<-simdata.type1(n=230,wna=1,pc1=0.4,pa=0.6,pn=0.4,NImargin=0.12,nIterations=1000)#0.08
t9<-simdata.type1(n=230,wna=1,pc1=0.4,pa=0.7,pn=0.4,NImargin=0.12,nIterations=1000)#0.18
t10<-simdata.type1(n=230,wna=1,pc1=0.4,pa=0.4,pn=0.2,NImargin=0.12,nIterations=1000)#0.08
t11<-simdata.type1(n=230,wna=1,pc1=0.4,pa=0.4,pn=0.5,NImargin=0.12,nIterations=1000)#-0.22

ggarrange(t4, t5, t6, t7,
          t8, t9, t10,t11,
          ncol = 4, nrow = 2)

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
  
  #Lower critical value given p1 (0ne sided of 95% CI = alpha=0.025)
  lower<-prop.test(x=c(pc1*n,(pc1-NImargin)*n),n=c(n,n), correct=FALSE) 
  lower.value<-lower$conf.int [1] #95% CI of NI margin 
  
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
  power.iv<-power.itt<- power.at<-power.pp<-c() 
  eff.iv<-eff.itt<- eff.at<-eff.pp<-c() 
  
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
        p11.vector<- (filter(simdata,intervention==1, simdata$randomisation==1))$outcome; p11.value<-mean(p11.vector)
        p00.vector<- (filter(simdata,intervention==0, simdata$randomisation==0))$outcome; p00.value<-mean(p00.vector)
        pz1.vector<- (filter(simdata, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
        pz0.vector<- (filter(simdata, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
        pd1.vector<- (filter(simdata, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
        pd0.vector<- (filter(simdata, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
        
        #estimate treatment effects 
        eff.itt[l]<- pz1.value-pz0.value            #intention to treat 
        eff.pp [l]<- p11.value-p00.value            #per protocol
        eff.at[l] <- pd1.value-pd0.value            #as treated 
        ivmodel=ivreg(simdata$outcome ~ simdata$intervention, ~ simdata$randomisation, x=TRUE, data=simdata)
        eff.iv[l]<-ivmodel$coef[2]                 #iv with 2 stage regression
      
      }
      #power 
      power.iv[i]<- mean(eff.iv<lower.value)
      power.itt[i]<- mean(eff.itt<lower.value)
      power.at[i]<- mean(eff.at<lower.value)
      power.pp[i]<- mean(eff.pp<lower.value)
    }
    return(data.frame(wc, power.iv, power.itt, power.pp, power.at))
  }
  
  sim.df<-sim()
  
  sim.df$pred.p.iv<- predict(lm(power.iv~log(wc), data=sim.df))
  sim.df$pred.p.itt<- predict(lm(power.itt~log(wc), data=sim.df))
  sim.df$pred.p.at<- predict(lm(power.at~log(wc), data=sim.df))
  sim.df$pred.p.pp<- predict(lm(power.pp~log(wc), data=sim.df))
  
  plot <- ggplot(sim.df, aes(sim.df[1]))+ 
    geom_point(aes(y=sim.df[2], colour="IV Power")) + 
    geom_point(aes(y=sim.df[3], colour="ITT Power")) + 
    geom_point(aes(y=sim.df[4], colour="PP Power")) + 
    geom_point(aes(y=sim.df[5], colour="AT Power")) + 
    geom_line(aes(y = pred.p.iv, colour="IV Power" ))+
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

ggarrange(p1, p2, p3, ncol = 3, nrow = 1)

p4<- p1 
p5<-simdata.power(n=230,wna=1,pc1=0.4,pa=0.4,pc0=0.3,pn=0.4,NImargin=-0.12, nIterations=1000) #pc1-pc0=0.1                                               
p6<-simdata.power(n=230,wna=1,pc1=0.4,pa=0.4,pc0=0.5,pn=0.4,NImargin=-0.12, nIterations=1000) #pc1-pc0= -0.1
p7<-simdata.power(n=230,wna=1,pc1=0.4,pa=0.5,pc0=0.4,pn=0.3,NImargin=-0.12, nIterations=1000) #pc1-pc0=0
p8<-simdata.power(n=230,wna=1,pc1=0.5,pa=0.3,pc0=0.5,pn=0.5,NImargin=-0.12, nIterations=1000) #pc1-pc0=0
p9<- simdata.power(n=10000,wna=1,pc1=0.4,pa=0.4,pc0=0.5,pn=0.4,NImargin=-0.12, nIterations=10) #pc1-pc0=-0.1, can be overcome with larger sample size

ggarrange(p4, p5, p6, p7, p8, p9,
          ncol = 3, nrow = 2)

##########GROUP SEQUENTIAL#########
###################################
simdata.gs<- function(n, wna, pc1, pc0, pa, pn, NImargin, nIterations){  
  # n: number of participants per group:
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
  
  #make empty vectors for iterations 
  z.itt1<-z.pp1<-z.at1<-z.iv1<-z.itt2<-z.pp2<-z.at2<-z.iv2<-c()
  z.itt3<-z.pp3<-z.at3<-z.iv3<-z.itt4<-z.pp4<-z.at4<-z.iv4<-c()
  .z.itt<-.z.pp<-.z.at<-.z.iv<-c()

  #simulate and derive treatment effect 
  sim<- function() {
    for(i in 1:length(wc)) {  print (i)
      for(l in 1:nIterations) { print(l)
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
        
        #Shuffle rows 
        simdata <- simdata[sample(nrow(simdata)),]
        
        #data for each interim analysis 
        s1<-simdata[1:(2*n/4),]
        s2<-simdata[1:(2*n/4*2),]
        s3<-simdata[1:(2*n/4*3),]
        s4<-simdata[1:(2*n),]
        
        #Compute z values at each interim analysis (3 interim and 1 final) 
        #generate proportions in simulated data s1
        p11.vector<- (filter(s1,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
        p00.vector<- (filter(s1,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
        pz1.vector<- (filter(s1,randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
        pz0.vector<- (filter(s1,randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
        pd1.vector<- (filter(s1,intervention==1))$outcome;pd1.value<- mean(pd1.vector)
        pd0.vector<- (filter(s1,intervention==0))$outcome;pd0.value<- mean(pd0.vector)
        
        ###estimate treatment effects 
        eff.itt1[l]<- pz1.value-pz0.value            #intention to treat 
        eff.pp1[l] <- p11.value-p00.value            #per protocol
        eff.at1 [l]<- pd1.value-pd0.value            #as treated 
        ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s1)
        eff.iv1[l]<-ivmodel$coef[2]                 #iv with 2 stage regression
        
        #generate proportions in simulated data s2
        p11.vector<- (filter(s2,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
        p00.vector<- (filter(s2,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
        pz1.vector<- (filter(s2, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
        pz0.vector<- (filter(s2, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
        pd1.vector<- (filter(s2, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
        pd0.vector<- (filter(s2, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
        
        ###estimate treatment effects 
        eff.itt2[l]<- pz1.value-pz0.value            #intention to treat 
        eff.pp2[l] <- p11.value-p00.value            #per protocol
        eff.at2[l] <- pd1.value-pd0.value            #as treated 
        ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s2)
        eff.iv2[l]<-ivmodel$coef[2]                 #iv with 2 stage regression
        
        #generate proportions in simulated data s3
        p11.vector<- (filter(s3,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
        p00.vector<- (filter(s3,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
        pz1.vector<- (filter(s3, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
        pz0.vector<- (filter(s3, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
        pd1.vector<- (filter(s3, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
        pd0.vector<- (filter(s3, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
        
        ###estimate treatment effects 
        eff.itt3[l]<- pz1.value-pz0.value            #intention to treat 
        eff.pp3 [l]<- p11.value-p00.value            #per protocol
        eff.at3 [l]<- pd1.value-pd0.value            #as treated 
        ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s3)
        eff.iv3[l]<-ivmodel$coef[2]                 #iv with 2 stage regression
        
        #generate proportions in simulated data s4
        p11.vector<- (filter(s4,intervention==1, randomisation==1))$outcome; p11.value<-mean(p11.vector)
        p00.vector<- (filter(s4,intervention==0, randomisation==0))$outcome; p00.value<-mean(p00.vector)
        pz1.vector<- (filter(s4, randomisation==1))$outcome; pz1.value<-mean(pz1.vector)
        pz0.vector<- (filter(s4, randomisation==0))$outcome; pz0.value<-mean(pz0.vector) 
        pd1.vector<- (filter(s4, intervention==1))$outcome;pd1.value<- mean(pd1.vector)
        pd0.vector<- (filter(s4, intervention==0))$outcome;pd0.value<- mean(pd0.vector)
        
        ###estimate treatment effects 
        eff.itt4[l]<- pz1.value-pz0.value            #intention to treat 
        eff.pp4 [l]<- p11.value-p00.value            #per protocol
        eff.at4 [l]<- pd1.value-pd0.value            #as treated 
        ivmodel=ivreg(outcome ~ intervention, ~ randomisation, x=TRUE, data=s4)
        eff.iv4[l]<-ivmodel$coef[2]                 #iv with 2 stage regression
        
      }
      ###sds of effects s1
      sd.eff.itt1<-sd(eff.itt1)
      sd.eff.pp1<- sd(eff.pp1)
      sd.eff.at1<- sd(eff.at1)
      sd.eff.iv1<- sd(eff.iv1)
      
      ###z values s1
      z.itt1[i]<-(mean(eff.itt1)-NImargin)/sd.eff.itt1
      z.iv1[i]<-(mean(eff.iv1)-NImargin)/sd.eff.iv1
      z.pp1[i]<-(mean(eff.pp1)-NImargin)/sd.eff.pp1
      z.at1[i]<-(mean(eff.at1)-NImargin)/sd.eff.at1
      
      ###sds of effects s2
      sd.eff.itt2<-sd(eff.itt2)
      sd.eff.pp2<- sd(eff.pp2)
      sd.eff.at2<- sd(eff.at2)
      sd.eff.iv2<- sd(eff.iv2)
      
      ###z values s2
      z.itt2[i]<-(mean(eff.itt2)-NImargin)/sd.eff.itt2
      z.iv2[i]<-(mean(eff.iv2)-NImargin)/sd.eff.iv2
      z.pp2[i]<-(mean(eff.pp2)-NImargin)/sd.eff.pp2
      z.at2[i]<-(mean(eff.at2)-NImargin)/sd.eff.at2
      
      ###sds of effects s3
      sd.eff.itt3<-sd(eff.itt3)
      sd.eff.pp3<- sd(eff.pp3)
      sd.eff.at3<- sd(eff.at3)
      sd.eff.iv3<- sd(eff.iv3)
      
      ###z values s3
      z.itt3[i]<-(mean(eff.itt3)-NImargin)/sd.eff.itt3
      z.iv3[i]<-(mean(eff.iv3)-NImargin)/sd.eff.iv3
      z.pp3[i]<-(mean(eff.pp3)-NImargin)/sd.eff.pp3
      z.at3[i]<-(mean(eff.at3)-NImargin)/sd.eff.at3
      
      ###sds of effects s4
      sd.eff.itt4<-sd(eff.itt4)
      sd.eff.pp4<- sd(eff.pp4)
      sd.eff.at4<- sd(eff.at4)
      sd.eff.iv4<- sd(eff.iv4)
      
      ###z values s4
      z.itt4[i]<-(mean(eff.itt4)-NImargin)/sd.eff.itt4
      z.iv4[i]<-(mean(eff.iv4)-NImargin)/sd.eff.iv4
      z.pp4[i]<-(mean(eff.pp4)-NImargin)/sd.eff.pp4
      z.at4[i]<-(mean(eff.at4)-NImargin)/sd.eff.at4
      
      .z.itt<-list(z.itt1, z.itt2,z.itt3,z.itt4)
      .z.at<-list(z.at1, z.at2,z.at3,z.at4)
      .z.pp<-list(z.pp1,z.pp2,z.pp3,z.pp4)
      .z.iv<-list(z.iv1,z.iv2,z.iv3,z.iv4)
    } 
    
    z<-gsDesign(k=4, n.fix=2*n, delta0 = NImargin, n.I=c(nrow(s1),nrow(s2) ,nrow(s3),nrow(s4)), maxn.IPlan = 2*n, beta=0.2) #boundary 
    z<-z$upper$bound #4 boundaries 
    z.itt<-data.frame(.z.itt)
    z.at<-data.frame(.z.at)
    z.pp<-data.frame(.z.pp)
    z.iv<-data.frame(.z.iv)
    z.c<-data.frame(cbind(prop.c, z.itt, z.at, z.pp,z.iv, 
                          rep(z[1],length(wc)),
                          rep(z[2],length(wc)),
                          rep(z[3],length(wc)),
                          rep(z[4],length(wc))))
    names(z.c)<-c('c','itt1','itt2','itt3','itt4',
                  'at1','at2','at3','at4',
                  'pp1','pp2','pp3','pp4',
                  'iv1','iv2','iv3','iv4',
                  'z1','z2','z3','z4')
    
    return(z.c)
    
  }
  
  z.c<- sim()
  
  #plot Z values and boundaries against wc
  z.c$pred.z.itt1<- predict(lm(itt1~log(wc), data=z.c))
  z.c$pred.z.itt2<- predict(lm(itt2~log(wc), data=z.c))
  z.c$pred.z.itt3<- predict(lm(itt3~log(wc), data=z.c))
  z.c$pred.z.itt4<- predict(lm(itt4~log(wc), data=z.c))
  z.c$pred.z.at1<- predict(lm(at1~log(wc), data=z.c))
  z.c$pred.z.at2<- predict(lm(at2~log(wc), data=z.c))
  z.c$pred.z.at3<- predict(lm(at3~log(wc), data=z.c))
  z.c$pred.z.at4<- predict(lm(at4~log(wc), data=z.c))
  z.c$pred.z.pp1<- predict(lm(pp1~log(wc), data=z.c))
  z.c$pred.z.pp2<- predict(lm(pp2~log(wc), data=z.c))
  z.c$pred.z.pp3<- predict(lm(pp3~log(wc), data=z.c))
  z.c$pred.z.pp4<- predict(lm(pp4~log(wc), data=z.c))
  z.c$pred.z.iv1<- predict(lm(iv1~log(wc), data=z.c))
  z.c$pred.z.iv2<- predict(lm(iv2~log(wc), data=z.c))
  z.c$pred.z.iv3<- predict(lm(iv3~log(wc), data=z.c))
  z.c$pred.z.iv4<- predict(lm(iv4~log(wc), data=z.c))
  
  
  plot.z.itt<- ggplot(z.c, aes(z.c[1]))+ 
    geom_point(aes(y = z.c[2], colour="Z value ITT 1" ))+
    geom_point(aes(y = z.c[3], colour="Z value ITT 2"))+
    geom_point(aes(y = z.c[4], colour="Z value ITT 3"))+
    geom_point(aes(y = z.c[5], colour="Z value ITT 4"))+
    geom_line(aes(y = pred.z.itt1, colour="Z value ITT 1" ))+
    geom_line(aes(y = pred.z.itt2, colour="Z value ITT 2"))+
    geom_line(aes(y = pred.z.itt3, colour="Z value ITT 3" ))+
    geom_line(aes(y = pred.z.itt4, colour="Z value ITT 4"))+
    geom_line(aes(y = z1, colour="Z value ITT 1"),linetype="dotted")+
    geom_line(aes(y = z2, colour="Z value ITT 2"),linetype="dotted")+
    geom_line(aes(y = z3, colour="Z value ITT 3"),linetype="dotted")+
    geom_line(aes(y = z4, colour="Z value ITT 4"),linetype="dotted")+
    xlab("Proportion of compliers")+
    ylab("Z value")
  
  plot.z.at<- ggplot(z.c, aes(z.c[1]))+ 
    geom_point(aes(y = z.c[6], colour="Z value AT 1" ))+
    geom_point(aes(y = z.c[7], colour="Z value AT 2"))+
    geom_point(aes(y = z.c[8], colour="Z value AT 3"))+
    geom_point(aes(y = z.c[9], colour="Z value AT 4"))+
    geom_line(aes(y = pred.z.at1, colour="Z value AT 1" ))+
    geom_line(aes(y = pred.z.at2, colour="Z value AT 2"))+
    geom_line(aes(y = pred.z.at3, colour="Z value AT 3" ))+
    geom_line(aes(y = pred.z.at4, colour="Z value AT 4"))+
    geom_line(aes(y = z1, colour="Z value AT 1"),linetype="dotted")+
    geom_line(aes(y = z2, colour="Z value AT 2"),linetype="dotted")+
    geom_line(aes(y = z3, colour="Z value AT 3"),linetype="dotted")+
    geom_line(aes(y = z4, colour="Z value AT 4"),linetype="dotted")+
    xlab("Proportion of compliers")+
    ylab("Z value")
  
  plot.z.pp<- ggplot(z.c, aes(z.c[1]))+ 
    geom_point(aes(y = z.c[10], colour="Z value PP 1" ))+
    geom_point(aes(y = z.c[11], colour="Z value PP 2"))+
    geom_point(aes(y = z.c[12], colour="Z value PP 3"))+
    geom_point(aes(y = z.c[13], colour="Z value PP 4"))+
    geom_line(aes(y = pred.z.pp1, colour="Z value PP 1" ))+
    geom_line(aes(y = pred.z.pp2, colour="Z value PP 2"))+
    geom_line(aes(y = pred.z.pp3, colour="Z value PP 3" ))+
    geom_line(aes(y = pred.z.pp4, colour="Z value PP 4"))+
    geom_line(aes(y = z1, colour="Z value PP 1"),linetype="dotted")+
    geom_line(aes(y = z2, colour="Z value PP 2"),linetype="dotted")+
    geom_line(aes(y = z3, colour="Z value PP 3"),linetype="dotted")+
    geom_line(aes(y = z4, colour="Z value PP 4"),linetype="dotted")+
    xlab("Proportion of compliers")+
    ylab("Z value")
  
  plot.z.iv<- ggplot(z.c, aes(z.c[1]))+ 
    geom_point(aes(y = z.c[14], colour="Z value IV 1" ))+
    geom_point(aes(y = z.c[15], colour="Z value IV 2"))+
    geom_point(aes(y = z.c[16], colour="Z value IV 3"))+
    geom_point(aes(y = z.c[17], colour="Z value IV 4"))+
    geom_line(aes(y = pred.z.iv1, colour="Z value IV 1" ))+
    geom_line(aes(y = pred.z.iv2, colour="Z value IV 2"))+
    geom_line(aes(y = pred.z.iv3, colour="Z value IV 3" ))+
    geom_line(aes(y = pred.z.iv4, colour="Z value IV 4"))+
    geom_line(aes(y = z1, colour="Z value IV 1"),linetype="dotted")+
    geom_line(aes(y = z2, colour="Z value IV 2"),linetype="dotted")+
    geom_line(aes(y = z3, colour="Z value IV 3"),linetype="dotted")+
    geom_line(aes(y = z4, colour="Z value IV 4"),linetype="dotted")+
    xlab("Proportion of compliers")+
    ylab("Z value")
  
  plot<- ggarrange(plot.z.itt, plot.z.at, plot.z.pp, plot.z.iv, ncol = 2, nrow = 2)
  
  return(plot)
}

gs1<-simdata.gs(n=230,wna=1,pc0=0.6, pc1=0.6,pa=0.6,pn=0.6,nIterations=100,NImargin=-.12)  
gs2<-simdata.gs(n=300,wna=1,pc0=0.6, pc1=0.6,pa=0.6,pn=0.6,nIterations=100,NImargin=-.12) # increased sample size
gs3<-simdata.gs(n=230,wna=1,pc0=0.6, pc1=0.7,pa=0.6,pn=0.6,nIterations=100,NImargin=-.12) # increased effect of intervention 
gs4<-simdata.gs(n=230,wna=3,pc0=0.6, pc1=0.6,pa=0.6,pn=0.6,nIterations=100,NImargin=-.12) # increased effect wna ratio
gs5<-simdata.gs(n=230,wna=1,pc0=0.6, pc1=0.6,pa=0.8,pn=0.6,nIterations=100,NImargin=-.12) # increased pa
gs6<-simdata.gs(n=230,wna=1,pc0=0.6, pc1=0.6,pa=0.6,pn=0.8,nIterations=100,NImargin=-.12) # increased pn


#using gsDesign package to determine sample size for REGARD-VAP 
n.fix<- nBinomial (p1=0.6, p2=0.6, delta0=-0.12)
x<-gsDesign(n.fix=n.fix, k=4, beta=0.2, delta0=-0.12)
ceiling(x$n.I)
y<-gsProbability(theta=x$delta*seq(0,2,0.25),d=x)
plot(y, plottype=6, lty=2, lwd=3)
