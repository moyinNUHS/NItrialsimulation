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

######################################BIAS######################################
################################################################################

simdata.bias<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
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
    geom_point(aes(y=sim.df[2], colour="IV effect",  alpha=0.4)) + 
    geom_point(aes(y=sim.df[3], colour="ITT effect", alpha=0.4)) + 
    geom_point(aes(y=sim.df[4], colour="PP effect",  alpha=0.4)) + 
    geom_point(aes(y=sim.df[5], colour="AT effect",  alpha=0.4)) + 
    geom_point(aes(y=sim.df[6], colour="Modified PP effect", alpha=0.4)) + 
    geom_point(aes(y=sim.df[7], colour="Propensity score effect", alpha=0.4)) + 
    geom_smooth(aes(y=sim.df[2], colour="IV effect"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=sim.df[3], colour="ITT effect"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=sim.df[4], colour="PP effect"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=sim.df[5], colour="AT effect"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=sim.df[6], colour="Modified PP effect"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_smooth(aes(y=sim.df[7], colour="Propensity score effect"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    geom_line(aes(y=sim.df[interval,4], colour='True effect'),linetype="dotted") +
    scale_alpha_manual(values = c("MS"=0.2, "MT"=0.2, "TLR" = 1), guide = 'none') +
    xlab("Proportion of compliant participants")+
    ylab("Effect")+
    theme_bw()+
    scale_colour_manual(values=cbPalette)+
    theme(legend.position="bottom")+
    theme(legend.title=element_blank())+
    scale_x_continuous(limits=c(0.60, 1))+
    scale_y_continuous(limits=c(-0.2, 0))
  
  return(bias.plot)
} 

######BIAS: Varying effect of known confounder on outcome (without unknown confounder)########
#use scale -0.2 to 0 
bias1<-simdata.bias(n=300, p1=0.3, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
bias2<-simdata.bias(n=300, p1=0.3, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = -0.2, confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
bias3<-simdata.bias(n=300, p1=0.3, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 1000000,  confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
bias4<-simdata.bias(n=300, p1=0.3, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 1000000,  confounder.eff.o = -0.2, confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
bias5<-simdata.bias(n=300, p1=0.3, p0=0.4, i0=0.45, i1=0.55, confounder.eff.i= 0.01,     confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
bias6<-simdata.bias(n=300, p1=0.3, p0=0.4, i0=0.45, i1=0.55, confounder.eff.i= 0.01,     confounder.eff.o = -0.2, confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
biasfig1a<-ggarrange(bias1, bias3, ncol = 1, nrow = 2,labels = c(paste(LETTERS[1:2], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(biasfig1a, top = text_grob("Scenario 1: Non-compliance in group assigned to intervention ", face = "bold", size = 12))
biasfig1b<-ggarrange(bias2, bias4, ncol = 1, nrow = 2,labels = c(paste(LETTERS[1:2], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(biasfig1b, top = text_grob("Scenario 2: Non-compliance in group assigned to standard-of-care", face = "bold", size = 12))
biasfig1c<-ggarrange(bias5, bias6, ncol = 1, nrow = 2,labels = c(paste(LETTERS[1:2], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(biasfig1c, top = text_grob("Scenario 3: Non-compliance in both groups", face = "bold", size = 12))

######BIAS: Varying effect of unknown confounder on outcome########
bias7<-simdata.bias(n=300,  p1=0.3, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 0.001,  confounder.eff.o = 0.3,confounder.u.eff.i=0.001, confounder.u.eff.o=0.2,nIterations=1000)
bias8<-simdata.bias(n=300,  p1=0.3, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 1000,  confounder.eff.o = 0.3,confounder.u.eff.i=1000,  confounder.u.eff.o=0.2,nIterations=1000)
bias9<-simdata.bias(n=300,  p1=0.3, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 0.001,  confounder.eff.o = 0.5,confounder.u.eff.i=0.001,confounder.u.eff.o=-0.2,nIterations=1000)
bias10<-simdata.bias(n=300, p1=0.3, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 1000,   confounder.eff.o = 0.5,confounder.u.eff.i=1000, confounder.u.eff.o=-0.2,nIterations=1000)
biasfig2<-ggarrange(bias7, bias8, bias9, bias10, ncol = 2, nrow = 2,labels = c(paste(LETTERS[1:4], sep="")), common.legend = TRUE, legend = "bottom")
annotate_figure(biasfig2, top = text_grob("Figure 3b. Varying effect of unknown confounder", face = "bold", size = 12))

######BIAS: Varying event rates########
#change scale to -2 to +0.2
bias11<-bias1
bias12<-simdata.bias(n=300, p1=0.4, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
bias13<-simdata.bias(n=300, p1=0.5, p0=0.4, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
biasfig3<-ggarrange(bias11, bias12, bias13, ncol = 3, nrow = 1,labels = c(paste(LETTERS[1:3], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(biasfig3, top = text_grob("Figure 3c. Varying effect of event rates in intervention and control groups", face = "bold", size = 12))

##################################TYPE 1 ERROR##################################
################################################################################
simdata.type1<- function(n, i0, i1, NImargin, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations){  
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
  
  if (NImargin0) {
    
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
      geom_point(aes(y=sim.df[2], colour="IV type 1 error", alpha=0.4)) + 
      geom_point(aes(y=sim.df[3], colour="ITT type 1 error",alpha=0.4)) + 
      geom_point(aes(y=sim.df[4], colour="PP type 1 error", alpha=0.4)) + 
      geom_point(aes(y=sim.df[5], colour="AT type 1 error", alpha=0.4)) + 
      geom_point(aes(y=sim.df[6], colour="Modified PP type 1 error", alpha=0.4)) +
      geom_point(aes(y=sim.df[7], colour="Propensity score type 1 error", alpha=0.4)) + 
      geom_smooth(aes(y=sim.df[2], colour="IV type 1 error"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[3], colour="ITT type 1 error"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[4], colour="PP type 1 error"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[5], colour="AT type 1 error"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[6], colour="Modified PP type 1 error"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[7], colour="Propensity score type 1 error"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_line(aes(y=0.025, colour='True type 1 error'),linetype="dotted") +
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

###########Varying effect of known confounder########
t1<-simdata.type1(n=200, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t2<-simdata.type1(n=200, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = -0.2, confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t3<-simdata.type1(n=200, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000,  confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t4<-simdata.type1(n=200, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000,  confounder.eff.o = -0.2, confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t5<-simdata.type1(n=200, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.01,     confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t6<-simdata.type1(n=200, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.01,     confounder.eff.o = -0.2, confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
tfig1<-ggarrange(t1, t2, t3, t4,t5, t6, ncol = 3, nrow = 2,labels = c(paste(LETTERS[1:6], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(tfig1, top = text_grob("Figure 4a. Varying effect of known confounder on outcome", face = "bold", size = 12))

#########Varying effect of unknown confounder###########
t7<-simdata.type1(n=200, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.001, confounder.eff.o = 0.3,  confounder.u.eff.i=0.001, confounder.u.eff.o=0.2, nIterations=1000)
t8<-simdata.type1(n=200, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000,  confounder.eff.o = 0.3,  confounder.u.eff.i=1000,  confounder.u.eff.o=0.2, nIterations=1000)
t9<-simdata.type1(n=200, p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.001,  confounder.eff.o = 0.5,  confounder.u.eff.i=0.001, confounder.u.eff.o=-0.2,nIterations=1000)
t10<-simdata.type1(n=200,p1=0.4, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000,  confounder.eff.o = 0.5,  confounder.u.eff.i=1000,  confounder.u.eff.o=-0.2,nIterations=1000)
tfig2<-ggarrange(t7, t8, t9, t10, ncol = 2, nrow = 2,labels = c(paste(LETTERS[1:4], sep="")), common.legend = TRUE, legend = "bottom")
annotate_figure(tfig2, top = text_grob("Figure 4b. Varying effect of unknown confounder", face = "bold", size = 12))

###########Varying outcome rate###########
t11<-simdata.type1(n=200, p1=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1, confounder.u.eff.o=0,nIterations=1000)
t12<-t2
t13<-simdata.type1(n=200, p1=0.5, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1, confounder.u.eff.o=0,nIterations=1000)
tfig3<-ggarrange(t11, t12, t13, ncol = 3, nrow = 1,labels = c(paste(LETTERS[1:3], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(tfig3, top = text_grob("Figure 4c. Varying effect of event rates in intervention and control groups", face = "bold", size = 12))

#########Varying sample size ###############
t14<-simdata.type1(n=300, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.01,confounder.eff.o = 0.5,confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
t15<-t5
t16<-simdata.type1(n=400, p1=0.4, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.01,confounder.eff.o = 0.5,confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
tfig4<-ggarrange(t14, t15, t16, ncol = 3, nrow = 1,labels = c(paste(LETTERS[1:3], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(tfig4, top = text_grob("Figure 4d. Varying sample size", face = "bold", size = 12))

######################################POWER#####################################
################################################################################
simdata.power<- function(n, i0, i1, p1, p0, confounder.eff.o, confounder.eff.i, confounder.u.eff.i,confounder.u.eff.o, NImargin, nIterations){  
  
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
      geom_point(aes(y=sim.df[2], colour="IV power", alpha=0.4)) + 
      geom_point(aes(y=sim.df[3], colour="ITT power", alpha=0.4)) + 
      geom_point(aes(y=sim.df[4], colour="PP power", alpha=0.4)) + 
      geom_point(aes(y=sim.df[5], colour="AT power", alpha=0.4)) + 
      geom_point(aes(y=sim.df[6], colour="Modified PP power", alpha=0.4)) +
      geom_point(aes(y=sim.df[7], colour="Propensity score power", alpha=0.4)) + 
      geom_smooth(aes(y=sim.df[2], colour="IV power"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[3], colour="ITT power"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[4], colour="PP power"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[5], colour="AT power"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[6], colour="Modified PP power"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_smooth(aes(y=sim.df[7], colour="Propensity score power"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
      geom_line(aes(y=.8, colour='Target power'),linetype="dotted") +
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

############Varying effect of known confounder########
p1<-simdata.power(n=200, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
p2<-simdata.power(n=200, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = -0.2, confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
p3<-simdata.power(n=200, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000,  confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
p4<-simdata.power(n=200, p1=0.4, p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000000,  confounder.eff.o = -0.2, confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
p5<-simdata.power(n=200, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.01,     confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
p6<-simdata.power(n=200, p1=0.4, p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.01,     confounder.eff.o = -0.2, confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
pfig1<-ggarrange(p1, p2, p3, p4,p5, p6, ncol = 3, nrow = 2,labels = c(paste(LETTERS[1:6], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(pfig1, top = text_grob("Figure 5a. Varying effect of known confounder on outcome", face = "bold", size = 12))

#########Varying effect of unknown confounder###########
p7<-simdata.power(n=200, p1=0.4, p0=0.3,NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.0001, confounder.eff.o = 0.3,  confounder.u.eff.i=0.001, confounder.u.eff.o=0.2,nIterations=1000)
p8<-simdata.power(n=200, p1=0.4, p0=0.3,NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000,  confounder.eff.o = 0.3,  confounder.u.eff.i=1000,  confounder.u.eff.o=0.2,nIterations=1000)
p9<-simdata.power(n=200, p1=0.4, p0=0.3,NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.0001,   confounder.eff.o = 0.5,  confounder.u.eff.i=0.001,confounder.u.eff.o=-0.2,nIterations=1000)
p10<-simdata.power(n=200,p1=0.4, p0=0.3,NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 1000,   confounder.eff.o = 0.5,  confounder.u.eff.i=1000, confounder.u.eff.o=-0.2,nIterations=1000)
pfig2<-ggarrange(p7, p8, p9, p10, ncol = 2, nrow = 2,labels = c(paste(LETTERS[1:4], sep="")), common.legend = TRUE, legend = "bottom")
annotate_figure(pfig2, top = text_grob("Figure 5b. Varying effect of unknown confounder", face = "bold", size = 12))

###########Varying outcome rate###########
p11<-simdata.power(n=200, p1=0.3,p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1, confounder.u.eff.o=0,nIterations=1000)
p12<-t3
p13<-simdata.power(n=200, p1=0.5,p0=0.3, NImargin=0.12, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1, confounder.u.eff.o=0,nIterations=1000)
pfig3<-ggarrange(p11, p12, p13, ncol = 3, nrow = 1,labels = c(paste(LETTERS[1:3], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(pfig3, top = text_grob("Figure 5c. Varying effect of event rates in intervention and control groups", face = "bold", size = 12))

#########Varying sample size ###############
p14<-simdata.power(n=100, p1=0.4,p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.01, confounder.eff.o = 0.5,   confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
p15<-t5
p16<-simdata.power(n=300, p1=0.4,p0=0.3, NImargin=0.12, i0=0.45, i1=0.55, confounder.eff.i= 0.01, confounder.eff.o = 0.5,   confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000)
pfig4<-ggarrange(p14, p15, p16, ncol = 3, nrow = 1,labels = c(paste(LETTERS[1:3], sep="")), common.legend = TRUE, legend = "bottom" )
annotate_figure(pfig4, top = text_grob("Figure 5d. Varying sample size", face = "bold", size = 12))

######################################GROUP SEQUENTIAL#####################################
###########################################################################################
simdata.gs<- function(n, p1, p0, i0, i1, NImargin, confounder.eff.o, confounder.eff.i,confounder.u.eff.i,confounder.u.eff.o, nIterations){  

  if ((NImargin <0) & (p1-p0 > NImargin)) { #built with alternative hypothesis= pshort - plong > - NI 
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
            
            #shuffle rows in simulated data 
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
      
      plot.z.itt<- ggplot(z.c, aes(z.c[1]))+ 
        geom_point(aes(y = z.c[2], colour="Z value 1" , alpha=0.4))+
        geom_point(aes(y = z.c[3], colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = z.c[4], colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = z.c[5], colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y= z.c[2], colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[3], colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[4], colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[5], colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
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
        geom_point(aes(y = z.c[6], colour="Z value 1" , alpha=0.4))+
        geom_point(aes(y = z.c[7], colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = z.c[8], colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = z.c[9], colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y=z.c[6], colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[7], colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[8], colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[9], colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
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
        geom_point(aes(y = z.c[10], colour="Z value 1", alpha=0.4 ))+
        geom_point(aes(y = z.c[11], colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = z.c[12], colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = z.c[13], colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y=z.c[10], colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[11], colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[12], colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[13], colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
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
        geom_point(aes(y = z.c[14], colour="Z value 1" , alpha=0.4))+
        geom_point(aes(y = z.c[15], colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = z.c[16], colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = z.c[17], colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y=z.c[14], colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[15], colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[16], colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[17], colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
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
        geom_point(aes(y = z.c[18], colour="Z value 1" , alpha=0.4))+
        geom_point(aes(y = z.c[19], colour="Z value 2", alpha=0.4))+
        geom_point(aes(y = z.c[20], colour="Z value 3", alpha=0.4))+
        geom_point(aes(y = z.c[21], colour="Z value 4", alpha=0.4))+
        geom_smooth(aes(y=z.c[18], colour="Z value 1"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[19], colour="Z value 2"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[20], colour="Z value 3"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
        geom_smooth(aes(y=z.c[21], colour="Z value 4"), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
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
        scale_x_continuous(limits=c(0.55, 1))+
        scale_y_continuous(limits=c(0.5, 4))
      
      plot<- ggarrange(plot.z.itt, plot.z.at, plot.z.pp, plot.z.mpp,plot.z.ps,plot.z.iv, ncol = 3, nrow = 2,
                       labels = c("Intention to treat", "As treated","Per protocol","Modified PP","Propensity score","Instrumental variable"), common.legend = TRUE, legend = "bottom")
      
      return(plot)
    }else (print ("Number of iterations must be more than 1."))
  } else (print ("NI margin must be negative, and p1-p0 (in terms of unfavorable outcomes) must be more than NI margin in this simulation"))
}

gs1<-simdata.gs(n=200,p0=0.3, p1=0.4, i0=0.01, i1=0.99, confounder.eff.i= 0.000001, confounder.eff.o = 0.5,  confounder.u.eff.i=1,confounder.u.eff.o=0,nIterations=1000, NImargin=-0.12)
annotate_figure(gs1, top = text_grob("Figure 1. Varying effect of unknown confounder", face = "bold", size = 12))

gs2<-simdata.gs(n=200,p0=0.3, p1=0.4,i0=0.01, i1=0.99,confounder.eff.o=-0.2,nIterations=1000,confounder.u.eff.i= 0.05 ,confounder.u.eff.o=0.05,NImargin=-.12) 
annotate_figure(gs2, top = text_grob("Figure 2. Varying effect of unknown confounder", face = "bold", size = 12))

gs3<-simdata.gs(n=200,p0=0.6, p1=0.6,confounder.eff.o=10,nIterations=250,confounder.u.eff.i= 0.05 ,confounder.u.eff.o=0.05,NImargin=-.12) 
annotate_figure(gs3, top = text_grob("Figure 3. Varying effect of unknown confounder", face = "bold", size = 12))


gs4<-simdata.gs(n=300,p0=0.6, p1=0.6,confounder.eff.o=20,nIterations=250,confounder.u.eff.i= 0.05 ,confounder.u.eff.o=0.05,NImargin=-.12) 
annotate_figure(gs4, top = text_grob("Figure 4. Varying effect of unknown confounder", face = "bold", size = 12))

gs5<-simdata.gs(n=300,p0=0.4, p1=0.6,confounder.eff.o=20,nIterations=250,confounder.u.eff.i= 0.05 ,confounder.u.eff.o=0.05,NImargin=-.12) 
annotate_figure(gs5, top = text_grob("Figure 5. Varying effect of unknown confounder", face = "bold", size = 12))

gs6<-simdata.gs(n=300,p0=0.6, p1=0.6,confounder.eff.o=20,nIterations=250,confounder.u.eff.i= 0.05 ,confounder.u.eff.o=0.05,NImargin=-.12) 
annotate_figure(gs6, top = text_grob("Figure 6. Varying effect of unknown confounder", face = "bold", size = 12))

gs7<-simdata.gs(n=400,p0=0.4, p1=0.6,confounder.eff.o=20,nIterations=250,confounder.u.eff.i= 0.05 ,confounder.u.eff.o=0.05,NImargin=-.12) 
annotate_figure(gs7, top = text_grob("Figure 7. Varying effect of unknown confounder", face = "bold", size = 12))

################################################################
#using gsDesign package to determine sample size for REGARD-VAP 
n.fix<- nBinomial (p1=0.6, p2=0.6, delta0=-0.12)
x<-gsDesign(n.fix=n.fix, k=4, beta=0.2, delta0=-0.12)
ceiling(x$n.I)
y<-gsProbability(theta=x$delta*seq(0,2,0.25),d=x)
plot(y, plottype=6, lty=2, lwd=3)

################################################################
#plot type 1 and type 2 errors
library(grDevices)
par(mfrow=c(1,2))

#type 1 error 
x <- seq(-4, 4, length=1000) # x-axis 
hx <- dnorm(x, mean=2, sd=1) # normal distribution of null hypotehsis 
plot(x, hx, type="n", xlim=c(-4, 4), ylim=c(0, 0.5), #plot null hypothesis first 
     ylab = "",
     xlab = "",
     main= expression(paste("a. Type I (", alpha, ") error")), axes=FALSE)
axis(1, at = c(-4,-qnorm(.025), 0.12,4), #label x axis 
     labels = expression("","NI margin","p-value",""))

# Print null hypothesis area
col_null = "#AAAAAA"
polygon(c(min(x), x,max(x)), c(0,hx,0), col=col_null, lwd=2, density=c(10, 40), angle=-45, border=0)
lines(x, hx, lwd=2, lty="dashed", col=col_null)

## The area where the alpha is
col3 = adjustcolor("#D55E00", alpha.f=0.2)
rect(-4, 0,0.1, max(hx),border = NA, col=col3)

abline(v=0.1, lwd=2, col="#000088", lty="dashed")

legend("top", inset=.015, title="Legend",
       c("Null hypothesis", "Type I error"), fill=c(col_null, col3), 
       angle=-45,
       density=c(20, 1000, 1000), horiz=FALSE)


#type 2 error 
x <- seq(-6, 4, length=1000) # x-axis 
hx <- dnorm(x, mean=0, sd=1) # normal distribution of null hypotehsis 
plot(x, hx, type="n", xlim=c(-6, 4), ylim=c(0, 0.5), #plot null hypothesis first 
     ylab = "",
     xlab = "",
     main= expression(paste("b. Type II (", beta, ") error")), axes=FALSE)
axis(1, at = c(-6, qnorm(.025), 0, 4), #label x axis 
     labels = expression(-infinity,"p-value","NI margin" , infinity))

# Print null hypothesis area
col_null = "#AAAAAA"
polygon(c(min(x), x,max(x)), c(0,hx,0), col=col_null, lwd=2, density=c(10, 40), angle=-45, border=0)
lines(x, hx, lwd=2, lty="dashed", col=col_null)

#The alternative hypothesis area
shift = qnorm(0.025, mean=0, sd=1)*1.5 #shift alternative hypothesis graph to the left 
xfit2 <- x + shift
yfit2 <- dnorm(xfit2, mean=shift, sd=1)

## The underpowered area
lb <- min(xfit2)
ub <- round(qnorm(0.025),2)
col1 = adjustcolor("#0072B2", alpha.f=0.2)
i <- xfit2 >= lb & xfit2 <= ub
polygon(c(lb,xfit2[i],ub), c(0,yfit2[i],0), col=col1)

## The area where the power is
i <- xfit2 >= ub
col2 = adjustcolor("#F0E442", alpha.f=0.2)
polygon(c(ub,xfit2[i],max(xfit2)), c(0,yfit2[i],0), col=col2)

# Outline the alternative hypothesis
lines(xfit2, yfit2, lwd=2)

abline(v=ub, lwd=2, col="#000088", lty="dashed")

legend("top", inset=.015, title="Legend",
       c("Null hypothesis","Power", "Type II error"), fill=c(col_null, col1, col2), 
       angle=-45,
       density=c(20, 1000, 1000), horiz=FALSE)

################################################################
#for checking proportions
simdata<- function(r0, r1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o){  
  
  #make up vectors for simulations 
  .i1<-.i0<-.pp1<-.pp0<-.pp<-c()
  i1p<-i0p<-pp1p<-pp0p<-ppp<-c()
  .pp.out1<-.pp.out0<-.at.out1<-.at.out0<-.itt.out1<-.ittout0<-c()
  ppout1<-ppout0<-atout0<- atout1<- ittout0<- ittout1<-c()
  
  #number of data points in simulation 
  interval<-4 
  n<-300
  
  #simulate and derive treatment effect 
  sim<- function() { 
    for(i in 1:interval) { print(paste ("interval value",i))
      for(l in 1:100) { 
        #simulate data frame 
        
        #RANDOMISATION ratio 1:1
        simdata<- data.frame ("randomisation" = c(rep(1,n), rep(0,n)))     #half randomised to 1, half randomised to 0 
        simdata$id<- seq(1,(2*n), by=1) #create participant id 
        
        #confounder normal distribution ranging 0-1 
        simdata$confounder<- rep(NA,2*n) 
        simdata$confounder<- rnorm(2*n) # assume counfounder to be SOFA(severity), normally distributed 
        simdata$confounder<- sample((simdata$confounder - min(simdata$confounder))/(max(simdata$confounder) - min(simdata$confounder)))
        
        
        #unknown confounder normal distribution ranging 0-1 
        simdata$confounder.u<- rep(NA,2*n) 
        simdata$confounder.u<- rnorm(2*n) # assume counfounder to be normally distributed 
        simdata$confounder.u<- sample((simdata$confounder.u - min(simdata$confounder.u))/(max(simdata$confounder.u) - min(simdata$confounder.u)))
        
        #INTERVENTION dependent on randomisation and confounders
        simdata$intervention <- rep(NA,2*n)
        
        odd0<-r0/(1-r0)                                       # Number of intervention=1 when C=0 R=0/ Number of intervention=0 when C=0 R=0
        
        odd1<-confounder.eff.i                                # (Number of intervention=1 when C=1 / Number of intervention=0 when C=1) OVER
        # (Number of of intervention =1 when C=0 / Number of intervention =0 when C=0 ) (controlled for R)
        # if confounder has no negative effect on compliance - (1/1)/(1/1)
        # if confounder has strong negative effect on compliance - (1/10)/(10/1)
        odd2<-confounder.u.eff.i                              # (Number of intervention=1 when UC=1 / Number of intervention=0 when UC=1) OVER
        # (Number of of intervention =1 when UC=0 / Number of intervention=0 when UC=0 ) (controlled for R)
        # Strength of unknown confounder on intervention
        odd3<-(r1/(1-r1))/(r0/(1-r0))                         # (Number of intervention=1 when R=1 / Number of intervention=0 when R=1) OVER
        # (Number of of intervention =1 when R=0 / Number of of intervention =0 when R=0 ) (controlled for C)
        # randomisation is assumed to have more effect on compliance than confounder 
        
        .b0<-log (odd0)
        b0<- seq(.b0,log(0.01/(1-0.01)), length.out = interval)
        .b1<-log (odd1)
        b1<- seq(.b1,0, length.out = interval)
        .b2<-log (odd2)
        b2<- seq(.b2,0, length.out = interval)
        .b3<-log (odd3)
        b3<- seq(.b3,log((0.99/(1-0.99))/(0.01/(1-0.01))), length.out = interval)
        
        if (b1[i]==0) { pi<- simdata$randomisation} else {
          logit.pi<- b0[i]+b1[i]*simdata$confounder+b2[i]*simdata$confounder.u+b3[i]*simdata$randomisation
          pi<-exp(logit.pi)/(1+exp(logit.pi))
        }
        simdata$intervention <- rbinom(2*n,1,prob=pi)
        
        #OUTCOME dependent on confounders and intervention
        simdata$outcome <- rep(NA,2*n)
        
        odd0<-p0/(1-p0)                                       # Number of outcome=1 when C=0 I=0/ Number of outcome=0 when C=0 I=0 (mortality and VAP recurrences in ICU)
        odd1<-confounder.eff.o                                # (Number of outcome=1 when C=1 / Number of outcome=0 when C=1) OVER
        # (Number of of outcome =1 when C=0 / Number of outcome =0 when C=0 ) (controlled for R)
        # Strength of confounder on outcome
        odd2<-confounder.u.eff.o                              # (Number of outcome=1 when UC=1 / Number of outcome=0 when UC=1) OVER
        # (Number of of outcome =1 when UC=0 / Number of outcome=0 when UC=0 ) (controlled for R)
        # Strength of unknown confounder on outcome
        odd3<-(p1/(1-p1))/(p0/(1-p0))                         # (Number of outcome=1 when I=1 / Number of outcome=0 when I=1) OVER
        # (Number of of outcome =1 when I=0 / Number of outcome=0 when I=0 ) (controlled for C)
        
        b0<-log (odd0)
        b1<-log (odd1)
        b2<-log (odd2)
        b3<-log (odd3)
        
        logit.pi<- b0+b1*simdata$confounder+b2*simdata$confounder.u+b3*simdata$intervention
        pi<-exp(logit.pi)/(1+exp(logit.pi))
        simdata$outcome <- rbinom(2*n,1,prob=pi)
        
        #proportion of patients who were compliant to protocol 
        i1<- filter(simdata, intervention==1)
        i0<- filter(simdata, intervention == 0)
        pp1<- filter(simdata, randomisation==1 & intervention==1)
        pp0<-filter(simdata, randomisation==0 & intervention==0)
        pp<- rbind(pp1,pp0)
        
        .i1[l]<- nrow(i1)/(2*n)
        .i0[l]<- nrow(i0)/(2*n)
        .pp1[l]<- nrow(pp1)/(n)
        .pp0[l]<-nrow(pp0)/(n)
        .pp[l]<- nrow(pp)/(2*n)
        
        pp.out0<- filter(simdata, randomisation==0 & intervention == 0 & outcome==1)
        pp.out1<- filter(simdata, randomisation==1 & intervention == 1 & outcome==1)
        at.out0<- filter(simdata, intervention == 0 & outcome==1)
        at.out1<- filter(simdata, intervention == 1 & outcome==1)
        itt.out0<- filter(simdata, randomisation == 0 & outcome==1)
        itt.out1<- filter(simdata, randomisation == 1 & outcome==1)
        
        .pp.out0<- nrow(pp.out0)/nrow(pp0)
        .pp.out1<- nrow(pp.out1)/nrow(pp1)
        .at.out0<- nrow(at.out0)/nrow(i0)
        .at.out1<- nrow(at.out1)/nrow(i1)
        .itt.out0<- nrow(itt.out0)/n
        .itt.out1<- nrow(itt.out1)/n
        
      }
      
      #proportion of patients who were compliant to protocol 
      i1p[i]<- mean(.i1)
      i0p[i]<- mean(.i0)
      pp1p[i]<- mean(.pp1)
      pp0p[i]<- mean(.pp0)
      ppp [i]<- mean(.pp)
      
      ppout0[i]<- mean(.pp.out0)
      ppout1[i]<- mean(.pp.out1)
      atout0[i]<- mean(.at.out0)
      atout1[i]<- mean(.at.out1)
      ittout0[i]<- mean(.itt.out0)
      ittout1[i]<- mean(.itt.out1)
    }
    return(data.frame(i0p, i1p, pp0p, pp1p, ppp, atout0, atout1,ppout0, ppout1,ittout0, ittout1))
  }
  sim.df<-sim()
  return(sim.df)
} 

simdata(r0=0.2, r1=0.25, p0=0.4, p1=0.3, confounder.eff.i=0.001,confounder.eff.o=5, confounder.u.eff.i=1,confounder.u.eff.o=1)
