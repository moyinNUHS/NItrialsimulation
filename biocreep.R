
################################################################################################################
###################Using causal inference to address non-Adherence in non inferiority trials###################
################################################################################################################

######## Set up #########
rm(list=ls()) # Clean working environment 

# Required libraries 
library(plotrix)
library(dplyr) #for filter
library(gtools) #for logit function
library(parallel) #use multiple cores
library(forcats) #ridge plot 
library(ggridges)
library(ggplot2)
library(gridExtra)

setwd("/Users/moyin/Desktop/angelsfly/NItrialsimulation/codes")

cl<-makeCluster(detectCores()-2)

#set seed for simulation 
set.seed(1234)

#Define parameters 
# n:                  number of participants per group; 
# i0:                 probability of intervention = 1 when Randomization=0
# i1:                 probability of intervention = 1 when Randomization=1
# p1:                 mortality in population when intervention =1
# p0:                 mortality in population when intervention =0 
# confounder.eff.o:   effect of confounder on outcome (odds ratio)
# confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio)
# confounder.eff.i:   effect of confounder on intervention (odds ratio)
# confounder.u.eff.o: effect of unknown confounder on outcome (odds ratio)
# nIterations:        number of iterations 

simdata.bias.2<- function(n, i0, i1, p0, p1, confounder.eff.i,confounder.eff.o, confounder.u.eff.i,confounder.u.eff.o, nIterations, interval){  
  
  #make up vectors for simulations 
  eff.itt2<-.eff.itt2<-c() 
  eff.itt<-.eff.itt<-c() 
  prop.c<-.prop.c<-c()
  
  #alpha error and critical value 
  z = qnorm(1-0.05/2) #alpha=0.025 (one sided)
  
  #true effect 
  true.eff<-p1-p0 
  
  #simulate and derive treatment effect 
    for(i in 1:interval) { print(paste ("interval value",i))
      for(l in 1:nIterations) { 
        tryCatch({ #allow the function to run in case of errors
          #simulate trial data frame 
          
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
          
          b1<-log (confounder.eff.o) 
          b2<-log (confounder.u.eff.o)
          
          if (true.eff>=0) { # ensure all probabilities < 1
            logit.prob.outcome.1<-  b1*simdata$confounder+ b2*simdata$confounder.u
            prob.outcome.1<-exp(logit.prob.outcome.1)/(1+exp(logit.prob.outcome.1))      #probability of outcome if intervention = 1
            simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1) 
            
            prob.outcome.0<-    prob.outcome.1 - true.eff        #probability of outcome if intervention = 0
            simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0)
          } else {
            logit.prob.outcome.0<-  b1*simdata$confounder+ b2*simdata$confounder.u
            prob.outcome.0<-exp(logit.prob.outcome.0)/(1+exp(logit.prob.outcome.0))      #probability of outcome if intervention = 0
            simdata$outcome0 <- rbinom(2*n,1,prob=prob.outcome.0) 
            
            prob.outcome.1<-    prob.outcome.0 + true.eff        #probability of outcome if intervention = 1
            simdata$outcome1 <- rbinom(2*n,1,prob=prob.outcome.1)
          }
          
          #INTERVENTION dependent on randomisation and confounders by logistric regression - creating increasing Adherence
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
          pz1.vector<- simdata[which(simdata$randomisation==1),]$outcome
          pz1.value<- mean(pz1.vector)  
          pz0.vector<- simdata[which(simdata$randomisation==0),]$outcome
          pz0.value<- mean(pz0.vector)
          .eff.itt[l]<- pz1.value-pz0.value    
          var.eff.itt<- pz1.value*(1-pz1.value)/length(pz1.vector) + pz0.value*(1-pz0.value)/length(pz0.vector)
          CI.itt<- .eff.itt[l] + c(-1,1)*z*sqrt(var.eff.itt)
          .eff.itt2[l]<-CI.itt[2]
          
          ## per protocol 
          pp<- rbind(filter(simdata, randomisation==1 & intervention==1), (filter(simdata, randomisation==0 & intervention==0))) # perprotocol population
          .prop.c[l]<- nrow(pp)/(2*n) #proportion of patients who were compliant to protocol
          
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #receive error message if there is an error 
        
      }
      
      # mean proportion of patients who were compliant to protocol 
      prop.c [[i]]<- .prop.c
      eff.itt[[i]]<-.eff.itt
      eff.itt2[[i]]<-.eff.itt2
      
    }
    
    # make a dataframe of the outcomes of the iterations 
    eff.itt.df<- cbind(unlist(prop.c),unlist(eff.itt),unlist(eff.itt2))

    return(eff.itt.df)
  
} 

trial1<-simdata.bias.2(n=500, p1=0.4, p0=0.3, i0=0.45, i1=0.55,confounder.eff.i= 0.001,confounder.eff.o = 1000,  confounder.u.eff.i=1,confounder.u.eff.o=1,nIterations=2000,interval = 50)
trial2<-simdata.bias.2(n=500, p1=0.5, p0=0.4, i0=0.45, i1=0.55,confounder.eff.i= 0.001,confounder.eff.o = 1000,  confounder.u.eff.i=1,confounder.u.eff.o=1,nIterations=2000,interval = 50)
trial3<-simdata.bias.2(n=500, p1=0.6, p0=0.5, i0=0.45, i1=0.55,confounder.eff.i= 0.001,confounder.eff.o = 1000,  confounder.u.eff.i=1,confounder.u.eff.o=1,nIterations=2000,interval = 50)
trial4<-simdata.bias.2(n=500, p1=0.7, p0=0.6, i0=0.45, i1=0.55,confounder.eff.i= 0.001,confounder.eff.o = 1000,  confounder.u.eff.i=1,confounder.u.eff.o=1,nIterations=2000,interval = 50)

stopCluster(cl)

fullset<-list(trial1, trial2, trial3, trial4)
save(fullset, file="fullset.Rdata")

#format data for plot 
load("fullset.RData")

NImargin<-0.1
upperbound70<- 0.72
lowerbound70<-0.68
upperbound80<- 0.82
lowerbound80<-0.78
upperbound90<- 0.92
lowerbound90<-0.88
upperbound100<- 1
lowerbound100<-0.98

subset<- function(fullset) { 
  
  set100<-set90<-set80<-set70<-c()
  
  for (i in 1:4) {
    
    .fullset<-fullset[[i]] #Take every set in trial1, trial2, trial3, trial4 
    
    set100[[i]]<-.fullset[,3][which(.fullset[,1]>lowerbound100 & .fullset[,1]<upperbound100)] #subset those rows with 100% Adherence and end up with a list of data with 100% Adherence from 4 trials
    set90[[i]]<-.fullset[,3][which(.fullset[,1]>lowerbound90 & .fullset[,1]<upperbound90)]
    set80[[i]]<-.fullset[,3][which(.fullset[,1]>lowerbound80 & .fullset[,1]<upperbound80)]
    set70[[i]]<-.fullset[,3][which(.fullset[,1]>lowerbound70 & .fullset[,1]<upperbound70)]
  }
  
  return(list(set100, set90, set80, set70))
}

formatforplot<- function(data) {
  
  forplot<-c()
  
  for (i in 1:4) { #every i is a Adherence level 
  
  set1<-data[[i]][[1]]#trial1
  set2<-data[[i]][[2]]#trial2
  set3<-data[[i]][[3]]#trial3
  set4<-data[[i]][[4]]#trial4
  
  false1<- sum(set1<=0.1)/length(set1)
  false2<- sum(set2<=0.1)/length(set2)*false1
  false3<- sum(set3<=0.1)/length(set3)*false2
  false4<- sum(set4<=0.1)/length(set4)*false3

  set1.df<- as.data.frame(cbind(rep(1, length(set1)), set1))
  colnames(set1.df)<-c('y','x')
  
  .set2.plot<- set2+0.1
  set2.diff<- 0.2-max(.set2.plot[{q<-rank(.set2.plot)/length(.set2.plot);q<false2}])
  set2.plot<- .set2.plot+set2.diff
  set2.df<- as.data.frame(cbind(rep(2, length(set2.plot)), set2.plot))
  colnames(set2.df)<-c('y','x')
  
  .set3.plot<- set3+0.2
  set3.diff<- 0.3-max(.set3.plot[{q<-rank(.set3.plot)/length(.set3.plot);q<false3}])
  if (set3.diff==Inf) {set3.diff=0.3-min(.set3.plot)}
  set3.plot<- .set3.plot+set3.diff
  set3.df<- as.data.frame(cbind(rep(3, length(set3.plot)), set3.plot))
  colnames(set3.df)<-c('y','x')
  
  .set4.plot<- set4+0.3
  set4.diff<- 0.4-max(.set4.plot[{q<-rank(.set4.plot)/length(.set4.plot);q<false4}])
  if (set4.diff==Inf) {set4.diff=0.4-min(.set4.plot)}
  set4.plot<- .set4.plot+set4.diff
  set4.df<- as.data.frame(cbind(rep(4, length(set4.plot)), set4.plot))
  colnames(set4.df)<-c('y','x')
    
  set.df<- as.data.frame(rbind(set1.df, set2.df, set3.df,set4.df))
  set.df$y<-fct_rev(as.factor(set.df$y))
  
  
  
  forplot[[i]]<-set.df
  
  }
  
  return(forplot)
  
}

plot<-function(set.df,Adherence){
  
  gg<- ggplot(set.df, aes(x = x, y = y, group = y),height = ..count..) +
    geom_density_ridges(aes(x = x, fill = paste(y)), 
                        alpha = .8, color = "grey", from = 0, to = 0.55) +
    labs(
      subtitle = paste(Adherence, "%","", "Adherence"),
      y= "", 
      x=""
    ) +
    scale_fill_cyclical(
      values = c("#0868ac", "#43a2ca","#7bccc4","#bae4bc")
    ) +
    scale_x_continuous(limits=c(-0.02,0.55), breaks = c(0,0.1, 0.2, 0.3,0.4)) +
    scale_y_discrete(labels=c("C vs D","B vs C","A vs B", "Standard vs A"))+
    theme_ridges(grid = F) +theme(axis.text=element_text(size=10))
  
  d <- ggplot_build(gg)$data[[1]]
  .d1<- d[which(d$x <=0.1),]
  d1<- .d1[which(.d1$ymin ==4),]
  .d2<- d[which(d$x <=0.2),]
  d2<- .d2[which(.d2$ymin ==3),]
  .d3<- d[which(d$x <=0.3),]
  d3<- .d3[which(.d3$ymin ==2),]
  .d4<- d[which(d$x <=0.4),]
  d4<- .d4[which(.d4$ymin ==1),]
  
  # Add geom_ribbon for shaded area
  plot<- gg +
    geom_ribbon(
      data = transform(d1, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d2, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d3, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d4, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")
  
  return(plot)
}

plot1<-function(set.df,Adherence){
  
  gg<-ggplot(set.df, aes(x = x, y = y, group = y),height = ..count..) +
    geom_density_ridges(aes(x = x, fill = paste(y)), 
                        alpha = .8, color = "grey", from = 0, to = 0.55) +
    labs(
      subtitle = paste(Adherence, "%","", "Adherence"),
      y= "", 
      x="Decrease in treatment efficacy with respect to standard-of-care"
    ) +
    scale_fill_cyclical(
      values = c("#0868ac", "#43a2ca","#7bccc4","#bae4bc")
    ) +
    scale_x_continuous(limits=c(-0.02,0.55), breaks = c(0,0.1, 0.2, 0.3,0.4)) +
    scale_y_discrete(labels=c("C vs D","B vs C","A vs B", "Standard vs A"))+ 
    theme_ridges(grid=F) +theme(axis.text=element_text(size=10), axis.title.x = element_text(hjust =0.5, size=10))
  
  d <- ggplot_build(gg)$data[[1]]
  .d1<- d[which(d$x <=0.1),]
  d1<- .d1[which(.d1$ymin ==4),]
  .d2<- d[which(d$x <=0.2),]
  d2<- .d2[which(.d2$ymin ==3),]
  .d3<- d[which(d$x <=0.3),]
  d3<- .d3[which(.d3$ymin ==2),]
  .d4<- d[which(d$x <=0.4),]
  d4<- .d4[which(.d4$ymin ==1),]
  
  # Add geom_ribbon for shaded area
  plot<- gg +
    geom_ribbon(
      data = transform(d1, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d2, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d3, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")+
    geom_ribbon(
      data = transform(d4, y = group),
      aes(x, ymin = ymin, ymax = ymax, group = group),
      fill = "#df65b0")
  
  return(plot)
  
}


###########################################
###########################################
data<-subset(fullset)
forplot<-formatforplot(data)
plot100<- plot(forplot[[1]],100)
plot90<- plot(forplot[[2]],90)
plot80<- plot(forplot[[3]],80)
plot70<- plot1(forplot[[4]],70)

grid.arrange(plot100, plot90,plot80,plot70,nrow = 4)



